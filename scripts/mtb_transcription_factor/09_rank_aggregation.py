"""
Script for performing rank aggregation across a variety of methods
in order to investigate metabolic perturbations caused by the
transcrition factors
"""

# Setup
# Imports
# Standard Library Imports
import pathlib
from typing import Hashable

# External Imports
import cobra  # type:ignore
import numpy as np  # type:ignore
import pandas as pd  # type:ignore
from sklearn.preprocessing import StandardScaler  # type:ignore

# Local Imports
from metabolic_modeling_utils.rra import robust_rank_aggreg  # type:ignore

# Path Setup
try:
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
except NameError:
    BASE_PATH = pathlib.Path(".").absolute()
CACHE_PATH = BASE_PATH / "cache"
DATA_PATH = BASE_PATH / "data"
RESULTS_PATH = BASE_PATH / "results" / "transcription_factors"

# Create Directories if needed
CACHE_PATH.mkdir(parents=True, exist_ok=True)
RESULTS_PATH.mkdir(parents=True, exist_ok=True)
#######################
# Script Parameters ###
#######################
cobra.Configuration().solver = "hybrid"
scaler = StandardScaler()
scaler.set_output(transform="pandas")
PHENOTYPE_RESISTANT_STD_CUTOFF = (
    1.5  # Cutoff for z-score to be considered resistant
)
PHENOTYPE_SUSCEPTIBLE_STD_CUTOFF = (
    -1.5
)  # Cutoff for z-score to be considered susceptible

########################################
# SECTION: Read in Model Information ###
########################################
reaction_information_df = pd.read_csv(
    CACHE_PATH / "model_information" / "reaction_information.csv"
)
metabolite_information_df = pd.read_csv(
    CACHE_PATH / "model_information" / "metabolite_information.csv"
)


###########################
# SECTION: Read in data ###
###########################
# SUBSECTION: Divergence
divergence_df = scaler.fit_transform(
    pd.read_csv(RESULTS_PATH / "divergence_results.csv", index_col=0)
    .replace([np.inf, -np.inf], np.nan)
    .clip(lower=0.0)
)  # Read in and scale the divergence data to account for different sized rxn groups

# Split the divergence df into reaction, subsystem, metabolite_synthesis, and
# metabolite consumptions
subsystem_divergence_df = divergence_df.loc[
    :, divergence_df.columns.str.startswith("subsystem__")
]
subsystem_divergence_df.columns = subsystem_divergence_df.columns.str.replace(
    "^subsystem__", "", regex=True
)
kegg_divergence_df = divergence_df.loc[
    :, divergence_df.columns.str.startswith("kegg__")
]
kegg_divergence_df.columns = kegg_divergence_df.columns.str.replace(
    "^kegg__", "", regex=True
)
metabolite_synthesis_divergence_df = divergence_df.loc[
    :, divergence_df.columns.str.startswith("metabolite_synthesis__")
]
metabolite_synthesis_divergence_df.columns = (
    metabolite_synthesis_divergence_df.columns.str.replace(
        "^metabolite_synthesis__", "", regex=True
    )
)
metabolite_consumption_divergence_df = divergence_df.loc[
    :, divergence_df.columns.str.startswith("metabolite_consumption__")
]
metabolite_consumption_divergence_df.columns = (
    metabolite_consumption_divergence_df.columns.str.replace(
        "^metabolite_consumption__", "", regex=True
    )
)
reaction_divergence_df = divergence_df.loc[
    :, divergence_df.columns.str.startswith("reaction__")
]
reaction_divergence_df.columns = reaction_divergence_df.columns.str.replace(
    "^reaction__", "", regex=True
)

# SUBSECTION: GSVA
gsva_df = pd.read_csv(RESULTS_PATH / "metabolite_gsva.csv", index_col=0)
# Split the GSVA into consumption and synthesis
gsva_synthesis_df = gsva_df.loc[
    :, gsva_df.columns.str.endswith("_synthesis_network")
]
gsva_synthesis_df.columns = gsva_synthesis_df.columns.str.replace(
    "_synthesis_network$", "", regex=True
)
gsva_consumption_df = gsva_df.loc[
    :, gsva_df.columns.str.endswith("_consumption_network")
]
gsva_consumption_df.columns = gsva_consumption_df.columns.str.replace(
    "_consumption_network$", "", regex=True
)

# SUBSECTION: Density
density_df = pd.read_csv(RESULTS_PATH / "tf_target_density.csv", index_col=0).T

# SUBSECTION: Metabolite Enrichment
enrichment_df = (
    pd.read_csv(
        RESULTS_PATH / "tf_target_metabolite_network_enrichment.csv",
        usecols=[
            "tf",
            "metabolite",
            "metabolite network direction",
            "p-value",
        ],
    )
    .query("`metabolite network direction` == 'synthesis'")
    .pivot(index="metabolite", columns="tf", values="p-value")
)


# SUBSECTION: Drug Phenotype
# NOTE: These are the diff, not the z-scores
inh_broth_df = pd.read_excel(
    DATA_PATH / "transcription_factors" / "trip_phenotype_inh.xlsx",
    sheet_name="TableS1-TRIPAbundanceFoldChange",
    index_col=0,
    usecols="A,D,E",
    skiprows=3,
    names=["TF", "Broth", "INH"],  # type: ignore
)
inh_broth_df["INH_DIFF"] = inh_broth_df["INH"] - inh_broth_df["Broth"]


# type: ignore
other_drug_df = pd.read_excel(  # type: ignore
    DATA_PATH / "transcription_factors" / "all_drug.xlsx",
    sheet_name="All",
    index_col=0,
    usecols="A,L:O",
    names=["TF", "BDQ_DIFF", "OFX_DIFF", "PBTZ_DIFF", "RIF_DIFF"],  # type: ignore
)
phenotype_df = (
    pd.merge(
        other_drug_df,
        inh_broth_df[["INH_DIFF"]],
        how="inner",
        left_index=True,
        right_index=True,
    )
    .sort_index()
    .drop("Empty_Plasmid", axis=0)
)


#####################################################
# SECTION: Synthesis/Consumption Rank Aggregation ###
#####################################################
# NOTE: For each TF, perform rank aggregation between the synthesis and consumption
# metabolite networks

# SUBSECTION: Divergence
divergence_synth_cons_df_list: list[pd.DataFrame] = []
num_metabolites = max(
    len(metabolite_synthesis_divergence_df.columns),
    len(metabolite_consumption_divergence_df.columns),
)
for tf in divergence_df.index:
    # Extract the synthesis and consumption series for this TF
    synth_series = metabolite_synthesis_divergence_df.loc[tf]
    cons_series = metabolite_consumption_divergence_df.loc[tf]
    # Convert these into the rank dictionary
    div_rank_dict = {
        "synthesis": list(synth_series.sort_values(ascending=False).index),
        "consumption": list(cons_series.sort_values(ascending=False).index),
    }
    # Perform the rank aggregation
    rra_res = (
        robust_rank_aggreg(
            rank_list_dict=div_rank_dict,
            num_ranked_elems=num_metabolites,
            full=True,
            exact=True,
        )
        .reset_index(drop=True)
        .reset_index(drop=False, names="rank")
    )
    # Add a column describing which TF this was performed for
    rra_res["tf"] = tf
    # Append to the results list
    divergence_synth_cons_df_list.append(rra_res)
# Combine the results for all the TFs
divergence_synth_cons_df = pd.concat(
    divergence_synth_cons_df_list, axis=0, ignore_index=True
)
# Add in the Metabolite Information
divergence_synth_cons_df = pd.merge(
    divergence_synth_cons_df,
    metabolite_information_df,
    how="left",
    left_on="Name",
    right_on="id",
)
# Save the results
divergence_synth_cons_df.to_csv(
    RESULTS_PATH / "metabolite_synthesis_consumption_divergence_rank_agg.csv",
    index=False,
)

# SUBSECTION: GSVA
gsva_synth_cons_df_list: list[pd.DataFrame] = []
num_metabolites = max(
    len(gsva_synthesis_df.columns),
    len(gsva_consumption_df.columns),
)
for tf in gsva_df.index:
    # Extract the synthesis and consumption series for this TF
    synth_series = gsva_synthesis_df.loc[tf]
    cons_series = gsva_consumption_df.loc[tf]
    # Convert these into the rank dictionary
    gsva_rank_dict = {
        "synthesis": list(synth_series.sort_values(ascending=False).index),
        "consumption": list(cons_series.sort_values(ascending=False).index),
    }
    # Perform the rank aggregation
    rra_res = (
        robust_rank_aggreg(
            rank_list_dict=gsva_rank_dict,
            num_ranked_elems=num_metabolites,
            full=True,
            exact=True,
        )
        .reset_index(drop=True)
        .reset_index(drop=False, names="rank")
    )
    # Add a column describing which TF this was performed for
    rra_res["tf"] = tf
    # Append to the results list
    gsva_synth_cons_df_list.append(rra_res)
# Combine the results for all the TFs
gsva_synth_cons_df = pd.concat(
    gsva_synth_cons_df_list, axis=0, ignore_index=True
)
# Add in the Metabolite Information
gsva_synth_cons_df = pd.merge(
    gsva_synth_cons_df,
    metabolite_information_df,
    how="left",
    left_on="Name",
    right_on="id",
)
# Save the results
gsva_synth_cons_df.to_csv(
    RESULTS_PATH / "metabolite_synthesis_consumption_gsva_rank_agg.csv",
    index=False,
)

################################################
## SECTION: Rank Aggregation for Metabolites ###
################################################

# For each TF, aggregate across the different methods for
# evaluating metabolite perturbations
# - GSVA (2-directional?)
# - Enrichment
# - Divergence

# GSVA synthesis df index is tfs
# Enrichment Df columns are tf
# Divergence metabolite synthesis divergence index is tfs

tf_list = sorted(
    set(gsva_synthesis_df.index)
    & set(enrichment_df.columns)
    & set(metabolite_synthesis_divergence_df.index)
)

num_metabolites = max(
    len(set(gsva_synthesis_df.columns)),
    len(set(enrichment_df.index)),
    len(set(metabolite_synthesis_divergence_df.columns)),
)

# NOTE: Filtering down for only enrichment with p-values below 0.1

metabolite_agg_df_list: list[pd.DataFrame] = []

for tf in tf_list:
    rank_dict = {
        "gsva_synth_descending": list(
            gsva_synthesis_df.loc[tf].sort_values(ascending=False).index
        ),
        "gsva_synth_ascending": list(
            gsva_synthesis_df.loc[tf].sort_values(ascending=True).index
        ),
        "enrichment": list(
            enrichment_df[tf]
            .loc[lambda x: x <= 0.10]
            .sort_values(ascending=True)
            .index
        ),
        "divergence": list(
            metabolite_synthesis_divergence_df.loc[tf]
            .sort_values(ascending=False)
            .index
        ),
    }
    rank_agg_df = (
        robust_rank_aggreg(
            rank_list_dict=rank_dict,
            num_ranked_elems=num_metabolites,
            full=False,
            exact=True,
        )
        .reset_index(drop=True)
        .reset_index(drop=False, names="rank")
    )
    rank_agg_df["tf"] = tf
    metabolite_agg_df_list.append(rank_agg_df)

# Combine the results
metabolite_agg_df = pd.concat(metabolite_agg_df_list, axis=0)

# Save these results after adding in metabolite information
metabolite_agg_df.merge(
    metabolite_information_df, left_on="Name", right_on="id", how="left"
).to_csv(RESULTS_PATH / "tf_metabolite_rank_aggregation.csv", index=False)

#############################################
# SECTION: Rank Aggregation for Reactions ###
#############################################

# Divergence
# Density

tf_list = sorted(set(reaction_divergence_df.index) & set(density_df.index))
rxn_count = max(
    len(set(reaction_divergence_df.columns)), len(set(density_df.columns))
)

reaction_agg_df_list: list[pd.DataFrame] = []

for tf in tf_list:
    rank_dict = {
        "density": list(
            density_df.loc[tf]
            .loc[lambda x: x > 0.0]
            .sort_values(ascending=False)
            .index
        ),
        "divergence": list(
            reaction_divergence_df.loc[tf]
            .dropna()
            .sort_values(ascending=False)
            .index
        ),
    }
    rank_agg_df = (
        robust_rank_aggreg(
            rank_list_dict=rank_dict,
            num_ranked_elems=rxn_count,
            full=False,
            exact=True,
        )
        .reset_index(drop=True)
        .reset_index(drop=False, names="rank")
    )
    rank_agg_df["tf"] = tf
    reaction_agg_df_list.append(rank_agg_df)
# Combine the results
reaction_agg_df = pd.concat(reaction_agg_df_list, axis=0)

# Save the results after adding the reaction information
reaction_agg_df.merge(
    reaction_information_df, left_on="Name", right_on="id", how="left"
).to_csv(RESULTS_PATH / "tf_reaction_rank_aggregation.csv", index=False)


###################################################
# SECTION: Rank Aggregation for Drug Phenotypes ###
###################################################
# NOTE: Perform aggregation for the various metrics across
# the resistant and susceptible strains


# Helper function for converting a dataframe to a list of ranks
def df_to_rank_dict(
    input_df: pd.DataFrame, ascending: bool = False, axis: str | int = 0
) -> dict[Hashable, list[Hashable]]:
    """
    Function to convert a dataframe into a dictionary of rank lists

    Parameters
    ----------
    input_df : pd.DataFrame
        Dataframe to convert to a dict of rank lists
    ascending : bool
        Whether the values should be sorted in ascending order or not
    axis : str or int
        Which axis should be used as the keys of the dictionary
        (columns or 1 will use the column headers as keys for example)

    Returns
    -------
    dict of Hashable to lists of Hashable
        The dict of rank lists
    """
    match axis:
        case 0 | "index":
            df_iter = input_df.iterrows()
        case 1 | "columns":
            df_iter = input_df.items()
        case _:
            raise ValueError("Invalid axis specifier")
    rank_dict = {}
    for key, value_series in df_iter:
        rank_dict[key] = list(
            value_series.dropna().sort_values(ascending=ascending).index
        )
    return rank_dict


# Create dataframes from the aggregated metabolites and TFs
rxn_agg_df = reaction_agg_df.pivot(index="Name", columns="tf", values="Score")
met_agg_df = metabolite_agg_df.pivot(
    index="Name", columns="tf", values="Score"
)


drug_phenotype_rra_metabolite_res_list: list[pd.DataFrame] = []
drug_phenotype_rra_reaction_res_list: list[pd.DataFrame] = []
drug_phenotype_rra_subsystem_res_list: list[pd.DataFrame] = []
drug_phenotype_rra_kegg_res_list: list[pd.DataFrame] = []
for drug, phenotype_series in phenotype_df.items():
    for phenotype_direction in ["resistant", "susceptible"]:
        if phenotype_direction == "resistant":
            tf_set = set(
                phenotype_series[
                    phenotype_series >= PHENOTYPE_RESISTANT_STD_CUTOFF
                ].index
            )
        elif phenotype_direction == "susceptible":
            tf_set = set(
                phenotype_series[
                    phenotype_series <= PHENOTYPE_SUSCEPTIBLE_STD_CUTOFF
                ].index
            )
        else:
            assert False, (
                "Iterator over known list resulted in unknown element...somehow"
            )
        # SUBSECTION: Metabolite Network Divergence
        for metabolite_network_direction, met_net_div_df in zip(
            ["synthesis", "consumption"],
            [
                metabolite_synthesis_divergence_df,
                metabolite_consumption_divergence_df,
            ],
        ):
            met_net_div_tf_list = list(
                set(met_net_div_df.index).intersection(tf_set)
            )
            if len(met_net_div_tf_list) <= 3:
                continue
            met_net_div_rank_dict = df_to_rank_dict(
                met_net_div_df.loc[met_net_div_tf_list],
                ascending=False,
                axis="index",
            )
            met_net_div_rra_res = (
                robust_rank_aggreg(
                    rank_list_dict=met_net_div_rank_dict,  # type:ignore
                    num_ranked_elems=len(met_net_div_df.columns),
                    full=True,
                    exact=True,
                )
                .reset_index(drop=True)
                .reset_index(drop=False, names="rank")
            )
            met_net_div_rra_res["drug"] = drug
            met_net_div_rra_res["phenotype direction"] = phenotype_direction
            met_net_div_rra_res["method"] = (
                f"divergence metabolite {metabolite_network_direction}"
            )
            drug_phenotype_rra_metabolite_res_list.append(met_net_div_rra_res)

        # SUBSECTION: Aggregation for GSVA metabolite results
        for metabolite_network_direction, met_net_gsva_df in zip(
            ["synthesis", "consumption"],
            [gsva_synthesis_df, gsva_consumption_df],
        ):
            for sort_direction in ["ascending", "descending"]:
                gsva_tf_list = list(
                    set(met_net_gsva_df.index).intersection(tf_set)
                )
                if len(gsva_tf_list) <= 5:
                    continue
                met_net_gsva_rank_dict = df_to_rank_dict(
                    input_df=met_net_gsva_df.loc[gsva_tf_list],
                    ascending=(sort_direction == "ascending"),
                    axis="index",
                )
                met_net_gsva_rra_res = (
                    robust_rank_aggreg(
                        rank_list_dict=met_net_gsva_rank_dict,  # type:ignore
                        num_ranked_elems=len(met_net_gsva_df.columns),
                        full=True,
                        exact=True,
                    )
                    .reset_index(drop=True)
                    .reset_index(drop=False, names="rank")
                )
                met_net_gsva_rra_res["drug"] = drug
                met_net_gsva_rra_res["phenotype direction"] = (
                    phenotype_direction
                )
                met_net_gsva_rra_res["method"] = (
                    f"GSVA metabolite {metabolite_network_direction} {sort_direction}"
                )
                drug_phenotype_rra_metabolite_res_list.append(
                    met_net_gsva_rra_res
                )
        # SUBSECTION: Aggregation for Reaction Divergence
        rxn_div_tf_list = list(
            set(reaction_divergence_df.index).intersection(tf_set)
        )
        rxn_div_rank_dict = df_to_rank_dict(
            input_df=reaction_divergence_df.loc[rxn_div_tf_list],
            ascending=False,
            axis="index",
        )
        rxn_div_rra_res = (
            robust_rank_aggreg(
                rank_list_dict=rxn_div_rank_dict,  # type:ignore
                num_ranked_elems=len(
                    reaction_divergence_df.dropna(axis=1, how="all").columns
                ),
            )
            .reset_index(drop=True)
            .reset_index(drop=False, names="rank")
        )
        rxn_div_rra_res["drug"] = drug
        rxn_div_rra_res["phenotype direction"] = phenotype_direction
        rxn_div_rra_res["method"] = "reaction divergence"
        drug_phenotype_rra_reaction_res_list.append(rxn_div_rra_res)
        # SUBSECTION: Aggregation for reaction density
        density_tf_list = list(set(density_df.index).intersection(tf_set))
        density_rank_dict = {}
        for tf in density_tf_list:
            density_series = density_df.loc[tf]
            # Filter out any that are exactly 0
            density_series = density_series[density_series > 0.0]
            # add to the rank dict
            density_rank_dict[tf] = list(
                density_series.sort_values(ascending=False).index
            )
        density_rra_res = (
            robust_rank_aggreg(
                rank_list_dict=density_rank_dict,
                num_ranked_elems=len(density_df.columns),
                full=False,
                exact=True,
            )
            .reset_index(drop=True)
            .reset_index(drop=False, names="rank")
        )
        density_rra_res["drug"] = drug
        density_rra_res["phenotype direction"] = phenotype_direction
        density_rra_res["method"] = "target density"
        drug_phenotype_rra_reaction_res_list.append(density_rra_res)
        # SUBSECTION: Subsystem Divergence
        subsystem_div_tf_list = list(
            set(subsystem_divergence_df.index).intersection(tf_set)
        )
        subsystem_rank_dict = df_to_rank_dict(
            subsystem_divergence_df.loc[subsystem_div_tf_list],
            ascending=False,
            axis="index",
        )
        subsystem_rra_res = (
            robust_rank_aggreg(
                rank_list_dict=subsystem_rank_dict,  # type: ignore
                num_ranked_elems=len(subsystem_divergence_df.columns),
                full=True,
                exact=True,
            )
            .reset_index(drop=True)
            .reset_index(drop=False, names="rank")
        )
        subsystem_rra_res["drug"] = drug
        subsystem_rra_res["phenotype direction"] = phenotype_direction
        drug_phenotype_rra_subsystem_res_list.append(subsystem_rra_res)
        # SUBSECTION: KEGG Divergence
        kegg_div_tf_list = list(
            set(kegg_divergence_df.index).intersection(tf_set)
        )
        kegg_rank_dict = df_to_rank_dict(
            kegg_divergence_df.loc[kegg_div_tf_list],
            ascending=False,
            axis="index",
        )
        kegg_rra_res = (
            robust_rank_aggreg(
                rank_list_dict=kegg_rank_dict,  # type: ignore
                num_ranked_elems=len(kegg_divergence_df.columns),
                full=True,
                exact=True,
            )
            .reset_index(drop=True)
            .reset_index(drop=False, names="rank")
        )
        kegg_rra_res["drug"] = drug
        kegg_rra_res["phenotype direction"] = phenotype_direction
        drug_phenotype_rra_kegg_res_list.append(kegg_rra_res)
        # SUBSECTION: Aggregate TF reaction results
        rxn_agg_tf_list = list(set(rxn_agg_df.columns).intersection(tf_set))
        rxn_agg_rank_dict = df_to_rank_dict(
            input_df=rxn_agg_df[rxn_agg_tf_list],
            ascending=True,
            axis="columns",
        )
        rxn_agg_rra_res = (
            robust_rank_aggreg(
                rank_list_dict=rxn_agg_rank_dict,  # type:ignore
                num_ranked_elems=len(rxn_agg_df.index),
                full=False,
            )
            .reset_index(drop=True)
            .reset_index(drop=False, names="rank")
        )
        rxn_agg_rra_res["drug"] = drug
        rxn_agg_rra_res["phenotype direction"] = phenotype_direction
        rxn_agg_rra_res["method"] = "aggregated"
        drug_phenotype_rra_reaction_res_list.append(rxn_agg_rra_res)
        # SUBSECTION: Aggregate TF metabolite results
        met_agg_tf_list = list(set(met_agg_df.columns).intersection(tf_set))
        met_agg_rank_dict = df_to_rank_dict(
            input_df=met_agg_df[met_agg_tf_list],
            ascending=True,
            axis="columns",
        )
        met_agg_rra_res = (
            robust_rank_aggreg(
                rank_list_dict=met_agg_rank_dict,  # type:ignore
                num_ranked_elems=len(met_agg_df.index),
                full=False,
            )
            .reset_index(drop=True)
            .reset_index(drop=False, names="rank")
        )
        met_agg_rra_res["drug"] = drug
        met_agg_rra_res["phenotype direction"] = phenotype_direction
        met_agg_rra_res["method"] = "aggregated"
        drug_phenotype_rra_metabolite_res_list.append(met_agg_rra_res)


# SUBSECTION: Combine together results dataframes, add in reaction and metabolite Information
metabolite_rra_res = pd.merge(
    pd.concat(
        drug_phenotype_rra_metabolite_res_list, axis=0, ignore_index=True
    ),
    metabolite_information_df,
    how="left",
    left_on="Name",
    right_on="id",
)
reaction_rra_res = pd.merge(
    pd.concat(drug_phenotype_rra_reaction_res_list, axis=0, ignore_index=True),
    reaction_information_df,
    how="left",
    left_on="Name",
    right_on="id",
)
subsystem_rra_res = pd.concat(
    drug_phenotype_rra_subsystem_res_list, axis=0, ignore_index=True
)
kegg_rra_res = pd.concat(
    drug_phenotype_rra_kegg_res_list, axis=0, ignore_index=True
)

# SUBSECTION: Save the results
metabolite_rra_res.to_csv(RESULTS_PATH / "metabolite_drug_phenotype_rra.csv")
reaction_rra_res.to_csv(RESULTS_PATH / "reaction_drug_phenotype_rra.csv")
subsystem_rra_res.to_csv(RESULTS_PATH / "subsystem_drug_phenotype_rra.csv")
kegg_rra_res.to_csv(RESULTS_PATH / "kegg_drug_phenotype_rra.csv")
