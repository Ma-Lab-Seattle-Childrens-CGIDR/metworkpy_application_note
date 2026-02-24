"""
Script to evaluate the impact TF gene targets
on the metabolic network using ko-divergence
"""

# Setup
# Imports
# Standard Library Imports
import logging
import pathlib
import sys
import tomllib
from typing import cast
import warnings

# External Imports
import cobra
import metworkpy
from metabolic_modeling_utils.false_discovery_control import fdr_with_nan
import numpy as np
import pandas as pd
from scipy import stats

# Local Imports
from common_functions import get_metabolite_network

# Path Setup
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".").absolute()  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
DATA_PATH = BASE_PATH / "data"
RESULTS_PATH = BASE_PATH / "results" / "mtb_transcription_factors"
CACHE_PATH = BASE_PATH / "cache"
METABOLITE_NETWORKS_PATH = CACHE_PATH / "metabolite_networks" / "7H9_ADC"
GENE_KO_DIVERGENCE_PATH = CACHE_PATH / "gene_ko_divergence" / "7h9_adc"
BASE_MODEL_PATH = BASE_PATH / "models" / "iEK1011_v2_7H9_ADC_glycerol.json"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"

# Create Directories if needed
CACHE_PATH.mkdir(parents=True, exist_ok=True)
GENE_KO_DIVERGENCE_PATH.mkdir(parents=True, exist_ok=True)
METABOLITE_NETWORKS_PATH.mkdir(parents=True, exist_ok=True)
RESULTS_PATH.mkdir(parents=True, exist_ok=True)
LOG_PATH.mkdir(parents=True, exist_ok=True)


logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=LOG_PATH / "04_tf_ko_divergence.log",
    filemode="w",
    level=logging.INFO,
)

# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)

# Script Parameters
cobra.Configuration().solver = CONFIG["cobra"]["solver"]

# Read in the Base Model
logger.info("Reading in the base model")
BASE_MODEL = metworkpy.read_model(BASE_MODEL_PATH)

# Create a list of reactions to remove from the metabolite networks
reactions_to_remove: list[str] = []
subsystems_to_ignore: list[str] = [
    "Intracellular demand",
    "Biomass and maintenance functions",
    "Extracellular exchange",
]
for rxn in BASE_MODEL.reactions:
    if rxn.subsystem in subsystems_to_ignore:
        reactions_to_remove.append(rxn.id)
metabolite_synthesis_rxn_net_df_path = (
    METABOLITE_NETWORKS_PATH / "metabolite_synthesis_reaction_network.csv"
)
metabolite_synthesis_rxn_network_df = get_metabolite_network(
    out_path=metabolite_synthesis_rxn_net_df_path,
    model=BASE_MODEL,
    rxns_to_remove=reactions_to_remove,
    proportion=CONFIG["mtb_tf"]["divergence"]["essential-proportion"],
    synthesis=True,
    processes=CONFIG["processes"],
)
# Start by dropping any columns that are all 0
metabolite_synthesis_rxn_network_df = metabolite_synthesis_rxn_network_df.drop(
    metabolite_synthesis_rxn_network_df.loc[
        :, metabolite_synthesis_rxn_network_df.sum() == 0
    ].columns,
    axis=1,
)

# Find correlation between the metabolite networks
# and create a list of representative metabolites
corr_df = metabolite_synthesis_rxn_network_df.corr()
# Create an edge list for finding connected components
edge_list: list[tuple[str, str]] = (
    corr_df[corr_df >= 0.75].stack().index.to_list()
)
node_list = corr_df.columns.tolist()
# Find the representative nodes (with the representative node being the one with the highest degree)
representative_metabolite_dict: dict[str, set[str]] = {
    k: v - {k}
    for k, v in metworkpy.utils.find_representative_nodes(
        node_list=node_list, edge_list=edge_list
    ).items()
}

metabolites_to_drop: set[str] = set()
for mets in representative_metabolite_dict.values():
    metabolites_to_drop |= mets
metabolite_synthesis_rxn_network_df = metabolite_synthesis_rxn_network_df.drop(
    list(metabolites_to_drop), axis=1
)
# Create a dict of metabolite: reaction set for each network still in the dataframe
metabolite_network_dict: dict[str, list[str]] = {}
for met, rxns in metabolite_synthesis_rxn_network_df.items():
    metabolite_network_dict[met] = list(rxns[rxns].index)  # type: ignore

# Create a dict of reactions to find the divergence for
divergence_targets: dict[str, list[str]] = {}

subsystems_to_ignore_reactions: list[str] = [
    "Intracellular demand",
    "Extracellular exchange",
]

for reaction in BASE_MODEL.reactions:
    if reaction.subsystem in subsystems_to_ignore:
        continue
    divergence_targets[f"{reaction.id}__reaction"]

# Add in the metabolite networks
for met, net in metabolite_network_dict.items():
    if len(net) == 0:
        warnings.warn(
            f"Found empty network for metabolite: {met}"
        )  # This will print to slurm output
        logger.warning(f"Found empty network for metabolite: {met}")
        continue
    divergence_targets[f"{met}__metabolite"] = net

# Perform the KO divergence calculations
logger.info("Reading in or generating the ko divergence results")
ko_divergence_df_path = (
    GENE_KO_DIVERGENCE_PATH / "gene_ko_divergence_results.csv"
)
if ko_divergence_df_path.exists():
    logger.info("Reading in the ko divergence dataframe")
    ko_divergence_df = pd.read_csv(ko_divergence_df_path, index_col=0)
else:
    logger.info("Performing ko divergence analysis")
    ko_divergence_df = cast(
        pd.DataFrame,
        metworkpy.divergence.ko_divergence(
            model=BASE_MODEL,
            genes_to_ko=BASE_MODEL.genes.list_attr("id"),
            target_networks=divergence_targets,
            divergence_type="kl",
            n_neighbors=CONFIG["mtb_tf"]["ko_divergence"]["n-neighbors"],
            sample_count=CONFIG["mtb_tf"]["ko_divergence"]["sample-count"],
            processes=CONFIG["processes"],
            sampler_seed=1618,
        ),
    ).clip(lower=0.0)
    logger.info("Saving the ko divergence results")
    ko_divergence_df.to_csv(ko_divergence_df_path, index=True)
# Filter out columns that are all nan
ko_divergence_df = ko_divergence_df.dropna(axis="columns")

# Get a list of genes in the model
model_gene_set: set[str] = set(BASE_MODEL.genes.list_attr("id"))

# Find the TF targets
logger.info("Finding the TF targets")
tf_fc_df = pd.read_excel(
    DATA_PATH / "mtb_transcription_factors" / "tfoe_targets.xlsx",
    sheet_name="SupplementaryTableS2",
    skiprows=list(range(8)) + [9],
    usecols="A,E:HB",
    index_col=0,
)
tf_pval_df = pd.read_excel(
    DATA_PATH / "mtb_transcription_factors" / "tfoe_targets.xlsx",
    sheet_name="SupplementaryTableS2",
    skiprows=list(range(8)) + [9],
    usecols="A,HC:OZ",
    index_col=0,
)
tf_pval_df.columns = tf_pval_df.columns.str.replace(".1", "")

tf_target_df = (
    tf_fc_df.abs() >= CONFIG["mtb_tf"]["ko_divergence"]["target-fc-cutoff"]
) & (tf_pval_df <= CONFIG["mtb_tf"]["ko_divergence"]["target-pval-cutoff"])

# Create a dictionary of TF: reaction targets
logger.info("Creating a reaction target dictionary for all the TFs")
tf_target_dict: dict[str, list[str]] = {}
for tf, target_series in tf_target_df.items():
    tf_target_dict[str(tf)] = [
        str(g)
        for g in target_series[target_series].index
        if g in model_gene_set
    ]

# Filter for only TFs which have at least 5 targets in the model
tf_target_dict = {
    k: v
    for k, v in tf_target_dict.items()
    if len(v) >= CONFIG["mtb_tf"]["ko_divergence"]["min-target-count"]
}

# Now, for each TF test if it targets genes which cause
# relatively large perturbations in each of the various metabolite synthesis networks
logger.info("Performing tests for TF target ko-divergence")
results_df_list: list[pd.DataFrame] = []

# Filter the ko_divergence_df for only metabolites for the TF testing
ko_divergence_df = ko_divergence_df.loc[
    :, ko_divergence_df.columns.str.endswith("__metabolite")
]
ko_divergence_df.columns = ko_divergence_df.columns.str.replace(
    "__metabolite$", "", regex=True
)

for tf, target_list in tf_target_dict.items():
    logger.info(f"Starting tests for TF: {tf}")
    tf_res_df = pd.DataFrame(
        np.nan,
        index=ko_divergence_df.columns,
        columns=pd.Index(
            [
                "Mann-Whitney U1",
                "Mann-Whitney U2",
                "rho",
                "p-value",
                "adj p-value",
            ]
        ),
    )
    for metabolite, ko_divergence_series in ko_divergence_df.items():
        # Drop na values and infs
        ko_divergence_series = ko_divergence_series.replace(
            [np.inf, -np.inf], np.nan
        ).dropna()
        # Seperate the targeted and non-targeted genes
        targeted_divergence = ko_divergence_series[
            ko_divergence_series.index.isin(target_list)
        ]
        non_targeted_divergence = ko_divergence_series[
            ~ko_divergence_series.index.isin(target_list)
        ]
        if (
            len(targeted_divergence)
            < CONFIG["mtb_tf"]["ko_divergence"]["min-target-count"]
        ) or (
            len(non_targeted_divergence)
            < CONFIG["mtb_tf"]["ko_divergence"]["min-target-count"]
        ):
            # Don't calculate for TFs which target too many, or too few genes
            continue
        # Calculate the results
        # Calculate the Mann-Whitney test results
        mann_whitney_res = stats.mannwhitneyu(
            x=targeted_divergence,
            y=non_targeted_divergence,
            alternative="greater",
        )
        u1 = mann_whitney_res.statistic
        # Calculate u2 based on Scipy Stats documnetation
        u2 = len(targeted_divergence) * len(non_targeted_divergence) - u1
        # Calculate the rho (aka AUC)
        rho = u1 / (len(targeted_divergence) * len(non_targeted_divergence))
        pval = mann_whitney_res.pvalue
        tf_res_df.loc[metabolite, "Mann-Whitney U1"] = u1
        tf_res_df.loc[metabolite, "Mann-Whitney U2"] = u2
        tf_res_df.loc[metabolite, "rho"] = rho
        tf_res_df.loc[metabolite, "p-value"] = pval
    logger.info(f"Finished performing tests for {tf}")
    # Correct for false discovery rate
    tf_res_df["adj p-value"] = fdr_with_nan(tf_res_df["p-value"])
    # Reset the index to get a metabolite column
    tf_res_df = tf_res_df.reset_index(drop=False, names="metabolite")
    tf_res_df["tf"] = tf
    # Add the results dataframe to the overall results
    results_df_list.append(tf_res_df)

# Combine the results for all the TFs
logger.info("Finished performing results, combining across TFs")
tf_ko_divergence_res_df = pd.concat(results_df_list, axis=0)

# Add in a column describing which metabolites are represented by each representative
represented_met_df = pd.Series(
    {k: str(v) for k, v in representative_metabolite_dict.items()}
).to_frame(name="represented metabolites")
tf_ko_divergence_res_df = tf_ko_divergence_res_df.merge(
    right=represented_met_df,
    left_on="metabolite",
    right_index=True,
    how="left",
)


# Add in metabolite information
metabolite_info_df = pd.read_csv(
    CACHE_PATH / "model_information" / "metabolite_information.csv"
)

tf_ko_divergence_res_df = pd.merge(
    tf_ko_divergence_res_df,
    metabolite_info_df,
    how="left",
    left_on="metabolite",
    right_on="id",
)
# Save the final results
logger.info("Saving the final results")
tf_ko_divergence_res_df.to_csv(
    RESULTS_PATH / "ko_divergence_tf_target_analysis.csv",
    index=False,
)
logger.info("Finished performing KO-divergence analysis for the TFs! ;)")
