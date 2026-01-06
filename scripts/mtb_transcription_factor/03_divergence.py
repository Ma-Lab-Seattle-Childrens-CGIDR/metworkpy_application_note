"""
Calculate the divergence between each of the TFOE models and a wildtype reference,
for all of the reactions, subsystems, and KEGG pathways
"""

# Setup
# Imports
# Standard Library Imports
from collections import defaultdict
import logging
import pathlib
import warnings
from typing import Any

# External Imports
import cobra  # type:ignore
from metabolic_modeling_utils import kegg_rest
import metworkpy  # type:ignore
from metworkpy.divergence import calculate_divergence_grouped  # type:ignore
from metworkpy.metabolites import metabolite_network  # type:ignore
import numpy as np
import pandas as pd
from tqdm import tqdm  # type:ignore

# Local Imports


# Path setup
try:
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
except NameError:
    BASE_PATH = pathlib.Path(".").absolute()
DATA_PATH = BASE_PATH / "data"
CACHE_PATH = BASE_PATH / "cache"
METABOLITE_NETWORKS_PATH = CACHE_PATH / "metabolite_networks" / "7H9_ADC"
RESULTS_PATH = BASE_PATH / "results" / "transcription_factors"
FLUX_SAMPLES_PATH = CACHE_PATH / "tf_model_flux_samples"
MODELS_PATH = BASE_PATH / "models"

# Make required directories if they don't exist
CACHE_PATH.mkdir(parents=True, exist_ok=True)
METABOLITE_NETWORKS_PATH.mkdir(parents=True, exist_ok=True)
RESULTS_PATH.mkdir(parents=True, exist_ok=True)
FLUX_SAMPLES_PATH.mkdir(parents=True, exist_ok=True)

# Read in the base model for finding subsystems
cobra.Configuration().solver = "hybrid"
BASE_MODEL = metworkpy.read_model(
    MODELS_PATH / "iEK1011_v2_7H9_ADC_glycerol.json"
)

# Run Parameters
PROCESSES = 12
N_NEIGHBORS = 5
METRIC = 2.0
ESSENTIAL_PROPORTION = 0.1
REACTION_PROPORTION = 0.1
PROGRESS_BAR = True

# Setup logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=BASE_PATH
    / "logs"
    / "transcription_factors"
    / "03_divergence.log",
    filemode="w",
    level=logging.INFO,
)

# Determine the metabolite reaction synthesis/consumption networks
# NOTE: The network dataframes are indexed by reactions, with columns for each metabolite

# Create a list of reactions to remove from the metabolite networks
reactions_to_remove = []
subsystems_to_ignore = [
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
if metabolite_synthesis_rxn_net_df_path.exists():
    metabolite_synthesis_rxn_network_df = pd.read_csv(
        metabolite_synthesis_rxn_net_df_path, index_col=0
    )
else:
    # Create the metabolite network
    metabolite_synthesis_rxn_network_df = (
        metabolite_network.find_metabolite_synthesis_network_reactions(
            model=BASE_MODEL.copy(),
            method="essential",
            essential_proportion=ESSENTIAL_PROPORTION,
            processes=PROCESSES,
        )
    )
    # Drop the reactions to ignore
    metabolite_synthesis_rxn_network_df = (
        metabolite_synthesis_rxn_network_df.drop(reactions_to_remove, axis=0)
    )
    # Save the results
    metabolite_synthesis_rxn_network_df.to_csv(
        metabolite_synthesis_rxn_net_df_path, index=True
    )

# Repeat for the consuming networks
metabolite_consumption_rxn_network_df_path = (
    METABOLITE_NETWORKS_PATH / "metabolite_consumption_reaction_network.csv"
)
if metabolite_consumption_rxn_network_df_path.exists():
    metabolite_consumption_rxn_network_df = pd.read_csv(
        metabolite_consumption_rxn_network_df_path, index_col=0
    )
else:
    # NOTE: Indexed by reactions, columns are metabolites
    metabolite_consumption_rxn_network_df = (
        metabolite_network.find_metabolite_consuming_network_reactions(
            BASE_MODEL,
            reaction_proportion=REACTION_PROPORTION,
            progress_bar=False,
            processes=PROCESSES,
        )
    )
    # Drop the reactions being ignored
    metabolite_consumption_rxn_network_df = (
        metabolite_consumption_rxn_network_df.drop(reactions_to_remove, axis=0)
    )
    metabolite_consumption_rxn_network_df.to_csv(
        metabolite_consumption_rxn_network_df_path, index=True
    )
# Get the Kegg pathways
kegg_path_df = kegg_rest.get_kegg_pathway_genes("mtu")
kegg_desc_df = kegg_rest.get_kegg_pathway_descriptions(
    "mtu", remove_str=" - Mycobacterium tuberculosis H37Rv"
)
kegg_path_df = pd.merge(
    kegg_path_df,
    kegg_desc_df,
    how="left",
    left_on="pathway",
    right_on="pathway",
)

# Get a dict of reaction_group:reaction, this will include
# subsystems in the model, metabolite reaction networks, and
# KEGG pathways

# Create a dict to rename errorneous subsystems
subsys_rename_dict = {
    "Arabinogalactan bioynthesis": "Arabinogalactan biosynthesis",
    "Propanoate metabolism": "Propanoate Metabolism",
}

divergence_rxn_group_dict: dict[str, list[Any]] = defaultdict(list)
for rxn in BASE_MODEL.reactions:
    subsys = (
        rxn.subsystem
        if rxn.subsystem not in subsys_rename_dict
        else subsys_rename_dict[rxn.subsystem]
    )
    if subsys and (subsys not in subsystems_to_ignore):
        divergence_rxn_group_dict[f"subsystem__{subsys}"].append(rxn.id)
    if subsys not in subsystems_to_ignore:
        divergence_rxn_group_dict["subsystem__Whole Metabolism"].append(rxn.id)
        # Add a 'subsystem' containing a single rxn for all reactions
        # which are not in excluded subsystems
        divergence_rxn_group_dict[f"reaction__{rxn.id}"].append(rxn.id)
# Add in reaction groups for each metabolite
for metabolite in BASE_MODEL.metabolites:
    divergence_rxn_group_dict[f"metabolite_synthesis__{metabolite.id}"] = (
        sorted(
            set(
                metabolite_synthesis_rxn_network_df[
                    metabolite_synthesis_rxn_network_df[metabolite.id]
                ].index
            )
        )
    )
    divergence_rxn_group_dict[f"metabolite_consumption__{metabolite.id}"] = (
        sorted(
            set(
                metabolite_consumption_rxn_network_df[
                    metabolite_consumption_rxn_network_df[metabolite.id]
                ].index
            )
        )
    )
# Add in reaction groups for each kegg path
model_genes = set(BASE_MODEL.genes.list_attr("id"))
for path_desc, df in kegg_path_df.groupby("description"):
    gene_list = [g for g in list(df["gene"]) if g in model_genes]
    if len(gene_list) == 0:
        continue
    rxn_list = sorted(
        set(
            metworkpy.gene_to_reaction_list(
                model=BASE_MODEL, gene_list=gene_list
            )
        )
    )
    divergence_rxn_group_dict[f"kegg__{path_desc}"] = rxn_list

# Filter the divergence groups down so none are empty
divergence_rxn_group_dict = {
    k: v for k, v in divergence_rxn_group_dict.items() if len(v) > 0
}


sample_list = [
    s.name.split(".")[0] for s in (FLUX_SAMPLES_PATH).glob("*.parquet")
]
sample_list.remove("wildtype")

# Perform the divergence analysis for the diff models
wt_sample = pd.read_parquet(FLUX_SAMPLES_PATH / "wildtype.parquet")

# Create the results dataframe or read it in if it exists
divergence_results_path = RESULTS_PATH / "divergence_results.csv"
divergence_results_df = pd.DataFrame(
    columns=pd.Index(divergence_rxn_group_dict.keys())
)

# Iterate through all the samples, finding the divergence for all the subsystems
for sample in tqdm(sample_list, disable=not PROGRESS_BAR):
    sample_name = sample
    logger.info(f"Starting calculating divergence for {sample_name}")
    # Read in the sample flux distribution
    sample_flux_dist = pd.read_parquet(FLUX_SAMPLES_PATH / f"{sample}.parquet")
    # Find the divergence for all the divergence groups
    with warnings.catch_warnings(action="ignore"):
        div_res = calculate_divergence_grouped(
            wt_sample,
            sample_flux_dist,
            divergence_groups=divergence_rxn_group_dict,
            divergence_type="kl",
            processes=PROCESSES,
            n_neighbors=N_NEIGHBORS,
            distance_metric=METRIC,
            discrete=False,
            jitter=None,
        )
    # Append the results to the output dataframe
    divergence_results_df.loc[sample_name] = div_res

divergence_results_df = divergence_results_df.copy()
# Clip the values in the divergence results
divergence_results_df = divergence_results_df.clip(
    lower=0,
    upper=divergence_results_df.replace(np.inf, np.nan).max(
        axis=0, skipna=True
    ),  # type:ignore # Can be arraylike
    axis="columns",
)
# Save the clipped divergence
divergence_results_df.to_csv(
    divergence_results_path,
)

print("Finished calculating divergence :)")
