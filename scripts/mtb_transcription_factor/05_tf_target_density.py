"""
Script for calculating the target density of the transcription
factor targets in the metabolic model
"""

# Setup
# Imports
# Standard Library Imports
from collections import defaultdict
import logging
import pathlib
import sys
import tomllib

# External Imports
import cobra  # type: ignore
import metworkpy  # type:ignore
import pandas as pd
import scipy.stats as stats

# Local Imports

# Path Setup
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".").absolute()  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
DATA_PATH = BASE_PATH / "data"
CACHE_PATH = BASE_PATH / "cache"
RESULTS_PATH = BASE_PATH / "results" / "mtb_transcription_factors"
BASE_MODEL_PATH = BASE_PATH / "models" / "iEK1011_v2_7H9_ADC_glycerol.json"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"

# Create Directories if needed
CACHE_PATH.mkdir(parents=True, exist_ok=True)
RESULTS_PATH.mkdir(parents=True, exist_ok=True)
LOG_PATH.mkdir(parents=True, exist_ok=True)

# Logging Config
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=LOG_PATH / "05_tf_target_density.log",
    filemode="w",
    level=logging.INFO,
)

# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)

# Read in the metabolic model
logger.info("Reading in the base model")
cobra.Configuration().solver = CONFIG["cobra"]["solver"]
BASE_MODEL = metworkpy.read_model(BASE_MODEL_PATH)

# SECTION: Construct the reaction network
# Find reactions to exlude
logger.info("Finding reactions to exlude")
reactions_to_exclude = []
SUBSYSTEMS_TO_EXCLUDE = [
    "Biomass and maintenance functions",
    "Intracellular demand",
    "Extracellular exchange",
]
metabolite_rxn_count_dict: dict[str, int] = defaultdict(int)
model_reactions = set()
for rxn in BASE_MODEL.reactions:
    if rxn.subsystem in SUBSYSTEMS_TO_EXCLUDE:
        reactions_to_exclude.append(rxn.id)
        continue
    model_reactions.add(rxn.id)
    for metabolite in rxn.metabolites:
        metabolite_rxn_count_dict[metabolite.id] += 1
REACTIONS_TO_EXCLUDE = sorted(set(reactions_to_exclude))


# Also exclude highly connected metabolites
logger.info("Finding reactions to exclude")
METABOLITES_TO_EXCLUDE = list(
    set(
        map(
            lambda t: t[0],
            sorted(
                metabolite_rxn_count_dict.items(),
                key=lambda i: i[1],
                reverse=True,
            )[: CONFIG["mtb_tf"]["target_density"]["exclude-top-n-met"]],
        )
    )
)

# Create a list of nodes to exclude
nodes_to_exclude = REACTIONS_TO_EXCLUDE + METABOLITES_TO_EXCLUDE

# Create the actual metabolic network
logger.info("Create metabolic network")
metabolic_network = metworkpy.network.create_metabolic_network(
    model=BASE_MODEL,
    weighted=False,
    directed=False,
    nodes_to_remove=nodes_to_exclude,
)

# Project the metabolic network onto only the reactions
logger.info("Projecting metabolic network onto reactions only")
reaction_network = metworkpy.network.bipartite_project(
    network=metabolic_network,
    node_set=list(model_reactions),
    directed=False,
    weight=None,
)

# SECTION: TF targets
logger.info("Reading in TF target information")
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

# Create the TF target dataframe
logger.info("Finding TF targets")
tf_target_df = (
    tf_fc_df.abs() >= CONFIG["mtb_tf"]["target_density"]["target-fc-cutoff"]
) & (
    tf_pval_df <= CONFIG["mtb_tf"]["target_density"]["target-pval-cutoff"]
).astype("float")

# Read in the reaction information dataframe
rxn_info_df = pd.read_csv(
    CACHE_PATH / "model_information" / "reaction_information.csv"
)

# SECTION: Target Density
logger.info("Finding TF target density")
tf_target_density_list = []
for tf, target_series in tf_target_df.items():
    logger.info(f"Finding target density for {tf}")
    target_density_series = metworkpy.network.gene_target_density(
        metabolic_network=reaction_network,
        metabolic_model=BASE_MODEL,
        gene_labels=target_series,
        radius=CONFIG["mtb_tf"]["target_density"]["radius"],
    )
    target_density_series.name = tf
    tf_target_density_list.append(target_density_series)
# Combine all the target densities into a single dataframe
logger.info("Found all target density, combining results")
tf_density_df = pd.concat(tf_target_density_list, axis=1)

# Add in the reaction information
tf_density_df = tf_density_df.merge(
    rxn_info_df, how="left", left_index=True, right_on="id"
)

# SECTION: Save the final results
logger.info("Saving final results")
tf_density_df.to_csv(RESULTS_PATH / "tf_target_density.csv", index=True)
logger.info("Finished finding TF target density! ;)")

# SECTION: Target Enrichment
logger.info("Finding TF target enrichment")
tf_target_enrichment_list: list[pd.Series] = []
for tf, target_series in tf_target_df.items():
    logger.info(f"Finding target enrichment for {tf}")
    target_enrichment_series = metworkpy.network.gene_target_enrichment(
        metabolic_network=reaction_network,
        metabolic_model=BASE_MODEL,
        gene_targets=list(target_series[target_series].index),
        alternative="greater",
        metric="p-value",
        radius=CONFIG["mtb_tf"]["target_density"]["radius"],
    )
    target_enrichment_series = pd.Series(
        stats.false_discovery_control(target_enrichment_series),
        index=target_enrichment_series.index,
    )
    target_enrichment_series.name = tf
    tf_target_enrichment_list.append(target_enrichment_series)
# Combine all the target enrichments into a single dataframe
logger.info("Found all target densities, combining results")
tf_enrichment_df = pd.concat(tf_target_enrichment_list, axis=1)

# Add the reaction information to the tf_enrichment _df
tf_enrichment_df = tf_enrichment_df.merge(
    rxn_info_df, how="left", left_index=True, right_on="id"
)

# Save the final dataframe
tf_enrichment_df.to_csv(RESULTS_PATH / "tf_target_enrichment.csv", index=False)
