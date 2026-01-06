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

# External Imports
import cobra  # type: ignore
import metworkpy  # type:ignore
import pandas as pd

# Local Imports
from metabolic_modeling_utils import gene_target_density

# Path Setup
try:
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
except NameError:
    BASE_PATH = pathlib.Path(".").absolute()
DATA_PATH = BASE_PATH / "data"
CACHE_PATH = BASE_PATH / "cache"
RESULTS_PATH = BASE_PATH / "results" / "transcription_factors"
BASE_MODEL_PATH = BASE_PATH / "models" / "iEK1011_v2_7H9_ADC_glycerol.json"

# Create Directories if needed
CACHE_PATH.mkdir(parents=True, exist_ok=True)
RESULTS_PATH.mkdir(parents=True, exist_ok=True)

# Logging Config
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=BASE_PATH
    / "logs"
    / "transcription_factors"
    / "05_tf_target_density.log",
    filemode="w",
    level=logging.INFO,
)

# Script Parameters
EXCLUDE_TOP_N_MET = 10  # Number of highly connected metabolites to not include (h20, atp, etc.)
DENSITY_RADIUS = 1  # Radius to use when calculating density
TF_TARGET_FOLD_CHANGE_CUTOFF = (
    1.0  # log2(fold-change) cutoff for considering a gene a TF target
)
TF_TARGET_PVALUE_CUTOFF = (
    0.05  # p-value cutoff for considering a gene a TF target
)

# Read in the metabolic model
logger.info("Reading in the base model")
cobra.Configuration().solver = "hybrid"
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
            )[:EXCLUDE_TOP_N_MET],
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
    DATA_PATH / "transcription_factors" / "tfoe_targets.xlsx",
    sheet_name="SupplementaryTableS2",
    skiprows=list(range(8)) + [9],
    usecols="A,E:HB",
    index_col=0,
)
tf_pval_df = pd.read_excel(
    DATA_PATH / "transcription_factors" / "tfoe_targets.xlsx",
    sheet_name="SupplementaryTableS2",
    skiprows=list(range(8)) + [9],
    usecols="A,HC:OZ",
    index_col=0,
)
tf_pval_df.columns = tf_pval_df.columns.str.replace(".1", "")

# Create the TF target dataframe
logger.info("Finding TF targets")
tf_target_df = (tf_fc_df.abs() >= TF_TARGET_FOLD_CHANGE_CUTOFF) & (
    tf_pval_df <= TF_TARGET_PVALUE_CUTOFF
).astype("float")

# SECTION: Target Density
logger.info("Finding TF activating target density")
tf_target_density_list = []
for tf, target_series in tf_target_df.items():
    logger.info(f"Finding target density for {tf}")
    target_density_series = gene_target_density.gene_target_density(
        metabolic_network=reaction_network,
        metabolic_model=BASE_MODEL,
        gene_labels=target_series,
        radius=DENSITY_RADIUS,
    )
    target_density_series.name = tf
    tf_target_density_list.append(target_density_series)
# Combine all the target densities into a single dataframe
logger.info("Found all target density, combining results")
tf_density_df = pd.concat(tf_target_density_list, axis=1)

# SECTION: Save the final results
logger.info("Saving final results")
tf_density_df.to_csv(RESULTS_PATH / "tf_target_density.csv", index=True)
logger.info("Finished finding TF target density! ;)")
