"""
Perform GSVA for the transcription factor expression data in the
metabolite networks (consuming and synthesis)
"""

# Setup
# Imports
# Standard Library Imports
import logging
import pathlib
from typing import cast
import sys
import tomllib

# External Imports
import cobra  # type:ignore
import decoupler as dc  # type: ignore
import metworkpy  # type:ignore
from metworkpy.metabolites import metabolite_network  # type:ignore
import pandas as pd

# Local Imports
from common_functions import get_metabolite_network


# Path setup
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".")  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
DATA_PATH = BASE_PATH / "data"
CACHE_PATH = BASE_PATH / "cache"
METABOLITE_NETWORKS_PATH = CACHE_PATH / "metabolite_networks" / "7H9_ADC"
RESULTS_PATH = BASE_PATH / "results" / "mtb_transcription_factors"
MODELS_PATH = BASE_PATH / "models"
LOG_PATH = BASE_PATH / "mtb_transcription_factors"

# Make required directories if they don't exist
CACHE_PATH.mkdir(parents=True, exist_ok=True)
METABOLITE_NETWORKS_PATH.mkdir(parents=True, exist_ok=True)
RESULTS_PATH.mkdir(parents=True, exist_ok=True)
LOG_PATH.mkdir(parents=True, exist_ok=True)

# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)

# Read in the base model for finding subsystems
cobra.Configuration().solver = CONFIG["cobra"]["solver"]
BASE_MODEL = metworkpy.read_model(
    MODELS_PATH / "iEK1011_v2_7H9_ADC_glycerol.json"
)

# Setup logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=LOG_PATH / "07_metabolite_gsva.log",
    filemode="w",
    level=logging.INFO,
)
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
metabolite_synthesis_rxn_network_df = get_metabolite_network(
    metabolite_synthesis_rxn_net_df_path,
    model=BASE_MODEL,
    rxns_to_remove=reactions_to_remove,
    proportion=CONFIG["mtb_tf"]["metabolite_gsva"]["essential-proportion"],
    synthesis=True,
    processes=CONFIG["processes"],
)


# Repeat for the consuming networks
metabolite_consumption_rxn_network_df_path = (
    METABOLITE_NETWORKS_PATH / "metabolite_consumption_reaction_network.csv"
)
metabolite_consumption_rxn_network_df = get_metabolite_network(
    metabolite_consumption_rxn_network_df_path,
    model=BASE_MODEL,
    rxns_to_remove=reactions_to_remove,
    proportion=CONFIG["mtb_tf"]["metabolite_gsva"]["reaction-proportion"],
    synthesis=False,
    processes=CONFIG["processes"],
)


# Create a rxn to gene dataframe to translate reaction networks
# into gene networks
rxn_list = []
gene_list = []
for rxn in BASE_MODEL.reactions:
    for gene in rxn.genes:
        rxn_list.append(rxn.id)
        gene_list.append(gene.id)
rxn_to_gene_df = pd.DataFrame({"reaction": rxn_list, "gene": gene_list})

# Pivot the network dfs longer, so one row per reaction metabolite pair
metabolite_synthesis_network_long = pd.merge(
    metabolite_synthesis_rxn_network_df.reset_index(
        drop=False, names="reaction"
    ).melt(id_vars="reaction", var_name="metabolite", value_name="to_keep"),
    rxn_to_gene_df,
    how="left",
    left_on="reaction",
    right_on="reaction",
)
metabolite_synthesis_network_long = cast(
    pd.DataFrame,
    metabolite_synthesis_network_long[
        metabolite_synthesis_network_long["to_keep"]
    ][["metabolite", "gene"]],
)
# Add a suffix to the networks to keep track of synthesis vs consumption
metabolite_synthesis_network_long["metabolite"] = (
    metabolite_synthesis_network_long["metabolite"] + "_synthesis_network"  # type: ignore
)
# Repeat above for the consumption network
metabolite_consumption_network_long = pd.merge(
    metabolite_consumption_rxn_network_df.reset_index(
        drop=False, names="reaction"
    ).melt(id_vars="reaction", var_name="metabolite", value_name="to_keep"),
    rxn_to_gene_df,
    how="left",
    left_on="reaction",
    right_on="reaction",
)
metabolite_consumption_network_long = cast(
    pd.DataFrame,
    metabolite_consumption_network_long[
        metabolite_consumption_network_long["to_keep"]
    ][["metabolite", "gene"]],
)
# Add a suffix to the networks to keep track of whether it is for synthesis or
# consumption
metabolite_consumption_network_long["metabolite"] = (
    # type: ignore
    metabolite_consumption_network_long["metabolite"] + "_consumption_network"
)
# Join the two metabolite network dataframes
metabolite_network_combined = pd.concat(
    [metabolite_synthesis_network_long, metabolite_consumption_network_long],
    axis=0,
    ignore_index=True,
)
metabolite_network_combined = metabolite_network_combined.drop_duplicates()
# Change the column names to the expected names for decoupler
metabolite_network_combined = metabolite_network_combined.rename(
    {"metabolite": "source", "gene": "target"}, axis=1
)

# Read in the log2(fold-change) data
tfoe_l2fc = (
    pd.read_excel(
        DATA_PATH / "mtb_transcription_factors" / "tfoe_targets.xlsx",
        sheet_name="SupplementaryTableS2",
        skiprows=list(range(8)) + [9],
        usecols="A,E:HB",
        index_col=0,
    )
    .T.rename(
        {
            "Rv0061": "Rv0061c",
            "Rv2427Ac": "Rv2427A",
        },
        axis=1,
    )
    .drop(
        [
            "Rv1784",
            "Rvns01",
            "Rvns02",
            "Rvnt01",
            "Rvnt02",
            "Rvnt03",
            "Rvnt05",
            "Rvnt06",
            "Rvnt07",
            "Rvnt08",
            "Rvnt11",
            "Rvnt12",
            "Rvnt13",
            "Rvnt15",
            "Rvnt17",
            "Rvnt19",
            "Rvnt21",
            "Rvnt22",
            "Rvnt24",
            "Rvnt27",
            "Rvnt28",
            "Rvnt29",
            "Rvnt30",
            "Rvnt32",
            "Rvnt33",
            "Rvnt34",
            "Rvnt40",
            "Rvnt41",
        ],
        axis=1,
    )
)

# Perform the GSVA analysis
gsva_res, _ = dc.mt.gsva(  # type:ignore
    data=tfoe_l2fc, net=metabolite_network_combined, kcdf="gaussian"
)
# gsva_res is a dataframe since that is the type of data
gsva_res = cast(pd.DataFrame, gsva_res)

# Save the GSVA results
gsva_res.to_csv(RESULTS_PATH / "metabolite_gsva.csv", index=True)
