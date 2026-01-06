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

# External Imports
import cobra  # type:ignore
import decoupler as dc
import metworkpy  # type:ignore
from metworkpy.metabolites import metabolite_network  # type:ignore
import pandas as pd

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
MODELS_PATH = BASE_PATH / "models"

# Make required directories if they don't exist
CACHE_PATH.mkdir(parents=True, exist_ok=True)
METABOLITE_NETWORKS_PATH.mkdir(parents=True, exist_ok=True)
RESULTS_PATH.mkdir(parents=True, exist_ok=True)

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
    metabolite_synthesis_network_long["metabolite"] + "_synthesis_network"
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
# Add a suffix to the networks to keep track of whether it is for synthesis or consumption
metabolite_consumption_network_long["metabolite"] = (
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
        DATA_PATH / "transcription_factors" / "tfoe_targets.xlsx",
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
