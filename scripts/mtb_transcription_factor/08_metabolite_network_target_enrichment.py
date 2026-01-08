"""
Script to analyze enrichment of metabolite targets in
different metabolite networks
"""

# Setup
# Imports
# Standard Library Imports
import logging
import pathlib
import sys
import tomllib

# External Imports
import cobra  # type:ignore
import metworkpy  # type:ignore
from metabolic_modeling_utils.false_discovery_control import fdr_with_nan
import numpy as np
import pandas as pd
from scipy import stats  # type:ignore

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

# Setup logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=LOG_PATH / "08_metabolite_network_target_enrichment.log",
    filemode="w",
    level=logging.INFO,
)

# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)

# Read in the base model for finding metabolite networks
cobra.Configuration().solver = CONFIG["cobra"]["solver"]
BASE_MODEL = metworkpy.read_model(
    MODELS_PATH / "iEK1011_v2_7H9_ADC_glycerol.json"
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
metabolite_synthesis_rxn_network_df = get_metabolite_network(
    metabolite_synthesis_rxn_net_df_path,
    model=BASE_MODEL,
    rxns_to_remove=reactions_to_remove,
    proportion=CONFIG["mtb_tf"]["target_enrichment"]["essential-proportion"],
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
    proportion=CONFIG["mtb_tf"]["target_enrichment"]["reaction-proportion"],
    synthesis=False,
    processes=CONFIG["processes"],
)


# Determine the TF targets

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
    tf_fc_df.abs() >= CONFIG["mtb_tf"]["target_enrichment"]["target-fc-cutoff"]
) & (tf_pval_df <= CONFIG["mtb_tf"]["target_enrichment"]["target-pval-cutoff"])

# Get a set of genes in the model
model_gene_set = set(BASE_MODEL.genes.list_attr("id"))

# Now, for each type of metabolite network (synthesis, consumption)
# Identify any significant overlaps
target_enrichment_res_list: list[pd.DataFrame] = []
for metabolite_network_direction, metabolite_network_df in zip(
    ["synthesis", "consumption"],
    [
        metabolite_synthesis_rxn_network_df,
        metabolite_consumption_rxn_network_df,
    ],
):
    for tf, tf_target_series in tf_target_df.items():
        # Find the TF targets which are in the model
        tf_target_set = set(
            tf_target_series[tf_target_series].index
        ).intersection(model_gene_set)
        if len(tf_target_set) <= 3:
            continue
        tf_enrichment_df = pd.DataFrame(
            np.nan,
            columns=pd.Index(
                [
                    "metabolite network direction",
                    "metabolite network size",
                    "tf target count",
                    "tf target-metabolite network overlap",
                    "total genes",
                    "odds-ratio",
                    "p-value",
                    "adj p-value",
                ]
            ),
            index=metabolite_network_df.columns,
        )
        tf_enrichment_df["metabolite network direction"] = (
            metabolite_network_direction
        )
        tf_enrichment_df["tf target count"] = len(tf_target_set)
        tf_enrichment_df["total genes"] = len(model_gene_set)
        for metabolite, metabolite_net_series in metabolite_network_df.items():
            # Create the contingency table
            contingency_table = pd.DataFrame(
                np.nan,
                index=pd.Index(["in met net", "out met net"]),
                columns=pd.Index(["tf target", "not tf target"]),
            )
            # Find the metabolite network gene set
            met_net_gene_set = set(
                metworkpy.reaction_to_gene_list(
                    model=BASE_MODEL,
                    reaction_list=list(
                        metabolite_net_series[metabolite_net_series].index
                    ),
                    essential=False,
                )
            )
            if len(met_net_gene_set) <= 3:
                continue
            # Fill out the contingency table
            contingency_table.loc["in met net", "tf target"] = len(
                tf_target_set & met_net_gene_set
            )
            contingency_table.loc["in met net", "not tf target"] = len(
                met_net_gene_set - tf_target_set
            )
            contingency_table.loc["out met net", "tf target"] = len(
                tf_target_set - met_net_gene_set
            )
            contingency_table.loc["out met net", "not tf target"] = len(
                model_gene_set - (tf_target_set | met_net_gene_set)
            )
            # Perform the enrichment test
            fisher_res = stats.fisher_exact(
                table=contingency_table.to_numpy(),
                alternative="greater",
            )
            # Fill out the results table for this row
            tf_enrichment_df.loc[metabolite, "metabolite network size"] = len(  # type:ignore
                met_net_gene_set
            )
            tf_enrichment_df.loc[
                metabolite, "tf target-metabolite network overlap"  # type:ignore
            ] = contingency_table.loc["in met net", "tf target"]
            tf_enrichment_df.loc[metabolite, "odds-ratio"] = (
                fisher_res.statistic
            )  # type:ignore
            tf_enrichment_df.loc[metabolite, "p-value"] = fisher_res.pvalue  # type:ignore
        tf_enrichment_df = tf_enrichment_df.reset_index(
            drop=False, names="metabolite"
        )
        tf_enrichment_df["tf"] = tf
        tf_enrichment_df["adj p-value"] = fdr_with_nan(
            tf_enrichment_df["p-value"]
        )
        # Drop rows which still have NaN
        tf_enrichment_df = tf_enrichment_df.dropna(axis="index")
        target_enrichment_res_list.append(tf_enrichment_df)
# Combine all of the target enrichment
target_enrichment_res_df = pd.concat(target_enrichment_res_list, axis=0)

# Read in the metabolite information dataframe
metabolite_info_df = pd.read_csv(
    CACHE_PATH / "model_information" / "metabolite_information.csv"
)

# Join the metabolite informatino to the target enrichment results
target_enrichment_res_df = pd.merge(
    target_enrichment_res_df,
    metabolite_info_df,
    how="left",
    left_on="metabolite",
    right_on="id",
)

# Save the results
target_enrichment_res_df.to_csv(
    RESULTS_PATH / "tf_target_metabolite_network_enrichment.csv", index=False
)
