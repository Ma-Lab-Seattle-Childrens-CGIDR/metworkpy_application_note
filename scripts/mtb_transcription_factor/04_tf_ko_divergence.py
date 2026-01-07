"""
Script to evaluate the impact TF gene targets
on the metabolic network using ko-divergence
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
import numpy as np  # type:ignore
import pandas as pd  # type:ignore
from scipy import stats  # type:ignore

# Local Imports
from metabolic_modeling_utils.false_discovery_control import fdr_with_nan

# Path Setup
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".")  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
DATA_PATH = BASE_PATH / "data"
RESULTS_PATH = BASE_PATH / "results" / "mtb_transcription_factors"
CACHE_PATH = BASE_PATH / "cache"
GENE_KO_DIVERGENCE_PATH = CACHE_PATH / "gene_ko_divergence" / "7h9_adc"
BASE_MODEL_PATH = BASE_PATH / "models" / "iEK1011_v2_7H9_ADC_glycerol.json"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"

# Create Directories if needed
CACHE_PATH.mkdir(parents=True, exist_ok=True)
GENE_KO_DIVERGENCE_PATH.mkdir(parents=True, exist_ok=True)
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
SUBSYSTEMS_TO_IGNORE = {
    "Biomass and maintenance functions",
    "Extracellular exchange",
    "Intracellular demand",
}
SUBSYSTEM_RENAME_DICT = {
    "Propanoate metabolism": "Propanoate Metabolism",
    "Arabinogalactan bioynthesis": "Arabinogalactan biosynthesis",
}

# Read in the Base Model
logger.info("Reading in the base model")
BASE_MODEL = metworkpy.read_model(BASE_MODEL_PATH)

# Create a dictionary of subsystem to reactions
logger.info("Creating subsystem dictionary")
subsystem_to_rxn_dict: dict[str, list[str]] = defaultdict(list)
for rxn in BASE_MODEL.reactions:
    if rxn.subsystem in SUBSYSTEMS_TO_IGNORE:
        continue
    subsys = (
        SUBSYSTEM_RENAME_DICT[rxn.subsystem]
        if rxn.subsystem in SUBSYSTEM_RENAME_DICT
        else rxn.subsystem
    )
    subsystem_to_rxn_dict[subsys].append(rxn.id)
    subsystem_to_rxn_dict["whole_metabolism"].append(rxn.id)

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
    ko_divergence_df = metworkpy.divergence.ko_divergence(
        model=BASE_MODEL,
        genes_to_ko=BASE_MODEL.genes.list_attr("id"),
        target_networks=subsystem_to_rxn_dict,
        divergence_metric="kl",
        n_neighbors=CONFIG["mtb"]["ko_divergence"]["n-neighbors"],
        sample_count=CONFIG["mtb"]["ko_divergence"]["sample-count"],
        processes=CONFIG["processes"],
    ).clip(lower=0.0)
    logger.info("Saving the ko divergence results")
    ko_divergence_df.to_csv(ko_divergence_df_path, index=True)
# Filter out columns that are all nan
ko_divergence_df = ko_divergence_df.dropna(axis="columns")

# Get a list of genes in the model
model_gene_set = set(BASE_MODEL.genes.list_attr("id"))

# Find the TF targets
logger.info("Finding the TF targets")
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

tf_target_df = (
    tf_fc_df.abs() >= CONFIG["mtb"]["ko_divergence"]["target-fc-cutoff"]
) & (tf_pval_df <= CONFIG["mtb"]["ko_divergence"]["target-pval-cutoff"])

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
    if len(v) >= CONFIG["mtb"]["ko_divergence"]["min-target-count"]
}

# Now, for each TF test if it targets genes which cause
# relatively large perturbations in each of the various subsystems
logger.info("Performing tests for TF target ko-divergence")
results_df_list: list[pd.DataFrame] = []

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
    for subsystem, ko_divergence_series in ko_divergence_df.items():
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
            < CONFIG["mtb"]["ko_divergence"]["min-target-count"]
            or len(non_targeted_divergence)
            < CONFIG["mtb"]["ko_divergence"]["min-target-count"]
        ):
            continue  # Don't calculate for TFs which target too many, or too few genes
        # Calculate the results
        # Calculate the Mann-Whitney test results
        mann_whitney_res = stats.mannwhitneyu(
            x=targeted_divergence, y=non_targeted_divergence
        )
        u1 = mann_whitney_res.statistic
        u2 = len(targeted_divergence) * len(non_targeted_divergence) - u1
        rho = u1 / (len(targeted_divergence) * len(non_targeted_divergence))
        pval = mann_whitney_res.pvalue
        tf_res_df.loc[subsystem, "Mann-Whitney U1"] = u1  # type:ignore
        tf_res_df.loc[subsystem, "Mann-Whitney U2"] = u2  # type:ignore
        tf_res_df.loc[subsystem, "rho"] = rho  # type:ignore
        tf_res_df.loc[subsystem, "p-value"] = pval  # type:ignore
    logger.info(f"Finished performing tests for {tf}")
    # Correct for false discovery rate
    tf_res_df["adj p-value"] = fdr_with_nan(tf_res_df["p-value"])
    # Reset the index to get a subsystem column
    tf_res_df = tf_res_df.reset_index(drop=False, names="subsystem")
    tf_res_df["tf"] = tf
    # Add the results dataframe to the overall results
    results_df_list.append(tf_res_df)

# Combine the results for all the TFs
logger.info("Finished performing results, combining across TFs")
tf_ko_divergence_res_df = pd.concat(results_df_list, axis=0)

# Save the final results
logger.info("Saving the final results")
tf_ko_divergence_res_df.to_csv(
    RESULTS_PATH / "ko_divergence_tf_target_analysis.csv",
    index=False,
)
logger.info("Finished performing KO-divergence analysis for the TFs! ;)")
