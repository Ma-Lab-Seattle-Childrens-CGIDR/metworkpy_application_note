"""
Script to compare divergence approach to
iMAT plus FBA
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
import pandas as pd  # type:ignore

# Local Imports

# Setup Path
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".").absolute()  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
CACHE_PATH = BASE_PATH / "cache"
DATA_PATH = BASE_PATH / "data"
MODEL_PATH = BASE_PATH / "models"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"
RESULTS_PATH = BASE_PATH / "results" / "mtb_transcription_factors"

# Create directories if needed
LOG_PATH.mkdir(parents=True, exist_ok=True)
RESULTS_PATH.mkdir(parents=True, exist_ok=True)

# Setup Logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=LOG_PATH / "11_imat_test.log",
    filemode="w",
    level=logging.INFO,
)

# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)

# Change the default cobra solver
cobra.Configuration().solver = CONFIG["cobra"]["solver"]

# Read in base model
logger.info("Reading in m7H9 model")
BASE_MODEL = metworkpy.read_model(
    MODEL_PATH / "iEK1011_v2_7H9_ADC_glycerol.json"
)
# Read in the transcription factor differential expression data
logger.info("Reading in the TFOE differential expression data")

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
).T

# Get the gene expression for ArgR (Rv1657)
gene_expr = tfoe_l2fc["Rv1657"]

# Create a gene weights series
gene_weights = pd.Series(0.0, index=gene_expr.index)
gene_weights[gene_expr >= CONFIG["mtb_tf"]["imat"]["pos-fold-change"]] = 1.0
gene_weights[gene_expr <= CONFIG["mtb_tf"]["imat"]["neg-fold-change"]] = -1.0

# Convert the gene weights into reactions weights
rxn_weights = metworkpy.gene_to_rxn_weights(
    model=BASE_MODEL,
    gene_weights=gene_weights,
)

# Find the fluxes for the base model (pFBA)
base_fluxes = cobra.flux_analysis.pfba(
    model=BASE_MODEL, fraction_of_optimum=0.95
).fluxes
if not isinstance(base_fluxes, pd.Series):
    raise ValueError("Couldn't get fluxes from base model")
base_fluxes.name = "pFBA fluxes"

# Find the fluxes for the IMAT directly
imat_fluxes = metworkpy.imat.imat(
    model=BASE_MODEL,
    rxn_weights=rxn_weights,
    epsilon=CONFIG["mtb_tf"]["imat"]["epsilon"],
    threshold=CONFIG["mtb_tf"]["imat"]["threshold"],
).fluxes
if not isinstance(imat_fluxes, pd.Series):
    raise ValueError("Couldn't get fluxes from imat")
imat_fluxes.name = "IMAT fluxes"

# Create an FVA model, then perform pFBA
imat_model = metworkpy.imat.generate_model(
    model=BASE_MODEL,
    rxn_weights=rxn_weights,
    method="fva",
    epsilon=CONFIG["mtb_tf"]["imat"]["epsilon"],
    threshold=CONFIG["mtb_tf"]["imat"]["threshold"],
    loopless=False,
    processes=CONFIG["processes"],
    objective_tolerance=CONFIG["mtb_tf"]["imat"]["objective-tolerance"],
)
fva_model_fluxes = cobra.flux_analysis.pfba(
    imat_model, fraction_of_optimum=0.95
).fluxes
if not isinstance(fva_model_fluxes, pd.Series):
    raise ValueError("Couldn't get fluxes from fva imat model")
fva_model_fluxes.name = "FVA IMAT pFBA fluxes"

# Combine the fluxes together
results_df = pd.concat([base_fluxes, imat_fluxes, fva_model_fluxes], axis=1)

# Calculate the differences in fluxes
results_df["diff imat"] = results_df["IMAT fluxes"] - results_df["pFBA fluxes"]
results_df["diff fva imat"] = (
    results_df["FVA IMAT pFBA fluxes"] - results_df["pFBA fluxes"]
)

# Read in the reaction information dataframe
rxn_info_df = pd.read_csv(
    CACHE_PATH / "model_information" / "reaction_information.csv"
)

# Join the reaction information dataframe to the results dataframe
results_df = pd.merge(
    results_df, rxn_info_df, how="left", left_index=True, right_on="id"
)

# Save the results dataframe
results_df.to_csv(RESULTS_PATH / "imat_compare.csv")
