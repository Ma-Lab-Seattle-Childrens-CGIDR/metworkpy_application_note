"""
Script to generate condition specific models for all of the transcription factor over expression strains
"""

# Setup
# Imports
# Standard Library Imports
import logging
import pathlib
import warnings

# External Imports
import cobra  # type:ignore
import numpy as np
import metworkpy  # type:ignore
import pandas as pd  # type:ignore

# Local Imports

# Setup Path
try:
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
except NameError:
    BASE_PATH = pathlib.Path(".").absolute()
DATA_PATH = BASE_PATH / "data"
CACHE_PATH = BASE_PATH / "cache"
MODEL_PATH = BASE_PATH / "models"
MODEL_OUT_PATH = CACHE_PATH / "tf_models"

# Create directories if needed
CACHE_PATH.mkdir(parents=True, exist_ok=True)
MODEL_OUT_PATH.mkdir(parents=True, exist_ok=True)

# Setup Logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=BASE_PATH
    / "logs"
    / "transcription_factors"
    / "01_model_generation.log",
    filemode="w",
    level=logging.INFO,
)

# Change the default cobra solver
cobra.Configuration().solver = "hybrid"

# Read in base model
logger.info("Reading in m7H9 model")
BASE_MODEL = metworkpy.read_model(
    MODEL_PATH / "iEK1011_v2_7H9_ADC_glycerol.json"
)

# Run Parameters
PROCESSES = 12
EPSILON = 1.0
THRESHOLD = 0.01
NEG_FOLD_CHANGE = -1
POS_FOLD_CHANGE = 1
OBJECTIVE_TOLERANCE = 5e-2

# Read in the transcription factor differential expression data
logger.info("Reading in the TFOE differential expression data")

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

# Next Generate the differential expression models
for transcription_factor in tfoe_l2fc.index:
    logger.info(f"Starting transcription factor: {transcription_factor}")
    model_out_path = MODEL_OUT_PATH / f"{transcription_factor}.json"
    if model_out_path.exists():
        continue  # Model already generated
    # Convert log2fc to gene weights
    logger.info("Converting log2fc to gene weights")
    sample_l2fc = tfoe_l2fc.loc[transcription_factor]
    gene_weights = np.zeros(sample_l2fc.size)
    gene_weights[sample_l2fc <= NEG_FOLD_CHANGE] = -1
    gene_weights[sample_l2fc >= POS_FOLD_CHANGE] = 1
    gene_weights = pd.Series(gene_weights, index=sample_l2fc.index)  # type:ignore
    # Convert gene weights to reaction weights
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"Genes .* are in model but not.*",
            category=UserWarning,
            module=r".*metworkpy.*",
        )
        rxn_weights = metworkpy.gpr.gene_to_rxn_weights(
            model=BASE_MODEL,
            gene_weights=gene_weights,
            fn_dict=metworkpy.gpr.gpr_functions.IMAT_FUNC_DICT,
            fill_val=0,
        )
    # Generate iMAT Model
    logger.info("Generating diff model")
    imat_model = metworkpy.imat.generate_model(
        model=BASE_MODEL,
        rxn_weights=rxn_weights,
        method="fva",
        loopless=False,
        epsilon=EPSILON,
        threshold=THRESHOLD,
        objective_tolerance=OBJECTIVE_TOLERANCE,
    )
    # Save the model
    logger.info("Saving model")
    metworkpy.write_model(model=imat_model, model_path=model_out_path)
# For the diff models, there is not WT version other than the base model,
# so that will be used as the WT
metworkpy.write_model(
    model=BASE_MODEL, model_path=MODEL_OUT_PATH / "wildtype.json"
)
