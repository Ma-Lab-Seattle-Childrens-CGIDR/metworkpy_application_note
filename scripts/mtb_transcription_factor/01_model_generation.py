"""
Script to generate condition specific models for all of the transcription
factor over expression strains
"""

# Setup
# Imports
# Standard Library Imports
import logging
import pathlib
import warnings
import sys
import tomllib

# External Imports
import cobra
import metworkpy
import pandas as pd

# Local Imports

# Setup Path
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".").absolute()  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
DATA_PATH = BASE_PATH / "data"
CACHE_PATH = BASE_PATH / "cache"
MODEL_PATH = BASE_PATH / "models"
MODEL_OUT_PATH = CACHE_PATH / "tf_models"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"

# Create directories if needed
CACHE_PATH.mkdir(parents=True, exist_ok=True)
MODEL_OUT_PATH.mkdir(parents=True, exist_ok=True)
LOG_PATH.mkdir(parents=True, exist_ok=True)

# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)


def create_imat_model(
    base_model: cobra.Model, tf: str, expr_series: pd.Series, subset
):
    logger.info(f"*************Creating iMAT model for {tf}*************")
    model_out_path = MODEL_OUT_PATH / f"{tf}.json"
    if model_out_path.exists():
        return  # Model already generated
    logger.info("Finding gene weights")
    gene_weights = metworkpy.utils.expr_to_imat_gene_weights(
        expr_series,
        quantile=(
            CONFIG["mtb_tf"]["imat"]["lower-quantile"],
            CONFIG["mtb_tf"]["imat"]["upper-quantile"],
        ),
        subset=subset,
    )
    logger.info("Finding reaction weights")
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"Genes .* are in model but not.*",
            category=UserWarning,
            module=r".*metworkpy.*",
        )
        rxn_weights = metworkpy.gpr.gene_to_rxn_weights(
            model=base_model,
            gene_weights=gene_weights,
            fn_dict=metworkpy.gpr.gpr_functions.IMAT_FUNC_DICT,
            fill_val=0,
        )
    logger.info("Generating iMAT model")
    imat_model = metworkpy.imat.generate_model(
        model=BASE_MODEL,
        rxn_weights=rxn_weights,
        method=CONFIG["mtb_tf"]["imat"]["method"],
        loopless=None,
        epsilon=CONFIG["mtb_tf"]["imat"]["epsilon"],
        threshold=CONFIG["mtb_tf"]["imat"]["threshold"],
        objective_tolerance=CONFIG["mtb_tf"]["imat"]["objective-tolerance"],
        processes=1,  # Solver parallelism works better, and more stably
    )
    logger.info("Saving the model")
    metworkpy.write_model(model=imat_model, model_path=model_out_path)


if __name__ == "__main__":
    # Setup Logging
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        filename=LOG_PATH / "01_model_generation.log",
        filemode="w",
        level=logging.INFO,
    )

    # Change the default cobra solver
    cobra.Configuration().solver = CONFIG["cobra"]["solver"]

    # Read in base model
    logger.info("Reading in m7H9 model")
    BASE_MODEL = metworkpy.read_model(
        MODEL_PATH / "iEK1011_v2_7H9_ADC_glycerol.json"
    )
    BASE_MODEL_GENES = BASE_MODEL.genes.list_attr("id")

    # Read in the Transcription factor expression data, and the
    # metadata
    logger.info("Reading in the TFOE expression data")
    tfoe_micro_df = pd.read_csv(
        DATA_PATH / "mtb_transcription_factors" / "trip_microarray.csv",
        index_col=0,
    )
    tfoe_micro_df.index = tfoe_micro_df.index.str.replace("\n", "")
    tfoe_meta_df = pd.read_csv(
        DATA_PATH / "mtb_transcription_factors" / "trip_microarray_meta.tsv",
        sep="\t",
    )
    tfoe_meta_df["locus"] = (
        "Rv"
        + tfoe_meta_df["Title"].str.extract(
            r"TFOE_\d{4}_(\d{4}[ABc]?)(_[ABC])?"
        )[0]
    )

    # Get the median expression across
    tfoe_micro_median = tfoe_micro_df.median(axis=0)

    # Add on the Locus tag information
    tfoe_micro_df = (
        tfoe_micro_df.merge(
            tfoe_meta_df[["Accession", "locus"]],
            how="left",
            left_index=True,
            right_on="Accession",
        )
        .set_index("Accession")
        .groupby("locus")
        .median()
    ).T.rename(
        {
            "Rv0061": "Rv0061c",
            "Rv2427Ac": "Rv2427A",
        },
        axis=1,
    )

    tfoe_micro_df["wildtype"] = tfoe_micro_median

    # For each TF/Sample, create an iMAT model
    for tf, expr_series in tfoe_micro_df.items():
        assert isinstance(tf, str)
        create_imat_model(
            base_model=BASE_MODEL,
            tf=tf,
            expr_series=expr_series,
            subset=BASE_MODEL_GENES
            if CONFIG["mtb_tf"]["imat"]["subset-genes"]
            else None,
        )
