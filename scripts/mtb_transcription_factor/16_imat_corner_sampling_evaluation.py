"""
Script to investigate the iMAT corner sampling based approach using ArgR
as a test case
"""

# Setup
# Imports
# Standard Library Imports
import logging
import pathlib
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
MODEL_PATH = BASE_PATH / "models"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"
RESULTS_PATH = BASE_PATH / "results" / "mtb_transcription_factors"

# Create directories if needed
LOG_PATH.mkdir(parents=True, exist_ok=True)
RESULTS_PATH.mkdir(parents=True, exist_ok=True)

# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)


if __name__ == "__main__":
    # Setup Logging
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        filename=LOG_PATH / "16_imat_corner_sampling_evaluation.log",
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

    tfoe_micro_df["median"] = tfoe_micro_median

    # Generate the FVA model for ArgR, and median, then sample from
    # them and compute the divergence
    # First the gene weights for both of them
    argr_gene_weights = metworkpy.utils.expr_to_imat_gene_weights(
        tfoe_micro_df["Rv1657"],
        quantile=(
            CONFIG["mtb_tf"]["imat"]["lower-quantile"],
            CONFIG["mtb_tf"]["imat"]["upper-quantile"],
        ),
        subset=BASE_MODEL_GENES
        if CONFIG["mtb_tf"]["imat"]["subset-genes"]
        else None,
    )
    median_gene_weights = metworkpy.utils.expr_to_imat_gene_weights(
        tfoe_micro_df["median"],
        quantile=(
            CONFIG["mtb_tf"]["imat"]["lower-quantile"],
            CONFIG["mtb_tf"]["imat"]["upper-quantile"],
        ),
        subset=BASE_MODEL_GENES
        if CONFIG["mtb_tf"]["imat"]["subset-genes"]
        else None,
    )
    # Convert these to reaction weights
    argr_rxn_weights = metworkpy.gpr.gene_to_rxn_weights(
        BASE_MODEL, argr_gene_weights
    )
    median_rxn_weights = metworkpy.gpr.gene_to_rxn_weights(
        BASE_MODEL, median_gene_weights
    )

    # Generate iMAT models
    argr_fva_imat = metworkpy.imat.generate_model(
        model=BASE_MODEL,
        rxn_weights=argr_rxn_weights,
        method="fva",
        epsilon=CONFIG["mtb_tf"]["imat"]["epsilon"],
        threshold=CONFIG["mtb_tf"]["imat"]["threshold"],
        objective_tolerance=CONFIG["mtb_tf"]["imat"]["objective-tolerance"],
        processes=1,  # Solver parallelism works better, and more stably
    )
    median_fva_imat = metworkpy.imat.generate_model(
        model=BASE_MODEL,
        rxn_weights=median_rxn_weights,
        method="fva",
        epsilon=CONFIG["mtb_tf"]["imat"]["epsilon"],
        threshold=CONFIG["mtb_tf"]["imat"]["threshold"],
        objective_tolerance=CONFIG["mtb_tf"]["imat"]["objective-tolerance"],
        processes=1,  # Solver parallelism works better, and more stably
    )

    # Get samples from each of the models
    argr_fva_sampler = cobra.sampling.OptGPSampler(
        model=argr_fva_imat,
        thinning=CONFIG["mtb_tf"]["sampling"]["thinning"],
        processes=CONFIG["processes"],
    )
    median_fva_sampler = cobra.sampling.OptGPSampler(
        model=argr_fva_imat,
        thinning=CONFIG["mtb_tf"]["sampling"]["thinning"],
        processes=CONFIG["processes"],
    )
    argr_fva_samples = argr_fva_sampler.sample(
        CONFIG["mtb_tf"]["sampling"]["num-samples"]
    )
    argr_fva_samples = argr_fva_samples[
        argr_fva_sampler.validate(argr_fva_samples) == "v"  # type:ignore
    ]
    median_fva_samples = median_fva_sampler.sample(
        CONFIG["mtb_tf"]["sampling"]["num-samples"]
    )
    median_fva_samples = median_fva_samples[
        median_fva_sampler.validate(median_fva_samples) == "v"  # type:ignore
    ]

    # Compute the divergence between argr_fva and median
    argr_fva_divergence = (
        metworkpy.divergence.kl_divergence_functions.kl_divergence_array(
            median_fva_samples,
            argr_fva_samples,
            processes=CONFIG["processes"],
            n_neighbors=CONFIG["mtb_tf"]["divergence"]["n-neighbors"],
            distance_metric=CONFIG["mtb_tf"]["divergence"]["metric"],
            discrete=False,
            jitter=1e-12,
            jitter_seed=230472309487,
        )
    )
    assert isinstance(argr_fva_divergence, pd.Series)

    # Now, use the iMAT corner sampling based approach
    argr_corner_samples = metworkpy.imat.imat_sampling(
        model=BASE_MODEL,
        rxn_weights=argr_rxn_weights,
        epsilon=CONFIG["mtb_tf"]["imat"]["epsilon"],
        threshold=CONFIG["mtb_tf"]["imat"]["threshold"],
        n_samples=CONFIG["mtb_tf"]["sampling"]["num-samples"],
        processes=CONFIG["processes"],
        seed=237928734,
        fva_kwargs={"loopless": None},
    )
    median_corner_samples = metworkpy.imat.imat_sampling(
        model=BASE_MODEL,
        rxn_weights=argr_rxn_weights,
        epsilon=CONFIG["mtb_tf"]["imat"]["epsilon"],
        threshold=CONFIG["mtb_tf"]["imat"]["threshold"],
        n_samples=CONFIG["mtb_tf"]["sampling"]["num-samples"],
        processes=CONFIG["processes"],
        seed=528736283764,
        fva_kwargs={"loopless": None},
    )
    argr_corner_divergence = (
        metworkpy.divergence.kl_divergence_functions.kl_divergence_array(
            median_corner_samples,
            argr_corner_samples,
            processes=CONFIG["processes"],
            n_neighbors=CONFIG["mtb_tf"]["divergence"]["n-neighbors"],
            distance_metric=CONFIG["mtb_tf"]["divergence"]["metric"],
            discrete=False,
            jitter=1e-12,
            jitter_seed=230472309487,
        )
    )

    # Combine the two methods of getting the divergence into a dataframe and save it
    results_df = pd.DataFrame(
        {
            "FVA iMAT Divergence": argr_fva_divergence,
            "Corner Sampling Divergence": argr_corner_divergence,
        }
    )

    results_df.to_csv(
        RESULTS_PATH / "corner_sampling_divergence_comparison.csv",
        index=True,
    )
