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
import cobra
import metworkpy
import numpy as np
import pandas as pd
from scipy import stats

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


if __name__ == "__main__":
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
    BASE_MODEL_GENES = BASE_MODEL.genes.list_attr("id")
    # Read in the transcription factor differential expression data
    logger.info("Reading in the TFOE microarray expression data")

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

    # Get the gene expression for ArgR (Rv1657)
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

    # Convert the gene weights into reactions weights
    argr_rxn_weights = metworkpy.gene_to_rxn_weights(
        model=BASE_MODEL,
        gene_weights=argr_gene_weights,
    )
    median_rxn_weights = metworkpy.gene_to_rxn_weights(
        model=BASE_MODEL,
        gene_weights=argr_gene_weights,
    )

    # Find the fluxes for the base model (pFBA)
    base_fluxes = cobra.flux_analysis.pfba(
        model=BASE_MODEL, fraction_of_optimum=0.95
    ).fluxes
    if not isinstance(base_fluxes, pd.Series):
        raise ValueError("Couldn't get fluxes from base model")
    base_fluxes.name = "Unconstrained pFBA fluxes"

    # Find the fluxes for the IMAT directly
    argr_imat_fluxes = metworkpy.imat.imat(
        model=BASE_MODEL,
        rxn_weights=argr_rxn_weights,
        epsilon=CONFIG["mtb_tf"]["imat"]["epsilon"],
        threshold=CONFIG["mtb_tf"]["imat"]["threshold"],
    ).fluxes
    if not isinstance(argr_imat_fluxes, pd.Series):
        raise ValueError("Couldn't get fluxes from imat")
    argr_imat_fluxes.name = "ArgR iMAT fluxes"
    median_imat_fluxes = metworkpy.imat.imat(
        model=BASE_MODEL,
        rxn_weights=median_rxn_weights,
        epsilon=CONFIG["mtb_tf"]["imat"]["epsilon"],
        threshold=CONFIG["mtb_tf"]["imat"]["threshold"],
    ).fluxes
    if not isinstance(median_imat_fluxes, pd.Series):
        raise ValueError("Couldn't get fluxes from imat")
    median_imat_fluxes.name = "Median iMAT fluxes"

    # Create an FVA model, then perform pFBA
    argr_imat_model = metworkpy.imat.generate_model(
        model=BASE_MODEL,
        rxn_weights=argr_rxn_weights,
        method="fva",
        epsilon=CONFIG["mtb_tf"]["imat"]["epsilon"],
        threshold=CONFIG["mtb_tf"]["imat"]["threshold"],
        loopless=False,
        processes=CONFIG["processes"],
        objective_tolerance=CONFIG["mtb_tf"]["imat"]["objective-tolerance"],
    )
    median_imat_model = metworkpy.imat.generate_model(
        model=BASE_MODEL,
        rxn_weights=median_rxn_weights,
        method="fva",
        epsilon=CONFIG["mtb_tf"]["imat"]["epsilon"],
        threshold=CONFIG["mtb_tf"]["imat"]["threshold"],
        loopless=False,
        processes=CONFIG["processes"],
        objective_tolerance=CONFIG["mtb_tf"]["imat"]["objective-tolerance"],
    )
    argr_fva_model_fluxes = cobra.flux_analysis.pfba(
        argr_imat_model, fraction_of_optimum=0.95
    ).fluxes
    if not isinstance(argr_fva_model_fluxes, pd.Series):
        raise ValueError("Couldn't get fluxes from fva imat model")
    argr_fva_model_fluxes.name = "ArgR FVA IMAT pFBA fluxes"
    median_fva_model_fluxes = cobra.flux_analysis.pfba(
        median_imat_model, fraction_of_optimum=0.95
    ).fluxes
    if not isinstance(median_fva_model_fluxes, pd.Series):
        raise ValueError("Couldn't get fluxes from fva imat model")
    median_fva_model_fluxes.name = "Median FVA IMAT pFBA fluxes"

    # Combine the fluxes together
    results_df = pd.concat(
        [
            base_fluxes,
            argr_imat_fluxes,
            median_imat_fluxes,
            argr_fva_model_fluxes,
            median_fva_model_fluxes,
        ],
        axis=1,
    )

    # Calculate the differences in fluxes
    results_df["diff imat"] = (
        results_df["ArgR iMAT fluxes"] - results_df["Median iMAT fluxes"]
    )
    results_df["diff fva imat"] = (
        results_df["ArgR FVA IMAT pFBA fluxes"]
        - results_df["Median FVA IMAT pFBA fluxes"]
    )

    # Get samples from the Median Model and the iMAT FVA model
    median_imat_samples = cobra.sampling.sample(
        model=median_imat_model,
        n=CONFIG["mtb_tf"]["sampling"]["num-samples"],
        method="optgp",
        thinning=CONFIG["mtb_tf"]["sampling"]["thinning"],
        processes=CONFIG["processes"],
        seed=812309712,
    )
    argr_imat_samples = cobra.sampling.sample(
        model=argr_imat_model,
        n=CONFIG["mtb_tf"]["sampling"]["num-samples"],
        method="optgp",
        thinning=CONFIG["mtb_tf"]["sampling"]["thinning"],
        processes=CONFIG["processes"],
        seed=923875928735,
    )
    # For each reaction, calculate t-tests and ks-tests
    stat_test_res = pd.DataFrame(
        np.nan,
        index=median_imat_samples.columns,
        columns=pd.Index(
            [
                "t-stat",
                "t p-value",
                "ks-stat",
                "ks p-value",
                "kruskal stat",
                "kruskal p-value",
                "mannwhitneyu stat",
                "mannwhitneyu p-value",
            ]
        ),
    )

    # Perform statistical tests for the samples
    for rxn in stat_test_res.index:
        ttest = stats.ttest_ind(
            median_imat_samples[rxn],
            argr_imat_samples[rxn],
            alternative="two-sided",
            equal_var=False,
        )
        kstest = stats.ks_2samp(
            median_imat_samples[rxn],
            argr_imat_samples[rxn],
            alternative="two-sided",
        )
        kruskal = stats.kruskal(
            median_imat_samples[rxn],
            argr_imat_samples[rxn],
        )
        mannu = stats.mannwhitneyu(
            median_imat_samples[rxn],
            argr_imat_samples[rxn],
            alternative="two-sided",
        )
        stat_test_res.loc[rxn, "t-stat"] = ttest.statistic
        stat_test_res.loc[rxn, "t p-value"] = ttest.pvalue
        stat_test_res.loc[rxn, "ks-stat"] = kstest.statistic
        stat_test_res.loc[rxn, "ks p-value"] = kstest.pvalue
        stat_test_res.loc[rxn, "kruskal stat"] = kruskal.statistic
        stat_test_res.loc[rxn, "kruskal p-value"] = kruskal.pvalue
        stat_test_res.loc[rxn, "mannwhitneyu stat"] = mannu.statistic
        stat_test_res.loc[rxn, "mannwhitneyu p-value"] = mannu.pvalue

    # Add on the statistical test comparison
    results_df = results_df.merge(
        stat_test_res, how="left", left_index=True, right_index=True
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
