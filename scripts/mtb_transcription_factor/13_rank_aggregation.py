"""
Script for performing rank aggregation for the TF Divergence
and the TF Target density
"""

# Setup
# Imports
# Standard Library Imports
import logging
import pathlib
import sys
import tomllib

# External Imports
import numpy as np
import pandas as pd
from robustrankaggregpy.aggregate_ranks import (
    rank_matrix_from_df,
    aggregate_ranks,
)

# Path Setup
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".").absolute()  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
DATA_PATH = BASE_PATH / "data"
CACHE_PATH = BASE_PATH / "cache"
RESULTS_PATH = BASE_PATH / "results" / "mtb_transcription_factors"
BASE_MODEL_PATH = BASE_PATH / "models" / "iEK1011_v2_7H9_ADC_glycerol.json"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"


if __name__ == "__main__":
    # Create Directories if needed
    CACHE_PATH.mkdir(parents=True, exist_ok=True)
    RESULTS_PATH.mkdir(parents=True, exist_ok=True)
    LOG_PATH.mkdir(parents=True, exist_ok=True)
    # Logging Config
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        filename=LOG_PATH / "14_rank_aggregation.log",
        filemode="w",
        level=logging.INFO,
    )
    # Read in the configuration file
    with open(BASE_PATH / "config.toml", "rb") as f:
        CONFIG = tomllib.load(f)

    # Read in the Density results, and set the reaction id as the index
    density_df = pd.read_csv(
        RESULTS_PATH / "tf_target_density.csv", index_col=0
    ).set_index("id")

    # Read in the divergence results
    divergence_df = pd.read_csv(
        RESULTS_PATH / "divergence_results.csv", index_col=0
    )
    divergence_df = divergence_df.loc[
        :, divergence_df.columns.str.startswith("reaction__")
    ]
    divergence_df.columns = divergence_df.columns.str.replace(
        "^reaction__", "", regex=True
    )
    # Replace infs with nan
    divergence_df = divergence_df.replace([np.inf, -np.inf], np.nan)
    # Drop any columns which are all nan
    divergence_df = divergence_df.dropna(axis=1, how="all")
    # Transpose to match the density df
    divergence_df = divergence_df.T

    # Find the set of reactions in both
    density_rxn_set = set(density_df.index)
    divergence_rxn_set = set(divergence_df.index)
    common_reactions = sorted(density_rxn_set & divergence_rxn_set)
    # Find the set of transcription factors in both
    density_tf_set = set(density_df.columns)
    divergence_tf_set = set(divergence_df.columns)
    common_tfs = sorted(density_tf_set & divergence_tf_set)

    # Filter for only the common reactions/tfs
    density_df = density_df.loc[common_reactions, common_tfs]
    divergence_df = divergence_df.loc[common_reactions, common_tfs]

    # Create the results dataframe
    results_df = pd.DataFrame(
        np.nan, index=pd.Index(common_reactions), columns=pd.Index(common_tfs)
    )

    # Iterate through each TF, and get the density and perform RRA
    for tf in results_df.columns:
        # Construct the dataframe for performing the ranking
        density_series = density_df[tf]
        divergence_series = divergence_df[tf]
        comb_df = pd.DataFrame(
            {"density": density_series, "divergence": divergence_series}
        )
        # Rank this matrix
        rank_matrix = rank_matrix_from_df(
            data=comb_df,
            full=True,  # True means that NaN are preserved
            ascending=False,
            rank_method="max",
        )
        # Perform the rank aggregation
        agg_series = aggregate_ranks(
            rank_matrix=rank_matrix,
            method="rra",
            full=True,
        )
        # Add this to the results series
        results_df.loc[agg_series.index, tf] = agg_series

    # Add in the reaction information
    rxn_info_df = pd.read_csv(
        CACHE_PATH / "model_information" / "reaction_information.csv"
    )
    results_df = results_df.merge(
        rxn_info_df, how="left", left_index=True, right_on="id"
    )

    # Save the results
    results_df.to_csv(RESULTS_PATH / "tf_rra.csv", index=False)
