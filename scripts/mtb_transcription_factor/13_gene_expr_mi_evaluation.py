"""
Script for evaluating the Mutual information network
constructed from the RNAseq gene expression compendia from
Yoo et al.
"""

# Setup
# Imports
# Standard Library Imports
import logging
import pathlib
import sys
import tomllib

# External Imports
import metworkpy
import networkx as nx
import pandas as pd
from scipy import stats

# Local Imports


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
MI_NETWORK_ADJACENCY_PATH = CACHE_PATH / "mtb_mutual_information"
RESULTS_PATH = BASE_PATH / "results" / "mtb_transcription_factors"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"

if __name__ == "__main__":
    # Setup Logging
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        filename=LOG_PATH / "13_gene_expr_mi_evaluation.log",
        filemode="w",
        level=logging.INFO,
    )
    # Read in the configuration file
    with open(BASE_PATH / "config.toml", "rb") as f:
        CONFIG = tomllib.load(f)

    # Read in the gene expression data
    rnaseq_compendia = pd.read_csv(
        DATA_PATH / "rnaseq_compendia" / "log_tpm_norm.csv", index_col=0
    ).T

    # Find the mutual information network (if required)
    mi_adj_out_path = MI_NETWORK_ADJACENCY_PATH / "gene_expr_mi_adjacency.csv"
    if mi_adj_out_path.exists():
        mi_adj_df = pd.read_csv(mi_adj_out_path, index_col=0)
    else:
        mi_adj_df = metworkpy.information.mi_pairwise(
            rnaseq_compendia,
            calculate_pvalue=False,
            processes=CONFIG["processes"],
            n_neighbors=CONFIG["mtb_tf"]["mutual-information"]["n-neighbors"],
            metric_x=CONFIG["mtb_tf"]["mutual-information"]["x-metric"],
            metric_y=CONFIG["mtb_tf"]["mutual-information"]["y-metric"],
            clip=True,
        ).fillna(0.0)  # type: ignore  # Fill any NA values with 0
        mi_adj_df.to_csv(mi_adj_out_path, index=True)

    # Convert from the mutual information adjacency matrix to a networkx network
    mi_network: nx.Graph = nx.from_pandas_adjacency(
        mi_adj_df, create_using=nx.Graph
    )

    # Evaluate the relationship between the mutual information centrality and
    # essentiality
    results_dict: dict[str, float] = {}
    # Essentiality/Vulnerability index
    vi_df = pd.read_excel(
        DATA_PATH / "gene_info" / "bosch_vi.xlsx",
        sheet_name="(1) Mtb H37Rv",
        index_col=0,
        usecols="A,B,D,E,Z",
    )
    vi_df.index = vi_df.index.str.replace("^RVBD", "Rv", regex=True)
    compendia_genes = set(mi_adj_df.index)
    vi_genes = set(vi_df.index)
    common_genes = sorted(compendia_genes & vi_genes)
    # Get the tnseq essentiality series
    tnseq_ess = vi_df.loc[common_genes, "tnseq_ess"]
    # Get the VI series
    vi_series = vi_df.loc[common_genes, "Vulnerability Index"]

    # Compute the eigenvector centrality of the networkx graph
    gene_centrality_series = pd.Series(
        nx.eigenvector_centrality(mi_network, weight="weight")
    )[common_genes]
    gene_centrality_series.name = "Gene Expr MI Eigenvector Centrality"
    gene_centrality_series.to_csv(
        RESULTS_PATH / "gene_expr_mi_gene_centrality.csv"
    )

    # Start by comparing the mutual information centrality of essential genes
    # and non-essential genes
    mannu_res = stats.mannwhitneyu(
        gene_centrality_series[tnseq_ess == "Essential"],
        gene_centrality_series[tnseq_ess != "Essential"],
        alternative="greater",  # The distribution underlying x is stochastically greater than y
    )
    u1 = mannu_res.statistic
    u2 = (
        gene_centrality_series[tnseq_ess == "Essential"].shape[0]
        * gene_centrality_series[tnseq_ess != "Essential"].shape[0]
        - u1
    )
    auc = u1 / (u1 + u2)
    pvalue = mannu_res.pvalue
    results_dict["TNSeq essential vs not Mann-Whitney U-test u1"] = u1
    results_dict["TNSeq essential vs not Mann-Whitney U-test u2"] = u2
    results_dict["TNSeq essential vs not Mann-Whitney U-test auc"] = auc
    results_dict["TNSeq essential vs not Mann-Whitney U-test p-value"] = pvalue

    # Evaluate the correlation between the Mutual Information Centrality and the Vulnerability Index
    # Using Kendall-Tau
    # NOTE: More negative VI values are MORE vulnerable
    kendall_res = stats.kendalltau(
        gene_centrality_series, vi_series, alternative="less"
    )
    spearman_res = stats.spearmanrho(
        gene_centrality_series, vi_series, alternative="less"
    )
    pearson_res = stats.pearsonr(
        gene_centrality_series, vi_series, alternative="less"
    )
    results_dict[
        "Vulnerability Index Correlation MI Centrality Kendall-Tau statistic"
    ] = kendall_res.statistic
    results_dict[
        "Vulnerability Index Correlation MI Centrality Kendall-Tau p-value"
    ] = kendall_res.pvalue
    results_dict[
        "Vulnerability Index Correlation MI Centrality Spearman's Rho statistic"
    ] = spearman_res.statistic
    results_dict[
        "Vulnerability Index Correlation MI Centrality Spearman's Rho p-value"
    ] = spearman_res.pvalue
    results_dict[
        "Vulnerability Index Correlation MI Centrality Pearson's R statistic"
    ] = pearson_res.statistic
    results_dict[
        "Vulnerability Index Correlation MI Centrality Pearson's R p-value"
    ] = pearson_res.pvalue

    # Create a series to write this to a csv
    essentiality_res_series = pd.Series(
        results_dict, name="Mutual Information Essentiality Statistical Tests"
    )
    # Save the results
    essentiality_res_series.to_csv(
        RESULTS_PATH
        / "gene_expr_mutual_information_vi_essentiality_statistics.csv",
        index=True,
    )
