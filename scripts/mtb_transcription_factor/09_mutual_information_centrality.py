"""
Script to analyze the mutual information centrality of the
iEK1011_v2 model, and the tendency of the transcription factors to
target more or less central reactions.
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
from metabolic_modeling_utils.false_discovery_control import fdr_with_nan
import networkx as nx
import numpy as np
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
FLUX_SAMPLES_PATH = CACHE_PATH / "tf_model_flux_samples"
RESULTS_PATH = BASE_PATH / "results" / "mtb_transcription_factors"
MODELS_PATH = BASE_PATH / "models"
LOG_PATH = BASE_PATH / "mtb_transcription_factors"

# Create directories if needed
CACHE_PATH.mkdir(parents=True, exist_ok=True)
MI_NETWORK_ADJACENCY_PATH.mkdir(parents=True, exist_ok=True)
LOG_PATH.mkdir(parents=True, exist_ok=True)

# Setup Logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=LOG_PATH / "08_metabolite_network_target_enrichment.log",
    filemode="w",
    level=logging.INFO,
)
# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)

# Read in the base model for finding mutual information centrality
cobra.Configuration().solver = CONFIG["cobra"]["solver"]
BASE_MODEL = metworkpy.read_model(
    MODELS_PATH / "iEK1011_v2_7H9_ADC_glycerol.json"
)

# Find reactions which should be ignored
reactions_to_remove: list[str] = []
subsystems_to_ignore: list[str] = [
    "Intracellular demand",
    "Biomass and maintenance functions",
    "Extracellular exchange",
]
for rxn in BASE_MODEL.reactions:
    if rxn.subsystem in subsystems_to_ignore:
        reactions_to_remove.append(rxn.id)

# If needed get flux samples for the iEK1011_v2 model
logger.info("Generating samples if needed")
wt_flux_sample_path = FLUX_SAMPLES_PATH / "wildtype.parquet"
if wt_flux_sample_path.exists():
    logger.info("Reading in previous sampling")
    flux_sample_df = pd.read_parquet(wt_flux_sample_path)
else:
    logger.info("Starting sampling")
    sampler = cobra.sampling.OptGPSampler(
        model=BASE_MODEL,
        thinning=CONFIG["mtb_tf"]["sampling"]["thinning"],
        processes=CONFIG["processes"],
    )
    samples = flux_sample_df = sampler.sample(
        CONFIG["mtb_tf"]["sampling"]["num-samples"]
    )
    valid_samples = samples[sampler.validate(samples) == "v"]  # type: ignore
    logger.info(
        f"Validated samples, {
            len(valid_samples)
            / CONFIG['mtb_tf']['sampling']['num-samples']:.2%} valid"
    )
    valid_samples.to_parquet(wt_flux_sample_path, index=False)
    flux_sample_df = valid_samples

# Find the mutual information network (if required)
mi_adj_out_path = MI_NETWORK_ADJACENCY_PATH / "metabolic_mi_adjacency.csv"
if mi_adj_out_path.exists():
    mi_adj_df = pd.read_csv(mi_adj_out_path, index_col=0)
else:
    mi_adj_df = (
        metworkpy.information.mi_network_adjacency_matrix(
            samples=flux_sample_df.drop(
                reactions_to_remove, axis=1
            ),  # Drop reactions to remove
            processes=CONFIG["processes"],
            progress_bar=False,
        )
        .fillna(0.0)
        .clip(0.0)
    )  # Fill any NA values with 0, remove all negative values
    mi_adj_df.to_csv(mi_adj_out_path, index=True)


# Convert from the mutual information adjacency matrix to a networkx network
mi_network: nx.Graph = nx.from_pandas_adjacency(
    mi_adj_df, create_using=nx.Graph
)


# Compute the eigenvector centrality of the networkx graph
eigenvector_centrality = pd.Series(
    nx.eigenvector_centrality(mi_network, weight="weight")
)

# Read in the targets for all of the transcription factors
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


# For each of the transcription factors, evaluate if it has a
# significantly higher MI centrality than background
results_df = pd.DataFrame(
    np.nan,
    index=tf_target_df.columns,
    columns=pd.Index(["u1", "u2", "auc", "p-value"]),
)
logger.info("Evaluating Eigenvector centrality significance for TF targets")
for tf, tf_target_series in tf_target_df.items():
    target_gene_set = (
        set(tf_target_series[tf_target_series].index) & model_gene_set
    )
    # Translate the target gene set into a reaction list
    target_reaction_list = metworkpy.utils.gene_to_reaction_list(
        model=BASE_MODEL, gene_list=list(target_gene_set)
    )
    # Get the eigenvector centrality of these reactions
    targeted_centrality = eigenvector_centrality[target_reaction_list]
    non_targeted_centrality = eigenvector_centrality.drop(target_reaction_list)
    if (
        len(targeted_centrality)
        < CONFIG["mtb_tf"]["mutual_information"]["min-target-count"]
    ) or (
        len(non_targeted_centrality)
        < CONFIG["mtb_tf"]["mutual_information"]["min-target-count"]
    ):
        continue
    # Perform the Mann-Whitney U-test
    mannu_res = stats.mannwhitneyu(
        targeted_centrality,
        non_targeted_centrality,
        alternative="greater",
        nan_policy="omit",
    )
    # Get u1 and u2
    u1 = mannu_res.statistic
    u2 = len(targeted_centrality) * len(non_targeted_centrality) - u1
    # Get the AUC
    auc = u1 / (len(targeted_centrality) * len(non_targeted_centrality))
    # Get the p-value
    pvalue = mannu_res.pvalue
    # Save results to the dataframe
    results_df.loc[tf, "u1"] = u1
    results_df.loc[tf, "u2"] = u2
    results_df.loc[tf, "auc"] = auc
    results_df.loc[tf, "p-value"] = pvalue
# Perform false discovery correction
logger.info("Finished finding significance, adjusting p-values")
results_df["adj p-value"] = fdr_with_nan(results_df["p-value"])


# Save the final results
logger.info("Saving final results")
results_df.to_csv(RESULTS_PATH / "mutual_information_tf_target_centrality.csv")


# Evaluate the relationship between the mutual information centrality and
results_dict: dict[str, float] = {}
# Essentiality/Vulnerability index
vi_df = pd.read_excel(
    DATA_PATH / "gene_info" / "bosch_vi.xlsx",
    sheet_name="(1) Mtb H37Rv",
    index_col=0,
    usecols="A,B,D,E,Z",
)
vi_df.index = vi_df.index.str.replace("^RVBD", "Rv", regex=True)
# Find the overlap of genes in the model and in the essentiality data
model_genes = set(BASE_MODEL.genes.list_attr("id"))
vi_genes = set(vi_df.index)
common_genes = sorted(model_genes & vi_genes)
# Get the tnseq essentiality series
tnseq_ess = vi_df.loc[common_genes, "tnseq_ess"]
# Get the VI series
vi_series = vi_df.loc[common_genes, "Vulnerability Index"]

# Convert from reaction centrality to gene centrality
gene_centrality_series = pd.Series(np.nan, index=pd.Index(model_genes))
for gene in gene_centrality_series.index:
    centrality_list = []
    for rxn in BASE_MODEL.genes.get_by_id(gene).reactions:
        if rxn.id not in eigenvector_centrality.index:
            continue
        centrality_list.append(eigenvector_centrality[rxn.id])
    gene_centrality_series[gene] = max(centrality_list)
gene_centrality_series = gene_centrality_series[common_genes]
# Save the gene centrality series
gene_centrality_series.name = "Eigenvector Centrality"
gene_centrality_series.to_csv(RESULTS_PATH / "flux_mi_gene_centrality.csv")

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
    RESULTS_PATH / "mutual_information_vi_essentiality_statistics.csv",
    index=True,
)
