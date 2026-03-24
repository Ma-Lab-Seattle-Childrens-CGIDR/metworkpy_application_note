"""
Script to examine the centrality in the reaction metabolic network
of transcription factor targets
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
import cobra
import metworkpy
import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

# Local Imports
from common_functions import fdr_with_nan, bootstrap

# Path Setup
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".")  # Use current dir as base path
    PROGRESS_BAR = True
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
    PROGRESS_BAR = False  # Don't use progress bar when run as script
DATA_PATH = BASE_PATH / "data"
RESULTS_PATH = BASE_PATH / "results" / "mtb_transcription_factors"
BASE_MODEL_PATH = BASE_PATH / "models" / "iEK1011_v2_7H9_ADC_glycerol.json"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"

# Create Directories if needed
RESULTS_PATH.mkdir(parents=True, exist_ok=True)
LOG_PATH.mkdir(parents=True, exist_ok=True)

# Logging Setup
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=LOG_PATH / "06_reaction_network_centrality.log",
    filemode="w",
    level=logging.INFO,
)

# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)

cobra.Configuration().solver = CONFIG["cobra"]["solver"]

# Read in the base model
logger.info("Reading in base model")
BASE_MODEL = metworkpy.read_model(BASE_MODEL_PATH)

# Get a list of genes in the model
model_genes = BASE_MODEL.genes.list_attr("id")

# Find reactions and metabolites to be removed
logger.info("Finding reactions/metabolites to remove")
subsystems_to_ignore = {
    "Intracellular demand",
    "Biomass and maintenance functions",
    "Extracellular exchange",
}
reactions_to_remove = set()
for rxn in BASE_MODEL.reactions:
    if rxn.subsystem in subsystems_to_ignore:
        reactions_to_remove.add(rxn.id)

# Find how many reactions each metabolite is involved in
metabolite_counts: dict[str, int] = defaultdict(int)
for rxn in BASE_MODEL.reactions:
    for metabolite in rxn.metabolites:
        metabolite_counts[metabolite.id] += 1
metabolites_to_exclude = set(
    map(
        lambda t: t[0],
        sorted(metabolite_counts.items(), key=lambda i: i[1], reverse=True)[
            : CONFIG["mtb_tf"]["network_properties"]["exclude-top-n-met"]
        ],
    )
)

# Create a list of nodes to exclude from the network
nodes_to_exclude = list(reactions_to_remove | metabolites_to_exclude)

# Construct the metabolic network
logger.info("Constructing the metabolic network")
metabolic_network = metworkpy.network.create_metabolic_network(
    model=BASE_MODEL,
    weighted=False,
    directed=CONFIG["mtb_tf"]["network_properties"]["directed"],
    nodes_to_remove=nodes_to_exclude,
)

# Project the network onto only the reactions
logger.info("Projecting the bipartite network onto the reactions only")
metabolic_reaction_network_nodes = [
    rxn
    for rxn in BASE_MODEL.reactions.list_attr("id")
    if rxn not in reactions_to_remove
]
metabolic_reaction_network = metworkpy.bipartite_project(
    metabolic_network,
    node_set=metabolic_reaction_network_nodes,
    weight=None,
)

# Find the closeness centrality of reactions in the network
logger.info("Calculating closeness centrality")
closeness_centrality_series = pd.Series(
    nx.closeness_centrality(metabolic_reaction_network)
)
# Find the betweeness centrality of reaction in the network
logger.info("Calculating betweenness centrality")
betweenness_centrality_series = pd.Series(
    nx.betweenness_centrality(metabolic_reaction_network, normalized=True)
)
# Save these results
pd.DataFrame(
    {
        "Closeness": closeness_centrality_series,
        "Betweenness": betweenness_centrality_series,
    }
).to_csv(RESULTS_PATH / "metabolic_reaction_network_centrality.csv")


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
    tf_fc_df.abs()
    >= CONFIG["mtb_tf"]["network_properties"]["target-fc-cutoff"]
) & (
    tf_pval_df <= CONFIG["mtb_tf"]["network_properties"]["target-pval-cutoff"]
)

# Create a dictionary of TF: reaction targets
logger.info("Creating a reaction target dictionary for all the TFs")
tf_target_dict: dict[str, list[str]] = {}
for tf, target_series in tf_target_df.items():
    tf_target_dict[str(tf)] = metworkpy.utils.gene_to_reaction_list(
        model=BASE_MODEL,
        gene_list=[
            str(g)
            for g in target_series[target_series].index
            if g in model_genes
        ],
    )

# Find which reactions in the model are associated with genes
model_reactions_assoc_with_genes = [
    r.id for r in BASE_MODEL.reactions if len(r.genes) > 0
]

# For each TF, test if it targets more central reactions
tf_centrality_results_df = pd.DataFrame(
    np.nan,
    index=tf_target_df.columns,
    columns=pd.Index(
        [
            "mean closeness",
            "std closeness",
            "mean betweenness",
            "std betweenness",
            "closeness u1",
            "closeness u2",
            "closeness rho",
            "closeness p-value",
            "closeness adj p-value",
            "closeness bootstrap p-value",
            "closeness bootstrap adj p-value",
            "betweenness u1",
            "betweenness u2",
            "betweenness rho",
            "betweenness p-value",
            "betweenness adj p-value",
            "betweenness bootstrap p-value",
            "betweenness bootstrap adj p-value",
        ]
    ),
)

logger.info("Testing centrality for all of the TFs")
for tf, tf_target_rxns in tqdm(
    tf_target_dict.items(), disable=not PROGRESS_BAR
):
    logger.info(f"Starting on {tf}")
    non_targeted_rxns = [  # Really inefficient, but not rate limiting
        r for r in model_reactions_assoc_with_genes if r not in tf_target_rxns
    ]
    if len(tf_target_rxns) < 5 or len(non_targeted_rxns) < 5:
        continue  # Exclude TFs which either target too few/many reactions
    # Start by finding the mean/std for closeness and betweenness
    tf_centrality_results_df.loc[tf, "mean closeness"] = (
        closeness_centrality_series[tf_target_rxns].mean()
    )
    tf_centrality_results_df.loc[tf, "std closeness"] = (
        closeness_centrality_series[tf_target_rxns].std()
    )
    tf_centrality_results_df.loc[tf, "mean betweenness"] = (
        betweenness_centrality_series[tf_target_rxns].mean()
    )
    tf_centrality_results_df.loc[tf, "std betweenness"] = (
        betweenness_centrality_series[tf_target_rxns].std()
    )
    # Perform the Mann-Whitney U-test for the closeness centrality
    targeted_rxn_closeness = closeness_centrality_series[tf_target_rxns]
    non_targeted_rxn_closeness = closeness_centrality_series[non_targeted_rxns]
    mannu_res = stats.mannwhitneyu(
        x=targeted_rxn_closeness,
        y=non_targeted_rxn_closeness,
        alternative="greater",
    )
    closeness_u1 = mannu_res.statistic
    closeness_u2 = (
        len(targeted_rxn_closeness) * len(non_targeted_rxn_closeness)
        - closeness_u1
    )
    closeness_rho = closeness_u1 / (closeness_u2 + closeness_u1)
    closeness_pval = mannu_res.pvalue

    tf_centrality_results_df.loc[tf, "closeness u1"] = closeness_u1
    tf_centrality_results_df.loc[tf, "closeness u2"] = closeness_u2
    tf_centrality_results_df.loc[tf, "closeness rho"] = closeness_rho
    tf_centrality_results_df.loc[tf, "closeness p-value"] = closeness_pval

    # Perform the Mann-Whitney U-test for the betweenness centrality
    targeted_rxn_betweenness = betweenness_centrality_series[tf_target_rxns]
    non_targeted_rxn_betweenness = betweenness_centrality_series[
        non_targeted_rxns
    ]
    mannu_res = stats.mannwhitneyu(
        x=targeted_rxn_betweenness,
        y=non_targeted_rxn_betweenness,
        alternative="greater",
    )
    betweenness_u1 = mannu_res.statistic
    betweenness_u2 = (
        len(targeted_rxn_betweenness) * len(non_targeted_rxn_betweenness)
        - betweenness_u1
    )
    betweenness_rho = betweenness_u1 / (betweenness_u2 + betweenness_u1)
    betweenness_pval = mannu_res.pvalue

    tf_centrality_results_df.loc[tf, "betweenness u1"] = betweenness_u1
    tf_centrality_results_df.loc[tf, "betweenness u2"] = betweenness_u2
    tf_centrality_results_df.loc[tf, "betweenness rho"] = betweenness_rho
    tf_centrality_results_df.loc[tf, "betweenness p-value"] = betweenness_pval

    # Perform Bootstrap p-value tests for mean betweeness and closeness
    _, closeness_centrality_bootstrap_pvalue = bootstrap(
        closeness_centrality_series[model_reactions_assoc_with_genes],
        sample_set=tf_target_rxns,
        statistic=np.mean,
        iterations=CONFIG["mtb_tf"]["network_properties"][
            "bootstrap-iterations"
        ],
        alternative=CONFIG["mtb_tf"]["network_properties"][
            "bootstrap-alternative"
        ],
        seed=324,
    )
    tf_centrality_results_df.loc[tf, "closeness bootstrap p-value"] = (
        closeness_centrality_bootstrap_pvalue
    )
    _, betweenness_centrality_bootstrap_pvalue = bootstrap(
        betweenness_centrality_series[model_reactions_assoc_with_genes],
        sample_set=tf_target_rxns,
        statistic=np.mean,
        iterations=CONFIG["mtb_tf"]["network_properties"][
            "bootstrap-iterations"
        ],
        alternative=CONFIG["mtb_tf"]["network_properties"][
            "bootstrap-alternative"
        ],
        seed=324,
    )
    tf_centrality_results_df.loc[tf, "betweenness bootstrap p-value"] = (
        betweenness_centrality_bootstrap_pvalue
    )

# The bootstrap can generate p-values slightly greater than 1.0
# due to floating point error
tf_centrality_results_df["closeness bootstrap p-value"] = (
    tf_centrality_results_df["closeness bootstrap p-value"].clip(
        lower=0.0, upper=1.0
    )
)
tf_centrality_results_df["betweenness bootstrap p-value"] = (
    tf_centrality_results_df["betweenness bootstrap p-value"].clip(
        lower=0.0, upper=1.0
    )
)

# Perform the false discovery correction
tf_centrality_results_df["closeness adj p-value"] = fdr_with_nan(
    tf_centrality_results_df["closeness p-value"]
)
tf_centrality_results_df["closeness bootstrap adj p-value"] = fdr_with_nan(
    tf_centrality_results_df["closeness bootstrap p-value"]
)
tf_centrality_results_df["betweenness adj p-value"] = fdr_with_nan(
    tf_centrality_results_df["betweenness p-value"]
)
tf_centrality_results_df["betweenness bootstrap adj p-value"] = fdr_with_nan(
    tf_centrality_results_df["betweenness bootstrap p-value"]
)

# Read in the gene information from mycobrowser
gene_info_df = pd.read_csv(
    DATA_PATH / "gene_info" / "Mycobacterium_tuberculosis_H37Rv_txt_v5.txt",
    sep="\t",
)
# Join this with the centrality results
tf_centrality_results_df = pd.merge(
    tf_centrality_results_df,
    gene_info_df,
    how="left",
    left_index=True,
    right_on="Locus",
).set_index("Locus")

# Save the results
logger.info("Finished analyzing centrality, saving results")
tf_centrality_results_df.to_csv(
    RESULTS_PATH / "metabolic_reaction_network_centrality_analysis.csv",
    index=True,
)
logger.info("Finished performing centrality analysis! ;)")
