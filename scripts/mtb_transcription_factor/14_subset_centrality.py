"""
Script for evaluating the centrality of reactions/metabolites
using subset centrality measures
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
import networkx as nx
import pandas as pd

# Local Imports

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
CACHE_PATH = BASE_PATH / "cache"
RESULTS_PATH = BASE_PATH / "results" / "mtb_transcription_factors"
BASE_MODEL_PATH = BASE_PATH / "models" / "iEK1011_v2_7H9_ADC_glycerol.json"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"

if __name__ == "__main__":
    # Create Directories if needed
    RESULTS_PATH.mkdir(parents=True, exist_ok=True)
    LOG_PATH.mkdir(parents=True, exist_ok=True)

    # Logging Setup
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        filename=LOG_PATH / "14_subset_centrality.log",
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
    metabolites_to_remove = set(
        metworkpy.network.get_top_metabolites(
            BASE_MODEL,
            n=CONFIG["mtb_tf"]["network_properties"]["exclude-top-n-met"],
            type="substrate",
        )
    )
    nodes_to_exclude = list(reactions_to_remove | metabolites_to_remove)

    # Construct the metabolic network
    logger.info("Constructing the metabolic network")
    metabolic_network = metworkpy.network.create_metabolic_network(
        model=BASE_MODEL,
        weighted=False,
        directed=CONFIG["mtb_tf"]["network_properties"]["directed"],
        nodes_to_remove=nodes_to_exclude,
    )
    # Determine the nodes in the network which are reactions
    # (not just all reactions since some reactions are removed)
    metabolic_network_rxn_nodes = {
        n
        for n in metabolic_network.nodes
        if n in BASE_MODEL.reactions.list_attr("id")
    }

    # Evaluate the closeness of nodes in the network
    base_closeness = pd.Series(nx.closeness_centrality(metabolic_network))
    # Evaluate the base betweenness of nodes in the network
    base_betweenness = pd.Series(
        metworkpy.network.centrality.betweenness_centrality_bipartite_subset(
            metabolic_network,
            metabolic_network_rxn_nodes,
            metabolic_network_rxn_nodes,
        )
    )

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
        tf_pval_df
        <= CONFIG["mtb_tf"]["network_properties"]["target-pval-cutoff"]
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

    subset_betweenness_df = pd.DataFrame(
        index=pd.Index(metabolic_network.nodes)
    )
    subset_closeness_df = pd.DataFrame(index=pd.Index(metabolic_network.nodes))
    subset_betweenness_df["FULL"] = base_betweenness
    subset_closeness_df["FULL"] = base_closeness

    # For each TF, calculate the subset betweenness/closeness, subtract off the
    # base centrality, and add it to the dataframe
    for tf, reaction_targets in tf_target_dict.items():
        subset_closeness_df[tf] = (
            pd.Series(
                metworkpy.network.centrality.closeness_centrality_subset(
                    metabolic_network, targets=reaction_targets
                )
            )
            - base_closeness
        )
        subset_betweenness_df[tf] = (
            pd.Series(
                metworkpy.network.centrality.betweenness_centrality_bipartite_subset(
                    metabolic_network,
                    node_partition=metabolic_network_rxn_nodes,
                    targets=reaction_targets,
                )
            )
            - base_betweenness
        )
    subset_closeness_df = subset_closeness_df.reset_index(
        drop=False, names="Node ID"
    )
    subset_betweenness_df = subset_betweenness_df.reset_index(
        drop=False, names="Node ID"
    )

    # Add in reaction and metabolite information
    rxn_info_df = pd.read_csv(
        CACHE_PATH / "model_information" / "reaction_information.csv"
    )
    metabolite_info_df = pd.read_csv(
        CACHE_PATH / "model_information" / "metabolite_information.csv"
    )

    # Join both of these
    subset_closeness_df = (
        subset_closeness_df.merge(
            rxn_info_df, how="left", left_on="Node ID", right_on="id"
        )
        .rename({"id": "reaction id"}, axis=1)
        .merge(
            metabolite_info_df, how="left", left_on="Node ID", right_on="id"
        )
        .rename({"id": "merabolite id"}, axis=1)
    )
    subset_betweenness_df = (
        subset_betweenness_df.merge(
            rxn_info_df, how="left", left_on="Node ID", right_on="id"
        )
        .rename({"id": "reaction id"}, axis=1)
        .merge(
            metabolite_info_df, how="left", left_on="Node ID", right_on="id"
        )
        .rename({"id": "merabolite id"}, axis=1)
    )

    # Save the results
    subset_closeness_df.to_csv(
        RESULTS_PATH / "tf_subset_closeness.csv", index=False
    )
    subset_betweenness_df.to_csv(
        RESULTS_PATH / "tf_subset_betweenness.csv", index=False
    )
