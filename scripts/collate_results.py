"""
Script to combine the results from the various analyses into
a single file
"""

# Setup
# Imports
# Standard Library Imports
import json
import pathlib
import sys
import tomllib
from typing import cast

# External Imports
import cobra
import pandas as pd
import metworkpy
import networkx as nx
from sklearn.preprocessing import StandardScaler

# Path Setup
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = (
        pathlib.Path(".").absolute().parent
    )  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent
MODEL_PATH = BASE_PATH / "models"
RESULTS_PATH = BASE_PATH / "results"
SIMULATION_RESULTS_PATH = RESULTS_PATH / "simulation"
MTB_TF_RESULTS_PATH = RESULTS_PATH / "mtb_transcription_factors"

# Script Configuration
# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)
cobra.Configuration().solver = CONFIG["cobra"]["solver"]


##########################
### Simulation Results ###
##########################
# Use a function to name space the simulation result collation
def collate_simulation_results():
    # -----------------------
    # -- Model Information --
    # -----------------------
    simulation_model = metworkpy.read_model(
        MODEL_PATH / "simulation_model.json"
    )
    model_df = pd.DataFrame(
        "",
        index=pd.Index(simulation_model.reactions.list_attr("id")),
        columns=pd.Index(["Equation", "Gene-Reaction Rule", "Subsystem"]),
    )
    for rxn in simulation_model.reactions:
        rxn = cast(cobra.Reaction, rxn)
        model_df.loc[rxn.id, "Equation"] = rxn.build_reaction_string()
        model_df.loc[rxn.id, "Gene-Reaction Rule"] = rxn.gene_reaction_rule
        model_df.loc[rxn.id, "Subsystem"] = rxn.subsystem
    model_df = model_df.reset_index(drop=False, names="Reaction")

    # ----------------------
    # -- Metabolic Graphs --
    # ----------------------
    def to_pandas_adjacency(graph: nx.Graph) -> pd.DataFrame:
        """Convert a networkx graph to a pandas dataframe"""
        idx = pd.Index(graph.nodes)
        return pd.DataFrame(nx.to_numpy_array(graph), index=idx, columns=idx)

    def read_graph(path: pathlib.Path) -> nx.Graph:
        """
        Read a networkx graph from a json file
        """
        with open(path, "r") as f:
            graph = nx.node_link_graph(
                json.load(f), multigraph=False, directed=False, edges="edges"
            )
        return graph

    # Read in the metabolic, reaction, and metabolite graphs
    metabolic_connectivity_adjacency = to_pandas_adjacency(
        read_graph(
            SIMULATION_RESULTS_PATH
            / "metabolic_networks"
            / "metabolic_network.json"
        )
    )
    reaction_connectivity_adjacency = to_pandas_adjacency(
        read_graph(
            SIMULATION_RESULTS_PATH
            / "metabolic_networks"
            / "metabolic_reaction_network.json"
        )
    )
    metabolite_connectivity_adjacency = to_pandas_adjacency(
        read_graph(
            SIMULATION_RESULTS_PATH
            / "metabolic_networks"
            / "metabolic_metabolite_network.json"
        )
    )

    metabolite_centrality = pd.read_csv(
        SIMULATION_RESULTS_PATH
        / "metabolic_networks"
        / "metabolite_centrality.csv",
        index_col=0,
    ).reset_index(drop=False, names="Metabolite")

    reaction_centrality = pd.read_csv(
        SIMULATION_RESULTS_PATH
        / "metabolic_networks"
        / "reaction_centrality.csv",
        index_col=0,
    ).reset_index(drop=False, names="Reaction")

    # --------------------
    # -- Flux MI Graphs --
    # --------------------
    mi_adjacency = pd.read_csv(
        SIMULATION_RESULTS_PATH / "mutual_information" / "mi_adjacency.csv",
        index_col=0,
    )
    mi_centrality = pd.read_csv(
        SIMULATION_RESULTS_PATH / "mutual_information" / "mi_centrality.csv",
        index_col=0,
    ).reset_index(drop=False, names="Reaction")
    # ----------------------------
    # -- Metabolite Subnetworks --
    # ----------------------------
    metabolite_consuming_network = pd.read_csv(
        SIMULATION_RESULTS_PATH
        / "metabolite_networks"
        / "metabolite_consuming_network.csv",
        index_col=0,
    )
    metabolite_synthesis_network = pd.read_csv(
        SIMULATION_RESULTS_PATH
        / "metabolite_networks"
        / "metabolite_synthesis_network.csv",
        index_col=0,
    )
    # ----------------------------------------------
    # -- Reaction Neighborhood Enrichment/Density --
    # ----------------------------------------------
    gene_target_density = pd.read_csv(
        SIMULATION_RESULTS_PATH / "target_density" / "gene_target_density.csv",
        index_col=0,
    ).reset_index(names="Reaction")
    gene_target_enrichment = pd.read_csv(
        SIMULATION_RESULTS_PATH
        / "target_density"
        / "gene_target_enrichment_pval.csv",
        index_col=0,
    ).reset_index(names="Reaction")

    # -------------------
    # -- KO Divergence --
    # -------------------
    ko_divergence = pd.read_csv(
        SIMULATION_RESULTS_PATH
        / "ko_divergence"
        / "ko_divergence_results.csv",
        index_col=0,
    ).reset_index(names="Gene")
    ko_divergence = ko_divergence.melt(
        id_vars="Gene", var_name="Divergence Group", value_name="Divergence"
    )
    ko_divergence["Divergence Group Type"] = (
        ko_divergence["Divergence Group"].str.split("__").str[0]
    )
    ko_divergence["Divergence Group"] = (
        ko_divergence["Divergence Group"]
        .str.split("__")
        .str[1:]
        .str.join("__")
    )
    ko_divergence = ko_divergence[
        ["Gene", "Divergence Group Type", "Divergence Group", "Divergence"]
    ]

    # ---------------------
    # -- IMAT Divergence --
    # ---------------------
    imat_div = pd.read_csv(
        SIMULATION_RESULTS_PATH / "iMAT" / "imat_divergence.csv", index_col=0
    ).reset_index(names="Divergence Group")
    imat_div["Divergence Group Type"] = (
        imat_div["Divergence Group"].str.split("__").str[0]
    )
    imat_div["Divergence Group"] = (
        imat_div["Divergence Group"].str.split("__").str[1:].str.join("__")
    )
    imat_div = imat_div[
        ["Divergence Group Type", "Divergence Group", "IMAT Divergence"]
    ]

    # ---------------------------
    # -- Save Combined Results --
    # ---------------------------
    with pd.ExcelWriter(
        RESULTS_PATH / "simulation_model_results.xlsx"
    ) as writer:
        # Model Information
        model_df.to_excel(writer, sheet_name="Model Information", index=False)
        # Metabolic Connectivity Networks
        metabolic_connectivity_adjacency.to_excel(
            writer, sheet_name="Metabolic SCN"
        )
        reaction_connectivity_adjacency.to_excel(
            writer, sheet_name="Reaction SCN", index=True
        )
        metabolite_connectivity_adjacency.to_excel(
            writer, sheet_name="Metabolite SCN", index=True
        )
        metabolite_centrality.to_excel(
            writer, sheet_name="SCN Metabolite Centrality", index=False
        )
        reaction_centrality.to_excel(
            writer, sheet_name="SCN Reaction Centrality", index=False
        )
        # Mutual information networks
        mi_adjacency.to_excel(
            writer, sheet_name="Flux MI Adjacency Matrix", index=True
        )
        mi_centrality.to_excel(
            writer, sheet_name="Flux MI Centrality", index=False
        )
        # Metabolite Networks
        metabolite_synthesis_network.to_excel(
            writer, sheet_name="Metabolite synthesis Networks", index=True
        )
        metabolite_consuming_network.to_excel(
            writer, sheet_name="Metabolite consuming Networks", index=True
        )
        # Reaction neighborhood enrichment/density
        gene_target_density.to_excel(
            writer, sheet_name="Gene Target Density", index=False
        )
        gene_target_enrichment.to_excel(
            writer, sheet_name="Gene Target Enrichment", index=False
        )
        # KO Divergence
        ko_divergence.to_excel(writer, sheet_name="KO Divergence", index=False)
        # IMAT Divergence
        imat_div.to_excel(writer, sheet_name="IMAT Divergence", index=False)


collate_simulation_results()
