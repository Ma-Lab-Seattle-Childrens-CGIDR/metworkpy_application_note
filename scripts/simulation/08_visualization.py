"""
Script to generate visualizations for the previous
results on the simulation model
"""

# Setup
# Imports
# Standard Library Imports
from collections import defaultdict
import json
import pathlib
import sys
import tomllib

# External Imports
import cobra  # type:ignore
import iplotx as ipx  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
from metabolic_modeling_utils import escher_maps  # type: ignore
import metworkpy
import networkx as nx
import numpy as np
import pandas as pd

# Local Imports

# Path Setup
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".").absolute()  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
MODEL_PATH = BASE_PATH / "models"
RESULTS_PATH = BASE_PATH / "results"
ESCHER_MAP_PATH = BASE_PATH / "escher_maps" / "simulation_model.json"

# Script Configuration
# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)
cobra.Configuration().solver = CONFIG["cobra"]["solver"]

# Read in the simulation model
sim_model = metworkpy.read_model(MODEL_PATH / "simulation_model.json")

#########################
# Metabolite Networks ###
#########################
# Read in the metabolite network dataframe
metabolite_synthesis_networks = pd.read_csv(
    RESULTS_PATH / "metabolite_networks" / "metabolite_synthesis_network.csv",
    index_col=0,
)
metabolite_consuming_networks = pd.read_csv(
    RESULTS_PATH / "metabolite_networks" / "metabolite_consuming_network.csv",
    index_col=0,
)

metabolite_net_escher_out_dir = (
    RESULTS_PATH / "metabolite_networks" / "escher_maps"
)
metabolite_net_escher_out_dir.mkdir(parents=True, exist_ok=True)

for net_direction, metabolite_network_df in zip(
    ["synthesis", "consuming"],
    [metabolite_synthesis_networks, metabolite_consuming_networks],
):
    for metabolite, network_series in metabolite_network_df.items():
        # Convert the series to float to work with the json export
        network_series = network_series.astype("float")
        # Add the metabolite network data to the map
        escher_maps.escher_map_add_data(
            input_map=ESCHER_MAP_PATH,
            output_dir=metabolite_net_escher_out_dir,
            output_prefix=f"{metabolite}_{net_direction}_",
            reaction_data=network_series,
            reaction_scale=[
                {"type": "value", "value": 0.0, "color": "black", "size": 25},
                {"type": "value", "value": 1.0, "color": "red", "size": 25},
            ],
        )


#################################
# Reaction Network Centrality ###
#################################

# Read in the reaction centrality results
rxn_centrality_df = pd.read_csv(
    RESULTS_PATH / "metabolic_networks" / "reaction_centrality.csv",
    index_col=0,
)
met_centrality_df = pd.read_csv(
    RESULTS_PATH / "metabolic_networks" / "metabolite_centrality.csv",
    index_col=0,
)

centrality_escher_out_dir = RESULTS_PATH / "metabolic_networks" / "escher_maps"
centrality_escher_out_dir.mkdir(parents=True, exist_ok=True)

for cent_measure in ["betweenness", "closeness"]:
    escher_maps.escher_map_add_data(
        input_map=ESCHER_MAP_PATH,
        output_dir=centrality_escher_out_dir,
        output_prefix=f"{cent_measure}_",
        reaction_data=rxn_centrality_df[cent_measure],
        metabolite_data=met_centrality_df[cent_measure],
        reaction_scale=[
            {"type": "min", "color": "blue", "size": 10},
            {"type": "max", "color": "red", "size": 30},
        ],
        metabolite_scale=[
            {"type": "min", "color": "blue", "size": 10},
            {"type": "max", "color": "red", "size": 30},
        ],
    )

###################
# KO Divergence ###
###################
# Read in the KO divergence resutlts
ko_div_res = pd.read_csv(
    RESULTS_PATH / "ko_divergence" / "ko_divergence_results.csv", index_col=0
)
# Extract the reaction divergence
ko_rxn_div_res = ko_div_res.loc[
    :, ko_div_res.columns.str.startswith("reaction__")
]
ko_rxn_div_res.columns = ko_rxn_div_res.columns.str.replace("reaction__", "")
# Extract the metabolite divergence
ko_met_div_res = ko_div_res.loc[
    :, ko_div_res.columns.str.startswith("metabolite_synthesis__")
]
ko_met_div_res.columns = ko_met_div_res.columns.str.replace(
    "metabolite_synthesis__", ""
)

# Set the output directory
ko_div_map_out_dir = RESULTS_PATH / "ko_divergence" / "escher_maps"
ko_div_map_out_dir.mkdir(parents=True, exist_ok=True)

for gene in ko_div_res.index:
    rxn_div_series = ko_rxn_div_res.loc[gene]
    met_div_series = ko_met_div_res.loc[gene]
    escher_maps.escher_map_add_data(
        input_map=ESCHER_MAP_PATH,
        output_dir=ko_div_map_out_dir,
        output_prefix=f"{gene}_ko_div_",
        reaction_data=rxn_div_series.replace(
            [np.inf, -np.inf], np.nan
        ).replace(np.nan, 0.0),
        metabolite_data=met_div_series.replace(
            [np.inf, -np.inf], np.nan
        ).replace(np.nan, 0.0),
        reaction_scale=[
            {
                "type": "min",
                "color": "blue",
                "size": 10,
            },
            {
                "type": "max",
                "color": "red",
                "size": 30,
            },
        ],
        metabolite_scale=[
            {
                "type": "min",
                "color": "blue",
                "size": 10,
            },
            {
                "type": "max",
                "color": "red",
                "size": 30,
            },
        ],
    )


########################
# Mutual Information ###
########################
# Read in the mutual information centrality results
mi_cent_df = pd.read_csv(
    RESULTS_PATH / "mutual_information" / "mi_centrality.csv", index_col=0
)

mi_cent_out_dir = RESULTS_PATH / "mutual_information" / "escher_maps"
mi_cent_out_dir.mkdir(parents=True, exist_ok=True)

escher_maps.escher_map_add_data(
    input_map=ESCHER_MAP_PATH,
    output_dir=mi_cent_out_dir,
    output_prefix="mi_centrality_",
    reaction_data=mi_cent_df["eigenvector"],
    reaction_scale=[
        {"type": "min", "color": "blue", "size": 10},
        {"type": "max", "color": "red", "size": 30},
    ],
)

#############
# Density ###
#############
# Read in the target density results
target_density_df = pd.read_csv(
    RESULTS_PATH / "target_density" / "gene_target_density.csv", index_col=0
)
target_density_df.columns = target_density_df.columns.str.replace(
    "^Radius: ", "", regex=True
)
density_escher_out_dir = RESULTS_PATH / "target_density" / "escher_maps"
density_escher_out_dir.mkdir(parents=True, exist_ok=True)

for radius, density_series in target_density_df.items():
    escher_maps.escher_map_add_data(
        input_map=ESCHER_MAP_PATH,
        output_dir=density_escher_out_dir,
        output_prefix=f"density_r{radius}_",
        reaction_data=density_series,
        reaction_scale=[
            {"type": "min", "color": "black", "size": 25},
            {"type": "max", "color": "red", "size": 25},
        ],
    )

# Read in the target enrichment p-value results
target_enrichment_pval_df = pd.read_csv(
    RESULTS_PATH / "target_density" / "gene_target_enrichment_pval.csv",
    index_col=0,
)
target_enrichment_pval_df.columns = (
    target_enrichment_pval_df.columns.str.replace("^Radius: ", "", regex=True)
)
enrichment_escher_out_dir = RESULTS_PATH / "target_density" / "escher_maps"
enrichment_escher_out_dir.mkdir(parents=True, exist_ok=True)

for radius, enrichment_series in target_enrichment_pval_df.items():
    escher_maps.escher_map_add_data(
        input_map=ESCHER_MAP_PATH,
        output_dir=enrichment_escher_out_dir,
        output_prefix=f"enrichment_pval_r{radius}_",
        reaction_data=enrichment_series,
        reaction_scale=[
            {"type": "value", "value": 1.0, "color": "black", "size": 25},
            {"type": "value", "value": 0.0, "color": "red", "size": 25},
        ],
    )

# Read in the target enrichment odds results
target_enrichment_odds_df = pd.read_csv(
    RESULTS_PATH / "target_density" / "gene_target_enrichment_odds.csv",
    index_col=0,
)
target_enrichment_odds_df.columns = (
    target_enrichment_odds_df.columns.str.replace("^Radius: ", "", regex=True)
)
enrichment_escher_out_dir = RESULTS_PATH / "target_density" / "escher_maps"
enrichment_escher_out_dir.mkdir(parents=True, exist_ok=True)

for radius, enrichment_series in target_enrichment_odds_df.items():
    # Clip the odds ratio to max instead of inf
    enrichment_series = enrichment_series.clip(
        upper=enrichment_series.replace(np.inf, np.nan).max()
    )
    escher_maps.escher_map_add_data(
        input_map=ESCHER_MAP_PATH,
        output_dir=enrichment_escher_out_dir,
        output_prefix=f"enrichment_odds_r{radius}_",
        reaction_data=enrichment_series,
        reaction_scale=[
            {"type": "value", "value": 0.0, "color": "black", "size": 25},
            {"type": "max", "color": "red", "size": 25},
        ],
    )

##########
# IMAT ###
##########
# Read in the IMAT activity results
imat_activity_series = pd.read_csv(
    RESULTS_PATH / "iMAT" / "imat_activity.csv", index_col=0
)["IMAT Activity"]
# Create directory to save the escher maps
imat_escher_map_out_dir = RESULTS_PATH / "iMAT" / "escher_maps"
imat_escher_map_out_dir.mkdir(parents=True, exist_ok=True)


escher_maps.escher_map_add_data(
    input_map=ESCHER_MAP_PATH,
    output_dir=imat_escher_map_out_dir,
    output_prefix="imat_solution_",
    reaction_data=imat_activity_series,
    reaction_scale=[
        {"type": "value", "value": -1.0, "color": "blue", "size": 10},
        {"type": "value", "value": 0.0, "color": "grey", "size": 20},
        {"type": "value", "value": 1.0, "color": "red", "size": 30},
    ],
)


# Read in the IMAT divergence results
imat_div_res = pd.read_csv(
    RESULTS_PATH / "iMAT" / "imat_divergence.csv", index_col=0
)
imat_rxn_div = imat_div_res[imat_div_res.index.str.startswith("reaction__")]
imat_rxn_div.index = imat_rxn_div.index.str.replace("reaction__", "")
imat_met_div = imat_div_res[imat_div_res.index.str.startswith("metabolite__")]
imat_met_div.index = imat_met_div.index.str.replace("metabolite__", "")

escher_maps.escher_map_add_data(
    input_map=ESCHER_MAP_PATH,
    output_dir=imat_escher_map_out_dir,
    output_prefix="imat_divergence_",
    reaction_data=imat_rxn_div["IMAT Divergence"],
    reaction_scale=[
        {"type": "min", "color": "blue", "size": 10},
        {"type": "max", "color": "red", "size": 30},
    ],
    metabolite_data=imat_met_div["IMAT Divergence"],
    metabolite_scale=[
        {"type": "min", "color": "blue", "size": 10},
        {"type": "max", "color": "red", "size": 30},
    ],
)

##################################
# Reaction Network Visualization #
##################################
REACTION_NODE_COLOR = "tab:blue"
METABOLITE_NODE_COLOR = "tab:red"
NODE_SIZE = 20
EDGE_ALPHA = 0.5
EDGE_WIDTH = 2
EDGE_COLOR = "k"
PLOT_MARGIN = 0.2
IMG_FORMAT = "svg"
FIG_WIDTH = 30
FIG_HEIGHT = 30
LAYOUT_SEED = 314
# Create (if needed) directory to save images
network_graph_viz_path = (
    RESULTS_PATH / "metabolic_networks" / "graph_visualization"
)
network_graph_viz_path.mkdir(parents=True, exist_ok=True)
sim_model_rxn_list = sim_model.reactions.list_attr("id")

# Create a function to create a graph representation


def draw_graph(
    network: nx.DiGraph | nx.Graph,
    figure_size: tuple[float, float],
    node_colors: dict,
    out_path: pathlib.Path,
    seed: int | np.random.RandomState | None = None,
    **kwargs,
):
    """
    Create a image representing a graph

    Parameters
    ----------
    network : nx.Graph or nx.DiGraph
        Network to create an image of
    figure_size : tuple of float
        Tuple describing size of the figure in inches as (width, height)
    node_colors : dict of str to list of str
        Dictionary describing node colors, the keys should be the colors,
        and the values should be lists of nodes which will have that color
    seed : int, RanomState instance or None, default=None
    kwargs
        Any additional arguments are passed to the iplotx network function
    """
    # Draw the metabolic network
    fig, ax = plt.subplots()
    fig.set_size_inches(*figure_size)
    # Add color attributes to the graph
    # Create a list of the colors for the nodes
    node_color_list = []
    for n in network.nodes:
        for color, node_set in node_colors.items():
            if n in node_set:
                node_color_list.append(color)
                break
    # Plot the network with iplotx
    ipx.network(
        network,
        ax=ax,
        layout=nx.spring_layout(metabolic_network, seed=seed),
        vertex_marker="r",
        vertex_labels=True,
        vertex_facecolor=node_color_list,
        style="hollow",
        **kwargs,
    )
    fig.savefig(out_path)
    plt.close()


# Read in the reaction networks
with open(  # type: ignore
    RESULTS_PATH / "metabolic_networks" / "metabolic_network.json", "r"
) as f:
    metabolic_network = nx.node_link_graph(json.load(f))
with open(  # type: ignore
    RESULTS_PATH / "metabolic_networks" / "metabolic_reaction_network.json",
    "r",
) as f:
    metabolic_rxn_network = nx.node_link_graph(json.load(f))
with open(  # type: ignore
    RESULTS_PATH / "metabolic_networks" / "metabolic_metabolite_network.json",
    "r",
) as f:
    metabolic_metabolite_network = nx.node_link_graph(json.load(f))

# Find colors for the nodes in the graphs
graph_node_colors = defaultdict(list)
for n in metabolic_network.nodes:
    graph_node_colors[
        REACTION_NODE_COLOR
        if n in sim_model_rxn_list
        else METABOLITE_NODE_COLOR
    ].append(n)

# Draw the metabolic network
draw_graph(
    metabolic_network,
    figure_size=(30, 30),
    node_colors=graph_node_colors,
    out_path=network_graph_viz_path / f"metabolic_network.{IMG_FORMAT}",
    seed=LAYOUT_SEED,
)

# Draw the metabolic reaction network
draw_graph(
    metabolic_rxn_network,
    figure_size=(20, 20),
    node_colors=graph_node_colors,
    out_path=network_graph_viz_path
    / f"metabolic_reaction_network.{IMG_FORMAT}",
    seed=LAYOUT_SEED,
)

# Draw the metabolic metabolite network
draw_graph(
    metabolic_metabolite_network,
    figure_size=(15, 15),
    node_colors=graph_node_colors,
    out_path=network_graph_viz_path
    / f"metabolic_metabolite_network.{IMG_FORMAT}",
    seed=LAYOUT_SEED,
)
