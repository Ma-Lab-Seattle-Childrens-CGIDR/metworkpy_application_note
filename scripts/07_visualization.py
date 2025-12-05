"""
Script to generate visualizations for the previous
results on the simulation model
"""

# Setup
# Imports
# Standard Library Imports
import json
import pathlib
import sys
import tomllib

# External Imports
import cobra  # type:ignore
from metabolic_modeling_utils import escher_maps
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
    BASE_PATH = pathlib.Path(__file__).parent.parent
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
    for metabolite, network_series in metabolite_network_df.iterrows():
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
density_escher_out_dir = RESULTS_PATH / "target_density" / "escher_maps"
density_escher_out_dir.mkdir(parents=True, exist_ok=True)

for radius, density_series in target_density_df.items():
    escher_maps.escher_map_add_data(
        input_map=ESCHER_MAP_PATH,
        output_dir=density_escher_out_dir,
        output_prefix=f"density_r{radius}",
        reaction_data=density_series,
        reaction_scale=[
            {"type": "min", "color": "black", "size": 25},
            {"type": "max", "color": "red", "size": 25},
        ],
    )
