"""
Script to generate visualizations for the previous
results on the simulation model
"""

# Setup
# Imports
# Standard Library Imports
import pathlib
import sys
import tomllib

# External Imports
import cobra  # type:ignore
from metabolic_modeling_utils.escher_maps import escher_map_add_data
import metworkpy
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

# Script Configuration
# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)
cobra.Configuration().solver = CONFIG["cobra"]["solver"]

# Read in the simulation model
sim_model = metworkpy.read_model(MODEL_PATH / "simulation_model.json")


escher_maps.escher_map_add_data(
    map_in_path,
    ESCHER_MAPS_OUT_DIR,
    output_prefix=f"{tf}_",
    reaction_data=reaction_divergence_df.loc[tf]
    .replace([np.inf, -np.inf], np.nan)
    .dropna(),
    metabolite_data=metabolite_divergence_df.loc[tf]
    .replace([np.inf, -np.inf], np.nan)
    .dropna(),
    reaction_scale=[
        {"type": "max", "color": "red", "size": 50},
        {"type": "min", "color": "blue", "size": 50},
        {"type": "value", "value": 0, "color": "grey", "size": 10},
    ],
    metabolite_scale=[
        {"type": "max", "color": "red", "size": 50},
        {"type": "min", "color": "blue", "size": 50},
        {"type": "value", "value": 0, "color": "grey", "size": 10},
    ],
)
