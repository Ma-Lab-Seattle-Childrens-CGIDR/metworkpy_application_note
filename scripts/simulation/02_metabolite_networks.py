"""
Script to generate metabolite networks from the simulation model
"""

# Setup
# Imports
# Standard Library Imports
import pathlib
import sys
import tomllib

# External Imports
import cobra  # type: ignore
import metworkpy
from metworkpy.metabolites.metabolite_network import (
    find_metabolite_synthesis_network_reactions,
    find_metabolite_consuming_network_reactions,
)

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
RESULTS_PATH = BASE_PATH / "results" / "metabolite_networks"

# Create directories if needed
RESULTS_PATH.mkdir(parents=True, exist_ok=True)

# Script Configuration
# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)
cobra.Configuration().solver = CONFIG["cobra"]["solver"]


# Read in the Simulation Model
sim_model = metworkpy.read_model(MODEL_PATH / "simulation_model.json")

# Create a list of reactions which will be removed from the metabolite network
# (These will be the exchange pseudo reactions, and the biomass reaction )
rxns_to_ignore_set = set()
for rxn in sim_model.reactions:
    if rxn.subsystem in CONFIG["to-ignore"]["subsystems"]:
        rxns_to_ignore_set.add(rxn.id)
rxns_to_ignore = list(rxns_to_ignore_set)  # Pandas wants a list for drop

# Generate the metabolite synthesis network dataframe
metabolite_synthesis_network = find_metabolite_synthesis_network_reactions(
    model=sim_model,
    method="essential",
    essential_proportion=CONFIG["metabolite-networks"]["essential-proportion"],
    progress_bar=False,
).drop(rxns_to_ignore, axis=0)

# Save the metabolite synthesis network dataframe
metabolite_synthesis_network.to_csv(
    RESULTS_PATH / "metabolite_synthesis_network.csv", index=True
)

# Generate the metabolite consuming network dataframe
metabolite_consuming_network = find_metabolite_consuming_network_reactions(
    model=sim_model,
    reaction_proportion=CONFIG["metabolite-networks"]["reaction-proportion"],
    check_reverse=True,
    progress_bar=False,
).drop(rxns_to_ignore, axis=0)

# Save the metabolite consuming network dataframe
metabolite_consuming_network.to_csv(
    RESULTS_PATH / "metabolite_consuming_network.csv"
)
