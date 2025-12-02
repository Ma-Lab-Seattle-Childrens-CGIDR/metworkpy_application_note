"""
Script to generate metabolic reaction network, and perform
some centrality analysis
"""

# Setup
# Imports
# Standard Library Imports
import json
import pathlib
import sys
import tomllib

# External Imports
import cobra  # type: ignore
import metworkpy
import networkx as nx
import pandas as pd


# Local Imports

# Path Setup
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".")  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent
MODEL_PATH = BASE_PATH / "models"
RESULTS_PATH = BASE_PATH / "results" / "metabolic_networks"

# Create directories if needed
RESULTS_PATH.mkdir(parents=True, exist_ok=True)

# Script Configuration
# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)
cobra.Configuration().solver = CONFIG["cobra"]["solver"]

# Read in the simulation Model
sim_model = metworkpy.read_model(MODEL_PATH / "simulation_model.json")

# Create a list of reaction to remove
rxns_to_ignore_set = set()
for rxn in sim_model.reactions:
    if rxn.subsystem in CONFIG["to-ignore"]["subsystems"]:
        rxns_to_ignore_set.add(rxn.id)

# Additionally ignore all extracellular metabolites
met_to_ignore_set = set()
for met in sim_model.metabolites:
    if met.compartment in CONFIG["to-ignore"]["compartments"]:
        met_to_ignore_set.add(met.id)

# Combine the reactions/metabolites to ignore into one list
nodes_to_ignore = sorted(rxns_to_ignore_set | met_to_ignore_set)


# Generate the bipartite reaction model
metabolic_network = metworkpy.create_metabolic_network(
    model=sim_model,
    weighted=False,
    directed=CONFIG["metabolic-network"]["directed"],
    nodes_to_remove=nodes_to_ignore,
)

# Save the metabolic network
with open(RESULTS_PATH / "metabolic_network.json", "w") as f:
    json.dump(nx.node_link_data(metabolic_network), f)


# Project the network onto the reactions only
model_reactions = set(sim_model.reactions.list_attr("id")) - rxns_to_ignore_set

reaction_network = metworkpy.bipartite_project(
    metabolic_network,
    node_set=model_reactions,
    directed=CONFIG["metabolic-network"]["directed"],
)

# Save the reaction network
with open(RESULTS_PATH / "metabolic_reaction_network.json", "w") as f:
    json.dump(nx.node_link_data(reaction_network), f)

# Project the network onto metabolites only
model_metabolites = (
    set(sim_model.metabolites.list_attr("id")) - met_to_ignore_set
)

metabolite_network = metworkpy.bipartite_project(
    metabolic_network,
    node_set=model_metabolites,
    directed=CONFIG["metabolic-network"]["directed"],
)

# Save the metabolite network
with open(RESULTS_PATH / "metabolic_metabolite_network.json", "w") as f:
    json.dump(nx.node_link_data(metabolite_network), f)


# Perform centrality analysis
reaction_closeness_centrality = pd.Series(
    nx.closeness_centrality(reaction_network)
)
reaction_betweenness_centrality = pd.Series(
    nx.betweenness_centrality(reaction_network)
)
metabolite_closeness_centrality = pd.Series(
    nx.closeness_centrality(metabolite_network)
)
metabolite_betweenness_centrality = pd.Series(
    nx.betweenness_centrality(metabolite_network)
)

# Create dataframes to save results
reaction_centrality = pd.DataFrame(
    {
        "closeness": reaction_closeness_centrality,
        "betweenness": reaction_betweenness_centrality,
    }
)

metabolite_centrality = pd.DataFrame(
    {
        "closeness": metabolite_closeness_centrality,
        "betweenness": metabolite_betweenness_centrality,
    }
)

# Save the results
reaction_centrality.to_csv(
    RESULTS_PATH / "reaction_centrality.csv", index=True
)
metabolite_centrality.to_csv(
    RESULTS_PATH / "metabolite_centrality.csv", index=True
)
