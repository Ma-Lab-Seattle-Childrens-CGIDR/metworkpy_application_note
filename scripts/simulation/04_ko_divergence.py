"""
Script to perform knock-out divergence analysis on the simulation model
"""

# Setup
# Imports
# Standard Library Imports
from collections import defaultdict
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
    BASE_PATH = pathlib.Path(".")  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
MODEL_PATH = BASE_PATH / "models"
RESULTS_PATH = BASE_PATH / "results" / "ko_divergence"

# Make directories if needed
RESULTS_PATH.mkdir(parents=True, exist_ok=True)

# Script Configuration
# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)
cobra.Configuration().solver = CONFIG["cobra"]["solver"]

# Read in the simulation model
sim_model = metworkpy.read_model(MODEL_PATH / "simulation_model.json")

# Find which reactions to not include in the metabolite networks
rxns_to_ignore_set = set()
for rxn in sim_model.reactions:
    if rxn.subsystem in CONFIG["to-ignore"]["subsystems"]:
        rxns_to_ignore_set.add(rxn.id)
rxns_to_ignore = list(rxns_to_ignore_set)  # Pandas wants a list for drop

# Generate the metabolite synthesis network
metabolite_synthesis_network = find_metabolite_synthesis_network_reactions(
    model=sim_model,
    method="essential",
    essential_proportion=CONFIG["metabolite-networks"]["essential-proportion"],
    progress_bar=False,
).drop(rxns_to_ignore, axis=0)

# Generate the metabolite consuming network dataframe
metabolite_consuming_network = find_metabolite_consuming_network_reactions(
    model=sim_model,
    reaction_proportion=CONFIG["metabolite-networks"]["reaction-proportion"],
    check_reverse=True,
    progress_bar=False,
).drop(rxns_to_ignore, axis=0)

# Create a dictionary of divergence targets to evaluate
divergence_targets: dict[str, list[str]] = defaultdict(list)

# Iterate through the reactions to create the divergence targets for subsystems
# and reactions
for rxn in sim_model.reactions:
    if rxn.subsystem in CONFIG["to-ignore"]["subsystems"]:
        continue
    divergence_targets[f"subsystem__{rxn.subsystem}"].append(rxn.id)
    divergence_targets["subsystem__whole_metabolism"].append(rxn.id)
    divergence_targets[f"reaction__{rxn.id}"].append(rxn.id)

# Iterate through the metabolite networks to create divergence targets for
# metabolite networks
for met, rxn_series in metabolite_synthesis_network.items():
    divergence_targets[f"metabolite_synthesis__{met}"] = list(
        rxn_series[rxn_series].index
    )
for met, rxn_series in metabolite_consuming_network.items():
    divergence_targets[f"metabolite_consuming__{met}"] = list(
        rxn_series[rxn_series].index
    )

# Filter the divergence targets for any that are of size 0
# (metabolite A for example has a synthesis network of size 0)
divergence_targets = {
    k: v for k, v in divergence_targets.items() if len(v) > 0
}

# Perform the KO divergence analysis
ko_divergence_df = metworkpy.divergence.ko_divergence(
    model=sim_model,
    genes_to_ko=sorted(sim_model.genes.list_attr("id")),
    target_networks=divergence_targets,
    divergence_metric=CONFIG["divergence"]["type"],
    n_neighbors=CONFIG["divergence"]["n-neighbors"],
).clip(lower=0)  # Divergence should be >0, but
# this is an estimate, so it can be slightly negative
# Clipping to correct this somewhat


# Save the results of the KO divergence
ko_divergence_df.to_csv(RESULTS_PATH / "ko_divergence_results.csv", index=True)
