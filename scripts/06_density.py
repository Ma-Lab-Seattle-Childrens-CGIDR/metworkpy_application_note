"""
Script to simulate gene target density on the simulation
model
"""

# Setup
# Imports
# Standard Library Imports
import pathlib
import sys
import tomllib

# External Imports
import cobra  # type:ignore
import metworkpy
from metworkpy.network.density import (
    gene_target_density,
    gene_target_enrichment,
)
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
RESULTS_PATH = BASE_PATH / "results" / "target_density"

# Create directories if needed
RESULTS_PATH.mkdir(parents=True, exist_ok=True)

# Script Configuration
# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)
cobra.Configuration().solver = CONFIG["cobra"]["solver"]

# Read in the simulation model
sim_model = metworkpy.read_model(MODEL_PATH / "simulation_model.json")

# Construct the reaction network from the model
metabolic_network = metworkpy.create_metabolic_network(
    model=sim_model,
    weighted=False,
    directed=CONFIG["metabolic-network"]["directed"],
)
reaction_network = metworkpy.network.bipartite_project(
    metabolic_network,
    directed=CONFIG["metabolic-network"]["directed"],
    node_set=sim_model.reactions.list_attr("id"),
)


# Iterate through the radius list and find the density for each radius
density_series_list: list[pd.Series] = []
for radius in CONFIG["target-density"]["radius-list"]:
    # Create the series of target density
    target_density_series = gene_target_density(
        metabolic_network=reaction_network,
        metabolic_model=sim_model,
        gene_labels=CONFIG["target-density"]["targeted-genes"],
        radius=radius,
    )
    target_density_series.name = f"Radius: {radius}"
    density_series_list.append(target_density_series)
# Concatenate together the different series
density_df = pd.concat(density_series_list, axis=1)

# Iterate through the radius list and find the enrichment p-value
# for each radius
enrichment_pval_series_list: list[pd.Series] = []
for radius in CONFIG["target-density"]["radius-list"]:
    # Create the series of target density
    target_enrichment_pval_series = gene_target_enrichment(
        metabolic_network=reaction_network,
        metabolic_model=sim_model,
        gene_targets=CONFIG["target-density"]["targeted-genes"],
        metric="p-value",
        alternative="greater",
        radius=radius,
    )
    target_enrichment_pval_series.name = f"Radius: {radius}"
    enrichment_pval_series_list.append(target_enrichment_pval_series)
# Concatenate together the different series
enrichment_pval_df = pd.concat(enrichment_pval_series_list, axis=1)

# Iterate through the radius list and find the enrichment odds-ratio for
# each radius
enrichment_odds_series_list: list[pd.Series] = []
for radius in CONFIG["target-density"]["radius-list"]:
    print(f"radius: {radius}")
    # Create the series of target density
    target_enrichment_odds_series = gene_target_enrichment(
        metabolic_network=reaction_network,
        metabolic_model=sim_model,
        gene_targets=CONFIG["target-density"]["targeted-genes"],
        metric="odds-ratio",
        alternative="greater",
        radius=radius,
    )
    target_enrichment_odds_series.name = f"Radius: {radius}"
    enrichment_odds_series_list.append(target_enrichment_odds_series)
# Concatenate together the different series
enrichment_odds_df = pd.concat(enrichment_odds_series_list, axis=1)

# Save the results
density_df.to_csv(RESULTS_PATH / "gene_target_density.csv")
enrichment_pval_df.to_csv(RESULTS_PATH / "gene_target_enrichment_pval.csv")
enrichment_odds_df.to_csv(RESULTS_PATH / "gene_target_enrichment_odds.csv")
