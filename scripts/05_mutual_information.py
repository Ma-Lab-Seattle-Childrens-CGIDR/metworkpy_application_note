"""
Script to determine flux sample based mutual information between the reactions
of the simulation model, and then determine the eigenvalue centrality
of this mutual information network
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
import networkx as nx
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
CACHE_PATH = BASE_PATH / "cache"
RESULTS_PATH = BASE_PATH / "results" / "mutual_information"
FLUX_SAMPLE_PATH = CACHE_PATH / "flux_samples"

# Make directories if needed
RESULTS_PATH.mkdir(parents=True, exist_ok=True)
CACHE_PATH.mkdir(parents=True, exist_ok=True)
FLUX_SAMPLE_PATH.mkdir(parents=True, exist_ok=True)

# Script Configuration
# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)
cobra.Configuration().solver = CONFIG["cobra"]["solver"]

# Read in the simulation model
sim_model = metworkpy.read_model(MODEL_PATH / "simulation_model.json")

# Create a list of reactions to not include in the centrality calculations
rxns_to_ignore_set = set()
for rxn in sim_model.reactions:
    if rxn.subsystem in CONFIG["to-ignore"]["subsystems"]:
        rxns_to_ignore_set.add(rxn.id)
rxns_to_ignore: list[str] = list(rxns_to_ignore_set)

# Perform flux sampling (if needed, otherwise just read in the samples)
flux_sample_cache_path = FLUX_SAMPLE_PATH / "flux_samples.parquet"
if flux_sample_cache_path.exists():
    flux_samples: pd.DataFrame = pd.read_parquet(flux_sample_cache_path)
else:
    sampler = cobra.sampling.OptGPSampler(
        model=sim_model,
        thinning=CONFIG["flux-sampling"]["thinning"],
        processes=CONFIG["processes"],
    )
    # Generate samples
    flux_samples = sampler.sample(CONFIG["flux-sampling"]["num-samples"])
    # Validate samples
    flux_samples = flux_samples[sampler.validate(flux_samples) == "v"]
    # Save the flux samples
    flux_samples.to_parquet(
        flux_sample_cache_path,
    )

# Calculate the mutual information network (if needed)
mi_adj_mat_out_path = RESULTS_PATH / "mi_adjacency.csv"
if mi_adj_mat_out_path.exists():
    mi_adj_mat = pd.read_csv(mi_adj_mat_out_path, index_col=0)
else:
    mi_adj_mat = metworkpy.information.mi_pairwise(
        flux_samples,
        processes=CONFIG["processes"],
        progress_bar=False,
        n_neighbors=CONFIG["mutual-information"]["n-neighbors"],
        metric_x=CONFIG["mutual-information"]["x-metric"],
        metric_y=CONFIG["mutual_information"]["y-metric"],
        truncate=True,
    )
    mi_adj_mat.to_csv(mi_adj_mat_out_path, index=True)

# Create the mutual information network
mi_adj_mat = mi_adj_mat.drop(rxns_to_ignore, axis=0).drop(
    rxns_to_ignore, axis=1
)
mi_network = nx.from_pandas_adjacency(mi_adj_mat)

# Find the pagerank and eigenvalue centrality
eigenvector_centrality = pd.Series(
    nx.eigenvector_centrality(mi_network, weight="weight")
)
pagerank_centrality = pd.Series(nx.pagerank(mi_network, weight="weight"))

# Create a dataframe for saving the centrality values
mi_centrality_df = pd.DataFrame(
    {"eigenvector": eigenvector_centrality, "pagerank": pagerank_centrality}
)

# Save the centrality dataframe
mi_centrality_df.to_csv(RESULTS_PATH / "mi_centrality.csv", index=True)
