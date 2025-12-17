"""
Script to simulate divergence following creation of an iMAT model
"""

# Setup
# Imports
# Standard Library Imports
import pathlib
import sys
import tomllib
from typing import Hashable

# External Imports
import cobra  # type:ignore
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
RESULTS_PATH = BASE_PATH / "results" / "iMAT"

# Create directories if needed
RESULTS_PATH.mkdir(parents=True, exist_ok=True)

# Script Configuration
# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)
cobra.Configuration().solver = CONFIG["cobra"]["solver"]

# Read in the simulation model
sim_model = metworkpy.read_model(MODEL_PATH / "simulation_model.json")

# Create a gene weights series
gene_weights = pd.Series(0.0, index=pd.Index(sim_model.genes.list_attr("id")))
gene_weights[CONFIG["imat"]["up-regulated-genes"]] = 1.0
gene_weights[CONFIG["imat"]["down-regulated-genes"]] = -1.0

# Convert the gene weights into reaction weights
reaction_weights = metworkpy.gpr.gene_to_rxn_weights(
    model=sim_model,
    gene_weights=gene_weights,
    fill_val=0.0,
    # This is also the default, just being explicit
    fn_dict=metworkpy.gpr.gpr_functions.IMAT_FUNC_DICT,
)

# Run the basic iMAT algorithm
simple_imat_model = metworkpy.imat.add_imat_constraints(
    model=sim_model,
    rxn_weights=reaction_weights,
    epsilon=CONFIG["imat"]["epsilon"],
    threshold=CONFIG["imat"]["threshold"],
)
metworkpy.imat.add_imat_objective_(
    model=simple_imat_model, rxn_weights=reaction_weights
)
imat_solution: cobra.core.Solution = simple_imat_model.optimize()
# Create a series with 1.0 for active reactions, and -1.0 for inactive
# reactions
imat_activity_series = pd.Series(
    0.0, index=pd.Index(sim_model.reactions.list_attr("id"))
)
imat_fluxes_pos_weight = imat_solution.fluxes[reaction_weights > 0.0]
imat_fluxes_neg_weight = imat_solution.fluxes[reaction_weights < 0.0]
imat_activity_series[
    imat_fluxes_pos_weight[
        (imat_fluxes_pos_weight.abs() >= CONFIG["imat"]["epsilon"])
    ].index
] = 1.0
imat_activity_series[
    imat_fluxes_neg_weight[
        (imat_fluxes_neg_weight.abs() <= CONFIG["imat"]["threshold"])
    ].index
] = -1.0
imat_activity_series.name = "IMAT Activity"
imat_activity_series.to_csv(RESULTS_PATH / "imat_activity.csv")


# Generate an iMAT model based on the reaction weights
imat_model = metworkpy.imat.generate_model(
    model=sim_model,
    rxn_weights=reaction_weights,
    method="fva",
    epsilon=CONFIG["imat"]["epsilon"],
    threshold=CONFIG["imat"]["threshold"],
    objective_tolerance=CONFIG["imat"]["objective-tolerance"],
    loopless=False,  # This can error on larger models
    processes=CONFIG["processes"],
)

# Create the sampler for the both the base model and the imat model
base_sampler = cobra.sampling.OptGPSampler(
    model=sim_model,
    thinning=CONFIG["flux-sampling"]["thinning"],
    processes=CONFIG["processes"],
)
imat_sampler = cobra.sampling.OptGPSampler(
    model=imat_model,
    thinning=CONFIG["flux-sampling"]["thinning"],
    processes=CONFIG["processes"],
)

# Sample using the sampler
base_samples = base_sampler.sample(
    CONFIG["flux-sampling"]["num-samples"], fluxes=True
)
imat_samples = imat_sampler.sample(
    CONFIG["flux-sampling"]["num-samples"], fluxes=True
)

# Validate the samples to ensure they all meet the model constraints
base_samples = base_samples[base_sampler.validate(base_samples) == "v"]
imat_samples = imat_samples[imat_sampler.validate(imat_samples) == "v"]

# Calculate the divergence between the imat model and the base model
# for all the reactions in the model, and all the metabolites in the model

# Read in the Metabolite networks for the model
metabolite_networks = pd.read_csv(
    RESULTS_PATH.parent
    / "metabolite_networks"
    / "metabolite_synthesis_network.csv",
    index_col=0,
)

# Create a dictionary for divergence groups
divergence_groups: dict[str, list[Hashable]] = {}

# For each metabolite, create a divergence group for it
for metabolite, met_network in metabolite_networks.items():
    divergence_groups[f"metabolite__{metabolite}"] = list(
        met_network[met_network].index
    )

# Also, add in divergence groups for each reaction
for rxn in metabolite_networks.index:
    divergence_groups[f"reaction__{rxn}"] = [rxn]

# Filter out any zero-size divergence groups
divergence_groups = {k: v for k, v in divergence_groups.items() if len(v) > 0}


# Compute the divergence for all these groups
imat_divergence = (
    metworkpy.divergence.group_divergence.calculate_divergence_grouped(
        dataset1=base_samples,
        dataset2=imat_samples,
        divergence_groups=divergence_groups,
        divergence_type="kl",
        processes=CONFIG["processes"],
    )
).clip(lower=0.0)
imat_divergence.name = "IMAT Divergence"


# Save the results
imat_divergence.to_csv(RESULTS_PATH / "imat_divergence.csv", index=True)
