"""
Script to generate a cobra model for use in performing simulations for
the various methods of metworkpy
"""

# Setup
# Imports
# Standard Library Imports
import itertools
import pathlib  # Handle paths
from pprint import pprint
from string import ascii_uppercase
import sys  # Used to check if running in REPL or from file
import tempfile

# External Imports
from cobra import Model, Reaction, Metabolite, io  # type:ignore
import metworkpy

if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".")  # Use current dir as base path
else:
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent
MODEL_OUT_PATH = BASE_PATH / "models"

# Make directories if needed
MODEL_OUT_PATH.mkdir(parents=True, exist_ok=True)

# Define Global Variables
EXTERNAL_EXCHANGE_BOUNDS = (-100, 100)
INTERNAL_REVERSIBLE = (-50, 50)
INTERNAL_FORWARD = (0, 50)
IMPORT_REVERSIBLE = (-50, 50)
EXPORT_IRREVERSIBLE = (0, 50)

# Create the model
sim_model = Model("simulation_model")

# Create metabolites for the model
metabolite_dict: dict[str, Metabolite] = {}
# External Metabolites
for met in itertools.chain(
    ascii_uppercase[: ascii_uppercase.find("G") + 1],
    ["N", "R", "U", "V", "W", "X"],
):
    metabolite_dict[f"{met}_E"] = Metabolite(
        f"{met}_E",
        formula=f"{met}",
        name=f"Extracellular {met}",
        compartment="E",
        charge=0,
    )

# Internal Metabolites
for met in ascii_uppercase[: ascii_uppercase.find("X") + 1]:
    metabolite_dict[f"{met}_C"] = Metabolite(
        f"{met}_C",
        formula=f"{met}",
        name=f"Cytosolic {met}",
        compartment="C",
        charge=0,
    )

# Create reactions for the model
reaction_list: list[Reaction] = []

# External Exchange Reactions
for met in itertools.chain(
    ascii_uppercase[: ascii_uppercase.find("G") + 1],
    ["N", "R", "U", "V", "W", "X"],
):
    rxn = Reaction(
        id=f"{met}_E_ex",
        name=f"External {met} exchange",
        subsystem="External Exchange Reactions",
        lower_bound=EXTERNAL_EXCHANGE_BOUNDS[0],
        upper_bound=EXTERNAL_EXCHANGE_BOUNDS[1],
    )
    rxn.add_metabolites({metabolite_dict[f"{met}_E"]: -1.0})  # Met <->
    reaction_list.append(rxn)

# Reversible Import Reactions
for met in ascii_uppercase[: ascii_uppercase.find("G") + 1]:
    rxn = Reaction(
        id=f"{met}_import",
        name=f"External {met} import",
        subsystem="Transport",
        lower_bound=IMPORT_REVERSIBLE[0],
        upper_bound=IMPORT_REVERSIBLE[1],
    )
    rxn.add_metabolites(
        {
            metabolite_dict[f"{met}_E"]: -1.0,
            metabolite_dict[f"{met}_C"]: 1.0,
        }
    )  # Met External <-> Met Internal
    reaction_list.append(rxn)

# Irreversible Export Reactions
export_met_to_subsys = {}
for subsys, met_set in {
    "S3": {"N"},
    "S5": {"W", "R"},
    "S7": {"X", "V", "U"},
}.items():
    for m in met_set:
        export_met_to_subsys[m] = subsys


irreversible_export_dict: dict[str, Reaction] = {}
for met in ["N", "R", "U", "V", "W", "X"]:
    rxn = Reaction(
        id=f"{met}_exp",
        name=f"{met} export",
        subsystem=export_met_to_subsys[met],
        lower_bound=EXPORT_IRREVERSIBLE[0],
        upper_bound=EXPORT_IRREVERSIBLE[1],
    )
    rxn.add_metabolites(
        {
            metabolite_dict[f"{met}_C"]: -1.0,
            metabolite_dict[f"{met}_E"]: 1.0,
        }
    )
    irreversible_export_dict[f"{met}_exp"] = rxn

# Add in genes for some of the export reactions
irreversible_export_dict["W_exp"].gene_reaction_rule = "g017"
irreversible_export_dict["R_exp"].gene_reaction_rule = "g018"
irreversible_export_dict["X_exp"].gene_reaction_rule = "g021"
irreversible_export_dict["V_exp"].gene_reaction_rule = "g022"
irreversible_export_dict["U_exp"].gene_reaction_rule = "g023"

# Add all the export reactions to the reaction list
for rxn in irreversible_export_dict.values():
    reaction_list.append(rxn)


# Internal Subsystem 1
r_A_B__G_H = Reaction(
    id="R_A_B__G_H",
    name="Reaction A+B<->G+H",
    subsystem="S1",
    lower_bound=INTERNAL_REVERSIBLE[0],
    upper_bound=INTERNAL_REVERSIBLE[1],
)
r_A_B__G_H.add_metabolites(
    {
        metabolite_dict["A_C"]: -1.0,
        metabolite_dict["B_C"]: -1.0,
        metabolite_dict["G_C"]: 1.0,
        metabolite_dict["H_C"]: 1.0,
    }
)
r_A_B__G_H.gene_reaction_rule = "g001"
reaction_list.append(r_A_B__G_H)

# Internal Subsystem 2
# C + H <-> I
r_C_H__I = Reaction(
    id="R_C_H__I",
    name="Reaction C+H<->I",
    subsystem="S2",
    lower_bound=INTERNAL_REVERSIBLE[0],
    upper_bound=INTERNAL_REVERSIBLE[1],
)
r_C_H__I.add_metabolites(
    {
        metabolite_dict["C_C"]: -1.0,
        metabolite_dict["H_C"]: -1.0,
        metabolite_dict["I_C"]: 1.0,
    }
)
r_C_H__I.gene_reaction_rule = "( g002 or g003 )"
reaction_list.append(r_C_H__I)

# C + D <-> J
r_C_D__J = Reaction(
    id="R_C_D__J",
    name="Reaction C+D<->J",
    subsystem="S2",
    lower_bound=INTERNAL_REVERSIBLE[0],
    upper_bound=INTERNAL_REVERSIBLE[1],
)
r_C_D__J.add_metabolites(
    {
        metabolite_dict["C_C"]: -1.0,
        metabolite_dict["D_C"]: -1.0,
        metabolite_dict["J_C"]: 1.0,
    }
)
r_C_D__J.gene_reaction_rule = "g004"
reaction_list.append(r_C_D__J)

# I <-> P
r_I__P = Reaction(
    id="R_I__P",
    name="Reaction I<->P",
    subsystem="S2",
    lower_bound=INTERNAL_REVERSIBLE[0],
    upper_bound=INTERNAL_REVERSIBLE[1],
)
r_I__P.add_metabolites(
    {
        metabolite_dict["I_C"]: -1.0,
        metabolite_dict["P_C"]: 1.0,
    }
)
r_I__P.gene_reaction_rule = "g008"
reaction_list.append(r_I__P)


# J <-> Q
r_J__Q = Reaction(
    id="R_J__Q",
    name="Reaction J<->Q",
    subsystem="S2",
    lower_bound=INTERNAL_REVERSIBLE[0],
    upper_bound=INTERNAL_REVERSIBLE[1],
)
r_J__Q.add_metabolites(
    {
        metabolite_dict["J_C"]: -1.0,
        metabolite_dict["Q_C"]: 1.0,
    }
)
r_J__Q.gene_reaction_rule = "g010"
reaction_list.append(r_J__Q)


# Create a biomass reaction
r_biomass = Reaction(
    id="biomass",
    name="Biomass",
    subsystem="Biomass",
    lower_bound=INTERNAL_FORWARD[0],
    upper_bound=INTERNAL_FORWARD[1],
)
r_biomass.add_metabolites(
    {
        metabolite_dict["W_C"]: -1.0,
        metabolite_dict["X_C"]: -1.0,
        metabolite_dict["U_C"]: -1.0,
    }
)
reaction_list.append(r_biomass)

# Add all the reactions to the model
sim_model.add_reactions(reaction_list=reaction_list)

# Set the objective to be the biomass reaction
sim_model.objective = "biomass"

# Write the model to a temporary file
# in order to validate the sbml
with tempfile.NamedTemporaryFile(suffix=".xml") as f_sbml:
    io.write_sbml_model(sim_model, f_sbml.name)
    _, report = io.validate_sbml_model(f_sbml.name)
for err_type, err_list in report.items():
    if len(err_list) > 0:
        raise RuntimeError(
            f"Error in model-- Type: {err_type}, Errors: {err_list}"
        )


# Save the model in various formats
metworkpy.write_model(sim_model, MODEL_OUT_PATH / "simulation_model.json")
metworkpy.write_model(sim_model, MODEL_OUT_PATH / "simulation_model.xml")
metworkpy.write_model(sim_model, MODEL_OUT_PATH / "simulation_model.mat")
metworkpy.write_model(sim_model, MODEL_OUT_PATH / "simulation_model.yaml")
