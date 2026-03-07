"""
Script to generate a cobra model for use in performing simulations for
the various methods of metworkpy
"""

# Setup
# Imports
# Standard Library Imports
import itertools  # Used for chaining iterators
import pathlib  # Handle paths
from string import ascii_uppercase  # Used for iteration to create metabolites
import sys  # Used to check if running in REPL or from file
import tempfile  # Used to create temporary file for model validation
import tomllib
from typing import cast

# External Imports
from cobra import Configuration, Model, Reaction, Metabolite, io
import metworkpy  # Used for convienience function for writing the model
import pandas as pd

if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".")  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
MODEL_OUT_PATH = BASE_PATH / "models"

# Make directories if needed
MODEL_OUT_PATH.mkdir(parents=True, exist_ok=True)

# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)
# Set the cobra solver
Configuration.solver = CONFIG["cobra"]["solver"]

######################
# Create the model ###
######################
sim_model = Model("simulation_model")

########################
# Create Metabolites ###
########################
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

######################
# Create Reactions ###
######################
reaction_list: list[Reaction] = []

#################################
# External Exchange Reactions ###
#################################

for met in itertools.chain(
    ascii_uppercase[: ascii_uppercase.find("G") + 1],
    ["N", "R", "U", "V", "W", "X"],
):
    rxn = Reaction(
        id=f"{met}_E_ex",
        name=f"External {met} exchange",
        subsystem="External Exchange Reactions",
        lower_bound=CONFIG["simulation"]["reaction-bounds"][
            "external-exchange"
        ][0],
        upper_bound=CONFIG["simulation"]["reaction-bounds"][
            "external-exchange"
        ][1],
    )
    rxn.add_metabolites({metabolite_dict[f"{met}_E"]: -1.0})  # Met <->
    reaction_list.append(rxn)

#################################
# Reversible Import Reactions ###
#################################
for met in ascii_uppercase[: ascii_uppercase.find("G") + 1]:
    rxn = Reaction(
        id=f"{met}_import",
        name=f"External {met} import",
        subsystem="Transport",
        lower_bound=CONFIG["simulation"]["reaction-bounds"][
            "import-reversible"
        ][0],
        upper_bound=CONFIG["simulation"]["reaction-bounds"][
            "import-reversible"
        ][1],
    )
    rxn.add_metabolites(
        {
            metabolite_dict[f"{met}_E"]: -1.0,
            metabolite_dict[f"{met}_C"]: 1.0,
        }
    )  # Met External <-> Met Internal
    reaction_list.append(rxn)

###################################
# Irreversible Export Reactions ###
###################################
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
        lower_bound=CONFIG["simulation"]["reaction-bounds"][
            "export-irreversible"
        ][0],
        upper_bound=CONFIG["simulation"]["reaction-bounds"][
            "export-irreversible"
        ][1],
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


##########################
# Internal Subsystem 1 ###
##########################
r_A_B__G_H = Reaction(
    id="R_A_B__G_H",
    name="Reaction A+B<->G+H",
    subsystem="S1",
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        0
    ],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        1
    ],
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

##########################
# Internal Subsystem 2 ###
##########################
# C + H <-> I
r_C_H__I = Reaction(
    id="R_C_H__I",
    name="Reaction C+H<->I",
    subsystem="S2",
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        0
    ],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        1
    ],
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
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        0
    ],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        1
    ],
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
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        0
    ],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        1
    ],
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
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        0
    ],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        1
    ],
)
r_J__Q.add_metabolites(
    {
        metabolite_dict["J_C"]: -1.0,
        metabolite_dict["Q_C"]: 1.0,
    }
)
r_J__Q.gene_reaction_rule = "g010"
reaction_list.append(r_J__Q)

##########################
# Internal Subsystem 3 ###
##########################
# E <-> M
r_E__M = Reaction(
    id="R_E__M",
    name="Reaction E<->M",
    subsystem="S3",
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        0
    ],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        1
    ],
)
r_E__M.add_metabolites(
    {
        metabolite_dict["E_C"]: -1.0,
        metabolite_dict["M_C"]: 1.0,
    }
)
r_E__M.gene_reaction_rule = "g005"
reaction_list.append(r_E__M)

# F + M -> N
r_F_M__N = Reaction(
    id="R_F_M__N",
    name="Reaction F+M->N",
    subsystem="S3",
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][0],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][1],
)
r_F_M__N.add_metabolites(
    {
        metabolite_dict["F_C"]: -1.0,
        metabolite_dict["M_C"]: -1.0,
        metabolite_dict["N_C"]: 1.0,
    }
)
r_F_M__N.gene_reaction_rule = "( g006 and g007 )"
reaction_list.append(r_F_M__N)

# N <-> T
r_N__T = Reaction(
    id="R_N__T",
    name="Reaction N<->T",
    subsystem="S3",
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        0
    ],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        1
    ],
)
r_N__T.add_metabolites(
    {
        metabolite_dict["N_C"]: -1.0,
        metabolite_dict["T_C"]: 1.0,
    }
)
r_N__T.gene_reaction_rule = "( g011 or g012 )"
reaction_list.append(r_N__T)

# N -> U
r_N__U = Reaction(
    id="R_N__U",
    name="Reaction N->U",
    subsystem="S3",
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][0],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][1],
)
r_N__U.add_metabolites(
    {
        metabolite_dict["N_C"]: -1.0,
        metabolite_dict["U_C"]: 1.0,
    }
)
r_N__U.gene_reaction_rule = "( g011 and g013 )"
reaction_list.append(r_N__U)

##########################
# Internal Subsystem 4 ###
##########################
# H <-> K
r_H__K = Reaction(
    id="R_H__K",
    name="Reaction H<->K",
    subsystem="S4",
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        0
    ],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        1
    ],
)
r_H__K.add_metabolites(
    {
        metabolite_dict["H_C"]: -1.0,
        metabolite_dict["K_C"]: 1.0,
    }
)
r_H__K.gene_reaction_rule = "g009"
reaction_list.append(r_H__K)

# G + K <-> L
r_G_K__L = Reaction(
    id="R_G_K__L",
    name="Reaction G+K<->L",
    subsystem="S4",
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        0
    ],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-reversible"][
        1
    ],
)
r_G_K__L.add_metabolites(
    {
        metabolite_dict["G_C"]: -1.0,
        metabolite_dict["K_C"]: -1.0,
        metabolite_dict["L_C"]: 1.0,
    }
)
r_G_K__L.gene_reaction_rule = "g008"
reaction_list.append(r_G_K__L)

# K -> O
r_K__O = Reaction(
    id="R_K__O",
    name="Reaction K->O",
    subsystem="S4",
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][0],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][1],
)
r_K__O.add_metabolites(
    {
        metabolite_dict["K_C"]: -1.0,
        metabolite_dict["O_C"]: 1.0,
    }
)
r_K__O.gene_reaction_rule = "g024"
reaction_list.append(r_K__O)

# L -> W
r_L__W = Reaction(
    id="R_L__W",
    name="Reaction L->W",
    subsystem="S4",
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][0],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][1],
)
r_L__W.add_metabolites(
    {
        metabolite_dict["L_C"]: -1.0,
        metabolite_dict["W_C"]: 1.0,
    }
)
r_L__W.gene_reaction_rule = "( g014 and g020 )"
reaction_list.append(r_L__W)

##########################
# Internal Subsystem 5 ###
##########################
# O + P -> R
r_O_P__R = Reaction(
    id="R_O_P__R",
    name="Reaction O+P->R",
    subsystem="S4",
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][0],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][1],
)
r_O_P__R.add_metabolites(
    {
        metabolite_dict["O_C"]: -1.0,
        metabolite_dict["P_C"]: -1.0,
        metabolite_dict["R_C"]: 1.0,
    }
)
r_O_P__R.gene_reaction_rule = "( g015 and g020 )"
reaction_list.append(r_O_P__R)


##########################
# Internal Subsystem 6 ###
##########################
# P + Q -> S
r_P_Q__S = Reaction(
    id="R_P_Q__S",
    name="Reaction P+Q->S",
    subsystem="S6",
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][0],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][1],
)
r_P_Q__S.add_metabolites(
    {
        metabolite_dict["P_C"]: -1.0,
        metabolite_dict["Q_C"]: -1.0,
        metabolite_dict["S_C"]: 1.0,
    }
)
r_P_Q__S.gene_reaction_rule = "g016"
reaction_list.append(r_P_Q__S)

# S + T -> V + X
r_S_T__V_X = Reaction(
    id="R_S_T__V_X",
    name="Reaction S+T->V+X",
    subsystem="S6",
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][0],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][1],
)
r_S_T__V_X.add_metabolites(
    {
        metabolite_dict["S_C"]: -1.0,
        metabolite_dict["T_C"]: -1.0,
        metabolite_dict["V_C"]: 1.0,
        metabolite_dict["X_C"]: 1.0,
    }
)
r_S_T__V_X.gene_reaction_rule = "g019"
reaction_list.append(r_S_T__V_X)


######################
# Biomass Reaction ###
######################
r_biomass = Reaction(
    id="biomass",
    name="Biomass",
    subsystem="Biomass",
    lower_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][0],
    upper_bound=CONFIG["simulation"]["reaction-bounds"]["internal-forward"][1],
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

######################
# Model Validation ###
######################

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
# Check that the optimal biomass generation is near 25.0
if abs(sim_model.slim_optimize() - 25.0) > 1e-7:
    raise ValueError("Model has incorrect biomass objective")

####################
# Save the Model ###
####################
metworkpy.write_model(sim_model, MODEL_OUT_PATH / "simulation_model.json")
metworkpy.write_model(sim_model, MODEL_OUT_PATH / "simulation_model.xml")
metworkpy.write_model(sim_model, MODEL_OUT_PATH / "simulation_model.mat")
metworkpy.write_model(sim_model, MODEL_OUT_PATH / "simulation_model.yaml")

###############################################
# Create a Table of Reaction, Equation, GPR ###
###############################################
model_df = pd.DataFrame(
    "",
    index=pd.Index(sim_model.reactions.list_attr("id")),
    columns=pd.Index(["Equation", "Gene-Reaction Rule", "Subsystem"]),
)
for rxn in sim_model.reactions:
    rxn = cast(Reaction, rxn)
    model_df.loc[rxn.id, "Equation"] = rxn.build_reaction_string()
    model_df.loc[rxn.id, "Gene-Reaction Rule"] = rxn.gene_reaction_rule
    model_df.loc[rxn.id, "Subsystem"] = rxn.subsystem

# Save the datafame to the models folder
model_df.reset_index(drop=False, names="Reaction").to_csv(
    MODEL_OUT_PATH / "simulation_model_table.csv", index=False
)
