"""
Script to get information about reactions from the GSMM
"""

# Standard Library Imports
import logging
import pathlib
import sys
import tomllib
import warnings

# External Imports
import cobra  # type:ignore
import metworkpy  # type:ignore
import pandas as pd

# Local Imports


# Path SEtup
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".")  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
DATA_PATH = BASE_PATH / "data"
CACHE_PATH = BASE_PATH / "cache" / "model_information"
BASE_MODEL_PATH = BASE_PATH / "models" / "iEK1011_v2_7H9_ADC_glycerol.json"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"

# Create directories if needed
CACHE_PATH.mkdir(parents=True, exist_ok=True)
LOG_PATH.mkdir(parents=True, exist_ok=True)

# Logging setup
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=LOG_PATH / "00_reaction_info.log",
    filemode="w",
    level=logging.INFO,
)

# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)

# Script Parameters
cobra.Configuration().solver = CONFIG["cobra"]["solver"]


BASE_MODEL = metworkpy.read_model(BASE_MODEL_PATH)
reaction_info_df = pd.DataFrame(
    "",
    index=pd.Index(BASE_MODEL.reactions.list_attr("id")),
    columns=pd.Index(
        [
            "name",
            "annotation",
            "expression",
            "gpr",
            "notes",
            "id",
            "subsystem",
            "genes",
        ]
    ),
)
for rxn in BASE_MODEL.reactions:
    reaction_info_df.loc[rxn.id, "name"] = rxn.name
    reaction_info_df.loc[rxn.id, "annotation"] = str(rxn.annotation)
    reaction_info_df.loc[rxn.id, "expression"] = rxn.build_reaction_string()
    reaction_info_df.loc[rxn.id, "gpr"] = rxn.gene_reaction_rule
    reaction_info_df.loc[rxn.id, "notes"] = str(rxn.notes)
    reaction_info_df.loc[rxn.id, "id"] = rxn.id
    reaction_info_df.loc[rxn.id, "subsystem"] = rxn.subsystem
    reaction_info_df.loc[rxn.id, "genes"] = ",".join([g.id for g in rxn.genes])
reaction_info_df = reaction_info_df.reset_index(drop=True)

# Save the reaction information
logger.info("Saving reaction information")
reaction_info_df.to_csv(CACHE_PATH / "reaction_information.csv", index=False)
