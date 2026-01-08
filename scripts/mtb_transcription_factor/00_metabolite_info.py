"""
Script to scrape BiGG and other databases for metabolite information
"""

# Standard Library Imports
import logging
import pathlib
import re
import sys
import tomllib
from typing import Optional
import warnings

# External Imports
import cobra  # type: ignore
import metworkpy
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
    filename=LOG_PATH / "00_metabolite_info.log",
    filemode="w",
    level=logging.INFO,
)

# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)

# Script Parameters
KEGG_PATHWAY_DESCRIPTIONS = pd.read_csv(
    "https://rest.kegg.jp/list/pathway",
    sep="\t",
    names=["pathway", "description"],
)

# Read in the base model to get information about the metabolites
logger.info("Reading in base model")
cobra.Configuration().solver = CONFIG["cobra"]["solver"]
BASE_MODEL = metworkpy.read_model(BASE_MODEL_PATH)
# Create a metabolite information dataframe using the BASE_MODEL
metabolite_info_df = pd.DataFrame(
    "",
    index=pd.Index(BASE_MODEL.metabolites.list_attr("id")),
    columns=pd.Index(["name", "compartment", "formula", "id"]),
)
logger.info("Creating initial metabolite information dataframe")
for met in BASE_MODEL.metabolites:
    metabolite_info_df.loc[met.id, "name"] = met.name
    metabolite_info_df.loc[met.id, "compartment"] = met.compartment
    metabolite_info_df.loc[met.id, "formula"] = met.formula
    metabolite_info_df.loc[met.id, "id"] = met.id

metabolite_info_df = metabolite_info_df.reset_index(drop=True)


# Read in the bigg_metabolite_df
logger.info("Reading in BiGG metabolite information")
BIGG_METABOLITE_DF = pd.read_csv(
    DATA_PATH / "bigg_info" / "bigg_models_metabolites.txt", sep="\t"
).fillna("")

# Create a dict to map from current metabolite ids to Bigg IDS
logger.info("Creating dict of metabolite to BiGG id")
metabolite_to_bigg_id_dict = {}

# Some metabolites don't exist in the BiGG database
BIGG_MISSING_METABOLITES = {
    "5ohhiccoa[c]",
    "hieccoa[c]",
    "cocheacoa[c]",
    "moodacoa[c]",
    "mbkacoa[c]",
}

MET_REGEX = re.compile(r"^([\w_]+)\[(\w)\]$")
for metabolite in metabolite_info_df["id"]:
    # If the metabolite is already in the dictionary
    # i.e. it has been hand translated, skip it
    if metabolite in metabolite_to_bigg_id_dict:
        continue
    # If the metabolite is known not to be in BiGG, Skip it
    if metabolite in BIGG_MISSING_METABOLITES:
        continue
    # Split the name and the compartment
    met_name, met_compartment = MET_REGEX.findall(metabolite)[0]
    # Create a potential ID
    potential_bigg_id = f"{met_name}_{met_compartment}"
    # Check if the metabolite is a correct BiGG id
    if (
        len(
            BIGG_METABOLITE_DF.loc[
                BIGG_METABOLITE_DF["bigg_id"] == potential_bigg_id
            ]
        )
        == 1
    ):
        metabolite_to_bigg_id_dict[metabolite] = BIGG_METABOLITE_DF.loc[
            BIGG_METABOLITE_DF["bigg_id"] == potential_bigg_id, "bigg_id"
        ].iloc[0]
        continue
    # Check if the metabolite is an old name
    if (
        len(
            BIGG_METABOLITE_DF[
                BIGG_METABOLITE_DF["old_bigg_ids"]
                .str.contains(metabolite)
                .astype("bool")
            ]
            | BIGG_METABOLITE_DF["old_bigg_ids"]
            .str.contains(potential_bigg_id)
            .astype("bool")
        )
        == 1
    ):
        new_id = (
            BIGG_METABOLITE_DF[
                BIGG_METABOLITE_DF["old_bigg_ids"].str.contains(metabolite)
            ]
            | BIGG_METABOLITE_DF["old_bigg_ids"]
            .str.contains(potential_bigg_id)["bigg_id"]
            .iloc[0]  # type: ignore
        )
        metabolite_to_bigg_id_dict[metabolite] = new_id
        continue
    warnings.warn(f"Couldn't find {metabolite}")


# Create Regexes for the different metabolite database links
METANETX_REGEX = re.compile(
    r"MetaNetX \(MNX\) Chemical: http://identifiers\.org/metanetx\.chemical/(\w+)"
)
KEGG_REGEX = re.compile(
    r"KEGG Compound: http://identifiers\.org/kegg\.compound/(C\d+)"
)
CHEBI_REGEX = re.compile(r"CHEBI: http://identifiers\.org/chebi/CHEBI:(\d+)")
SEED_REGEX = re.compile(
    r"SEED Compound: http://identifiers\.org/seed\.compound/(cpd\d+)"
)

DATABASE_REGEX_DICT = {
    "MetaNetX": METANETX_REGEX,
    "KEGG": KEGG_REGEX,
    "CHEBI": CHEBI_REGEX,
    "SEED": SEED_REGEX,
}


# Create a function to determine database link information
def find_database_ids(
    identifier: str,
    bigg_metabolite_df: pd.DataFrame,
    db_regex_dict: dict[str, re.Pattern],
) -> Optional[dict[str, list[str]]]:
    """
    Extract database identifier information from the bigg_metabolite_df
    """
    # Get the database links string
    bigg_met_df_filtered = bigg_metabolite_df.loc[
        bigg_metabolite_df["bigg_id"] == identifier, ["database_links"]
    ]
    if len(bigg_met_df_filtered) == 0:
        return None  # Couldn't find id
    db_link_str = bigg_met_df_filtered.iloc[0, 0]
    # Find the various database identifiers of interest
    db_id_dict = {}
    for db_name, db_regex in db_regex_dict.items():
        db_id_dict[db_name] = db_regex.findall(db_link_str)
    return db_id_dict


# Add the database id columns to the metabolite_info_df
metabolite_info_df["bigg_id"] = ""
metabolite_info_df["MetaNetX id"] = ""
metabolite_info_df["KEGG id"] = ""
metabolite_info_df["SEED id"] = ""
metabolite_info_df["CHEBI id"] = ""

logger.info("Getting database links for each metabolite")

for metabolite in metabolite_info_df["id"]:
    # Get the bigg id
    if metabolite not in metabolite_to_bigg_id_dict:
        continue
    bigg_id = metabolite_to_bigg_id_dict[metabolite]
    metabolite_info_df.loc[
        metabolite_info_df["id"] == metabolite, "bigg_id"
    ] = bigg_id
    # Get database links for the id
    db_link_ids = find_database_ids(
        identifier=bigg_id,
        bigg_metabolite_df=BIGG_METABOLITE_DF,
        db_regex_dict=DATABASE_REGEX_DICT,
    )
    if db_link_ids is None:
        continue
    metabolite_info_df.loc[
        metabolite_info_df["id"] == metabolite, "MetaNetX id"
    ] = ";".join(db_link_ids["MetaNetX"])
    metabolite_info_df.loc[
        metabolite_info_df["id"] == metabolite, "KEGG id"
    ] = ";".join(db_link_ids["KEGG"])
    metabolite_info_df.loc[
        metabolite_info_df["id"] == metabolite, "SEED id"
    ] = ";".join(db_link_ids["SEED"])
    metabolite_info_df.loc[
        metabolite_info_df["id"] == metabolite, "CHEBI id"
    ] = ";".join(db_link_ids["CHEBI"])

KEGG_COMPOUND_TO_PATHWAY = pd.read_csv(
    "https://rest.kegg.jp/link/compound/pathway",
    sep="\t",
    names=["pathway", "compound"],
)
KEGG_COMPOUND_TO_PATHWAY["compound"] = KEGG_COMPOUND_TO_PATHWAY[
    "compound"
].str.replace("^cpd:", "", regex=True)
KEGG_COMPOUND_TO_PATHWAY["pathway"] = KEGG_COMPOUND_TO_PATHWAY[
    "pathway"
].str.replace("^path:", "", regex=True)
# Join this with the pathway descriptions
KEGG_COMPOUND_TO_PATHWAY = pd.merge(
    KEGG_COMPOUND_TO_PATHWAY,
    KEGG_PATHWAY_DESCRIPTIONS,
    how="left",
    left_on="pathway",
    right_on="pathway",
)


# Get KEGG pathway information
def find_kegg_pathways(kegg_id: str) -> Optional[str]:
    """
    Find KEGG pathways associated with a metabolite
    """
    compound_df = KEGG_COMPOUND_TO_PATHWAY[
        KEGG_COMPOUND_TO_PATHWAY["compound"] == kegg_id
    ]
    return ";".join(compound_df["description"])


# Get KEGG module information
kegg_compound_to_module = pd.read_csv(
    "https://rest.kegg.jp/link/compound/module",
    sep="\t",
    names=["module", "compound"],
)
kegg_compound_to_module["module"] = kegg_compound_to_module[
    "module"
].str.replace("^md:", "", regex=True)
kegg_compound_to_module["compound"] = kegg_compound_to_module[
    "compound"
].str.replace("^cpd:", "", regex=True)
kegg_module_descriptions = pd.read_csv(
    "https://rest.kegg.jp/list/module",
    sep="\t",
    names=["module", "description"],
)
kegg_compound_to_module = pd.merge(
    kegg_compound_to_module,
    kegg_module_descriptions,
    how="left",
    left_on="module",
    right_on="module",
)


def find_kegg_modules(kegg_id: str) -> Optional[str]:
    """
    Find KEGG modules associated with a metabolite
    """
    compound_df = kegg_compound_to_module[
        kegg_compound_to_module["compound"] == kegg_id
    ]
    return ";".join(compound_df["description"])


# Add the KEGG pathway information to the metabolite info dataframe
logger.info("Starting to add kegg information to metabolite dataframe")
metabolite_info_df["KEGG pathways"] = ""
metabolite_info_df["KEGG modules"] = ""
for metabolite in metabolite_info_df["id"]:
    logger.info(f"Getting KEGG pathways for {metabolite}")
    kegg_id = metabolite_info_df.loc[
        metabolite_info_df["id"] == metabolite, "KEGG id"
    ].iloc[0]  # type:ignore
    # If the string is not empty, find the associated pathways
    kegg_pathway_list: list[str] = []
    kegg_module_list: list[str] = []
    for id in kegg_id.split(";"):
        k_paths = find_kegg_pathways(id)
        k_modules = find_kegg_modules(id)
        if k_paths is not None:
            kegg_pathway_list.append(k_paths)
        if k_modules is not None:
            kegg_module_list.append(k_modules)
    metabolite_info_df.loc[
        metabolite_info_df["id"] == metabolite, "KEGG pathways"
    ] = ";".join(kegg_pathway_list)
    metabolite_info_df.loc[
        metabolite_info_df["id"] == metabolite, "KEGG modules"
    ] = ";".join(kegg_module_list)

# Save the updated metbolite information dataframe
logger.info("Saving metabolite information")
metabolite_info_df.to_csv(
    CACHE_PATH / "metabolite_information.csv",
    index=False,
)
