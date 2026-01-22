"""
Script to generate visualizations for the Mtb Transcription
Factor Analyses
"""

# Setup
# Imports
# Standard Library Imports
import json
import pathlib
import logging
import sys
import tomllib
import warnings

# External Imports
import cobra
import matplotlib.pyplot as plt
from metabolic_modeling_utils import escher_maps
import metworkpy
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

# Local Imports

# Path Setup# Path Setup
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".").absolute()  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
CACHE_PATH = BASE_PATH / "cache"
MODEL_PATH = BASE_PATH / "models"
RESULTS_PATH = BASE_PATH / "results" / "mtb_transcription_factors"
ESCHER_MAPS_PATH = BASE_PATH / "escher_maps" / "mtb"
ESCHER_MAPS_OUT_DIR = RESULTS_PATH / "escher_maps"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"

# Made directory if required
LOG_PATH.mkdir(exist_ok=True, parents=True)
ESCHER_MAPS_OUT_DIR.mkdir(exist_ok=True, parents=True)
# Setup Logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=LOG_PATH / "10_mtb_tf_visualization.log",
    filemode="w",
    level=logging.INFO,
)

# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)
cobra.Configuration().solver = CONFIG["cobra"]["solver"]

# Create a list of the escher maps
ESCHER_MAPS_INPUT_LIST = [
    ESCHER_MAPS_PATH / "mtb" / f"{map}.json"
    for map in [
        "Arabinogalactan_peptidoglycan_complex",
        "central_carbon",
        "sulfur_and_folate",
        "nitrogen",
    ]
]

# Create the sklearn scaler
SCALER = StandardScaler()
SCALER.set_output(transform="pandas")

# Read in the iEk1011 and iEK1011_v2 base models
model_v2 = metworkpy.read_model(
    MODEL_PATH / "iEK1011_v2_7H9_ADC_glycerol.json"
)

model_v1 = metworkpy.read_model(MODEL_PATH / "iEK1011_7H9_ADC_glycerol.json")

# Create a dictionary to map reaction and metabolite names
# from iEK1011_v1 to iEK1011_v2
# Create a rename dictionary for reactions from iEK1011_v2 to iEK1011
model_v1_rxn_ids = model_v1.reactions.list_attr("id")

rxn_rename_dict = {
    "CITL": "CITLr",
    "ICDHyr": "ICDHy",
    "PDHcr": "PDHc",
    "SUCD1": "SUCD1i",
    "CYRDAR2": "CYRDAR",
    "FOMETRi": "FOMETR",
    "HEMEAS": "HEMEAS_1",
    "ADMDC": "ADMDCr",
    "SHSL1": "SHSL1r",
    "OCOAT1": "OCOAT1r",
    "DHNPA2r": "DHNPA2",
    "GHMT2r": "GHMT2",
    "PRMICI": "PRMICIi",
    "3M2OPLOXRD": "BKDA1",
    "PPTT": "PPTTm",
    "PAPSR": "PAPR",
    "NAPRT": "NAPRTr",
    "NNATr": "NNAT",
    "PDX5POi": "PDX5PO",
    "G6PDH2r": "G6PDH",
    "G6PDH3": "G6PDH",
    "ACALD": "ACALDi",
    "L_LACDcm_1": "L_LACD",
    "ALDD19xr": "ALDD19x",
    "KARA1": "KARA1i",
    "KARA2": "KARA2i",
    "GLYDHDA_copy2": "GXRA",
    "PYDXS": "PDBL_1",
    "PLPS": "PDBL_2",
    "CPMPS": "MBP_1",
    "APPLDHr": "APRO",
    "G3PD5": "GLPD2",
    "G6PP": "PPGKr",
    "ALCD2ir": "ALCD2x",
    "MI3PS": "INO1",
    "ACOADH1": "OIVD3",
    "MICITDr": "2MID",
    "MCITL2": "2MIL",
    "DHORDi": "PYRD1",
    "MGSA": "SR2",
    "LCADi": "LALDD",
}
for rxn in model_v2.reactions:
    if rxn.subsystem in {
        "Extracellular exchange",
        "Biomass and maintenance functions",
        "Intracellular demand",
        "Transport",
    }:
        # These types of reaction aren't represented on the Escher map anyway
        continue
    if rxn.id in {
        "HIPEOX",
        "IPDCF",
        "ECHA20_2",
        "FADA6",
        "FADE31",
        "MBKTHIO",
        "SDH1",
        "ATPM",
        "BiomassGrowth_sMtb",
        "BIOMASS__2",
        "BIOMASS__2.1",
    }:
        # reactions added by v2 update
        continue
    if rxn.id in rxn_rename_dict:
        continue  # hand mapped already
    if rxn.id in {"GLYDHDA_copy1", "GF4GL_7", "IPDF"}:
        continue  # Unmappable
    try:
        # Try just using the same id
        new_id = model_v1.reactions.get_by_id(rxn.id).id
        rxn_rename_dict[rxn.id] = new_id
    except KeyError:
        try:
            # Try finding a reaction with the same name
            new_id = model_v1.reactions.get_by_any(rxn.name)[0].id
            rxn_rename_dict[rxn.id] = new_id
        except KeyError:
            # If this is a single gene reaction,
            # try mapping through that
            if len(rxn.genes) == 1:
                rxn_gene = rxn.genes.__iter__().__next__()
                try:
                    model_v1_gene = model_v1.genes.get_by_id(rxn_gene.id)
                    if len(model_v1_gene.reactions) == 1:
                        new_id = (
                            model_v1_gene.reactions.__iter__().__next__().id
                        )
                        rxn_rename_dict[rxn.id] = new_id
                        continue
                except KeyError:
                    warnings.warn(f"Couldn't find mapping for {rxn.id}")
            try:
                split_id = "".join(rxn.id.split("_")[:-1])
                new_id = model_v1.reactions.get_by_id(split_id).id
                rxn_rename_dict[rxn.id] = new_id
            except KeyError:
                for potential_new_id in model_v1_rxn_ids:
                    if (
                        rxn.id in potential_new_id
                        or potential_new_id in rxn.id
                    ):
                        warnings.warn(
                            f"Mapping {rxn.id} to {potential_new_id}"
                        )
                        rxn_rename_dict[rxn.id] = potential_new_id
                        break
                else:
                    warnings.warn(f"Couldn't find mapping for {rxn.id}")

# Create a dict for renaming metabolites as well
# Use the metabolite id, which has the BiGG IDs
met_rename_dict = (
    pd.read_csv(
        CACHE_PATH / "model_information" / "metabolite_information.csv",
    )
    .set_index("id")["bigg_id"]
    .to_dict()
)
met_rename_dict = {
    k: v for k, v in met_rename_dict.items() if not isinstance(v, float)
}


##################
### Divergence ###
##################
# Read in the divergence results
divergence_df = SCALER.fit_transform(
    pd.read_csv(RESULTS_PATH / "divergence_results.csv", index_col=0)
)
# Split reaction and metabolite divergence
rxn_div_df = divergence_df.loc[
    :, divergence_df.columns.str.startswith("reaction__")
]
rxn_div_df.columns = rxn_div_df.columns.str.replace(
    "^reaction__", "", regex=True
)

met_div_df = divergence_df.loc[
    :, divergence_df.columns.str.startswith("reaction__")
]
met_div_df.columns = met_div_df.columns.str.replace(
    "^reaction__", "", regex=True
)

# Extract the ArgR values
# Renaming them to match the iEK1011 identifiers
argr_rxn_divergence = rxn_div_df.loc["Rv1657"].rename(rxn_rename_dict)
argr_met_divergence = met_div_df.loc["Rv1657"].rename(met_rename_dict)

# Show the divergence for all of the escher maps
for input_map in ESCHER_MAPS_INPUT_LIST:
    escher_maps.escher_map_add_data(
        input_map=input_map,
        output_dir=ESCHER_MAPS_OUT_DIR,
        output_prefix="ArgR_divergence_",
        reaction_data=argr_rxn_divergence.replace(
            [np.inf, -np.inf], np.nan
        ).dropna(),
        metabolite_data=argr_met_divergence.replace(
            [np.inf, -np.inf], np.nan
        ).dropna(),
        reaction_scale=[
            {"type": "min", "color": "blue", "size": 10},
            {"type": "max", "color": "red", "size": 30},
        ],
        metabolite_scale=[
            {"type": "min", "color": "blue", "size": 10},
            {"type": "max", "color": "red", "size": 30},
        ],
    )

###############
### Density ###
###############
# Read in the density
density_df = pd.read_csv(RESULTS_PATH / "tf_target_density.csv", index_col=0)

argr_density = density_df["Rv1657"].rename(rxn_rename_dict)

for input_map in ESCHER_MAPS_INPUT_LIST:
    escher_maps.escher_map_add_data(
        input_map=input_map,
        output_dir=ESCHER_MAPS_OUT_DIR,
        output_prefix="ArgR_density_",
        reaction_data=argr_density.dropna(),
        reaction_scale=[
            {"type": "max", "color": "red", "size": 30},
            {"type": "value", "value": 0.0, "color": "grey", "size": 10},
        ],
    )
