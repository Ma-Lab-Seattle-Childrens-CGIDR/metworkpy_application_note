"""
Script to analyze enrichment of metabolite targets in
different metabolite networks
"""

# Setup
# Imports
# Standard Library Imports
import logging
import pathlib
import sys
from collections import defaultdict
from typing import Iterable, Mapping, cast

# External Imports
import cobra
import metworkpy
import numpy as np
import pandas as pd
import tomllib

# Local Imports
from common_functions import (
    fdr_with_nan,
    get_kegg_pathway_descriptions,
    get_kegg_pathway_genes,
    get_metabolite_network,
)
from scipy import stats

# Path setup
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = (
        pathlib.Path(".").absolute().parent.parent
    )  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
DATA_PATH = BASE_PATH / "data"
CACHE_PATH = BASE_PATH / "cache"
METABOLITE_NETWORKS_PATH = CACHE_PATH / "metabolite_networks" / "7H9_ADC"
RESULTS_PATH = BASE_PATH / "results" / "mtb_transcription_factors"
MODELS_PATH = BASE_PATH / "models"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"


def find_overlap_enrichment(
    target_dict: Mapping[str, Iterable[str]],
    gene_set_dict: Mapping[str, Iterable[str]],
    population_genes: Iterable[str],
) -> pd.DataFrame:
    # Convert the dicts to use sets for values to make the overlaps/unions easier
    target_dict: dict[str, set[str]] = {
        k: set(v) for k, v in target_dict.items()
    }
    gene_set_dict: dict[str, set[str]] = {
        k: set(v) for k, v in gene_set_dict.items()
    }
    population_genes = set(population_genes)

    enrich_results_list: list[pd.DataFrame] = []
    for regulator, target_set in target_dict.items():
        # Filter target set to be within population genes
        target_set = target_set & population_genes
        reg_enrich_df = pd.DataFrame(
            0.0,
            index=pd.Index(gene_set_dict.keys()),
            columns=pd.Index(
                [
                    "Gene Set Size",
                    "Regulator Target Set Size",
                    "Overlap",
                    "Total Genes",
                    "odds-ratio",
                    "p-value",
                ]
            ),
        )
        reg_enrich_df["Regulator Name"] = regulator
        reg_enrich_df["Total Genes"] = len(population_genes)
        for gene_set_name, gene_set in gene_set_dict.items():
            # Filter gene_set to be within population genes
            gene_set = gene_set & population_genes
            # Setup the contigency table
            # +---------------------------------------------------+
            # |                | In Target Set | Not In Target Set|
            # +================+===============+==================+
            # | In Pathway     |               |                  |
            # +----------------+---------------+------------------+
            # | Not in Pathway |               |                  |
            # +----------------+---------------+------------------+
            cont_table = np.zeros((2, 2))
            cont_table[0, 0] = len(target_set & gene_set)
            cont_table[0, 1] = len(gene_set - target_set)
            cont_table[1, 0] = len(target_set - gene_set)
            cont_table[1, 1] = len(population_genes) - len(
                gene_set | target_set
            )
            fisher_res = stats.fisher_exact(cont_table, alternative="greater")
            reg_enrich_df.loc[
                gene_set_name,
                [
                    "Gene Set Size",
                    "Regulator Target Set Size",
                    "Overlap",
                    "odds-ratio",
                    "p-value",
                ],
            ] = (
                len(gene_set),
                len(target_set),
                len(gene_set & target_set),
                fisher_res.statistic,
                fisher_res.pvalue,
            )

        # Adjust the p-values
        reg_enrich_df["Adjusted p-value"] = fdr_with_nan(
            reg_enrich_df["p-value"]
        )
        # Set the Gene Set names to be a new column
        reg_enrich_df = reg_enrich_df.reset_index(
            drop=False, names="Gene Set Name"
        )[
            [
                "Regulator Name",
                "Gene Set Name",
                "Regulator Target Set Size",
                "Gene Set Size",
                "Overlap",
                "Total Genes",
                "odds-ratio",
                "p-value",
                "Adjusted p-value",
            ]
        ]
        enrich_results_list.append(reg_enrich_df)
    return pd.concat(enrich_results_list, axis=0)


if __name__ == "__main__":
    # Make required directories if they don't exist
    CACHE_PATH.mkdir(parents=True, exist_ok=True)
    METABOLITE_NETWORKS_PATH.mkdir(parents=True, exist_ok=True)
    RESULTS_PATH.mkdir(parents=True, exist_ok=True)
    LOG_PATH.mkdir(parents=True, exist_ok=True)

    # Setup logging
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        filename=LOG_PATH / "08_metabolite_network_target_enrichment.log",
        filemode="w",
        level=logging.INFO,
    )

    # Read in the configuration file
    with open(BASE_PATH / "config.toml", "rb") as f:
        CONFIG = tomllib.load(f)

    # Read in the base model for finding metabolite networks
    cobra.Configuration().solver = CONFIG["cobra"]["solver"]
    BASE_MODEL = metworkpy.read_model(
        MODELS_PATH / "iEK1011_v2_7H9_ADC_glycerol.json"
    )

    # Determine the metabolite reaction synthesis/consumption networks
    # NOTE: The network dataframes are indexed by reactions, with columns for each metabolite
    # Create a list of reactions to remove from the metabolite networks
    reactions_to_remove = []
    subsystems_to_ignore = [
        "Intracellular demand",
        "Biomass and maintenance functions",
        "Extracellular exchange",
    ]
    for rxn in BASE_MODEL.reactions:
        if rxn.subsystem in subsystems_to_ignore:
            reactions_to_remove.append(rxn.id)
    metabolite_synthesis_rxn_net_df_path = (
        METABOLITE_NETWORKS_PATH / "metabolite_synthesis_reaction_network.csv"
    )
    metabolite_synthesis_rxn_network_df = get_metabolite_network(
        metabolite_synthesis_rxn_net_df_path,
        model=BASE_MODEL,
        rxns_to_remove=reactions_to_remove,
        proportion=CONFIG["mtb_tf"]["target_enrichment"][
            "essential-proportion"
        ],
        synthesis=True,
        processes=CONFIG["processes"],
    )

    # Repeat for the consuming networks
    metabolite_consumption_rxn_network_df_path = (
        METABOLITE_NETWORKS_PATH
        / "metabolite_consumption_reaction_network.csv"
    )
    metabolite_consumption_rxn_network_df = get_metabolite_network(
        metabolite_consumption_rxn_network_df_path,
        model=BASE_MODEL,
        rxns_to_remove=reactions_to_remove,
        proportion=CONFIG["mtb_tf"]["target_enrichment"][
            "reaction-proportion"
        ],
        synthesis=False,
        processes=CONFIG["processes"],
    )

    # Determine the TF targets

    # Find the TF targets
    logger.info("Finding the TF targets")
    tf_fc_df = pd.read_excel(
        DATA_PATH / "mtb_transcription_factors" / "tfoe_targets.xlsx",
        sheet_name="SupplementaryTableS2",
        skiprows=list(range(8)) + [9],
        usecols="A,E:HB",
        index_col=0,
    )
    tf_pval_df = pd.read_excel(
        DATA_PATH / "mtb_transcription_factors" / "tfoe_targets.xlsx",
        sheet_name="SupplementaryTableS2",
        skiprows=list(range(8)) + [9],
        usecols="A,HC:OZ",
        index_col=0,
    )
    tf_pval_df.columns = tf_pval_df.columns.str.replace(".1", "")

    tf_target_df = (
        tf_fc_df.abs()
        >= CONFIG["mtb_tf"]["target_enrichment"]["target-fc-cutoff"]
    ) & (
        tf_pval_df
        <= CONFIG["mtb_tf"]["target_enrichment"]["target-pval-cutoff"]
    )

    # Get a set of genes in the model
    model_gene_set = set(BASE_MODEL.genes.list_attr("id"))

    # Convert the TF target sets and Metabolite Networks into
    # Dictionaries

    # Metabolite Networks
    metabolite_synthesis_gene_network_dict: dict[str, set[str]] = {}
    for metabolite, rxn_series in metabolite_synthesis_rxn_network_df.items():
        met_net_gene_set = set(
            metworkpy.reaction_to_gene_list(
                model=BASE_MODEL,
                reaction_list=list(rxn_series[rxn_series].index),
                essential=False,
            )
        )
        if len(met_net_gene_set) <= 3:
            continue
        metabolite = cast(str, metabolite)
        metabolite_synthesis_gene_network_dict[metabolite] = met_net_gene_set
    metabolite_consuming_gene_network_dict: dict[str, set[str]] = {}
    for (
        metabolite,
        rxn_series,
    ) in metabolite_consumption_rxn_network_df.items():
        met_net_gene_set = set(
            metworkpy.reaction_to_gene_list(
                model=BASE_MODEL,
                reaction_list=list(rxn_series[rxn_series].index),
                essential=False,
            )
        )
        if len(met_net_gene_set) <= 3:
            continue
        metabolite = cast(str, metabolite)
        metabolite_consuming_gene_network_dict[metabolite] = met_net_gene_set

    # TF Target Sets
    tf_target_dict: dict[str, set[str]] = {}
    for tf, tf_target_series in tf_target_df.items():
        tf = cast(str, tf)
        tf_target_dict[tf] = set(tf_target_series[tf_target_series].index)

    # Now for synthesis and consumption networks, identify significant overlaps
    metabolite_enrichment_list = []
    for metabolite_network_direction, metabolite_gene_network_dict in zip(
        ["synthesis", "consuming"],
        [
            metabolite_synthesis_gene_network_dict,
            metabolite_consuming_gene_network_dict,
        ],
    ):
        enrichment_df = find_overlap_enrichment(
            tf_target_dict,
            gene_set_dict=metabolite_gene_network_dict,
            population_genes=model_gene_set,
        )
        enrichment_df["Metabolite Network Direction"] = (
            metabolite_network_direction
        )
    metabolite_network_enrichment = pd.concat(
        metabolite_enrichment_list, axis=0
    )[
        [
            "Regulator Name",
            "Gene Set Name",
            "Metabolite Network Direction",
            "Total Genes",
            "Gene Set Size",
            "Regulator Target Set Size",
            "Overlap",
            "Total Genes",
            "odds-ratio",
            "p-value",
            "Adjusted p-value",
        ]
    ]

    # Save the results
    metabolite_network_enrichment.to_csv(
        RESULTS_PATH / "tf_target_metabolite_network_enrichment.csv",
        index=False,
    )

    # Repeat for the Subsystems
    # Get a dict of subsystem to sets of genes
    subsystem_to_gene_dict: dict[str, set[str]] = defaultdict(set)
    for rxn in BASE_MODEL.reactions:
        if rxn.subsystem in subsystems_to_ignore:
            continue
        for gene in rxn.genes:
            subsystem_to_gene_dict[rxn.subsystem].add(gene.id)

    subsystem_enrichment_df = find_overlap_enrichment(
        tf_target_dict,
        gene_set_dict=subsystem_to_gene_dict,
        population_genes=model_gene_set,
    )[
        [
            "Regulator Name",
            "Gene Set Name",
            "Total Genes",
            "Gene Set Size",
            "Regulator Target Set Size",
            "Overlap",
            "Total Genes",
            "odds-ratio",
            "p-value",
            "Adjusted p-value",
        ]
    ]

    subsystem_enrichment_df.to_csv(
        RESULTS_PATH / "tf_target_subsystem_network_enrichment.csv",
        index=False,
    )

    # Finally for the KEGG Pathways
    kegg_path_df = get_kegg_pathway_genes("mtu")
    kegg_desc_df = get_kegg_pathway_descriptions(
        "mtu", remove_str=" - Mycobacterium tuberculosis H37Rv"
    )

    # Convert the kegg pathway genes into a dict, and get a set of all genes in KEGG pathways
    kegg_path_dict: dict[str, set[str]] = cast(
        dict[str, set[str]],
        {
            path: set(df["gene"])
            for path, df in kegg_path_df.groupby("pathway")
        },
    )
    kegg_gene_set = set(kegg_path_df["gene"])
    kegg_enrichment_df = (
        find_overlap_enrichment(
            tf_target_dict,
            gene_set_dict=kegg_path_dict,
            population_genes=kegg_gene_set,
        )
        .merge(
            kegg_desc_df,
            how="left",
            left_on="Gene Set Name",
            right_on="pathway",
        )[
            [
                "Regulator Name",
                "Gene Set Name",
                "description",
                "Total Genes",
                "Gene Set Size",
                "Regulator Target Set Size",
                "Overlap",
                "Total Genes",
                "odds-ratio",
                "p-value",
                "Adjusted p-value",
            ]
        ]
        .rename({"description": "KEGG Pathway Description"}, axis=1)
    )

    # Save results
    kegg_enrichment_df.to_csv(
        RESULTS_PATH / "tf_target_kegg_enrichment.csv", index=False
    )
