"""
Script to combine the results from the various analyses into
a single file
"""

# Setup
# Imports
# Standard Library Imports
import json
import pathlib
import sys
import tomllib
from typing import cast

# External Imports
import cobra
import pandas as pd
import metworkpy
import networkx as nx
import numpy as np
from sklearn.preprocessing import StandardScaler

# Path Setup
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".").absolute()  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent
MODEL_PATH = BASE_PATH / "models"
RESULTS_PATH = BASE_PATH / "results"
SIMULATION_RESULTS_PATH = RESULTS_PATH / "simulation"
MTB_TF_RESULTS_PATH = RESULTS_PATH / "mtb_transcription_factors"
CACHE_PATH = BASE_PATH / "cache"

# Script Configuration
# Read in the configuration file
with open(BASE_PATH / "config.toml", "rb") as f:
    CONFIG = tomllib.load(f)
cobra.Configuration().solver = CONFIG["cobra"]["solver"]


##########################
### Simulation Results ###
##########################
# Use a function to name space the simulation result collation
def collate_simulation_results():
    # -----------------------
    # -- Model Information --
    # -----------------------
    simulation_model = metworkpy.read_model(
        MODEL_PATH / "simulation_model.json"
    )
    model_df = pd.DataFrame(
        "",
        index=pd.Index(simulation_model.reactions.list_attr("id")),
        columns=pd.Index(["Equation", "Gene-Reaction Rule", "Subsystem"]),
    )
    for rxn in simulation_model.reactions:
        rxn = cast(cobra.Reaction, rxn)
        model_df.loc[rxn.id, "Equation"] = rxn.build_reaction_string()
        model_df.loc[rxn.id, "Gene-Reaction Rule"] = rxn.gene_reaction_rule
        model_df.loc[rxn.id, "Subsystem"] = rxn.subsystem
    model_df = model_df.reset_index(drop=False, names="Reaction")

    # ----------------------
    # -- Metabolic Graphs --
    # ----------------------
    def to_pandas_adjacency(graph: nx.Graph) -> pd.DataFrame:
        """Convert a networkx graph to a pandas dataframe"""
        idx = pd.Index(graph.nodes)
        return pd.DataFrame(nx.to_numpy_array(graph), index=idx, columns=idx)

    def read_graph(path: pathlib.Path) -> nx.Graph:
        """
        Read a networkx graph from a json file
        """
        with open(path, "r") as f:
            graph = nx.node_link_graph(
                json.load(f), multigraph=False, directed=False, edges="edges"
            )
        return graph

    # Read in the metabolic, reaction, and metabolite graphs
    metabolic_connectivity_adjacency = to_pandas_adjacency(
        read_graph(
            SIMULATION_RESULTS_PATH
            / "metabolic_networks"
            / "metabolic_network.json"
        )
    )
    reaction_connectivity_adjacency = to_pandas_adjacency(
        read_graph(
            SIMULATION_RESULTS_PATH
            / "metabolic_networks"
            / "metabolic_reaction_network.json"
        )
    )
    metabolite_connectivity_adjacency = to_pandas_adjacency(
        read_graph(
            SIMULATION_RESULTS_PATH
            / "metabolic_networks"
            / "metabolic_metabolite_network.json"
        )
    )

    metabolite_centrality = pd.read_csv(
        SIMULATION_RESULTS_PATH
        / "metabolic_networks"
        / "metabolite_centrality.csv",
        index_col=0,
    ).reset_index(drop=False, names="Metabolite")

    reaction_centrality = pd.read_csv(
        SIMULATION_RESULTS_PATH
        / "metabolic_networks"
        / "reaction_centrality.csv",
        index_col=0,
    ).reset_index(drop=False, names="Reaction")

    # --------------------
    # -- Flux MI Graphs --
    # --------------------
    mi_adjacency = pd.read_csv(
        SIMULATION_RESULTS_PATH / "mutual_information" / "mi_adjacency.csv",
        index_col=0,
    )
    mi_centrality = pd.read_csv(
        SIMULATION_RESULTS_PATH / "mutual_information" / "mi_centrality.csv",
        index_col=0,
    ).reset_index(drop=False, names="Reaction")
    # ----------------------------
    # -- Metabolite Subnetworks --
    # ----------------------------
    metabolite_consuming_network = pd.read_csv(
        SIMULATION_RESULTS_PATH
        / "metabolite_networks"
        / "metabolite_consuming_network.csv",
        index_col=0,
    )
    metabolite_synthesis_network = pd.read_csv(
        SIMULATION_RESULTS_PATH
        / "metabolite_networks"
        / "metabolite_synthesis_network.csv",
        index_col=0,
    )
    # ----------------------------------------------
    # -- Reaction Neighborhood Enrichment/Density --
    # ----------------------------------------------
    gene_target_density = pd.read_csv(
        SIMULATION_RESULTS_PATH / "target_density" / "gene_target_density.csv",
        index_col=0,
    ).reset_index(names="Reaction")
    gene_target_enrichment = pd.read_csv(
        SIMULATION_RESULTS_PATH
        / "target_density"
        / "gene_target_enrichment_pval.csv",
        index_col=0,
    ).reset_index(names="Reaction")

    # -------------------
    # -- KO Divergence --
    # -------------------
    ko_divergence = pd.read_csv(
        SIMULATION_RESULTS_PATH
        / "ko_divergence"
        / "ko_divergence_results.csv",
        index_col=0,
    ).reset_index(names="Gene")
    ko_divergence = ko_divergence.melt(
        id_vars="Gene", var_name="Divergence Group", value_name="Divergence"
    )
    ko_divergence["Divergence Group Type"] = (
        ko_divergence["Divergence Group"].str.split("__").str[0]
    )
    ko_divergence["Divergence Group"] = (
        ko_divergence["Divergence Group"]
        .str.split("__")
        .str[1:]
        .str.join("__")
    )
    ko_divergence = ko_divergence[
        ["Gene", "Divergence Group Type", "Divergence Group", "Divergence"]
    ]

    # ---------------------
    # -- IMAT Divergence --
    # ---------------------
    imat_div = pd.read_csv(
        SIMULATION_RESULTS_PATH / "iMAT" / "imat_divergence.csv", index_col=0
    ).reset_index(names="Divergence Group")
    imat_div["Divergence Group Type"] = (
        imat_div["Divergence Group"].str.split("__").str[0]
    )
    imat_div["Divergence Group"] = (
        imat_div["Divergence Group"].str.split("__").str[1:].str.join("__")
    )
    imat_div = imat_div[
        ["Divergence Group Type", "Divergence Group", "IMAT Divergence"]
    ]

    # ---------------------------
    # -- Save Combined Results --
    # ---------------------------
    with pd.ExcelWriter(
        RESULTS_PATH / "simulation_model_results.xlsx"
    ) as writer:
        # Model Information
        model_df.to_excel(writer, sheet_name="Model Information", index=False)
        # Metabolic Connectivity Networks
        metabolic_connectivity_adjacency.to_excel(
            writer, sheet_name="Metabolic SCN"
        )
        reaction_connectivity_adjacency.to_excel(
            writer, sheet_name="Reaction SCN", index=True
        )
        metabolite_connectivity_adjacency.to_excel(
            writer, sheet_name="Metabolite SCN", index=True
        )
        metabolite_centrality.to_excel(
            writer, sheet_name="SCN Metabolite Centrality", index=False
        )
        reaction_centrality.to_excel(
            writer, sheet_name="SCN Reaction Centrality", index=False
        )
        # Mutual information networks
        mi_adjacency.to_excel(
            writer, sheet_name="Flux MI Adjacency Matrix", index=True
        )
        mi_centrality.to_excel(
            writer, sheet_name="Flux MI Centrality", index=False
        )
        # Metabolite Networks
        metabolite_synthesis_network.to_excel(
            writer, sheet_name="Metabolite synthesis Networks", index=True
        )
        metabolite_consuming_network.to_excel(
            writer, sheet_name="Metabolite consuming Networks", index=True
        )
        # Reaction neighborhood enrichment/density
        gene_target_density.to_excel(
            writer, sheet_name="Gene Target Density", index=False
        )
        gene_target_enrichment.to_excel(
            writer, sheet_name="Gene Target Enrichment", index=False
        )
        # KO Divergence
        ko_divergence.to_excel(writer, sheet_name="KO Divergence", index=False)
        # IMAT Divergence
        imat_div.to_excel(writer, sheet_name="IMAT Divergence", index=False)


######################
### Mtb TF Results ###
######################
def collate_mtb_tf_results():
    """
    Combine the results of the analysis of the Mtb Transcription Factors
    """
    # -----------------------
    # -- Model Information --
    # -----------------------
    reaction_info = pd.read_csv(
        CACHE_PATH / "model_information" / "reaction_information.csv"
    )
    metabolite_info = pd.read_csv(
        CACHE_PATH / "model_information" / "metabolite_information.csv"
    )

    # ----------------------
    # -- Metabolic Graphs --
    # ----------------------
    rxn_centrality = pd.read_csv(
        MTB_TF_RESULTS_PATH / "metabolic_reaction_network_centrality.csv",
        index_col=0,
    ).reset_index(drop=False, names="Reaction")
    rxn_centrality_analysis = (
        pd.read_csv(
            MTB_TF_RESULTS_PATH
            / "metabolic_reaction_network_centrality_analysis.csv",
            index_col=0,
        )
        .reset_index(names="Transcription Factor")
        .loc[:, "Transcription Factor":"betweenness bootstrap adj p-value"]
    ).dropna(how="any", axis="index")

    # ------------------------
    # -- Mutual Information --
    # ------------------------
    flux_mi_gene_centrality = pd.read_csv(
        MTB_TF_RESULTS_PATH / "flux_mi_gene_centrality.csv", index_col=0
    ).reset_index(names="Gene")
    flux_mi_ess_vi_stats = pd.read_csv(
        MTB_TF_RESULTS_PATH
        / "mutual_information_vi_essentiality_statistics.csv",
        index_col=0,
    )

    # ----------------------------
    # -- Metabolite Subnetworks --
    # ----------------------------
    tf_metabolite_network_enrichment = pd.read_csv(
        MTB_TF_RESULTS_PATH / "tf_target_metabolite_network_enrichment.csv",
    ).loc[:, "metabolite":"tf"][
        [
            "tf",
            "metabolite",
            "metabolite network direction",
            "metabolite network size",
            "tf target count",
            "tf target-metabolite network overlap",
            "total genes",
            "odds-ratio",
            "p-value",
            "adj p-value",
        ]
    ]
    tf_subsystem_network_enrichment = pd.read_csv(
        MTB_TF_RESULTS_PATH / "tf_target_subsystem_network_enrichment.csv"
    )[
        [
            "tf",
            "subsystem",
            "subsystem size",
            "tf target count",
            "tf target-subsystem network overlap",
            "total genes",
            "odds-ratio",
            "p-value",
            "adj p-value",
        ]
    ]
    tf_metabolite_gsva = (
        pd.read_csv(MTB_TF_RESULTS_PATH / "metabolite_gsva.csv", index_col=0)
        .reset_index(names="TF")
        .melt(id_vars="TF", var_name="Metabolite", value_name="GSVA")
    )
    tf_metabolite_gsva["Network Direction"] = (
        tf_metabolite_gsva["Metabolite"].str.split("_").str[-2]
    )
    tf_metabolite_gsva["Metabolite"] = (
        tf_metabolite_gsva["Metabolite"].str.split("_").str[:-2].str.join("_")
    )

    # ----------------------------------------------
    # -- Reaction Neighborhood Enrichment/Density --
    # ----------------------------------------------
    tf_target_density = (
        pd.read_csv(MTB_TF_RESULTS_PATH / "tf_target_density.csv")
        .set_index("id")["Rv1657"]
        .to_frame(name="ArgR Neighborhood Target Density")
        .reset_index(names="Reaction ID")
    )
    tf_rxn_neghborhood_enrichment = (
        pd.read_csv(MTB_TF_RESULTS_PATH / "tf_target_enrichment.csv")
        .set_index("id")["Rv1657"]
        .to_frame(name="ArgR Neighborhood Target Enrichment")
        .reset_index(names="Reaction ID")
    )

    # -------------------
    # -- KO Divergence --
    # -------------------
    ko_divergence_df = pd.read_csv(
        (
            CACHE_PATH
            / "gene_ko_divergence"
            / "7h9_adc"
            / "gene_ko_divergence_results.csv"
        ),
        index_col=0,
    ).clip(lower=0.0)
    ko_divergence_df = (
        ko_divergence_df.loc[
            :,
            (
                (ko_divergence_df.columns.str.endswith("__metabolite"))
                | (ko_divergence_df.columns == "BIOMASS__2__reaction")
            ),
        ]
        .reset_index(names="Gene")
        .melt(
            id_vars="Gene",
            var_name="Divergence Group",
            value_name="Divergence",
        )
    )
    ko_divergence_df["Divergence Group Type"] = (
        ko_divergence_df["Divergence Group"].str.split("__").str[-1]
    )
    ko_divergence_df["Divergence Group"] = (
        ko_divergence_df["Divergence Group"]
        .str.replace("__reaction", "")
        .str.replace("__metabolite", "")
    )

    tf_target_ko_divergence = pd.read_csv(
        MTB_TF_RESULTS_PATH / "ko_divergence_tf_target_analysis.csv"
    ).rename({"rho": "AUC-ROC"}, axis=1)[
        [
            "tf",
            "metabolite",
            "represented metabolites",
            "Mann-Whitney U1",
            "Mann-Whitney U2",
            "AUC-ROC",
            "p-value",
            "adj p-value",
        ]
    ]
    tf_target_ko_divergence["represented metabolites"] = (
        tf_target_ko_divergence["represented metabolites"].str.replace(
            "set()", ""
        )
    )

    # Add the represented metabolites column to the ko_divergence_df
    represented_metabolites = (
        tf_target_ko_divergence.set_index("metabolite")[
            "represented metabolites"
        ]
        .to_frame("Represented Metabolites")
        .reset_index(names="Metabolite")
        .drop_duplicates()
    )
    ko_divergence_df = ko_divergence_df.merge(
        represented_metabolites,
        how="left",
        left_on="Divergence Group",
        right_on="Metabolite",
    )[
        [
            "Gene",
            "Divergence Group",
            "Divergence Group Type",
            "Represented Metabolites",
            "Divergence",
        ]
    ]
    # ---------------------
    # -- IMAT Divergence --
    # ---------------------
    scaler = StandardScaler()
    scaler.set_output(transform="pandas")
    imat_divergence = scaler.fit_transform(
        pd.read_csv(
            MTB_TF_RESULTS_PATH / "divergence_results.csv", index_col=0
        ).replace([np.inf, -np.inf], np.nan)
    ).dropna(axis="columns", how="all")
    imat_divergence.index.name = "TF"
    imat_reaction_divergence = imat_divergence.loc[
        :, imat_divergence.columns.str.startswith("reaction__")
    ]
    imat_reaction_divergence.columns = (
        imat_reaction_divergence.columns.str.replace("reaction__", "")
    )

    imat_metabolite_synthesis_divergence = imat_divergence.loc[
        :, imat_divergence.columns.str.startswith("metabolite_synthesis__")
    ]
    imat_metabolite_synthesis_divergence.columns = (
        imat_metabolite_synthesis_divergence.columns.str.replace(
            "metabolite_synthesis__", ""
        )
    )
    # ------------------
    # -- IMAT Compare --
    # ------------------
    imat_compare_df = (
        pd.read_csv(MTB_TF_RESULTS_PATH / "imat_compare.csv")
        .set_index("id")
        .reset_index(names="Reaction")
    ).rename(
        {
            "pFBA fluxes": "pFBA fluxes",
            "IMAT fluxes": "IMAT solution fluxes",
            "FVA IMAT pFBA fluxes": "FVA IMAT Model pFBA fluxes",
            "diff imat": "IMAT solution Fluxes - pFBA fluxes",
            "diff fva imat": "FVA IMAT Model pFBA fluxes - pFBA fluxes",
        },
        axis=1,
    )[
        [
            "Reaction",
            "pFBA fluxes",
            "IMAT solution fluxes",
            "FVA IMAT Model pFBA fluxes",
            "IMAT solution Fluxes - pFBA fluxes",
            "FVA IMAT Model pFBA fluxes - pFBA fluxes",
        ]
    ]
    # Add in the subsystem columns
    imat_compare_df = imat_compare_df.set_index("Reaction")
    imat_compare_df["subsystem"] = reaction_info.set_index("id")["subsystem"]

    with pd.ExcelWriter(RESULTS_PATH / "mtb_tf_results.xlsx") as writer:
        # Model Information
        reaction_info.to_excel(
            writer, sheet_name="Reaction Information", index=False
        )
        metabolite_info.to_excel(
            writer, sheet_name="Metabolite Information", index=False
        )
        # Metabolic Graphs
        rxn_centrality.to_excel(
            writer,
            sheet_name="Reaction SCN Centrality",
            index=False,
        )
        rxn_centrality_analysis.to_excel(
            writer, sheet_name="TF Target Centrality", index=False
        )
        # Flux MI Centrality
        flux_mi_gene_centrality.to_excel(
            writer, sheet_name="Flux MI Network Centrality", index=False
        )
        flux_mi_ess_vi_stats.to_excel(
            writer, sheet_name="MI Centrality vs Essentialiy", index=True
        )
        # Metabolite Enrichment
        tf_metabolite_network_enrichment.to_excel(
            writer, sheet_name="TF Metabolite Enrichment", index=False
        )
        tf_subsystem_network_enrichment.to_excel(
            writer,
            sheet_name="TF Subsystem Enrichment",
            index=False,
        )
        tf_metabolite_gsva.to_excel(
            writer, sheet_name="TF Metabolite GSVA", index=False
        )
        # Target density/enrichment
        tf_target_density.to_excel(
            writer, sheet_name="ArgR Rxn Neighborhood Density", index=False
        )
        tf_rxn_neghborhood_enrichment.to_excel(
            writer, sheet_name="ArgR Rxn Neighbor Enrichment", index=False
        )
        # KO Divergence
        ko_divergence_df.to_excel(
            writer, sheet_name="Gene KO Divergence", index=False
        )
        tf_target_ko_divergence.to_excel(
            writer, sheet_name="TF Target KO Divergence", index=False
        )
        # IMAT Divergence
        imat_reaction_divergence.to_excel(
            writer, sheet_name="Normalized IMAT Reaction Div", index=True
        )
        imat_metabolite_synthesis_divergence.to_excel(
            writer, sheet_name="Normalized IMAT Metabolite Div", index=True
        )
        # IMAT Compare
        imat_compare_df.to_excel(
            writer, sheet_name="ArgR IMAT Fluxes", index=True
        )


if __name__ == "__main__":
    collate_simulation_results()
    collate_mtb_tf_results()
