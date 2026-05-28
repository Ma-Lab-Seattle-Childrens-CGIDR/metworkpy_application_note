"""
Script for evaluating the impacts of parameters of network generation
"""

# Setup
# Imports
# Standard Library Imports
import logging
import pathlib
import sys
import tomllib

# External Imports
import cobra
import metworkpy
import networkx as nx
import pandas as pd

# Local Imports

# Setup Path
if hasattr(sys, "ps1"):
    # Running in a REPL
    BASE_PATH = pathlib.Path(".").absolute()  # Use current dir as base path
else:
    # Running as a file
    # Use file path to find root
    BASE_PATH = pathlib.Path(__file__).parent.parent.parent
CACHE_PATH = BASE_PATH / "cache"
DATA_PATH = BASE_PATH / "data"
MODEL_PATH = BASE_PATH / "models"
LOG_PATH = BASE_PATH / "logs" / "mtb_transcription_factors"
RESULTS_PATH = BASE_PATH / "results" / "mtb_transcription_factors"

# Create directories if needed
LOG_PATH.mkdir(parents=True, exist_ok=True)
RESULTS_PATH.mkdir(parents=True, exist_ok=True)


def evaluate_network(
    model: cobra.Model,
    reactions_to_remove: list[str],
    metabolite_counts: pd.Series,
    n: int = 10,
) -> pd.Series:
    """
    Constructs a network after removing the top `n` most connected metabolites,
    then evaluates various properties of the network
    """
    # Construct the network after removing the metabolites
    logger.info("Constructing the network")
    rxn_network = metworkpy.network.create_reaction_network(
        model=model,
        weighted=False,
        directed=False,
        nodes_to_remove=reactions_to_remove
        + list(metabolite_counts.sort_values(ascending=False).iloc[:n]),
    )
    # Create a dict to hold the results
    results_dict = {}

    # Calculate the number of components
    logger.info("Finding the connected components")
    results_dict["component count"] = nx.number_connected_components(
        rxn_network
    )
    # Get only the largest component for the remaining evaluations
    rxn_network = rxn_network.subgraph(
        max(nx.connected_components(rxn_network), key=len)
    )
    results_dict["largest component size"] = len(rxn_network.nodes)

    # Find the algebraic connectivity of the graph
    logger.info("Finding algebraic connectivity")
    results_dict["algebraic connectivity"] = nx.algebraic_connectivity(
        rxn_network
    )

    # Find the small-world coefficients of a graph
    logger.info("Finding small-world coefficients")
    results_dict["small-world sigma"] = nx.sigma(rxn_network)
    results_dict["small-world omega"] = nx.omega(rxn_network)

    # Estimate the average clustering of the graph
    logger.info("Finding the average clustering")
    results_dict["average clustering estimate"] = (
        nx.approximation.average_clustering(rxn_network)
    )

    # For several centrality measures, find the min/max/avg/stdev
    for measure, centrality_function in zip(
        ["closeness", "betweenness", "degree"],
        [
            nx.closeness_centrality,
            nx.betweenness_centrality,
            nx.degree_centrality,
        ],
    ):
        logger.info(f"Finding the {measure} centrality")
        cent_series = pd.Series(centrality_function(rxn_network))
        results_dict[f"{measure} centrality max"] = cent_series.max()
        results_dict[f"{measure} centrality min"] = cent_series.min()
        results_dict[f"{measure} centrality mean"] = cent_series.mean()
        results_dict[f"{measure} centrality std"] = cent_series.std()

    # Return the results dict as a series
    logger.info("Returning results")
    res_series = pd.Series(results_dict)
    res_series.name = str(n)
    return res_series


if __name__ == "__main__":
    # Setup Logging
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        filename=LOG_PATH / "12_network_evaluation.log",
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Read in the configuration file
    with open(BASE_PATH / "config.toml", "rb") as f:
        CONFIG = tomllib.load(f)

    # Change the default cobra solver
    cobra.Configuration().solver = CONFIG["cobra"]["solver"]

    # Read in base model
    logger.info("Reading in m7H9 model")
    BASE_MODEL = metworkpy.read_model(
        MODEL_PATH / "iEK1011_v2_7H9_ADC_glycerol.json"
    )

    logger.info("Finding reactions to remove")
    subsystems_to_ignore = {
        "Intracellular demand",
        "Biomass and maintenance functions",
        "Extracellular exchange",
    }
    reactions_to_remove = set()
    for rxn in BASE_MODEL.reactions:
        if rxn.subsystem in subsystems_to_ignore:
            reactions_to_remove.add(rxn.id)

    # Get a count of the number of reactions each metabolite is involved in
    logger.info("Getting counts of metabolite connectivity")
    stoich_mat = cobra.util.create_stoichiometric_matrix(
        model=BASE_MODEL, array_type="DataFrame"
    )
    assert isinstance(stoich_mat, pd.DataFrame), (
        "COBRA returned wrong stoichiometric matrix type"
    )
    metabolite_counts = (stoich_mat.abs() > 0.0).sum(axis=1)

    res_list = []
    for n in [0, 1, 5, 10, 15, 20, 30, 40, 50]:
        logger.info(f"Working on n={n}")
        res_list.append(
            evaluate_network(
                model=BASE_MODEL,
                reactions_to_remove=reactions_to_remove,
                metabolite_counts=metabolite_counts,
                n=n,
            )
        )
    logger.info("Finished evaluating the networks, saving results")
    res_df = pd.concat(res_list, axis=1).T
    res_df.to_csv(RESULTS_PATH / "network_metabolite_removal_evaluation.csv")
