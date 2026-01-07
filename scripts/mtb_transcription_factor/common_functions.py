"""
Some common functionality used across multiple scripts
"""

import pathlib

import cobra  # type: ignore
import metworkpy
import pandas as pd


def get_metabolite_network(
    out_path: pathlib.Path,
    model: cobra.Model,
    rxns_to_remove: list[str] | None = None,
    proportion: float = 0.10,
    synthesis: bool = True,
    processes: int = -1,
) -> pd.DataFrame:
    """
    Read in or generate a metabolite network

    Parameters
    ----------
    out_path : pathlib.Path
        Path to read from or Save generated network dataframe to
    model : cobra.Model
        Model to use for generating the metabolite networks if needed
    rxns_to_remove : list of str
        Reactions to remove from the metabolite networks
    proportion : float
        Proportion used in finding the metabolite networks, either the
        essential proportion of the synthesis method, or the reaction
        proportion of the consuming method
    synthesis : bool
        Whether to find the synthesis network, if false finds the consuming
        network instead
    processes : int
        Number of processes to use when finding the metabolite network

    Return
    ------
    metabolite_network : pd.DataFrame
        Dataframe describing the metabolite networks with
        metabolites as columns, and reactions as rows
    """
    if out_path.exists():
        return pd.read_csv(out_path, index_col=0)
    if synthesis:
        metabolite_network = (
            metworkpy.metabolites.find_metabolite_synthesis_network_reactions(
                model=model,
                method="essential",
                essential_proportion=proportion,
                progress_bar=False,
                processes=processes,
            )
        )
    else:
        metabolite_network = (
            metworkpy.metabolites.find_metabolite_consuming_network_reactions(
                model=model,
                reaction_proportion=proportion,
                progress_bar=False,
                processes=processes,
                loopless=False,
            )
        )
    metabolite_network = metabolite_network.drop(rxns_to_remove, axis=0)
    metabolite_network.to_csv(out_path, index=True)
    return metabolite_network
