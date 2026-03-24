"""
Some common functionality used across multiple scripts
"""

import os
import pathlib
from collections.abc import Callable
from typing import (
    Any,
    AnyStr,
    Dict,
    Hashable,
    List,
    Literal,
    Optional,
    Union,
    cast,
)

import cobra

# External Imports
import escher
import metworkpy
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import MinMaxScaler, StandardScaler


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


KEGG_REST_PREFIX = r"https://rest.kegg.jp/"


def get_kegg_pathway_genes(
    organism_code: str, strip: bool = True
) -> pd.DataFrame:
    """
    Get a DataFrame describing the genes in the KEGG pathways for
    the desired organism

    Parameters
    ----------
    organism_code : str
        The KEGG organism code for the organism of interest,
        for example code for *Mycobacterium tuberculosis*  is
        'mtu'
    strip : bool, default=True
        Whether to remove common prefixes from the
        gene and pathway identifiers, specifically
        removes 'path:' from pathway identifiers and
        '<organism code>:' from gene identifiers

    Returns
    -------
    kegg_pathway_df : pd.DataFrame
        Dataframe describing the KEGG pathways, with a column
        (pathway) with the pathway identifier, and a column
        (gene) with the gene identifier.
    """
    pathway_df = pd.read_csv(
        f"{KEGG_REST_PREFIX}link/pathway/{organism_code}",
        sep="\t",
        names=["gene", "pathway"],
    )
    if strip:
        pathway_df["gene"] = pathway_df["gene"].str.replace(
            f"{organism_code}:", ""
        )
        pathway_df["pathway"] = pathway_df["pathway"].str.replace("path:", "")
    return pathway_df


def get_kegg_pathway_descriptions(
    organism_code: str, remove_str: str
) -> pd.DataFrame:
    """
    Get a DataFrame describing the KEGG pathways for
    the desired organism

    Parameters
    ----------
    organism_code : str
        The KEGG organism code for the organism of interest,
        for example the code for *Mycobacterium tuberculosis*
        is mtu
    remove_str : str
        A string to remove from the pathway descriptions,
        will be treated as a regex, see
        `pandas documentation<https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.replace.html>`_
        for replace. The descriptions are then stripped of whitespace at the
        start and end of the string.

    Returns
    -------
    kegg_pathway_descriptions : pd.DataFrame
        A dataframe with the KEGG pathway descriptions,
        with a column (pathway) with the KEGG
        pathway identified, and another column
        (description) with the description
        of the pathway
    """
    pathway_description_df = pd.read_csv(
        f"{KEGG_REST_PREFIX}list/pathway/{organism_code}",
        sep="\t",
        names=["pathway", "description"],
    )
    pathway_description_df["description"] = (
        pathway_description_df["description"]
        .str.replace(remove_str, "")
        .str.strip()
    )
    return pathway_description_df


def fdr_with_nan(ps: pd.Series, **kwargs) -> pd.Series:
    """
    Perform false discovery correction on a series of p-values which may contain
    missing values

    Parameters
    ----------
    ps : pd.Series
        Pandas series containing p-values to correct for false discovery
    kwargs
        Key word arguments passed to scipy.stats.false_discovery_control

    Returns
    -------
    corrected_ps : pd.Series
        Series with corrected ps, with the same index as the input series
    """
    corrected_ps = pd.Series(np.nan, index=ps.index)
    no_nan_series = ps.dropna(inplace=False)
    corrected_ps[no_nan_series.index] = stats.false_discovery_control(
        no_nan_series, **kwargs
    )
    return corrected_ps


def bootstrap(
    data: pd.Series,
    sample_set: list[Hashable],
    statistic: Callable[[pd.Series], float],
    iterations: int = 1_000,
    method: Literal["kernel", "empirical"] = "kernel",
    alternative: Literal["two-sided", "less", "greater"] = "two-sided",
    seed: Optional[int] = None,
) -> tuple[float, float]:
    """
    Perform a bootstrap significance test to determine if the statistic of the
    sample set is significantly greater than expected by chance

    Parameters
    ----------
    data : pd.Series
        The data used for testing, should include all the values in the sample set
        as well as all the additional data not in the sample set
    sample_set : list of Hashable
        The samples included in the sample set, should be the same type as the
        index of the data series, and all values should be included in the
        data series index
    statistic : function taking a pandas series and returning a float
        The statistic to bootstrap the significance for, must accept a pandas Series,
        and return a float which is the value of the statistic
    iterations : int, default=1000
        The number of bootstrap iterations to perform
    method : "kernel" or "empirical"
        Method to use when estimating the p-value, kernel uses Gaussian Kernel
        Density to estimate the p-value, while emprical uses an emprical CDF
    alternative : "two-sided", "less", or "greater"
        What should be used as an alternative, less tests whether the statistic
        is less than expected by chance, greater tests whether the statistic is
        greater than expected by chance, and two-sided tests if it is either

    Returns
    -------
    statistic : float
        The value of the statistic for the sample set
    p-value : float
        The estimated p-value for the statistic
    """
    sample_statistic = statistic(data[sample_set])  # type: ignore
    data_idx = list(data.index)
    # Create the random number generator
    rng = np.random.default_rng(seed=seed)
    # Create the bootstrap distribution
    bootstrap_statistic_distribution = np.zeros((iterations,), dtype=float)
    for iteration in range(iterations):
        bootstrap_statistic_distribution[iteration] = statistic(
            data[
                list(rng.choice(data_idx, size=len(sample_set), replace=True))
            ]
        )
    # Create the desired distribution
    if method == "kernel":
        distribution = stats.gaussian_kde(bootstrap_statistic_distribution)
        prob_less = cast(
            float, distribution.integrate_box_1d(-np.inf, sample_statistic)
        )
        prob_greater = cast(
            float, distribution.integrate_box_1d(sample_statistic, np.inf)
        )
    elif method == "empirical":
        # Apply an adjustment based on 'Permutation p-values should never be zero:
        # calculating exact p-values when permutations are randomly drawn'
        # First find eps, which is floating point tolerance
        # Also based on Scipy's permutation test implementation
        eps = (
            0
            if not np.isdtype(
                bootstrap_statistic_distribution.dtype, "real floating"
            )
            else np.finfo(bootstrap_statistic_distribution.dtype).eps * 100
        )
        gamma = np.abs(eps * sample_statistic)
        prob_less = float(
            np.count_nonzero(
                bootstrap_statistic_distribution <= sample_statistic + gamma
            )
            + 1
        ) / float(iterations + 1)
        prob_greater = float(
            np.count_nonzero(
                bootstrap_statistic_distribution >= sample_statistic - gamma
            )
            + 1
        ) / float(iterations + 1)
        empirical_distribution = stats.ecdf(bootstrap_statistic_distribution)
        prob_less: float = empirical_distribution.cdf.evaluate(
            sample_statistic
        )
        prob_greater: float = empirical_distribution.sf.evaluate(
            sample_statistic
        )
    else:
        raise ValueError(
            f"method must be 'kernel' or 'emprical' but received {method}"
        )
    if alternative == "less":
        return sample_statistic, prob_less
    elif alternative == "greater":
        return sample_statistic, prob_greater
    elif alternative == "two-sided":
        return sample_statistic, 2 * min(prob_less, prob_greater)
    else:
        raise ValueError(
            f"alternative must be either 'less', 'greater', or 'two-sided', but received {alternative}"
        )


# Local Imports

# Create a pathlike type
PathLike = Union[str, os.PathLike, pathlib.Path]


# Function to read in an escher map, add reaction
# and metabolite data to it, then save the map
def escher_map_add_data(
    input_map: PathLike,
    output_dir: PathLike,
    output_prefix: AnyStr,
    reaction_data: Optional[pd.Series] = None,
    reaction_data_scaling: Optional[Literal["minmax", "standard"]] = None,
    metabolite_data: Optional[pd.Series] = None,
    metabolite_data_scaling: Optional[Literal["minmax", "standard"]] = None,
    reaction_scale: Optional[List[Dict[AnyStr, Any]]] = None,
    metabolite_scale: Optional[List[Dict[AnyStr, Any]]] = None,
):
    """
    Add reaction and metabolite data to an Escher map and save the resulting HTML

    Parameters
    ----------
    input_map : PathLike
        Path to the map to update with reaction and metabolite data
    reaction_data : pd.Series, optional
        Series containing reaction data to add to the Escher map, whose index
        matches the reaction IDS found in the map
    reaction_data_scaling : 'minmax' or 'standard'
        Whether to use MinMaxScaler or StandardScaler from scikit-learn
        on the reaction data, default is to not perform scaling
    metabolite_data : pd.Series, optional
        Series containing metabolite data to add to the Escher map, whose index
        matches the metabolite IDS found in the map
    metabolite_data_scaling : 'minmax' or 'standard'
        Whether to use MinMaxScaler or StandardScaler from scikit-learn
        on the metabolite data, default is to not perform scaling
    reaction_scale : list of dict from str to Any
        Used to update the reaction scale of the builder if not None, see
        `Escher documentation<https://escher.readthedocs.io/en/latest/escher-python.html>`_
        for more details
    metabolite_scale : list of dict from str to Any
        Used to update the metabolite scale of the builder if not None, see
        `Escher documentation<https://escher.readthedocs.io/en/latest/escher-python.html>`_
        for more details
    """
    # Convert the input map path into a Pathlib Path
    input_map_path = pathlib.Path(input_map)
    map_name = input_map_path.stem
    # Scale the reaction/metabolite data as needed
    if reaction_data_scaling is not None and reaction_data is not None:
        if reaction_data_scaling == "minmax":
            scaler = MinMaxScaler()
        elif reaction_data_scaling == "standard":
            scaler = StandardScaler()
        else:
            raise ValueError(
                f"Invalid reaction scaling choice, must be either minmax or standard but received {reaction_data_scaling}"
            )
        scaler.set_output(transform="pandas")
        reaction_data = scaler.fit_transform(
            reaction_data.to_frame("rxn_data")
        )["rxn_data"]
    if metabolite_data_scaling is not None and metabolite_data is not None:
        if metabolite_data_scaling == "minmax":
            scaler = MinMaxScaler()
        elif metabolite_data_scaling == "standard":
            scaler = StandardScaler()
        else:
            raise ValueError(
                f"Invalid metabolite scaling choice, must be either minmax or standard but received {metabolite_data_scaling}"
            )
        scaler.set_output(transform="pandas")
        metabolite_data = scaler.fit_transform(
            metabolite_data.to_frame("met_data")
        )["met_data"]
    # Create the building loading in the map
    builder = escher.Builder(map_json=str(input_map_path))
    # Change the scroll behaviour
    builder.scroll_behavior = "zoom"
    # Stop the absolute value being used for visualizaiton
    builder.reaction_styles = ["color", "size", "text"]
    # Add the reaction data and metabolite data to the map
    if reaction_data is not None:
        builder.reaction_data = reaction_data.to_dict()
    if metabolite_data is not None:
        builder.metabolite_data = metabolite_data.to_dict()
    # Add the reaction and metabolite scales if provided
    if reaction_scale:
        builder.reaction_scale = reaction_scale
    if metabolite_scale:
        builder.metabolite_scale = metabolite_scale
    # Save the map in html format
    output_dir = pathlib.Path(output_dir)
    builder.save_html(str(output_dir / f"{output_prefix}{map_name}.html"))
