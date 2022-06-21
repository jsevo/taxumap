# Authors: Jonas Schluter <jonas.schluter@nyulangone.org>, Grant Hussey <grant.hussey@nyulangone.org>
# License: MIT

import os
import sys
import warnings

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.spatial.distance as ssd
from sklearn.preprocessing import MinMaxScaler

# from taxumap.custom_logging import setup_logger

# logger_tools = setup_logger("tools", verbose=False, debug=False)


def tax_agg(rel_abundances, taxonomy, agg_levels, distance_metric, weights, low_precision = False):
    """Generates a distance matrix aggregated on each designated taxon

    Args:
        rel_abundances (Pandas df): Relative abundance df with row-wise compositional data, row: sample, columns: OTU/ASV label
        taxonomy (Pandas df): Row: OTU/ASV label, columns: hierarchy of taxonomy for that ASV/OTU
        agg_levels (list of str): Taxons to aggregate
        distance_metric (str): String to pass to ssd.cdist()
        weights (list of int): Weights of the non-ASV/OTU taxons
        low_precision(float): default False, set threshold for minimum mean relative abundances below which columns are ignored

    Returns:
        pandas df: distance table, row and columns are sample ids
    """

    _X = rel_abundances.copy()
    # remove columns that are always zero
    _X = _X.loc[:, (_X != 0).any(axis=0)]
    if low_precision:
        _threshold = low_precision*_X.shape[0]
        _Xs = _Xs.sum(axis=0)>_threshold
        _X = _X.loc[:,_Xs.values]
    # do a PCA first? memoize?
    Xdist = ssd.cdist(_X, _X, distance_metric)
    Xdist = pd.DataFrame(Xdist, index=_X.index, columns=_X.index)

    for agg_level, weight in zip(agg_levels, weights):
        warnings.warn("aggregating on %s" % agg_level)
        Xagg = aggregate_at_taxlevel(_X, taxonomy, agg_level)
        Xagg = ssd.cdist(Xagg, Xagg, distance_metric)
        Xagg = pd.DataFrame(Xagg, index=_X.index, columns=_X.index)
        Xagg = Xagg * weight

        Xdist = Xdist + Xagg

    return Xdist

def aggregate_at_taxlevel(X, tax, level):
    """Sum relative abundances of all members of a taxon at a given taxonomic level"""
    _X_agg = X.copy()
    _X_agg.columns = [tax.loc[x][level] for x in _X_agg.columns]
    _X_agg = _X_agg.groupby(_X_agg.columns, axis=1).sum()
    try:
        assert np.allclose(
            _X_agg.sum(axis=1), 1.0
        ), "At taxonomic aggregation level %s, rows do not sum to 1."
    except AssertionError:
        print("At taxonomic aggregation level %s, rows do not sum to 1. Moving on anyway")


    return _X_agg


def scale(X, scaler=MinMaxScaler(), remove_rare_asv_level=0):
    """Min max scaling of relative abundances to ensure that different taxonomic levels have comparable dynamic ranges.
    Params
    ===============
    X: ASV table
    scaler: one of the sklearn.preprocessing scalers, defaults to MinMaxScaler

    Returns
    ===============
    Xscaled: scaled ASV table
    """
    # X_sum = X.sum()
    X_stats = X.apply(["max"]).T

    if remove_rare_asv_level > 0:
        # if an ASV has never reached at least `remove_rare_asv_level` threshold, ignore.
        X_consider = X_stats.applymap(lambda v: v > remove_rare_asv_level).apply(np.any, axis=1)
        X_consider = X_consider[X_consider.values].index
    else:
        X_consider = X.columns

    Xscaled = scaler.fit_transform(X[X_consider])

    return Xscaled
