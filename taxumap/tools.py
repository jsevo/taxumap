import os

#!/usr/bin/env python
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.spatial.distance as ssd
from sklearn.preprocessing import MinMaxScaler


def aggregate_at_taxlevel(X, tax, level):
    """Helper function. For a given taxonomic level, aggregate relative abundances by summing all members of corresponding taxon."""
    _X_agg = X.copy()
    _X_agg.columns = [tax.loc[x][level] for x in _X_agg.columns]
    _X_agg = _X_agg.groupby(_X_agg.columns, axis=1).sum()
    try:
        assert np.allclose(
            _X_agg.sum(axis=1), 1.0
        ), "At taxonomic aggregation level %s, rows do not sum to 1."
    except AssertionError:
        print("moving on anyway")

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
    X_sum = X.sum()
    X_stats = X.apply(["max"]).T

    if remove_rare_asv_level > 0:
        # if an ASV has never reached at least `remove_rare_asv_level` threshold, ignore.
        X_consider = X_stats.applymap(lambda v: v > remove_rare_asv_level).apply(
            np.any, axis=1
        )
        X_consider = X_consider[X_consider.values].index
    else:
        X_consider = X.columns

    Xscaled = scaler.fit_transform(X[X_consider])

    return Xscaled
