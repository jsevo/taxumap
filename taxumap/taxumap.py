# Authors: Jonas Schluter <jonas.schluter@nyulangone.org>
# License: BSD 3 clause
__author__ = "Jonas Schluter"
__copyright__ = "Copyright 2020, MIT License"

import os

#!/usr/bin/env python
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.spatial.distance as ssd
from sklearn.preprocessing import MinMaxScaler

# Ideally should have these installed
# from hctmicrobiomemskcc.tools.microbiotatools import fill_taxonomy_table


def aggregate_at_taxlevel(X, tax, level):
    """Helper function. For a given taxonomic level, aggregate relative abundances by summing all members of corresponding taxon."""
    _X_agg = X.copy()
    _X_agg.columns = [tax.loc[x][level] for x in _X_agg.columns]
    _X_agg = _X_agg.groupby(_X_agg.columns, axis=1).sum()
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


def parse_microbiome_data(fp, idx_col="index_column", idx_dtype=str):
    """Load the microbiota data"""

    fp = Path(fp)

    # There's probably a better way of doing this test -
    # fp.resolve(strict=True) will error if file does not exist?
    try:
        f = fp.resolve(strict=True)

    except FileNotFoundError as fe:
        print("{0}".format(fe))
        print(
            "The microbiota composition table should be located in the data/ subfolder and named microbiota_table.csv"
        )
        sys.exit(2)

    # Why can't the above statement just be folded into this?
    if fp.is_file():
        try:

            # This must be implemented in a two-liner to make sure dtype
            # of index is str
            X = pd.read_csv(fp, dtype={idx_col: idx_dtype})
            X.set_index(idx_col, inplace=True)
            X = X.astype(np.float64)

            if not np.allclose(X.sum(axis=1), 1):
                warnings.warn("rows do not sum to 1. Is this intentional?")

            return X.fillna(0)

        except ValueError as ve:
            print("{0}".format(ve))
            print(
                "Please make sure the microbiota_table has a column labeled 'index_column' which contains the sample IDs"
            )
        except:
            print(
                "An unknown error occurred during microbiota_table parsing. Please see the instructions for how to run taxumap."
            )


def parse_taxonomy_data(fp):
    """Load the taxonomy data."""
    """Todo: Make it so this fxn takes 'idx_col' parameter instead of forcing
             users to set the ASV/OTU col at the end"""

    fp = Path(fp)

    try:

        f = fp.resolve(strict=True)

    except FileNotFoundError as fe:
        print("{0}".format(fe))
        print(
            "The taxonomy table should be located in the data/ subfolder and named taxonomy.csv"
        )

        sys.exit(2)

    if f.is_file():

        try:
            tax = pd.read_csv(fp,)

            print(
                "Reading taxonomy table. Assuming columns are ordered by phylogeny with in descending order of hierarchy."
            )
            print("e.g. Kingdom, Phylum, ... , Genus, Species, etc")

            print("Make sure that the OTU/ASV is the LAST COLUMN")

            _low = tax.columns[-1]
            try:
                assert _low.lower() in [
                    "asv",
                    "otu",
                ], "%s not in ['ASV', 'OTU'], moving on"

            except AssertionError as ae:
                print(ae)

            print(
                "Setting %s as index, THIS ASSUMES %s IS THE LOWEST TAXONOMIC LEVEL"
                % (_low, _low)
            )

            tax = tax.set_index(_low)
            if np.any(tax.isna()):
                warnings.warn(
                    "Missing values (NaN) found for some taxonomy levels, you should consider filling with higher taxonomic level names /n \
                    Please consult the documentation for best way to move forward"
                )

                # tax = fill_taxonomy_table(tax)

            return tax

        except ValueError as ve:
            print("{0}".format(ve))
            # print(
            #     "Please make sure the taxonomy has a column labeled 'index_column' which contains the sample IDs"
            # )
            print("Error - please make sure that the OTU/ASV is the LAST COLUMN")

        except:
            print(
                "An unknown error occurred during taxonomy table parsing. Please see the instructions for how to run taxumap."
            )


def parse_asvcolor_data(fp):
    """Load the taxonomy data."""
    """Todo: This function only works rn with 'ASV' as the 
             index. This should become a parameter"""

    if type(fp) is str:

        fp = Path(fp)
        try:
            f = fp.resolve(strict=True)
        except FileNotFoundError as fe:
            print("{0}".format(fe))
            print(
                "The color per ASV table should be located in the data/ subfolder and named asvcolors.csv"
            )
            sys.exit(2)
        if f.is_file():
            try:
                taxcolors = pd.read_csv(fp)
                print("Reading color per ASV table.")
                try:
                    assert taxcolors.columns[[0, 1]].to_list() == ["ASV", "HexColor"]
                except AssertionError:
                    print(
                        'Column names should be:  ["ASV", "HexColor"]. Choosing colors automatically.'
                    )
                    return ()

                taxcolors = taxcolors.set_index("ASV")
                if np.any(taxcolors.isna()):
                    warnings.warn(
                        "Missing values (NaN) found for some taxcolors. Filling with 'grey'"
                    )
                    taxcolors = taxcolors.fillna("grey")
                return taxcolors
            except ValueError as ve:
                print("{0}".format(ve))
                print(
                    "Please make sure the taxcolors has columns labeled ['ASV','HexColor'], and contain as values the ASV labels as strings and valid hex color stings"
                )
            except:
                print(
                    "Please make sure the taxcolors has columns labeled ['ASV','HexColor'], and contain as values the ASV labels as strings and valid hex color stings"
                )
    elif type(fp) is pd.DataFrame:
        print("using provided taxcolors")
        taxcolors = fp
        return taxcolors


def taxonomic_aggregation(
    X, tax, agg_levels, distanceperlevel=False, distancemetric="braycurtis"
):
    """
    Expand the original microbiota sample data by aggregating taxa at different higher taxonomic levels.

    Parameters
    ----------

    X : pandas.DataFrame, microbiota data table at lowers taxonomic level

    tax : pandas.DataFrame, taxonomy table

    agg_level : list; which taxonomic levels to include. Defaults to the second highest (usually phylum) and third lowest (usually family).

    distanceperlevel: bool, should sample-to-sample distances be calculated per tax level

    Returns 
    ----------

    X : Expanded microbiota table.

    """
    assert (type(agg_levels) == list) | (
        agg_levels == None
    ), "Aborting: agg_levels should be a list of taxonomic levels or explicitly `None`"
    import numpy as np

    if agg_levels == None:
        try:
            assert (
                len(tax.columns) > 3
            ), "the taxonomy table has very few columns. Cannot aggregate taxonomic levels. Reverting to regular UMAP"
            agg_levels = np.array(tax.columns)[[1, -2]]
            print(agg_levels)
            if agg_levels[0] == agg_levels[1]:
                # in case the taxonomy table is weird and has few columns
                agg_levels = agg_levels[0]
        except AssertionError as ae:
            print("{0}".format(ae))

    _X = X.copy()
    if distanceperlevel:
        Xdist = ssd.cdist(_X, _X)
        Xdist = pd.DataFrame(Xdist, index=_X.index, columns=_X.index)
        for l in agg_levels:
            print("aggregating on %s" % l)
            Xagg = aggregate_at_taxlevel(_X, tax, l)
            if distanceperlevel:
                Xagg = ssd.cdist(Xagg, Xagg, distancemetric)
                Xagg = pd.DataFrame(Xagg, index=_X.index, columns=_X.index)
                Xdist = Xdist + Xagg
        X = Xdist
    else:
        for l in agg_levels:
            print("aggregating on %s" % l)
            Xagg = aggregate_at_taxlevel(_X, tax, l)
            X = X.join(Xagg, lsuffix="_r")
        assert np.allclose(_X.sum(axis=1), X.sum(axis=1) / (len(agg_levels) + 1)), (
            "During aggregation, the sum of relative abundances is not equal to %d-times the original relative abundances. This would have been expected due to aggregating and joining"
            % ((len(agg_levels) + 1))
        )
    return X


def taxumap(
    agg_levels,
    withscaling,
    distanceperlevel,
    distancemetric,
    print_figure,
    print_with_diversity,
    X=None,
    tax=None,
    asv_colors=None,
    loadembedding=False,
    withusercolors=False,
    bgcolor="white",
    fpx="data/microbiota_table.csv",
    fpt="data/taxonomy.csv",
    fpc="data/asvcolors.csv",
    transform_seed=None,
    debug=False,
    save_embedding=False,
    neigh=120,
    min_dist=0.2,
):
    """
    taxumap embedding of microbiota composition data. Data is expected to be compositional, i.e. each row sums to 1. Two tables are required: the microbiota data and a taxonomy table.

    The microbiota data file ('data/microbiota_table.csv') must have a column with sample indeces, labeled 'index_column'. The remaining columns are expected to be the lowest level taxa (OTU/ASV/...):

    index_column, ASV1, ASV2
    'sample1', 0.5, 0.5
    'sample2', 0.2, 0.8

    The taxonomy table ('data/taxonomy.csv') is expected to resolve higher taxonomic groups for the columns in the microbiota table. The columns of this table should contain taxonomic levels. They should be ordered from left to right in decreasing taxonomic hierarchy, e.g.

    kingdom, phylum, ..., ASV
    'Bacteria', 'Firmicutes', ..., 'ASV1'

    The data is expected to be located in the data/ folder. Results will be written to the results/ folder

    Parameters
    ----------
    agg_level : list; which taxonomic levels to include. Defaults to the second highest (usually phylum) and third lowest (usually family).

    print_figure : bool, flag if pretty scatter plot of taxumap embedded data should be produced.

    print_with_diversity : bool, flag if the scatter plot should have a background indicating alpha-diversity of samples in the region (inverse Simpson).

    fpx : File path to microbiota compositional data table in csv table format. Defaults to ./data/microbiota_table.csv.

    fpt : File path to the taxonomy table in csv table format. Defaults to ./data/taxonomy.csv. 

    fpc : (optional) File path to the color per ASV table in csv table format. Defaults to ./data/asvcolors.csv. 

    Returns
    ----------

    TAXUMAP: UMAP object fit to data

    X_embedded : pandas.DataFrame of 2-D coordinates with indeces from the microbiota_table.l

    """

    if X is None:
        X = parse_microbiome_data(fpx)

    # _cols = X.sum(axis=0).sort_values().tail(100).index
    # X = X[_cols]
    if tax is None:
        tax = parse_taxonomy_data(fpt)

    if asv_colors is None:
        if withusercolors:
            asv_colors = parse_asvcolor_data(fpc)
        else:
            asv_colors = None

    # aggregate taxonomic levels
    Xagg = taxonomic_aggregation(X, tax, agg_levels, distanceperlevel=distanceperlevel)

    # scale
    if withscaling:
        Xscaled = scale(Xagg)
    else:
        print("not scaling")
        Xscaled = Xagg

    from umap import UMAP

    if loadembedding:
        X_embedded = pd.read_csv("results/embedded.csv", index_col="index_column")
        embedding = X_embedded.values
        TAXUMAP = {}

    else:
        if withscaling:
            TAXUMAP = UMAP(
                n_neighbors=neigh, min_dist=min_dist, metric="manhattan"
            ).fit(Xscaled)

        elif distanceperlevel:
            TAXUMAP = UMAP(
                n_neighbors=neigh, min_dist=min_dist, metric="precomputed"
            ).fit(Xscaled)
        else:
            TAXUMAP = UMAP(
                n_neighbors=neigh, min_dist=min_dist, metric=distancemetric
            ).fit(Xscaled)

        embedding = TAXUMAP.transform(Xscaled)

        if transform_seed is not None:
            TAXUMAP.transform_seed = transform_seed

        X_embedded = pd.DataFrame(
            embedding, index=X.index, columns=["phyUMAP-1", "phyUMAP-2"]
        )

        if save_embedding:

            if not os.path.isdir("results"):
                print()
                print("No results folder in current working directory.")
                print("Will create one.")
                os.mkdir("results")

            X_embedded.to_csv("results/embedded.csv")

    if print_figure:
        # diversity per sample
        # if not np.all(X.sum(axis=1) > 0):
        #     warnings.warn("some rows have zero sum. adding tiny value (1e-16)")
        #     X = X + 1e-16
        if print_with_diversity:
            ivs = X.apply(lambda r: 1 / np.sum(r ** 2), axis=1)
        else:
            ivs = X.iloc[:, 1]
        pretty_print(
            X,
            embedding,
            ivs,
            tax,
            with_diversity_background=print_with_diversity,
            usercolors=asv_colors,
            bgcolor=bgcolor,
        )

    if debug:
        return (TAXUMAP, X_embedded, Xscaled, X)
    else:
        return (TAXUMAP, X_embedded)


if __name__ == "__main__":
    import sys, getopt

    """
    use with options:

    -e OR WITH --scalingpertaxlevel
    -p OR WITH --print_figure
    -s OR WITH --with_diversity_background
    -a OR --asvcolors 
    -c OR --scatterbgcolor : a valid matplotlib color for the background of the scatter plot, default: white
    -l OR --loadembedding : load already calculated embedding from result/embedding.csv
    """
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(
            argv,
            "epsalc:",
            [
                "scalingpertaxlevel",
                "print_figure",
                "with_diversity_background",
                "asvcolors",
                "loadembedding",
                "scatterbgcolor=",
            ],
        )

    except getopt.GetoptError:
        sys.exit("unknown options")

    agg_levels = None
    withscaling = False
    distanceperlevel = False
    distancemetric = "braycurtis"
    print_figure = False
    distanceperlevel = False
    withusercolors = False
    with_diversity_background = False
    loadembedding = False
    bgcolor = "white"

    for opt, arg in opts:
        if opt in ("-e", "--scalingpertaxlevel"):
            withscaling = True
        elif opt in ("-p", "--print_figure"):
            print_figure = True
        elif opt in ("-s", "--with_diversity_background"):
            with_diversity_background = True
        elif opt in ("-c", "--scatterbgcolor"):
            bgcolor = str(arg)
        elif opt in ("-a", "--asvcolors"):
            withusercolors = True
        elif opt in ("-l", "--loadembedding"):
            loadembedding = True

    _ = taxumap(
        agg_levels,
        withscaling,
        distanceperlevel,
        distancemetric,
        print_figure,
        withusercolors=withusercolors,
        print_with_diversity=with_diversity_background,
        loadembedding=loadembedding,
        bgcolor=bgcolor,
        save_embedding=True,
    )
