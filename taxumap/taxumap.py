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

import taxumap.dataloading as parse
import taxumap.tools as tls
import taxumap.visualizations as viz

from umap import UMAP


class Taxumap:
    def __init__(
        self,
        agg_levels=["Phylum", "Family"],
        weight=None,
        taxonomy=None,
        taxonomy_meta=None,
        fpt=None,
        fpx=None,
    ):

        self.agg_levels = agg_levels

        self.fpt = "/Users/granthussey/github/taxumap/taxumap/data/taxonomy.csv"
        self.fpx = "/Users/granthussey/github/taxumap/taxumap/data/microbiota_table.csv"

        # I am pretty sure that my use of If...Else below violates the EAFP principles - should use Try..Except instead

        # Todo: Add in self.fpx is None check, etc.

        if weight is None:
            self.weight = [1] * len(agg_levels)
        else:
            # TODO: Input checks
            self.weight = weight

        if taxonomy is None:
            self.taxonomy = parse.parse_microbiome_data(self.fpx)
        else:
            # TODO: Input checks
            self.taxonomy = taxonomy

        if taxonomy_meta is None:
            self.taxonomy_meta = parse.parse_taxonomy_data(self.fpt)
        else:
            # TODO: Input checks
            self.taxonomy_meta = taxonomy_meta

    def transform(self, debug=False):

        # I hard-coded distanceperlevel for now.
        Xagg = taxonomic_aggregation(
            self.taxonomy, self.taxonomy_meta, self.agg_levels, distanceperlevel=False
        )

        # I hard-coded scaling for now.
        print("Not scaling")
        Xscaled = Xagg

        # Default parameters from old legacy file
        neigh = 120
        min_dist = 0.2
        distance_metric = "braycurtis"

        self.taxumap = UMAP(
            n_neighbors=neigh, min_dist=min_dist, metric=distance_metric
        ).fit(Xscaled)
        self.embedding = self.taxumap.transform(Xscaled)

        if debug:
            self.Xscaled = Xscaled
            self.Xagg = Xagg

        return self

    def __repr__(self):
        return "Taxumap('{}', '{}', '{}', '{}')".format(
            self.agg_levels, self.weight, self.fpx, self.fpt
        )

    def __str__(self):

        # create an if...elif blcok for if fpt exists or not
        # this is a part of the package where I build in the ability to
        # create from file or create from pandas df.

        message = (
            "Taxumap with agg_levels = {} and weights = {}. \n \n".format(
                self.agg_levels, self.weight
            )
            + "The taxonomy and taxonomy_meta data generated from files located \
         at {} and {}, respectively".format(
                self.fpx, self.fpt
            )
        )
        return message

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
            Xagg = tls.aggregate_at_taxlevel(_X, tax, l)
            if distanceperlevel:
                Xagg = ssd.cdist(Xagg, Xagg, distancemetric)
                Xagg = pd.DataFrame(Xagg, index=_X.index, columns=_X.index)
                Xdist = Xdist + Xagg
        X = Xdist
    else:
        for l in agg_levels:
            print("aggregating on %s" % l)
            Xagg = tls.aggregate_at_taxlevel(_X, tax, l)
            X = X.join(Xagg, lsuffix="_r")
        assert np.allclose(_X.sum(axis=1), X.sum(axis=1) / (len(agg_levels) + 1)), (
            "During aggregation, the sum of relative abundances is not equal to %d-times the original relative abundances. This would have been expected due to aggregating and joining"
            % ((len(agg_levels) + 1))
        )
    return X


def taxumap_legacy(
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
        X = parse.parse_microbiome_data(fpx)

    # _cols = X.sum(axis=0).sort_values().tail(100).index
    # X = X[_cols]
    if tax is None:
        tax = parse.parse_taxonomy_data(fpt)

    if asv_colors is None:
        if withusercolors:
            asv_colors = parse.parse_asvcolor_data(fpc)
        else:
            asv_colors = None

    # aggregate taxonomic levels
    Xagg = taxonomic_aggregation(X, tax, agg_levels, distanceperlevel=distanceperlevel)

    # scale
    if withscaling:
        Xscaled = tls.scale(Xagg)
    else:
        print("not scaling")
        Xscaled = Xagg

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
        viz.pretty_print(
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

    _ = taxumap_legacy(
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
