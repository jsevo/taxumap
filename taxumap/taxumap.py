# Authors: Jonas Schluter <jonas.schluter@nyulangone.org>, Grant Hussey <grant.hussey@nyulangone.org>
# License: BSD 3 clause # Which license?

__author__ = ["Jonas Schluter", "Grant Hussey"]
__copyright__ = "Copyright 2020, MIT License"

import os
import sys
import warnings
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
import scipy.spatial.distance as ssd
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from umap import UMAP

import taxumap.dataloading as parse
import taxumap.tools as tls
import taxumap.visualizations as viz


class Taxumap:
    """Taxumap object for running taxUMAP algorithm"""

    def __init__(
        self,
        agg_levels=["Phylum", "Family"],
        weight=None,
        rel_abundances=None,
        taxonomy=None,
        fpt=None,
        fpx=None,
        name=None,
    ):

        """Constructor method for the Taxumap object

        Args:
            agg_levels (list, optional): Determines which taxonomic levels to aggregate on. See taxUMAP documentation. Defaults to ["Phylum", "Family"].
            weight (list of int, optional): Determines the weight of each of the agg_levels. Length of this list must match length of agg_levels. Defaults to None.
            rel_abundances (dataframe, optional): Compositional data (relative counts) of abundance of ASV/OTU. Defaults to None.
            taxonomy (dataframe, optional): Dataframe with index of ASV/OTU and columns representing decsending taxonomic levels. Defaults to None.
            fpt (str, optional): Filepath to the rel_abundances meta dataframe, if saved on disk. Defaults to None.
            fpx (str, optional): Filepath to the rel_abundances dataframe, if saved on disk. Defaults to None.
            name (str, optional): A useful name for the project. Used in graphing and saving methods. Defaults to None.
        """

        self.agg_levels = list(map(lambda x: x.capitalize(), agg_levels))

        # I am pretty sure that my use of If...Else below violates the EAFP principles - should use Try..Except instead

        if weight is None:
            self.weight = [1] * len(agg_levels)
        else:
            if len(weight) != len(agg_levels):
                raise ValueError(
                    "The length of the weight must match that of agg_levels"
                )
            else:
                self.weight = weight

        # Set rel_abundances df
        try:
            if rel_abundances is None and fpx is None:
                raise ValueError
            elif isinstance(rel_abundances, pd.DataFrame):
                print("Recognized `rel_abundances` parameter as Pandas DataFrame")
                self.fpx = None
                self.rel_abundances = rel_abundances
                parse.check_if_compositional(
                    rel_abundaces, name="locally-supplied microbiota (rel_abundances)"
                )
            elif isinstance(fpx, str):
                self.fpx = fpx
                self.rel_abundances = parse.parse_microbiome_data(fpx)
            else:
                raise NameError
        except (ValueError, NameError) as e:
            _name_value_error(e, "fpx", "rel_abundances")
            raise

        # Set taxonomy df
        try:
            if taxonomy is None and fpt is None:
                raise ValueError
            elif isinstance(taxonomy, pd.DataFrame):
                print("Recognized `taxonomy` parameter as Pandas DataFrame")
                self.fpt = None
                self.taxonomy = taxonomy
                self.taxonomy.columns = map(str.capitalize, self.taxonomy.columns)
                parse.check_tax_is_consistent(self.taxonomy)
            elif isinstance(fpx, str):
                self.fpt = fpt
                self.taxonomy = parse.parse_taxonomy_data(fpt)
            else:
                raise NameError
        except (ValueError, NameError) as e:
            _name_value_error(e, "fpt", "taxonomy")
            raise

        # Set name attribute
        if isinstance(name, str):
            self.name = name
        else:
            self.name = None

    @property
    def _exist_tax_meta_df(self):
        return isinstance(self.taxonomy, pd.DataFrame)

    @property
    def _exist_tax_meta_fp(self):
        """Returns true if the `fpt` parameter is a valid filepath"""
        try:
            logic = Path(self.fpt).is_file()
        except TypeError:
            return False
        else:
            return logic

    @property
    def _exist_tax_df(self):
        return isinstance(self.rel_abundances, pd.DataFrame)

    @property
    def _exist_tax_fp(self):
        """Returns true if the `fpx` parameter is a valid filepath"""
        try:
            logic = Path(self.fpx).is_file()
        except TypeError:
            return False
        else:
            return logic

    @property
    def _is_df_loaded(self):
        """Returns true if the object was initialized with df objects"""
        return all((self._exist_tax_df, self._exist_tax_meta_df)) and not all(
            (self._exist_tax_fp, self._exist_tax_meta_fp)
        )

    @property
    def _is_fp_loaded(self):
        """Returns true if the object was initialized with filepaths, or
        if all such pre-requisites exist"""
        return all(
            (
                self._exist_tax_df,
                self._exist_tax_fp,
                self._exist_tax_meta_fp,
                self._exist_tax_meta_df,
            )
        )

    @property
    def _embedded_csv_name(self):
        """Filename for saving self.embedded to a csv file. Uses self.name if available."""
        if isinstance(self.name, str):
            return "_".join([self.name, "embedded.csv"])
        else:
            return "embedded.csv"

    @property
    def _embedded_pickle_name(self):
        """Filename for saving self to a pickle. Uses self.name if available."""
        if isinstance(self.name, str):
            return "_".join([self.name, "taxumap_pickle.pickle"])
        else:
            return "taxumap_pickle.pickle"

    @property
    def taxumap1(self):
        """Get a label for first dimension of embedded space"""
        first_letters = [agg_level[0] for agg_level in self.agg_levels]
        initials = "".join(first_letters)
        return "taxumap-{}-1".format(initials)

    @property
    def taxumap2(self):
        """Get a label for second dimension of embedded space"""
        first_letters = [agg_level[0] for agg_level in self.agg_levels]
        initials = "".join(first_letters)
        return "taxumap-{}-2".format(initials)

    @property
    def _is_transformed(self):
        try:
            if isinstance(self.embedding, np.ndarray):
                return True
            else:
                warnings.warn(
                    "taxumap.embedding is not an ndarray, something went wrong"
                )
        except AttributeError:
            return False

    @property
    def df_embedding(self):

        if self._is_transformed:
            return pd.DataFrame(
                self.embedding,
                columns=[self.taxumap1, self.taxumap2],
                index=self.index,
            )
        else:
            raise AttributeError(
                "Please run Taxumap.transform_self() to generate your embedding"
            )

    @property
    def df_dominant_taxon(self):

        # DataFrame where each row is the maximum taxon corresponding to the
        # shared index in the column index_column
        df_dom_tax_per_sample = (
            pd.DataFrame(
                self.rel_abundances.idxmax(axis="columns"), columns=["max_tax"]
            )
            .reset_index()
            .set_index("max_tax")
        )

        # This is where a merge is done to add in the full taxonomical
        # data for each of the "max_tax" dominant taxon
        prelim_df_table = df_dom_tax_per_sample.merge(
            self.taxonomy, right_index=True, left_index=True
        )

        # all below here, I am cleaning the prelim_df_table for final return
        tax_hierarchy = prelim_df_table.columns[
            prelim_df_table.columns != "index_column"
        ]

        new_tax_hierarchy = [
            "dom_" + each_level.lower() for each_level in tax_hierarchy
        ]

        change_labels = dict(zip(tax_hierarchy, new_tax_hierarchy))

        df_final = prelim_df_table.rename(columns=change_labels)
        df_final.index.name = "max_tax"
        df_final = df_final.reset_index().set_index("index_column")

        return df_final

    def transform_self(
        self, scale=False, debug=False, save=False, outdir=None, pickle=False
    ):
        """If rel_abundances and taxonomy dataframes are available, will run the taxUMAP transformation.

        Args:
            debug (bool, optional): If True, self will be given X and Xscaled variables (debug only). Defaults to False.
            save (bool, optional): If True, will attempt to save the resulting embedded as a csv file. Defaults to False.
            outdir (str, optional): Path to where files should be saved, if save=True. Defaults to None.
            pickle (bool, optional): If True, will save self object as a pickle. Defaults to False.
        """

        # Shouldn't need `try...except` because any Taxumap object should have proper attributes
        Xagg = taxonomic_aggregation(
            self.rel_abundances, self.taxonomy, self.agg_levels, distanceperlevel=False,
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
        # self._is_transformed = True
        self.index = Xscaled.index

        if debug:
            self.Xscaled = Xscaled
            self.Xagg = Xagg

        if save:
            try:
                outdir = Path(outdir).resolve(strict=True)
            except (FileNotFoundError, TypeError) as e:
                print(e)
                print(
                    '\nNo valid outdir was declared.\nSaving data into "./results" folder.\n'
                )
                outdir = Path("./results").resolve()
            except Exception as e:
                throw_unknown_save_error(e)

            if outdir.is_dir():
                self.save_it(outdir, pickle)
            else:
                try:
                    os.mkdir(outdir)
                except FileExistsError as e:
                    print(e)
                    print(
                        '\nUnable to automatically save embedding, most likely do to a file called "results" in your PWD.\n'
                    )
                    sys.exit(2)
                else:
                    self.save_it(outdir, pickle)
        return self

    def save_it(self, outdir=None, pickle=False):

        if outdir is None:
            print("Saving to ./results")
            outdir = Path("./results").resolve()
            try:
                os.mkdir(outdir)
                print("Making ./results folder...")
            except FileExistsError:
                print("./results folder already exists")

        try:
            pd.DataFrame(
                self.embedding, columns=["taxUMAP1", "taxUMAP2"], index=self.index
            ).to_csv(os.path.join(outdir, self._embedded_csv_name))
        except AttributeError as e:
            print(e)
            print(
                "\nEmbedding not currently populated. Please run taxumap.Taxumap.transform_self(save=True).\n"
            )
        except Exception as e:
            throw_unknown_save_error(e)
        else:
            if pickle:
                import pickle

                with open(os.path.join(outdir, self._embedded_pickle_name), "wb") as f:
                    pickle.dump(self, f)

    def scatter(self, figsize=(16, 10), save=False, **kwargs):

        if not self._is_transformed:
            raise AttributeError(
                "Your Taxumap has yet to be transformed. Run Taxumap.transform_self() first."
            )

        fig, ax = plt.subplots(figsize=figsize)

        ax.scatter(
            self.df_embedding[self.df_embedding.columns[0]],
            self.df_embedding[self.df_embedding.columns[1]],
            **kwargs
        )
        ax.set_xlabel(self.taxumap1)
        ax.set_ylabel(self.taxumap2)

        sns.despine()

        return fig, ax

    def __repr__(self):

        if self._is_df_loaded:
            return "Taxumap(agg_levels = {}, weight = {}, rel_abundances = '{}', taxonomy = '{}')".format(
                self.agg_levels,
                self.weight,
                "loaded from local scope",
                "loaded from local scope",
            )

        elif self._is_fp_loaded:
            return "Taxumap(agg_levels = {}, weight = {}, fpx = '{}', fpt = '{}')".format(
                self.agg_levels, self.weight, self.fpx, self.fpt
            )

        else:
            return "Taxumap(agg_levels = {}, weight = {}, fpx = '{}', fpt = '{}')".format(
                self.agg_levels, self.weight, self.fpx, self.fpt
            )

    def __str__(self):

        # create an if...elif blcok for if fpt exists or not
        # this is a part of the package where I build in the ability to
        # create from file or create from pandas df.

        messages = [
            "Taxumap with agg_levels = {} and weights = {}.".format(
                self.agg_levels, self.weight
            )
        ]

        if self._is_df_loaded:
            messages.append(
                "The `rel_abundances` and `taxonomy` dataframes were passed in from local variables."
            )
            return "\n \n".join(messages)

        elif self._is_fp_loaded:
            messages.append(
                "The rel_abundances and taxonomy dataframes were generated from files located at\n'{}'\nand\n'{}',\nrespectively".format(
                    self.fpx, self.fpt
                )
            )
            return "\n \n".join(messages)

        else:
            return repr(self)


def throw_unknown_save_error(e):
    print(e)
    print("\nUnknown error has occured. Cannot save embedding as instructed.\n")
    sys.exit(2)


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

    if agg_levels == None:
        try:
            assert (
                len(tax.columns) < 3
            ), "the taxonomy table has very few columns. Cannot aggregate taxonomic levels. Reverting to regular UMAP"

            # Automaticaly get the second and second-to-last agg_levels
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

        # Check to see if data was aggregated properly
        assert np.allclose(_X.sum(axis=1), X.sum(axis=1) / (len(agg_levels) + 1)), (
            "During aggregation, the sum of relative abundances is not equal to %d-times the original relative abundances. This would have been expected due to aggregating and joining"
            % ((len(agg_levels) + 1))
        )
    return X


def _name_value_error(e, fp_param, df_param):
    print(e)
    print()
    print(
        "Please provide the constructor with one of the following: a filepath to your {} file via parameter `{}`, or with an initialized variable for your {} dataframe via the parameter `{}`".format(
            df_param, fp_param, df_param, df_param
        )
    )


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
    """Please use a '/' delimiter for weight and agg_levels"""
    import argparse

    parser = argparse.ArgumentParser(description="Get options for taxumap run")
    parser.add_argument(
        "-m", "--microbiota_data", help="Microbiota (rel_abundances) Table"
    )
    parser.add_argument("-t", "--taxonomy", help="Taxonomy Table")
    parser.add_argument("-w", "--weight", help="Weights")
    parser.add_argument("-a", "--agg_levels", help="Aggregation Levels")

    args = parser.parse_args()

    # try:
    #     agg_levels = list(map(lambda x: x.capitalize(), args.agg_levels.split("/")))
    # except AttributeError:
    #     agg_levels = args.agg_levels

    # AGG_LEVELS
    if args.agg_levels is not None:
        agg_levels = list(map(lambda x: x.capitalize(), args.agg_levels.split("/")))
    else:
        agg_levels = args.agg_levels

    # WEIGHT
    if args.weight is not None:
        weights = args.weight.split("/")
        weights = list(map(lambda x: int(x), weights))
    else:
        weights = args.weight

    # TAX META
    if args.taxonomy is not None:
        fpt = args.taxonomy
    else:
        fpt = "./data/taxonomy.csv"

    # rel_abundances
    if args.microbiota_data is not None:
        fpx = args.microbiota_data
    else:
        fpx = "./data/microbiota_table.csv"

    t = Taxumap(agg_levels=agg_levels, weight=weights, fpt=fpt, fpx=fpx)

    # TODO: opt to ask for saving, default yes
    t.transform_self()
