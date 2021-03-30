# Authors: Jonas Schluter <jonas.schluter@nyulangone.org>, Grant Hussey <grant.hussey@nyulangone.org>
# License: BSD 3 clause # Which license?

__author__ = ["Jonas Schluter", "Grant Hussey"]
__copyright__ = "Copyright 2020, MIT License"

import os
import sys
import warnings
import logging
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt
from umap import UMAP

import taxumap.dataloading as dataloading
from taxumap.tools import *
from taxumap.input_validation import validate_inputs
import taxumap.visualizations as viz
from taxumap.custom_logging import setup_logger
from taxumap._taxumap import TaxumapMixin, _save
from taxumap.errors import throw_unknown_save_error, _name_value_error

logger_taxumap = setup_logger("taxumap", verbose=False, debug=False)



class Taxumap:
    """Taxumap object for running TaxUMAP algorithm"""

    def __init__(
        self,
        agg_levels=["Phylum", "Family"],
        weights=None,
        rel_abundances=None,
        taxonomy=None,
        name=None,
        random_state=42,
    ):

        """Constructor method for the Taxumap object

        Args:
            agg_levels (list, optional): Determines which taxonomic levels to aggregate on. See taxUMAP documentation. Defaults to ["Phylum", "Family"].
            weights (list of int, optional): Determines the weights of each of the agg_levels. Length of this list must match length of agg_levels. Defaults to None.
            rel_abundances (dataframe, optional): Compositional data (relative counts) of abundance of ASV/OTU. Defaults to None.
            taxonomy (dataframe, optional): Dataframe with index of ASV/OTU and columns representing descending taxonomic levels. Defaults to None.
            fpt (str, optional): Filepath to the rel_abundances meta dataframe, if saved on disk. Defaults to None.
            fpx (str, optional): Filepath to the rel_abundances dataframe, if saved on disk. Defaults to None.
            name (str, optional): A useful name for the project. Used in graphing and saving methods. Defaults to None.
        """
        self.random_state = random_state
        self._is_transformed =False
        self.agg_levels = list(map(lambda x: x.capitalize(), agg_levels))


        weights, rel_abundances, taxonomy = validate_inputs(
            weights, rel_abundances, taxonomy, agg_levels, logger_taxumap
        )
        self.weights = weights
        self.rel_abundances = rel_abundances
        self.taxonomy = taxonomy

        # Set name attribute
        if isinstance(name, str):
            self.name = name
        else:
            self.name = None

    def transform_self(
        self, scale=False, debug=False, save=False, outdir=None, **kwargs
    ):
        """If rel_abundances and taxonomy dataframes are available, will run the taxUMAP transformation.

        Args:
            debug (bool, optional): If True, self will be given X and Xagg variables (debug only). Defaults to False.
            save (bool, optional): If True, will attempt to save the resulting embedded as a csv file. Defaults to False.
            outdir (str, optional): Path to where files should be saved, if save=True. Defaults to None.
            pickle (bool, optional): If True, will save self object as a pickle. Defaults to False.
        """

        # Maybe better way of implementing this

        if "neigh" not in kwargs:
            logger_taxumap.warning(
                "Please set neigh parameter to approx. the size of individals in the dataset. See documentation."
            )
            neigh = (
                120 if len(self.rel_abundances) > 120 else len(self.rel_abundances) - 1
            )
        else:
            neigh = kwargs["neigh"]

        if "min_dist" not in kwargs:
            logger_taxumap.info("Setting min_dist to 0.05/sum(weights)")
            min_dist = 0.05 / np.sum(self.weights)
        else:
            min_dist = kwargs["min_dist"]

        if "epochs" not in kwargs:
            epochs = (
                5000
                if neigh < 120
                else (1000 if len(self.rel_abundances) < 5000 else 1000)
            )
            logger_taxumap.info("Setting epochs to %d" % epochs)
        else:
            epochs = kwargs["epochs"]

        distance_metric = "braycurtis"

        # Shouldn't need `try...except` because any Taxumap object should have proper attributes
        Xagg = tax_agg(
            self.rel_abundances,
            self.taxonomy,
            self.agg_levels,
            distance_metric,
            self.weights,
        )

        rs = np.random.RandomState(seed=self.random_state)

        if self._is_transformed:
            print(
                "TaxUMAP has already been fit. Re-running could yield a different embedding due to random state changes betweeen first and second run. Re-starting RandomState."
            )
            rs = np.random.RandomState(seed=self.random_state)

        self.taxumap = UMAP(
            n_neighbors=neigh,
            min_dist=min_dist,
            n_epochs=epochs,
            metric="precomputed",
#            transform_seed=1,
            random_state=rs,
        ).fit(Xagg)
        self._is_transformed = True

        self.embedding = self.taxumap.transform(Xagg)
        self.index = Xagg.index

        if debug:
            self.Xagg = Xagg

        if save:
            self.save_embedding(outdir=outdir)

        return self

    @property
    def df_embedding(self):
        """Creates a dataframe from the embedding generated by Taxumap.transform_self()"""
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

    @classmethod
    def from_pickle(cls, fp):
        import pickle

        try:
            with open(fp, "rb") as f:
                data = pickle.load(f)
        except Exception:
            logger_taxumap.exception("Something went wrong loading from pickle")
        else:
            logger_taxumap.info("Successfully located pickle file")

        return data

    def to_pickle(self, outdir=".", name=None):

        import pickle

        if name is None:
            name = self._embedded_pickle_name

        with open(os.path.join(outdir, name), "wb") as f:
            pickle.dump(self, f)

        return self

    def save_embedding(self, outdir=None):

        try:
            _save(self.df_embedding.to_csv, outdir, self._embedded_csv_name)
        except AttributeError as e:
            logger_taxumap.warning(
                "\nEmbedding not currently populated. Please run taxumap.Taxumap.transform_self(save=True).\n"
            )
        except Exception as e:
            throw_unknown_save_error(e)

        return self

    def scatter(
        self, figsize=(16, 10), save=False, outdir=None, ax=None, fig=None, **kwargs
    ):
        # TODO I would like this removed. There is a visualizations module which has already got appropriate functions. should not be a method of this class.

        if not self._is_transformed:
            raise AttributeError(
                "Your Taxumap has yet to be transformed. Run Taxumap.transform_self() first."
            )

        if (fig is None) or (ax is None):
            fig, ax = plt.subplots(figsize=figsize)

        ax.scatter(
            self.df_embedding[self.df_embedding.columns[0]],
            self.df_embedding[self.df_embedding.columns[1]],
            **kwargs,
        )
        ax.set_xlabel(self.taxumap1)
        ax.set_ylabel(self.taxumap2)

        ax.set_title(self.name)

        sns.despine(trim=True, offset=5)

        if save:
            _save(fig.savefig, outdir, self._plot_name)

        return fig, ax

    def __repr__(self):
        return f"Taxumap(agg_levels = {self.agg_levels}, weights = {self.weights})"

    def __str__(self):
        return (
            f"Taxumap with agg_levels = {self.agg_levels} and weights = {self.weights}."
        )



def _save(fxn, outdir, filename, **kwargs):
    """[summary]

    Args:
        fxn (function handle): f(Path-like object or str)
        outdir ([type]): [description]
        filename ([type]): [description]

    Raises:
        TypeError: [description]
    """

    if not callable(fxn):
        raise TypeError("'fxn' passed is not callable")

    if outdir is not None:
        try:
            outdir = Path(outdir).resolve(strict=True)
        except (FileNotFoundError, TypeError) as e:
            logger_taxumap.warning(
                '\nNo valid outdir was declared.\nSaving data into "./results" folder.\n'
            )

    elif outdir is None:
        outdir = Path("./results").resolve()
        try:
            os.mkdir(outdir)
            logger_taxumap.info("Making ./results folder...")
        except FileExistsError:
            logger_taxumap.info("./results folder already exists")
        except Exception as e:
            throw_unknown_save_error(e)
            sys.exit(2)

    try:
        fxn(os.path.join(outdir, filename), **kwargs)
    except Exception as e:
        throw_unknown_save_error(e)
    else:
        logger_taxumap.info("Save successful")


def throw_unknown_save_error(e):
    logger_taxumap.exception(
        "\nUnknown error has occured. Cannot save embedding as instructed."
    )
    sys.exit(2)
