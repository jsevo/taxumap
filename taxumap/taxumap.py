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
import taxumap.tools as tools
import taxumap.visualizations as viz
from taxumap.custom_logging import setup_logger
from taxumap._taxumap import TaxumapMixin, _save
from taxumap.errors import throw_unknown_save_error, _name_value_error

logger_taxumap = setup_logger("taxumap", verbose=False, debug=False)


class Taxumap(TaxumapMixin):
    """Taxumap object for running taxUMAP algorithm"""

    def __init__(
        self,
        agg_levels=["Phylum", "Family"],
        weights=None,
        rel_abundances=None,
        taxonomy=None,
        fpt=None,
        fpx=None,
        name=None,
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

        self.agg_levels = list(map(lambda x: x.capitalize(), agg_levels))

        # I am pretty sure that my use of If...Else below violates the EAFP principles - should use Try..Except instead

        if weights is None:
            self.weights = [1] * len(agg_levels)
        else:
            if len(weights) != len(agg_levels):
                raise ValueError(
                    "The length of the weights must match that of agg_levels"
                )
            else:
                self.weights = weights

        # Set rel_abundances df
        try:
            if rel_abundances is None and fpx is None:
                raise ValueError
            elif isinstance(rel_abundances, pd.DataFrame):
                logger_taxumap.info(
                    "Recognized `rel_abundances` parameter as Pandas DataFrame"
                )
                self.fpx = None
                self.rel_abundances = rel_abundances

                # Validate the rel_abundances dataset
                try:
                    self.rel_abundances.set_index("index_column")
                except KeyError:
                    # if it can't set the index to 'index_column', maybe it's already set
                    # let's check to see if that is not the case
                    if self.rel_abundances.index.name != "index_column":
                        logger_taxumap.exception(
                            "Your rel_abundances df needs to contain 'index_column' containing sample identifiers of each row containing relative OTU/ASV abundances"
                        )
                        sys.exit(2)
                    else:
                        logger_taxumap.info(
                            "Your rel_abundances dataframe already had 'index_column' as the index."
                        )

                else:
                    logger_taxumap.info(
                        "index_column has been set as index for self.rel_abundances"
                    )

                # lastly, quickly check if the data is compositional.
                dataloading.check_if_compositional(
                    self.rel_abundances,
                    name="locally-supplied microbiota (rel_abundances)",
                )

            elif isinstance(fpx, str):
                self.fpx = fpx
                self.rel_abundances = dataloading.parse_microbiome_data(fpx)
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
                logger_taxumap.info(
                    "Recognized `taxonomy` parameter as Pandas DataFrame"
                )
                self.fpt = None
                self.taxonomy = taxonomy
                self.taxonomy.columns = map(str.capitalize, self.taxonomy.columns)
                dataloading.check_tax_is_consistent(self.taxonomy)
            elif isinstance(fpx, str):
                self.fpt = fpt
                self.taxonomy = dataloading.parse_taxonomy_data(fpt)
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
        if "neigh" in kwargs:
            neigh = kwargs["neigh"]
        else:
            # TODO: Add in documentation guideance on neigh
            logger_taxumap.warning(
                "Please set neigh parameter to approx. the size of individals in the dataset. See documentation."
            )
            neigh = 120 if len(self.rel_abundances) > 120 else len(self.rel_abundances)

        if "min_dist" in kwargs:
            min_dist = kwargs["min_dist"]
        else:
            logger_taxumap.info("Setting min_dist to 0.05/sum(weights)")
            min_dist = 0.05 / np.sum(self.weights)

        if "epochs" in kwargs:
            epochs = kwargs["epochs"]
        else:
            epochs = (
                5000
                if neigh < 120
                else (1000 if len(self.rel_abundances) < 5000 else 1000)
            )
            logger_taxumap.info("Setting epochs to %d" % epochs)

        distance_metric = "braycurtis"

        # Shouldn't need `try...except` because any Taxumap object should have proper attributes
        Xagg = tools.tax_agg(
            self.rel_abundances,
            self.taxonomy,
            self.agg_levels,
            distance_metric,
            self.weights,
        )

        self.taxumap = UMAP(
            n_neighbors=neigh, min_dist=min_dist, n_epochs=epochs, metric="precomputed"
        ).fit(Xagg)

        self.embedding = self.taxumap.transform(Xagg)
        # self._is_transformed = True

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

    def scatter(
        self, figsize=(16, 10), save=False, outdir=None, ax=None, fig=None, **kwargs
    ):

        if not self._is_transformed:
            raise AttributeError(
                "Your Taxumap has yet to be transformed. Run Taxumap.transform_self() first."
            )

        if fig or ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        ax.scatter(
            self.df_embedding[self.df_embedding.columns[0]],
            self.df_embedding[self.df_embedding.columns[1]],
            **kwargs
        )
        ax.set_xlabel(self.taxumap1)
        ax.set_ylabel(self.taxumap2)

        ax.set_title(self.name)

        sns.despine(trim=True, offset=5)

        if save:
            _save(fig.savefig, outdir, self._plot_name)

        return fig, ax

    def __repr__(self):

        if self._is_df_loaded:
            return "Taxumap(agg_levels = {}, weights = {}, rel_abundances = '{}', taxonomy = '{}')".format(
                self.agg_levels,
                self.weights,
                "loaded from local scope",
                "loaded from local scope",
            )

        elif self._is_fp_loaded:
            return "Taxumap(agg_levels = {}, weights = {}, fpx = '{}', fpt = '{}')".format(
                self.agg_levels, self.weights, self.fpx, self.fpt
            )

        else:
            return "Taxumap(agg_levels = {}, weights = {}, fpx = '{}', fpt = '{}')".format(
                self.agg_levels, self.weights, self.fpx, self.fpt
            )

    def __str__(self):

        # create an if...elif blcok for if fpt exists or not
        # this is a part of the package where I build in the ability to
        # create from file or create from pandas df.

        messages = [
            "Taxumap with agg_levels = {} and weights = {}.".format(
                self.agg_levels, self.weights
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
