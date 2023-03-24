# Authors: Jonas Schluter <jonas.schluter@nyulangone.org>, Grant Hussey <grant.hussey@nyulangone.org>
# License: MIT

__author__ = ["Jonas Schluter", "Grant Hussey"]
__copyright__ = "Copyright 2020, MIT License"

import warnings
import os

# import logging

import matplotlib as mpl

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from umap import UMAP

from taxumap.tools import tax_agg
from taxumap.input_validation import validate_inputs

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


class Taxumap:
    """Taxumap object for running TaxUMAP algorithm."""

    def __init__(
        self,
        agg_levels=["Phylum", "Family"],
        weights=None,
        microbiota_data=None,
        taxonomy=None,
        name=None,
        random_state=42,
    ):
        """Instantiate the Taxumap object.

        Args:
            agg_levels (list, optional): Determines which taxonomic levels to aggregate on. See taxUMAP documentation. Defaults to ["Phylum", "Family"].
            weights (list of int, optional): Determines the weights of each of the agg_levels. Length of this list must match length of agg_levels. Defaults to None.
            microbiota_data (dataframe, optional): Compositional data (relative counts) of abundance of ASV/OTU. Defaults to None.
            taxonomy (dataframe, optional): Dataframe with index of ASV/OTU and columns representing descending taxonomic levels. Defaults to None.
            fpt (str, optional): Filepath to the rel_abundances meta dataframe, if saved on disk. Defaults to None.
            fpx (str, optional): Filepath to the rel_abundances dataframe, if saved on disk. Defaults to None.
            name (str, optional): A useful name for the project. Used in graphing and saving methods. Defaults to None.
        """
        self.random_state = random_state
        self._is_transformed = False
        self.agg_levels = list(map(lambda x: x.capitalize(), agg_levels))

        weights, microbiota_data, taxonomy = validate_inputs(
            weights, microbiota_data, taxonomy, agg_levels)
        print("post validate inputs main", taxonomy)
        self.weights = weights
        self.microbiota_data = microbiota_data
        self.taxonomy = taxonomy

        # Set name attribute
        if isinstance(name, str):
            self.name = name
        else:
            self.name = None

    def transform_self(self,
                       scale=False,
                       debug=False,
                       distance_metric="braycurtis",
                       **kwargs):
        """Run the taxUMAP transformation.

        Requires microbiota_data and taxonomy dataframes.

        Args:
            debug (bool, optional): If True, self will be given X and Xagg variables (debug only). Defaults to False.
        """
        if ("neigh" not in kwargs) or (kwargs["neigh"] is None):
            warnings.warn(
                "Please set neigh parameter to approx. the size of individals in the dataset. See documentation."
            )
            neigh = (120 if len(self.microbiota_data) > 120 else
                     len(self.microbiota_data) - 1)
        else:
            neigh = kwargs["neigh"]

        if ("min_dist" not in kwargs) or (kwargs["min_dist"] is None):
            warnings.warn("Setting min_dist to 0.05/sum(weights)")
            min_dist = 0.05 / np.sum(self.weights)
        else:
            min_dist = kwargs["min_dist"]

        if ("epochs" not in kwargs) or (kwargs["epochs"] is None):
            epochs = (5000 if neigh < 120 else
                      (1000 if len(self.microbiota_data) < 5000 else 1000))
            warnings.warn("Setting epochs to %d" % epochs)
        else:
            epochs = kwargs["epochs"]

        if "low_precision" not in kwargs:
            low_precision = False
        else:
            low_precision = np.float(kwargs["low_precision"])

        # Shouldn't need `try...except` because any Taxumap object should have proper attributes
        Xagg = tax_agg(
            self.microbiota_data,
            self.taxonomy,
            self.agg_levels,
            distance_metric,
            self.weights,
            low_precision=low_precision,
        )

        rs = np.random.RandomState(seed=self.random_state)

        if self._is_transformed:
            print(
                "TaxUMAP has already been fit. Re-running could yield a different embedding due to random state changes betweeen first and second run. Re-seeding random number generator RandomState."
            )

            rs = np.random.RandomState(seed=self.random_state)

        with warnings.catch_warnings():
            # this will suppress all warnings in this block
            warnings.simplefilter("ignore")

            self.taxumap = UMAP(
                n_neighbors=neigh,
                min_dist=min_dist,
                n_epochs=epochs,
                metric="precomputed",
                #            transform_seed=1,
                random_state=rs,
            ).fit(Xagg)

        # fit embedding as numpy array
        self.embedding = self.taxumap.transform(Xagg)
        # raw data sample index as taxumap property
        self.index = Xagg.index
        # fit embedding as data frame
        self.df_embedding = pd.DataFrame(self.embedding,
                                         columns=["taxumap1", "taxumap2"],
                                         index=self.index)
        self._is_transformed = True

        if debug:
            self.Xagg = Xagg

        return self

    @property
    def df_dominant_taxon(self):
        """Return the taxon with max abundance in sample."""
        # DataFrame where each row is the maximum taxon corresponding to the
        # shared index in the column index_column
        df_dom_tax_per_sample = (pd.DataFrame(
            self.microbiota_data.idxmax(axis="columns"),
            columns=["max_tax"]).reset_index().set_index("max_tax"))

        # This is where a merge is done to add in the full taxonomical
        # data for each of the "max_tax" dominant taxon
        prelim_df_table = df_dom_tax_per_sample.merge(self.taxonomy,
                                                      right_index=True,
                                                      left_index=True)

        # all below here, I am cleaning the prelim_df_table for final return
        tax_hierarchy = prelim_df_table.columns[
            prelim_df_table.columns != "index_column"]

        new_tax_hierarchy = [
            "dom_" + each_level.lower() for each_level in tax_hierarchy
        ]

        change_labels = dict(zip(tax_hierarchy, new_tax_hierarchy))

        df_final = prelim_df_table.rename(columns=change_labels)
        df_final.index.name = "max_tax"
        df_final = df_final.reset_index().set_index("index_column")

        return df_final

    def save_embedding(self, path=None):
        """Write embedding to file."""
        if path is None:
            path = "taxumap_embedding.csv"

        try:
            self.df_embedding.to_csv(path)
        except Exception as e:
            print(e)

        return self

    def scatter(
        self,
        figsize=(6, 4),
        save=False,
        outdir=None,
        plotname="taxumap_scatterplot.pdf",
        **kwargs,
    ):
        """Scatter plot."""
        if not self._is_transformed:
            raise AttributeError(
                "Taxumap has yet to be transformed. Run Taxumap.transform_self"
            )

        fig, ax = plt.subplots(figsize=figsize)

        ax.scatter(
            self.df_embedding[self.df_embedding.columns[0]],
            self.df_embedding[self.df_embedding.columns[1]],
            **kwargs,
        )

        ax.set_xlabel("TaxUMAP-1")
        ax.set_ylabel("TaxUMAP-2")

        ax.set_title(self.name)

        if save:
            fig.savefig(os.path.join(outdir, plotname))

        return fig, ax

    def __repr__(self):
        return f"Taxumap(agg_levels = {self.agg_levels}, weights = {self.weights})"

    def __str__(self):
        return (
            f"Taxumap with agg_levels = {self.agg_levels} and weights = {self.weights}."
        )
