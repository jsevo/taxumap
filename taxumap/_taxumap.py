from pathlib import Path
import pandas as pd
import numpy as np

import sys
import os

from taxumap.custom_logging import setup_logger
from taxumap.errors import throw_unknown_save_error, _name_value_error


logger_taxumapmixins = setup_logger("taxumap_mixin", verbose=False, debug=False)


class TaxumapMixin:
    @property
    def _exist_tax_meta_df(self):
        return isinstance(self.taxonomy, pd.DataFrame)

    @property
    def _exist_tax_df(self):
        return isinstance(self.rel_abundances, pd.DataFrame)

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
    def _plot_name(self, ext="png"):
        """Filename for saving a plot to file. Uses self.name if available."""
        if isinstance(self.name, str):
            return "_".join([self.name, "plot." + ext])
        else:
            return "plot." + ext

    @classmethod
    def from_pickle(cls, fp):
        import pickle

        try:
            with open(fp, "rb") as f:
                data = pickle.load(f)
        except Exception:
            logger_taxumapmixins.exception("Something went wrong loading from pickle")
        else:
            logger_taxumapmixins.info("Successfully located pickle file")

        return data

    def to_pickle(self, outdir=".", name=None):

        import pickle

        if name is None:
            name = self._embedded_pickle_name

        with open(os.path.join(outdir, name), "wb") as f:
            pickle.dump(self, f)

        return self

    @property
    def taxumap1(self):
        """Get a label for first dimension of taxUMAP embedded space"""
        first_letters = [agg_level[0] for agg_level in self.agg_levels]
        initials = "".join(first_letters)
        return "taxumap-{}-1".format(initials)

    @property
    def taxumap2(self):
        """Get a label for second dimension of taxUMAP embedded space"""
        first_letters = [agg_level[0] for agg_level in self.agg_levels]
        initials = "".join(first_letters)
        return "taxumap-{}-2".format(initials)

    @property
    def _is_transformed(self):
        """Returns True if Taxumap.transform_self() has been previously run"""
        try:
            if isinstance(self.embedding, np.ndarray):
                return True
            else:
                logger_taxumapmixins.warning(
                    "taxumap.embedding is not an ndarray, something went wrong"
                )
        except AttributeError:
            return False

    def save_embedding(self, outdir=None):

        try:
            _save(self.df_embedding.to_csv, outdir, self._embedded_csv_name)
        except AttributeError as e:
            logger_taxumapmixins.warning(
                "\nEmbedding not currently populated. Please run taxumap.Taxumap.transform_self(save=True).\n"
            )
        except Exception as e:
            throw_unknown_save_error(e)

        return self


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
            logger_taxumapmixins.warning(
                '\nNo valid outdir was declared.\nSaving data into "./results" folder.\n'
            )

    elif outdir is None:
        outdir = Path("./results").resolve()
        try:
            os.mkdir(outdir)
            logger_taxumapmixins.info("Making ./results folder...")
        except FileExistsError:
            logger_taxumapmixins.info("./results folder already exists")
        except Exception as e:
            throw_unknown_save_error(e)
            sys.exit(2)

    try:
        fxn(os.path.join(outdir, filename), **kwargs)
    except Exception as e:
        throw_unknown_save_error(e)
    else:
        logger_taxumapmixins.info("Save successful")

