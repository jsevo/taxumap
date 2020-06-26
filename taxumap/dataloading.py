import os

#!/usr/bin/env python
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.spatial.distance as ssd
from sklearn.preprocessing import MinMaxScaler


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

            check_if_compositional(X, name="microbiota (rel_abundances)")

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


def check_if_compositional(X, name=""):
    if not np.allclose(X.sum(axis=1), 1):
        warnings.warn(
            "Rows in the {} dataframe do not sum to 1. Is this intentional?".format(
                name
            )
        )


def parse_taxonomy_data(fp, idx_col=["ASV", "OTU"]):
    """Load the taxonomy data."""

    fp = Path(fp)

    # This tests specifically if a path exists, regardless if the path points to a file or a directory.

    # To follow EAFP guidelines, we will run this and then attempt to simply read in the file via Pandas.

    try:

        try:
            fp = fp.resolve(strict=True)

        except FileNotFoundError as fe:
            print("{0}".format(fe))
            print(
                "If using default settings, please make sure /n \
                that the taxonomy table is located in the 'data/' subfolder and named 'taxonomy.csv /n'"
            )
            print(
                "Otherwise, please make sure that you are directing taxumap \
                    to the proper location of your taxonomy file."
            )
            sys.exit(2)

        try:
            tax = pd.read_csv(fp)

            print()
            print(
                "Reading taxonomy table. Assuming columns are ordered by phylogeny with in descending order of hierarchy:"
            )
            print("e.g. Kingdom, Phylum, ... , Genus, Species, etc.")
            print()
            print("The OTU or ASV column must be labeled as 'OTU' or 'ASV'.")

        except (IsADirectoryError, pd.errors.ParserError, FileNotFoundError) as e:
            print(e)
            print(
                "Your file path may be pointing to a file that doesn't exist, or \
                simply to a directory. Please re-check your file path."
            )
            sys.exit(2)

        try:

            idx_lowest_level = tax.columns[
                tax.columns.str.lower().isin([x.lower() for x in idx_col])
            ]

            tax.set_index(idx_lowest_level.to_list(), inplace=True)

        except ValueError as e:
            print(e)
            print("ASV/OTU not found in your columns")
            sys.exit(2)

        else:

            if np.any(tax.isna()):
                warnings.warn(
                    "Missing values (NaN) found for some taxonomy levels, you should consider filling with higher taxonomic level names. Please consult the documentation for best way to move forward."
                )

            # should we just read_csv() and only allow certain datatypes (i.e. str, obj?)
            if any(tax.dtypes == (int, float)):
                warnings.warn(
                    "Your taxonomy table contains columns may contain numerical data. Please consult documentation because you may have incorrectly formatted your dataframe."
                )

            return tax

    except:
        print(
            "An unknown error occurred during taxonomy table parsing. Please see the documentation for how to run taxumap."
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
