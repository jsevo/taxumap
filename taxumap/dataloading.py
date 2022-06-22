import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# from taxumap.custom_logging import setup_logger
# logger_dataloading = setup_logger("dataloading")


def parse_microbiome_data(fp, idx_col="index_column", idx_dtype=str):
    """Attempt to load a `microbiota_table.csv` file from specified file location.

    Args:
        fp (str): Path to the microbiota_table.csv file
        idx_col (str, optional): The column in microbiota_table.csv to be used as the index. Should be either ASV/OTU. Defaults to "index_column".
        idx_dtype (dtype, optional): Data type for the index. For consistency, str is a good choice. Defaults to str.
    Returns:
        Pandas df: a properly-formatted microbiota_data df
    """
    fp = Path(fp)
    # There's probably a better way of doing this test -
    # fp.resolve(strict=True) will error if file does not exist?
    try:
        _ = fp.resolve(strict=True)

    except FileNotFoundError:
        warnings.warn(
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

            check_if_compositional(X, name="microbiota (microbiota_data)")

            return X.fillna(0)

        except ValueError:
            warnings.warn(
                "Please make sure the microbiota_table has a column labeled 'index_column' which contains the sample IDs"
            )
        except:
            warnings.warn(
                "An unknown error occurred during microbiota_table parsing. Please see the instructions for how to run taxumap."
            )


def parse_taxonomy_data(fp, idx_col=["ASV", "OTU"], idx_dtype=str):
    """Load the taxonomy data."""

    fp = Path(fp)

    # This tests specifically if a path exists, regardless if the path points to a file or a directory.

    # To follow EAFP guidelines, we will run this and then attempt to simply read in the file via Pandas.

    try:

        # See if supplied file exists
        try:
            fp = fp.resolve(strict=True)

        except FileNotFoundError:
            warnings.warn(
                "If using default settings, please make sure\n \
                that the taxonomy table is located in the 'data/' subfolder and named 'taxonomy.csv\n \
                Otherwise, please make sure that you are directing taxumap\n\
                    to the proper location of your taxonomy file."
            )
            sys.exit(2)

        # Read file into Pandas df
        try:
            tax = pd.read_csv(fp)

            warnings.warn(
                "Reading taxonomy table. Assuming columns are ordered by phylogeny with in descending order of hierarchy:\n \
                e.g. Kingdom, Phylum, ... , Genus, Species, etc.\n \
                Additionally, the OTU or ASV column must be labeled as 'OTU' or 'ASV' unless otherwise specified"
            )

        except (IsADirectoryError, pd.errors.ParserError, FileNotFoundError):
            warnings.warn(
                "Your file path may be pointing to a file that doesn't exist, or \
                simply to a directory. Please re-check your file path."
            )

            sys.exit(2)

        # Set index to ASV/OTU
        try:
            if "ASV" in tax.columns:
                tax.set_index("ASV", inplace=True)
                # tax.index.name = "SEQ"

            else:
                tax.set_index("OTU", inplace=True)

                # TODO: Change this to "SEQ" after I make sure this wont break anything (i dont think it will)

                tax.index.name = "ASV"
            tax.columns = map(str.capitalize, tax.columns)

        except ValueError:

            warnings.warn(
                "Neither ASV nor OTU not found in your columns. This is required as base level of Taxonomy. A fix can be to rename your lowest level to ASV."
            )

        else:
            check_tax_is_consistent(tax)
            tax.index = tax.index.astype(idx_dtype)
            return tax

    except:
        warnings.warn(
            "An unknown error occurred during taxonomy table parsing. Please see the documentation for how to run taxumap."
        )


def check_if_compositional(X, name=""):
    """Check if all rows sum to 1."""
    if not np.allclose(X.sum(axis=1), 1):
        warnings.warn(
            "Rows in the {} dataframe do not sum to 1. Is this intentional?".format(name)
        )


def check_tax_is_consistent(df):
    """Check if taxonomy tbl has nans, datatypes; warns if non str found."""
    if np.any(df.isna()):
        warnings.warn(
            "Missing values (NaN) found for some taxonomy levels, you should consider filling with higher taxonomic level names. Please consult the documentation for best way to move forward."
        )
    # should we just read_csv() & only allow certain datatypes (i.e. str, obj?)
    if any(df.dtypes == (int, float)):
        warnings.warn(
            "Your taxonomy table contains columns may contain numerical data. Please consult documentation because you may have incorrectly formatted your dataframe."
        )


def parse_asvcolor_data(fp):
    """Load the taxonomy data."""
    """Todo: This function only works rn with 'ASV' as the 
             index. This should become a parameter"""

    if type(fp) is str:

        fp = Path(fp)
        try:
            f = fp.resolve(strict=True)
        except FileNotFoundError:
            warnings.warn(
                "The color per ASV table should be located in the data/ subfolder and named asvcolors.csv"
            )
            sys.exit(2)
        if f.is_file():
            try:
                taxcolors = pd.read_csv(fp)
                warnings.warn("Reading color per ASV table.")
                try:
                    assert taxcolors.columns[[0, 1]].to_list() == ["ASV", "HexColor"]
                except AssertionError:
                    warnings.warn(
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
            except ValueError:
                warnings.warn(
                    "Please make sure the taxcolors has columns labeled ['ASV','HexColor'], and contain as values the ASV labels as strings and valid hex color stings"
                )
            except:
                warnings.warn(
                    "Please make sure the taxcolors has columns labeled ['ASV','HexColor'], and contain as values the ASV labels as strings and valid hex color stings"
                )
    elif type(fp) is pd.DataFrame:
        warnings.warn("using provided taxcolors")
        taxcolors = fp
        return taxcolors
