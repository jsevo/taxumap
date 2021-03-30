# Authors: Jonas Schluter <jonas.schluter@nyulangone.org>
# License: BSD 3 clause
__author__ = "Jonas Schluter"
__copyright__ = "Copyright 2020, MIT License"


import unittest
import pandas as pd
import numpy as np
from hctmicrobiomemskcc.tools.microbiotatools import fill_taxonomy_table
from .taxumap import Taxumap


def test_reproducibility_of_embedding(n=500, d=15):
    """Does repeated fitting change the result? """
    rX = (
        pd.DataFrame(np.random.standard_normal((n, d)))
        .applymap(np.abs)
        .apply(lambda r: r / r.sum(), axis=1)
    )
    rX.rename(columns={0: "index_column"}, inplace=True)
    tax = pd.DataFrame(
        [np.ones(d - 1), np.arange(d - 1)],
        index=["Phylum", "Family"],
        columns=np.arange(d - 1) + 1,
    ).T
    tax.iloc[10::, 0] = 2
    tax["Class"] = 1
    tax["Genus"] = 0
    tax = tax.astype(str)

    # tumap 1, transform twice
    tumap = Taxumap(rel_abundances=rX, taxonomy=tax)
    tumap.transform_self()
    em1 = tumap.embedding.copy()
    tumap.transform_self()
    em2 = tumap.embedding.copy()

    # tumap 2, transform twice
    tumap2 = Taxumap(rel_abundances=rX, taxonomy=tax)
    tumap2.transform_self()
    em3 = tumap2.embedding.copy()

    tumap2.transform_self()
    em4 = tumap2.embedding.copy()

    # all embeddings must be close
    return np.allclose(em1, em2, em3, em4)


class TestFillTaxonomyTable(unittest.TestCase):
    """
    basic test class
    """

    def test_fill_taxonomy_table(self):
        """
        The actual test.
        Any method which starts with ``test_`` will considered as a test case.
        """
        taxtable = pd.DataFrame(
            [
                ["Bacteria", "Firmicutes", np.nan, np.nan],
                ["Bacteria", "Firmicutes", "SomeClass", "SomeOrder"],
                ["Bacteria", np.nan, "SomeClass2", "SomeOrder2"],
                ["Bacteria", "Bacteroides", "SomeClass3", np.nan],
                [np.nan, "Bacteroides", "SomeClass3", np.nan],
            ],
            columns=["Kingdom", "Phylum", "Class", "Order"],
            index=["ASV1", "ASV2", "ASV3", "ASV4", "ASV5"],
        )
        filled_taxtable = fill_taxonomy_table(taxtable)
        x1 = filled_taxtable.loc["ASV1"]["Class"]
        x2 = filled_taxtable.loc["ASV1"]["Order"]
        x3 = filled_taxtable.loc["ASV2"]["Class"]
        x4 = filled_taxtable.loc["ASV3"]["Phylum"]
        x5 = filled_taxtable.loc["ASV4"]["Order"]
        x6 = filled_taxtable.loc["ASV5"]["Kingdom"]
        x7 = filled_taxtable.loc["ASV5"]["Order"]
        self.assertEqual(x1, "Phylum_Firmicutes_of_unknown_Class____ASV1")
        self.assertEqual(x2, "Phylum_Firmicutes_of_unknown_Order____ASV1")
        self.assertEqual(x3, "SomeClass")
        self.assertEqual(x4, "Kingdom_Bacteria_of_unknown_Phylum____ASV3")
        self.assertEqual(x5, "Class_SomeClass3_of_unknown_Order____ASV4")
        # a stupid scenario
        self.assertEqual(x6, "unknown_Kingdom____ASV5")

        self.assertEqual(x7, "Kingdom_unknown_Kingdom_of_unknown_Order____ASV5")

    def test_reproducibility(self):
        test_reproducibility_of_embedding(n=1500, d=70)


if __name__ == "__main__":
    unittest.main()
