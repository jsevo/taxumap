# Authors: Jonas Schluter <jonas.schluter@nyulangone.org>
# License: BSD 3 clause
__author__ = "Jonas Schluter"
__copyright__ = "Copyright 2020, MIT License"


import unittest
import pandas as pd
import numpy as np
from taxumap import fill_taxonomy_table


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

            ],
            columns=["Kingdom", "Phylum", "Class", "Order"],
            index=["ASV1", "ASV2", "ASV3", "ASV4"],
        )
        filled_taxtable = fill_taxonomy_table(taxtable)
        x1 = filled_taxtable.loc["ASV1"]["Class"]
        x2 = filled_taxtable.loc["ASV1"]["Order"]
        x3 = filled_taxtable.loc["ASV2"]["Class"]
        x4 = filled_taxtable.loc["ASV3"]["Phylum"]

        self.assertEqual(x1,
                         "Phylum_Firmicutes_of_unknown_Class____ASV1")
        self.assertEqual(x2,
                         "Phylum_Firmicutes_of_unknown_Order____ASV1")
        self.assertEqual(x3,
                         "SomeClass")
        self.assertEqual(x4,
                         "Kingdom_Bacteria_of_unknown_Phylum____ASV3")






if __name__ == "__main__":
    unittest.main()
