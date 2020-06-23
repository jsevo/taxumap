import os
import shutil
import unittest
from pathlib import Path

import taxumap.taxumap as t


class TestTaxumap(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Makes temp directory to save everything"""
        unique_folder_name = Path("./temp").resolve()  # maybe beef this up later lol
        print("Making temp directory at {}.".format(unique_folder_name))

        try:
            os.mkdir(unique_folder_name)
        except FileExistsError:
            try:
                os.rmdir(unique_folder_name)
            except OSError:
                print("Pease delete the temp directory in taxumap/taxumap/unittests")

    @classmethod
    def tearDownTest(cls):
        """"Deletes temp directory"""
        unique_folder_name = Path("./temp").resolve()  # maybe beef this up later lol
        shutil.rmtree(unique_folder_name)
        pass

    def setUp(self):
        self.broke_name = t.Taxumap(
            fpt="/Users/granthussey/github/taxumap/taxumap/unittests/taxonomy.csv",
            fpx="/Users/granthussey/github/taxumap/taxumap/unittests/microbiota_table.csv",
            name=3465,
        )

        self.t1 = t.Taxumap(
            fpt="/Users/granthussey/github/taxumap/taxumap/unittests/taxonomy.csv",
            fpx="/Users/granthussey/github/taxumap/taxumap/unittests/microbiota_table.csv",
            name="test1",
        )

        self.t2 = t.Taxumap(
            fpt="/Users/granthussey/github/taxumap/taxumap/unittests/taxonomy.csv",
            fpx="/Users/granthussey/github/taxumap/taxumap/unittests/microbiota_table.csv",
            name="test2_transformed",
        )

        self.t2.transform()

    def test_interactive_loading(self):

        with self.assertRaises(NameError):
            t.Taxumap(taxonomy={}, taxonomy_meta={})

        with self.assertRaises(ValueError):
            t.Taxumap(taxonomy=None, taxonomy_meta=None)

        self.assertIsNone(self.broke_name.name)

    def test_fp_loading(self):
        pass

    def test_is_loaded_from_logic(self):
        pass


if __name__ == "__main__":
    unittest.main()
