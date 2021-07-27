# Authors: Jonas Schluter <jonas.schluter@nyulangone.org>, Grant Hussey <grant.hussey@nyulangone.org>
# License: MIT

import os
import shutil
import unittest
from pathlib import Path
import logging

import taxumap.taxumap as t

# The data for unittests is from Olin and Axel:
# Olin, Axel (2018), “Stereotypic Immune System Development in Newborn Children”, Mendeley Data, v1


# setup logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
stream_format = logging.Formatter("%(funcName)s:%(levelname)s\n%(message)s\n")
sh = logging.StreamHandler()
sh.setFormatter(stream_format)
sh.setLevel(logging.WARNING)  # this level and ABOVE


class TestTaxumap(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Makes temp directory to save everything"""
        unique_folder_name = Path("./temp").resolve()  # maybe beef this up later lol
        logger.warning("Making temp directory at {}.".format(unique_folder_name))

        try:
            os.mkdir(unique_folder_name)
        except FileExistsError:
            try:
                os.rmdir(unique_folder_name)
            except OSError:
                logger.critical(
                    "Pease delete the temp directory in taxumap/taxumap/unittests"
                )

    @classmethod
    def tearDownTest(cls):
        """"Deletes temp directory"""
        unique_folder_name = Path("./temp").resolve()  # maybe beef this up later lol
        shutil.rmtree(unique_folder_name)
        pass

    def setUp(self):
        self.broke_name = t.Taxumap(
            fpt="taxonomy.csv", fpx="microbiota_table.csv", name=3465,
        )

        self.t1 = t.Taxumap(
            fpt="taxonomy.csv", fpx="microbiota_table.csv", name="test1",
        )

        self.t2 = t.Taxumap(
            fpt="taxonomy.csv", fpx="microbiota_table.csv", name="test2_transformed",
        )

        self.t2.transform_self()

    def test_interactive_loading(self):

        with self.assertRaises(NameError):
            t.Taxumap(rel_abundances={}, taxonomy={})

        with self.assertRaises(ValueError):
            t.Taxumap(rel_abundances=None, taxonomy=None)

        self.assertIsNone(self.broke_name.name)

    def test_fp_loading(self):
        pass

    def test_is_loaded_from_logic(self):
        pass


if __name__ == "__main__":
    unittest.main()
