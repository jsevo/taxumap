# Authors: Jonas Schluter <jonas.schluter@nyulangone.org>, Grant Hussey <grant.hussey@nyulangone.org>
# License: MIT
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="taxumap",
    version="0.1",
    description="UMAP visualization for microbiota compositions with taxonomic structure.",
    url="http://github.com/jsevo/taxumap",
    author="Jonas Schluter",
    author_email="jonas.schluter+github@gmail.com",
    license="MIT License",
    packages=["taxumap"],
    install_requires=[
        "matplotlib",
        "pandas",
        "seaborn",
        "numpy",
        "scipy",
        "numba",
        "umap-learn",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    scripts=['run_taxumap.py'],
    zip_safe=False,
)
