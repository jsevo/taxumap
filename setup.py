from setuptools import setup

setup(
    name="taxumap",
    version="0.1",
    description="UMAP visualization for microbiota compositions with taxonomic structure.",
    url="http://github.com/jsevo/taxumap",
    author="Jonas Schluter",
    author_email="jonas.schluter+github@gmail.com",
    license="Apache License 2.0",
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
    zip_safe=False,
)
