from setuptools import setup

setup(
    name='phylo-umap',
    version='0.1',
    description='UMAP visualization for microbiota compositions with phylogenetic structure.',
    url='http://github.com/jsevo/phylo-umap',
    author='Jonas Schluter',
    author_email='jonas.schluter+github@gmail.com',
    license='Apache License 2.0',
    packages=['phylo-umap'],
    install_requires=[
        'matplotlib',
        'pandas',
        'seaborn',
        'numpy',
        'scipy',
        'numba',
        'umap-learn',
    ],
    zip_safe=False)
