# TaxUMAP

Visualize structure in large microbiome datasets. Implements a microbiome-research-focused extension of the Uniform Manifold Approximation and Projection (UMAP) by calculating microbiome sample-by-sample distances at different taxonomic aggregations of taxon abundances, and allowing for custom weighting of aggregates.

## Installation

> *Notice:* TaxUMAP will be made available on both PyPi and Bioconda for installation via pip and conda. But until then, please use `pip install -e .` as described below to install in developer mode.


```
git clone https://github.com/jsevo/taxumap.git
pip install -e .
```

## Quickstart: Notebook Example

* An interactive Jupyter notebook file is provided in `examples/taxumap_example.ipynb`. ***This demonstrates TaxUMAP on an example dataset we provide.*** See [below](#example_data) for more information.

* Additionally, two other notebooks are provided to cover both adjusting the TaxUMAP default `agg_levels` and `weights` parameters (`examples/adjusting_taxumap_parameters.ipynb`), as well as cleaning the taxonomy table (`examples/cleaning_taxonomy_table.ipynb`).


## Data required
Two tables are required: the **microbiota data** and a **taxonomy table**.

The ***microbiota data file*** (e.g., `examples/example_data/microbiota_table.csv`) must have a column with sample indices labeled 'index_column'. The remaining columns are expected to be amplicon sequence variant (ASV) labels or operational taxonomic unit (OTU) labels, i.e., the lowest level of taxonomy:

| index_column | ASV1 | ASV2 | ... | ASV500 |
| :--- | :---: | :---: | :---: | :---: |
|'sample1'| 0.5| 0.4| ... | 0.1 |
|'sample2'|0.2| 0.6| ... | 0.2 |
| ... | ... | ... | ... | ... |
|'sample3'|0.1| 0.4| ... | 0.5 |


The ***taxonomy table*** (e.g., `examples/example_data/taxonomy.csv`) is expected to resolve higher taxonomic groups for each ASV/OTU in the microbiota table. The **index** of the taxonomy table should be **ASV/OTU labels**, while the **columns** of the taxonomy table should be **higher taxonomic categories** (e.g., kingdom, phylum, etc.). The columns must be ordered from left to right in decreasing taxonomic hierarchy, e.g.:


| ASV | Kingdom    | Phylum       | ...   | Genus    | Species |
| :---: | :---:       | :---:        | :---: | :---:  |:---:  |
| 'ASV1' | 'Bacteria' | 'Firmicutes' | ...   | 'Staphylococcus' | 'aureus' |
| 'ASV2' | 'Bacteria' | 'Bacillota' | ...   | '[Ruminococcus]' | 'gnavus' |
| ... | ... | ... | ...   | ... | ... |
| 'ASV500' | 'Bacteria' | 'Verrucomicrobia' | ...  | 'Akkermansia' | 'muciniphila' |

In the above tables, the ``''`` designates strings. **Any UNKNOWN taxonomic levels (e.g., 'unknown species') should be set to np.nan or the string 'nan'.** For more information on how to properly resolve unknown taxonomic levels for TaxUMAP, **please see the notebook `examples/cleaning_taxonomy_table.ipynb`**. Finally, the taxonomy table should be *monophyletic*.


---

## Usage

### Command line:

```bash
run_taxumap.py -t examples/example_data/taxonomy.csv -m examples/example_data/microbiota_table.csv
```
The embedding will be saved in the current working folder, or to a location with the `-o path/to/folder/` flag. Additionally, for best results, the flag `-n` should be folllowed by the number of unique patients in your dataset (see **Optional** flag information below for more details).


### Python:
```python
from taxumap.taxumap_base import Taxumap

##### Initialize Taxumap object #####

# From file
tu = Taxumap(taxonomy='path/to/taxonomy.csv',
            microbiome_data='path/to/microbiota_table.csv')


##### Run the transformation and look at the results #####

# Transform the data (an inplace function)
tu.transform_self()

# Raw embedding dataframe
tu.df_embedding

# "Which taxon dominate each sample?" dataframe
tu.df_dominant_taxon

# Visualize the embedding (will save to present working directory as "taxumap_scatterplot.pdf")
tu.scatter(save=True)

# Save the embedding (will save to present working directory as "taxumap_embedding.csv" if no parameter passed)
tu.save_embedding('path/to/embedding.csv')

```

---

## Flags for Command-line Interface (CLI)

### Required

* `-t` or `--taxonomy`: filepath to your `taxonomy.csv` file
* `-m` or `--microbiota`: filepath to your `microbiota_table.csv` file

### Optional, but recommended

* `-a` or `--agg_levels`: Which taxonomic levels to aggregate, in the form of a `/`-delimined string (e.g. `Phylum/Family`, `Family`, `Phylum/Order/Genus`). Defaults to `Phylum/Family`.
* `-w` or `--weights`: Weights to give to each taxonomic level defined in `--agg_levels`, in the form of a `/`-delimined string (e.g. `5/6`, `0.5/2`, `6/2/1`, etc). Defaults to 1 for each specified aggregation level.
* `-n` or `--neigh`: Change the `n_neighbors` parameter passed to the UMAP algorithm. See documentation [here](https://umap-learn.readthedocs.io/en/latest/parameters.html?highlight=n_neighbors#n-neighbors).


### Optional, change default behavior

* `-o` or `--outdir`: Where to save embedding. Defaults to the present working directory.
* `-s` or `--save`: Set to False to not save the embedding. Defaults to True.
* `-b` or `--min_dist`: Change the `min_dist` parameter passed to the UMAP algorithm. See documentation [here](https://umap-learn.readthedocs.io/en/latest/parameters.html?highlight=min_dist#min-dist).


---

## Example notebook (with example data) <a name="example_data"></a>

A dataset provided by Olin et al. can be used to try out  features of TaxUMAP and how to format new data properly.

* [Link to original publication](https://pubmed.ncbi.nlm.nih.gov/30142345/)
* [Link to the dataset](http://dx.doi.org/10.17632/ynhdrcxtcc.1)

Publication
> Olin A, Henckel E, Chen Y, et al. Stereotypic Immune System Development in Newborn Children. Cell. 2018;174(5):1277-1292.e14. doi:10.1016/j.cell.2018.06.045

Dataset
> Olin, Axel (2018), “Stereotypic Immune System Development in Newborn Children”, Mendeley Data, v1

For convenience, we are providing in the `taxumap/examples/example_data` directory a pre-cleaned version of this dataset, as allowed under the `CC BY 4.0` license. An accompanying Jupyter Notebook is provided to demonstrated how the data was cleaned.

### generated summary

This code defines a class called `Taxumap`, which is used for running the taxUMAP algorithm. The class constructor takes several arguments, including `agg_levels`, `weights`, `microbiota_data`, `taxonomy`, `name`, and `random_state`. These arguments are used to initialize attributes of the `Taxumap` object.

Some of the methods of the `Taxumap` class include `transform_self()`, `scatter()`, `save_embedding()`, and `df_dominant_taxon()`. These methods are used to perform the taxUMAP transformation, generate a scatter plot of the transformed data, save the embedding to a file, and get the dominant taxonomic group based on the maximum abundance in each sample, respectively.

Overall, this code provides a framework for running the taxUMAP algorithm and visualizing the results.



## License

[MIT](https://choosealicense.com/licenses/mit/)
