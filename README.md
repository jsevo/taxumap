# TaxUMAP

Visualize structure in large microbiome datasets. Implements a microbiome-research-focused extension of the Uniform Manifold Approximation and Projection (UMAP) by calculating microbiome sample-by-sample distances at different taxonomic aggregations of taxon abundances, and allowing for custom weighting of aggregates.

## Installation

> *Notice:* TaxUMAP will be made available on both PyPi and Bioconda for installation via pip and conda. But until then, please use `pip install -e .` as described below to install in developer mode.


```
git clone https://github.com/jsevo/taxumap.git
pip install -e .
```

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
| ... | ... | ... | ...   | ... |
| 'ASV500' | 'Bacteria' | 'Verrucomicrobia' | ...  | 'Akkermansia' | 'muciniphila' |

In the above tables, the ``''`` designates strings. **Any UNKNOWN taxonomic levels (e.g., 'unknown species') should be set to np.nan or the string 'nan'.** For more information on how to properly resolve unknown taxonomic levels for TaxUMAP, **please see the notebook `examples/cleaning_taxonomy_table.ipynb`**. Finally, the taxonomy table should be *monophyletic*.


---

## Quickstart

### Command line:

```bash
python taxumap/run_taxumap.py -t taxonomy.csv -m microbiota_table.csv
```
Your embedding will be saved in your current folder, or you can provide a location with the `-o path/to/folder/` flag. Additionally, for best results, the flag `-n` should be folllowed by the number of unique patients in your dataset (see flag information below for more details).


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

### Notebook:

* An interactive Jupyter notebook file is provided in `examples/taxumap_example.ipynb`. ***This demonstrates TaxUMAP on an example dataset we provide.*** See [below](#example_data) for more information.

* Additionally, two other notebooks are provided to cover both adjusting the TaxUMAP default `agg_levels` and `weights` parameters (`examples/adjusting_taxumap_parameters.ipynb`), as well as cleaning the taxonomy table (`examples/cleaning_taxonomy_table.ipynb`).


---

## Flags for `run_taxumap.py`

### Required

* `-t` or `--taxonomy`: filepath to your `taxonomy.csv` file
* `-m` or `--microbiota`: filepath to your `microbiota_table.csv` file

### Optional, but recommended

* `-n` or `--neigh`: number of unique microbiota hosts in datasafe, if applicable
* `-a` or `--agg_levels`: Which taxonomic levels to aggregate, in the form of a `/`-delimined string (e.g. `Phylum/Family`, `Family`, `Phylum/Order/Genus`). Defaults to `Phylum/Family`.
* `-w` or `--weights`: Weights to give to each taxonomic level defined in `--agg_levels`, in the form of a `/`-delimined string (e.g. `5/6`, `0.5/2`, `6/2/1`, etc). Defaults to 1 for each specified aggregation level.
*

### Optional, change default behavior

* `-o` or `--outdir`: Where to save embedding. Defaults to the present working directory.
* `-s` or `--save`: Set to False to not save the embedding. Defaults to True.
* `-b` or `--min_dist`: Change the `min_dist` parameter passed to the UMAP algorithm. See documentation [here](https://umap-learn.readthedocs.io/en/latest/parameters.html?highlight=min_dist#min-dist).


## Documentation

 * To-do: Explain weights
 * To-do: Explain agg_levels
 * These concepts are illustrated in the notebook found at `examples/adjusting_taxumap_parameters.ipynb`.

---

## Details

## Roadmap

[1] As presented here, taxUMAP is used for visualizing large microbiota datasets. However, it has much broader applications to improving upon UMAP by informing the algorithm of the hierarchical structure of data.

We will be updating this package to include examples and adaptations needed for such use cases.

[2] We will be updating general user issues. Please submit a GitHub issue for any comments, questions, concerns.

---

## Example notebook (with example data) <a name="example_data"></a>

A dataset provided by Olin et al. works well for those wanting to try out the features of taxUMAP or to better understand how to format your own data properly.

* [Link to original publication](https://pubmed.ncbi.nlm.nih.gov/30142345/)
* [Link to the dataset](http://dx.doi.org/10.17632/ynhdrcxtcc.1)

Publication
> Olin A, Henckel E, Chen Y, et al. Stereotypic Immune System Development in Newborn Children. Cell. 2018;174(5):1277-1292.e14. doi:10.1016/j.cell.2018.06.045

Dataset
> Olin, Axel (2018), “Stereotypic Immune System Development in Newborn Children”, Mendeley Data, v1

Solely for convenience, I am providing in the `taxumap/examples/example_data` directory a pre-cleaned version of this dataset, as allowed under the `CC BY 4.0` license. I also provide a Jupyter Notebook to see how the data was cleaned.



## License

[MIT](https://choosealicense.com/licenses/mit/)
