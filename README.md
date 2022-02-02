# TaxUMAP

Visualize structure in large microbiome datasets. Implements a microbiome research-focused extension of the Uniform Manifold Approximation and Projection (UMAP) by calculating microbiome sample by sample distances at different taxonomic aggregations of taxon abundances, and allowing for custom weighting of aggregates.

## Installation

> *Notice:* TaxUMAP will be made available on both PyPi and Bioconda for installation via pip and conda. But until then, please use `pip install -e` as described below to install in developer mode.


## 
```
git clone https://github.com/jsevo/taxumap.git
pip install -e .
```

## Data required
Two tables are required: the microbiota data and a taxonomy table.

The ***microbiota data file*** (`microbiota_table.csv`) must have a column with sample indices labeled 'index_column'. The remaining columns are expected to be the lowest level taxa (OTU/ASV/...):

| index_column | ASV1 | ASV2 |
| :--- | :---: | :---: |
|'sample1'| 0.5| 0.5|
|'sample2'|0.2| 0.8|

You can see that the data is *compositional*, or that each row sums to 1. (This kind of table may also be referred to as a *relative abundance* or `rel_abundances` for that reason.) 

The ***taxonomy table*** (`taxonomy.csv`) is expected to resolve higher taxonomic groups for the columns in the microbiota table. The columns of this table should contain taxonomic levels. They should be ordered from left to right in decreasing taxonomic hierarchy, e.g.

| kingdom    | phylum       | ...   | ASV    |
| :---       | :---:        | :---: | :---:  |
| 'Bacteria' | 'Firmicutes' | ...   | 'ASV1' |

Unless designated by the `-t` and `-m` flags, the data is expected to be within the `data/` folder. Results are written to the `taxumap/results/` folder.

---

## Usage

### Command line:

```bash
python taxumap/run_taxumap.py -t taxonomy.csv -m microbiota_table.csv -n 15
```
Your embedding will be saved in the `taxumap/results/` folder. 


### Python:
```python
from taxumap.taxumap_base import Taxumap

##### Initialize Taxumap object #####

# From file
t = Taxumap(taxonomy='path/to/taxonomy.csv', 
            rel_abundances='path/to/microbiota_table.csv')

# OR #

# From local variable scope
## df_taxonomy :: pd.DataFrame
## rel_abundance :: pd.DataFrame

t = Taxumap(taxonomy=df_taxonomy, 
            rel_abundances=df_rel_abundances)


##### Run the transformation and look at the results #####

# Transform the data (an inplace function)
t.transform_self(neigh=13)

# Raw embedding dataframe
t.df_embedding

# "Which taxon dominate each sample?" dataframe
t.df_dominant_taxon

# Visualize the embedding
t.scatter()

# Save the embedding
t.save_embedding() 

```
---

## Flags for `run_taxumap.py`

### Required

* `-t` or `--taxonomy`: filepath to your `taxonomy.csv` file
* `-m` or `--microbiota`: filepath to your `microbiota_table.csv` file
* `-n` or `--neigh`: number of patients

### Optional, but recommended

* `-a` or `--agg_levels`: Which taxonomic levels to aggregate, in the form of a `/`-delimined string (e.g. `Phylum/Family`)
* `-w` or `--weights`: Weights to give to each taxonomic level defined in `--agg_levels`, in the form of a `/`-delimined string (e.g. `5/6`, `0.5/2`, `6/2/1`, etc). Defaults to 1 for each.
* 

### Optional, change default behavior

* `-o` or `--outdir`: Where to save embedding. Defaults to `taxumap/results`.
* `-v` or `--verbose`: Add flag to log INFO-level information.
* `-d` or `--debug`: Add flag to log DEBUG-level information.
* `-s` or `--save`: Set to False to not save the embedding. Defaults to True.
* `-b` or `--min_dist`: Change the `min_dist` parameter passed to the UMAP algorithm. See documentation [here](https://umap-learn.readthedocs.io/en/latest/parameters.html?highlight=min_dist#min-dist). 


## Documentation for Taxumap as a Package

See this link. (In progress)

---

## Details

## Roadmap

[1] As presented here, taxUMAP is used for visualizing large microbiota datasets. However, it has much broader applications to improving upon UMAP by informing the algorithm of the hierarchical structure of data.

We will be updating this package to include examples and adaptations needed for such use cases.

[2] We will be updating general user issues. Please submit a 

---

## Example data

A dataset provided by Axel Olin works well for those wanting to try out the features of taxUMAP or to better understand how to format your own data properly.

* [Link to original publication](https://pubmed.ncbi.nlm.nih.gov/30142345/)
* [Link to the dataset](http://dx.doi.org/10.17632/ynhdrcxtcc.1)

Publication
> Olin A, Henckel E, Chen Y, et al. Stereotypic Immune System Development in Newborn Children. Cell. 2018;174(5):1277-1292.e14. doi:10.1016/j.cell.2018.06.045

Dataset
> Olin, Axel (2018), “Stereotypic Immune System Development in Newborn Children”, Mendeley Data, v1

Solely for convenience, I am providing in the `taxumap/example_data` directory a pre-cleaned version of this dataset, as allowed under the `CC BY 4.0` license. I also provide a Jupyter Notebook to see how the data was cleaned.

## License

[MIT](https://choosealicense.com/licenses/mit/)
