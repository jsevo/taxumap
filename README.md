# taxumap

**Visualize structure in large microbiome data sets. Implements Uniform Manifold Approximation and Projection (UMAP) with phylogenetic hierarchy.**
## Installation

First, clone this repository. Install it. Copy the data to the expected location. Run.
## 
    git clone https://github.com/jsevo/taxumap.git
    cd taxumap
    pip install -e .
    cd taxumap
    cp /path/to/microbiota_table.csv data/
    cp /path/to/taxonomy.csv data/
    python3 taxumap

## Details
taxumap embedding of microbiota composition data. Data is expected to be compositional, i.e. each row sums to 1. Two tables are required: the microbiota data and a taxonomy table.

The microbiota data file (`data/microbiota_table.csv`) must have a column with sample indeces, labeled 'index_column'. The remaining columns are expected to be the lowest level taxa (OTU/ASV/...):

| index_column | ASV1 | ASV2 |
| :--- | :---: | :---: |
|'sample1'| 0.5| 0.5|
|'sample2'|0.2| 0.8|

The taxonomy table (`data/taxonomy.csv`) is expected to resolve higher taxonomic groups for the columns in the microbiota table. The columns of this table should contain taxonomic levels. They should be ordered from left to right in decreasing taxonomic hierarchy, e.g.

| kingdom    | phylum       | ...   | ASV    |
| :---       | :---:        | :---: | :---:  |
| 'Bacteria' | 'Firmicutes' | ...   | 'ASV1' |

The data is expected to be located in the `data/` folder. Results will be written to the `results/` folder

