#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import warnings
from sklearn.preprocessing import MinMaxScaler


def _fill_taxonomy_table(tax):
    """
    Helper function that fills nan's in a taxonomy table. Such gaps are filled 'from the left' with the next higher non-nan taxonomy level and the lowest level (e.g. OTU# or ASV#) appended.
    """
    taxlevels = list(tax.columns[1::])
    root_level = tax.columns[0]
    tax[root_level] = tax[root_level].fillna('unknown_%s' % root_level)
    for i, level in enumerate(taxlevels):
        _missing_l = tax[level].isna()
        tax.loc[_missing_l, level] = [
            'unknown_%s_of_' % level + str(x)
            for x in tax.loc[_missing_l][taxlevels[i - 1]]
        ]
    for i, (ix, c) in enumerate(tax.iteritems()):
        tax.loc[ix, c] = tax[c].astype(str)+'____'+tax[taxlevels[-1]].astype(str)
    tax = tax.applymap(lambda v: v if 'unknown' in v else v.split('____')[0])

    return (tax)


def _aggregate_at_phylolevel(X, tax, level):
    """Helper function. For a given taxonomic level, aggregate relative abundances by summing all members of corresponding taxon."""
    _X_agg = X.copy()
    _X_agg.columns = [tax.loc[x][level] for x in _X_agg.columns]
    _X_agg = _X_agg.groupby(_X_agg.columns, axis=1).sum()
    return (_X_agg)


def _scale(X, scaler=MinMaxScaler()):
    """Min max scaling of relative abundances to ensure that different taxonomic levels have comparable dynamic ranges."""
    scaler = MinMaxScaler()
    X_sum = X.sum()
    X_stats = X.apply(['mean', 'median', 'max']).T
    X_consider = X_stats.applymap(lambda v: v > 0.001).apply(np.any, axis=1)
    X_consider = X_consider[X_consider.values].index

    X_scaled = scaler.fit_transform(X[X_consider])
    return (X_scaled)


def parse_microbiome_data(fp='data/microbiota_table.csv',
                          idx_col='index_column'):
    """Load the microbiota data."""

    from pathlib import Path
    fp = Path(fp)
    try:
        f = fp.resolve(strict=True)
    except FileNotFoundError as fe:
        print("{0}".format(fe))
        print(
            "The microbiota composition table should be located in the data/ subfolder and named microbiota_table.csv"
        )
        pass
    if fp.is_file():
        try:
            X = pd.read_csv(
                fp,
                index_col=idx_col,
            )
            X = X.astype(np.float64)
            if not np.allclose(X.sum(axis=1), 1):
                warnings.warn("rows do not sum to 1. Is this intentional?")
            return (X.fillna(0))
        except ValueError as ve:
            print("{0}".format(ve))
            print(
                "Please make sure the microbiota_table has a column labeled 'index_column' which contains the sample IDs"
            )
        except:
            print(
                "An unknown error occurred during microbiota_table parsing. Please see the instructions for how to run phylo-umap."
            )


def parse_taxonomy_data(fp='data/taxonomy.csv'):
    """Load the taxonomy data."""
    from pathlib import Path
    fp = Path(fp)
    try:
        f = fp.resolve(strict=True)
    except FileNotFoundError as fe:
        print("{0}".format(fe))
        print(
            "The taxonomy table should be located in the data/ subfolder and named taxonomy.csv"
        )
        pass
    if f.is_file():
        try:
            tax = pd.read_csv(fp)
            print(
                "Reading taxonomy table. Assuming columns are ordered by phylogeny with in ascending order of hierarchy."
            )
            tax = tax.set_index(tax.columns[-1])
            if np.any(tax.isna()):
                warnings.warn(
                    "Missing values (NaN) found for some taxonomy levels, filling with higher taxonomic level names"
                )
                tax = _fill_taxonomy_table(tax)
            return (tax)
        except ValueError as ve:
            print("{0}".format(ve))
            print(
                "Please make sure the taxonomy has a column labeled 'index_column' which contains the sample IDs"
            )
        except:
            print(
                "An unknown error occurred during taxonomy table parsing. Please see the instructions for how to run phylo-umap."
            )


def phylogenetic_aggregation(X, tax, agg_levels):
    """
    Expand the original microbiota sample data by aggregating taxa at different higher phylogenetic levels.

    Parameters
    ----------

    X : pandas.DataFrame, microbiota data table at lowers taxonomic level

    tax : pandas.DataFrame, taxonomy table

    agg_level : list; which phylogenetic levels to include. Defaults to the second highest (usually phylum) and third lowest (usually family).

    Returns 
    ----------

    X : Expanded microbiota table.

    """
    if agg_levels == None:
        try:
            assert len(
                tax.columns
            ) > 3, "the taxonomy table has very few columns. Cannot aggregate phylogenetic levels. Reverting to regular UMAP"
            agg_levels = np.array(tax.columns)[[1, -3]]
            print(agg_levels)
            if agg_levels[0] == agg_levels[1]:
                # in case the taxonomy table is weird and has few columns
                agg_levels = agg_levels[0]
            _X = X.copy()
            for l in agg_levels:
                print("aggregating on %s" % l)
                Xagg = _aggregate_at_phylolevel(_X, tax, l)
                X = X.join(Xagg)
            assert np.allclose(
                _X.sum(axis=1), X.sum(axis=1) / (len(agg_levels) + 1)
            ), "During aggregation, the sum of relative abundances is not equal to %d-times the original relative abundances. This would have been expected due to aggregating and joining" % (
                (len(agg_levels) + 1))
            return (X)
        except AssertionError as ae:
            print("{0}".format(ae))


def pretty_print(X, embedding, ivs, with_diversity_background=True):

    """Make a scatter plot of phyloUMAP-embedded microbiota data. Samples are colored by their dominant taxon at lowest taxonomic levels. The top 15 most abundant taxa have a unique color, all other taxa are grey. Optionally, interpolate the diversity of samples in the local region of the embedded space and color the background accordingly, with darker shades indicating higher diversity."""

    import seaborn as sns
    import matplotlib.pyplot as plt

    from sklearn.preprocessing import LabelEncoder
    dominant_taxon_name = X.idxmax(axis=1)
    lenc = LabelEncoder().fit(dominant_taxon_name)
    _t = lenc.transform(dominant_taxon_name)
    dominant_taxon = pd.Series(_t, index=X.index)
    top_15_taxa = dominant_taxon.value_counts().sort_values(ascending=False).head(15)
    top_15_taxa_labels = dominant_taxon_name.value_counts().sort_values(ascending=False).head(15) 
    dominant_taxon = dominant_taxon.apply(
        lambda v: v if v in top_15_taxa else 15)


    from matplotlib import cm
    _ncolors = len(top_15_taxa)
    _ncolors = _ncolors if _ncolors <= 15 else 16

    cmap = cm.get_cmap('tab20c', _ncolors)
    embedding_colors = [cmap(x) for x in dominant_taxon]
    embedding_labels = [
        lenc.inverse_transform([x])[0] if x != 16 else 'other'
        for x in dominant_taxon
    ]

    ##set up figure
    plt.close('all')
    fig, ax = plt.subplots(figsize=(5, 5))
    with_diversity_background = True
    if with_diversity_background:
        ## heatmap as background indicateing interpolated diversity in that region
        cmap = sns.dark_palette(color='white', as_cmap=True, reverse=True)
        from scipy.interpolate import griddata
        xmin, xmax = np.floor(min(embedding[:, 0])), np.ceil(
            max(embedding[:, 0]))
        ymin, ymax = np.floor(min(embedding[:, 1])), np.ceil(
            max(embedding[:, 1]))
        grid_x, grid_y = np.mgrid[xmin:xmax:15j, ymin:ymax:15j]
        grid_z1 = griddata(
            embedding,
            ivs, (grid_x, grid_y),
            method='linear',
            fill_value=np.nan)
        # plot heatmap
        ax.imshow(
            np.flipud(grid_z1.T),
            extent=(xmin, xmax, ymin, ymax),
            cmap=cmap,
            vmin=1,
            vmax=15,
            alpha=0.25)

    #ax.set_aspect('equal',adjustable='box')
    ## tsne scatter
    ax.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c=embedding_colors,
        s=3,
        alpha=1,
        marker='o',
        rasterized=True)

    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0],
                              marker='o',
                              linestyle='',
                              alpha=1,
                              color=c,#cmap(c),
                              label=n) for (n, c)
                       in set(zip(embedding_labels, embedding_colors))]
    ax.legend(handles=legend_elements, loc=(1.1,.01))
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel('phyUMAP-2')
    ax.set_xlabel('phyUMAP-1')
    sns.despine()
    plt.gcf().savefig('results/projection.pdf', dpi=250, bbox_inches='tight')
    plt.axis('off')
    ax.legend().remove()
    plt.gcf().savefig(
        'results/no_axes_projection.png',
        dpi=250,
    )


def phylo_umap(fpx='data/microbiota_table.csv',
               fpt='data/taxonomy.csv',
               agg_levels=None,
               print_figure=False,
               print_with_diversity=False):
    """
    phylo-UMAP embedding of microbiota composition data. Data is expected to be compositional, i.e. each row sums to 1. Two tables are required: the microbiota data and a taxonomy table.

    The microbiota data file ('data/microbiota_table.csv') must have a column with sample indeces, labeled 'index_column'. The remaining columns are expected to be the lowest level taxa (OTU/ASV/...):

    index_column, ASV1, ASV2
    'sample1', 0.5, 0.5
    'sample2', 0.2, 0.8

    The taxonomy table ('data/taxonomy.csv') is expected to resolve higher taxonomic groups for the columns in the microbiota table. The columns of this table should contain taxonomic levels. They should be ordered from left to right in decreasing taxonomic hierarchy, e.g.

    kingdom, phylum, ..., ASV
    'Bacteria', 'Firmicutes', ..., 'ASV1'

    The data is expected to be located in the data/ folder. Results will be written to the results/ folder

    Parameters
    ----------
    fpx : File path to microbiota compositional data table in csv table format. Defaults to ./data/microbiota_table.csv.

    fpx : File path to the taxonomy table in csv table format. Defaults to ./data/taxonomy.csv. 

    agg_level : list; which phylogenetic levels to include. Defaults to the second highest (usually phylum) and third lowest (usually family).

    print_figure : bool, flag if pretty scatter plot of phylo-UMAP embedded data should be produced.

    print_with_diversity : bool, flag if the scatter plot should have a background indicating alpha-diversity of samples in the region (inverse Simpson).


    Returns 
    ----------

    phyloUMAP : UMAP object fit to data

    X_embedded : pandas.DataFrame of 2-D coordinates with indeces from the microbiota_table.

    """
    X = parse_microbiome_data(fpx)
    tax = parse_taxonomy_data(fpt)
    # aggregate phylogenetic levels
    Xagg = phylogenetic_aggregation(X, tax, agg_levels)
    # scale
    X_scaled = _scale(Xagg)

    from umap import UMAP
    neigh = 5 if int(X_scaled.shape[0] / 100) < 5 else int(
        X_scaled.shape[0] / 100)

    phyloUMAP = UMAP(
        n_neighbors=neigh, min_dist=0.25, metric='euclidean').fit(X_scaled)

    embedding = phyloUMAP.transform(X_scaled)
    X_embedded = pd.DataFrame(
        embedding, index=X.index, columns=['phyUMAP-1', 'phyUMAP-2'])
    X_embedded.to_csv('results/embedded.csv')
    if print_figure:
        # diversity per sample
        if not np.all(X.sum(axis=1) > 0):
            warnings.warn("some rows have zero sum. adding tiny value (1e-16)")
            X = X + 1e-16
        if print_with_diversity:
            ivs = X.apply(lambda r: 1 / np.sum(r**2), axis=1)
        else:
            ivs = X.iloc[:, 1]
        pretty_print(
            X, embedding, ivs, with_diversity_background=print_with_diversity)

    return (phyloUMAP, X_embedded)


if __name__ == "__main__":
    import sys, getopt
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(
            argv, "hp:s:", ["print_figure=", "with_diversity_background="])
    except getopt.GetoptError:
        sys.exit(2)

    print_figure = False
    with_diversity_background = False

    for opt, arg in opts:
        if opt in ("-p", "print_figure"):
            print_figure = True if int(arg) == 1 else False
            print(opt, arg)
        elif opt in ("-s", "--with_diversity_background"):
            with_diversity_background = True if int(arg) == 1 else False
            print(opt, arg)

    _ = phylo_umap(
        print_figure=print_figure,
        print_with_diversity=with_diversity_background)
