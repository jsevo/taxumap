#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import warnings
import scipy.spatial.distance as ssd
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
        tax.loc[
            ix,
            c] = tax[c].astype(str) + '____' + tax[taxlevels[-1]].astype(str)
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


def parse_microbiome_data(fp, idx_col='index_column'):
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
        sys.exit(2)
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


def parse_taxonomy_data(fp):
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
        sys.exit(2)
    if f.is_file():
        try:
            tax = pd.read_csv(fp, )
            print(
                "Reading taxonomy table. Assuming columns are ordered by phylogeny with in descending order of hierarchy."
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


def parse_asvcolor_data(fp):
    """Load the taxonomy data."""
    from pathlib import Path
    fp = Path(fp)
    try:
        f = fp.resolve(strict=True)
    except FileNotFoundError as fe:
        print("{0}".format(fe))
        print(
            "The color per ASV table should be located in the data/ subfolder and named asvcolors.csv"
        )
        sys.exit(2)
    if f.is_file():
        try:
            taxcolors = pd.read_csv(fp)
            print(
                "Reading color per ASV table."
            )
            try:
                assert taxcolors.columns[[0,1]].to_list()== ["ASV", "HexColor"]
            except AssertionError:
                print( 'Column names should be:  ["ASV", "HexColor"]. Choosing colors automatically.')
                return()

            taxcolors = taxcolors.set_index('ASV')
            if np.any(taxcolors.isna()):
                warnings.warn(
                    "Missing values (NaN) found for some taxcolors. Filling with 'grey'"
                )
                taxcolors = taxcolors.fillna('grey')
            return (taxcolors)
        except ValueError as ve:
            print("{0}".format(ve))
            print(
                "Please make sure the taxcolors has columns labeled ['ASV','HexColor'], and contain as values the ASV labels as strings and valid hex color stings"
            )
        except:
            print(
                "Please make sure the taxcolors has columns labeled ['ASV','HexColor'], and contain as values the ASV labels as strings and valid hex color stings"
            )


def phylogenetic_aggregation(X,
                             tax,
                             agg_levels,
                             distanceperlevel=False,
                             distancemetric='braycurtis'):
    """
    Expand the original microbiota sample data by aggregating taxa at different higher phylogenetic levels.

    Parameters
    ----------

    X : pandas.DataFrame, microbiota data table at lowers taxonomic level

    tax : pandas.DataFrame, taxonomy table

    agg_level : list; which phylogenetic levels to include. Defaults to the second highest (usually phylum) and third lowest (usually family).

    distanceperlevel: bool, should sample-to-sample distances be calculated per tax level

    Returns 
    ----------

    X : Expanded microbiota table.

    """
    import numpy as np
    if agg_levels == None:
        try:
            assert len(
                tax.columns
            ) > 3, "the taxonomy table has very few columns. Cannot aggregate phylogenetic levels. Reverting to regular UMAP"
            agg_levels = np.array(tax.columns)[[2, -2]]
            print(agg_levels)
            if agg_levels[0] == agg_levels[1]:
                # in case the taxonomy table is weird and has few columns
                agg_levels = agg_levels[0]
            _X = X.copy()
            if distanceperlevel:
                Xdist = ssd.cdist(_X, _X)
                Xdist = pd.DataFrame(Xdist, index=_X.index, columns=_X.index)
                for l in agg_levels:
                    print("aggregating on %s" % l)
                    Xagg = _aggregate_at_phylolevel(_X, tax, l)
                    if distanceperlevel:
                        Xagg = ssd.cdist(Xagg, Xagg, distancemetric)
                        Xagg = pd.DataFrame(Xagg,
                                            index=_X.index,
                                            columns=_X.index)
                        Xdist = Xdist + Xagg
                X = Xdist
            else:
                for l in agg_levels:
                    print("aggregating on %s" % l)
                    Xagg = _aggregate_at_phylolevel(_X, tax, l)
                    X = X.join(Xagg, lsuffix="_r")
                assert np.allclose(
                    _X.sum(axis=1),
                    X.sum(axis=1) / (len(agg_levels) + 1)
                ), "During aggregation, the sum of relative abundances is not equal to %d-times the original relative abundances. This would have been expected due to aggregating and joining" % (
                    (len(agg_levels) + 1))
            return (X)
        except AssertionError as ae:
            print("{0}".format(ae))


def pretty_print(X,
                 embedding,
                 ivs,
                 tax,
                 usercolors=None,
                 with_diversity_background=True,
                 bgcolor='white'):
    """Make a scatter plot of phyloUMAP-embedded microbiota data. Samples are colored by their dominant Genus. The top 15 most abundant genera have a unique color, all other taxa are grey. Optionally, interpolate the diversity of samples in the local region of the embedded space and color the background accordingly, with darker shades indicating higher diversity."""

    import seaborn as sns
    import matplotlib.pyplot as plt

    from sklearn.preprocessing import LabelEncoder
    dominant_taxon_name = X.idxmax(axis=1)
    dominant_taxon_name = dominant_taxon_name.apply(
        lambda v: tax.loc[v]['Genus'])
    dominant_taxon = dominant_taxon_name.copy()

    top_15_taxa = dominant_taxon.value_counts().sort_values(
        ascending=False).head(15)
    top_15_taxa_labels = top_15_taxa.index

    dominant_taxon = dominant_taxon.apply(lambda v: v
                                          if v in top_15_taxa else '-1')

    lenc = LabelEncoder().fit(dominant_taxon)
    _t = lenc.transform(dominant_taxon)
    dominant_taxon = pd.Series(_t, index=X.index)

    from matplotlib import cm
    _ncolors = len(top_15_taxa)
    _ncolors = _ncolors if _ncolors <= 15 else 16

    cmap = cm.get_cmap('tab20c', _ncolors)
    embedding_colors = [
        cmap(x) if x != 0 else 'whitesmoke' for x in dominant_taxon
    ]
    embedding_labels = [
        lenc.inverse_transform([x])[0]
        if lenc.inverse_transform([x])[0] != '-1' else 'other'
        for x in dominant_taxon
    ]

    ##set up figure
    plt.close('all')
    fig, ax = plt.subplots(figsize=(5, 5))
    if with_diversity_background:
        ## heatmap as background indicateing interpolated diversity in that region
        cmap = sns.dark_palette(color='white', as_cmap=True, reverse=True)
        from scipy.interpolate import griddata
        xmin, xmax = np.floor(min(embedding[:,
                                            0])), np.ceil(max(embedding[:, 0]))
        ymin, ymax = np.floor(min(embedding[:,
                                            1])), np.ceil(max(embedding[:, 1]))
        grid_x, grid_y = np.mgrid[xmin:xmax:15j, ymin:ymax:15j]
        grid_z1 = griddata(embedding,
                           ivs, (grid_x, grid_y),
                           method='linear',
                           fill_value=np.nan)
        # plot heatmap
        ax.imshow(np.flipud(grid_z1.T),
                  extent=(xmin, xmax, ymin, ymax),
                  cmap=cmap,
                  vmin=1,
                  vmax=15,
                  alpha=0.25)
        ax.set_facecolor(bgcolor)

    #ax.set_aspect('equal',adjustable='box')
    ## phyloUMAP scatter
    if usercolors is None:
        noncolored_idx = list(map(lambda x: x == 'whitesmoke', embedding_colors))
        ax.scatter(embedding[noncolored_idx, 0],
                   embedding[noncolored_idx, 1],
                   c=np.array(embedding_colors)[noncolored_idx],
                   s=3,
                   alpha=1,
                   marker='o',
                   rasterized=True)
        colored_idx = list(map(lambda x: x != 'whitesmoke', embedding_colors))
        ax.scatter(embedding[colored_idx, 0],
                   embedding[colored_idx, 1],
                   c=np.array(embedding_colors)[colored_idx],
                   s=3,
                   alpha=1,
                   marker='o',
                   rasterized=True)
        ax.scatter(embedding[:, 0],
                   embedding[:, 1],
                   facecolor='none',
                   edgecolor='k',
                   linewidth=.1,
                   s=3,
                   alpha=1,
                   marker='o',
                   rasterized=True)
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D(
                [0],
                [0],
                marker='o',
                linestyle='',
                alpha=1,
                color=c,  #cmap(c),
                label=n) for (n, c) in set(zip(embedding_labels, embedding_colors))
        ]
        ax.legend(handles=legend_elements, loc=(1.1, .01))


    else:
        dominant_asv = X.idxmax(axis=1)
        dominant_asv_rel = X.max(axis=1)
        embedding_colors = ['whitesmoke' if dominant_asv_rel[i]<0.3 else usercolors.loc[x].values[0] for i,x in dominant_asv.iteritems()]
        noncolored_idx = list(map(lambda x: x == 'whitesmoke', embedding_colors))
        ax.scatter(embedding[noncolored_idx, 0],
                   embedding[noncolored_idx, 1],
                   c=np.array(embedding_colors)[noncolored_idx],
                   s=3,
                   alpha=1,
                   linewidth=0.1,
                   marker='o',
                   rasterized=True)
        colored_idx = list(map(lambda x: x != 'whitesmoke', embedding_colors))
        ax.scatter(embedding[colored_idx, 0],
                   embedding[colored_idx, 1],
                   c=np.array(embedding_colors)[colored_idx],
                   s=3,
                   alpha=1,
                   linewidth=0.1,
                   marker='o',
                   rasterized=True)
       
        from matplotlib.lines import Line2D
        most_dominating = dominant_asv.loc[dominant_asv_rel>=0.3].apply(lambda v: tax.loc[v]['Genus']).value_counts().sort_values(ascending=False).head(30)

        most_dominating_color = dominant_asv.loc[dominant_asv_rel>=0.3].apply(lambda v: usercolors.loc[v].values[0]).value_counts().sort_values(ascending=False).head(30).index

        legend_names = np.array(list(map(lambda v: tax.loc[v].Genus.values, [dominant_asv[dominant_asv_rel>0.3].value_counts().head(30).index.to_list()]))).reshape(-1)
        legend_colors = np.array(list(map(lambda v: usercolors.loc[v].values, [dominant_asv[dominant_asv_rel>0.3].value_counts().head(30).index.to_list()]))).reshape(-1)


        legend_elements = [
            Line2D(
                [0],
                [0],
                marker='o',
                linestyle='',
                alpha=1,
                color=c,  #cmap(c),
                label=n) for (n, c) in set(zip(legend_names,
                                               legend_colors))]
        ax.legend(handles=legend_elements, loc=(1.1, .01))




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


def phylo_umap(agg_levels,
               withscaling,
               distanceperlevel,
               distancemetric,
               print_figure,
               print_with_diversity,
               loadembedding=False,
               withusercolors=False,
               bgcolor='white',
               fpx='data/microbiota_table.csv',
               fpt='data/taxonomy.csv',
               fpc='data/asvcolors.csv'):
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
    agg_level : list; which phylogenetic levels to include. Defaults to the second highest (usually phylum) and third lowest (usually family).

    print_figure : bool, flag if pretty scatter plot of phylo-UMAP embedded data should be produced.

    print_with_diversity : bool, flag if the scatter plot should have a background indicating alpha-diversity of samples in the region (inverse Simpson).

    fpx : File path to microbiota compositional data table in csv table format. Defaults to ./data/microbiota_table.csv.

    fpt : File path to the taxonomy table in csv table format. Defaults to ./data/taxonomy.csv. 

    fpc : (optional) File path to the color per ASV table in csv table format. Defaults to ./data/asvcolors.csv. 

    Returns
    ----------

    phyloUMAP : UMAP object fit to data

    X_embedded : pandas.DataFrame of 2-D coordinates with indeces from the microbiota_table.l

    """
    X = parse_microbiome_data(fpx)
    # _cols = X.sum(axis=0).sort_values().tail(100).index
    # X = X[_cols]
    tax = parse_taxonomy_data(fpt)
    if withusercolors:
        asv_colors = parse_asvcolor_data(fpc)
    else:
        asv_colors = None
    # aggregate phylogenetic levels
    Xagg = phylogenetic_aggregation(X,
                                    tax,
                                    agg_levels,
                                    distanceperlevel=distanceperlevel)

    # scale
    if withscaling:
        X_scaled = _scale(Xagg)
    else:
        print('not scaling')
        X_scaled = Xagg

    from umap import UMAP
    # neigh = 5 if int(X_scaled.shape[0] / 100) < 5 else int(X_scaled.shape[0] /
    #                                                        100)
    neigh = 120
    min_dist = 0.2

    if loadembedding:
        X_embedded = pd.read_csv('results/embedded.csv', index_col = "index_column" )
        embedding = X_embedded.values
        phyloUMAP = {}

    else:
        if withscaling:
            phyloUMAP = UMAP(n_neighbors=neigh,
                            min_dist=min_dist,
                            metric='manhattan').fit(X_scaled)

        elif distanceperlevel:
            phyloUMAP = UMAP(n_neighbors=neigh,
                            min_dist=min_dist,
                            metric='precomputed').fit(X_scaled)
        else:
            phyloUMAP = UMAP(n_neighbors=neigh,
                            min_dist=min_dist,
                            metric=distancemetric).fit(X_scaled)

        embedding = phyloUMAP.transform(X_scaled)
        X_embedded = pd.DataFrame(embedding,
                                index=X.index,
                                columns=['phyUMAP-1', 'phyUMAP-2'])
        X_embedded.to_csv('results/embedded.csv')
    if print_figure:
        # diversity per sample
        # if not np.all(X.sum(axis=1) > 0):
        #     warnings.warn("some rows have zero sum. adding tiny value (1e-16)")
        #     X = X + 1e-16
        if print_with_diversity:
            ivs = X.apply(lambda r: 1 / np.sum(r**2), axis=1)
        else:
            ivs = X.iloc[:, 1]
        pretty_print(X,
                     embedding,
                     ivs,
                     tax,
                     with_diversity_background=print_with_diversity,
                     usercolors = asv_colors,
                     bgcolor=bgcolor)

    return (phyloUMAP, X_embedded)


if __name__ == "__main__":
    import sys, getopt
    """
    use with options:

    -e OR WITH --scalingpertaxlevel
    -p OR WITH --print_figure
    -s OR WITH --with_diversity_background
    -a OR --asvcolors 
    -c OR --scatterbgcolor : a valid matplotlib color for the background of the scatter plot, default: white
    -l OR --loadembedding : load already calculated embedding from result/embedding.csv
    """
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "epsalc:", [
            "scalingpertaxlevel", "print_figure", "with_diversity_background",
            "asvcolors","loadembedding",
            "scatterbgcolor="
        ])
    except getopt.GetoptError:
        sys.exit("unknown options")

    agg_levels = None
    withscaling = False
    distanceperlevel =False  
    distancemetric = 'braycurtis'
    print_figure = False
    distanceperlevel = False
    withusercolors=False
    with_diversity_background = False
    loadembedding = False

    for opt, arg in opts:
        if opt in ("-e", "--scalingpertaxlevel"):
            withscaling = True
        elif opt in ("-p", "--print_figure"):
            print_figure = True
        elif opt in ("-s", "--with_diversity_background"):
            with_diversity_background = True
        elif opt in ("-c", "--scatterbgcolor"):
            bgcolor = str(arg)
        elif opt in ("-a", "--asvcolors"):
            withusercolors = True 
        elif opt in ("-l", "--loadembedding"):
            loadembedding = True 



    _ = phylo_umap(agg_levels,
                   withscaling,
                   distanceperlevel,
                   distancemetric,
                   print_figure,
                   withusercolors=withusercolors,
                   print_with_diversity=with_diversity_background,
                   loadembedding=loadembedding,
                   bgcolor=bgcolor)
