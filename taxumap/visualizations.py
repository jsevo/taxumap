# Authors: Jonas Schluter <jonas.schluter@nyulangone.org>, Grant Hussey <grant.hussey@nyulangone.org>
# License: MIT

import os

#!/usr/bin/env python
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.spatial.distance as ssd
from sklearn.preprocessing import MinMaxScaler


def pretty_print(
    X,
    embedding,
    ivs,
    tax,
    usercolors=None,
    with_diversity_background=True,
    bgcolor="white",
):
    """Make a scatter plot of taxumap-embedded microbiota data. Samples are colored by their dominant Genus. The top 15 most abundant genera have a unique color, all other taxa are grey. Optionally, interpolate the diversity of samples in the local region of the embedded space and color the background accordingly, with darker shades indicating higher diversity."""

    import seaborn as sns
    import matplotlib.pyplot as plt

    from sklearn.preprocessing import LabelEncoder

    dominant_taxon_name = X.idxmax(axis=1)
    dominant_taxon_name = dominant_taxon_name.apply(lambda v: tax.loc[v]["Genus"])
    dominant_taxon = dominant_taxon_name.copy()

    top_15_taxa = dominant_taxon.value_counts().sort_values(ascending=False).head(15)
    top_15_taxa_labels = top_15_taxa.index

    dominant_taxon = dominant_taxon.apply(lambda v: v if v in top_15_taxa else "-1")

    lenc = LabelEncoder().fit(dominant_taxon)
    _t = lenc.transform(dominant_taxon)
    dominant_taxon = pd.Series(_t, index=X.index)

    from matplotlib import cm

    _ncolors = len(top_15_taxa)
    _ncolors = _ncolors if _ncolors <= 15 else 16

    cmap = cm.get_cmap("tab20c", _ncolors)
    embedding_colors = [cmap(x) if x != 0 else "whitesmoke" for x in dominant_taxon]
    embedding_labels = [
        lenc.inverse_transform([x])[0]
        if lenc.inverse_transform([x])[0] != "-1"
        else "other"
        for x in dominant_taxon
    ]

    ##set up figure
    plt.close("all")
    fig, ax = plt.subplots(figsize=(5, 5))
    if with_diversity_background:
        ## heatmap as background indicateing interpolated diversity in that region
        cmap = sns.dark_palette(color="white", as_cmap=True, reverse=True)
        from scipy.interpolate import griddata

        xmin, xmax = np.floor(min(embedding[:, 0])), np.ceil(max(embedding[:, 0]))
        ymin, ymax = np.floor(min(embedding[:, 1])), np.ceil(max(embedding[:, 1]))
        grid_x, grid_y = np.mgrid[xmin:xmax:15j, ymin:ymax:15j]
        grid_z1 = griddata(
            embedding, ivs, (grid_x, grid_y), method="linear", fill_value=np.nan
        )
        # plot heatmap
        ax.imshow(
            np.flipud(grid_z1.T),
            extent=(xmin, xmax, ymin, ymax),
            cmap=cmap,
            vmin=1,
            vmax=15,
            alpha=0.25,
        )
        ax.set_facecolor(bgcolor)

    # ax.set_aspect('equal',adjustable='box')
    ## taxumap scatter
    if usercolors is None:
        noncolored_idx = list(map(lambda x: x == "whitesmoke", embedding_colors))
        ax.scatter(
            embedding[noncolored_idx, 0],
            embedding[noncolored_idx, 1],
            c=np.array(embedding_colors)[noncolored_idx],
            s=3,
            alpha=1,
            marker="o",
            rasterized=True,
        )
        colored_idx = list(map(lambda x: x != "whitesmoke", embedding_colors))
        ax.scatter(
            embedding[colored_idx, 0],
            embedding[colored_idx, 1],
            c=np.array(embedding_colors)[colored_idx],
            s=3,
            alpha=1,
            marker="o",
            rasterized=True,
        )
        ax.scatter(
            embedding[:, 0],
            embedding[:, 1],
            facecolor="none",
            edgecolor="k",
            linewidth=0.1,
            s=3,
            alpha=1,
            marker="o",
            rasterized=True,
        )
        from matplotlib.lines import Line2D

        legend_elements = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                alpha=1,
                color=c,  # cmap(c),
                label=n,
            )
            for (n, c) in set(zip(embedding_labels, embedding_colors))
        ]
        ax.legend(handles=legend_elements, loc=(1.1, 0.01))

    else:
        dominant_asv = X.idxmax(axis=1)
        dominant_asv_rel = X.max(axis=1)
        embedding_colors = [
            "whitesmoke" if dominant_asv_rel[i] < 0.3 else usercolors.loc[x].values[0]
            for i, x in dominant_asv.iteritems()
        ]
        noncolored_idx = list(map(lambda x: x == "whitesmoke", embedding_colors))
        ax.scatter(
            embedding[noncolored_idx, 0],
            embedding[noncolored_idx, 1],
            c=np.array(embedding_colors)[noncolored_idx],
            s=3,
            alpha=1,
            linewidth=0.1,
            marker="o",
            rasterized=True,
        )
        colored_idx = list(map(lambda x: x != "whitesmoke", embedding_colors))
        ax.scatter(
            embedding[colored_idx, 0],
            embedding[colored_idx, 1],
            c=np.array(embedding_colors)[colored_idx],
            s=3,
            alpha=1,
            linewidth=0.1,
            marker="o",
            rasterized=True,
        )

        from matplotlib.lines import Line2D

        most_dominating = (
            dominant_asv.loc[dominant_asv_rel >= 0.3]
            .apply(lambda v: tax.loc[v]["Genus"])
            .value_counts()
            .sort_values(ascending=False)
            .head(30)
        )

        most_dominating_color = (
            dominant_asv.loc[dominant_asv_rel >= 0.3]
            .apply(lambda v: usercolors.loc[v].values[0])
            .value_counts()
            .sort_values(ascending=False)
            .head(30)
            .index
        )

        legend_names = np.array(
            list(
                map(
                    lambda v: tax.loc[v].Genus.values,
                    [
                        dominant_asv[dominant_asv_rel > 0.3]
                        .value_counts()
                        .head(30)
                        .index.to_list()
                    ],
                )
            )
        ).reshape(-1)
        legend_colors = np.array(
            list(
                map(
                    lambda v: usercolors.loc[v].values,
                    [
                        dominant_asv[dominant_asv_rel > 0.3]
                        .value_counts()
                        .head(30)
                        .index.to_list()
                    ],
                )
            )
        ).reshape(-1)

        legend_elements = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                alpha=1,
                color=c,  # cmap(c),
                label=n,
            )
            for (n, c) in set(zip(legend_names, legend_colors))
        ]
        ax.legend(handles=legend_elements, loc=(1.1, 0.01))

    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel("TaxUMAP-2")
    ax.set_xlabel("TaxUMAP-1")
    sns.despine()
    try:
        plt.gcf().savefig("results/projection.pdf", dpi=250, bbox_inches="tight")
        plt.axis("off")
        ax.legend().remove()
        plt.gcf().savefig(
            "results/no_axes_projection.png", dpi=250,
        )
        return fig

    except FileNotFoundError as fe:
        print(fe)
        if not os.path.isdir("results"):
            print()
            print("No results folder in current working directory.")
            print("Will create one.")
            os.mkdir("results")
    except:
        print("Not saving figure. Unknown error.")
        pass


def scatter_plot_with_colors(
    embedding,
    colors,
    legend_items=[],
    bgcolor="white",
    scatter_alpha=1.0,
    s=1,
    fig=None,
    ax=None,
    ax2=None,
    figsize=(5, 5),
    fontsize=8,
    marker="o",
):
    ##set up figure
    import matplotlib.pyplot as plt
    import seaborn as sns
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    if fig is None:
        plt.close("all")
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(10, 1)
        ax = fig.add_subplot(gs[:-1, :])
        # ax2 = inset_axes(ax, width="20%", height="10%", loc=5)
        ax2 = fig.add_subplot(gs[-1, :])

        ax.scatter(
            embedding[:, 0],
            embedding[:, 1],
            edgecolors="none",
            facecolors=colors,
            alpha=scatter_alpha,
            s=s,
            marker=marker,
            zorder=0,
        )
    else:

        ax.scatter(
            embedding[:, 0],
            embedding[:, 1],
            edgecolors="none",
            facecolors=colors,
            alpha=scatter_alpha,
            s=s,
            marker=marker,
        )

    ax.axis("equal")
    ax.set_ylabel("TaxUMAP-2", fontsize=fontsize)
    ax.set_xlabel("TaxUMAP-1", fontsize=fontsize)
    ax.set_facecolor(bgcolor)
    # ax.secondary_xaxis('top', functions=(deg2rad, rad2deg))

    #    ax.spines["top"].set_position(ax.spines['bottom'].get_position())
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    sns.despine(trim=True, top=False, bottom=False, ax=ax, offset=5)

    ax.spines["bottom"].set_visible(False)

    ax2.set_yticks([])

    return (fig, ax, ax2)
