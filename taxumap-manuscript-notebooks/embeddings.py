from sklearn.manifold import TSNE 
from sklearn.decomposition import PCA, KernelPCA
from umap import UMAP
from sklearn.preprocessing import MinMaxScaler

RUNEMBEDDINGS = True 
if RUNEMBEDDINGS:
    #simple PCA
    pcaembedding = PCA(n_components=2).fit_transform(XASV.fillna(0))
    
    #base embedding (kernel pca)
    kernelpcaembedding = KernelPCA(n_components=2).fit_transform(XASV.fillna(0))
    
    # non-phylo umap
    embedding_non_phylo_unscaled = UMAP(n_neighbors=120,min_dist=0.2, metric="manhattan").fit_transform(XASV)
    
    
    # embedding_non_phylo_scaled = UMAP(n_neighbors=120,min_dist=0.2, metric="manhattan").fit_transform(MinMaxScaler().fit_transform(XASV))


RUNTAXUMAPS = True 
if RUNTAXUMAPS: 
    from taxumap.taxumap import taxumap
    agg_levels = ["Phylum", "Family"]
    withscaling = False # do not scale the columns of X
    distanceperlevel = False # do not calculate a separate distance matrix at each phylogenetic level because we are using the manhattan distance 
    distancemetric = "manhattan"
    printfigure=False
    printwithdiversity=False #dont plot the average diversity in the background of the scatter plot
    X_in = XASV
    tax = taxonomy
    withusercolors=taxonomy_meta[["HexColor"]]


    TAXUMAP, X_embedded, taxumap_Xscaled, taxumap_X = taxumap(agg_levels,
                     withscaling,
                     distanceperlevel,
                     distancemetric,
                     printfigure,
                     printwithdiversity,
                     X_in,
                     tax,
                     withusercolors,
                     debug=True, #return tables
                     save_embedding=True#save xy coordinates
#                                                              );
    
    TAXUMAP_alllevels, X_embedded_alllevels, taxumap_Xscaled_alllevels, taxumap_X_alllevels = taxumap(["Phylum", "Class", "Order", "Family", "Genus"],
                     withscaling,
                     distanceperlevel,
                     distancemetric,
                     printfigure,
                     printwithdiversity,
                     X_in,
                     tax,
                     withusercolors,
                     debug=True, #return tables
                     save_embedding=False #save xy coordinates
                                                             );

#     TAXUMAPSCALED, X_embedded_scaled, taxumap_Xscaled_scaled, taxumap_X_scaled = taxumap(
#                      agg_levels,
#                      True,
#                      False,
#                      "euclidean",
#                      printfigure,
#                      printwithdiversity,
#                      X_in,
#                      tax,
#                      withusercolors,
#                      debug=True, #return tables
#                      save_embedding=True#save xy coordinates
#                                                              );

#     TAXUMAPSCALEDeuclidean, X_embedded_scaledeuclidean, taxumap_Xscaled_scaledeuclidean, taxumap_X_scaledeuclidean = taxumap(
#                      agg_levels,
#                      True,
#                      False,
#                      "euclidean",
#                      printfigure,
#                      printwithdiversity,
#                      X_in,
#                      tax,
#                      withusercolors,
#                      debug=True, #return tables
#                      save_embedding=True#save xy coordinates
#                                                              );
LOADPCoAS = False
if LOADPCoAS:
    pcoa_embedding_unweighted_unifrac = PCA(n_components=2).fit_transform(unweighted_unifrac.set_index("SampleID"))
    #Weighted Unifrac
    pcoa_embedding_weighted_unifrac = PCA(n_components=2).fit_transform(weighted_unifrac.set_index("SampleID"))

    
del unweighted_unifrac
del weighted_unifrac
#del TAXUMAPSCALED, taxumap_Xscaled_scaled, taxumap_X_scaled
#del TAXUMAPSCALEDeuclidean, taxumap_Xscaled_scaledeuclidean, taxumap_X_scaledeuclidean
del TAXUMAP_alllevels, taxumap_Xscaled_alllevels, taxumap_X_alllevels

write_now=False
if write_now:
    for (em,n) in zip(
        [pcaembedding,
         pcoa_embedding_unweighted_unifract[:,0:2], 
         pcoa_embedding_weighted_unifract, 
         embedding_non_phylo_unscaled,
         X_embedded_alllevels.values,
         X_embedded.values],
        ["pcaembedding",
        "pcoa_unweighted_unifrac_embedding", 
         "pcoa_weighted_unifrac_embedding",
        "embedding_nontax_umap_unscaled",
        "taxumap_alllevels",
        "current_taxumap_embedding"]):
        pd.DataFrame(em, index=XASV.index).to_csv("results/%s.csv"%n)