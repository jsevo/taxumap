# Authors: Jonas Schluter <jonas.schluter@nyulangone.org>
# License: BSD 3 clause
__author__ = "Jonas Schluter"
__copyright__ = "Copyright 2020, MIT License"


### prepare tables
alldata = load_all_data()
counts = alldata["counts"]
taxonomy = alldata["taxonomy"]
# make generically labeled Genera unique
taxonomy.loc[taxonomy.Genus.str.contains("metagenome") & ~taxonomy.Genus.str.contains("ASV"), 'Genus'] = taxonomy.loc[taxonomy.Genus.str.contains("metagenome") & ~taxonomy.Genus.str.contains("ASV"), 'Genus'] + "___" + taxonomy.loc[taxonomy.Genus.str.contains("metagenome") & ~taxonomy.Genus.str.contains("ASV"), 'Genus'].index
taxonomy_meta = alldata["taxonomy_meta"]
samples = alldata["samples"]
qpcr = alldata["qpcr"]
qpcr = (
    qpcr.groupby("SampleID")
    .agg(lambda v: np.power(np.product(v), 1 / len(v)))
    .reset_index()
)
hct_meta = alldata["metadata"].copy()

asv_table = get_composition_at_taxlevel(counts, taxonomy.reset_index(), "ASV")
family_table = get_composition_at_taxlevel(counts, taxonomy.reset_index(), "Family")

# some samples recorded but not yet sequenced/sequencing failed. reindex
samples = samples.loc[samples.SampleID.isin(asv_table.index)].copy()
samples["Timepoint"] = samples["Timepoint"].astype(np.float)
hct_meta["DayOfTransplant"] = hct_meta["TimepointOfTransplant"].astype(np.float)
hct_meta["TimepointOfTransplant"] = hct_meta["DayOfTransplant"].astype(np.float)

asv_table_pp = asv_table.join(samples.assign(hctday = lambda r: r.DayRelativeToNearestHCT).set_index("SampleID"))
asv_table_pp = asv_table_pp.dropna()
asv_table_pp = asv_table_pp.sort_values(["PatientID", "hctday"])

X = pd.merge_asof(
    asv_table_pp.reset_index().sort_values("hctday"),
    hct_meta.sort_values("DayOfTransplant")[["DayOfTransplant", "PatientID",]],
    left_on="hctday",
    right_on="DayOfTransplant",
    by="PatientID",
    direction="nearest",
)
# select the ASV relative abundance columns
XASV = X.set_index("SampleID")[[x for x in X.columns if "ASV_" in x]]


########################################
scatter_meta_data = pd.read_csv(
    "/Users/schluj05/data/MSKCCDATA/meta_data_for_scatter_jitter_volatility_clean.csv",
    index_col="SampleID",
)
scatter_meta_data.index.name = "SampleID"

scatter_meta_data = scatter_meta_data.drop_duplicates()
samples['Day'] = samples['DayRelativeToNearestHCT']
scatter_meta_data['Day'] = scatter_meta_data['DayRelativeToNearestHCT']
########################################

_tax = taxonomy.loc[XASV.columns[0:100]]
_tax.index.name = "ASV"
# _tax = _tax.reset_index()
_tax = _tax[["Kingdom", "Phylum", "Genus"]]
lineages = {str(ix): r.to_list() for i, (ix, r) in enumerate(_tax.iterrows())}
tree = TreeNode.from_taxonomy(lineages.items())
tree_rooted = tree.root_at(tree)

for n in tree.traverse(include_self=False):
    n.length = 1.0
for n in tree_rooted.traverse(include_self=False):
    n.length = 1.0



load_distances = True   

if load_distances:

    unweightedunifrac_distances = pd.read_csv(
        "results/unweighted_unifrac_dist.csv", index_col="SampleID"
    )

    weightedunifrac_distances = pd.read_csv(
        "results/weighted_unifrac_dist.csv", index_col="SampleID"
    )


else:
    unweightedunifrac_distances = beta_diversity(
        "unweighted_unifrac",
        XASV.iloc[:, 0:100].values * 10000,
        XASV.index.to_list(),
        otu_ids=XASV.columns[0:100].to_list(),
        tree=tree_rooted,
        validate=True,
    )

    weightedunifrac_distances = beta_diversity(
        "weighted_unifrac",
        XASV.iloc[:, 0:100].values * 10000,
        XASV.index.to_list(),
        otu_ids=XASV.columns[0:100].to_list(),
        tree=tree_rooted,
        validate=True,
    )

    pd.DataFrame(
        unweightedunifrac_distances.data, index=XASV.index, columns=XASV.index
    ).to_csv("results/unweighted_unifrac_dist.csv")

    pd.DataFrame(
        weightedunifrac_distances.data, index=XASV.index, columns=XASV.index
    ).to_csv("results/weighted_unifrac_dist.csv")
