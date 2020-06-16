# Authors: Jonas Schluter <jonas.schluter@nyulangone.org>
# License: BSD 3 clause
__author__ = "Jonas Schluter"
__copyright__ = "Copyright 2020, MIT License"

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from matplotlib.collections import PatchCollection
import matplotlib.cm as cm
import numpy as np
import seaborn as sns
from hctmicrobiomemskcc.dataloading.dataloading import (
    load_microbiome_tables,
    load_metadata,
    load_all_data,
)
from hctmicrobiomemskcc.dataloading.hct_datawrangling import calculate_hct_day
from hctmicrobiomemskcc.tools.microbiotatools import (
    calculate_relative_counts,
    get_composition_at_taxlevel,
)
from hctmicrobiomemskcc.tools.tools import isfloat
import matplotlib
import pandas as pd

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

### prepare tables
alldata = load_all_data()
counts = alldata["counts"]
taxonomy = alldata["taxonomy"]
taxonomy_meta = alldata["taxonomy_meta"]
samples = alldata["samples"]
qpcr = alldata["qpcr"]
qpcr = (
    qpcr.groupby("SampleID")
    .agg(lambda v: np.power(np.product(v), 1 / len(v)))
    .reset_index()
)
hct_meta = alldata["metadata"]
hct_meta["DayOfTransplant"] = hct_meta["DayOfTransplant"].astype(np.float)

asv_table = get_composition_at_taxlevel(counts, taxonomy.reset_index(), "ASV")


# some samples recorded but not yet sequenced/sequencing failed. reindex
samples = samples.loc[samples.SampleID.isin(asv_table.index)].copy()
samples["Timepoint"] = samples["Timepoint"].astype(np.float)
hct_meta["DayOfTransplant"] = hct_meta["DayOfTransplant"].astype(np.float)
hct_day = calculate_hct_day(samples.set_index("SampleID"), hct_meta)

# ## quick check: some samples have no hct day / no pid. are they contained in the qpcr? No:
# no_pid_or_day = (
#     (genus_table_pp.isna().sum(axis=1))
#     .loc[(genus_table_pp.isna().sum(axis=1)) == 1]
#     .index
# )
# print(len(qpcr.reindex(no_pid_or_day).dropna()))

asv_table_pp = asv_table.join(samples.set_index("SampleID").join(hct_day))
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


