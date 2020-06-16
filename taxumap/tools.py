
import numpy as np
import pandas as pd
from hctmicrobiomemskcc.tools.microbiotatools import fill_taxonomy_table

# def fill_taxonomy_table(tax):
#     """
#     Helper function that fills nan's in a taxonomy table. Such gaps are filled 'from the left' with the next higher non-nan taxonomy level and the lowest level (e.g. OTU# or ASV#) appended.
#     TODO make less ugly
#     """
#     taxlevels = list(tax.columns[0::])
#     root_level = tax.columns[0]

#     # root level: should be Kingdom
#     if "kingdom" not in root_level.lower():
#         print(
#             "the highest taxonomic level found is %s, not kingdom as expected. beware"
#             % root_level
#         )

#     # Fills in unknown Kingdom
#     tax[root_level] = tax[root_level].fillna("unknown_%s" % root_level)

#     for i, level in enumerate(taxlevels[1::]):
#         print(level)
#         # indexes of rows where level is missing
#         _missing_l = tax[level].isna()
#         # fill with next higher level
#         tax.loc[_missing_l, level] = "unknown_%s" % level

#     tax_mask = tax.applymap(lambda v: "unknown" in v)
#     tax_fill = tax.copy()
#     tax_fill[tax_mask] = np.nan

#     # lookup table for shifted tax level, e.g. "Class" -> "Phylum")
#     taxlevelshifted = pd.Series(taxlevels[:-1], index=taxlevels[1::])
#     taxlevelshifted.loc["Kingdom"] = "Kingdom"

#     # series with the higher level per ASV/OTU found
#     _highest_level = (tax_fill.isna()).idxmax(axis=1)
#     _highest_level = _highest_level.apply(lambda v: taxlevelshifted.loc[v])

#     # convert taxfill into empty string data frame except where tax is missing a level
#     for cn, c in tax_fill.iteritems():
#         tax_fill[cn] = tax_fill[cn].fillna(_highest_level)
#     tax_fill[~tax_mask] = ""
#     for ix, r in tax_fill.iterrows():
#         whatsit = r.apply(lambda v: "" if v == "" else "_" + tax.loc[ix, v])
#         tax_fill.loc[ix] += whatsit
#     tax_fill += "_of_"

#     # empty strings where tax is not missing values
#     tax_fill[~tax_mask] = ""

#     # pre pend the missing taxlevel to the tax table where tax table is missing
#     tax = tax_fill + tax

#     # if Kingdom was missing:
#     tax.loc[tax.Kingdom.str.contains("unknown"), "Kingdom"] = "unknown_Kingdom"

#     # append the unique sequence id from the index (e.g. ASV_X) to filled values so as to avoid aggregating on filled values.
#     for ix, r in tax.iterrows():
#         for c, v in r.items():
#             if "unknown" in v:
#                 tax.loc[ix, c] = v + "____" + str(ix)

#     return tax
