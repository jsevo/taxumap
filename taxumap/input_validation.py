# Authors: Grant Hussey <grant.hussey@nyulangone.org>
# License: MIT

import pandas as pd
import taxumap.dataloading as dataloading
import numpy as np
import warnings


def validate_weights(weights, agg_levels):
    if weights is None:
        weights = [1] * len(agg_levels)
    else:
        if len(weights) != len(agg_levels):
            raise ValueError("The length of the weights must match that of agg_levels")
    weights = weights.copy()
    return weights


def validate_microbiome_data_frame(microbiota_data):
    # Validate the microbiota_data dataset
    try:
        microbiota_data.set_index("index_column", inplace=True)
    except KeyError:
        # if it can't set the index to 'index_column', maybe it's already set
        # let's check to see if that is not the case
        if microbiota_data.index.name != "index_column":
            raise ValueError(
                "Your microbiota_data df needs to contain 'index_column' containing sample identifiers of each row containing relative OTU/ASV abundances"
            )

    # lastly, quickly check if the data is compositional.
    dataloading.check_if_compositional(
        microbiota_data,
        name="locally-supplied microbiota (microbiota_data)",
    )
    microbiota_data = microbiota_data.copy()

    return microbiota_data


def validate_microbiome_data(microbiota_data):
    if isinstance(microbiota_data, pd.DataFrame):
        warnings.warn("Recognized `microbiota_data` parameter as Pandas DataFrame")
    elif isinstance(microbiota_data, str):
        microbiota_data = dataloading.parse_microbiome_data(microbiota_data)
    else:
        raise ValueError(
            "The microbiome data should be provided as either a pd.DataFrame or path to csv file"
        )
    microbiota_data = validate_microbiome_data_frame(microbiota_data)

    return microbiota_data


# def fill_taxonomy_table(tax):
#     """
#     Helper function that fills nan's in a taxonomy table. Such gaps are filled 'from the left' with the next higher non-nan taxonomy level and the lowest level (e.g. OTU# or ASV#) appended.
#     TODO make less ugly
#     """
#     if "kingdom" not in root_level.lower():
#         warnings.warn(
#             f"the highest taxonomic level found is {root_level}, not kingdom as expected. not attempting to fill gaps in taxonomy"
#         )
#         return tax
#     taxlevels = list(tax.columns[0::])
#     root_level = tax.columns[0]
#     # root level: should be Kingdom
#     if "kingdom" not in root_level.lower():
#         warnings.warn(
#             f"the highest taxonomic level found is {root_level}, not kingdom as expected. beware"
#         )
#     tax[root_level] = tax[root_level].fillna(f"unknown_{root_level}")
#     for _, level in enumerate(taxlevels[1::]):
#         # indexes of rows where level is missing
#         _missing_l = tax[level].isna()
#         # fill with next higher level
#         tax.loc[_missing_l, level] = f"unknown_{level}"
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


def fill_tax_table(tax):
    """Fills missing values in the taxonomy table. Will recognize only 'np.nan' data types as empty values.

    Args:
        tax (pd.DataFrame): Dataframe with index of ASV/OTU and columns of left -> right increasing specificity in taxonomy (e.g., Kingdom -> Species)

    Output:
        new_tax (pd.DataFrame): Properly-filled taxonomy dataframe
    """

    if len(tax.index) != len(tax.index.unique()):
        print(
            "Repeated OTUs/ASVs in the taxonomy index. Check to make sure there is only _one_ entry per OTU in taxonomy table."
        )

    table_name = tax.index.name  # Important - don't remove this and its corresponding stpe below.

    container = []
    for otu in tax.index:

        series = tax.loc[otu]

        # If there are no NaNs in the OTU, don't do anything.
        if (~series.isna()).all():
            container.append(series)

        # However, if NaNs do exist, fill the taxonomy "from-the-left"
        else:

            # First, initialize the first "last_not_nan" in case the highest level is missing.
            last_not_nan_index = series.index[0]
            last_not_nan_value = f"unk_{last_not_nan_index}"

            for i in range(len(series)):

                # If the label is nan, fill it
                if series.isna().iloc[i]:
                    series.iloc[
                        i
                    ] = f"unk_{series.index[i]}_of_{last_not_nan_index}_{last_not_nan_value}__{series.name}"

                # Else, keep that non-nan value in case the next label is nan
                else:
                    last_not_nan_index = series.index[i]
                    last_not_nan_value = series.iloc[i]

            container.append(series)

    new_tax = pd.concat(container, axis=1).T
    new_tax.index.name = table_name

    return new_tax


def _fill_tax_table_old(tax):
    """Fills missing values in the taxonomy table. Will recognize only 'np.nan' data types as empty values.

    Args:
        tax (pd.DataFrame): Dataframe with index of ASV/OTU and columns of left -> right increasing specificity in taxonomy (e.g., Kingdom -> Species)

    Output:
        new_tax (pd.DataFrame): Properly-filled taxonomy dataframe
    """
    if len(tax.index) != len(tax.index.unique()):
        print(
            "Repeated OTUs/ASVs in the taxonomy index. Check to make sure there is only _one_ entry per OTU in taxonomy table."
        )

    # MUST be in increasing specificity order (Kingdom -> Species)
    # OTU/ASV must be the INDEX.
    tax_labels = tax.columns
    table_name = tax.index.name  # Important - don't remove this and its corresponding stpe below.

    # Gather all OTUs to iterate over
    otus = tax.index.unique()

    new_tax = []  # Collector for new taxonomy pd.Series
    for otu in otus:

        series = tax.loc[otu]

        # If there are no NaNs in the OTU, don't do anything.
        if (~series.isna()).all():
            new_tax.append(series)

        # However, if NaNs do exist, fill the taxonomy "from-the-left"
        else:
            first_nan = np.argwhere(series.isna().values == True)[0][0]

            # In case "Kingdom" is NaN (or other highest level taxa)
            if first_nan == 0:
                last_not_nan = first_nan
            else:
                last_not_nan = first_nan - 1

            ##### Below commented-out code I'm saving here, ignore #####
            # for i in range(first_nan, len(series)):
            #     series.iloc[i] = f'unk_{series.index[i]}_of_{series.index[i-1]}_{series.iloc[i-1]}'
            #####                                                  #####

            # Perform "fill-from-the-left"
            # For each and every NaN, fill it with the last non-NaN taxonomy, and append the ASV/OTU name at the end as well.
            for i in range(first_nan, len(series)):

                # In case "Kingdom" is NaN (or other highest level taxa)
                if i == 0:
                    series.iloc[i] = f"unk_{series.index[i]}"
                else:
                    series.iloc[
                        i
                    ] = f"unk_{series.index[i]}_of_{series.index[last_not_nan]}_{series.iloc[last_not_nan]}"

            # Add in the ASV/OTU name to the end of every unknown
            for i in range(first_nan, len(series)):
                series.iloc[i] = f"{series.iloc[i]}__{otu}"

            new_tax.append(series)

    new_tax = pd.concat(new_tax, axis=1).T

    # This name gets erased in the above transformation, so return it.
    new_tax.index.name = table_name

    return new_tax


def ensure_monophyletic_for_hct_dataset(taxonomy):
    taxonomy.loc[
        taxonomy.Genus.str.contains("metagenome") & ~taxonomy.Genus.str.contains("ASV"),
        "Genus",
    ] = (
        taxonomy.loc[
            taxonomy.Genus.str.contains("metagenome")
            & ~taxonomy.Genus.str.contains("ASV"),
            "Genus",
        ]
        + "___"
        + taxonomy.loc[
            taxonomy.Genus.str.contains("metagenome")
            & ~taxonomy.Genus.str.contains("ASV"),
            "Genus",
        ].index
    )

    taxonomy.loc[taxonomy.Family == "Family XI", "Family"] = (
        taxonomy.loc[taxonomy.Family == "Family XI", "Family"]
        + "_order_"
        + taxonomy.loc[taxonomy.Family == "Family XI", "Order"]
    )
    taxonomy.loc[taxonomy.Family == "Family gut metagenome", "Family"] = (
        taxonomy.loc[taxonomy.Family == "Family gut metagenome", "Family"]
        + taxonomy.loc[taxonomy.Family == "Family gut metagenome"].index
    )

    return taxonomy

def test_taxonomy_is_monophyletic(taxonomy):
    try:
      for (p,c) in zip(taxonomy.columns[0:-1],taxonomy.columns[1::]):
          gc = taxonomy[[p,c]].groupby(c).agg(lambda x: len(np.unique(x))<=1)
          if np.all(gc.values):
              pass
          else:
              print(p,c)
              print('not monophyletic')
    except:
        print("Could not ensure monophyly")

    


def normalize_taxonomy(taxonomy):
    try:
        taxonomy = fill_tax_table(taxonomy)
    except:
        warnings.warn("taxonomy table gap filling failed")
    try:
        ensure_monophyletic_for_hct_dataset(taxonomy)
    except:
        warnings.warn("ensure_monophyletic_for_hct_dataset failed. likely not hct data set.")
    return taxonomy


def validate_taxonomy(taxonomy):
    # Set taxonomy df
    try:
        if taxonomy is None:
            raise ValueError(
                "taxonomy is None. Provide a pd.DataFrame taxonomy table or path to a .csv"
            )
        elif isinstance(taxonomy, pd.DataFrame):
            taxonomy = taxonomy.copy()
            taxonomy.columns = map(str.capitalize, taxonomy.columns)
            dataloading.check_tax_is_consistent(taxonomy)
        elif isinstance(taxonomy, str):
            taxonomy = dataloading.parse_taxonomy_data(taxonomy)
    except:
        warnings.warn("Taxonomy validation failed")
    try:
        taxonomy = normalize_taxonomy(taxonomy)

    except:
        warnings.warn(
            "An exception occurred when harmonizing taxonomy date (filling gaps, testing if all monophyletic)"
        )
    try:
        test_taxonomy_is_monophyletic(taxonomy)
    except:
        warnings.warn("Review taxonomy table, it appears not to be monophyletic.")  
    return taxonomy


def validate_inputs(weights, microbiota_data, taxonomy, agg_levels):
    weights = validate_weights(weights, agg_levels)
    microbiota_data = validate_microbiome_data(microbiota_data)
    taxonomy = validate_taxonomy(taxonomy)
    return (weights, microbiota_data, taxonomy)
