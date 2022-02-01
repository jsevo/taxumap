# Authors: Grant Hussey <grant.hussey@nyulangone.org>
# License: MIT

import pandas as pd
import taxumap.dataloading as dataloading


def validate_weights(weights, agg_levels):
    if weights is None:
        weights = [1] * len(agg_levels)
    else:
        if len(weights) != len(agg_levels):
            raise ValueError("The length of the weights must match that of agg_levels")
    return weights


def validate_microbiome_data_frame(rel_abundances, logger):
    # Validate the rel_abundances dataset
    try:
        rel_abundances.set_index("index_column", inplace=True)
    except KeyError:
        # if it can't set the index to 'index_column', maybe it's already set
        # let's check to see if that is not the case
        if rel_abundances.index.name != "index_column":
            raise ValueError(
                "Your rel_abundances df needs to contain 'index_column' containing sample identifiers of each row containing relative OTU/ASV abundances"
            )

    # lastly, quickly check if the data is compositional.
    dataloading.check_if_compositional(
        rel_abundances, name="locally-supplied microbiota (rel_abundances)",
    )

    return rel_abundances


def validate_microbiome_data(rel_abundances, logger):
    if isinstance(rel_abundances, pd.DataFrame):
        logger.info("Recognized `rel_abundances` parameter as Pandas DataFrame")
    elif isinstance(rel_abundances, str):
        rel_abundances = dataloading.parse_microbiome_data(rel_abundances)
    else:
        raise ValueError(
            "The microbiome data should be provided as either a pd.DataFrame or path to csv file"
        )
    rel_abundances = validate_microbiome_data_frame(rel_abundances, logger)

    return rel_abundances


def fill_taxonomy_table(tax, logger):
    """
    Helper function that fills nan's in a taxonomy table. Such gaps are filled 'from the left' with the next higher non-nan taxonomy level and the lowest level (e.g. OTU# or ASV#) appended.
    TODO make less ugly
    """
    if "kingdom" not in root_level.lower():
        logger.info(
            f"the highest taxonomic level found is {root_level}, not kingdom as expected. not attempting to fill gaps in taxonomy"
        )
        return tax
    taxlevels = list(tax.columns[0::])
    root_level = tax.columns[0]
    # root level: should be Kingdom
    if "kingdom" not in root_level.lower():
        logger.info(
            f"the highest taxonomic level found is {root_level}, not kingdom as expected. beware"
        )
    tax[root_level] = tax[root_level].fillna(f"unknown_{root_level}")
    for _, level in enumerate(taxlevels[1::]):
        # indexes of rows where level is missing
        _missing_l = tax[level].isna()
        # fill with next higher level
        tax.loc[_missing_l, level] = f"unknown_{level}"
    tax_mask = tax.applymap(lambda v: "unknown" in v)
    tax_fill = tax.copy()
    tax_fill[tax_mask] = np.nan
    # lookup table for shifted tax level, e.g. "Class" -> "Phylum")
    taxlevelshifted = pd.Series(taxlevels[:-1], index=taxlevels[1::])
    taxlevelshifted.loc["Kingdom"] = "Kingdom"
    # series with the higher level per ASV/OTU found
    _highest_level = (tax_fill.isna()).idxmax(axis=1)
    _highest_level = _highest_level.apply(lambda v: taxlevelshifted.loc[v])
    # convert taxfill into empty string data frame except where tax is missing a level
    for cn, c in tax_fill.iteritems():
        tax_fill[cn] = tax_fill[cn].fillna(_highest_level)
    tax_fill[~tax_mask] = ""
    for ix, r in tax_fill.iterrows():
        whatsit = r.apply(lambda v: "" if v == "" else "_" + tax.loc[ix, v])
        tax_fill.loc[ix] += whatsit
    tax_fill += "_of_"
    # empty strings where tax is not missing values
    tax_fill[~tax_mask] = ""
    # pre pend the missing taxlevel to the tax table where tax table is missing
    tax = tax_fill + tax

    # if Kingdom was missing:
    tax.loc[tax.Kingdom.str.contains("unknown"), "Kingdom"] = "unknown_Kingdom"

    # append the unique sequence id from the index (e.g. ASV_X) to filled values so as to avoid aggregating on filled values.
    for ix, r in tax.iterrows():
        for c, v in r.items():
            if "unknown" in v:
                tax.loc[ix, c] = v + "____" + str(ix)
    return tax


def ensure_monophyletic_for_hct_dataset(taxonomy):
    taxonomy = taxonomy.loc[
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
    return taxonomy


def normalize_taxonomy(taxonomy, logger):
    try:
        taxonomy = fill_taxonomy_table(taxonomy)
    except:
        logger.info("taxonomy table gap filling failed")
    try:
        taxonomy = ensure_monophyletic_for_hct_dataset(taxonomy)
    except:
        logger.info(
            "ensure_monophyletic_for_hct_dataset failed. likely not hct data set."
        )
    return taxonomy


def validate_taxonomy(taxonomy, logger):
    # Set taxonomy df
    try:
        if taxonomy is None:
            raise ValueError(
                "taxonomy is None. Provide a pd.DataFrame taxonomy table or path to a .csv"
            )
        elif isinstance(taxonomy, pd.DataFrame):
            taxonomy.columns = map(str.capitalize, taxonomy.columns)
            dataloading.check_tax_is_consistent(taxonomy)
        elif isinstance(taxonomy, str):
            taxonomy = dataloading.parse_taxonomy_data(taxonomy)
    except:
        logger.critical("Taxonomy validation failed")
    try:
        normalize_taxonomy(taxonomy, logger)
    except:
        logger.info(
            "An exception occurred when harmonizing taxonomy date (filling gaps, testing if all monophyletic)"
        )
    return taxonomy


def validate_inputs(weights, rel_abundances, taxonomy, agg_levels, logger):

    weights = validate_weights(weights, agg_levels)
    rel_abundances = validate_microbiome_data(rel_abundances, logger)
    taxonomy = validate_taxonomy(taxonomy, logger)

    return (weights, rel_abundances, taxonomy)
