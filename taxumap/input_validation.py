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
    return taxonomy


def validate_inputs(weights, rel_abundances, taxonomy, agg_levels, logger):

    weights = validate_weights(weights, agg_levels)
    rel_abundances = validate_microbiome_data(rel_abundances, logger)
    taxonomy = validate_taxonomy(taxonomy, logger)

    return (weights, rel_abundances, taxonomy)
