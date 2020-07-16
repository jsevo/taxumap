import sys

from taxumap.custom_logging import setup_logger

logger_errors = setup_logger("errors", verbose=False, debug=False)


def throw_unknown_save_error(e):
    logger_errors.exception(
        "\nUnknown error has occured. Cannot save embedding as instructed."
    )
    sys.exit(2)


def _name_value_error(e, fp_param, df_param):
    logger_errors.exception(
        "Please provide the constructor with one of the following: a filepath to your {} file via parameter `{}`, or with an initialized variable for your {} dataframe via the parameter `{}`".format(
            df_param, fp_param, df_param, df_param
        )
    )
