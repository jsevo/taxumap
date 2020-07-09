import logging
import os
from pathlib import Path

# TODO: The pwd default thing isnt working, probably os.path.join() doesnt work with log_dir = "./".
# Use getcwd next time


def setup_logger(name, verbose=False, debug=False, base_dir=None):

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    general_format = logging.Formatter(
        "%(asctime)s:%(funcName)s:%(levelname)s - %(message)s\n"
    )

    stream_format = logging.Formatter("%(funcName)s:%(levelname)s\n%(message)s\n")

    # save to log folder unless it doesn't exist
    try:
        log_dir = Path(os.path.join(base_dir, "logs")).resolve(strict=True)
    except (TypeError, FileNotFoundError):
        logger.info("No directory log found in pwd. Saving to pwd.")
        log_dir = "./"

    # file handler for warning and above
    fh_warning = logging.FileHandler(os.path.join(log_dir, name + "_warnings.log"))
    fh_warning.setLevel(logging.WARNING)
    fh_warning.setFormatter(general_format)

    # anything warning and above will be printed to console
    sh = logging.StreamHandler()
    sh.setFormatter(stream_format)
    sh.setLevel(logging.WARNING)  # this level and ABOVE

    logger.addHandler(fh_warning)
    logger.addHandler(sh)

    if verbose:
        # set file handler for info and above
        fh_info = logging.FileHandler(os.path.join(log_dir, name + "_info.log"))
        fh_info.setLevel(logging.INFO)
        fh_info.setFormatter(general_format)
        logger.addHandler(fh_info)

    if debug:
        # set file handler for debug and above
        fh_debug = logging.FileHandler(os.path.join(log_dir, name + "_debug.log"))
        fh_debug.setLevel(logging.DEBUG)
        fh_debug.setFormatter(general_format)
        logger.addHandler(fh_debug)

    return logger
