import logging


def setup_logger(name, verbose=False, debug=False):

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    general_format = logging.Formatter(
        "%(asctime)s:%(funcName)s:%(levelname)s - %(message)s\n"
    )

    stream_format = logging.Formatter("%(funcName)s:%(levelname)s\n%(message)s\n")

    # file handler for warning and above
    fh_warning = logging.FileHandler(name + "_warnings.log")
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
        fh_info = logging.FileHandler(name + "_info.log")
        fh_info.setLevel(logging.INFO)
        fh_info.setFormatter(general_format)
        logger.addHandler(fh_info)

    if debug:
        # set file handler for debug and above
        fh_debug = logging.FileHandler(name + "_debug.log")
        fh_debug.setLevel(logging.DEBUG)
        fh_debug.setFormatter(general_format)
        logger.addHandler(fh_debug)

    return logger
