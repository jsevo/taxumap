import logging


def setup_logger(name):

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    general_format = logging.Formatter(
        "%(asctime)s:%(funcName)s:%(levelname)s - %(message)s\n"
    )
    stream_format = logging.Formatter("%(funcName)s:%(levelname)s\n%(message)s\n")

    # set file handler for info and above
    fh_info = logging.FileHandler(name + "_debug.log")
    fh_info.setLevel(logging.INFO)
    fh_info.setFormatter(general_format)

    # file handler for warning and above
    fh_warning = logging.FileHandler(name + "_warnings.log")
    fh_warning.setLevel(logging.WARNING)
    fh_warning.setFormatter(general_format)

    # anything warning and above will be printed to console
    sh = logging.StreamHandler()
    sh.setFormatter(stream_format)
    sh.setLevel(logging.WARNING)  # this level and ABOVE

    logger.addHandler(fh_info)
    logger.addHandler(fh_warning)
    logger.addHandler(sh)

    return logger
