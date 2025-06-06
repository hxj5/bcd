# xlog.py - logging


import logging
import sys


class XFormatter(logging.Formatter):
    """Logger enables customized format of record."""

    def __init__(self, datefmt = None):
        """
        Parameters
        ----------
        datefmt : str or None, default None
            The format of date.
        """
        super(XFormatter, self).__init__(fmt = None, datefmt = datefmt)
    
    def format(self, record):
        record.message = record.getMessage()
        s = ""
        s += "["

        if record.levelno == logging.DEBUG:
            s += "D"
        elif record.levelno == logging.INFO:
            s += "I"
        elif record.levelno == logging.WARNING:
            s += "W"
        elif record.levelno == logging.ERROR:
            s += "E"
        elif record.levelno == logging.CRITICAL:
            s += "C"
        else:         # unknown
            s += "U"

        if record.module:
            s += "::%s" % record.module
        if record.funcName:
            s += "::%s" % record.funcName
        s += "]"

        if self.datefmt:
            s += "[%s]" % self.formatTime(record, self.datefmt)
        s += " "

        if record.message:
            s += str(record.message)
       
        if record.exc_info:
            # Cache the traceback text to avoid converting it multiple times
            # (it's constant anyway)
            if not record.exc_text:
                record.exc_text = self.formatException(record.exc_info)
        if record.exc_text:
            if s[-1:] != "\n":
                s = s + "\n"
            s = s + record.exc_text
        if record.stack_info:
            if s[-1:] != "\n":
                s = s + "\n"
            s = s + self.formatStack(record.stack_info)
        return s
 


def init_logging(log_file = None, stream = None, 
    fh_level = logging.DEBUG, fh_datefmt = "%Y-%m-%d %H:%M:%S",
    ch_level = logging.INFO,  ch_datefmt = "%Y-%m-%d %H:%M:%S"
):
    """Initialize logging configuration.

    Internally, it calls :func:`~logging.basicConfig` to initialize the 
    logging configuration.

    Parameters
    ----------
    log_file : str or None, default None
        Path to the output logging file.
        If None, do not use it for logging.
    stream
        Steam (file object) used for logging.
        If None, do not use it for logging.
    fh_level : int, default logging.DEBUG
        Logging levels for `log_file`.
        One of {NOTSET, DEBUG, INFO, WARNING, ERROR, CRITICAL} in the
        `logging` module.
    fh_datefmt : str or None, default "%Y-%m-%d %H:%M:%S"
        The date format for `log_file`.
    ch_level : int, default logging.INFO
        Logging levels for `stream`.
        One of {NOTSET, DEBUG, INFO, WARNING, ERROR, CRITICAL} in the
        `logging` module.
    ch_datefmt : str or None, default "%Y-%m-%d %H:%M:%S"
        The date format for `stream`.

    Returns
    -------
    Void.

    Notes
    -----
    This function should be called only once, at the very beginning of the
    running program.
    """
    func = "init_logging"
    if log_file is None and stream is None:
        raise ValueError(
            "[E::%s] at least one of 'log_file' and 'stream' should not be None." % \
            func)
    
    fh = ch = None
    handlers = []
    if log_file:
        fh = logging.FileHandler(log_file, mode = "w")
        fh.setLevel(fh_level)
        fh.setFormatter(XFormatter(datefmt = fh_datefmt))
        handlers.append(fh)

    if stream:
        ch = logging.StreamHandler(stream = stream)
        ch.setLevel(ch_level)
        ch.setFormatter(XFormatter(datefmt = ch_datefmt))
        handlers.append(ch)

    if int(sys.version_info[0]) >= 3 and int(sys.version_info[1]) >= 3:
        logging.basicConfig(
            level = logging.DEBUG,
            handlers = handlers
        )
    else:
        logging.basicConfig(
            level = logging.DEBUG
        )
        if log_file:
            logging.getLogger('').addHandler(fh)
        if stream:
            logging.getLogger('').addHandler(ch)
