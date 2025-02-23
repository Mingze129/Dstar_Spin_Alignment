'''
Author: Mingze Li
Email: mingze.li@cern.ch
Date: 2025-02-20
'''
__version__ = "2.0.0"
__author__ = "Mingze Li"
__email__ = "mingze.li@cern.ch"

from .Formula import *
from .DataOps import DataOps
from .FitOps import FitOps
from .SpinOps import SpinOps
from .FracOps import FracOps
from .Logger import logger_config

__all__ = [
    "BetheBlochAlephNP",
    "gausscale",
    "linear",
    "gausgauslin",
    "gauslin",
    "gauss",
    "constant_func",
    "DataOps",
    "FitOps",
    "SpinOps",
    "FracOps",
    "logger_config"
]