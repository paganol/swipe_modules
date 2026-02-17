"""SWIPE modules: sky scanning strategies for LSPE-SWIPE."""

from .common import data_directory
from .raster_scanning_strategy import SwipeRasterScanningStrategy
from .spin_scanning_strategy import SwipeSpinScanningStrategy

__version__ = "0.1.0"
__author__ = "Luca Pagano"
__author_email__ = "luca.pagano@unife.it"

__all__ = [
    "SwipeSpinScanningStrategy",
    "SwipeRasterScanningStrategy",
    "data_directory",
]
