from litebird_sim import ScanningStrategy
from numpy import np

class SwipeScanningStrategy(ScanningStrategy):
    """A class containing the parameters of the sky scanning strategy 
    for SWIPE

    """
    def __init__(
        self,
        spin_sun_angle_rad,
        spin_rate_hz,
        start_time=astropy.time.Time("2023-01-01", scale="tdb"),
    ):


