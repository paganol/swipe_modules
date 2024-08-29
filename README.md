# SWIPE Modules

This repository contains routines that adapt [litebird_sim](https://pypi.org/project/litebird-sim/) to LSPE-SWIPE.

It sould contain:
* Modules which takes care of the scanning strategy for a Balloon
* A number of modules that can ingest systematic effects in the timelines 


## Installation

Set up the enviroment (if necessary)
```
# Create a conda environment
conda create -n swipe_env python=3.8

# Activate the environment
conda activate swipe_env
```
Obviously a similar structure works also for ``virtualenv``


Install the package
```
git clone https://github.com/paganol/swipe_modules.git
cd swipe_modules
pip install .
```
It can be installed in editable mode with ``pip install -e .``


## Install IMO

Download the ``IMO`` from the github [repository](https://github.com/paganol/swipe_imo/) and run
```
python -m litebird_sim.install_imo
```
Select option ``3`` and enter the directory ``swipe_imo/IMO`` insert a descriptive name, something like SWIPE_IMO, and save with ``s``

## Usage

```python
from swipe_modules import SwipeSpinScanningStrategy
import litebird_sim as lbs
import healpy as hp
import astropy

# Initialize the simulation
sim = lbs.Simulation(
    start_time=astropy.time.Time("2026-09-01T16:00"),
    duration_s=36*3600, # 36 hrs of mission
    description="SWIPE simulation",
    random_seed=12345,
)

# Scanning strategy
# delta_time_s is important, it sets the resolution of the computed pointing
# (then pointing gets interpolated at the sampling frequency)
# litebird_sim's default is 60sec, for SWIPE a more frequent sampling, 1sec is reasonably good
scanning = SwipeSpinScanningStrategy()
sim.set_scanning_strategy(scanning,delta_time_s=1)

# Initialize the instrument, this might be done using the IMo
# Here spin_boresight_angle_rad is practically ``pi/2 - elevation``
instr = lbs.InstrumentInfo(
    name="SWIPE",
    spin_boresight_angle_rad=np.deg2rad(90.0-51.5),
)

det = lbs.DetectorInfo(name="foo", sampling_rate_hz=100)

# Observation and pointing
sim.set_instrument(instr)
obs, = sim.create_observations(detectors=[det])
sim.prepare_pointings()
sim.precompute_pointings()

# Plot the scan
nside = 128
npix = hp.nside2npix(nside)
h = np.zeros(npix)
for idet in [0]:
    pixidx = hp.ang2pix(nside, obs.pointing_matrix[idet, :, 0], obs.pointing_matrix[idet, :, 1])
    pixel_occurrences = np.bincount(pixidx)
    h[0:len(pixel_occurrences)] += pixel_occurrences
hp.mollview(h)
```
<img src="https://user-images.githubusercontent.com/5398538/146160724-a04d6117-e39b-4690-a247-9a7bfdeb6ba5.png" width="500">


## Parameters

The relevant class for the spin scanning strategy is ``SwipeSpinScanningStrategy``.
It takes the following parameters:

    - `site_latitude_deg`: latitude of the launching site in deg

    - `site_longitude_deg`: longitude of the launching site in deg

    - `longitude_speed_deg_per_sec`: longitude speed in deg/sec

    - `spin_rate_rmp`: the number of rotations per minute (RPM) around
    the spin axis

    - `start_time`: an ``astropy.time.Time`` object representing the
      start of the observation. It's currently unused, but it is meant
      to represent the time when the rotation starts (i.e., the angle
      ωt is zero).

    - `balloon_latitude_deg`: latitude of a tabulated trajectory

    - `balloon_longitude_deg`: longitude of a tabulated trajectory

    - `balloon_time`: list of astropy.time.Time

For a constant elevation scan in azimuth the class to use is ``SwipeRasterScanningStrategy``.
It takes the following parameters:

    - `site_latitude_deg`: latitude of the launching site in deg

    - `site_longitude_deg`: longitude of the launching site in deg

    - `longitude_speed_deg_per_sec`: longitude speed in deg/sec

    - `azimuth_start_deg`: start azimuth for constant elevation scan  

    - `azimuth_amplitude_deg`: amplitude in azimuth for constant elevation scan

    - `azimuth_scan_speed_deg_per_s`: scan speed in azimuth

    - `start_time`: an ``astropy.time.Time`` object representing the
      start of the observation. It's currently unused, but it is meant
      to represent the time when the rotation starts (i.e., the angle
      ωt is zero).

    - `balloon_latitude_deg`: latitude of a tabulated trajectory

    - `balloon_longitude_deg`: longitude of a tabulated trajectory

    - `balloon_time`: list of astropy.time.Time

In both cases if ``balloon_latitude_deg``, ``balloon_longitude_deg`` and ``balloon_time``
are provided a tabulated trajectory will be used.
