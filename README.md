# SWIPE Modules

This repository contains routines that adapt [litebird_sim](https://pypi.org/project/litebird-sim/) to LSPE-SWIPE.

It sould contain:
* A module which takes care of the scanning strategy for a Balloon
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
from swipe_modules import scanning_strategy
import litebird_sim as lbs
import healpy as hp

# Initialize the simulation
sim = lbs.Simulation(
    start_time=0,
    duration_s=36*3600, # 36 hrs of mission
    description="SWIPE simulation",
)

# Scanning strategy 
# delta_time_s is important, it sets the resolution of the computed pointing
# (then pointing gets interpolated at the sampling frequency)
# litebird_sim's default is 60sec, for SWIPE a more frequent sampling, 1sec is reasonably good
sim.generate_spin2ecl_quaternions(scanning_strategy.SwipeScanningStrategy(),delta_time_s=1)

# Initialize the instrument, this might be done using the IMo
# Here spin_boresight_angle_rad is practically ``pi/2 - elevation``
instr = lbs.InstrumentInfo(
    name="swipe",
    spin_boresight_angle_rad=np.deg2rad(40),
)
det = lbs.DetectorInfo(name="foo", sampling_rate_hz=10)

# Observation and pointing
obs, = sim.create_observations(detectors=[det])
pointings = lbs.get_pointings(
    obs,
    sim.spin2ecliptic_quats,
    detector_quats=[det.quat],
    bore2spin_quat=instr.bore2spin_quat,
)

# Plot the scan
nside = 256
npix = hp.nside2npix(nside)
h = np.zeros(npix)
for idet in [0]:
    pixidx = hp.ang2pix(nside, pointings[idet, :, 0], pointings[idet, :, 1])
    pixel_occurrences = np.bincount(pixidx)
    h[0:len(pixel_occurrences)] += pixel_occurrences
hp.mollview(h)
```
<img src="https://user-images.githubusercontent.com/5398538/146160724-a04d6117-e39b-4690-a247-9a7bfdeb6ba5.png" width="500">


## Parameters

The relevant class for the scanning strategy is ``SwipeScanningStrategy``.
It takes the following parameters:

    - ``site_latitude_deg``: Latitude of the launching site

    - ``site_longitude_deg``: Longitude of the launching site

    - ``longitude_speed_deg_per_sec``: Longitude speed of balloon 

    - ``spin_rate_rmp``: Number of rotations around the spin axis per minute

    - ``start_time``: Launching time 
