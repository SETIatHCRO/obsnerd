import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
import time


hcro = EarthLocation(lat=40.8178049*u.deg, lon=-121.4695413*u.deg, height=986*u.m)

tz = -8.0*u.hour
obstime = Time('2023-12-31T12:00:00') - tz

use_b = 0.0 * u.deg
gp_l = np.arange(0, 359) * u.deg
gp_b = np.ones(len(gp_l)) * use_b


gal = SkyCoord(frame='galactic', l=gp_l, b=gp_b)
radec = gal.transform_to('icrs')
azel = radec.transform_to(AltAz(obstime=obstime, location=hcro))

el_starting = 30.0
above_horizon = np.where(azel.alt.value > el_starting)
start_horizon = above_horizon[0][0]
starting_l = gal.l.value[start_horizon]

plt.figure('At T=0')
plt.plot(radec.ra.value[above_horizon], radec.dec.value[above_horizon], '.')
plt.plot(azel.az.value[above_horizon], azel.alt.value[above_horizon], '.')

print(f"UTC {obstime}")
print(f"Start at l = {starting_l}, b={gal.b.value[start_horizon]}") 
print(f"         RA={radec.ra.value[start_horizon]}, Dec={radec.dec.value[start_horizon]}")
print(f"         Az={azel.az.value[start_horizon]}, El={azel.alt.value[start_horizon]}")

time_to_track = 20.0  # minutes
track_times = obstime + np.arange(0.0, time_to_track * 60.0, 1.0) * u.second
lstep = 0.1*u.deg  #deg/sec

ts = time.mktime(obstime.datetime.timetuple())
az = []
el = []
with open('track.txt', 'w') as fp:
    for i in range(len(track_times)):
        this_l = starting_l*u.deg  + i*lstep
        if this_l > 360.0*u.deg:
            this_l = this_l - 360.0*u.deg
        gal = SkyCoord(frame='galactic', l=this_l, b=use_b)
        radec = gal.transform_to('icrs')
        azel = radec.transform_to(AltAz(obstime=track_times[i], location=hcro))
        az.append(azel.az.value)
        el.append(azel.alt.value)
        print(f"{int(ts*1E9) + i*1000000000},{azel.az.value},{azel.alt.value}", file=fp)
plt.plot(az, el, lw=4)


