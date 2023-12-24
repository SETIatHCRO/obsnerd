Initial observing recipe

When you are ready to go for the session claim the antennas:
./obsnerd.py start <INITIALS HERE>

Use S.O.P.P. to pick a satellite and a frequency:
 - for initial testing, GOES-16 as a geostationary at around 1680 MHz at az,el = 121.958,23.603 is a good start (these are the defaults)
 - ideally however we want to point in front of an LEO and measure it passing through

Note, there are two gnuradio companion scripts:  nrdz_use.py and nrdz_nofiles.py
If you don't want to record the data, use the nrdz_nofiles.py option and skip the file conversion/viewing steps below. 

For the 'freq' and 'move' commands, if you leave off the argument it will use the defaults.  The examples below show explicit arguments for them.

After you've chosen the target satellite and frequency set the frequency and move the antenna (and add a note):
./obsnerd.py note 'GOES 16'
./obsnerd.py freq 1680
./obsnerd.py move 121.958,23.603

When the antenna has gotten to target, in another window start the USRP:
./nrdz_use.py

When you are done with the observation, use the 'X' in the corner to stop the USRP.  This writes a file called nrdz

If you have recorded data, you should convert the raw nrdz datafile to HDF5 setting a useful output filename, e.g. <satname>_YYMMDD_HHMM.h5
./f2h5.py goes16_231223_1345.h5

To look at the data, you can use:
./view_file.py goes16_231223_1345.h5

The goal is to find the satellite frequency and crossing time, so you can zoom around the datafile to find the time and frequency.  view_file.py is pretty rudimentary, so if you zoom in you'll likely need to add y and/or x ticks
./view_file.py goes16_231223_1345.h5 -y 20

We'll want to capture what you've learned (TBD), but for now add to the obslog:
./obsnerd.py note 'Found the millionth starlink at 12:35:23 instead of the expected 12:35:50 yada yada'

If you want to keep it, move the hdf5 file somewhere (TBD).

When you are done, release the antennas:
./obsnerd.py end
