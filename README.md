Initial observing recipe

Log into gnuradio and cd into the obsnerd directory

Note that for (hopefully) lack of confusion, all times should be in UTC (which is what the gnuradio server is set to).  I make this goof all of the time.

When you are ready to go for the session, claim the antennas (if this fails, run 'obsnerd.py end' first):
obsnerd.py start <INITIALS HERE>

Use S.O.P.P. to pick a satellite and a frequency:
 - I've adapted the S.O.P.P. example.py as find_sats.py.  To get GPS satellites, can type
find_sats.py -s NAVS
 - ideally however we want to point in front of an LEO and measure it passing through

Note, there are two gnuradio companion scripts:  nrdz_use.py and nrdz_nofiles.py
If you don't want to record the data, use the nrdz_nofiles.py option and skip the file conversion/viewing steps below. 

After you've chosen the target satellite and frequency set the frequency, move the antenna and add the source information
obsnerd.py source nav72  2023-12-24T18:01:36
obsnerd.py freq 1575
obsnerd.py move 238.275,51.614

When the antenna has gotten to target, in another window start the USRP:
nrdz_use.py

When you are done with the observation, use the 'X' in the corner to stop the USRP.  This writes a file called nrdz

If you have recorded data, you should convert the raw nrdz datafile to HDF5 setting a useful output filename: <satname>_YYMMDD_HHMMSS.h5
This is done automatically (make sure you run the obsnerd.py source command)
f2h5.py

To look at the data, you can use:
onview.py nav72_231224_180136.h5 -p [wf/series/spectra]
the -p options are
    wf - waterfall
    series - all time series, which also plots a line at the time specified in the filename
    spectra - all spectra
if you want log use -l

The goal is to find the satellite frequency and crossing time, for which the series option is probably best.

To look around the whole datafile, the waterfall is good where you can zoom around the datafile to find the time and frequency.  
onview.py is pretty rudimentary, so if you zoom in you'll likely need to add y and/or x ticks
onview.py nav72_231224_180136.h5 -y 20 -d

Add the results to
https://docs.google.com/spreadsheets/d/1z4b35RkIdoUa2pU_zBnD5AoqiXg8NwgIvBci6oa-Tmo/edit#gid=0

If you want to keep it, move the hdf5 file somewhere (TBD).

When you are done, release the antennas:
obsnerd.py end

Offline, make a .txt file for each .5h file and add the information from the onlog.log and the googlesheet, and then describe the RFI in the field.

================COMMANDS ONLY==================
- updatetle.py
- obsnerd.py start 'initials'
- ... loop for multiple satellites
- obsnerd.py source nav72 2023-12-24T18:01:36
- obsnerd.py freq 1575
- obsnerd.py move 238.275,51.614
- nrdz_use.py
- <END OBS X>
- f2h5.py
- ... observe more with above loop
- obsnerd.py end

- ...analyze
- onview nav72_231224_180136.h5 -p series -d
...etc
...update googlesheet
...make .txt file (nav72_231224_180136.txt)
