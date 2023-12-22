Initial observing recipe

Using S.O.P.P. pick a satellite and a frequency
 - for initial testing, GOES-16 as a geostationary at around 1680 MHz at az/el -> 121.958/23.603 is a good start (these are the defaults)
 - ideally however we want to point in front of an LEO and measure it passing through

Add to the _top_ of the "Observing" file
`./obsnerd.py start` to get the antenna into the desired group
`./obsnerd.py freq -f 1680` to set the frequency
`./nrdz_nofiles.py` to get the spectrum analyzer going without writing files
`./obsnerd.py move --az 121.958 --el 23.603
 
Explore a bit.  When you are ready to take data, stop nrdz_nofiles.py and start
`./nrdz_use.py` to write data to a file called 'nrdz'
this will keep writing data, so the file may get large

WHen you are done, convert the output files to hdf5 using
`./f2h5.py nrdz -o <fn>.h5`

Move the hdf5 file somewhere (TBD) and delete both files to save room.
