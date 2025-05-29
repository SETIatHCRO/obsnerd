#Terminology
- **source** - a unique source name.
- **obsid** - a unique observation identifier:  `<SOURCE>_<MJD{.5f}>`
- **obsrec** - a unique observation/hardware identifier:  `<OBSID>_<LO>_<CNODE>`
- **obsrec file** - the filename holding the obsrec information (generally `<OBSREC>.npz` ***note***)
- **experiment** - a session looking at sources (typically within an MJD day or two)
- **obsinfo** - a file containing information on obsids for a given experiment:  `obsinfo_<MJD>.json`
- **ods** - operational data sharing

***note:*** uvh5 files are put into directories that map to an obsrec:<br>
`/mnt/primary/ata/projects/pID#/YYYY-MM-DD-HH:MM:SS/uvh5...../Lo[A/B/C/D].C[####]/uvh5......uvh5`

#Observing Recipe
The focus right now is to observe the Starlink Direct-to-Cell (DTC) satellites.

##Find satellites##
Need to find and set up the satellites.  This is currently clumsily interactive.  The main problem is that THERE ARE SO MANY!

Set up the planning tool from with ipython and get 'tracks'

    from obsnerd import onp_plan
    plan = onp_plan.Plan()
    plan.get_tracks(satname='some-search-string', start='now+30m', duration=number_of_minutes)

For example `plan.get_tracks(satname='*', start='now+30m', duration=2*60)`

You can get more tracks using `get_tracks` again

If this looks like a good set, choose the ones you want (there is an auto-choose mode by default):

    plan.choose_tracks(auto=True)

If you want to continue with that set,

    plan.proc_tracks()

This will write two files **ods_\<MJD>.json** and **obsinfo_\<MJD>.json**

Check that they both look reasonable (primarily the ods file)

##Observe##
1. on the VNC open a terminal
2. put the file `ods_<MJD>.json` onto obs-node1 as `~/rfsoc_obs_scripts/p054/ods_rados.json`
3. type `on_obs_prep.py --add-to-calendar`
4. type `aoctkuser.py --enable-rados`
5. hit the *Observe* button and if you are confident select **yes** twice.
6. sit back and watch the action

##Process the data##
The data now sits in that deeply nested directory structure in the ***note*** above.  You now want to dump the autocorrelations for the antennas you want to much smalller 'npz' files, which means you need to find them, generate a bash script, and run the bash script.

Do this while logged into **obs-node1** and in the **~/rfsoc\_obs\_scripts/p054** directory

1. Find them by typing `on_gen_dump.py ?` at the terminal
2. Generate the scripts by typing `on_gen_dump.py <date-from-above>` (use `on_gen_dump.py -h` for options).  This generates two scripts:  *dump\_autos.sh* and *download\_files.sh*.
3. Run the generated script `bash dump_autos.sh` (this will take a very long time so do it under a screen).  This generates npz files containing the specified autocorrelations per antenna/polarization/cnode.  The name of the file is the **obsrec** defined above
4. Move the *download\_files.sh* to your local machine and run to download the files if you wish to process them locally.

##Generate the dashboard##
In the directory with the obsinfo file and the npz files in the appropriate directory specified therein (if not the default **data** directory), you can view the data and generate the dashboard.

To view a dashboard of a source from the command line: `on_look.py <SOURCE>`

For options, type `on_look.py -h`

E.g. `on_look.py STARLINK11139DTC2_A_C0352 --lo A --dash`

If you want to generate and save all of the data in the obsinfo file, type `on_gen_dash.py <MJD>`, which will write a bash script file (default *dash.sh*).

Running the 'dash.sh' script won't display any data, but rather save them to *png* files.  You can concatenate them and write them to a pdf file by using a image tool like magick:  `magick *.png out.pdf`
