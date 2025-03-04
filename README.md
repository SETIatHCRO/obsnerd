#Terminology
- source - a unique source name.
- obsid - a unique observation identifier:  `source_MJD{.5f}`
- obsrec - a unique observation/hardware identifier:  `obsid_LO_CNODE`
- obsrec file - the filename holding the obsrec information (generally `obsrec.npz` ***note***)
- experiment - a session looking at sources (typically within an MJD day or two)
- obsinfo - a file containing information on obsids for a given experiment:  `obsinfo_MJD.json`

***note:*** uvh5 files are put into directories that map to an obsrec:

`/mnt/primary/ata/projects/pID#/YYYY-MM-DD-HH:MM:SS/uvh5...../LoA.C0352/uvh5......uvh5`

#Observing Recipe
The focus right now is to observe the Starlink Direct-to-Cell (DTC) satellites.

##Find satellites##
Need to find and set up the satellites.  This is currently interactive.  The main problem is that THERE ARE SO MANY!

Set up the planning tool from with ipython and get 'tracks'

    from obsnerd import onp_plan
    plan = onp_plan.Plan()
    plan.get_tracks(satname='some-search-string', start='now', duration=number_of_minutes)

For example `plan.get_tracks(satname='1112', start='now', duration=8*60)`

You can get more tracks using `get_tracks` again

If this looks like a good set, choose the ones you want (choose track number and +/- if you want an active ODS):

    plan.choose_tracks()

If you want to continue,

    plan.proc_tracks()

This will write two files **ods.json** and **obsinfo_MJD.json**

Check that they both look reasonable (primarily the ods file)

##Observe##
1. put that ods file into '/opt/mnt/share/ods_rados/ods_rados.json'
2. on the VNC, open a terminal and type `aoctkuser.py --enable-rados`
3. enable the google calendar and make a calendar entry (could have done this earlier)
4. hit the *Observe* button and if you are confident select **yes**
5. sit back and watch the action

##Process the data for analysis##
The data now sits in that deeply nested directory structure in the ***note*** above.  You now want to dump the autocorrelations for the antennas you want to much smalller 'npz' files, which means you need to find them, generate a bash script, and run the bash script.

Do this will logged into **obs-node1** and in the **~/rfsoc\_obs\_scripts/p054** directory

1. Find them by typing `on_gen_dump.py ?` at the terminal
2. Generate the scripts by typing `on_gen_dump.py <date-from-above>` (use `on_gen_dump.py -h` for options).  This generates three scripts:  *copy\_files.sh*, *dump\_autos.sh*, *download\_files.sh*.
3. Run the script `bash copy_files.sh` to copy the files over to the local drive (currently seems faster to copy, but will still take a very long time)
4. Run the generated script `bash dump_autos.sh`
5. Move the *download\_files.sh* to your local machine and run to download the files.
