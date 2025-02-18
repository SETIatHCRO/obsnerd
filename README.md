Terminology:
  source - a unique source name.  For starlink that is S#####
  obsid - a unique observation identifier:  source_MJD{.4f}
  obsrec - a unique observation/hardware identifier:  <obsid>_LO_CNODE
  Obsrec file - the filename holding the obsrec information (generally <obsrec>.npz)
  experiment - a session looking at sources (typically within an MJD day or two)
  obsinfo - a file containing information on obsids for a given experiment:  obsinfo_MJD.json
  obsmeta - a json file containing index to the obsinfo files: obsmeta.json


Observing Recipe (most ingredients missing at the moment...)

The focus right now is to observe the Starlink Direct-to-Cell (DTC) satellites.

- 1 Need to find and set up the satellites.  This is currently interactive.  The main problem is that there ARE SO MANY!

Set up the planning tool from with ipython
> from obsnerd import onp_plan
> plan = onp_plan.Plan()
The line below will look for satellite names containing 'satname' and find all tracks starting now and lasting 8 hours.
Try satname='1112'
> plan.get_tracks(satname='some-search-string', start='now', duration=8 * 60)
> plan.choose_tracks()
Is you want to continue,
> plan.proc_tracks()

Now put that ods file into '/opt/mnt/share/ods_rados/ods_rados.json'

- 2 Run the observation: 
    aoctkuser.py --enable-rados
    then hit the Observe button

- 3 Copy the data into local npz files along with the obsinfo file
    $ on_gen_dump.py ?
    $ on_gen_dump.py <date-from-above-look>