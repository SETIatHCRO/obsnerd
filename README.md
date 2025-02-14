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