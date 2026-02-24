# from .on_observation import Observation
from param_track import Parameters

print("I THINK on_obsinfo.py IS DEPRECATABLE")
class Obsinfo(Parameters):
    def __init__(self, filename=None, **kwargs):
        ptinit = kwargs.pop('ptinit', None)
        super().__init__(ptnote="Obsinfo parameters", ptstrict=True, pterr=False, ptverbose=False, ptinit=ptinit,
                         pttype=False, pttypeerr=False, ptsetunits=False, **kwargs)
        self.ptadd(filename=filename)
        if hasattr(self, 'observations') and isinstance(self.observations, dict):
            self.proc_observations()
        if hasattr(self, 'filters') and isinstance(self.filters, list):
            self.proc_filters()

    def proc_observations(self, sources=None):
        """ Process sources from configuration. """
        if sources is not None:
            self.observations = {key: {} for key in sources}
        observations = {}
        for key, val in self.observations.items():
            #observations[key] = Observation(source=key, ptinit='parameters.yaml:observations')
            if val:
                observations[key].update(**val)
        self.observations = observations

    def proc_filters(self, check =['LOA', 'LOB', 'LOC', 'LOD']):
        """ Process filters from configuration. """
        for chk in check:
            print(f"Checking filters for {chk}...")
            if hasattr(self, chk) and getattr(self, chk) is not None:
                for filt in self.filters:
                    b = filt.get('band', None)
                    if b is not None:
                        newb = []
                        for b in filt['band']:
                            newb.append(b * getattr(self, chk).unit)
                    filt['band'] = newb
            return
