#! /usr/bin/env python


file_list = 'ataalarms,atafx.orig,ataj2000togalactic,atapmetergetspeed,atasetskyfreq,ataant2miriad,atafxsri,atapmetersetaverage,ataant.csh,'\
            'atafxsrimark2,atajdtomjd,atapmetersetspeed,atasettsys,ataantennastate.csh,atagalactictoj2000,atapointer,atasetvelocity,ataants,'\
            'atalistantennalocs,atapp,atageo,ataproblem,atashell.csh,ataantsyscntl,atalistants,ataprog,atasinglepoint,ataantsyscntl~,atagetalarm,'\
            'atasolarsysephem,ataantsyscntl.old,atagetarraypos,atalistarrayorigin,ataps,ataasciistatus,ataquerydb,atasprmonitor,ataastro,atagetazel,'\
            'atalistcatalog,atastatus,ataradecephem,ataautotune,atagetcryopower,atalistdeltat,atastopants,ataautotune~,atagetcryorejtemp,atarebootant,'\
            'ataazeltoradec,atagetcryotemp,atalisteop,atarecoverants,atasys,atabasketweave,atagetdet,ataref,atabcknds,atagetdetdbm,atalistgpib,atasyslog,'\
            'atagetfocus,atareloadpm,atabimaweather,atagetholdtime,atalistitems,ataresetantenna,atataitoiso8601,atagetitem,ataresetserver,atacam,atalistsats,'\
            'atatempref,atacamgetptz,atagetlnas,ataresetservos,atalistservers,atatestants,atacamsetpan,atagetlo,atasacapture,atatestpointer,atalistsolarsys,'\
            'atacamsettilt,atagetmetadataitem,atasacontrol,atatime,atagetpams,atalisttles,atasaspectrum,atacamsetzoom,atatrack,atagetpmoffsets,atalnaoff,'\
            'atasatephem,atacatalogephem,atatrackephem,atagetrack,atalnaon,atasatinfo,atacheck,atagetradec,atasats,atatransitscan,atagetrpa,atalockserver,'\
            'ataunlockserver,atacircularizedtrackephem,atasbcoffon,atacmds,atagetrpagain,atamainsystime,atasefd,ataupdatecatalog,ataservostatus,atacollision,'\
            'atagetsampleperiod,atamenu,ataupdatetles,atadbampconfig,atagetservocurrent,atasetalarm,atautil,atadbmref,atagetskyfreq,atamisc,atasetantsoff,'\
            'atavideocapture,atadddmmsstodegrees,atagettle,atamjdtoiso8601,atasetantson,atagettsys,atavideoclient,atadegreestodddmmss,atagetwrappot,atamjdtojd,'\
            'atasetazel,atagpibcmnd,atawaituntilazel,atadescribeephem,atagui,atasetfocus,atamovie,atasetholdtime,atawaituntiltime,atadish,atahelpers.rb,atasetlna,'\
            'atahelpers.rb~,atamq,atasetlo,atawaituntiltracking,ataemergency,atahhmmsstohours,atamultiataps,atanvsslist,atasetmetadataitem,atawatchdog,ataephem,'\
            'atahotcold,ataobsinfo,atasetpams,atawave,ataofflimit,atasetpams6670.csh,ataweather,atafield,atahourstohhmmss,ataonoff,atasetpmoffsets,ataweather.csh,'\
            'atafindencoders,ataifp,ataweathermon.csh,ataopticalpointer,atasetrack,atawindsockcontrol,atafixedephem,ataigsephem,atasetrpa,atawindsockstatus,'\
            'ataiso8601toatatai,atapeakup,atawrapephem,atafx,atapmeterconfigure,atasetrpagain,atawtr.csh,atafx_launch,ataiso8601tomjd,atapmetergetaverage,atafxnohistory,'\
            'atapmetergetpower,atasetsampleperiod'


class ATAControl:
    def __init__(self, from_file=False):
        self.files = []
        if from_file:
            with open('atacontrol.txt', 'r') as fp:
                for line in fp:
                    data = line.split()
                    self.files += data
        else:
            self.files = file_list.split(',')

    def show(self, expr):
        for line in self.files:
            if expr in line and not line.lower().endswith('.bat'):
                print(line)

    def get_file_list(self):
        trunc = []
        for line in self.files:
            if not line.lower().endswith('.bat'):
                trunc.append(line)
        return trunc
    

if __name__ == '__main__':
    import sys
    ata_control = ATAControl()
    ata_control.show(sys.argv[1])