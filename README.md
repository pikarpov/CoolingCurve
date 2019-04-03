# CoolingCurve

Import, run, and produce specific cooling curves in CIE, which includes the correction for density. 
Electron number dencity is temperature dependent, calculated from CIE ionization fractions of Bryans et al. 2008.
Ion-by-ion tables are taken from Gnat et al. 2011, which were calculated with Cloudy.

cooling.py: the main script containing CoolingCurve class.

run.py: a few examples of different situations to test + function descriptions.

SpeciesList.txt:  default script containing element abundance fractions relative to H.

SpeciesList_gnatSolar.txt: solar abundance fraction tables, taken from Gnat+11; provided for convenience.

IonByIon: ion by ion cooling efficiencies taken from http://wise-obs.tau.ac.il/~orlyg/ion_by_ion/

IonFraction: CIE ionization fractions from Bryans+08 http://adsabs.harvard.edu/abs/2008AAS...212.0302B
