"""
Import, run, and produce specific cooling curves in CIE, which includes the correction for density. 
Electron number dencity is temperature dependent, calculated from Bryans et al. 2008.
Ion-by-ion tables are taken from Gnat et al. 2011, which were calculated with Cloudy.

CoolingCurve():

DEFAULT INPUTS:
    z = 1: fraction of the metallicity
    name = 'cooling_table_test.txt': how to save the calculated cooling function
    speclist = 'SpeciesList.txt': which SpeciesList to use, which includes
                                  element abundance fractions relative to H

OUTPUT: n_e*n_H*lambda/rho^2

==============================        

CoolingCurve.Plot():

The function is used for plotting the most recently generate cooling function
DEFAULT INPUTS:
    savename = 'None': if specified, then the plot will be saved

-pikarpov
"""

from cooling import CoolingCurve

#Simplest single cooling curve generation & plotting
cc = CoolingCurve(speclist = 'SpeciesList_gnatSolar.txt')
cc.ProduceCurve()
cc.Plot()


#multifunction generation in a loop, with the results plotted in the same figure
zvals = ['1e-3', '1e-2', '5e-2', '1e-1', '5e-1','1e0', '1e1', '1e2']
cc = CoolingCurve()
cc.speclist = 'SpeciesList_gnatSolar.txt'
for z in zvals:
    cc.z = float(z)
    cc.name = 'cooling_table_%s.txt'%z
    cc.ProduceCurve()
    cc.Plot()
