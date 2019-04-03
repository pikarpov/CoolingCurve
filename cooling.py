'''
Generate a specific cooling curve for a system in CIE,
which the correction for density.

Class Parameters:
    z = 1: fraction of the metallicity
    name='cooling_table_test.txt': how to save the calculated cooling function
    speclist='SpeciesList.txt': which SpeciesList to use, which includes
                                element abundance fractions relative to H

Output: n_e*n_H*lambda/rho^2

-pikarpov
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from scipy import interpolate

#my preferred plotting parameters, which you should feel free to change as you please
params = {'axes.labelsize': 25, 'legend.fontsize': 21, 'xtick.labelsize': 21,'ytick.labelsize': 21,
          'axes.linewidth': 1, 'lines.linewidth': 1.5,
          'xtick.major.width': 1,'ytick.major.width': 1,'xtick.minor.width': 1,'ytick.minor.width': 1,
          'xtick.major.size': 4,'ytick.major.size': 4,'xtick.minor.size': 3,'ytick.minor.size': 3, 
          'text.usetex': True}
mpl.rcParams.update(params)

class CoolingCurve:
    
    #set some default quantities, that can be altered at the function call
    def __init__(self, z=1, name='cooling_table_test.txt', speclist='SpeciesList.txt'): 
        self.z = z
        self.name = name
        self.speclist = speclist
        self.x = []
        self.comb = []
        self.ax = 'initialize'
        
    def ProduceCurve(self):
        #===============================
        #The function calculates the specific cooling curve
        #===============================
        
        #determine which elements are in SpeciesList
        names = np.moveaxis(np.genfromtxt(self.speclist, dtype=str), -1,0).tolist()
        elements = ['h']+names[0]
        elements = [i.capitalize() for i in elements]
        
        #get the abundance ratios, and realtive H & e number densities
        ratio, nh, ne = self.Ratios()
        for i in range(len(elements)):
            #import cooling efficiency of the element's ions
            ion_filename = r'IonByIon/%s.txt'%elements[i]  
            n = self.Count(ion_filename)
            vals = np.moveaxis(np.genfromtxt(ion_filename, skip_header=n), -1,0)
            
            #import CIE ionization fractions for each ion of the looped element
            abund = np.moveaxis(np.genfromtxt(r'IonFraction/%s.csv'%elements[i], delimiter=',', skip_header=1), -1,0)

            #convert from log format
            abund[0] = [10**j for j in abund[0]]
            self.x = vals[0]
            
            cie = np.zeros([len(self.x)])
            nefrac = np.zeros([len(self.x)])
            
            #iterate over all ions
            for j in range(1,len(vals)-1):
                #spline interpolation of the fractional abundance table
                spl = interpolate.splrep(abund[0], abund[j]) 
                
                #generate new tables with corresponding fractional abundances
                abundy = interpolate.splev(self.x, spl) 
                
                #add to the cooling curve of the element
                cie += 10**(-abundy)*vals[j]
                
                #add the number of electrons contributed by the ion
                nefrac += 10**(-abundy)*(j-1)
            
            #correction to the total number of electrons in the system 
            #to get the number of free electrons at each T bin
            
            nefrac = [k/(len(vals)-3) for k in nefrac]
            
            #Note: if one wants to get pure cooling functions, then uncomment bellow
            #nefrac = np.full(np.size(nefrac), 1); nh = 1; ne = 1
            
            if i==0: self.comb = ratio[i]*cie*nh*ne*nefrac
            else: self.comb += ratio[i]*cie*nh*ne*nefrac

            self.write_output()
        return


    def Ratios(self):
        #===============================
        #The function reads in the fractional abundances from SpeciesList
        #and calculate the relative H number density (nh)
        #and relative electron number density (ne)
        #===============================
        
        #read-in the fractions from the SpeciesList
        vals = np.moveaxis(np.genfromtxt(self.speclist), -1,0)
        
        #apply metallicity (z) correction if needed; default z=1
        vals[-1][1:] = vals[-1][1:]*self.z
        
        #add H contribution, which is always 1, since everything is expressed as relative to H
        ratio = np.concatenate(([1],vals[-1]))
        
        #calculate the relative H number density
        u = 1.6605e-24
        ma = np.concatenate(([1],vals[2]))
        frac = 0
        for i in range(len(ratio)):
            frac += ma[i]*u*ratio[i]
        nh = 1/frac #mass is simply set to 1g, to get the "relative" quantity with correct units
        
        #calculate the total number of electrons available (both bound and free)
        for i in range(len(ratio)):
            if i==0: ne = ratio[i]*nh
            else: ne += ratio[i]*nh*vals[1][i-1]

        return ratio, nh, ne
    
    
    def write_output(self):
        #===============================
        #The functions writes output to a specified file
        #Default: name = 'cooling_table_test.txt'
        #===============================
        
        f = open(self.name, 'w+')
        
        f.write('%d\n'%len(self.x))
        for i in range(len(self.x)):
            f.write('%.2e\t%.2e\n'%(self.x[i], self.comb[i]))
        f.close()    


    def Plot(self, savename=None):
        #===============================
        #The function plots the final effective cooling curve(s)
        #and saves the output if "savename" was specified upon calling Plot()
        #===============================
        
        #if it is the first plot, then create a figure; others will be overplotted
        if self.ax=='initialize': self.ax = plt.figure(figsize=(9,6)).add_subplot(1,1,1)
        
        #basic plotting parameters
        self.ax.loglog(self.x, self.comb, basex=10, basey=10, marker = 'None', linestyle='-')
        self.ax.grid(True)
        self.ax.set_xlabel(r'$\rm T\;[K]$')
        self.ax.set_ylabel(r'$\rm n_e n_H \lambda/\rho^2 \; [erg \; cm^{-3} \;s^{-1}]$')
        plt.tight_layout()
        
        #save the figure if savename was specified
        if savename!=None: plt.savefig(savename)

        
    def Count(self, filename):
        #===============================
        #The function counts the number of string lines to skip 
        #in the datafiles so the tables can be read. Since the tables
        #start at T=1e4, that's the first line to be imported.
        #===============================
        
        f = open(filename, 'r')
        n = 0
        for line in f:
            if line.find("1.00e+04")!=-1:
                break
            n+=1
        f.close()
        return n