#!/usr/bin/env python

"""
Version of the dengue 4 serotype model implemented for Stochpy
"""
import stochpy
import os
from itertools import cycle
import pylab as P
import numpy as np
from stochpy.modules import InterfaceCain

def plot(t,s,l):
    P.figure()
    s = np.vstack(s).T
    co = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    sy = cycle(['o', '^', '>', '<', 's', '*', '+', '1'])
    print s.shape, len(t)
    for i in range(s.shape[1]):
        P.plot(t, s[:,i], co.next()+sy.next()+'-')
    P.legend(l, loc=0)
    
if __name__=="__main__":
# Loading the model
    smod = stochpy.SSA()
    smod.model_dir = os.getcwd()
    smod.output_dir = os.getcwd()
    smod.Model('Dengue_full.psc')
    
    InterfaceCain.getCainInputfile(smod.SSA, 100, 100, 1)
    with open(os.path.join(stochpy.temp_dir,'cain_in.txt')) as f:
        with open('export.txt','w') as g:
            g.write(f.read())

    smod.DoStochSim(trajectories = 1,mode = 'time',end = 100,  method="direct", IsTrackPropensities=False)
    #smod.DoCainStochSim(endtime=100,frames=10000,trajectories=1, solver="HomogeneousDirect2DSearch",IsTrackPropensities=False)
    smod.GetRegularGrid()
    smod.PlotSpeciesTimeSeries()
    #smod.PlotAveragePropensitiesTimeSeries()
    #smod.PlotAverageSpeciesAutocorrelations()
    
#    t = smod.data_stochsim_grid.time
#    series = smod.data_stochsim_grid.species_means
#    sds = smod.data_stochsim_grid.species_standard_deviations
#    l = smod.data_stochsim.species_labels
#    plot(t,seriegs,l)
#    P.show()
