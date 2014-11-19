#!/usr/bin/env python

"""
Version of the dengue 4 serotype model implemented for Stochpy
"""
import stochpy
import os


smod = stochpy.SSA()
smod.model_dir = os.getcwd()
smod.output_dir = os.getcwd()
smod.Model('Dengue_full.psc')

#smod.DoStochSim(trajectories = 10,mode = 'time',end = 100, epsilon=0.03, method="TauLeaping", IsTrackPropensities=False)
smod.DoCainStochSim(endtime=100,frames=10000,trajectories=10, solver="HomogeneousDirect2DSearch",IsTrackPropensities=False)
smod.GetRegularGrid()
smod.PlotAverageSpeciesTimeSeries()
#smod.PlotAveragePropensitiesTimeSeries()