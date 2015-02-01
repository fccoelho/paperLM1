#!/usr/bin/env python
"""
Version of the dengue 4 serotype model implemented for Stochpy
"""
import stochpy
import os
from itertools import cycle
import pylab as P
import numpy as np
# from stochpy.modules import InterfaceCain


def plot(t, s, l):
    P.figure()
    s = np.vstack(s).T
    co = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    sy = cycle(['o', '^', '>', '<', 's', '*', '+', '1'])
    print s.shape, len(t)
    for i in range(s.shape[1]):
        P.plot(t, s[:, i], co.next() + sy.next() + '-')
    P.legend(l, loc=0)

def agg_by_type(s, l):
    """
    return series aggregated by serotype
    :param s: series
    :param l: labels
    """
    s = np.vstack(s).T
    i1 = [l.index(i) for i in l if i.endswith('1') and i.startswith('I')]
    s1 = np.zeros(s.shape[0])
    for i in i1:
        s1 += s[:, i]

    i2 = [l.index(i) for i in l if i.endswith('2') and i.startswith('I')]
    s2 = np.zeros(s.shape[0])
    for i in i2:
        s2 += s[:, i]

    i3 = [l.index(i) for i in l if i.endswith('3') and i.startswith('I')]
    # print "I3: ", i3
    s3 = np.zeros(s.shape[0])
    for i in i3:
        s3 += s[:, i]

    i4 = [l.index(i) for i in l if i.endswith('4') and i.startswith('I')]
    s4 = np.zeros(s.shape[0])
    # print "I4: ", i4
    for i in i4:
        s4 += s[:, i]

    return s1, s2, s3, s4


def plot_4_types(t, s, l):
    """
    Plot all infectives for each serotype
    :param t: Times
    :param s: series
    :param l: labels
    """
    # print l
    P.figure()

    s1, s2, s3, s4 = agg_by_type(s, l)

    P.plot(t, s1, label='1')

    P.plot(t, s2, label='2')

    P.plot(t, s3, label='3')

    P.plot(t, s4, label='4')
    P.legend(loc=0)
    P.xlabel('time')
    P.ylabel('individuals')
    P.savefig('4types.png', dpi=300)

def plot_xcorr(s1, s2):
    P.figure("Cross Correlation")


if __name__ == "__main__":
    # Loading the model
    smod = stochpy.SSA(IsInteractive=False)
    smod.model_dir = os.getcwd()
    smod.output_dir = os.getcwd()
    smod.Model('Dengue_full4.psc')


    # InterfaceCain.getCainInputfile(smod.SSA, 100, 100, 1)
    #    with open(os.path.join(stochpy.temp_dir,'cain_in.txt')) as f:
    #        with open('export.txt','w') as g:
    #            g.write(f.read())

    smod.DoStochSim(trajectories=1, mode='time', end=1000, method="direct", IsTrackPropensities=True)
    #smod.DoCainStochSim(endtime=200,frames=10000,trajectories=1, solver="HomogeneousDirect2DSearch",IsTrackPropensities=False)
    smod.GetRegularGrid()
    smod.PlotSpeciesTimeSeries()
    stochpy.plt.savefig(os.path.join(smod.output_dir, 'full4.png'), dpi=300)
    #smod.PlotAveragePropensitiesTimeSeries()
    #smod.PlotAverageSpeciesAutocorrelations()

    smod.PlotWaitingtimesDistributions('Birth', colors=['y'],
                                       linestyle='None', marker='v')
    smod.PlotWaitingtimesDistributions('Inf_I2341', colors=['r'],
                                       linestyle='None', marker='^')
    smod.PlotWaitingtimesDistributions('Inf_I1234', colors=['purple'],
                                       linestyle='None', marker='o')
    smod.PlotWaitingtimesDistributions('Inf_I3', colors=['b'],
                                       linestyle='None', marker='o')
    t = smod.data_stochsim_grid.time
    series = smod.data_stochsim_grid.species_means
    sds = smod.data_stochsim_grid.species_standard_deviations
    l = smod.data_stochsim.species_labels
    plot_4_types(t, series, l)
    #    plot(t,series,l)
    #    P.show()
    stochpy.plt.show()
