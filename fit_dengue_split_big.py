#
# Multi-year fit of influenza epidemics
#
import pyximport;

pyximport.install(pyimport=False)
from BIP.Bayes.Melding import FitModel
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import scipy.stats as st
import numpy as np
import pylab as P
import copy
from collections import defaultdict
import datetime
import sys
import pandas as pd
import numba
from numba import jit


beta = 1.0  # Transmission coefficient
b1 = 19.9
b0 = 1.5  # low beta during winters
eta = .0  # infectivity of asymptomatic infections relative to clinical ones. FIXED
epsilon = 2.8  # latency rate
mu = 0  # nat/mortality rate
m = 1e-6 # influx of cases
tau = 1  # recovery rate. FIXED

N = 1  # Population of Rio

inicio=0
#~ Ss = {0: s0, 1: s1, 2: s2}  # Multiple Ss map


# Initial conditions
inits = np.array([0.999, 0.001, 0.0])  # initial values for state variables.


#@jit#('f8[:](f8[:],f8,f8,f8)')
def sir(y, t, *pars):
    '''ODE model'''
    S, I, R = y
    s0, m, tau = pars
    beta = 0 if t > 728 else iRt(t) * tau/s0
    # print S, I, beta
    
    lamb = beta * (I+m) * S

    return np.array([-lamb,  #dS/dt
            lamb - tau * I,  #dI/dt
            tau * I,  #dR/dt
    ])

#@jit#('f8[:,:](f8[:],f8,f8,f8)')
def jac(y, t, *pars):
    S, I, R = y
    s0, m, tau = pars
    beta = 0 if t > 728 else iRt(t) * tau/s0
    return np.array([[-(I + m)*beta, -S*beta, 0],
           [(I + m)*beta, S*beta - tau, 0],
           [0, tau, 0]])

def model(theta):
    # setting parameters
    s0, m, tau = theta

    t0 = 0
    tf = fim
    Y = np.zeros((tf-t0, 3))

    
    # Initial conditions
    #print inits, m
    inits[0] = s0  # Define S0
    inits[-1] = N - sum(inits[:2])  # Define R(0)
    Y[t0:tf, :] = odeint(sir, inits, np.arange(t0, tf, 1), args=(s0, m, tau), Dfun=jac)  #,tcrit=tcrit)
    #inits = Y[-1, :]
    Y[t0:tf, 1] = Y[t0:tf, 1] / 1 # Adjusting the output to just the reported I to compare with data
    return Y


def plot_data(data, dates, incidence, rawprev):
    """
    Plot the data loaded from disk
    :param data:
    :param dates:
    :param incidence:
    :param rawprev:
    :return:
    """
    P.plot(dates, incidence, label='Incidence')
    P.plot(dates, rawprev, label='Prevalence')
    P.setp(P.gca().xaxis.get_majorticklabels(), rotation=45)
    P.grid()
    P.legend()
    P.figure()
    P.plot(dates, data.Rt, label=r'$R_t$')
    P.plot(dates, data.lwr, 'r-.')
    P.plot(dates, data.upr, 'r-.')
    P.setp(P.gca().xaxis.get_majorticklabels(), rotation=45)
    P.show()


def prepdata(fname, sday=0, eday=None, mao=7):
    """
    Prepare the data for the analysis.

    Parameters:
    file: path to data
    sday: Starting day of the inference
    eday: final day
    mao: Moving average's Order
    """
    data = pd.read_csv(fname, header=0, delimiter=',', skiprows=[1, 2, 3], parse_dates=True)
    # slicing to the desired period
    data = data[sday:eday]
    pop = pd.read_csv("pop_rio_1980-2012.csv", header=0, delimiter=',', index_col=0)

    dates = [datetime.datetime.strptime(d, "%Y-%m-%d") for d in data.start]
    pop_d = np.array([pop.loc[d.year] for d in dates])  # population on each date

    eday = len(df) if eday is None else eday
    # print data.dtype.names
    incidence = data.cases  # daily incidence
    # Converting incidence to Prevalence
    dur = 1. / tau  # infectious period
    rawprev = np.convolve(incidence, np.ones(dur), 'same')
    rawprev.shape = rawprev.size, 1
    rawprev /= pop_d

    # plot_data(data, dates, incidence, rawprev)
    # Doing moving average of order mao
    if mao > 1:
        sw = np.ones(mao, dtype=float) / mao  #smoothing window
        prev = np.convolve(rawprev, sw, 'same')  #Smoothing data (ma)
    else:
        prev = rawprev
    # sw = np.ones(6, dtype=float) / 6  #smoothing window
    # rt_smooth = np.convolve(data.Rt2, sw, 'same')
    Rt = fix_rt(data.Rt)
    prev.shape = (prev.size,)
    d = {'time': dates, 'I': np.nan_to_num(prev), 'Rt': Rt}
    return d

@np.vectorize
def fix_rt(rt):
    """
    Replace both NANs and INFs by zero
    :param rt: reproductive number, scalar
    :return: fixed rt
    """
    if np.isnan(rt):
        return 0
    elif np.isinf(rt):
        return 0
    else:
        return rt

# # running the analysys
if __name__ == "__main__":
    dt = prepdata('aux/data_Rt_dengue_big.csv', 0, 728, 1)

    
    # Defining start and end of the simulations
    t0s = [0,  # Start of the 1996 epidemic
           dt['time'].index(datetime.datetime(1997, 12, 15)),  # Start of the 1998 epidemic
           dt['time'].index(datetime.datetime(1998, 12, 21)),  # Start of the 1999 epidemic
           dt['time'].index(datetime.datetime(1999, 12, 13)),  # Start of the 2000 epidemic
           dt['time'].index(datetime.datetime(2000, 12, 18)),  # Start of the 2001 epidemic
           dt['time'].index(datetime.datetime(2001, 9, 10)),  # Start of the 2002 epidemic
           dt['time'].index(datetime.datetime(2005, 8, 15)),  # Start of the 2006 epidemic
           dt['time'].index(datetime.datetime(2006, 9, 25)),  # Start of the 2007 epidemic
           dt['time'].index(datetime.datetime(2007, 8, 27)),  # Start of the 2008 epidemic
           dt['time'].index(datetime.datetime(2008, 9, 1)),  # Start of the 2009 epidemic
    ]
    tfs = t0s[1:] + [len(dt['time'])]
    tfs = [ dt['time'].index(datetime.datetime(1996, 7, 29)),  # end of the 1996 epidemic
            dt['time'].index(datetime.datetime(1998, 10, 12)),  # end of the 1998 epidemic
            dt['time'].index(datetime.datetime(1999, 8, 23)),  # end of the 1999 epidemic
            dt['time'].index(datetime.datetime(2000, 10, 2)),  # end of the 2000 epidemic
            dt['time'].index(datetime.datetime(2001, 9, 10)),  # end of the 2001 epidemic
            dt['time'].index(datetime.datetime(2002, 9, 2)),  # end of the 2002 epidemic
            dt['time'].index(datetime.datetime(2006, 7, 31)),  # end of the 2006 epidemic
            dt['time'].index(datetime.datetime(2007, 8, 27)),  # end of the 2007 epidemic
            dt['time'].index(datetime.datetime(2008, 9, 1)),  # end of the 2008 epidemic
            dt['time'].index(datetime.datetime(2009, 9, 28)),  # end of the 2009 epidemic
    ]
    # print tfs
    # Interpolated Rt
    iRt = interp1d(np.arange(dt['Rt'].size), np.array(dt['Rt']), kind='linear', bounds_error=False, fill_value=0)

    P.plot(dt['Rt'],'*')
    P.plot(np.arange(0, 728, .2), [iRt(t) for t in np.arange(0, 728, .2)])
    #print type(dt['Rt'])
    #print [iRt(t) for t in np.arange(0, 728, .2)]


    pnames = ['S', 'I', 'R']

    wl = dt['I'].shape[0]  #window length
    nw = len(dt['time']) / wl  #number of windows
    tf = wl * nw  #total duration of simulation

    #~ print calc_R0s()
    inicio = t0s[0]
    fim = tfs[-1]
    #~ y = model([.999*N, 1e-6, 1])
    #print y
    #~ P.figure()
    #~ P.plot(dt['I'], '*')
    #~ P.plot(y[:, 1])
    #~ top = y[:, 1].max()
    #~ P.vlines(t0s,0,top, colors='g')
    #~ P.vlines(tfs,0,top, colors='r')
    #~ P.legend([pnames[1]])
    #~ P.show()
    #Priors and limits for all countries




    for inicio, fim in zip(t0s, tfs)[3:4]: # Slice to force start from a different point
        dt = prepdata('data_Rt_dengue_big.csv', inicio, fim, 1)
        # Interpolated Rt
        iRt = interp1d(np.arange(dt['Rt'].size), np.array(dt['Rt']), kind='linear', bounds_error=False, fill_value=0)
        ano = dt['time'][0].year
        mes = dt['time'][0].month
        modname = "DengueS{}_{}".format(ano, mes)
        tnames = ['s_{}_{}'.format(ano, mes), 'm', 'tau']
        

        nt = len(tnames)
        pnames = ['S', 'I', 'R']
        nph = len(pnames)
        wl = fim - inicio
        nw = 1

        tpars = [(1, 1), (0, 5e-6),(.9999,.0002)]
        tlims = [(0, 1), (0, 5e-6),(.9999,1.0001)]

        inits = [1 - dt['I'][0], dt['I'][0], 0]
        dt2 = copy.deepcopy(dt)
        #print inits

        F = FitModel(10000, model, inits, fim-inicio, tnames, pnames,
                     wl, nw, verbose=1, burnin=5000, constraints=[])
        F.set_priors(tdists=[st.beta, st.uniform, st.uniform],
                     tpars=tpars,
                     tlims=tlims,
                     pdists=[st.beta] * nph, ppars=[(1, 1)] * nph, plims=[(0, 1)] * nph)
        
        
        F.run(dt, 'DREAM', likvar=1e-10, pool=False, ew=0, adjinits=True, dbname=modname, monitor=['I', 'S'])
        #~ print F.AIC, F.BIC, F.DIC
        #print F.optimize(data=dt,p0=[s0,s1,s2], optimizer='scipy',tol=1e-55, verbose=1, plot=1)
        F.plot_results(['S', 'I'], dbname=modname, savefigs=1)
        P.clf()
        P.clf()
