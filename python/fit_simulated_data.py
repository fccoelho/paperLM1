"""
This script applies the inference to simulated data to validate
the methodology.
"""
import pyximport

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
from sqlalchemy import create_engine
import seaborn as sns
# import numba
# from numba import jit


beta = 1.0  # Transmission coefficient
b1 = 19.9
b0 = 1.5  # low beta during winters
eta = .0  # infectivity of asymptomatic infections relative to clinical ones. FIXED

mu = 0  # nat/mortality rate
m = 1e-6  # influx of cases
tau = 1.  # recovery rate. FIXED

N = 1  # Population of Rio

inicio = 0
# ~ Ss = {0: s0, 1: s1, 2: s2}  # Multiple Ss map


# Initial conditions
inits = np.array([0.999, 0.001, 0.0])  # initial values for state variables.


# @jit
def sir(y, t, *pars):
    '''ODE model'''
    S, I, R = y
    s0, m = pars
    beta = 0 if t > 728 else iRt(t)*tau / 0.0621

    lamb = beta * (I + m) * S

    return np.array([-lamb,  # dS/dt
                     lamb - tau * I,  # dI/dt
                     tau * I,  # dR/dt
                     ])


# @jit
def jac(y, t, *pars):
    """
    Jacobian of the model
    :param y: State of the system
    :param t: current time
    :param pars: tuple of parameters
    :return: the jacobian matrix at time t
    """
    S, I, R = y
    s0, m = pars
    beta = 0 if t > 728 else iRt(t)*tau / 0.0621
    return np.array([[-(I + m) * beta, -S * beta, 0],
                     [(I + m) * beta, S * beta - tau, 0],
                     [0, tau, 0]])


# @jit
def model(theta):
    """
    Wrap-up function which is called at every step of the MCMC by BIP's FitModel.
    :param theta: tuple of parameter values
    :return: the output of the model
    """
    # setting parameters
    s0, m = theta

    t0 = 0
    tf = fim - inicio
    Y = np.zeros((tf - t0, 3))

    # print(iRt(np.arange(t0, tf, 1)))
    # Initial conditions

    inits[0] = s0  # Define S0
    inits[-1] = N - sum(inits[:2])  # Define R(0)
    Y[t0:tf, :] = odeint(sir, inits, np.arange(t0, tf, 1), args=(s0, m), Dfun=jac)  # ,tcrit=tcrit)
    # inits = Y[-1, :]
    Y[t0:tf, 1] = Y[t0:tf, 1] / 1  # Adjusting the output to just the reported I to compare with data
    return Y


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
    pop = pd.read_csv("../DATA/pop_rio_1980-2012.csv", header=0, delimiter=',', index_col=0)

    dates = [datetime.datetime.strptime(d, "%Y-%m-%d") for d in data.start]
    pop_d = np.array([pop.loc[d.year] for d in dates])  # population on each date

    eday = len(data) if eday is None else eday

    incidence = data.cases  # daily incidence
    # Converting incidence to Prevalence
    dur = 1./ tau  # infectious period
    rawprev = np.convolve(incidence, np.ones(dur), 'same')
    rawprev.shape = rawprev.size, 1
    rawprev /= pop_d

    # plot_data(data, dates, incidence, rawprev)
    # Doing moving average of order mao
    if mao > 1:
        sw = np.ones(mao, dtype=float) / mao  # smoothing window
        prev = np.convolve(rawprev, sw, 'same')  # Smoothing data (ma)
    else:
        prev = rawprev
    # sw = np.ones(6, dtype=float) / 6  #smoothing window
    # rt_smooth = np.convolve(data.Rt2, sw, 'same')
    Rt = fix_rt(data.Rt)
    prev.shape = (prev.size,)
    d = {'time': dates, 'I': np.nan_to_num(prev), 'Rt': Rt}
    return d


def get_simulated_data( pars):
    y = model(pars)

    return y


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


def read_data(dbname):
    """
    Returns dataframes from tables of the sqlite database with the results of the inference
    :param dbname: sqlite3 database name
    :return: list of pandas dataframes
    """
    eng = create_engine('sqlite:///{}'.format(dbname))
    df_data = pd.read_sql_table('data', eng, index_col=['time'], parse_dates={'time': '%Y-%m-%d %H:%M:%S'})
    df_pt = pd.read_sql_table('post_theta', eng, index_col=['time'], parse_dates={'time': '%Y-%m-%d %H:%M:%S'})
    df_series = pd.read_sql_table('series', eng, index_col=['time'], parse_dates={'time': '%Y-%m-%d %H:%M:%S'})
    return df_data, df_pt, df_series

def plot_data(dbname):
    """
    Plot the data loaded from disk
    :param data:
    :param dates:
    :param incidence:
    :param rawprev:
    :return:
    """
    data, theta, series = read_data(dbname + '.sqlite')
    print(series.columns)
    P.plot(pd.to_datetime(series.index), series.I, alpha=0.7, label='Incidence')

    P.setp(P.gca().xaxis.get_majorticklabels(), rotation=45)
    P.grid()
    P.legend()
    # P.figure()
    # P.plot(dates, data.Rt, label=r'$R_t$')
    # P.plot(dates, data.lwr, 'r-.')
    # P.plot(dates, data.upr, 'r-.')
    # P.setp(P.gca().xaxis.get_majorticklabels(), rotation=45)
    P.show()


# # running the analysys
if __name__ == "__main__":
    # dt = prepdata('aux/data_Rt_dengue_big.csv', 0, 728, 1)
    dt = prepdata('../DATA/data_Rt_dengue_complete.csv', 0, 971, 1)


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

           dt['time'].index(datetime.datetime(2010, 10, 17)),  # Start of the 2011 epidemic
           dt['time'].index(datetime.datetime(2011, 8, 28)),  # Start of the 2012 epidemic
           dt['time'].index(datetime.datetime(2012, 11, 11)),  # Start of the 2013 epidemic
           ]
    tfs = t0s[1:] + [len(dt['time'])]
    tfs = [dt['time'].index(datetime.datetime(1996, 7, 29)),  # end of the 1996 epidemic
           dt['time'].index(datetime.datetime(1998, 10, 12)),  # end of the 1998 epidemic
           dt['time'].index(datetime.datetime(1999, 8, 23)),  # end of the 1999 epidemic
           dt['time'].index(datetime.datetime(2000, 10, 2)),  # end of the 2000 epidemic
           dt['time'].index(datetime.datetime(2001, 9, 10)),  # end of the 2001 epidemic
           dt['time'].index(datetime.datetime(2002, 9, 2)),  # end of the 2002 epidemic
           dt['time'].index(datetime.datetime(2006, 7, 31)),  # end of the 2006 epidemic
           dt['time'].index(datetime.datetime(2007, 8, 27)),  # end of the 2007 epidemic
           dt['time'].index(datetime.datetime(2008, 9, 1)),  # end of the 2008 epidemic
           dt['time'].index(datetime.datetime(2009, 9, 28)),  # end of the 2009 epidemic

           dt['time'].index(datetime.datetime(2011, 8, 28)),  # end of the 2011 epidemic
           dt['time'].index(datetime.datetime(2012, 11, 11)),  # end of the 2012 epidemic
           dt['time'].index(datetime.datetime(2013, 8, 25)),  # end of the 2013 epidemic
           ]

    # Interpolated Rt
    # iRt = interp1d(np.arange(dt['Rt'].size), np.array(dt['Rt']), kind='linear', bounds_error=False, fill_value=0)
    # print([iRt(i) for i in range(0, 971)])

    pnames = ['S', 'I', 'R']

    wl = dt['I'].shape[0]  # window length
    nw = len(dt['time']) / wl  # number of windows
    tf = wl * nw  # total duration of simulation

    for inicio, fim in list(zip(t0s, tfs))[12:13]:  # Optional Slice to allow the fitting of a subset of years
        print(inicio, fim)
        dt = prepdata('../DATA/data_Rt_dengue_complete.csv', inicio, fim, 1)
        # get simulated data
        inits = [1 - dt['I'][0], dt['I'][0], 0]
        # print(iRt(np.arange(0, fim-inicio, 1)))

        # Interpolated Rt
        iRt = interp1d(np.arange(dt['Rt'].size), np.array(dt['Rt']), kind='linear', bounds_error=False, fill_value=0)
        y = get_simulated_data([.0621, 1e-6])
        dt['I'] = y[:, 1]
        ano = dt['time'][0].year
        mes = dt['time'][0].month
        modname = "Sim_DengueS{}_{}".format(ano, mes)
        tnames = ['s_{}_{}'.format(ano, mes), 'm']

        nt = len(tnames)
        pnames = ['S', 'I', 'R']
        nph = len(pnames)
        wl = fim - inicio
        # Taking a a look at the simulated S curve before inference
        P.plot(dt['time'], y[:, 0], label='S')
        P.figure()
        P.plot(dt['time'], y[:, 1], label='I')
        P.plot(dt['time'], dt['I'], 'ro', label='data')
        # P.plot(dt['time'], dt['Rt'], 'r^', label='Rt')
        # print([iRt(t) for t in np.arange(0, fim-inicio)])
        # P.plot(dt['time'], [iRt(t) for t in np.arange(0, fim-inicio)], label='iRt')
        P.legend(loc=0)
        P.show()
        nw = 1

        tpars = [(1, 1), (1e-6, 1e-8)]
        tlims = [(0.03, .1), (0, 5e-6)]
        del dt['Rt']
        dt2 = copy.deepcopy(dt)
        # print inits

        F = FitModel(5000, model, inits, fim - inicio, tnames, pnames,
                     wl, nw, verbose=1, burnin=1100, constraints=[])
        F.set_priors(tdists=[st.beta, st.norm],
                     tpars=tpars,
                     tlims=tlims,
                     pdists=[st.beta] * nph, ppars=[(1, 1)] * nph, plims=[(0, 1)] * nph)

        F.run(dt, 'DREAM', likvar=1e-8, likfun='Normal', pool=False, ew=0, adjinits=True, dbname=modname, monitor=['I', 'S'], initheta=[0.0621, 1e-6])
        # ~ print(F.AIC, F.BIC, F.DIC)
        # print (F.optimize(data=dt2, p0=[.0621, 1e-6, 1.0], optimizer='scipy',tol=1e-55, verbose=1, plot=1))
        F.plot_results(['S', 'I'], dbname=modname, savefigs=1)
        P.figure()
        P.figure()
        P.plot(dt['time'], y[:, 1], label='I')
        plot_data(modname)
        P.clf()
        P.clf()
