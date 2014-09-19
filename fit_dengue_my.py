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


beta = 1.0  # Transmission coefficient
b1 = 19.9
b0 = 1.5  # low beta during winters
eta = .0  # infectivity of asymptomatic infections relative to clinical ones. FIXED
epsilon = 2.8  # latency rate
mu = 0  # nat/mortality rate

tau = 1  # recovery rate. FIXED

N = 1  # Population of Rio
s0 = 0.9999*N  # fraction of susceptibles at the beginning of first epidemic
s1 = 0.999*N;
s2 = 0.999*N;

Ss = {0: s0, 1: s1, 2: s2}  # Multiple Ss map


# Initial conditions
inits = np.array([s0, 100, 0.0])  # initial values for state variables.

def model(theta):
    # setting parameters
    s0, s1, s2 = theta
    Ss = {0: s0, 1: s1, 2: s2}  # Multiple Ss map

    # ~ b1=19.9
    def sir(y, t):
        '''ODE model'''
        S, I, R = y
        
        beta = 0 if t > 242 else  iRt(t) * tau/ N#(Ss[ycode[int(t)]])

        #~ if (calc_R0s()<1.1).any():
        #~ print "flat"
        lamb = (beta * I * S)
        return [-lamb,  #dS/dt
                lamb - tau * I,  #dI/dt
                tau * I,  #dR/dt
        ]

    Y = np.zeros((wl, 3))
    # Initial conditions

    for i in range(len(Ss)):
        t0 = t0s[i]
        tf = tfs[i]
        #print t0,tf
        if i>0:
            #~ inits[1] = Y[t0-1,1]
            #~ print inits
            inits[1] = dt['I'][t0-1] if N-Ss[i] > dt['I'][t0-1] else N-Ss[i]
            #~ print inits
        inits[0] = Ss[i];  # Define S0
        inits[-1] = N - sum(inits[:2])  # Define R(0)
        Y[t0:tf, :] = odeint(sir, inits, np.arange(t0, tf, 1))  #,tcrit=tcrit)
        #inits = Y[-1, :]

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
    data = pd.read_csv(fname, header=0, delimiter=',', skiprows=[1, 2 , 3], parse_dates=True)

    dates = [datetime.datetime.strptime(d, "%Y-%m-%d") for d in data.start]
    sday = sday  # starting day for the fitting
    eday = len(df) if eday is None else eday
    # print data.dtype.names
    incidence = data.cases[sday:eday]  # daily incidence
    # Converting incidence to Prevalence
    dur = 1. / tau  # infectious period
    rawprev = np.convolve(incidence, np.ones(dur), 'same')
    rawprev /= 6e6

    #~ P.plot(dates, incidence, label='Incidence')
    #~ P.plot(dates, rawprev, label='Prevalence')
    #~ P.setp(P.gca().xaxis.get_majorticklabels(), rotation=45)
    #~ P.grid()
    #~ P.legend()
    #~ P.figure()
    #~ P.plot(dates, data.Rt, label=r'$R_t$')
    #~ P.plot(dates, data.lwr, 'r-.')
    #~ P.plot(dates, data.upr, 'r-.')
    #~ P.setp(P.gca().xaxis.get_majorticklabels(), rotation=45)
    P.show()
    # Doing moving average of order mao
    if mao > 1:
        sw = np.ones(mao, dtype=float) / mao  #smoothing window
        prev = np.convolve(rawprev, sw, 'same')  #Smoothing data (ma)
    else:
        prev = rawprev

    # assert len(prev) == len(dates)
    d = {'time': dates, 'I': np.nan_to_num(prev), 'Rt': data.Rt}
    return d


def fill_in_missing_weeks(data):
    """
    Fills in the missing weeks in the data set
    missing prevalence is set to NaN 
    """
    d2 = defaultdict(lambda: [])
    for k, v in data.items():
        newd = []
        for i, l in enumerate(v):
            w, d, p = l
            if i > 0:
                ew = (newd[-1][0] % 52) + 1  # expected week number
                ed = newd[-1][1] + datetime.timedelta(7)  # expected date
                while w != ew:
                    newd.append([ew, ed, 1e-1])
                    ed = newd[-1][1] + datetime.timedelta(7)
                    ew = (newd[-1][0] % 52) + 1
                    p = np.nan
            newd.append([w, d, p])
        d2[k] = np.array(newd)

    return d2


def year_code(dates):
    """
    Returns a series whose values are the year numbers
    """
    a = 0
    yc = np.zeros(len(dates), dtype=int)
    year = dates[0].year
    for i, d in enumerate(dates):
        if d.year >= 2010:
            yc[i] = 0
        elif d.year >2010:
            yc[i] = 1
        else:
            a += 1
            year = d.year
        yc[i] = a

    return yc


# TODO: Fix this based on the actual dates
def beta_step(timeline):
    """
    Returns the date ranges where beta should be high, i.e., during summers.
    """
    bstep = np.zeros(len(timeline) + 200)
    for i, d in enumerate(timeline):
        if d.year == 2009:  # In 2009 h1n1 started in the summer
            if d.month != 6:
                bstep[i] = 1
        else:
            if d.month < 7 or d.month > 8:
                bstep[i] = 1
    return np.array(bstep)


def s_index(bstep):
    si = np.zeros(len(bstep))
    s = 0
    for i, b in enumerate(bstep):
        if i == 0: continue
        if b and not bstep[i - 1]:
            s += b
            si[i] = s
    return si

# # running the analysys
if __name__ == "__main__":
    dt = prepdata('data_Rt_dengue.csv', 0, 243, 1)
    modname = "Dengue_S0"
    # print dt['I'][:,1]

    bstep = beta_step(dt['time'])
    ycode = year_code(dt['time'])
    sindex = s_index(bstep)
    tcrit = [i for i in xrange(len(dt['time'])) if i]
    # Defining start and end of the simulations
    t0s = [0, dt['time'].index(datetime.datetime(2011, 8, 7)), dt['time'].index(datetime.datetime(2012, 10, 28))]
    tfs = t0s[1:] + [len(dt['time'])]
    print tfs
    # Interpolated Rt
    iRt = interp1d(np.arange(dt['Rt'].size), dt['Rt'], kind='cubic', bounds_error=False, fill_value=0)

    #~ P.plot(dt['Rt'],'*')
    #~ P.plot(np.arange(0,242,.2),[iRt(t) for t in np.arange(0,242,.2)])
    
    #~ tfs = np.array(tfs)-1
    #~ print len(sindex), len(dt['I']), len(dt['time'])
    #~ print t0s,tfs
    #~ P.figure()
    #~ P.gca().xaxis_date()
    #~ P.plot(dt['time'],dt['I'],'*-')
    #~ P.plot(dt['time'],dt['I'][:,1],'*-')
    #~ P.plot(dt['time'],bstep[:len(dt['time'])],'-*', label='bstep')
    #~ P.plot(dt['time'],ycode[:len(dt['time'])],'-+', label='ycode')
    #~ P.plot(dt['time'],.001*sindex[:len(dt['time'])],'-v', label='sindex')
    #~ P.legend()
    #~ P.gcf().autofmt_xdate()
    #~ P.show()

    tnames = ['s0', 's1', 's2']
    nt = len(tnames)
    pnames = ['S', 'I', 'R']
    nph = len(pnames)
    wl = dt['I'].shape[0]  #window length
    nw = len(dt['time']) / wl  #number of windows
    tf = wl * nw  #total duration of simulation
    inits[1] = dt['I'][0]
    print inits
    #~ print calc_R0s()
    y = model([s0,s1,s2])
    print s0, s1, s2
    #~ P.figure()
    P.plot(dt['I'],'*')
    P.plot(y[:,1])
    P.legend([pnames[1]])
    P.show()
    #Priors and limits for all countries
    tpars = [(9, 1), (9, 1), (9, 1)]
    tlims = [(0, 1), (0, 1), (0, 1)]

    F = FitModel(2000, model, inits, tf, tnames, pnames,
                 wl, nw, verbose=1, burnin=200, constraints=[])
    F.set_priors(tdists=nt * [st.beta],
                 tpars=tpars,
                 tlims=tlims,
                 pdists=[st.beta] * nph, ppars=[(1, 1)] * nph, plims=[(0, 1)] * nph)

    F.run(dt, 'DREAM', likvar=1e-7, pool=False, ew=0, adjinits=False, dbname=modname, monitor=['I'])
    #~ print F.AIC, F.BIC, F.DIC
    #print F.optimize(data=dt,p0=[s0,s1,s2], optimizer='scipy',tol=1e-55, verbose=1, plot=1)
    F.plot_results(['I'], dbname=modname, savefigs=1)
