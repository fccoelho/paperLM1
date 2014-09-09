#
# Multi-year fit of influenza epidemics
#
import pyximport;
pyximport.install(pyimport=True)
from BIP.Bayes.Melding import FitModel
from scipy.integrate import odeint
import scipy.stats as st
import numpy as np
import pylab as P
import copy
from collections import defaultdict
import datetime
import pandas as pd



beta = 1.0  #Transmission coefficient
b1 = 19.9
b0 = 1.5  #low beta during winters
eta = .0  #infectivity of asymptomatic infections relative to clinical ones. FIXED
epsilon = 2.8  #latency rate
mu = 0  #nat/mortality rate


tau = 1.4  #infectious rate. FIXED
nflu = 0  #.00142
s0 = .241  #fraction of susceptibles at the beginning of first epidemic
s1 = 0.699;
s2 = 0.699;
s3 = 0.692;
s4 = 0.691;
s5 = 0.655;
s6 = .233

Ss = {0: s0, 1: s1, 2: s2, 3: s3, 4: s4, 5: s5, 6: s6}  #Multiple Ss map
r0 = 1.2
# Initial conditions
inits = [.9999, .0001, 0]  #initial values for state variables.


def model(theta):
    #setting parameters
    s0, s1, s2, s3, s4, s5, s6, r0 = theta
    Ss = {0: s0, 1: s1, 2: s2, 3: s3, 4: s4, 5: s5, 6: s6}  #Multiple Ss map

    #~ b1=19.9
    def sir(y, t):
        '''ODE model'''
        S, I, R = y
        b1 = r0 * tau / (As[astep[int(t)]] * Ss[astep[int(t)]])
        beta = b1 if bstep[int(t)] else b0

        #~ if (calc_R0s()<1.1).any():
        #~ print "flat"
        lamb = beta * I * S
        return [-lamb,  #dS/dt
                lamb - tau * I,  #dI/dt
                tau * I,  #dR/dt
        ]

    Y = np.zeros((wl, 3))
    for i in range(len(Ss)):
        t0 = t0s[i]
        tf = tfs[i]
        #~ if i>0:
        #~ inits[1] = Y[t0-1,1]
        inits[0] = Ss[i];
        inits[-1] = 1 - sum(inits[:1])
        Y[t0:tf, :] = odeint(sir, inits, np.arange(t0, tf, 1))  #,tcrit=tcrit)
        Y[t0:tf, 1] = Y[t0:tf, 1] * As[
            astep[int(t0)]]  #Adjusting the output to just the symptomatic I to compare with data

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
    data = pd.read_csv(fname)
    if not eday:
        eday = data.cases.shape[0]
    dates = data.date[sday:eday]
    sday = sday  #starting day for the fitting
    eday = rawprev.size() if eday == None else eday
    #print data.dtype.names
    incidence = data.ili[sday:eday]  #daily incidence
    dur = 1  #./tau
    den = np.mean(data.active[sday:eday])
    rawprev = np.convolve(incidence, np.ones(dur), 'same') / den
    if weekly:
        dates = dates[6:eday:7]  #First date is the last day of first week
        prev = np.zeros(len(dates))
        for w in range(len(dates)):
            prev[w] = sum(rawprev[w * 7:w * 7 + 7])
        rawprev = prev
        print len(rawprev), len(dates)
    #~ P.plot(dates,incidence)
    #~ P.plot(dates,rawprev)
    #P.plot(dates,data.active[sday:eday])
    #~ P.show()
    #Doing moving average of order mao
    if mao > 1:
        sw = np.ones(mao, dtype=float) / mao  #smoothing window
        prev = np.convolve(rawprev, sw, 'same')  #Smoothing data (ma)
    else:
        prev = rawprev
    #prev = data.incidence[sday:eday]
    assert len(prev) == len(dates)
    d = {'time': dates, 'I': np.nan_to_num(prev)}
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
                ew = (newd[-1][0] % 52) + 1  #expected week number
                ed = newd[-1][1] + datetime.timedelta(7)  #expected date
                while w != ew:
                    newd.append([ew, ed, 1e-1])
                    ed = newd[-1][1] + datetime.timedelta(7)
                    ew = (newd[-1][0] % 52) + 1
                    p = np.nan
            newd.append([w, d, p])
        d2[k] = np.array(newd)

    return d2


def beta_step(timeline):
    """
    Returns the date ranges where beta should be low, i.e., during summers. 
    """
    bstep = np.zeros(len(timeline) + 200)
    for i, d in enumerate(timeline):
        if d.year == 2009:  #In 209 h1n1 started in the summer
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


def calc_R0s():
    r0s = []
    for i in range(7):
        r0s.append(b1 / tau * Ss[i])
    return np.array(r0s)


def constrain_Re(theta):
    res = []
    for i in range(7):
        res.append(theta[i] * theta[i + 7] * theta[-2] / tau > 1.1)
    return np.array(res).all()

## running the analysys

#weekly prevalence


dt = prepdata('data_Rt_dengue.csv', 0, 243, 1);
modname = "Dengue_S0"
#print dt['I'][:,1]

bstep = beta_step(dt['time'])
astep = alpha_step(bstep)
sindex = s_index(bstep)
tcrit = [i for i in xrange(len(dt['time'])) if i]
t0s = [0] + [i for i in range(len(dt['time'])) if sindex[i]]
tfs = t0s[1:] + [len(dt['time'])]
#~ tfs = np.array(tfs)-1
#~ print len(sindex), len(dt['I']), len(dt['time'])
#~ print t0s,tfs
#~ P.figure()
#~ P.gca().xaxis_date()
#~ P.plot(dt['time'],dt['I'],'*-')
#~ P.plot(dt['time'],dt['I'][:,1],'*-')
#~ P.plot(dt['time'],bstep[:len(dt['time'])],'-*', label='bstep')
#~ P.plot(dt['time'],astep[:len(dt['time'])],'-+', label='astep')
#~ P.plot(dt['time'],.001*sindex[:len(dt['time'])],'-v', label='sindex')
#~ P.legend()
#~ P.gcf().autofmt_xdate()
#~ P.show()

tnames = ['a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6',
          's0', 's1', 's2', 's3', 's4', 's5', 's6', 'r0', 'm']
nt = len(tnames)
pnames = ['S', 'I', 'R']
nph = len(pnames)
wl = dt['I'].shape[0]  #window length 446
nw = len(dt['time']) / wl  #number of windows
tf = wl * nw  #total duration of simulation
inits[1] = max(st.nanmean(dt['I'][0]), 0.0001)
#~ #print inits
#~ print calc_R0s()
#~ y = model2([a0,a1,a2,a3,a4,a5,a6,s0,s1,s2,s3,s4,s5,s6,r0,k])
#~ P.figure()
#~ P.plot(dt['I'],'*')
#~ P.plot(y[:,1])
#~ P.legend(['d','e']+[pnames[1]])
#~ P.show()
#Priors and limits for all countries
tpars_be = [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1)] + \
           [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1)] + [(1, .4), (0, 4e-6)]
tlims_be = [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1)] + \
           [(0.1, 1), (0.2, 1), (0.2, 1), (0.2, 1), (0.2, 1), (0.2, 1), (0.1, 1)] + [(1.01, 1.4),
                                                                                     (0, .00004)]  #beta and k
tpars_nl = [(0, .6), (0, .4), (0, .4), (0, .4), (0, .4), (0, .4), (0, .6)] + \
           [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1)] + [(1, .4), (0, 4e-6)]  #beta and k
tlims_nl = [(0, .35), (0, .4), (0, .4), (0, .4), (0, .4), (0, .3), (0, .5)] + \
           [(0.2, .8), (0.2, 1), (0.2, .9), (0.2, .8), (0.2, .9), (0.2, .8), (0.2, .8)] + [(1.01, 1.4),
                                                                                           (0, .00004)]  #beta and k
tpars_pt = [(0, .6), (0, .4), (0, .4), (0, .4), (0, .4), (0, .4), (0, .6)] + \
           [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1)] + [(1, .4), (0, 4e-6)]  #beta and k
tlims_pt = [(0, .5), (0, .4), (0, .30), (0, .4), (0, .30), (0, .3), (0, .5)] + \
           [(0.2, .8), (0.2, 1), (0.1, .9), (0.2, .8), (0.1, .9), (0.2, .8), (0.2, .8)] + [(1.1, 1.4),
                                                                                           (0, .00004)]  #beta and k

F = FitModel(2000, model, inits, tf, tnames, pnames,
             wl, nw, verbose=1, burnin=2000, constraints=[])
F.set_priors(tdists=nt * [st.uniform],
             tpars=tpars_be,
             tlims=tlims_be,
             pdists=[st.uniform] * nph, ppars=[(0, 1)] * nph, plims=[(0, 1)] * nph)

F.run(dt, 'DREAM', likvar=1e-6, pool=False, ew=0, adjinits=False, dbname=modname, monitor=['I'])
print F.AIC, F.BIC, F.DIC
#~ print F.optimize(data=dt,p0=[a0,a1,a2,a3,a4,a5,a6,s0,s1,s2,s3,s4,s5,s6,b1,k], optimizer='oo',tol=1e-55, verbose=1, plot=1)
F.plot_results(['I'], dbname=modname, savefigs=1)
