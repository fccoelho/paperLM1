# coding:utf8
from __future__ import division
from cgillespie import Model
import time
from numpy import array, genfromtxt
import pylab as P
import numpy as np

from itertools import cycle

vnames = ['S', 'I1', 'I2', 'I3', 'I4', 'R1', 'R2', 'R3', 'R4',
          'I12', 'I13', 'I14', 'I21', 'I23', 'I24', 'I31', 'I32',
            'I34', 'I41', 'I42', 'I43', 'R']

# transitions
# S->I1	S->I2	S->I3	S->I4	I1=>R1	I2=>R2	I3=>R3	I4=>R4	R1=>I12	R1=>I13	R1=>I14	R2=>I21	R2=>I23	R2=>I24	R3=>I31	R3=>I32	R3=>I34	R4=>I41	R4=>I42	R4=>I43	I12=R	I13=>R	I14=R	I21=>R	I23=>R	I24=R	I31=>R	I32=R	I34=R	I41=>R	I42=R	I43=>R

tm = genfromtxt('tmat.csv', delimiter=',', skip_header=2, usecols=range(1, 56), dtype=int)

# Initial Conditions
# Initial conditions for each of the four Serotype model

N = 5000

ini = (N, 50, 30, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

pars = (1 / 70.,  # mu
        400,  # beta
        1.,  # phi
        .1,  # sigma
        1., # gamma
        0.001, #delta : cross immunity protection
        )
        
        
# infectives for each serotype
inf_types = {1: [ini[1], ini[12], ini[15], ini[18]],
             2: [ini[2], ini[9], ini[16], ini[19]],
             3: [ini[3], ini[10], ini[13], ini[20]],
             4: [ini[4], ini[11], ini[14], ini[17]]
            }

# Propensity functions
def gen_prop_functs():
    nat_code = "def fnat(r, ini): return r[0]*sum(ini)\n"
    mu_code = "def fmu{}(r, ini): return r[0]*max(0,ini[{}])\n"
    inf_code = """def f{}(r, ini): return ini[0]* lamb({})\n"""
    rec_code = """def f{}(r, ini): return r[3]*max(0,ini[{}])\n"""
    inf2_code = """def f{}(r, ini): return r[4]*r[5]*ini[{}]*lamb({})\n"""

    nt = tm.shape[1]
    with open('propfun.py', 'w') as f:
        f.write("import numpy as np\n\n")
        
        f.write("pars = {}\nini = {}\nN = {}\n\n".format(pars, ini, N))
        f.write("# infectives for each serotype\n\n")
        f.write("""inf_types = {1: [ini[1], ini[12], ini[15], ini[18]],
             2: [ini[2], ini[9], ini[16], ini[19]],
             3: [ini[3], ini[10], ini[13], ini[20]],
             4: [ini[4], ini[11], ini[14], ini[17]]
            }\n\n""")
        #TODO: Fix the indexing below
        f.write("lamb = lambda i: {}*(ini[i+1] + sum(np.array(inf_types[i+1])*{}))\n".format(pars[1], pars[2]))
        # f.write(inf_code)

        f.write("\n#Natality\n\n")
        f.write(nat_code)

        f.write("\n#Infection rates\n\n")
        for i, n in enumerate(['S_I1', 'S_I2', 'S_I3', 'S_I4']):  # Infection rates
            f.write(inf_code.format(n, i))

        f.write("\n#Primary recovery rates\n\n")
        for i, n in enumerate(['I1_R1', 'I2_R2', 'I3_R3', 'I4_R4']):
            f.write(rec_code.format(n, i+1))

        f.write("\n#Secondary Infection rates\n\n")
        for i, n in enumerate(
                ['R1_I12', 'R1_I13', 'R1_I14', 'R2_I21', 'R2_I23', 'R2_I24', 'R3_I31',
                 'R3_I32', 'R3_I34', 'R4_I41', 'R4_I42', 'R4_I43']):
            type_i = [1,2,3,0,2,3,0,1,3,0,1,2] #serotype number for each infection
            f.write(inf2_code.format(n, i // 3 + 5, type_i[i]))

        f.write("\n#Secondary Recovery rates\n\n")
        for i, n in enumerate(
                ['I12_R', 'I13_R', 'I14_R', 'I21_R', 'I23_R', 'I24_R', 'I31_R',
                 'I32_R', 'I34_R', 'I41_R', 'I42_R', 'I43_R']):
            f.write(rec_code.format(n, i + 9))

        f.write("\n#Mortality rates\n\n")
        for i, n in enumerate(
                ['S', 'I1', 'I2', 'I3', 'I4', 'R1', 'R2', 'R3', 'R4', 'I12', 'I13', 'I14', 'I21', 'I23',
                 'I24', 'I31', 'I32', 'I34', 'I41', 'I42', 'I43', 'R']):  # Mortality rates
            f.write(mu_code.format(n, i))


gen_prop_functs()
from propfun import *

propensity = [fnat,                             # Natality
              fS_I1, fS_I2, fS_I3, fS_I4,       # Primary infections
              fI1_R1, fI2_R2, fI3_R3, fI4_R4,   # primary recovery
              fR1_I12, fR1_I13, fR1_I14,        #Secondary infections
              fR2_I21, fR2_I23, fR2_I24,        #Secondary infections
              fR3_I31, fR3_I32, fR3_I34,        #Secondary infections
              fR4_I41, fR4_I42, fR4_I43,        #Secondary infections
              fI12_R, fI13_R, fI14_R,           # Secondary recovery
              fI21_R, fI23_R, fI24_R,           # Secondary recovery
              fI31_R, fI32_R, fI34_R,           # Secondary recovery
              fI41_R, fI42_R, fI43_R,           # Secondary recovery
              fmuS,
              fmuI1, fmuI2, fmuI3, fmuI4,
              fmuR1, fmuR2, fmuR3, fmuR4,
              fmuI12, fmuI13, fmuI14,
              fmuI21, fmuI23, fmuI24,
              fmuI31, fmuI32, fmuI34,
              fmuI41, fmuI42, fmuI43,
              fmuR]

assert len(propensity) == tm.shape[1]

M = Model(vnames=vnames, rates=pars, inits=ini, tmat=tm, propensity=propensity)
t0 = time.time()
M.run(tmax=500, reps=1)
print 'total time: {} seconds'.format(time.time() - t0)
t, series, steps = M.getStats()

ser = series.mean(axis=2)

#  Calculatin prevalence by serotype

p1 = sum(ser[:,i] for i, n in enumerate(vnames) if (n.startswith('I') and '1' in n))
p2 = sum(ser[:,i] for i, n in enumerate(vnames) if (n.startswith('I') and '2' in n))
p3 = sum(ser[:,i] for i, n in enumerate(vnames) if (n.startswith('I') and '3' in n))
p4 = sum(ser[:,i] for i, n in enumerate(vnames) if (n.startswith('I') and '4' in n))
# print evts
co = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
sy = cycle(['o', '^', '>', '<', 's', '*', '+', '1'])
for s in range(ser.shape[1]):
    P.plot(t, ser[:,s], co.next()+sy.next()+'-')
P.legend(vnames, loc=0)
P.figure()
P.plot(t, p1, 'r-.', label='DENV1')
P.plot(t, p2, 'g-*', label='DENV2')
P.plot(t, p3, 'b-.', label='DENV3')
P.plot(t, p4, 'y-.', label='DENV4')
P.legend()
P.figure()
P.plot(t,ser.sum(axis=1), label="Total")
P.legend()
# P.plot(t, ser[:, 1::3], 'g-^')  # I plots
# P.plot(t, ser[:, 2::3], 'b-o')  # R plots
# P.legend(M.vn[0::3] + M.vn[1::3] + M.vn[2::3], loc=0)
P.show()

