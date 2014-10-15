# coding:utf8

from BIP.SDE.gillespie import Model
import time
from numpy import array, genfromtxt
import pylab as P

vnames = ['S', 'I1', 'I2', 'I3', 'I4', 'R1', 'R2', 'R3', 'R4', 'I11',
          'I12', 'I13', 'I14', 'I21', 'I22', 'I23', 'I24', 'I31', 'I32',
          'I33', 'I34', 'I41', 'I42', 'I43', 'I44', 'R']

# transitions
# S->I1	S->I2	S->I3	S->I4	I1=>R1	I2=>R2	I3=>R3	I4=>R4	R1=>I11	R1=>I12	R1=>I13	R1=>I14	R2=>I21	R2=>I22	R2=>I23	R2=>I24	R3=>I31	R3=>I32	R3=>I33	R3=>I34	R4=>I41	R4=>I42	R4=>I43	R4=>I44	I11=>R	I12=R	I13=>R	I14=R	I21=>R	I22=R	I23=>R	I24=R	I31=>R	I32=R	I33=>R	I34=R	I41=>R	I42=R	I43=>R	I44=R

tm = genfromtxt('tmat.csv', delimiter=',', skip_header=2, usecols=range(1, 67), dtype=int)

# Initial Conditions
# Initial conditions for each of the four Serotype model

N = 1000

ini = [N - 40, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

pars = (1 / 70.,  # mu
        400,  # beta
        1,  # phi
        100,  # sigma
        1)  # gamma

# Propensity functions
def gen_prop_functs():
    nat_code = "def fnat(r, ini): return r[0]*N\n"
    mu_code = "def fmu{}(r, ini): return r[0]*ini[{}]\n"
    inf_code = """def fS(r, ini): return ini[0]*sum([lamb(i) for i in range(4)])\n"""
    inf2_code = """def f{}(r, ini): return ini[0]* lamb({})\n"""
    rec_code = """def f{}(r, ini): return r[-2]*ini[{}]\n"""
    inf3_code = """def f{}(r, ini): return r[-1]*ini[{}]*lamb({})\n"""

    nt = tm.shape[1]
    with open('propfun.py', 'w') as f:
        f.write("pars = {}\nini = {}\nN = {}\n\n".format(pars, ini, N))
        f.write("lamb = lambda i: {}*(ini[i+1] + sum(ini[9+i : 9+12+i+1 : 4])-ini[9+4*i+1])\n".format(pars[1]))
        f.write(inf_code)

        f.write("\n#Infection rates\n\n")
        for i, n in enumerate(['S_I1', 'S_I2', 'S_I3', 'S_I4']):  # Infection rates
            f.write(inf2_code.format(n, i))

        f.write("\n#Primary recovery rates\n\n")
        for i, n in enumerate(['I1_R1', 'I2_R2', 'I3_R3', 'I4_R4']):
            f.write(rec_code.format(n, i + 5))

        f.write("\n#Secondary Infection rates\n\n")
        for i, n in enumerate(
                ['R1_I11', 'R1_I12', 'R1_I13', 'R1_I14', 'R2_I21', 'R2_I22', 'R2_I23', 'R2_I24', 'R3_I31',
                 'R3_I32', 'R3_I33', 'R3_I34', 'R4_I41', 'R4_I42', 'R4_I43', 'R4_I44']):
            f.write(inf3_code.format(n, i // 4 + 5, i % 4))

        f.write("\n#Secondary Recovery rates\n\n")
        for i, n in enumerate(
                ['I11_R', 'I12_R', 'I13_R', 'I14_R', 'I21_R', 'I22_R', 'I23_R', 'I24_R', 'I31_R',
                 'I32_R', 'I33_R', 'I34_R', 'I41_R', 'I42_R', 'I43_R', 'I44_R']):
            f.write(rec_code.format(n, i + 9))

        f.write("\n#Mortality rates\n\n")
        for i, n in enumerate(
                ['S', 'I1', 'I2', 'I3', 'I4', 'R1', 'R2', 'R3', 'R4', 'I11', 'I12', 'I13', 'I14', 'I21', 'I22', 'I23',
                 'I24', 'I31', 'I32', 'I33', 'I34', 'I41', 'I42', 'I43', 'I44', 'R']):  # Mortality rates
            f.write(mu_code.format(n, i))


gen_prop_functs()
from propfun import *

propensity = [fS, fS_I1, fS_I2, fS_I3, fS_I4, fI1_R1, fI2_R2, fI3_R3, fI4_R4, fR1_I11, fR1_I12, fR1_I13, fR1_I14,
              fR2_I21, fR2_I22, fR2_I23, fR2_I24, fR3_I31, fR3_I32, fR3_I33, fR3_I34, fR4_I41, fR4_I42, fR4_I43,
              fR4_I44, fI11_R, fI12_R, fI13_R, fI14_R, fI21_R, fI22_R, fI23_R, fI24_R, fI31_R, fI32_R, fI33_R, fI34_R,
              fI41_R, fI42_R, fI43_R, fI44_R, fmuS, fmuI1, fmuI2, fmuI3, fmuI4, fmuI11, fmuI12, fmuI13, fmuI14, fmuI21,
              fmuI22, fmuI23, fmuI24, fmuI31, fmuI32, fmuI33, fmuI34, fmuI41, fmuI42, fmuI43, fmuI44, fmuR]

M = Model(vnames=vnames, rates=pars, inits=ini, tmat=tm, propensity=propensity)
t0 = time.time()
M.run(tmax=100, reps=1, viz=0, serial=0)
print 'total time: {} seconds'.format(time.time() - t0)
t, series, steps, evts = M.getStats()
ser = series.mean(axis=0)
# print evts
P.plot(t, ser)
P.legend(M.vn, loc=0)
# P.plot(t, ser[:, 0::3], 'r-.')  # S plots
# P.plot(t, ser[:, 1::3], 'g-^')  # I plots
# P.plot(t, ser[:, 2::3], 'b-o')  # R plots
# P.legend(M.vn[0::3] + M.vn[1::3] + M.vn[2::3], loc=0)
P.show()

