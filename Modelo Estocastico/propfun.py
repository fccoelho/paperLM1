import numpy as np

pars = (0.02, 0.0004, 1.1, 0.1, 1.1)
ini = [5000, 150, 100, 75, 50, 15, 15, 15, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
N = 5000

lamb = lambda i: 0.0004*(ini[i+1] + sum(np.array(ini[9+i : 9+12+i+1 : 4])*1.1)-ini[9+4*i+1])

#Natality

def fnat(r, ini): return r[0]*sum(ini)

#Infection rates

def fS_I1(r, ini): return ini[0]* lamb(0)
def fS_I2(r, ini): return ini[0]* lamb(1)
def fS_I3(r, ini): return ini[0]* lamb(2)
def fS_I4(r, ini): return ini[0]* lamb(3)

#Primary recovery rates

def fI1_R1(r, ini): return r[3]*ini[1]
def fI2_R2(r, ini): return r[3]*ini[2]
def fI3_R3(r, ini): return r[3]*ini[3]
def fI4_R4(r, ini): return r[3]*ini[4]

#Secondary Infection rates

def fR1_I11(r, ini): return r[4]*ini[5]*lamb(0)
def fR1_I12(r, ini): return r[4]*ini[5]*lamb(1)
def fR1_I13(r, ini): return r[4]*ini[5]*lamb(2)
def fR1_I14(r, ini): return r[4]*ini[5]*lamb(3)
def fR2_I21(r, ini): return r[4]*ini[6]*lamb(0)
def fR2_I22(r, ini): return r[4]*ini[6]*lamb(1)
def fR2_I23(r, ini): return r[4]*ini[6]*lamb(2)
def fR2_I24(r, ini): return r[4]*ini[6]*lamb(3)
def fR3_I31(r, ini): return r[4]*ini[7]*lamb(0)
def fR3_I32(r, ini): return r[4]*ini[7]*lamb(1)
def fR3_I33(r, ini): return r[4]*ini[7]*lamb(2)
def fR3_I34(r, ini): return r[4]*ini[7]*lamb(3)
def fR4_I41(r, ini): return r[4]*ini[8]*lamb(0)
def fR4_I42(r, ini): return r[4]*ini[8]*lamb(1)
def fR4_I43(r, ini): return r[4]*ini[8]*lamb(2)
def fR4_I44(r, ini): return r[4]*ini[8]*lamb(3)

#Secondary Recovery rates

def fI11_R(r, ini): return r[3]*ini[9]
def fI12_R(r, ini): return r[3]*ini[10]
def fI13_R(r, ini): return r[3]*ini[11]
def fI14_R(r, ini): return r[3]*ini[12]
def fI21_R(r, ini): return r[3]*ini[13]
def fI22_R(r, ini): return r[3]*ini[14]
def fI23_R(r, ini): return r[3]*ini[15]
def fI24_R(r, ini): return r[3]*ini[16]
def fI31_R(r, ini): return r[3]*ini[17]
def fI32_R(r, ini): return r[3]*ini[18]
def fI33_R(r, ini): return r[3]*ini[19]
def fI34_R(r, ini): return r[3]*ini[20]
def fI41_R(r, ini): return r[3]*ini[21]
def fI42_R(r, ini): return r[3]*ini[22]
def fI43_R(r, ini): return r[3]*ini[23]
def fI44_R(r, ini): return r[3]*ini[24]

#Mortality rates

def fmuS(r, ini): return r[0]*ini[0]
def fmuI1(r, ini): return r[0]*ini[1]
def fmuI2(r, ini): return r[0]*ini[2]
def fmuI3(r, ini): return r[0]*ini[3]
def fmuI4(r, ini): return r[0]*ini[4]
def fmuR1(r, ini): return r[0]*ini[5]
def fmuR2(r, ini): return r[0]*ini[6]
def fmuR3(r, ini): return r[0]*ini[7]
def fmuR4(r, ini): return r[0]*ini[8]
def fmuI11(r, ini): return r[0]*ini[9]
def fmuI12(r, ini): return r[0]*ini[10]
def fmuI13(r, ini): return r[0]*ini[11]
def fmuI14(r, ini): return r[0]*ini[12]
def fmuI21(r, ini): return r[0]*ini[13]
def fmuI22(r, ini): return r[0]*ini[14]
def fmuI23(r, ini): return r[0]*ini[15]
def fmuI24(r, ini): return r[0]*ini[16]
def fmuI31(r, ini): return r[0]*ini[17]
def fmuI32(r, ini): return r[0]*ini[18]
def fmuI33(r, ini): return r[0]*ini[19]
def fmuI34(r, ini): return r[0]*ini[20]
def fmuI41(r, ini): return r[0]*ini[21]
def fmuI42(r, ini): return r[0]*ini[22]
def fmuI43(r, ini): return r[0]*ini[23]
def fmuI44(r, ini): return r[0]*ini[24]
def fmuR(r, ini): return r[0]*ini[25]
