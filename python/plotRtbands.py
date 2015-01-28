# -*- coding:utf-8 -*-
u"""
Created on 22/01/15
by fccoelho
license: GPL V3 or Later
"""


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as P
import datetime
import numpy as np




DF = pd.read_csv("../DATA/data_Rt_dengue_complete.csv", header=0, delimiter=',', skiprows=[1, 2, 3], parse_dates=True)
# Creating a date-time index  to facilitate date axes below.
DF.index = pd.to_datetime(DF.start)
# DF.set_index(['start'])
print DF.info()

# create two subplots with  shared x axes
fig, (ax1, ax2) = P.subplots(2, 1, sharex=True)
ax1.plot(DF.index, DF.cases,'-')
ax1.set_ylabel('Reported cases')
# print DF.start
ax2.fill_between(DF.index, DF.lwr, DF.upr, facecolor='blue', alpha=.2)
ax2.plot(DF.index, DF.Rt2, 'k-', alpha=.8, lw=.5)
ax2.set_ylabel('$R_t$')
ax2.set_ylim(0, 5)
fig.tight_layout()
fig.autofmt_xdate()
P.savefig("../plots/rt_series.png", dpi=400)
P.show()
