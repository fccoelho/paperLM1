""" 
Script to generate the plots for the
multi-year simulations
"""
import pickle as CP
from itertools import cycle
from scipy import stats
import glob
from numpy import *
from matplotlib import font_manager
from matplotlib.ticker import FormatStrFormatter
import matplotlib.dates as mdates
import matplotlib.pyplot as P
from scipy.stats import gaussian_kde
import datetime
import matplotlib.collections as collections
import seaborn as sns
import pandas as pd
from sqlalchemy import create_engine
from collections import OrderedDict

from BIP.Bayes.PlotMeld import violin_plot, pred_new_cases, peakdet, plot_par_violin



def read_data(dbname):
    """
    Returns dataframes from tables of the sqlite database
    :param dbname: sqlite3 database name
    :return: list of pandas dataframes
    """
    eng = create_engine('sqlite:///{}'.format(dbname))
    df_data = pd.read_sql_table('data', eng, index_col=['time'], parse_dates={'time': '%Y-%m-%d %H:%M:%S'})
    df_pt = pd.read_sql_table('post_theta', eng, index_col=['time'], parse_dates={'time': '%Y-%m-%d %H:%M:%S'})
    df_series = pd.read_sql_table('series', eng, index_col=['time'], parse_dates={'time': '%Y-%m-%d %H:%M:%S'})
    return df_data, df_pt, df_series


def get_ordered_series(dbs):
    """
    Return three ordered dicts with the data ordered by year
    :param dbs: list of filenames of the sqlite databases
    :return:
    """
    p = {}
    s = {}
    o = {}
    pts = OrderedDict()
    series = OrderedDict()
    obs = OrderedDict()
    for db in dbs:
        y = db.split('_')[0][-4:]
        data, theta, srs = read_data(db)
        # print data
        s[y] = srs
        o[y] = data.I
        p[y] = theta
    for y in sorted(s.keys()):
        series[y] = s[y]
        obs[y] = o[y]
        pts[y] = p[y]
    return obs, series, pts


def create_tex_table(dbs):
    """
    Create Latex table with the Attack ratios for each epidemic
    :return:
    """
    obs, series, pts = get_ordered_series(dbs)

    head = r"""\begin{center}
\begin{tabular}{l|c|c|c}
\hline
"""
    head += r"""Year & Cases & median Attack Ratio $ $S_0$ \\
\hline
"""
    bot = r"""
\hline
\end{tabular}
\end{center}
        """
    body = r""
    st = []
    # years = sorted(list(series.keys()))
    print (series.keys())
    for i, (Y, V) in enumerate(series.items()):
        cases = obs[Y].sum()
        first_week = V.index[0]
        s0 = array(series[Y].S.ix[first_week])
        try:
            ratio = 1.0*cases/s0
            body += Y + r" & {:.3} & {:.2} ({:.2}-{:.2}) & {:.3}({:.2}-{:.2})\\".format(cases*100, nanmedian(ratio),
                                                                                      stats.scoreatpercentile(ratio, 2.5),
                                                                                      stats.scoreatpercentile(ratio, 97.5),
                                                                                      nanmedian(s0)*100,
                                                                                      stats.scoreatpercentile(s0, 2.5)*100,
                                                                                      stats.scoreatpercentile(s0, 97.2)*100
                                                                                      )
            body += "\n"
        except KeyError as e:
            print (Y, first_week, e)
        except ValueError as e:
            print (s0, e)

    return head + body + bot


def plot_AR_S0(dbs):
    """
    Scatter plot of Attack ratio by S0 by year.
    :param ar:
    :param s0:
    :return:
    """
    obs, series, pts = get_ordered_series(dbs)

    f, axes = P.subplots(7, 2, sharex='col', sharey='row')
    for i, (Y, V) in enumerate(series.items()):
        ax = axes.ravel()[i]
        cases = obs[Y].sum()
        first_week = V.index[0]
        s0 = array(series[Y].S.ix[first_week])
        ratio = 1.0*cases/s0
        ax.set_title("${}$".format(Y))
        ax.scatter(s0, ratio)
    axes.ravel()[-2].set_xlabel("$S_0$")
    axes.ravel()[-1].set_xlabel("$S_0$")



def series(nam='Dengue_S0'):
    fig = P.figure()
    pt, series, predseries, obs, weeks = _read_results(nam)
    wl = len(array(obs['time'])) / weeks
    # ~ print pt
    ax = fig.add_subplot(1, 1, 1)
    ax.grid()
    ax.set_ylabel('Weekly Incidence')
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b-%d,`%y"))
    t = array(obs['time'])
    plot_series2(t, obs, series, names=['I'], title=nam)
    P.setp(ax.get_xticklabels(), rotation=30, fontsize=8)
    P.setp(ax.get_yticklabels(), fontsize=8)

    fig2 = P.figure()
    pnames = pt[0].dtype.names

    labs = [r'$S_{0,04}$', r'$S_{0,05}$', r'$S_{0,06}$', r'$S_{0,07}$', r'$S_{0,08}$', r'$S_{0,09}$', r'$S_{0,10}$',
            r'$r_e$', r'$m$']
    #print len (labs), len(pnames)
    b, g, r, p = sns.color_palette("muted", 4)
    for i, n in zip(range(1, len(pnames)+1), pnames):

        ax2 = fig2.add_subplot(4, 4, i)
        #~ ax2.grid()
        sns.distplot(pt[0][n], color=p, ax=ax2[1, 1])

        P.xlabel(labs[i - 1])
        #P.ylabel(str(var(pt[0][n])))
        # P.setp(ax2.get_xticklabels(), fontsize=8)
        # P.setp(ax2.get_yticklabels(), fontsize=8)
        # ax2.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        #print P.setp(ax2)

def plot_rt_beta():
    df = pd.read_csv('data_Rt_dengue_complete.csv', header=0, delimiter=',', skiprows=[1, 2, 3], parse_dates=True)
    df.Rt.plot()
    df.Rt2.plot()
    df.lwr.plot()
    df.upr.plot()

def plot_concat_series(dbs):
    """
    Plot concatenated time series for susceptibles and infectious
    :param dbs:
    """
    series = []
    obs = []
    for db in dbs:
        data, theta, srs = read_data(db)
        series.append(srs)
        obs.append(data)
    upr = lambda x: stats.scoreatpercentile(x, 97.5)
    lwr = lambda x: stats.scoreatpercentile(x, 2.5)
    c = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    ax1 = P.subplot(211)

    for srs in series:
        co = c.next()
        s_median = srs.S.groupby(level='time').median()
        s_median.plot(style='k-', label='Median S')
        s_upr = srs.S.groupby(level='time').aggregate(upr)
        s_lwr = srs.S.groupby(level='time').aggregate(lwr)
        # P.legend(loc=0)
        P.fill_between(s_median.index, s_lwr, s_upr, facecolor=co, alpha=.2)
    ax1.set_ylabel("Fraction of population")
    ax2 = P.subplot(212, sharex=ax1)

    for srs, da in zip(series, obs):
        co = c.next()
        i_median = srs.I.groupby(level='time').median()
        i_median.plot(style='k-', label="Median I")
        # da.I.plot(style='ro', alpha=.5, label='obs')
        i_upr = srs.I.groupby(level='time').aggregate(upr)
        i_lwr = srs.I.groupby(level='time').aggregate(lwr)
        P.fill_between(i_median.index, i_lwr, i_upr, facecolor=co, alpha=.2)
        da.I.plot(style='b.', alpha=0.3)
        # P.plot_date(pd.to_datetime(da.index), da.I, 'b.', alpha=0.8, label='Observations')
    # P.legend(['Median', 'obs'])
    P.tight_layout()
    ax2.set_ylabel("Fraction of population")
    ax2.set_xlabel("Time (weeks)")
    P.savefig('../plots/concat_SI.svg')
    P.savefig('../plots/concat_SI.png', dpi=400)


def plot_pt_corr(df):
    """
    plot the correlation matrix of the posteriors of the parameters
    """

    f, ax = P.subplots(figsize=(9, 9))
    cmap = sns.blend_palette(["#00008B", "#6A5ACD", "#F0F8FF",
                              "#FFE6F8", "#C71585", "#8B0000"], as_cmap=True)
    sns.corrplot(df, annot=True, sig_stars=True, method='spearman',
                 diag_names=True, cmap=cmap, ax=ax)
    f.tight_layout()

def calc_mse(dbs):
    obs, series, pts = get_ordered_series(dbs)
    mse = OrderedDict()
    mae = OrderedDict()
    mape = OrderedDict()
    for Y, srs in series.items():
        n = len(obs[Y])
        i_median = srs.I.groupby(level='time').median()
        mse[Y] = (1./n) * sum(((i_median-obs[Y])**2)/obs[Y]**2)
        mae[Y] = (1./n) * sum(abs(i_median-obs[Y]))
        mape[Y] = (1./n) * sum(abs((i_median-obs[Y])/obs[Y]))
    return mse, mae, mape

if __name__ == "__main__":
    sns.set(style="darkgrid", palette="Set2")
    font_manager.FontProperties().set_size('x-small')
    #~ series('Dengue_S0_big')

    dbs = glob.glob("../DengueS*.sqlite")
    #plot_concat_series(dbs)
    print(create_tex_table(dbs))
    # plot_AR_S0(dbs)
    print (calc_mse(dbs))


    P.show()

