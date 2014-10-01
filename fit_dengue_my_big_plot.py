""" 
Script to generate the plots for the
multi-year simulations
"""
import cPickle as CP
from itertools import cycle
from scipy import stats
import glob
from numpy import *
from matplotlib import font_manager
from matplotlib.ticker import FormatStrFormatter
import matplotlib.dates as mdates
import matplotlib.pyplot as P
import Gnuplot
from scipy.stats import gaussian_kde
import datetime
import matplotlib.collections as collections
import seaborn as sns
import pandas as pd
from sqlalchemy import create_engine

from BIP.Bayes.PlotMeld import violin_plot, pred_new_cases, peakdet, plot_par_violin


def _read_results(nam):
    """
    read results from disk
    """
    pt, series, predseries = [], [], []
    nf = len(glob.glob('%s*.pickle' % nam))
    for w in range(nf):
        fn = "%s_%s.pickle" % (nam, w)
        print fn
        with open(fn, 'r') as f:
            a, b, obs, pred, samples = CP.load(f)
        # f.close()
        pt.append(a)
        series.append(b)
        predseries.append(pred)

    return pt, series, predseries, obs, nf

def read_data(dbname):
    """
    Returns dataframes from tables of the sqlite database
    :param dbname: sqlite3 database name
    :return: list of pandas dataframes
    """
    eng = create_engine('sqlite:///{}.sqlite'.format(dbname))
    df_data = pd.read_sql_table('data', eng, index_col=['pk', 'time'], parse_dates={'time': '%Y-%m-%d %H:%M:%S'})
    df_pt = pd.read_sql_table('post_theta', eng, index_col=['pk', 'time'], parse_dates={'time': '%Y-%m-%d %H:%M:%S'})
    df_series = pd.read_sql_table('series', eng, index_col=['pk', 'time'], parse_dates={'time': '%Y-%m-%d %H:%M:%S'})
    return df_data, df_pt, df_series


def create_tex_table(nms):
    pts = {}
    for nm in nms:
        pt, series, predseries, obs, weeks = _read_results(nm)
        pts[nm.split("-")[0]] = pt
    head = r"""\begin{center}
            \begin{tabular}{c|c c c}
            \hline
            """
    head += r"""Name & Belgium & Netherlands & Portugal\\
            \hline
            """
    bot = r"""
        \hline
        \end{tabular}
        \end{center}
        """
    body = r""
    st = []
    for n in pts['be'][0].dtype.names:
        body += n + r" & %1.3G, " % mean(pts['be'][0][n]) + "%1.3G " % var(pts['be'][0][n]) + \
                r" & %1.3G, " % mean(pts['nl'][0][n]) + "%1.3G " % var(pts['nl'][0][n]) + \
                r" & %1.3G, " % mean(pts['pt'][0][n]) + r"%1.3G \\" % var(pts['pt'][0][n])

    return head + body + bot


def plot_series2(tim, obs, series, names=[], title='Simulated vs Observed series', wl=7, lag=False):
    ser2 = {}
    for n in series[0].dtype.names:
        ser2[n] = concatenate([s[n] for s in series], axis=1)
    ls = ser2[n].shape[1]
    tim = tim[:ls]
    # print type (series)#.I.shape
    #fig =P.gcf()
    if not names:
        names = series[0].dtype.names
    c = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    s = cycle(['o', '^', '>', '<', 's', '*', '+', '1'])
    if isinstance(tim[0], datetime.date):
        lag = datetime.timedelta(int(lag) * wl)
    else:
        lag = int(lag) * wl
    for i, n in enumerate(names):
        ax = P.gca()
        co = c.next()
        labs = ['Influenzanet', 'EISN']
        if n in obs:
            if len(obs[n].shape) > 1:
                for i in range(obs[n].shape[1]):
                    ax.plot(tim, obs[n][:len(tim), i], c.next() + s.next(), alpha=0.7, label=labs[i])#print len(tim),  ls
            else:
                ax.plot(tim, obs[n][:len(tim)], c.next() + s.next(), alpha=0.7, label=n+"_obs")
        ax.plot(array(tim) + lag, median(ser2[n], axis=0), 'k-', label=n)
        lower = [stats.scoreatpercentile(t, 2.5) for t in ser2[n].T]
        upper = [stats.scoreatpercentile(t, 97.5) for t in ser2[n].T]
        if len(series) > 1:  #in the case of iterative simulations
            dif = (array(upper) - array(lower))
            dif = dif / max(dif) * 10
            pe, va = peakdet(dif, 1)
            xp = [0] + pe[:, 0].tolist() + [len(lower) - 1]
            lower = interp(range(len(lower)), xp, array(lower)[xp])  # valley-to-valley interpolated band
            upper = interp(range(len(upper)), xp, array(upper)[xp])  #peak-to-peak interpolated band
        ax.fill_between(array(tim) + lag, lower, upper, facecolor=co, alpha=0.2)
        #ax.fill_between(array(tim)+lag,lower,upper,facecolor='k',alpha=0.1)
        if i < (len(names) - 1): ax.xaxis.set_ticklabels([])
        ax.legend()
        if i == 0:
            ax.set_title(title)
    #ax.xaxis.set_visible(True)
    #P.title(title)
    P.xlabel('Weeks')
    #~ if isinstance(tim[0], datetime.date):
    #~ P.gcf().autofmt_xdate()


def multiplot():
    fig = P.figure()

    wl = 7
    for i, n in zip([1, 2, 3, 4], ["nl05-06", "nl06-07", "nl07-08", "nl08-09"]):
        pt, series, predseries, obs, weeks = _read_results(n)
        pt = [s.beta for s in pt]

        ax = fig.add_subplot(2, 2, i)
        ax.grid()
        if n in ["nl05-06", "nl07-08"]:
            ax.set_ylabel('Weekly Incidence')
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d, %Y'))

        pred_new_cases(obs, predseries, weeks, names=['I'], title=n.replace('nl', '20'))
        P.setp(ax.get_xticklabels(), rotation=30, fontsize=8)
        P.setp(ax.get_yticklabels(), fontsize=8)
        # plot_series2(array(obs['time']),obs,series, names=['I'],title=n)
        ax1 = ax.twinx()
        violin_plot(ax1, pt, array(obs['time'])[wl - 1::wl], bp=True)
        ax1.set_ylabel(r'$\beta$')
        #ax1.format_xdata = mdates.DateFormatter('%m-%d')t_manager.FontProperties().set_size('x-small')
        #~ print
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
        P.setp(ax1.get_xticklabels(), rotation=30, fontsize=8)
        P.setp(ax1.get_yticklabels(), fontsize=8)

    fig2 = P.figure()
    for i, n in zip([1, 2, 3, 4], ["nl05-06", "nl06-07", "nl07-08", "nl08-09"]):
        pt, series, predseries, obs, weeks = _read_results(n)
        ax2 = fig2.add_subplot(2, 2, i)
        ax2.grid()
        if n in ["nl05-06", "nl07-08"]:
            ax2.set_ylabel('Incidence')
        plot_series2(array(obs['time']), obs, series, names=['I'], title=n.replace('nl', '20'))
        ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
        P.setp(ax2.get_xticklabels(), rotation=30, fontsize=8)
        P.setp(ax2.get_yticklabels(), fontsize=8)


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
if __name__ == "__main__":
    sns.set(style="darkgrid", palette="Set2")
    font_manager.FontProperties().set_size('x-small')
    #~ series('Dengue_S0_big')
    
    data, theta, srs = read_data('Dengue_S0_big')
    plot_pt_corr(theta)
    # sns.tsplot(series.S, time=series.time)

    # series.groupby(level='time').plot()
    P.show()
    # print create_tex_table(['be-multiyear', 'nl-multiyear', 'pt-multiyear'])
