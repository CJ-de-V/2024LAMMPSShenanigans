import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from distinctipy import distinctipy
from seaborn import scatterplot as scat

colors = distinctipy.get_colors(25)  # color coding for plots
number = 1


# Seaborn scatter plot with 2 x legends supported
def myplot(datfram, Xdata, Ydata, leghue=None, legstyle=None, filter=[], log=True):
    global number
    plt.figure(number)
    number += 1

    for fil in filter:
        datfram = datfram[datfram[fil[0]].isin(fil[1])]
        # datfram = datfram[datfram[filters[0]] == filters[1]]

    titlstr = "Plot of " + Xdata + " vs " + Ydata + " with restrictions: " + str(filter)
    f, ax = plt.subplots()
    if log:
        titlstr = 'log-log ' + titlstr
        plt.xlabel("Log(" + Xdata + ")")
        plt.ylabel("Log(" + Ydata + ")")
        ax.set_yscale("log", base=2)
        ax.set_xscale("log", base=2)
    else:
        plt.xlabel(Xdata)
        plt.ylabel(Ydata)
        ax.set_yscale("linear")
        ax.set_xscale("linear")
    plt.title(titlstr)
    scat(data=datfram, x=Xdata, y=Ydata, hue=leghue, style=legstyle, palette="deep")


# Regular plot with marker connections available
def connectedplot(datfram, Xdata, Ydata, legend, filter=[], xerror=0, yerror=0, marker='o', log=True):
    global number
    plt.figure(number)
    number += 1

    for filters in filter:
        datfram = datfram[datfram[filters[0]] == filters[1]]

    print(datfram)
    legset = datfram[legend].unique()  # list of all entries in legend
    for i in range(len(legset)):
        ndata = datfram.loc[df[legend] == legset[i]]
        x = np.array(ndata[Xdata].tolist())
        y = np.array(ndata[Ydata].tolist())
        if log:
            plt.loglog(x, y, marker, base=2, c=colors[i])  # . formerly o-
        else:
            plt.plot(x, y, marker, c=colors[i])

        if xerror != 0:
            plt.errorbar(x, y, xerr=np.array(ndata[xerror].tolist()), c=colors[i])
        if yerror != 0:
            plt.errorbar(x, y, yerr=np.array(ndata[yerror].tolist()), c=colors[i])
    titlstr = "Plot of " + Xdata + " vs " + Ydata + " for various " + legend + " with restrictions: " + str(filter)
    if log:
        titlstr = 'log-log ' + titlstr
    plt.title(titlstr)
    plt.xlabel("Log(" + Xdata + ")")
    plt.ylabel("Log(" + Ydata + ")")
    plt.legend(legset, title=legend)


def save_multi_image(filename):
    pp = PdfPages(filename)
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()


df = pd.read_csv('data.csv')
df['Rg/zd'] = df['Rg'] / df['zd']
# Np	Nt	Kp	Kt A	lp	Rg	zd	lb lt	sd_lp	sd_Rg	sd_zd  Rg/zd E2E dot
# data columns for faster copypasting

myplot(df, 'Rg', 'lp', 'Kp', 'A', filter=[['Nt', [8]], ['A', [30, 90]]])
# Notice the nice correlation, Rg (Np dependent) caps, lp (Kp dependent) doesn't

# linear correspondence between RG and zd found, holds for all Kp and A "roughly" should be possible to stick a line
#through it corresponding to the free area occupiable to the polymer at the given cone width
# has a geometrical explanation... well as zd gets bumped up your free space grows linearly
myplot(df, 'zd', 'Rg', 'A', 'Np', [['Kp', [64]], ['Nt', [8]]], log=False)

# as you confine it further Rg/zd shrinks. Reflects that Rg shrinks as you confine it
myplot(df, 'A', 'Rg/zd', 'Np', 'Nt', [['Kp',[ 64]], ['Nt', [8]]], log=False)

# This is to aid in showing the above: As you squish the cone Rg shrinks, but only above a certain length/Rg
# it has to be large enough before squishing so that it plays a roll, if its too small then squishing it will
# not have any visible impact
myplot(df, 'A', 'Rg', leghue='Np', filter=[['Nt', [8]], ['Kp', [64]], ['Np', [17, 65]]], log=False)

# Ties in to the above E2E is proportional to Rg proportional to Np, this just shows that squishing the cone does force
# the endpoints closer together, but only after a certain size prerequisite has been reached
myplot(df, 'A', 'E2E', 'Np', filter=[['Nt', [8]], ['Kp', [64]]], log=False)

# not dissimilar to the above: it is shown that angle causes bending for certain lengths beyond certain stiffnesses
myplot(df, 'A', 'Angle', 'Np', 'Kp', filter=[['Kp', [0.25, 64]], ['Nt', [8]], ['Np', [9, 129]]], log=False)

myplot(df, 'A', 'lt', 'Kp', 'Np', filter=[ ['Kp', [64]], ['Nt', [8]]], log=False)
#Filter ['A', [30, 40]], ['Kp', [64, 32]], ['Nt', [8, 16]]
#I do not believe there to be any significant alterations to the length of the tether bonds. I propose most of the change
#comes from straightening out the tether most likely, keeping the bonds fixed

myplot(df, 'A', 'Rg', 'Kp', 'Np', [['Kp', [0.25,64]], ['Nt', [8]]], log=False)

myplot(df, 'Nt', 'lt', 'A', 'Np', [['Kp', [64]], ['Np', [129]]], log=False)

# We need to see something sensible for this, what was not sensible was the fact that we were unable to find the
# change that Kp brings on to this.

save_multi_image("Results.pdf")
plt.show()
