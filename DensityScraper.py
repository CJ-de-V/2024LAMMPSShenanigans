import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt


def skiplines(f, n):
    for i in range(n):
        f.readline()


basedir = os.getcwd()
p1 = np.linspace(0., 0., 3)
p2 = np.linspace(0., 0., 3)
N = 0  # Number of particles
Nb = 0  # number of bonds
Nc = 0  # number of correlation functions
bondvectors = 0
avgbondlength = 0

Np = np.array([ 64]) + 1  # +1 for even sides
Nt = [64]  #
Kp = [64]  #
Kt = [0.25]
alpha = [30,40,50,60,70,80]  #

for nop in Np:
    for nt in Nt:
        if nt > nop:  # No Tether dominated dynamics
            break
        discriminator = "Np_" + str(nop) + "_Nt_" + str(nt)
        workingdir = basedir + "/output/" + discriminator
        os.chdir(workingdir)
        n2 = nop // 2  # number of bondvectors inthehalf
        for kp in Kp:
            for kt in Kt:
                for alp in alpha:  # alp in alpha:
                    discriminator = "Np_" + str(nop) + "_Nt_" + str(nt) + "_Kp_" + str(kp) + "_Kt_" + str(
                        kt) + "_Alpha_" + str(alp)
                    Nb = nop - 1
                    Nc = Nb - 1
                    bondvectors = np.zeros((Nb, 3))
                    correlationfunctions = np.zeros(Nc)
                    avgbondlength = 0
                    zvals =[]

                    # persistence and bond length calculation
                    with open('dump' + discriminator + '.dynamics') as datafile:
                        print('starting with analysis of: ' + discriminator)
                        # Index-1 of the above indicates how many atoms separate the bondvectors
                        skiplines(datafile, 9)  # skips first boilerplate
                        line = 'liney'

                        while line != '':
                            line = (datafile.readline()).split(" ")
                            p1 = np.array([float(line[2]), float(line[3]), float(line[4])])
                            zvals.append(p1[2])
                            # position 1 1 is now real
                            p2 = np.array([0.0, 0.0, 0.0])
                            # reads particle 2 to N in and finds the bondvectors
                            startpos = p1  # position of leftmost bond
                            for i in range(1, nop):  # reading past entry 1 which we already read
                                line = (datafile.readline()).split(" ")
                                p2 = np.array([float(line[2]), float(line[3]), float(line[4])])
                                zvals.append(p2[2])
                                p1 = p2

                            ### Tether Shenaigans
                            skiplines(datafile, nt)  # skips the tether's lines
                            # above disabled in favor of averaging tether bond vector length
                            # End of timestep
                            skiplines(datafile, 8)  # skips boilerplate
                            line = datafile.readline()  # reads in last line of boilerplate to confirm we're not at EOF
                            # finished reading in and processing one timestep
                        binboys = np.arange(np.min(zvals), np.max(zvals) + 1,
                                            1.0)  # Should be average bondlength but using 1 will suffice
                        plt.hist(zvals, bins=binboys, histtype=u'step',label=discriminator)
                    #end of file reached generate the plot
plt.legend()
plt.show()