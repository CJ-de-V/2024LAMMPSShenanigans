import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


# returns unit bondvectors
def normalize(v):
    norm = np.linalg.norm(v)
    global avgbondlength
    avgbondlength += norm
    if norm == 0:
        norm = np.finfo(v.dtype).eps
    return v / norm


def skiplines(f, n):
    for i in range(n):
        f.readline()


def fitfunc(x, b):
    return np.exp(-x / b)


np.set_printoptions(threshold=np.inf)
basedir = os.getcwd()
df = pd.DataFrame(
    columns=["Np", "Nt", "Kp", "Kt", "A", "lp", "Rg", "E2E", "zd", "Angle", "lb", "lt", "sd_lp", "sd_Rg", "sd_zd",
             "sd_E2E",
             "sd_dot"])
# ["Np", "Nt", "Kp", "Kt", "alp", "lp", "Rg", "E2E", "zd", "lb", "sd_lp", "sd_Rg", "sd_zd", "sd_E2E"]
p1 = np.linspace(0., 0., 3)
p2 = np.linspace(0., 0., 3)
N = 0  # Number of particles
Nb = 0  # number of bonds
Nc = 0  # number of correlation functions
bondvectors = 0
avgbondlength = 0

Np = np.array([8, 16, 32, 64, 128]) + 1  # +1 for even sides
Nt = [8, 16, 32, 64]  #
Kp = [0.25, 0.5, 1, 2, 4, 8, 16, 32, 64]  #
Kt = [0.25]
alpha = [80, 70, 60, 50, 40, 30]  #

for nop in Np:
    for nt in Nt:
        if nt > nop:  # No Tether dominated dynamics
            break
        for alp in alpha:  # alp in alpha:
            discriminator = "Np_" + str(nop) + "_Nt_" + str(nt) + "_Alpha_" + str(alp)
            workingdir = basedir + "/output/" + discriminator
            os.chdir(workingdir)
            n2 = nop // 2  # number of bondvectors inthehalf
            for kp in Kp:
                for kt in Kt:

                    discriminator = "Np_" + str(nop) + "_Nt_" + str(nt) + "_Kp_" + str(kp) + "_Kt_" + str(
                        kt) + "_Alpha_" + str(alp)

                    Nb = nop - 1
                    Nc = Nb - 1
                    bondvectors = np.zeros((Nb, 3))
                    correlationfunctions = np.zeros(Nc)
                    avgbondlength = 0
                    e2edist = np.array([])
                    bondot = np.array([])  # average distance corrleation of the two branches of the polymer
                    numavg = 0  # Number of timesteps we're averaging over, used for the normalization
                    ltet = 0  # average bond length for the tether
                    print('starting with analysis of: ' + discriminator)

                    # persistence and bond length calculation
                    with open('dump' + discriminator + '.dynamics') as datafile:
                        # Index-1 of the above indicates how many atoms separate the bondvectors
                        skiplines(datafile, 9)  # skips first boilerplate
                        line = 'liney'

                        while line != '':
                            line = (datafile.readline()).split(" ")
                            p1 = np.array([float(line[2]), float(line[3]), float(line[4])])
                            # position 1 1 is now real
                            p2 = np.array([0.0, 0.0, 0.0])
                            # reads particle 2 to N in and finds the bondvectors

                            startpos = p1  # position of leftmost bond
                            pm = -999999999
                            for i in range(1, nop):  # reading past entry 1 which we already read
                                line = (datafile.readline()).split(" ")
                                p2 = np.array([float(line[2]), float(line[3]), float(line[4])])
                                bondvectors[i - 1] = normalize(np.subtract(p2, p1))
                                p1 = p2
                                # if matches the bonding point record its Z position
                                if int(line[0]) == nop // 2 + 1:
                                    pm = p2

                            E2Evector = np.subtract(p2, startpos)
                            e2edist = np.append(e2edist, np.sqrt(np.dot(E2Evector, E2Evector)))
                            bondot = np.append(bondot, np.arccos(
                                np.dot((p2 - pm) / np.linalg.norm(p2 - pm),
                                       (-pm + startpos) / np.linalg.norm(-pm + startpos))))
                            # Tail to tail vectors since we don't need to make this harder than it has to be
                            for i in range(Nc):  # iterates over different spacings particles can have
                                runninavg = 0.0
                                for j in range(0, Nb - i):  # iterates over all legal bonds with i bonds between them
                                    runninavg += np.dot(bondvectors[j], bondvectors[j + i])  # Here be where we absed
                                correlationfunctions[i] += runninavg / (Nb - i)

                            ### Tether Shenaigans
                            # skiplines(datafile, nt)  # skips the tether's lines
                            # above disabled in favor of averaging tether bond vector length

                            line = (datafile.readline()).split(" ")
                            p1 = np.array([float(line[2]), float(line[3]), float(line[4])])
                            for i in range(1, nt):  # reading past entry 1 which we already read
                                line = (datafile.readline()).split(" ")
                                p2 = np.array([float(line[2]), float(line[3]), float(line[4])])  # next monomer position
                                # add bond length to total and carry on.
                                ltet += np.linalg.norm(p1 - p2)
                                p1 = p2

                            # End of timestep
                            skiplines(datafile, 8)  # skips boilerplate
                            line = datafile.readline()  # reads in last line of boilerplate to confirm we're not at EOF
                            numavg += 1
                            # finished reading in and processing one timestep

                        # Averaging over time/monomers
                        avgbondlength = avgbondlength / (Nb * numavg)
                        e2e = np.mean(e2edist)
                        sd_e2e = np.sqrt(np.var(e2edist))
                        dot = np.mean(bondot)
                        sd_dot = np.sqrt(np.var(bondot))
                        ltet = ltet / ((nt - 1) * numavg)

                        # fitting for lp
                        y = correlationfunctions / numavg
                        x = np.arange(len(y)) * avgbondlength
                        # x range in LJ units
                        weights = np.reciprocal(np.arange(Nc + 0.0, 0.0, -1.0))
                        param, param_cov = curve_fit(fitfunc, x, y, maxfev=1000, sigma=weights)
                        lp = param[0]
                        sd_lp = np.sqrt(param_cov[0][0])

                    # RG data
                    with open('radius_of_gyration' + discriminator + '.dat') as datafile:
                        fulldat = datafile.readlines()  # these files are short enough for us to just use readlines
                        fulldat.pop(0)
                        for i in range(len(fulldat)):
                            fulldat[i] = float(fulldat[i].split(' ')[1])
                        Rg = np.mean(fulldat)
                        sd_Rg = np.var(fulldat)

                    # Zdistance data

                    with open('dump' + discriminator + '.dynamics') as datafile:
                        skiplines(datafile, 7)
                        zpeak = float(datafile.readline().split(" ")[1]) - 1
                        # this part is the one that'll raise the most eyebrows, given that the -1 is practically a
                        # magic number finds walls Z pos, the -1 is to account for the fact that our wall is 1 unit
                        # away from the boundary

                    with open('pcom' + discriminator + '.dat') as datafile:
                        fulldat = datafile.readlines()[2:]  # pops 2 first lines
                        for i in range(len(fulldat)):
                            fulldat[i] = float(fulldat[i].split(' ')[3])
                        numzdata = np.array(fulldat)
                        numzdata = zpeak - numzdata
                        zdis = np.mean(numzdata)
                        sd_zd = np.sqrt(np.var(numzdata))

                    entry = [nop, nt, kp, kt, alp,
                             lp, Rg, e2e, zdis, dot,
                             avgbondlength, ltet, sd_lp, sd_Rg, sd_zd, sd_e2e, sd_dot]
                    datline = np.asarray(entry, dtype=float)
                    df.loc[len(df)] = datline
print(df, '\n')
os.chdir(basedir)
df.to_csv('data.csv')
# print you pandas here
