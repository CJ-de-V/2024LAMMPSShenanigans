# Batchrun.py: runs a batch of LAMMPS experiments for the tethered version
# arguments: NONE
# requirements: in.setup, in.continue and a.out should be present in the run directory

import os
from mpi4py import MPI
from lammps import lammps
from subprocess import call
import numpy as nu
import math

me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
basedir = os.getcwd()
Np = nu.array([8, 16, 32, 64, 128]) + 1  # +1 for even sides #
Nt = [8, 16, 32, 64, 128]
Kp = [64, 32, 16, 8, 4, 2, 1, 0.5, 0.25]  # [64,62,16,8,4,2,1,0.5,0.25]
Kt = [0.25]
Anglesd = nu.array([90, 80, 70, 60, 50, 40, 30])  # conic angles 80, 70, 60, 50, 40, 30
Anglesr = Anglesd * math.pi / 180
lb = 0.95  # Starting separation of the bonds
for np in Np:
    for nt in Nt:
        if nt > np:  # No Tether dominated dynamics, excluding the equal length case since it blows stuff up?
            break
        for i in range(len(Anglesr)):

            discriminator = "Np_" + str(np) + "_Nt_" + str(nt) + "_Alpha_" + str(Anglesd[i])
            workingdir = basedir + "/output/" + discriminator
            os.makedirs(workingdir, exist_ok=True)
            os.chdir(workingdir)
            call(["../../a.out", str(np), str(nt), str(Anglesd[i]), "config.tether"])  # Make Config File
            lmp = lammps()
            lines = open(basedir + '/' + "in.startup", 'r').readlines()
            for line in lines:
                lmp.command(line)

            spacing = 1.13  # spacing slightly greater than wall's interaction range to prevent it from touching at the bottom
            maxbl = 1.5  # Maximum Bondlength
            lenience = 1.0  # just to shift it so that the polymers aren't on top of the wall, also if they're too close the wall
            # at the concave part of the cone can sever the FENE bond
            ###CONE IS FIXED WITH ITS TIP AT THE TOP OF THE BOX
            l_max = ((np - 1) / 2 + (nt - 1)) * maxbl + spacing + lenience  # Height of Cone
            lmp.command("print " + "\">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" + str(l_max) + "\"")
            r_lo = l_max * math.tan(Anglesr[i])  # radius of bottom side of cone to get desired alpha
            zhi = (np - 2.0) / 2.0 + lenience
            # Lenience is to prevent the thing from blowing up
            zlo = zhi - l_max
            conecommands = [
                "region myCone cone z 0 0 " + str(r_lo) + " 0 " + str(zlo) + " " + str(zhi),
                "fix confinedCone all wall/region myCone lj126 1.0 1.0 1.122462048"]
            lmp.commands_list(conecommands)

            for kp in Kp:
                for kt in Kt:
                    # Equilibrate and adjust to new parameters
                    rundistinction = ["angle_coeff   1  " + str(kp),  # whatever distinguishes it from the next run.
                                      "angle_coeff   2  " + str(kt)]  # ,   "run 50000"]
                    # some unrecorded runs to let it adjust to the new run parameters should be here
                    lmp.commands_list(rundistinction)
                    discriminator = "Np_" + str(np) + "_Nt_" + str(nt) + "_Kp_" + str(kp) + "_Kt_" + str(
                        kt) + "_Alpha_" + str(Anglesd[i])
                    lmp.command("run 20000")

                    rundetails = (
                            "print " + "\"" + ">>>>>>>>>>>>>>>>DETAILS OF RUN<<<<<<<<<<<<<<<<" + "\n>>>>>>>>>>>>>>>>"
                            + discriminator + "<<<<<<<<<<<<<<<<\"")
                    lmp.command(rundetails)
                    # Commands needed for both setup and restart runs
                    primersetup = ['log log.' + discriminator,
                                   "fix mythermofile all print 10000 \"$t ${mytemp} ${myepair}\" file thermo_output" + discriminator + ".dat screen no",
                                   # thermodynamic data outputted to the file appropriately named (sort of)
                                   "fix myRG2file all print 10000 \"$t ${RG2}\" file radius_of_gyration" + discriminator + ".dat screen no",
                                   "fix comfix polymer ave/time 1 1 10000 c_mycomcompute[*] file pcom" + discriminator + ".dat",
                                   # every 10 000 timesteps spits out COM of polymer
                                   "dump dum2 all custom 10000 dump" + discriminator + ".dynamics id type x y z",
                                   "dump_modify dum2  sort id"]
                    lmp.commands_list(primersetup)

                    lmp.command("run 50000")

                    cleansetup = ["write_restart restart." + discriminator,
                                  "unfix mythermofile",
                                  "unfix myRG2file",
                                  "unfix comfix",
                                  "undump dum2",
                                  "log log.dummylog"
                                  # Since I can't actually turn off the logfile I'm just gonna have it write all the trash about resizing to this one, it'll overwrite itself and be small
                                  ]
                    lmp.commands_list(cleansetup)

            lmp.close()  # End of NxN run

print("Proc %d out of %d procs has" % (me, nprocs), lmp)
