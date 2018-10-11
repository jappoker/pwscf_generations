import re
import os
import subprocess
from latticepredict import return_popt


class tmcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class VC:
    ## for VCRELAX
    PAW_folder = "/rigel/home/jz2907/pseudo"
    # PAW = "Fe.KS.GGA-PBE-paw.UPF"
    # PAW = "Fe.JTH-GGA-PBE-paw.UPF"
    # PAW = "fe_GBRV_pbe_v1.5.uspp.F.UPF"
    PAW = "epaw1_KS-5-23.UPF"
    PREFIX = "hcpFe"
    WMASS = "0.01"
    NUMB_of_atoms = "2"
    NUM_of_atomic_types = "1"
    E_CUTOFF = "90.0"
    SMEARING = "fd"
    K_PO = " 10 10 10"
    SHIFT = "1 1 1"
    ATOMIC_POSITION = "Atomic positions: \nFe   0.333333333333333   0.666666666666667   0.250000000000000\nFe   0.666666666666667   0.333333333333333   0.750000000000000"


class opti_condition:

    def __init__(self,T):
        self.T = T
    
    def degauss(self):
        return 0.0380017024103152 * self.T / 6000

    def folderfullname(self):
        folder = "%s" % self.T
        fulldir = os.path.join(os.getcwd(), folder)
        return fulldir

    def foldername(self):
        return "%s" % self.T


def scpfiles(foldername):
    copied_folder = "tobecopy"
    subprocess.call("scp %s/* %s/"%(copied_folder,foldername),shell=True)
    print("basic files for "+ tmcolors.HEADER+"%s" % foldername+ tmcolors.ENDC + " is copied")
    return None

class generator:

    def __init__(self,T,degauss,foldername):
        self.T = T
        self.degauss = degauss
        self.foldername = foldername

    def vc_qe_input(self,filename = "qe_input_data"):
        fulldir = os.path.join(self.foldername , filename)
        infile = "pseudopotential directory:%s \n"%VC.PAW_folder
        infile += "%s \n"%VC.PAW
        infile += "prefix: %s\n"%VC.PREFIX
        infile += "wmass: %s\n"%VC.WMASS
        infile += "number of atoms: %s\n"%VC.NUMB_of_atoms
        infile += "number of atomic types: %s\n"% VC.NUM_of_atomic_types
        infile += "energy cutoff: %s\n" %VC.E_CUTOFF
        infile += "smearing: %s\n"%VC.SMEARING
        infile += "occupations: smearing \ndegauss: %s\n"%self.degauss
        infile += "k-points: %s\n"%VC.K_PO
        infile += "shift: %s\n" %VC.SHIFT
        infile += "scratch folder: %s/tmp\n" %self.foldername
        infile += "%s\n"%VC.ATOMIC_POSITION
        with open(fulldir,"w") as f:
            f.write(infile)
            f.close()
        
        return None

    def vc_job_header(self,filename = "job_header"):
        fulldir = os.path.join(self.foldername, filename)
        infile = "Scheduler: SLURM\n"
        infile += "Necessary modules to be loades(if any): module load intel-parallel-studio/2017:\n"
        infile += "Number of nodes: 1\n"
        infile += "Number of processors: 24\n"
        infile += "#!/bin/sh\n"
        infile += "#SBATCH -A mphys\n"
        infile += "#SBATCH -N 1\n"
        infile += "#SBATCH -c 24\n"
        infile += "#SBATCH -J %s_%s\n"%(VC.PREFIX,self.T)
        infile += "#SBATCH --time=12:00:00\n"
        infile += "#SBATCH --mail-type=ALL\n"
        infile += "#SBATCH --mail-user=jz2907@columbia.edu"
        with open(fulldir, "w") as f:
            f.write(infile)
            f.close()
    
    def vc_input(self,filename = "input.dat"):
        fulldir = os.path.join(self.foldername, filename)
        infile = "-10 0 50 100 150 200 250 300 320 340 360 380 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500\n"
        infile += return_popt(self.T) +"\n"
        infile += "  1.000000000000000   0.000000000000000   0.000000000000000\n"
        infile += " -0.500000000000000   0.8660254037844386  0.000000000000000\n"
        infile += "  0.000000000000000   0.000000000000000   1.604000000000000\n"

        with open(fulldir,"w") as f:
            f.write(infile)
            f.close



def main():
    # T_range = [300]
    T_range = [300,1000,2000,3000,4000,5000,6000,7000,8000]

    for T in T_range:
        opti = opti_condition(T)            #get optimized temperature to generate vcrelax workflow & getscf, ph
        
        if os.path.exists(opti.folderfullname()) == False:   # create folder for one Temperature
            subprocess.call("mkdir %s"%opti.foldername(),shell = True)

        scpfiles(opti.foldername())                         # copy basic files for workflow & phonons

        # gen = generator(T,opti.degauss(),opti.folderfullname())
        # gen.vc_qe_input()
        # gen.vc_job_header()
        # gen.vc_input()

    return None

if __name__ == "__main__":
    main()
