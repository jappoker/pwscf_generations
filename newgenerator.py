import re
import os
import subprocess
import sys

from generator import tmcolors,VC
# from lazy_property import LazyProperty

# T_range = [300, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]
T_range = [1000,2000]
P_range = [0, 50, 100, 150, 200, 250, 300, 320, 340, 360, 380, 400,
           500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]
# P_range = [0,50,100]

class vcrelax:
    def __init__(self,T,P):
        self.T = T
        self.P = P
        self.degauss = 0.0380017024103152 * self.T / 6000
        self.origin_dir = self.get_origin_dir
        self.target_dir = self.get_target_dir
        self.origin_output_fname = self.get_origin_output_fname
        self.target_input_fname = self.get_target_input_fname

    @property
    def get_full_dir(self):
        path = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), "%s" % self.T)
        return path

    @property
    def get_target_dir(self):
        path = self.get_full_dir
        path = os.path.join(path, self.get_vc_name)
        return path
    
    @property
    def get_vc_name(self):
        return "vc_%s.0"%self.P

    @property
    def get_origin_dir(self):
        path = self.get_full_dir
        path = os.path.join(path, self.get_vc_name)
        return path.replace("newpaw", "explanet")

    @property
    def get_origin_output_fname(self):
        fname = self.get_vc_name
        fullname = "%s/%s.out"%(self.origin_dir,fname)
        return fullname

    @property
    def get_target_input_fname(self):
        fname = self.get_vc_name
        fullname = "%s/%s.in" % (self.target_dir, fname)
        return fullname
    
    @property
    def get_alat_coord(self):
        alat = ""
        coordi = ""
        with open(self.get_origin_output_fname) as f:
            switch = 0
            data = []
            # get new coordinates and cell para
            for line in f:
                if line.find("Begin final coordinates") != -1:
                    switch = 1
                    continue
                if line.find("End final coordinates") != -1:
                    switch = 0
                if switch == 1:
                    data.append(line)

            alat = data[3].replace(
                "CELL_PARAMETERS (alat=", "").replace(")", "").replace("\n", "")
            coordi = data[4] + data[5] + data[6]
            f.close()

        return alat, coordi

def make_vc_file(vc,alat,coord):
    if vc.P <= 100:
        alat = float(alat) * 0.99945
    else:
        alat = float(alat) * (1- (8884/543 + vc.P**2/444/888)/10000)

    output = "&control\ncalculation = 'vc-relax',\n"
    output += "restart_mode = 'from_scratch',\n"
    output += "tstress = .true.\ntprnfor = .true.\nprefix = '%s'\n" % VC.PREFIX
    output += "nstep = 800\ndt = 15\n"
    output += "pseudo_dir = '%s'\n" % VC.PAW_folder
    output += "outdir = './tmp'\n"
    output += "etot_conv_thr = 1.0D-7\nforc_conv_thr = 1.0D-6\n"
    output += "/\n&system\nibrav = 0\nnbnd = 28\n"
    output += "celldm(1) =  " + str(alat) + "\n"
    output += "nat = 2\nntyp = 1\necutwfc = %s , ecutrho = 450\n"%VC.E_CUTOFF
    output += "occupations='smearing', smearing='fd', degauss= %s,\n" % vc.degauss
    output += "nr1 = 48\nnr2 = 48\nnr3 = 64\n"
    output += "/\n&electrons\ndiagonalization='david'\nmixing_mode = 'plain'\n"
    output += "electron_maxstep = 2222\n"
    output += "mixing_beta = 0.11\nconv_thr = 1.0d-12\n"
    output += "/\n&ions\nion_dynamics='damp'\n"
    output += "/\n&cell\ncell_dynamics = 'damp-w'\n"
    output += "press = " + str(vc.P*10) + "\n"
    output += "wmass = 0.01\n"
    output += "/\nCELL_PARAMETERS\n" + str(coord)
    output += "ATOMIC_SPECIES\n"
    output += "Fe   55.845   %s\n" % VC.PAW
    # 55.845
    atomic_position = VC.ATOMIC_POSITION.split("\n")

    output += "ATOMIC_POSITIONS (crystal)\n%s\n%s\n" % (
        atomic_position[1], atomic_position[2])
    output += "K_POINTS (Automatic)\n"
    output += "10  10  10   1   1   1"

    return output


def make_sbatch_files(prefix, pressure_range, vc, run_type):
    pressure_text = ""
    for p in pressure_range:
        pressure_text += "%s " % float(p)

    job_name = "%s%s_%s" % (run_type, prefix, vc.T)
    job_file_name = "%s.sh" % job_name

    core_per_job = int( 24 / len(pressure_range) )
    # core_per_job = 1

    if run_type == "vc":
        runx = "pw"

    job_text = "#!/bin/sh\n"
    job_text += "#SBATCH -A mphys\n"
    job_text += "#SBATCH -N 1\n"
    job_text += "#SBATCH -c 24\n"
    job_text += "#SBATCH -J %s\n" % job_name
    job_text += "#SBATCH --time=4:00:00\n"
    job_text += "#SBATCH --mail-type=ALL\n"
    job_text += "#SBATCH --mail-user=jz2907@columbia.edu\n"
    job_text += "module load intel-parallel-studio/2017\n"
    job_text += "for i in %s; do\n" % pressure_text
    job_text += "cd vc_$i\n"
    job_text += "mpirun -np %s %s.x -in %s_$i.in > %s_$i.out &\n" % (
        core_per_job, runx, run_type, run_type)
    job_text += "sleep 3\ncd ..\ndone\nwait"

    with open(os.path.join(vc.get_full_dir, job_file_name), "w+") as f:
        f.write(job_text)
        f.close()

    return None


def make_sbatch(vc, run_type):

    sbatch_list = {'1': [300, 320, 340, 360, 380, 400], '2': [500, 600, 700, 800, 900, 1000], '3': [
        1100, 1200, 1300, 1400, 1500], '4': [0, 50, 100, 150, 200, 250]}
    for i in sbatch_list.keys():
        make_sbatch_files(i, sbatch_list[i], vc, run_type)
    return None

if __name__ == "__main__":
    for T in T_range:
        for P in P_range:
            vc = vcrelax(T,P)
            if os.path.exists(vc.origin_output_fname):
                alat,coord = vc.get_alat_coord
                file_content = make_vc_file(vc,alat,coord)
                
                if os.path.exists(vc.target_dir) == False:
                    subprocess.call( "mkdir %s"%vc.target_dir ,shell=True)
                    print("%s%s%s is created!"%(tmcolors.HEADER,vc.target_dir,tmcolors.ENDC))
                
                if os.path.exists(vc.target_input_fname) == False:
                    with open(vc.target_input_fname,"w") as f:
                        f.write(file_content)
                        print("%s%s.in%s is created!" %
                            (tmcolors.HEADER, vc.get_vc_name, tmcolors.ENDC))
                        f.close()

        # make_sbatch(vc,"vc")
            
