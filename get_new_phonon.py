import re
import os
import subprocess
import sys
from generator import VC, tmcolors
from newcheck import color

'''
read input from vc_* to make scf_*
'''

T_range = [8000]
# T_range = [300,3000,4000,5000,6000,7000,8000]
scratch_folder = "/scratch/sciteam/zhuang1"
this_scratch = "epaw1"
# T_range = [300]

class scfph:
    def __init__(self,T,P):
        self.T = T
        self.P = P

        self.T_folder = "%s"%T
        self.vc_folder = "vc_%s.0"%self.P
        self.vc_fdir = os.path.join(self.T_folder,self.vc_folder)
        self.vc_input_name = "vc_%s.0.in" % self.P
        self.vc_indir = os.path.join(self.vc_fdir,self.vc_input_name)
        self.vc_input, 
        self.vc_output_name = "vc_%s.0.out" % self.P
        self.vc_outdir = os.path.join(self.vc_fdir, self.vc_output_name)
        self.vc_output,
        self.vc_alat,
        self.vc_coord,

        self.scf_folder = "scf_%s.0" % self.P
        self.scf_fdir = os.path.join(self.T_folder,self.scf_folder)
        self.scf_input_name = "scf_%s.0.in"%P
        self.scf_indir = os.path.join(self.scf_fdir, self.scf_input_name)
        self.scf_input, 


        self.ph_fdir = self.scf_fdir
        self.ph_input_name = "ph_%s.0.in" % self.P
        self.ph_indir = os.path.join(self.ph_fdir, self.ph_input_name)
        self.ph_input, 
        self.scratch_fdir = os.path.join(scratch_folder,this_scratch)
        self.scratch_dir = os.path.join(self.scratch_fdir, self.scf_fdir)

    @property
    def vc_input(self):
        '''
        return a list of lines
        '''
        lines = []
        if os.path.exists(self.vc_indir):
            with open(self.vc_indir,"r") as f:
                lines = f.readlines()
        else:
            print(color("No such vc!").warning)
        return lines

    @property
    def vc_output(self):
        '''
        return a list of lines
        '''
        lines = []
        if os.path.exists(self.vc_outdir):
            with open(self.vc_outdir, "r") as f:
                lines = f.readlines()
        else:
            print(color("No such vc output!").warning)
        return lines
    

    @property
    def vc_alat(self):
        alat,_ = self.vc_get_para()
        return alat

    @property
    def vc_coord(self):
        _,coord = self.vc_get_para()
        return coord

    def vc_get_para(self):
        file_lines = self.vc_output
        for i in range(len(file_lines)):
            line = file_lines[i]
            if re.findall("Begin final coordinates", line):
                for _ in range(1):  # volune
                    i += 1
                    line = file_lines[i]
                    volume = float(line.split()[4])
                for _ in range(2):
                    i += 1
                for _ in range(1):  # alat
                    i += 1
                    line = file_lines[i]
                    alat = (line.split()[2].rstrip(")"))

                for _ in range(1):  # coord
                    i += 1
                    line = file_lines[i]
                    coorda = float(line.split()[0])

                    i += 1
                    line = file_lines[i]
                    coordb1 = float(line.split()[0])
                    coordb2 = float(line.split()[1])

                    i += 1
                    line = file_lines[i]
                    coordc = float(line.split()[-1])
                break

        coord = "%14.9f%14.9f%14.9f\n" % (coorda, 0, 0)
        coord += "%14.9f%14.9f%14.9f\n" % (coordb1, coordb2, 0)
        coord += "%14.9f%14.9f%14.9f\n" % (0, 0, coordc)

        return alat , coord

    def make_tmp_folder(self):
        subprocess.call("FILE=%s;if [-f $FILE];then;echo 'File $FILE exists.';else;echo 'File $FILE does not exist.';mkdir $FILE;fi"%(self.scratch_dir), shell=True)
    
    @property
    def scf_input(self):
        raw = self.vc_input
        input = ""
        status = False
        
        for line in raw:
            if re.findall("calculation",line):
                input += "calculation = 'scf',\n"
            
            elif re.findall("pseudo",line):
                input += "pseudo_dir = '/u/sciteam/zhuang1/pseudo'\n"

            elif re.findall("tmp",line):
                input += "outdir = '%s'\n" % (self.scratch_dir)

            elif re.findall("celldm(1)", line):
                input += ("celldm(1) = %s\n" % round(self.vc_alat,7))

            elif re.findall("conv_thr = 1.0d-12", line):
                input += "conv_thr = 1.0D-10\n"
                # print("F@@@@@@@@@CK")
                
            elif re.findall("0.000000000", line) and re.findall("Fe",line) == [] and status == True:
                pass

            elif re.findall("electron_maxstep",line):
                pass

            elif re.findall("CELL_PARAMETERS", line):
                input += "CELL_PARAMETERS\n%s" % self.vc_coord
                status = True
            else:
                input += line
        return input

    @property
    def ph_input(self):
        nq1 = 4
        nq2 = 4
        nq3 = 4
        phonon = "phonons at on grid\n"
        phonon += "&inputph\n"
        phonon +="prefix='%s',\n"%VC.PREFIX
        phonon +="ldisp=.true.,\n"
        phonon +="iverbosity=1,\n"
        phonon +="epsil=.false.,\n"
        phonon +="recover=.true.,\n"
        phonon +="fildyn='dyn',\n"
        phonon +="amass(1)=55.845,\n"
        phonon += "outdir='%s',\n"%self.scratch_dir
        phonon +="tr2_ph=1.0d-14,\n"
        phonon +="alpha_mix=0.3\n"
        phonon += "nq1=%s\n" % nq1
        phonon += "nq2=%s\n" % nq2
        phonon += "nq3=%s\n" % nq3
        phonon += "/\n"
        return phonon

    
    def anpai(self):
        if os.path.exists(self.scf_fdir)==False:
            os.mkdir(self.scf_fdir)
        with open(self.scf_indir,"w+") as f:
            f.write(self.scf_input)
        with open(self.ph_indir,"w+")as f: 
            f.write(self.ph_input)
        self.make_tmp_folder()
        print("inputs in\t%s\tis created"%color("%s"%self.scf_fdir).header)


def get_p_from_folder(folder):
    pressure = float(folder.split("_")[-1])
    pressure = int(pressure)
    return pressure


def make_sbatch_files(prefix, pressure_range, scfph, run_type):
    pressure_text = ""
    for p in pressure_range:
        pressure_text += "%s " % float(p)

    job_name = "%s%s_%s" % (run_type, prefix, scfph.T)
    job_file_name = "%s.sh" % job_name

    core_per_job = int(24 / len(pressure_range))
    core_per_job =  32

    if run_type == "scf":
        runx = "pw"
    else:
        runx = "ph"





# cd / u/sciteam/zhuang1/work/10-10-try/try2/
# export OMP_NUM_THREADS = 2
# aprun - n 4 pw.x - in scf.in > scf.out

    job_text = "#!/bin/sh\n"
    job_text += "#PBS -j oe\n"
    job_text += "#PBS -l nodes=1:ppn=32:xe\n"
    job_text += "#PBS -q normal\n"
    job_text += "#PBS -N %s\n" % job_name
    job_text += "#PBS -l walltime=24:00:00\n"
    job_text += "#PBS -m bea\n"
    job_text += "#PBS -M jz2907@columbia.edu\n"
    # job_text += "module load intel-parallel-studio/2017\n"
    job_text += "cd /u/sciteam/zhuang1/work/epaw1/aaaaa/%s"%scfph.T_folder
    job_text += "export OMP_NUM_THREADS = 2"
    job_text += "for i in %s; do\n" % pressure_text
    job_text += "cd scf_$i\n"
    job_text += "aprun -n %s %s.x -in %s_$i.in > %s_$i.out \n" % (
        core_per_job, runx, run_type, run_type)
    job_text += "sleep 3\ncd ..\ndone\nwait"

    with open(os.path.join(scfph.T_folder, job_file_name), "w+") as f:
        f.write(job_text)
        f.close()

    return None


def make_sbatch(scfph, run_type):

    sbatch_list = {'1': [300, 320, 340, 360, 380, 400], '2': [500, 600, 700, 800, 900, 1000], '3': [
        1100, 1200, 1300, 1400, 1500], '4': [0, 50, 100, 150, 200, 250]}
    for i in sbatch_list.keys():
        make_sbatch_files(i, sbatch_list[i], scfph, run_type)
    return None


def make_sh_files(prefix, pressure_range, scfph, run_type):
    pressure_text = ""
    for p in pressure_range:
        pressure_text += "%s " % float(p)

    job_name = "%s%s_%s" % (run_type, prefix, scfph.T)
    job_file_name = "%s.sh" % job_name

    core_per_job = int(24 / len(pressure_range))

    job_text = "#!/bin/sh\n"
    job_text += "#SBATCH -A mphys\n"
    job_text += "#SBATCH -N 1\n"
    job_text += "#SBATCH -c 24\n"
    job_text += "#SBATCH -J %s\n" % job_name
    job_text += "#SBATCH --time=120:00:00\n"
    job_text += "#SBATCH --mail-type=ALL\n"
    job_text += "#SBATCH --mail-user=jz2907@columbia.edu\n"
    job_text += "module load intel-parallel-studio/2017\n"
    job_text += "for i in %s; do\n" % pressure_text
    job_text += "cd scf_$i\n"
    job_text += "mpirun -np %s %s.x -in \"$(ls ../%s.in)\" > %s.out &\n" % (
        core_per_job, run_type, run_type, run_type)
    job_text += "sleep 3\ncd ..\ndone\nwait"

    with open(os.path.join(scfph.T_folder, job_file_name), "w+") as f:
        f.write(job_text)
        f.close()

    return None


def make_sh(scfph, run_type):
    sbatch_list = {'1': [300, 320, 340, 360, 380, 400], '2': [500, 600, 700, 800, 900, 1000], '3': [
        1100, 1200, 1300, 1400, 1500], '4': [0, 50, 100, 150, 200, 250]}
    for i in sbatch_list.keys():
        make_sh_files(i, sbatch_list[i], scfph, run_type)
    return None

def anpaiyixia(now):
    now.anpai()
    # print(now.scf_input + "\n")


    return None

def main():
    for T in T_range:
        file_list = os.listdir("%s" % T)
        for file in file_list:
            if re.findall("vc_", file):
                P = get_p_from_folder(file)
                now = scfph(T, P)
                anpaiyixia(now)

        types = ["scf", "ph"]
        for run_type in types:
            make_sbatch(now, run_type)

        types_sh = ["q2r", "matdyn"]
        for type_sh in types_sh:
            make_sh(now, type_sh)


if __name__ == "__main__":
    main()
