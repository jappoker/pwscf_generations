import re
import os
import subprocess
import sys
sys.path.append("..")
from generator import VC,tmcolors,opti_condition,generator

class parameter:
    PAW_hgh = "Fe.pbe-sp-hgh.UPF"
    PAW_JTH = "Fe.JTH-GGA-PBE-paw.UPF"
    PAW_GBRV = "fe_GBRV_pbe_v1.5.uspp.F.UPF"
    PAW_KS = "Fe.KS.GGA-PBE-paw.UPF"

    # PAW = PAW_JTH
    PAW = "epaw1_KS-5-23.UPF"

class get_scf_ph:

    def __init__(self,T,degauss,foldername):
        self.T = T
        self.degauss = degauss
        self.foldername = foldername
    
    def get_pressure_index(self,filename = "input.dat"):
        inputfile = os.path.join(self.foldername,filename)
        with open(inputfile,"r") as f:
            pressure_line = f.readlines()[0:1][0].replace("\r\n", "")
            pressure_range = pressure_line.split()
            f.close()
        return pressure_range

    def get_coord(self,filename,pressure):
        alat = ""
        coordi = ""
        
        with open(filename ,"r") as f:
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

    def make_scf_input(self, pressure, alat, coordi, scf_folder):
        output = "&control\ncalculation = 'scf',\n"
        output += "restart_mode = 'from_scratch',\n"
        output += "tstress = .true.\ntprnfor = .true.\nprefix = '%s'\n"%VC.PREFIX
        output += "nstep = 400\n"
        output += "pseudo_dir = '%s'\n"%VC.PAW_folder
        output += "outdir = '%s/tmp'\n" % os.path.join(
            self.foldername, scf_folder)
        output += "etot_conv_thr = 1.0D-6\nforc_conv_thr = 1.0D-5\n"
        output += "/\n&system\nibrav = 0\n"
        output += "celldm(1) =  " + str(alat) + "\n"
        output += "nat = 2\nntyp = 1\necutwfc = 90.0\n"
        output += "occupations='smearing', smearing='fd', degauss= %s,\n"%self.degauss
        output += "/\n&electrons\ndiagonalization='david'\nmixing_mode = 'plain'\n"
        output += "mixing_beta = 0.3\nconv_thr = 1.0d-10\n"
        output += "/\n&ions\nion_dynamics='damp'\n"
        output += "/\n&cell\ncell_dynamics = 'damp-w'\n"
        output += "press = " + str(pressure) + "\n"
        output += "wmass = %s\n"%VC.WMASS
        output += "/\nCELL_PARAMETERS\n" + str(coordi)
        output += "ATOMIC_SPECIES\n"
        output += "Fe   55.845   %s\n"%parameter.PAW

        atomic_position = VC.ATOMIC_POSITION.split("\n")

        output += "ATOMIC_POSITIONS (crystal)\n%s\n%s\n" %(atomic_position[1],atomic_position[2])
        output += "K_POINTS (Automatic)\n"
        output += "10  10  10   1   1   1"

        return output   

    def make_ph_input(self,phonon_folder):
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
        phonon += "outdir='%s/tmp',\n" % phonon_folder
        phonon +="tr2_ph=1.0d-14,\n"
        phonon +="alpha_mix=0.3\n"
        phonon += "nq1=%s\n"%nq1
        phonon += "nq2=%s\n" % nq2
        phonon += "nq3=%s\n" % nq3
        phonon += "/\n"

        return phonon


def make_files(scfph , pressure):                   # get scf.in & ph.in
    vc_file = "vc_%s.0/vc_%s.0.out"%(pressure,pressure)
    vc_fulldir = os.path.join(scfph.foldername, vc_file)

    scf_folder = "scf_%s.0"%pressure
    scf_input_name = "scf_%s.0.in"%pressure

    scf_fullfolder = os.path.join(scfph.foldername, scf_folder)

    ph_input_name = "ph_%s.0.in" % pressure


    if os.path.exists(vc_fulldir):
        alat, coordi = scfph.get_coord(vc_fulldir, pressure)
        scf_input = scfph.make_scf_input(
            pressure, alat, coordi, scf_folder)       # get scf input

        if os.path.exists(scf_fullfolder) == False:
            subprocess.call("mkdir %s/ ; echo '%s%s%s is created'"%(scf_folder,tmcolors.HEADER, scf_folder, tmcolors.ENDC),shell=True)
        
        with open(os.path.join(scf_fullfolder,scf_input_name), "w+") as f:      # write the scf.in
            f.write(scf_input)
            print("%s%s%s is created!"%(tmcolors.OKGREEN,scf_input_name , tmcolors.ENDC))
            f.close()
        
        ph_input = scfph.make_ph_input(scf_fullfolder)                          # get phonon input
        with open(os.path.join(scf_fullfolder, ph_input_name), "w+") as f:       # write the ph.in
            f.write(ph_input)
            print("%s%s%s is created!" %
                  (tmcolors.OKBLUE, ph_input_name, tmcolors.ENDC))
            f.close()


    else :
        print("There is no pressure =%s %s %svcrelax folder!" %
              (tmcolors.HEADER, pressure, tmcolors.ENDC))


    

    return None


def make_sbatch_files(prefix, pressure_range, scfph, run_type):
    pressure_text = ""
    for p in pressure_range:
        pressure_text += "%s "%float(p)

    job_name = "%s%s_%s"%(run_type,prefix,scfph.T)
    job_file_name = "%s.sh"%job_name

    core_per_job = int(24 / len(pressure_range))

    if run_type == "scf":
        runx = "pw"
    else:
        runx = "ph"

    job_text = "#!/bin/sh\n"
    job_text += "#SBATCH -A mphys\n"
    job_text += "#SBATCH -N 1\n"
    job_text += "#SBATCH -c 24\n"
    job_text += "#SBATCH -J %s\n"%job_name
    job_text += "#SBATCH --time=120:00:00\n"
    job_text += "#SBATCH --mail-type=ALL\n"
    job_text += "#SBATCH --mail-user=jz2907@columbia.edu\n"
    job_text += "module load intel-parallel-studio/2017\n"
    job_text += "for i in %s; do\n"%pressure_text
    job_text += "cd scf_$i\n"
    job_text += "mpirun -np %s %s.x -in %s_$i.in > %s_$i.out &\n" % (core_per_job,runx,run_type, run_type)
    job_text += "sleep 3\ncd ..\ndone\nwait"

    with open(os.path.join(scfph.foldername,job_file_name),"w+") as f:
        f.write(job_text)
        f.close()

    return None


def make_sbatch(scfph,run_type):
    
    sbatch_list = { '1':[300,320,340,360,380,400],'2':[500,600,700,800,900,1000],'3':[1100,1200,1300,1400,1500],'4':[0,50,100,150,200,250] }
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

    with open(os.path.join(scfph.foldername, job_file_name), "w+") as f:
        f.write(job_text)
        f.close()

    return None

def make_sh(scfph, run_type):
    sbatch_list = {'1': [300, 320, 340, 360, 380, 400], '2': [500, 600, 700, 800, 900, 1000], '3': [
        1100, 1200, 1300, 1400, 1500], '4': [0, 50, 100, 150, 200, 250]}
    for i in sbatch_list.keys():
        make_sh_files(i, sbatch_list[i], scfph, run_type)
    return None

def get_main():
    this_folder = os.getcwd()
    T = int(this_folder.split("/")[-1])
    opti = opti_condition(T)
    scfph = get_scf_ph(T, opti.degauss(), this_folder)

    pressure_range = scfph.get_pressure_index()      # get folder index list

    for pressure in pressure_range :        # make scf_* folder and scf.in for each pressure
        make_files(scfph, pressure)

    types = ["scf","ph"]
    for run_type in types:
        make_sbatch(scfph, run_type)
    
    types_sh = ["q2r","matdyn"]
    for type_sh in types_sh:
        make_sh(scfph, type_sh)
    

if __name__ == "__main__":
    get_main()
