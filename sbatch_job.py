import re
import os
import subprocess

import sys

from generator import tmcolors


def show_help():

    job_type = "%sjob_type%s"%(tmcolors.HEADER,tmcolors.ENDC)
    temperature = "%stemperature%s" % (tmcolors.HEADER, tmcolors.ENDC)
    job_index = "%sjob_index%s" % (tmcolors.HEADER, tmcolors.ENDC)

    print("HELP!")
    print("argument: [job_type][temperature][job_index]")
    print("\t--%s\tscf, ph, q2r, matdyn, vc"%job_type)
    print("\t--%s\te.g. 300, 2000, .."%temperature)
    print("\t--%s"%job_index)
    print("\t\t 1 : [300,320,340,360,380,400]")
    print("\t\t 2 : [500,600,700,800,900,1000]")
    print("\t\t 3 : [1100,1200,1300,1400,1500]")
    print("\t\t 4 : [0,50,100,150,200,250]")
    return None

def argv_valid(arguments):
    status = False
    if len(arguments)!= 4:
        status = False
    elif arguments[1] not in ["scf","ph","matdyn","q2r","vc"]:
        status = False
    elif arguments[3] not in ["1","2","3","4","all"]:
        status = False
    else:
        status = True

    return status


def sbatch_file(arguments):
    file_name = ""
    temperature = arguments[2]
    this_folder = os.path.dirname(os.path.abspath(__file__))
    folder_name = os.path.join(this_folder,temperature)

    prefix = arguments[1]
    job_index = arguments[3]

    sbatch_list = {'1': [300, 320, 340, 360, 380, 400], '2': [500, 600, 700, 800, 900, 1000], '3': [
        1100, 1200, 1300, 1400, 1500], '4': [0, 50, 100, 150, 200, 250]}
    
    pressure_range = sbatch_list[job_index]
    print("Job type is %s . The Running Condition is: T = %s "%(prefix, temperature) + " Pressure Range = %s "%pressure_range)
    file_name = "%s%s_%s.sh"%(prefix,job_index,temperature)
    
    return folder_name, file_name

def sbatch_it(folder_name, file_name):
    subprocess.call("cd %s ;qsub %s" %
                    (folder_name, file_name), shell=True)
    # subprocess.call("cd %s ;cat %s" % (folder_name, file_name), shell=True)
    return None

def sh_it(folder_name, file_name):
    subprocess.call("cd %s ;sh %s" % (folder_name, file_name), shell=True)
    return None

def make_sure():
    if input("Print y to continue: ") == "y":
        status = True
    else:
        status = False
    return status

if __name__ == "__main__":
    last_argu = (sys.argv[-1])
    if last_argu == "-h" or not argv_valid(sys.argv) :
        show_help()
    else:
        if sys.argv[3] == "all":
            for i in [1,2,3,4]:
                new_argv = sys.argv[:-1]
                new_argv.append(str(i))
                folder_name, file_name = sbatch_file(new_argv)
                full_name = os.path.join(folder_name, file_name)
                if os.path.exists(full_name):
                    # if make_sure():
                    if True:
                        if new_argv[1] in ["ph", "scf", "vc"]:
                            sbatch_it(folder_name, file_name)
                        if new_argv[1] in ["matdyn", "q2r"]:
                            sh_it(folder_name, file_name)
                else:
                    print(tmcolors.WARNING+"SBATCH FILE DOES NOT EXIST!"+tmcolors.ENDC)

        else:
            folder_name, file_name = sbatch_file(sys.argv)
            full_name = os.path.join(folder_name, file_name)
            if os.path.exists(full_name):
                # if make_sure():
                if True:
                    if sys.argv[1] in ["ph","scf","vc"]:
                        sbatch_it(folder_name, file_name)
                    if sys.argv[1] in ["matdyn","q2r"]:
                        sh_it(folder_name, file_name)
            else:
                print(tmcolors.WARNING+"SBATCH FILE DOES NOT EXIST!"+tmcolors.ENDC)
