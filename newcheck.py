import re
import os
import subprocess
import sys

from generator import tmcolors
from typing import List

T_range = [300, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]
# T_range = [300, 1000, 2000, 3000, 4000]
P_range = [0, 50, 100, 150, 200, 250, 300, 320, 340, 360, 380, 400,
           500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]

class color:
    def __init__(self,content):
        self.content = content
    
    @property
    def header(self):
        return "%s%s%s"%(tmcolors.HEADER,self.content,tmcolors.ENDC)

    @property
    def warning(self):
        return "%s%s%s" % (tmcolors.WARNING, self.content, tmcolors.ENDC)

    @property
    def bold(self):
        return "%s%s%s" % (tmcolors.BOLD, self.content, tmcolors.ENDC)

    @property
    def green(self):
        return "%s%s%s" % (tmcolors.OKGREEN, self.content, tmcolors.ENDC)
    
    @property
    def blue(self):
        return "%s%s%s" % (tmcolors.OKBLUE, self.content, tmcolors.ENDC)

    @property
    def underline(self):
        return "%s%s%s" % (tmcolors.UNDERLINE, self.content, tmcolors.ENDC)


def info_build(prefix,status):
    if len(status)  < 1:
        out = ""
    else:
        out = "%s : %s" % (color(prefix).underline, status)
    # if re.findall("End",status):
        # out = ""
    return out

class check_folder:
    def __init__(self, T, P):
        self.T = T
        self.P = P
        self.cP = color(P).blue
        self.full_dir = self.get_full_dir
        self.cg_file = self.get_cg_file
        self.vc_file = self.get_vc_file
        self.scf_file = self.get_scf_file
        self.ph_folder = self.get_ph_folder

    
    @property
    def get_full_dir(self):
        path = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(path, "%s"%self.T)

    @property
    def get_cg_file(self):
        cg_file = "cg_%s.0.out"%self.P
        folder = os.path.join(self.full_dir, "cg_%s.0" % self.P)
        return os.path.join(folder,cg_file)

    @property
    def get_vc_file(self):
        vc_file = "vc_%s.0.out" % self.P
        folder = os.path.join(self.full_dir, "vc_%s.0" % self.P)
        return os.path.join(folder, vc_file)

    @property
    def get_scf_file(self):
        scf_file = "scf_%s.0.out" % self.P
        folder = os.path.join(self.full_dir, "scf_%s.0" % self.P)
        return os.path.join(folder, scf_file)

    @property
    def get_ph_folder(self):
        folder = os.path.join(self.full_dir, "scf_%s.0" % self.P)
        return folder

    @property
    def folder_report(self):
        full = "-- P = %s" % self.cP
        return full

    @property
    def check_cg(self):
        full = ""
        prefix = "cg"
        if os.path.exists(self.cg_file):
            with open(self.cg_file,"r") as f:
                for line in f:
                    if re.findall("convergence NOT achieved",line):
                        full = color("Not Converged").warning
                        break
                    if re.findall("JOB DONE", line):
                        full = color("End").bold
                        break
                f.close()
        return info_build(prefix, full)
    
    @property
    def ethr(self):
        tmp_ethr = ""
        tmp_ethr_p = []
        aaaa = ""
        regex = r"\s*\d.\d+E-\d+\s*"
        if os.path.exists(self.vc_file):
            with open(self.vc_file, "r") as f:
                fline = f.readlines()
                f.close()
                for line in fline:
                    if re.findall("ethr", line):
                        tmp_ethr = line.split()[2].rstrip(",")
                    if re.match(regex,tmp_ethr):
                        tmp_ethr_p.append  (int(tmp_ethr.split("-")[-1]))
                        # tmp_ethr_p.append(float(tmp_ethr))
            bbbb = "▁▂▃▄▅▆▇█▉▊▋▌▍▎▏"
            tol = len(tmp_ethr_p)
            m = int(tol / 20)
            for i in range(m):
                tttt = tmp_ethr_p[i*20]
                aaaa += "%s"%bbbb[tttt-1]
            tttt = tmp_ethr_p[-1]
            aaaa += "%s" % bbbb[tttt-1]
        return aaaa

    @property
    def check_vc(self):
        full = ""
        prefix = "vc"
        status = False
        tmp_pressure  = " "*7
        tmp_ethr = ""
        origin_volume = ""
        new_volume = ""
        tmp_ethr_p = ""
        regex = r"\s*\d.\d+E-\d+\s*"
        # max_iteration = 0
        if os.path.exists(self.vc_file):
            with open(self.vc_file,"r") as f:
                fline = f.readlines()
                f.close()
                for line in fline:
                    if re.findall("    unit-cell volume",line):
                        origin_volume = line.split()[-2]
                        break
                for line in fline[int(len(fline) *0.5) : ]:
                    if re.findall("Begin final coordinates",line):
                        status = True
                        full = color("End").bold
                    if re.findall("NOT achieved", line):
                        status = True
                        full = color("Not Converged").warning
                    if re.findall("%",line):
                        status = True
                        full = color("Crashed!").warning
                    if status == True:
                        break

                    # if re.findall("",line):



                    if re.findall("ethr",line):
                        tmp_ethr = line.split()[2].rstrip(",")

                    if re.findall("P=",line):
                        tmp_pressure = line.split()[-1].lstrip("P=")
                    if re.findall("new unit-cell volume",line) and len(line.split()) == 6:
                        new_volume = line.split()[-2]
            if status == False:
                if re.match(regex, tmp_ethr):
                    tmp_ethr_p = (int(tmp_ethr.split("-")[-1])) * "▶"
                else:
                    tmp_ethr_p = ""
                full = color("RUNNING").green + " Now P = " +color(tmp_pressure).header + "\tethr = " + color(tmp_ethr).blue +"\t%s"%tmp_ethr_p
            # if re.findall("End",full):
                # full += "\tdV = " + \
                    # color(round((float(new_volume)-float(origin_volume))
                                # float(origin_volume) * 100, 5)).blue + "%"
        return info_build(prefix,full)

    @property
    def check_scf(self):
        full = ""
        prefix = "scf"
        cell_volume = ""
        status = False

        if os.path.exists(self.scf_file):
            with open(self.scf_file,"r") as f:
                for line in f:
                    if re.findall("unit-cell volume",line):
                        cell_volume = color( line.split()[3]  ) .underline + " au^3"
                    if re.findall("JOB DONE",line):
                        status = True
                    if re.findall("%%%%",line):
                        full = color("CRASHED").warning
            if status == True and full =="":
                full = make_tab([cell_volume,color("End").bold])
            else:
                full = make_tab([cell_volume, color("RUNNING").warning])
        return info_build(prefix, full)

    @property
    def check_ph(self):
        prefix = "ph"
        full = ""
        dynmat = None
        switch = -1
        filesize = 0

        if os.path.exists(self.ph_folder):
            all_file_list = os.listdir(self.ph_folder)
            for file_name in all_file_list:
                if file_name == "CRASH":
                    full = color("CRASHED!").warning
                    break

                if file_name.find("dyn") != -1 and file_name.find("matdyn") == -1:
                    dyn_num = int(file_name.replace("dyn", ""))
                    if dyn_num > switch:
                        switch = dyn_num
                        filesize = os.path.getsize(
                            os.path.join(self.ph_folder, file_name))
            if full == "":
                if switch == -1:
                    full = color("None").warning
                else:
                    if filesize == 0:
                        full = "dyn %s" % color("%2.0f"%switch).blue
                        tmp_folder = os.path.join(
                            self.ph_folder, "tmp/_ph0/hcpFe.phsave")
                        tmp_list = os.listdir(tmp_folder)
                        dynmat = 0
                        for tmp in tmp_list:
                            if re.findall("dynmat",tmp):
                                dynmat_0 = float(tmp.lstrip("dynmat.").rstrip(".xml"))
                                if dynmat_0 > dynmat:
                                    dynmat = dynmat_0
                        dynmat_out = color(dynmat).header
                        full += info_build("\tdynmat",dynmat_out)
                    elif switch != 12:
                        full = color("Processing").header
                    else:
                        full = color("End").bold
        return info_build(prefix, full)

    @property
    def check_if_matdyn(self):
        status = False
        if os.path.exists(self.ph_folder):
            all_file_list = os.listdir(self.ph_folder)
            if "matdyn.out" in all_file_list or "disp.out" in all_file_list:
                status = True
        return status


def print_T(T):
    print("\nT = %s"%color(T).header)


def make_tab(things:List) -> str:
    output = ""
    for i in things:
        if i == "":
            pass
        output += "%s\t"%i
    return output.rstrip("\t")


def do_check(T,P):
    check = check_folder(T, P)
    if check.check_if_matdyn:
        info = make_tab([check.folder_report,color("ALL DONE").blue])
        print(info)
    else: 
        info = make_tab([check.folder_report,check.check_vc,check.check_scf,check.check_ph])
        if len(info) > len(check.folder_report):
            print(info)

def do_print_ethr(T,P):
    check = check_folder(T, P)
    if check.check_if_matdyn:
        info = make_tab([check.folder_report, color("ALL DONE").blue])
        print(info)
    else:
        info = check.ethr
        print(info)


def check_T(T):
    print_T(T)
    for P in P_range:
        do_check(T,P)
    return None

def check_TP(T,P):
    print_T(T)
    do_check(T,P)
    do_print_ethr(T,P)
    return None

def start():
    for T in T_range:
        print_T(T)
        for P in P_range:
            do_check(T,P)
    return None


if __name__ == "__main__":
    if len(sys.argv) == 2:
        check_T(sys.argv[1])
    elif len(sys.argv) == 3:
        check_TP(sys.argv[1], sys.argv[2])
    else:
        start()

    print("\n")
