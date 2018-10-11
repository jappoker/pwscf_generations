#! /usr/bin/env python3

import os
import re
import shelve
import subprocess
import time
from math import exp

import numpy as np
from scipy.optimize import curve_fit

import functions

# Reading input

erro = 'CRASH'

working_dir = os.getcwd()

print('Starting workflow calculations.\n')

# Read input files

# Input.dat has information on the calculation:
# The pressures to be calculated, the initial guess for V0, K0 and Kp
# and the lattice vectors, necessary to calculate the lattice parameter with the volume

with open('input.dat', 'r') as input_params:
    lines = input_params.readlines()

# job_param is a file with the parameters of the job to be submitted:
# walltime
# number of nodes
# number of processors

with open('job_header', 'r') as job_in:
    job_lines = job_in.readlines()
    for i in range(0, len(job_lines)):
        if re.findall('^Scheduler', job_lines[i]):
            schedul = job_lines[i].split(':')[1].rstrip().lstrip()
        if re.findall('^Necessary modules to be', job_lines[i]):
            modul = job_lines[i].split(':')
        if re.findall('Number of nodes', job_lines[i]):
            nodes = int(job_lines[i].split(':')[1].rstrip().lstrip())
        if re.findall('Number of processor', job_lines[i]):
            procs = int(job_lines[i].split(':')[1].rstrip().lstrip())


press = np.array(lines[0].split(), dtype=float)
eos_parameter = np.array(lines[1].split(), dtype=float)
param_aux = [float(x) for x in lines[1].split()]
npress = len(press)


if npress < 7:
    print('You have %d pressure points. At least 7 are necessary for a good quality EoS fitting, please, modify your input files.' % npress)
else:
    print('You have %d pressure points, more than 7! You are a smart user!' % npress)

print('Initial equation of state parameters: V0 = %5.4f   K0 = %5.2f    Kp = %4.2f\n' %
      (eos_parameter[0], eos_parameter[1], eos_parameter[2]))

# lattice vectors a1, a2, a3
a1 = np.array(lines[2].split(), dtype=float)
a2 = np.array(lines[3].split(), dtype=float)
a3 = np.array(lines[4].split(), dtype=float)

vecs = np.array([a1, a2, a3])

# a = (V/(a3.a1^a2))^(1/3)

vol_sc = abs(np.dot(a3, np.cross(a1, a2)))

print('You have %s pressures and %d processors' %
      (str(npress), (procs * nodes)))
div = (nodes * procs) / npress
if (nodes * procs) % npress == 0:
    print('Therefore, I will consider %d procs per pressure.' % div)
else:
    print('Number of processors not divided by number of pressures. Exiting')
    quit()

# Calculate volumes for each pressure

signal_cg = True

dirs_cg = os.listdir(working_dir)  # Current file folder

num_cg = 0
cg_found_list = []
cont_cg = 0

for j in press:
    if os.path.exists('cg_' + str(j)):
        # Checks if cg_ folders exist.
        print('Folders cg for P = %s found.' % str(j))
        cg_found_list = cg_found_list + [cont_cg]
        num_cg = num_cg + 1
    cont_cg = cont_cg + 1  # Looping counter

if num_cg == npress:
    print('All cg calculations exist. Proceeding to vc-relax calculation.')
    signal_cg = False
elif num_cg == 0:
    print('Cg calculations not found. Proceeding to crude guess calculation.')
    signal_cg = True
    press_cg = np.copy(press)
else:
    print('Some remaining pressures to be calculated.')
    press_cg = np.delete(press, cg_found_list)
    signal_cg = True

while signal_cg:  # If cg_ exist, skip this part and goes to vc-relax optimization
    functions.create_files(press_cg, eos_parameter, vecs, vol_sc, 'cg_', 'scf')

    functions.create_job(job_lines, press_cg, div,
                         'job-cg.sh', 'cg_', 'pw', modul)

    job_sub = subprocess.Popen(['sbatch', 'job-cg.sh'], stdout=subprocess.PIPE)
    output = job_sub.stdout.read().split()[-1]
    job_id = int(output)

    print('Job submitted. The Job ID is %d\n' % job_id)
    print('Waiting calculation. This can take a while, you can grab a coffee!')

    signal_job_cg = True

#	This part monitors the job execution.
    while signal_job_cg:
        # request SLURM queue information
        job_status = subprocess.Popen(
            ['squeue', '-j', str(job_id)], stdout=subprocess.PIPE)
        outqueue = job_status.stdout.read().split()
        if output not in outqueue:  # output is a string that contains the Job ID. If it is not in the queue information, the job is done
            print(
                "Crude guess is done! Continuing calculation... I hope your coffee was good!")
            signal_job_cg = False
        else:  # If it is, sleeps a litle bit and request queue information again.
            time.sleep(2)
            
    signal_cg = False


Paux = []
Vaux = []
cont = 0
index_list = []
Number_error_files = 0
for i in press:
    if erro in os.listdir('cg_' + str(i)):
        # If CRASH files exist, something went wrong. Stop the workflow
        print('Something went wrong, CRASH files found. Check your cg outputs.')
        quit()
    else:
        file_path = 'cg_' + str(i) + '/' + 'cg_' + str(i) + '.out'
        try:
            with open(file_path) as qe_out:
                qe_lines = qe_out.readlines()
        except FileNotFoundError:  # Check if output files exist.
            print('The output file for P = %d was not found. Removing it from list and continuing calculation.' % i)
            index_list = index_list + [cont]
            # Counts number of files that wasn't found
            Number_error_files = Number_error_files + 1
            pass
        else:
            for j in range(0, len(qe_lines)):
                if re.findall('P=', qe_lines[j]):
                    Paux = Paux + \
                        [float(qe_lines[j].split('=')[-1].rstrip().lstrip()) / 10]
                if re.findall('volume', qe_lines[j]):
                    Vaux = Vaux + \
                        [float(qe_lines[j].split()[-2].rstrip().lstrip())]
    cont = cont + 1

if Number_error_files > int(npress / 2):
    # If more that a half of the calculations are not found, stop the workflow
    print('More than a half of pressures were not calculated. Exiting')
    quit()

# deletes the pressures that weren't found
press = np.delete(press, index_list)

# If the number of volumes are not equal to the number of pressures, something is wrong, stop the workflow
if len(press) != len(Paux) or len(Paux) != len(Vaux):
    print("len of press = %s, len of Paux = %s"%(len(press),len(Paux)))
    print('ERROR: Some pressures were not found in the files. Please, review your data.')
    quit()

print('The initial volumes and pressures were written in the file Initial_PxV.dat\n')

pv = open('Initial_PxV.dat', 'w')
pv.write("Initial PxV\n")
pv.write("P (GPa)     V (au^3)\n")

for i in range(0, len(Paux)):
    pv.write("%7.2f   %7.2f\n" % (Paux[i], Vaux[i]))

P = np.array(Paux)
V = np.array(Vaux)
param1 = np.array(param_aux)

# Fits the Equation of State with the input parameters as a guess
eosPar, eosCov = curve_fit(functions.vinet, V, P, param1)
std = np.sqrt(np.diag(eosCov))

##### Chi^2 teste ######
chi2 = 0
for i in range(0, len(Paux)):
    chi2 = chi2 + (Paux[i] - functions.vinet(Vaux[i],
                                             eosPar[0], eosPar[1], eosPar[2]))**2

chi = np.sqrt(chi2)

print("EoS fitted, these are the standart deviations:\n")
print("V0 = %6.3f" % std[0])
print("K0 = %6.3f" % std[1])
print("Kp = %6.3f\n" % std[2])

print("chi = %7.4f" % chi)
pv.write("Results for a Vinet EoS fitting:\n")
pv.write("V0 = %6.4f    K0 = %4.2f    Kp = %4.2f\n" %
         (eosPar[0], eosPar[1], eosPar[2]))
pv.close()

# Generate new files for a vc-relax job

signal_vc = True

dirs_vc = os.listdir(working_dir)

# The routines below are the same as thosr for cg_ calculation

num_vc = 0
cont_vc = 0
vc_found_list = []

for j in press:
    if os.path.exists('vc_' + str(j)):
        # Checks if cg_ folders exist.
        print('Folders vc for P = %s found.' % str(j))
        vc_found_list = vc_found_list + [cont_vc]
        num_vc = num_vc + 1
    cont_vc = cont_vc + 1

if num_vc == npress:
    print('All vc calculations exist. Proceeding to EoS fitting.')
    signal_vc = False

elif num_vc == 0:
    print('Vc calculations not found. Proceeding to vc-relax calculation.')
    signal_vc = True
    press_vc = np.copy(press)
else:
    print('Some remaining pressures to be calculated.')
    press_vc = np.delete(press, vc_found_list)
    signal_vc = True

while signal_vc:
    print('Now I will generate the files for a vc-relax calculation')

    time.sleep(1)

    functions.create_files(press_vc, eosPar, vecs, vol_sc, 'vc_', 'vc-relax')

    functions.create_job(job_lines, press_vc, div,
                         'job-vc.sh', 'vc_', 'pw', modul)

    job_sub_vc = subprocess.Popen(
        ['sbatch', 'job-vc.sh'], stdout=subprocess.PIPE)
    output_vc = job_sub_vc.stdout.read().split()[-1]
    job_id_vc = int(output_vc)

    print('Job for variable cell relaxation submitted. The Job ID is %d\n' % job_id_vc)
    print('Waiting calculation. This can take a long time, you can have lunch! Or take a nap! :)')

    signal_job_vc = True

    while signal_job_vc:
        job_status_vc = subprocess.Popen(
            ['squeue', '-j', str(job_id_vc)], stdout=subprocess.PIPE)
        outqueue_vc = job_status_vc.stdout.read().split()
        if output_vc not in outqueue_vc:
            print("Congratulations, the variable cell relaxation is done!")
            print("Fitting the final Equation of State.")
            signal_job_vc = False
        else:
            time.sleep(2)
    signal_vc = False


eosInitPar = eosPar

del Paux[:]
del Vaux[:]

prov = []
Paux = []
Vaux = []
del index_list[:]
cont = 0
for i in press:
    if erro in os.listdir('vc_' + str(i)):
        print('CRASH file found. Something went wrong with the vc-relax calcultion, check your files.')
        quit()
    else:
        file_path = 'vc_' + str(i) + '/' + 'vc_' + str(i) + '.out'
        try:
            with open(file_path) as vc_out:
                vc_lines = vc_out.readlines()
        except FileNotFoundError:
            print(
                'Output file for pressure %d not found. Removing it from pressure list.' % i)
            index_list = index_list + [cont]
            pass
        else:
            for j in range(0, len(vc_lines)):
                if re.findall('P=', vc_lines[j]):
                    prov = prov + \
                        [float(vc_lines[j].split('=')[-1].rstrip().lstrip()) / 10]
                if re.findall('final unit-cell volume', vc_lines[j]):
                    Vaux = Vaux + \
                        [float(vc_lines[j].split()[-2].rstrip().lstrip())]
            Paux = Paux + [prov[-1]]
    cont = cont + 1

press = np.delete(press, index_list)

out = open('Out_file.dat', 'w')
out.write('P (GPa)    V (au^3)\n')
if len(Paux) != len(Vaux):
    print('Some volumes are mossing, check your vc-relax outputs.')
    quit()

for i in range(0, len(Paux)):
    out.write("%7.2f   %7.2f\n" % (Paux[i], Vaux[i]))

P1 = np.array(Paux)
V1 = np.array(Vaux)

eosFinalPar, eosFinalCov = curve_fit(functions.vinet, V1, P1, eosInitPar)
std = np.sqrt(np.diag(eosFinalCov))

chi2 = 0
for i in range(0, len(Paux)):
    chi2 = chi2 + (Paux[i] - functions.vinet(Vaux[i],
                                             eosFinalPar[0], eosFinalPar[1], eosFinalPar[2]))**2

chi = np.sqrt(chi2)
print("Ready, the fitting is done. The standart deviations are:\n")
print("V0 = %6.3f" % std[0])
print("K0 = %6.3f" % std[1])
print("Kp = %6.3f\n" % std[2])

print("chi = %7.4f\n" % chi)
print("The final EoS is written in the file Out_file.dat\n")
print("I hope you have a great day! See you next time.")

out.write("Final Results for a Vinet EoS fitting:\n")
out.write("V0 = %6.4f    K0 = %4.2f    Kp = %4.2f\n" %
          (eosFinalPar[0], eosFinalPar[1], eosFinalPar[2]))
out.close()


# Saving the calculation parameters for external use. (Eg. phonon calculation)
calc_par = shelve.open('calculation_parameters')
calc_par['pressures'] = press
