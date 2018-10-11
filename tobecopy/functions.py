import numpy as np
import re
import subprocess
import os
import time
from random import randint
import shutil
import atommass as at
from func_cij import Ekl, apply_strain

def create_qe_input(alat, lat_vec, calc_type, name, pres, tmp_folder, atom_pos=''):
	""" Function to create an input file for quantum espresso pw.x"""
	a = float(alat)
	fp = open(name,'w')
	with open('qe_input_data') as qe_in:
		lines = qe_in.readlines()
		for i in range(0,len(lines)): # Find the number of the line whe atomic positions begin
			if re.findall('^atomic positions',lines[i]):
				pos_at = i
			elif re.findall('^pseudopotential',lines[i]):
				pseudo = lines[i].split(':')[1].rstrip().lstrip()
				np = i + 1 # line where pseudo begins
			elif re.findall('^prefix',lines[i]):
				pre = lines[i].split(':')[1].rstrip().lstrip()
			elif re.findall('^number of atoms',lines[i]):
				nat = lines[i].split(':')[1].rstrip().lstrip()
			elif re.findall('^number of atomic',lines[i]):
				ntyp = lines[i].split(':')[1].rstrip().lstrip()
			elif re.findall('^wmass',lines[i]):
				wmass = lines[i].split(':')[1].rstrip().lstrip()
			elif re.findall('^energy cutoff', lines[i]):
				ecut = lines[i].split(':')[1].rstrip().lstrip()
			elif re.findall('^k-points', lines[i]):
				kp = lines[i].split(':')[1].rstrip().lstrip().split()
			elif re.findall('^Atomic positions', lines[i]):
				na = i + 1
			elif re.findall('^smearing',lines[i]):
				sm = lines[i].split(':')[1].rstrip().lstrip()
			elif re.findall('^shift', lines[i]):
                                sp = lines[i].split(':')[1].rstrip().lstrip().split()
		if (sm == 'gaussian' or sm == 'mp' or sm == 'fd'):
			for i in range(0,len(lines)):
				if re.findall('^occupations', lines[i]):
					occ = lines[i].split(':')[1].rstrip().lstrip()
				elif re.findall('^degauss', lines[i]):
					deg = lines[i].split(':')[1].rstrip().lstrip()
	fp.write("&control\n")
	fp.write("calculation = '%s',\n" % calc_type)
	fp.write("restart_mode = 'from_scratch',\n")
	fp.write("tstress = .true.\n")
	fp.write("tprnfor = .true.\n")
	fp.write("prefix = '%s'\n"  % pre)
	fp.write("pseudo_dir = '%s'\n" % pseudo)
	fp.write("outdir = './%s'\n" % tmp_folder)
	fp.write("nstep = 400\n")
	fp.write("dt = 15\n")
	fp.write("etot_conv_thr = 1.0d-7\n")
	fp.write("forc_conv_thr = 1.0d-6\n")
	fp.write("/\n")
	fp.write("&system\n")
	fp.write("ibrav = 0\n")
	fp.write("celldm(1) = %9.5f\n" % a)
	fp.write("nat = %d\n" % int(nat))
	fp.write("ntyp = %d\n" % int(ntyp))
	fp.write("ecutwfc = %4.1f\n" % float(ecut))
	if (sm == 'gaussian' or sm == 'mp' or sm == 'fd'):
		fp.write("occupations='%s', smearing='%s', degauss= %7.9f,\n" % (occ, sm, float(deg)))
	fp.write("nr1 = 48\n")
	fp.write("nr2 = 48\n")
	fp.write("nr3 = 64\n")
	fp.write("/\n")
	fp.write("&electrons\n")
	fp.write("diagonalization='david'\n")
	fp.write("mixing_mode = 'plain'\n")
	fp.write("mixing_beta = 0.3\n")
	fp.write("conv_thr = 1.0d-12\n")
	fp.write("/\n")
	fp.write("&ions\n")
	fp.write("ion_dynamics='damp'\n")
	fp.write("/\n")
	fp.write("&cell\n")
	fp.write("cell_dynamics = 'damp-w'\n")
	fp.write("press = %d\n" % int(pres*10))
	fp.write("wmass = %7.5f\n" % float(wmass))
	fp.write("/\n")
	fp.write("CELL_PARAMETERS\n")
	fp.write("%10.7f   %10.7f   %10.7f\n" % (lat_vec[0][0], lat_vec[0][1], lat_vec[0][2]))
	fp.write("%10.7f   %10.7f   %10.7f\n" % (lat_vec[1][0], lat_vec[1][1], lat_vec[1][2]))
	fp.write("%10.7f   %10.7f   %10.7f\n" % (lat_vec[2][0], lat_vec[2][1], lat_vec[2][2]))
	fp.write("ATOMIC_SPECIES\n")
	for i in range(0,int(ntyp)): # for each o the atomic types, sets th atom mass as 1.0 and writes the name of the pseudopotential file
		atom = re.split('[. -]',lines[np+i])[0] # Gets the atom
		atom = "Fe"
		fp.write("%3s   1.0   %s" % (atom, lines[np+i]))
	fp.write("ATOMIC_POSITIONS (crystal)\n")
	if atom_pos:
		for i in range(0,int(nat)):
			fp.write(atom_pos[i])
	else:
		for i in range(0,int(nat)):
			fp.write(lines[na + i])
	fp.write("K_POINTS (Automatic)\n")
	fp.write("%2d  %2d  %2d  %2d  %2d  %2d" % (int(kp[0]), int(kp[1]), int(kp[2]), int(sp[0]), int(sp[1]), int(sp[2])) )
	fp.close()


def create_job(job_head, pressures, d, job_name, dir_name, prog, modules=''):
	"""Function to create the job to be submitted"""
	jb=open(job_name,'w')
	for i in range(4, len(job_head)):
		jb.write('%s\n' % job_head[i])
	if modules:
		for i in range(1, len(modules)):
			jb.write('%s\n' % modules[i])
	jb.write('for i in ')
	for i in pressures:
		jb.write('%s ' % str(i))

	jb.write('; do\n')
	jb.write('cd %s$i\n' % dir_name)
	if prog == 'pw':
		jb.write('mpirun -np %d pw.x -in %s$i.in > %s$i.out &\n' % (d,dir_name,dir_name))
	if prog =='ph':
		jb.write('mpirun -np %d ph.x -in %s$i-ph.in > %s$i-ph.out &\n' % (d,dir_name,dir_name))
	jb.write('sleep 3\n')
	jb.write('cd ..\n')
	jb.write('done\n')
	jb.write('wait\n')
	jb.close()

def create_job_cij(job_head, pressures, d, job_name, dir_name, cclass, nstrain, dist, modules=''):
	strain_rel = relates_strain(cclass)
	jb=open(job_name,'w')
	for i in range(4, len(job_head)):
		jb.write('%s\n' % job_head[i])
	if modules:
		for i in range(1, len(modules)):
			jb.write('%s\n' % modules[i])
	jb.write('for i in ')
	for i in pressures:
		jb.write('%s ' % str(i))
	jb.write('; do\n')
	jb.write('cd %s$i/cij\n' % dir_name)
	jb.write('for j in ')
	for i in range(1, nstrain + 1):
		jb.write('e%d ' % strain_rel[str(i)])
	jb.write('; do\n')
	jb.write('cd cij_strain-$j\n')
	jb.write('for k in ')
	for i in range(1,len(dist)):
		jb.write('%s ' % dist[i])
	jb.write('; do\n')
	jb.write('mpirun -np %d pw.x -in cij_strain-$j$k.in > cij_strain-$j$k.out &\n' % d)
	jb.write('sleep 3\n')
	jb.write('done\n')
	jb.write('cd ..\n')
	jb.write('done\n')
	jb.write('cd ../..\n')
	jb.write('done\n')
	jb.write('wait\n')

#def create_job(nd, prc, d, time, pressures, job_name, dir_name, prog='pw'):
#	"""Funtion to create a job to be submit in a Slurm scheduler"""
## Construct qe input
#	jb = open(job_name,'w')
#	jb.write('#!/bin/sh\n')
#	jb.write('#SBATCH -A mphys\n')
#	jb.write('#SBATCH -N %d\n' % nd)
#	jb.write('#SBATCH --tasks-per-node=%d\n' % prc)
#	jb.write('#SBATCH -J job\n')
#	jb.write('#SBATCH --time=%s:%s:%s\n\n' % (time[0], time[1], time[2]))
#	jb.write('module load intel-parallel-studio/2017\n\n')
#	jb.write('for i in ')
#	for i in pressures:
#		jb.write('%s ' % str(i))
#
#	jb.write('; do\n')
#	jb.write('cd %s$i\n' % dir_name)
#	if prog == 'pw':
#		jb.write('mpirun -np %d /rigel/home/mld2189/bin/programs/qe-6.1/bin/pw.x -in %s$i.in > %s$i.out &\n' % (d,dir_name,dir_name))
#	if prog =='ph':
#		jb.write('mpirun -np %d /rigel/home/mld2189/bin/programs/qe-6.1/bin/ph.x -in %s$i-ph.in > %s$i-ph.out &\n' % (d,dir_name,dir_name))
##	jb.write('PW $SLURM_JOBID %d %s$i.in 2 &\n' % (d,dir_name))
#	jb.write('sleep 3\n')
#	jb.write('cd ..\n')
#	jb.write('done\n')
#	jb.write('wait\n')
#	jb.close()

#def create_job_cij(nd, prc, d, time, pressures, job_name, dir_name, cclass, nstrain, dist):
#	strain_rel = relates_strain(cclass)
#	jb = open(job_name,'w')
#	jb.write('#!/bin/sh\n')
#	jb.write('#SBATCH -A mphys\n')
#	jb.write('#SBATCH -N %d\n' % nd)
#	jb.write('#SBATCH --tasks-per-node=%d\n' % prc)
#	jb.write('#SBATCH -J job\n')
#	jb.write('#SBATCH --time=%s:%s:%s\n\n' % (time[0], time[1], time[2]))
#	jb.write('module load intel-parallel-studio/2017\n\n')
#	jb.write('for i in ')
#	for i in pressures:
#		jb.write('%s ' % str(i))
#	jb.write('; do\n')
#	jb.write('cd %s$i/cij\n' % dir_name)
#	jb.write('for j in ')
#	for i in range(1, nstrain + 1):
#		jb.write('e%d ' % strain_rel[str(i)])
#	jb.write('; do\n')
#	jb.write('cd cij_strain-$j\n')
#	jb.write('for k in ')
#	for i in range(1,len(dist)):
#		jb.write('%s ' % dist[i])
#	jb.write('; do\n')
#	jb.write('mpirun -np %d /rigel/home/mld2189/bin/programs/qe-6.1/bin/pw.x -in cij_strain-$j$k.in > cij_strain-$j$k.out\n' % d)
#	jb.write('sleep 3\n')
#	jb.write('done\n')
#	jb.write('cd ..\n')
#	jb.write('done\n')
#	jb.write('cd ../..\n')
#	jb.write('done\n')
#	jb.write('wait\n')

def create_files(press, params_eos, vectors, volsc, file_n, calc):
	""" This function creates the files at each pressure"""
	with open('qe_input_data') as qe_in:
		lines = qe_in.readlines()
		for i in range(0,len(lines)): # Find the number of the line with the scratch directory
			if re.findall('^scratch',lines[i]):
				scratch_folder = lines[i].split(':')[1].rstrip().lstrip()
	for i in press:
		data = newton(i, params_eos[0], 100, params_eos)
		alat = (data/volsc)**(1/3)
		file_name = file_n + str(i)
		create_qe_input(alat, vectors, calc, file_name + '.in', i, 'tmp')
		os.makedirs(file_name)
		shutil.move(file_name + '.in', file_name)
		create_tmp_folder(scratch_folder, file_name, 'tmp')
#		subprocess.call(["mkdir", file_name])
#		subprocess.call(["mv", file_name + '.in', file_name])

def create_scf_files(vectors, alat, dir_scf, press, atom_pos):
	"""Creates scf files after a full vc-relax calculation"""
	with open('qe_input_data') as qe_in:
		lines = qe_in.readlines()
		for i in range(0,len(lines)):
			if re.findall('^scratch', lines[i]):
				scratch_folder = lines[i].split(':')[1].rstrip().lstrip()
	file_scf = dir_scf + '.in'
	os.makedirs(dir_scf)
	create_qe_input(alat, vectors, 'scf', file_scf, press, 'tmp', atom_pos)
	shutil.move(file_scf, dir_scf)
	create_tmp_folder(scratch_folder, dir_scf,'tmp')

def create_cij_files(vectors, alat, dir_scf, dir_cij, press, atom_pos, dist, nstrain,cclass):
	""" Creates the distored structures and put each file in its place"""
	strain_rel = relates_strain(cclass)
	with open('qe_input_data') as qe_in:
		lines = qe_in.readlines()
		for i in range(0,len(lines)):
			if re.findall('^scratch', lines[i]):
				scratch_folder = lines[i].split(':')[1].rstrip().lstrip()
	os.makedirs(dir_scf + '/' + dir_cij)
	for i in range(1, nstrain + 1):
		strain_number  = strain_rel[str(i)]
		os.makedirs(dir_scf + '/' + dir_cij + '/' + dir_cij + '_strain-e' + str(strain_number))
		file_cij = dir_cij + '_strain-e' + str(strain_number)
		for j in range(1, len(dist)):
			tmp_folder = 'tmp' + dist[j]
			ap_st = apply_strain()
			new_vet = ap_st.e(vectors, float(dist[j]), strain_number)
			create_qe_input(alat, new_vet, 'relax', file_cij + dist[j] + '.in', press, tmp_folder, atom_pos)
			shutil.move(file_cij + dist[j] + '.in', dir_scf + '/' + dir_cij + '/' + dir_cij + '_strain-e' + str(strain_number))
			create_tmp_folder(scratch_folder,  dir_scf + '/' + dir_cij + '/' + dir_cij + '_strain-e' + str(strain_number), tmp_folder)

#		if cclass == '1' or cclass == '-1' or cclass == '2' or cclass == 'm' or cclass == '2/m' or cclass == '222' or cclass == 'mm2' or cclass == '2mm' or cclass == 'mm' or cclass == 'mmm':

def calculate_cij_tensor(Cij, dir_cij, cclass, dist, j):
	strain_rel = relates_strain(cclass)
	strain_number = strain_rel[str(j)]
	index = strain_number - 1
	dir_strain = dir_cij + '/' + 'cij_strain-e' + str(strain_number)
	stress_p_filename = dir_strain + '/' + 'cij_strain-e' + str(strain_number) + str(dist[2]) + '.out'
	stress_n_filename = dir_strain + '/' + 'cij_strain-e' + str(strain_number) + str(dist[1]) + '.out'
	stress_p = read_stress(stress_p_filename)
	stress_n = read_stress(stress_n_filename)
	if strain_number == 1:
#		C11
		Cij[index, 0] = calc_cij(stress_p, stress_n, 0, dist[1])
#		C12
		Cij[index, 1] = calc_cij(stress_p, stress_n, 1, dist[1])
		if cclass == '-62m':
			Cij[1, 1] = Cij[0,0]
#		C66 for some crystal classes can be calculated with C11 and C12
		if cclass == '3' or cclass == '-3' or cclass == '32' or cclass == '3m' or cclass == '-3m' or cclass == '-32/m' or cclass == '6' or cclass == '-6' or cclass == '622' or cclass == '62' or cclass == '6mm' or cclass == '6/m' or cclass == '6/mmm' or cclass == '-6m2':
			Cij[5,5] = 0.5*(Cij[0,0] - Cij[0,1])
			Cij[1, 1] = Cij[0,0]
#		For cubic systems, C13, C23, C22 and C33 can be calculated at this point
		if cclass == '23' or cclass == '-43m' or cclass == '432' or cclass =='43' or cclass == 'm3' or cclass == '2/m-3' or cclass == 'm3m' or cclass == 'm-3m':
			Cij[0,2] = Cij[0,1]
			Cij[1,2] = Cij[0,1]
			Cij[1,1] = Cij[0,0]
			Cij[2,2] = Cij[0,0]
#		C13
		else:
			Cij[index, 2] = calc_cij(stress_p, stress_n, 2, dist[1])
#		For tetragonal trigonal and hexagonal cells, C13 defines C23
		if cclass == '4' or cclass == '-4' or cclass == '4/m' or cclass == '4mm' or cclass == '422' or cclass == '42' or cclass == '-42m' or cclass == '4/mmm' or cclass == '3' or cclass == '-3' or cclass == '32' or cclass == '3m' or cclass == '-3m' or cclass == '-32/m' or cclass == '6' or cclass == '-6' or cclass == '622' or cclass == '62' or cclass == '6mm' or cclass == '6/m' or cclass == '6/mmm':
			Cij[1,2] = Cij[0,2]
#		C14
		if cclass == '1' or cclass == '-1' or cclass == '3' or cclass == '-3':
			Cij[index, 3] = calc_cij(stress_p, stress_n, 3, dist[1])
			if cclass == '3' or cclass == '-3':
				Cij[1,3] = -Cij[0,3]
				Cij[4,5] = 2*Cij[0,3]
		elif cclass == '32' or cclass == '3m' or cclass == '-3m' or cclass == '-32/m':
			Cij[index, 3] = -calc_cij(stress_p, stress_n, 3, dist[1])
			Cij[1,3] = -Cij[0,3]
			Cij[4,5] = 2*Cij[0,3]
		else:
			Cij[index, 3] = 0.0
#		C15
		if cclass == '1' or cclass == '-1' or cclass == '2' or cclass == 'm' or cclass == '2/m':
			Cij[index, 4] = calc_cij(stress_p, stress_n, 4, dist[1])
		elif cclass == '3' or cclass == '-3':
			Cij[index, 4] = -calc_cij(stress_p, stress_n, 4, dist[1])
			Cij[1,4] = -Cij[0,4]
			Cij[3,5] = 2*Cij[1,4]
		else:
			Cij[index, 4] = 0.0
#		C16
		if cclass == '1' or cclass == '-1' or cclass == '4' or cclass == '-4' or cclass == '4/m':
			Cij[index, 5] = calc_cij(stress_p, stress_n, 5, dist[1])
			if cclass == '4' or cclass == '-4' or cclass == '4/m':
				Cij[1,5] = Cij[0,5]
		else:
			Cij[index, 5] = 0.0

	if strain_number == 2:
#		C22 and C23
		if cclass == '1' or cclass == '-1' or cclass == '2' or cclass == 'm' or cclass == '2/m' or cclass == '222' or cclass == 'mm2' or cclass == '2mm' or cclass == 'mm' or cclass == 'mmm':
			Cij[index, 1] = calc_cij(stress_p, stress_n, 1, dist[1])
			Cij[index, 2] = calc_cij(stress_p, stress_n, 2, dist[1])
#		C24
		if cclass == '1' or cclass == '-1':
			Cij[index, 3] = calc_cij(stress_p, stress_n, 3, dist[1])
#		C25
		if cclass == '1' or cclass == '-1' or cclass == '2' or cclass == 'm' or cclass == '2/m':
			Cij[index, 4] = calc_cij(stress_p, stress_n, 4, dist[1])
#		C26
		if cclass == '1' or cclass == '-1':
			Cij[index, 5] = calc_cij(stress_p, stress_n, 5, dist[1])
	if strain_number == 3:
		Cij[index, 2] = calc_cij(stress_p, stress_n, 2, dist[1])
		if cclass == '-6m2':
			Cij[1,2] = calc_cij(stress_p, stress_n, 2, dist[1])
#		C34
		if cclass == '1' or cclass == '-1':
			Cij[index, 3] = calc_cij(stress_p, stress_n, 3, dist[1])
			Cij[index, 5] = calc_cij(stress_p, stress_n, 5, dist[1])
#		C35
		if cclass == '1' or cclass == '-1' or cclass == '2' or cclass == 'm' or cclass == '2/m':
			Cij[index, 4] = calc_cij(stress_p, stress_n, 4, dist[1])
	if strain_number == 4:
		Cij[index, 3] = 0.5*calc_cij(stress_p, stress_n, 3, dist[1])
		if cclass == '4' or cclass == '-4' or cclass == '4/m' or cclass == '4mm' or cclass == '422' or cclass == '42' or cclass == '-42m' or cclass == '4/mmm' or cclass == '3' or cclass == '-3' or cclass == '32' or cclass == '3m' or cclass == '-3m' or cclass == '-32/m' or cclass == '6' or cclass == '-6' or cclass == '622' or cclass == '62' or cclass == '6mm' or cclass == '6/m' or cclass == '6/mmm':
			Cij[4,4] = Cij[index,3]
		if cclass == '23' or cclass == '-43m' or cclass == '432' or cclass =='43' or cclass == 'm3' or cclass == '2/m-3' or cclass == 'm3m' or cclass == 'm-3m':
			Cij[4,4] = Cij[index,3]
			Cij[5,5] = Cij[index,3]
#		C45
		if cclass == '1' or cclass == '-1':
			Cij[index, 4] = 0.5*calc_cij(stress_p, stress_n, 4, dist[1])
		if cclass == '1' or cclass == '-1' or cclass == '2' or cclass == 'm' or cclass == '2/m':
			Cij[index, 5] = 0.5*calc_cij(stress_p, stress_n, 5, dist[1])
	if strain_number == 5:
#		C55
		if cclass == '1' or cclass == '-1' or cclass == '2' or cclass == 'm' or cclass == '2/m' or cclass == '222' or cclass == 'mm2' or cclass == '2mm' or cclass == 'mm' or cclass == 'mmm' or cclass == '-6m2':
			Cij[index, 4] = 0.5*calc_cij(stress_p, stress_n, 4, dist[1])
		if cclass == '1' or cclass == '-1':
			Cij[index, 5] = 0.5*calc_cij(stress_p, stress_n, 5, dist[1])
	if strain_number == 6:
#		C66
		if cclass == '1' or cclass == '-1' or cclass == '2' or cclass == 'm' or cclass == '2/m' or cclass == '222' or cclass == 'mm2' or cclass == '2mm' or cclass == 'mm' or cclass == 'mmm' or cclass == '4' or cclass == '-4' or cclass == '4/m' or cclass == '4mm' or cclass == '422' or cclass == '42' or cclass == '-42m' or cclass == '4/mmm':
			Cij[index, 5] = 0.5*calc_cij(stress_p, stress_n, 5, dist[1])

#	Completing the tensor
	k = 1
	for i in range(0,6):
		for j in range(k,6):
			Cij[j,i] = Cij[i,j]
		k = k + 1


def relates_strain(cclass):
	"""Function that relates a number between 1 and 6 to the necessary strains to be applied to each class"""
	if cclass == '1' or cclass == '-1' or cclass == '2' or cclass == 'm' or cclass == '2/m' or cclass == '222' or cclass == 'mm2' or cclass == '2mm' or cclass == 'mm' or cclass == 'mmm':
		dic_strain = {'1' : 1,  '2' : 2, '3' : 3, '4' : 4, '5' : 5, '6' : 6}
	elif cclass == '4' or cclass == '-4' or cclass == '4/m' or cclass == '4mm' or cclass == '422' or cclass == '42' or cclass == '-42m' or cclass == '4/mmm':
		dic_strain = {'1' : 1, '2' : 3, '3' : 4, '4' : 6}
	elif cclass == '3' or cclass == '-3' or cclass == '32' or cclass == '3m' or cclass == '-3m' or cclass == '-32/m' or cclass == '6' or cclass == '-6' or cclass == '622' or cclass == '62' or cclass == '6mm' or cclass == '6/m' or cclass == '6/mmm':
		dic_strain = {'1' : 1, '2' : 3, '3' : 4}
	elif cclass == '-6m2':
		dic_strain = {'1' : 1, '2' : 3, '3' : 5}
	elif cclass == '23' or cclass == '-43m' or cclass == '432' or cclass =='43' or cclass == 'm3' or cclass == '2/m-3' or cclass == 'm3m' or cclass == 'm-3m':
		dic_strain = {'1' : 1, '2' : 4}
	return dic_strain

def create_ph_input(name):
	""" Creates input files for phonons calculations."""
	with open('qe_input_data') as qe_in:
		lines_qe = qe_in.readlines()
		for i in range(0,len(lines_qe)):
			if lines_qe[i].startswith('number of atomic'):
				ntyp = lines_qe[i].split(':')[1].rstrip().lstrip()
			if lines_qe[i].startswith('pseudopotential'):
				np = i + 1
			if lines_qe[i].startswith('prefix'):
				pre = lines_qe[i].split(':')[1].rstrip().lstrip()

	with open('ph_input_data') as ph_in:
		lines_ph = ph_in.readlines()
		for i in range(0,len(lines_ph)):
			if lines_ph[i].startswith('epsil'):
				ep =  lines_ph[i].split(':')[1].rstrip().lstrip()
			if lines_ph[i].startswith('ldisp'):
				lp =  lines_ph[i].split(':')[1].rstrip().lstrip()
			if lines_ph[i].startswith('nq1'):
				nq1 =  lines_ph[i].split(':')[1].rstrip().lstrip()
			if lines_ph[i].startswith('nq2'):
				nq2 =  lines_ph[i].split(':')[1].rstrip().lstrip()
			if lines_ph[i].startswith('nq3'):
				nq3 =  lines_ph[i].split(':')[1].rstrip().lstrip()

	fh = open(name,'w')
	fh.write('Phonon\n')
	fh.write('&inputph\n')
	fh.write('prefix = %s\n' % pre)
	fh.write('tr2_ph=1.0d-14\n')
	fh.write('ldisp = .true.\n')
	fh.write("verbosity = 'high'\n")
	fh.write('epsil = %s\n' % ep)
	fh.write("fildyn = 'dyn'\n")
	fh.write("recover = .true.\n")
	fh.write("outdir = './tmp'\n")
	for i in range(0, int(ntyp)):
		atom = re.split('[. -]', lines_qe[np + i])[0] # Reads the atom
		fh.write('amass(%d) = %7.5f\n' % (i + 1, at.atWeight[atom]))
	fh.write('nq1 = %d\n' % int(nq1))
	fh.write('nq2 = %d\n' % int(nq2))
	fh.write('nq3 = %d\n' % int(nq3))
	fh.write('/')

def vinet(V, V0, K0, Kp):
	""" Vinet EoS functions """
	f = (V/V0)**(1/3)
	P = 3*K0*((1 - f)/f**2)*np.exp(1.5*(Kp - 1)*(1 - f))
	return P
def dvinet(V,paramd):
	""" First derivative ofo Vinet EoS. Necessary to find the volumes with the Newton-Raphson method"""
	f = (V/paramd[0])**(1/3)
	dP = -3*paramd[1]*(1/f**2 + (2*(1 - f))/f**3 + (1.5*(1-f)*(paramd[2] - 1))/f**2)*np.exp(1.5*(paramd[2] - 1)*(1 - f))*(1/(3*paramd[0]*f**2))
	return dP

def newton(P,x0, niter, paramn):
	""" VERY primitive implementation of Newton-Raphson method. Assumes guess is close to final results e always performs 100 iterations. """
	""" Needs improvement """
	# Find a new initial guess closer to the desired pressure
#	m = -(paramn[1]/paramn[0])
#	guess = (P/m) + x0
	x = x0 - 5
	while vinet(x,paramn[0], paramn[1], paramn[2]) < P:
		x0 = x0 - 5
		x = x0
	guess  = 0.5*(x + x0)
	for i in range(1,niter):
		x = guess - (vinet(guess,paramn[0], paramn[1], paramn[2]) - P)/dvinet(guess,paramn)
		guess = x
	return x
def create_tmp_folder(scratch, folder_name, tmp_folder_name):
	""" creates temporary folder in the system's scratch directory """
	tm = time.strftime("%d.%m.%Y.%M") # System's time
	random_number = randint(1,30000) # Rando number between 1 and 30000
	scratch_dir = scratch + '/' + 'tmp.' + tm + str(random_number) # Scratch folder string
	os.makedirs(scratch_dir) # Creates scradtch folder
	os.symlink(scratch_dir, folder_name + '/' + tmp_folder_name) # Creates a link called tmp at the current directory that points to the scratch folder

def read_stress(file_name):
	try:
		with open(file_name) as inp:
			stress_lines = inp.readlines()
			for j in range(0,len(stress_lines)):
				if re.findall('total   stress', stress_lines[j]):
					sigma11 = float(stress_lines[j+1].split()[-3].rstrip().lstrip())
					sigma12 = float(stress_lines[j+1].split()[-2].rstrip().lstrip())
					sigma13 = float(stress_lines[j+1].split()[-1].rstrip().lstrip())
					sigma22 = float(stress_lines[j+2].split()[-2].rstrip().lstrip())
					sigma23 = float(stress_lines[j+2].split()[-1].rstrip().lstrip())
					sigma33 = float(stress_lines[j+3].split()[-1].rstrip().lstrip())
			stress_tensor = np.array([sigma11, sigma22, sigma33, sigma23, sigma13, sigma12])
	except FileNotFoundError:
		print('Stress Output file not found (File name: %s' % file_name)
		quit()
#	stress_tensor = np.array([(sigma11, sigma12, sigma13), (sigma12, sigma22, sigma23), (sigma13, sigma23, sigma33)])
	return stress_tensor

def calc_cij(stress_p, stress_n, Cij_index, dist):
	i = Cij_index
	d = float(dist)
	Cij = abs((stress_p[i] - stress_n[i])/(2*d))/10
	return Cij

def print_cij(cij_file, Cij):
#	cij_f = open(fname, 'a')
	for i in range(0,6):
		for j in range(0,6):
			cij_file.write('%8.3f   ' % Cij[i,j])
		cij_file.write('\n')
#	cij_f.close()
def print_cij_thermo(fname, Cij, vol):
	cij_f = open(fname, 'a')
	cij_f.write('%8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f\n' % (Cij[0,0], Cij[1,1], Cij[2,2], Cij[0,1], Cij[0,2], Cij[1,2], Cij[3,3], Cij[4,4], Cij[5,5]))
	cij_f.close()

