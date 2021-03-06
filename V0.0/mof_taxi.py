import os
from subprocess import call, Popen, PIPE
from operator import itemgetter
from subprocess import call, Popen, PIPE
from Queue import Queue, Empty
from threading import Thread
import numpy as np

LAMMPS_EXEC = os.environ.get('LAMMPS_EXEC')
RASPA_EXEC = os.environ.get('RASPA_DIR')+'bin/simulate'
RASPA_DIR = os.environ.get('RASPA_DIR')
TAXI_DIR = os.environ.get('TAXI_DIR')

def lmps2cif(f = None):
	if f is None:
		print('No file was passed')
		return
	elif not isinstance(f, str):
		print('f needs to be the string name of the file')
		return
	try:
		fl = open(f, 'r')
	except IOError:
		print('The lammps file specified does not exist')
		return
	name = f.split('.')[0]
	atoms = []
	atom_types = {}
	start_parsing = False
	start_typing = False
	skip = False
	for line in fl:
		row = line.split()
		if skip == True:
			skip = False
			continue
		if start_typing == True and len(row) > 0:
			atom_types.update({int(row[0]):row[3]})
		elif start_typing == True and len(row) == 0:
			start_typing = False
			continue
		if row == ['Masses']:
			skip = True
			start_typing = True
        	if len(row) > 3 and row[2] == 'xlo':
            		xlo = float(row[0])
            		xhi = float(row[1])
            		dx = xhi - xlo
        	elif len(row) > 3 and row[2] == 'ylo':
            		ylo = float(row[0])
            		yhi = float(row[1])
            		dy = yhi - ylo
        	elif len(row) > 3 and row[2] == 'zlo':
            		zlo = float(row[0])
            		zhi = float(row[1])
            		dz = zhi - zlo
	        elif len(row) > 0 and row[0] == 'Atoms':
	            start_parsing = True
		if start_parsing == True and len(row) > 3:
			x = (float(row[4]) - xlo)/dx
			y = (float(row[5]) - ylo)/dy
			z = (float(row[6]) - zlo)/dz
			charge = float(row[3])
			a_ind = int(row[0])
			at = atom_types[int(row[2])]
			atoms.append([at, x, y, z, charge, a_ind])
		if start_parsing == True and row == ['Velocities'] or start_parsing == True and row == ['Bonds']:
			break
	atoms = sorted(atoms, key = itemgetter(5))
	fl.close()
	fl = open(name+'.cif', 'w')
	fl.write('#======================================================================\n')
	fl.write('# CRYSTAL DATA\n')
	fl.write('#----------------------------------------------------------------------\n')
	fl.write(' \n')
	fl.write('data_' + name + '\n')
	fl.write('\n\n')
	fl.write('_pd_phase_name                         MOF\n')
	fl.write('_cell_length_a                         ' + str(dx) + '\n')
        fl.write('_cell_length_b                         ' + str(dy) + '\n')
        fl.write('_cell_length_c                         ' + str(dz) + '\n')
	fl.write('_cell_angle_alpha                      90.00000\n')
        fl.write('_cell_angle_beta                       90.00000\n')
        fl.write('_cell_angle_gamma                      90.00000\n')
	fl.write('_symmetry_space_group_name_H-M         "P 1" \n')
	fl.write('_symmetry_Int_Tables_number            1\n')	
	fl.write('loop_\n')
	fl.write('_symmetry_equiv_pos_as_xyz\n')
	fl.write('   "x,y,z"\n\n\n')
	fl.write('loop_\n')
	fl.write('   _atom_site_label\n')
	fl.write('   _atom_site_fract_x\n')
	fl.write('   _atom_site_fract_y\n')
	fl.write('   _atom_site_fract_z\n')
	fl.write('   _atom_site_charge\n')
	for i in atoms:	
		fl.write(i[0] + '\t' + str(i[1]) + '\t' +str(i[2]) + '\t' + str(i[3]) + '\t' + str(i[4]) + '\n')
	
	fl.close()
	return

def write_pseudo_atoms(f = None, adsorbate = None, path = None):
        if f is None:
                print('No file was passed')
                return
        elif not isinstance(f, str):
                print('f needs to be the string name of the file')
                return
        try:
                fl = open(f, 'r')
        except IOError:
                print('The lammps file specified does not exist')
                return
        if f is None:
                print('No adsorbate name was passed')
                return
        elif not isinstance(f, str):
                print('adsorbate needs to be the string name of the molecule')
                return
	if path is None:
		print('using default hard coded path')
		path1 = RASPA_DIR + 'share/raspa/molecules/TraPPE-UA/'
		path2 = RASPA_DIR + 'share/raspa/forcefield/TraPPE/'
	else:
		print('the path must be specified in the python code and not passed to the function')
	'''
	Get the framework information
	'''

        start_typing = False
        skip = False
	masses = []
	pseudo_atoms = []
	chem = []	
        for line in fl:
                row = line.split()
                if skip == True:
                        skip = False
                        continue
                if start_typing == True and len(row) > 0:
			masses.append(row[1])
			pseudo_atoms.append(row[3])
			if row[3][0] == 'L':
				chem.append(row[3][1])
			else:
				chem.append(row[3][0])
                elif start_typing == True and len(row) == 0:
                        start_typing = False
                        break
                if row == ['Masses']:
                        skip = True
                        start_typing = True
	fl.close()

	'''
	Get the adsorbate information
	'''

	fl = open(path1 + adsorbate+'.def', 'r')
	start = False
	temp_names = []
	for line in fl:
		row = line.split()
		if row == ['#', 'atomic', 'positions']:
			start = True
			continue
		if start == True:
			if row[0] != '#':
				temp_names.append(row[1])
			else:
				start = False
				break
	fl.close()
	temp_names = list(set(temp_names))
	temp_output = []
	fl = open(path2 + 'pseudo_atoms.def', 'r')
	for line in fl:
		row = line.split()
		if row[0] in temp_names:
			temp_output.append(line)
	fl.close()	

	'''
	write the pseudo_atoms.def file
	'''

	fl = open('pseudo_atoms.def', 'w')
	fl.write('#number of pseudo atoms\n')
	fl.write(' ' + str(len(masses)+1+len(temp_output)) + '\n')
	fl.write('#type   print   as    chem  oxidation   mass        charge   polarization B-factor    radii  connectivity anisotropic anisotropic-type   tinker-type\n')
	for i in range(0, len(masses), 1):
		fl.write(pseudo_atoms[i]+'\tyes\t'+chem[i]+'\t'+chem[i]+'\t0\t'+masses[i]+'  \t0.0      0.0          1.0      1.0    0            0           relative           0\n')
	fl.write('He      yes     He      He      0       4.002602        0.0      0.0          1.0      1.0    0            0           relative           0\n')
	for i in temp_output:
		fl.write(i)
	fl.close()
	return
def write_force_field():
	f=open('force_field.def', 'w')
	f.write('# rules to overwrite\n')
	f.write('0\n')
	f.write('# number of defined interactions\n')
	f.write('0\n')
	f.write('# type type2 interaction\n\n')
	f.write('# mixing rules to overwrite\n')
	f.write('0\n')
	f.close()
	return
	
def write_mixing_rules(f = None, adsorbate = None, path = None):
        if f is None:
                print('No file was passed')
                return
        elif not isinstance(f, str):
                print('f needs to be the string name of the file')
                return
        try:
                fl = open(f, 'r')
        except IOError:
                print('The lammps file specified does not exist')
                return
        if path is None:
                print('using default hard coded path')
                path1 = RASPA_DIR + 'share/raspa/molecules/TraPPE-UA/'
                path2 = RASPA_DIR + 'share/raspa/forcefield/TraPPE/'
        else:
                print('the path must be specified in the python code and not passed to the function')

        start_typing = False
        skip = False
        sigmas = []
        epsilons = []
        at = []

	'''
	mixing rules for the framework
	'''

        for line in fl:
                row = line.split()
                if skip == True:
                        skip = False
                        continue
                if start_typing == True and len(row) > 0:
			sigmas.append(float(row[2]))
			epsilons.append(float(row[1])/0.0019872041)
			at.append(row[4])
                elif start_typing == True and len(row) == 0:
                        start_typing = False
                        break
                if row == ['Pair', 'Coeffs']:
                        skip = True
                        start_typing = True
	fl.close()

	'''
	mixing rules for the adsorbate
	'''

        fl = open(path1 + adsorbate+'.def', 'r')
        start = False
        temp_names = []
        for line in fl:
                row = line.split()
                if row == ['#', 'atomic', 'positions']:
                        start = True
                        continue
                if start == True:
                        if row[0] != '#':
                                temp_names.append(row[1])
                        else:
                                start = False
                                break
        fl.close()
        temp_names = list(set(temp_names))
        fl.close()
        temp_names = list(set(temp_names))
        temp_output = []
        fl = open(path2 + 'force_field_mixing_rules.def', 'r')
        for line in fl:
                row = line.split()
                if row[0] in temp_names:
                        temp_output.append(line)
        fl.close()

	'''
	write force_field_mixing_rules.def
	'''

        fl = open('force_field_mixing_rules.def', 'w')
	fl.write('# general rule for shifted vs truncated\n')
	fl.write('truncated\n')
	fl.write('# general rule tailcorrections\n')
	fl.write('yes\n')
	fl.write('# number of defined interactions \n')
	fl.write(str(len(sigmas)+1+len(temp_output)) + '\n')
	fl.write('# type interaction, parameters.    IMPORTANT: define shortest matches first, so that more specific ones overwrites these\n')
	for i in range(0, len(sigmas), 1):
		fl.write(at[i] + ' \tlennard-jones \t' + str(epsilons[i]) + ' \t' +str(sigmas[i]) + '\n')
	fl.write('He      lennard-jones   10.9            2.64\n')
        for i in temp_output:
                fl.write(i)
	fl.write('# general mixing rule for Lennard-Jones\n')
	fl.write('Lorentz-Berthelot')
	fl.close()
	return

def combine_lmps_RaspaRestart(f_lmps = None, working_dir = None, adsorbate = None, cyclenumber = None):
	if working_dir is None:
		working_dir = os.getcwd()
	if f_lmps is None:
		print('You must provide a string indicating the initial lmps file')
		return
	if adsorbate is None:
		print('You must provide an adsorbate name that matchs the share RASPA lammps files')
		return
	else:
		try:
			f2 = open(RASPA_DIR + 'share/raspa/lmps/' + adsorbate + '.lmps', 'r')
		except:
			print('Unable to open file, check the location or adsorbate name')
			return
	if cyclenumber == 0:
		raspa_file = os.listdir(working_dir+'/Restart/System_0/')[0]
		f = open(working_dir + '/Restart/System_0/' + raspa_file, 'r')
	else:
                lst = os.listdir(working_dir+'/Restart/System_0/')
                for i in lst:
                        if len(i.split('_')) > 1:
                                raspa_file = i
                f = open(working_dir + '/Restart/System_0/' + raspa_file, 'r')
	'''
	Extract data from  the adsorbate RASPA .lmps file
	'''

	adsorbate_masses = []
	adsorbate_pair = []
	adsorbate_bond_types = []
	adsorbate_angle_types = []
	adsorbate_dihedral_types = []
	adsorbate_improper_types = []
	adsorbate_bonds = []
	adsorbate_angles = []
	adsorbate_dihedrals = []
	adsorbate_impropers = []
	start_m = False
	start_p = False
	start_a = False
	start_b = False
	start_d = False
	start_i = False
	find_b = False
	find_a = False
	find_d = False
	find_i = False
	for line in f2:
		row = line.split()
		if not line.strip():
                        continue
                if len(row) > 1:
                        if row[1] == 'atoms':
                                p_per = int(row[0])
			elif row[1] == 'bonds':
				b_per = int(row[0])
			elif row[1] == 'angles':
				a_per = int(row[0])
			elif row[1] == 'dihedrals':
				d_per = int(row[0])
			elif row[1] == 'impropers':
				i_per = int(row[0])
                if len(row) > 2:
			if row[1] == 'atom':
				pt_per = int(row[0])
			elif row[1] == 'bond':
				bt_per = int(row[0])
			elif row[1] == 'angle':
				at_per = int(row[0])
			elif row[1] == 'dihedral':
				dt_per = int(row[0])
			elif row[1] == 'improper':
				it_per = int(row[0])
		if row[0] == 'Masses':
			start_m = True
			continue
		elif row[0] == 'Pair':
			start_m = False
			start_p = True
			start_b = False
			start_a = False
			start_d = False
			start_i = False
		        find_b = False
 	 	        find_a = False
        		find_d = False
        		find_i = False
			continue
		elif row[0] == 'Bond':
                        start_m = False
                        start_p = False
                        start_b = True
                        start_a = False
                        start_d = False
                        start_i = False
                        find_b = False
                        find_a = False
                        find_d = False
                        find_i = False
			continue
		elif row[0] == 'Angle':
                        start_m = False
                        start_p = False
                        start_b = False
                        start_a = True
                        start_d = False
                        start_i = False
                        find_b = False
                        find_a = False
                        find_d = False
                        find_i = False
			continue
		elif row[0] == 'Dihedral':
                        start_m = False
                        start_p = False
                        start_b = False
                        start_a = False
                        start_d = True
                        start_i = False
                        find_b = False
                        find_a = False
                        find_d = False
                        find_i = False
			continue
		elif row[0] == 'Improper':
                        start_m = False
                        start_p = False
                        start_b = False
                        start_a = False
                        start_d = False
                        start_i = True
                        find_b = False
                        find_a = False
                        find_d = False
                        find_i = False
			continue
		elif row[0] == 'Atoms':
                        start_m = False
                        start_p = False
                        start_b = False
                        start_a = False
                        start_d = False
                        start_i = False
                        find_b = False
                        find_a = False
                        find_d = False
                        find_i = False
			continue
		elif row[0] == 'Bonds':
                        start_m = False
                        start_p = False
                        start_b = False
                        start_a = False
                        start_d = False
                        start_i = False
                        find_b = True
                        find_a = False
                        find_d = False
                        find_i = False
			continue
		elif row[0] == 'Angles':
                        start_m = False
                        start_p = False
                        start_b = False
                        start_a = False
                        start_d = False
                        start_i = False
                        find_b = False
                        find_a = True
                        find_d = False
                        find_i = False
			continue
		elif row[0] == 'Dihedrals':
                        start_m = False
                        start_p = False
                        start_b = False
                        start_a = False
                        start_d = False
                        start_i = False
                        find_b = False
                        find_a = False
                        find_d = True
                        find_i = False
			continue
		elif row[0] == 'Impropers':
                        start_m = False
                        start_p = False
                        start_b = False
                        start_a = False
                        start_d = False
                        start_i = False
                        find_b = False
                        find_a = False
                        find_d = False
                        find_i = True
			continue

		if start_m == True:
			adsorbate_masses.append(row)
		elif start_p == True:
			adsorbate_pair.append(row)
		elif start_b == True:
			adsorbate_bond_types.append(row)
		elif start_a == True:
			adsorbate_angle_types.append(row)
		elif start_d == True:
			adsorbate_dihedral_types.append(row)
		elif start_i == True:
			adsorbate_improper_types.append(row)
		elif find_b == True:
			adsorbate_bonds.append(row)
		elif find_a == True:
			adsorbate_angles.append(row)
		elif find_d == True:
			adsorbate_dihedrals.append(row)
		elif find_i == True:
			adsorbate_impropers.append(row)
	f2.close()
	'''
	Get the relevant adsorbate information from the restart file
	'''
	adsorbate_type_and_positions = []
	adsorbate_charge = []
	for line in f:
		if not line.strip():
			continue
		row = line.split()
		if row[0] == 'Component:':
			number_of_adsorbate_molecules = int(row[3])
			continue
		if row[0] == 'Adsorbate-atom-position:':
			adsorbate_type_and_positions.append([int(row[2]), float(row[3]), float(row[4]), float(row[5])])
			continue
		if row[0] == 'Adsorbate-atom-charge:':
			adsorbate_charge.append(float(row[3]))
			continue
	number_of_adsorbates = len(adsorbate_type_and_positions)
	f.close()
	'''
	Create the adsorbate type mapping from adsorbate.def file
	'''
	raspa_mapping = {}
	f = open(os.environ.get('RASPA_DIR')+'share/raspa/molecules/TraPPE-UA/'+adsorbate+'.def', 'r')
	start = False
	for line in f:
		if not line.strip():
                        continue
		row = line.split()
		if len(row) > 1 and row[1] == 'atomic':
			start = True
			continue
		if start == True:
			if row[0] == '#':
				start = False
			else:
				raspa_mapping.update({row[0]:row[1]})				
	f.close()
	'''
	open the lammps file and include the adsorbate information
	'''
	if cyclenumber == 0:
		f = open(f_lmps, 'r')
		x_off, y_off, z_off, dx, dy, dz = extract_framework_info(f_lmps)
	else:
		f = open(working_dir + '/Archive/' + str(cyclenumber-1)+'_aftermd.lmps', 'r')
                x_off, y_off, z_off, dx, dy, dz = extract_framework_info(working_dir + '/Archive/' + str(cyclenumber-1)+'_aftermd.lmps')
	f_n_at, f_n_b, f_n_a, f_n_d, f_n_at_t, f_n_b_t, f_n_a_t, f_n_d_t = extract_adsorbent_info(f_lmps)
	fl = open('temp_lmps.lmps','w')
	check = False
	holder = 0
	lammps_mapping = {}
	add_masses = False
	add_ptypes = False
	add_btypes = False
	add_atypes = False
	add_dtypes = False
	add_itypes = False
	skip_velocity = False
	dihedrals_written = False
	count = 0
	temporary_bond_list, temporary_angle_list, temporary_dihedral_list, temporary_improper_list = [], [], [], []
	for line in f:
		row = line.split()
		if not line.strip():
			if add_masses == True:
				if count < 1:
					count += 1
				else:
					for i in adsorbate_masses:
						fl.write(str(int(i[0])+f_n_at_t)+' '+i[1]+'\n') 
					add_masses = False
					count = 0
			if add_ptypes == True:
				if count < 1:
					count += 1
				else:
					for i in adsorbate_pair:
						fl.write(str(int(i[0])+f_n_at_t)+' '+i[1]+'\t'+i[2]+'\n')
						lammps_mapping.update({i[4]:str(int(i[0])+f_n_at_t)})
					add_ptypes = False
					count = 0
			if add_btypes == True:
                                if count < 1:
                                        count += 1
                                else:
                                        for i in adsorbate_bond_types:
                                                fl.write(str(int(i[0])+f_n_b_t)+' '+i[1]+'\t'+i[2]+'\n')
                                        add_btypes = False
                                        count = 0
                        if add_atypes == True:
                                if count < 1:
                                        count += 1
                                else:
                                        for i in adsorbate_angle_types:
                                                fl.write(str(int(i[0])+f_n_a_t)+' '+i[1]+'\t'+i[2]+'\n')
                                        add_atypes = False
                                        count = 0
                        if add_dtypes == True:
                                if count < 1:
                                        count += 1
                                else:
                                        for i in adsorbate_dihedral_types:
                                                fl.write(str(int(i[0])+f_n_d_t)+' '+i[1]+'\t'+i[2]+' '+i[3]+'\n')
                                        add_dtypes = False
                                        count = 0
                        if add_itypes == True:
                                if count < 1:
                                        count += 1
                                else:
                                        for i in adsorbate_improper_types:
                                                fl.write(str(int(i[0])+n_improper_types)+' '+i[1]+'\t'+i[2]+' '+i[3]+'\n')
                                        add_itypes = False
                                        count = 0
			fl.write('\n')
			continue
		if row[0] == 'Velocities':
			skip_velocity = True
			continue
		if skip_velocity == True and row[0] != 'Bonds':
			continue
		else:
			skip_velocity = False	
		if len(row) > 1:
			if row[1] == 'atoms':
				n_atoms = int(row[0])
				fl.write(str(f_n_at+number_of_adsorbates)+ ' atoms\n')
				continue
			elif row[1] == 'bonds':
				n_bonds = int(row[0])
				fl.write(str(f_n_b+number_of_adsorbate_molecules*len(adsorbate_bonds))+' bonds\n')
				continue
			elif row[1] == 'angles':
				n_angles = int(row[0])
				fl.write(str(f_n_a+number_of_adsorbate_molecules*len(adsorbate_angles))+' angles\n')
				continue
			elif row[1] == 'dihedrals':
				n_dihedrals = int(row[0])
				fl.write(str(f_n_d+number_of_adsorbate_molecules*len(adsorbate_dihedrals))+' dihedrals\n')
				continue
			elif row[1] == 'impropers':
				n_impropers = int(row[0])
				fl.write(str(n_impropers+number_of_adsorbate_molecules*len(adsorbate_impropers))+' impropers\n')
				continue
		if len(row) > 2:
			if row[1] == 'atom':
				fl.write(str(f_n_at_t+len(adsorbate_masses))+ ' atom types\n')
				continue
			elif row[1] == 'bond':
                                fl.write(str(f_n_b_t+len(adsorbate_bond_types))+ ' bond types\n')
				continue
			elif row[1] == 'angle':
                                n_angle_types = int(row[0])
                                fl.write(str(f_n_a_t+len(adsorbate_angle_types))+ ' angle types\n')
				continue
			elif row[1] == 'dihedral':
				n_dihedral_types = int(row[0])
                                fl.write(str(f_n_d_t+len(adsorbate_dihedral_types))+ ' dihedral types\n')
				continue
			elif row[1] == 'improper':
                                n_improper_types = int(row[0])
                                fl.write(str(n_improper_types+len(adsorbate_improper_types))+ ' improper types\n')
				continue
		if row[0] == 'Masses':
			add_masses = True
			fl.write(line)
			continue
		elif row[0] == 'Pair':
			add_ptypes = True
			fl.write(line)
			continue
		elif row[0] == 'Bond':
                        add_btypes = True
			fl.write(line)
                        continue
		elif row[0] == 'Angle':
			add_atypes = True
                        fl.write(line)
                        continue
		elif row[0] == 'Dihedral':
			add_dtypes = True
                        fl.write(line)
                        continue
		elif row[0] == 'Improper':
			add_itypes = True
                        fl.write(line)
                        continue

		if len(row) == 7 or len(row) == 10:
			count = 0
                        if int(row[0]) > f_n_at:
				continue
			if int(row[1]) > holder:
				holder = int(row[1])
			if len(row) == 10:
				fl.write(str(row[0]) + '\t' + str(row[1]) + '\t' + str(row[2]) + '\t' + str(row[3]) + '\t' +str(row[4]) +'\t'+ str(row[5]) +'\t'+ str(row[6]) + '\n')
			else:
				fl.write(line)
			if int(row[0]) == f_n_at:
				for i in range(0, len(adsorbate_type_and_positions), 1):
					z_hold = 0
					y_hold = 0
					x_hold = 0
					fl.write(str(i+1+f_n_at) + ' ')
                                        fl.write(str(i/p_per+1 + holder) + ' ')
					fl.write(lammps_mapping[raspa_mapping[str(adsorbate_type_and_positions[i][0])]]+' ')
					fl.write(str(round(adsorbate_charge[i],4)) + ' ')
					if adsorbate_type_and_positions[i][1] + x_off < x_off:
						x_hold -= 1
					if adsorbate_type_and_positions[i][2] + y_off < y_off:
                                                y_hold -= 1
                                        if adsorbate_type_and_positions[i][3] + z_off < z_off:
                                                z_hold -= 1
					fl.write(str(round(adsorbate_type_and_positions[i][1],4)+x_off)+' '+str(round(adsorbate_type_and_positions[i][2],4)+y_off)+' '+str(round(adsorbate_type_and_positions[i][3],4)+z_off)+ '\n')	
			continue
		if len(row) == 4 and row[3].isdigit():
			if add_dtypes == True:
				fl.write(line)
				continue
			if int(row[0]) == 1:
				for i in range(0,number_of_adsorbate_molecules,1):
					for t in adsorbate_bonds:
						temp_line = str(int(t[0])+f_n_b+i*len(adsorbate_bonds))+'\t' + str(int(t[1])+f_n_b_t)+'\t' + str(int(t[2])+f_n_at+i*p_per)+'\t'+str(int(t[3])+f_n_at+i*p_per)+'\n'
						temporary_bond_list.append(temp_line)
                        if int(row[1]) > f_n_b_t:
                                continue
                        temporary_bond_list.append(line)
			continue
		if len(row) == 5 and row[3].isdigit():
			if int(row[0]) == 1:
				for i in range(0, number_of_adsorbate_molecules, 1):
					for t in adsorbate_angles:
						temp_line = str(int(t[0])+f_n_a+i*len(adsorbate_angles))+'\t' + str(int(t[1])+f_n_a_t)+'\t' + str(int(t[2])+f_n_at+i*p_per)+'\t'+str(int(t[3])+f_n_at+i*p_per)+'\t'+str(int(t[4])+f_n_at+i*p_per)+'\n'
						temporary_angle_list.append(temp_line)
                        if int(row[1]) > f_n_a_t:
                                continue
                        temporary_angle_list.append(line)
			continue
		if row[0] == 'Impropers':
			dihedrals_written = True
		if len(row) == 6 and row[4].isdigit() and dihedrals_written == False:
			if int(row[0]) == 1:
				for i in range(0, number_of_adsorbate_molecules, 1):
					for t in adsorbate_dihedrals:
						temp_line = str(int(t[0])+f_n_d+i*len(adsorbate_dihedrals))+'\t' + str(int(t[1])+f_n_d_t)+'\t' + str(int(t[2])+f_n_at+i*p_per)+'\t'+str(int(t[3])+f_n_at+i*p_per)+'\t'+str(int(t[4])+f_n_at+i*p_per)+'\t'+str(int(t[5])+f_n_at+i*p_per)+'\n'
						temporary_dihedral_list.append(temp_line)
                        if int(row[1]) > f_n_d_t:
                                continue
                        temporary_dihedral_list.append(line)
			continue
		elif len(row) == 6 and row[4].isdigit() and dihedrals_written == True:
                        fl.write(line)
                        if int(row[0]) == 1:
                                for i in range(0, number_of_adsorbate_molecules, 1):
                                        for t in adsorbate_impropers:
                                                fl.write(str(int(t[0])+n_impropers+i*len(adsorbate_impropers))+'\t')
                                                fl.write(str(int(t[1])+n_improper_types)+'\t')
                                                fl.write(str(int(t[2])+f_n_at+i*p_per)+'\t'+str(int(t[3])+f_n_at+i*p_per)+'\t'+str(int(t[4])+f_n_at+i*p_per)+'\t'+str(int(t[5])+f_n_at+i*p_per)+'\n')
                        continue
		if add_masses == True:
			if int(float(row[0])) > f_n_at_t:
				continue
		if add_ptypes == True:
			if int(float(row[0])) > f_n_at_t:
                                continue
                if add_btypes == True:
                        if int(float(row[0])) > f_n_b_t:
                                continue
                if add_atypes == True:
                        if int(float(row[0])) > f_n_a_t:
                                continue
                if add_dtypes == True:
                        if int(float(row[0])) > f_n_d_t:
                                continue
		if row[0] == 'Bonds' or row[0] == 'Angles' or row[0] == 'Dihedrals':
			continue
		fl.write(line)
	fl.write('\nBonds \n \n')
	count = 1
	print(number_of_adsorbate_molecules)
	for i in temporary_bond_list:
		row = i.split()
		fl.write(str(count) + '\t')
		fl.write(row[1] + '\t')
		fl.write(row[2] + '\t')
		fl.write(row[3] + '\n')	
		count += 1
        fl.write('\nAngles \n \n')
	count = 1
        for i in temporary_angle_list:
                row = i.split()
                fl.write(str(count) + '\t')
                fl.write(row[1] + '\t')
                fl.write(row[2] + '\t')
                fl.write(row[3] + '\t')
                fl.write(row[4] + '\n')
                count += 1
        fl.write('\nDihedrals \n \n')
        count = 1
        for i in temporary_dihedral_list:
                row = i.split()
                fl.write(str(count) + '\t')
                fl.write(row[1] + '\t')
                fl.write(row[2] + '\t')
                fl.write(row[3] + '\t')
                fl.write(row[4] + '\t')
                fl.write(row[5] + '\n')
                count += 1
	return

def updateRaspaRestart(f_lmps = None, num_framework_atoms = None, num_adsorbate_types = None, num_framework_types = None, working_dir = None, cyclenumber=None):
	if working_dir is None:
		working_dir = os.getcwd()
	if f_lmps is None:
		print('you must provide the name of the lammps file post MD')
		return
	if num_framework_atoms is None:
		print('you must specify how many atoms are in the adosrbent material')
		return
	if cyclenumber is None:
		print('you must specify the cyclenumber')
		return


	'''
	Get atoms from the LAMMPS file after the MD run
	'''

	f = open(f_lmps, 'r')
	atoms = []
	box_bounds = []
	store_atoms = False
	count = 0
	for line in f:
		if not line.strip():
			if store_atoms == True:
				if count < 1:
					count += 1
					continue
				else:
					store_atoms = False
					count = 0
					break
                        continue
		row = line.split()
		if row[0] == 'Atoms':
			store_atoms = True
			continue
		if store_atoms == True:
			if int(row[0]) <= num_framework_atoms:
				continue
			atoms.append([int(row[0])-num_framework_atoms-1, int(row[2])-num_framework_types-num_adsorbate_types, row[3], str(float(row[4])-xlo), str(float(row[5])-ylo), str(float(row[6])-zlo)])
		if len(row) > 2:
			if row[2] == 'xlo':
				box_bounds.append(float(row[1])-float(row[0]))
				xlo = float(row[0])
				xhi = float(row[1])
			elif row[2] == 'ylo':
				box_bounds.append(float(row[1])-float(row[0]))
				ylo = float(row[0])
				yhi = float(row[1])
			elif row[2] == 'zlo':
				box_bounds.append(float(row[1])-float(row[0]))
				zlo = float(row[0])
				zhi = float(row[1])
	atoms = sorted(atoms, key = itemgetter(0))
	f.close()
	'''
	open restart file and write a new one with  new adsorbate positions
	'''
	if cyclenumber == 0:
       		raspa_file = os.listdir(working_dir+'/Restart/System_0/')[0]
        	f = open(working_dir + '/Restart/System_0/' + raspa_file, 'r')
		fl = open(working_dir + '/Restart/System_0/temprestart.txt', 'w')
        else:
		lst = os.listdir(working_dir+'/Restart/System_0/')
		for i in lst:
			if len(i.split('_')) > 1:
        		        raspa_file = i
                f = open(working_dir + '/Restart/System_0/' + raspa_file, 'r')
                fl = open(working_dir + '/Restart/System_0/temprestart.txt', 'w')
	count = 0
	for line in f:
		if not line.strip():
			fl.write('\n')
			continue
		row = line.split()
		if row[0] == 'unit-cell-vector-a:' or row[0] == 'cell-vector-a:':
			fl.write(row[0]+'\t'+str(box_bounds[0])+'\t'+row[2]+'\t'+row[3]+'\n')
			continue
		elif row[0] == 'unit-cell-vector-b:' or row[0] == 'cell-vector-b:':
                        fl.write(row[0]+'\t'+row[1]+'\t'+str(box_bounds[1])+'\t'+row[3]+'\n')
                        continue
                elif row[0] == 'unit-cell-vector-c:' or row[0] == 'cell-vector-c:':
                        fl.write(row[0]+'\t'+row[1]+'\t'+row[2]+'\t'+str(box_bounds[2])+'\n')
                        continue
		elif row[0] == 'cell-lengths:':
			fl.write(row[0] + '\t'+str(box_bounds[0])+'\t'+str(box_bounds[1])+'\t'+str(box_bounds[2])+'\n')
			continue
		elif row[0] == 'Adsorbate-atom-position:':
			fl.write(row[0]+' '+row[1]+' '+row[2]+'\t'+str(atoms[count][3])+'\t'+str(atoms[count][4])+'\t'+str(atoms[count][5])+'\n')
			count += 1
			continue
		elif row[0] == 'Adsorbate-atom-force:':
			continue
		else:
			fl.write(line)
	fl.close()
	'''	
	rename the restart files appropriately
	'''

	if cyclenumber == 0:
		os.system('cp -avr '+working_dir+'/Restart/System_0/'+raspa_file+' '+working_dir+'/Archive/'+str(cyclenumber)+'.restart')
                os.system('cp -avr '+working_dir+'/Restart/System_0/temprestart.txt '+working_dir+'/Restart/System_0/'+raspa_file)
	else:
                os.system('cp -avr '+working_dir+'/RestartInitial/System_0/'+raspa_file+' '+working_dir+'/Archive/'+str(cyclenumber)+'.restart')
                os.system('cp -avr '+working_dir+'/Restart/System_0/temprestart.txt '+working_dir+'/RestartInitial/System_0/'+raspa_file)
	return

def updatecif(f_lmps=None, f_cif=None, num_framework_atoms=None, cyclenumber=None):
	if f_lmps is None:
		print('you must specify the .lmps file to update the cif')
		return
        if f_cif is None:
                print('you must specify the .cif file to be updated')
                return
	if num_framework_atoms is None:
		print('you must specify the number of framework atoms')
		return
	if cyclenumber is None:
		print('you must specify the cycle number')
		return
	f = open(f_lmps, 'r')
	fl = open(f_cif, 'r')
        atoms = []
        start_parsing = False
	for line in f:
                if not line.strip():
                        continue
		row = line.split()
                if len(row) > 3 and row[2] == 'xlo':
                        xlo = float(row[0])
                        xhi = float(row[1])
                        dx = xhi - xlo
                elif len(row) > 3 and row[2] == 'ylo':
                        ylo = float(row[0])
                        yhi = float(row[1])
                        dy = yhi - ylo
                elif len(row) > 3 and row[2] == 'zlo':
                        zlo = float(row[0])
                        zhi = float(row[1])
                        dz = zhi - zlo
                elif len(row) > 0 and row[0] == 'Atoms':
                    start_parsing = True
                if start_parsing == True and len(row) > 3 and int(row[0]) <= num_framework_atoms:
                        x = (float(row[4]) - xlo)/dx
                        y = (float(row[5]) - ylo)/dy
                        z = (float(row[6]) - zlo)/dz
                        charge = float(row[3])
                        a_ind = int(row[0])
                        atoms.append([x, y, z, charge, a_ind])
                if start_parsing == True and row == ['Velocities'] or start_parsing == True and row == ['Bonds']:
                        break
	f.close()
        atoms = sorted(atoms, key = itemgetter(4))
        f = open(f_cif, 'r')
	fl = open('temp_cif.cif', 'w')
	count = 0
	for line in f:
		if not line.strip():
			fl.write('\n')
                        continue
		row = line.split()
		if row[0] == '_cell_length_a':
		        fl.write('_cell_length_a                         ' + str(dx) + '\n')
			continue
		elif row[0] == '_cell_length_b':
                        fl.write('_cell_length_b                         ' + str(dy) + '\n')
			continue
                elif row[0] == '_cell_length_c':
                        fl.write('_cell_length_c                         ' + str(dz) + '\n')
                        continue
		elif len(row) == 5:
			fl.write(row[0]+'\t'+str(atoms[count][0])+'\t'+str(atoms[count][1])+'\t'+str(atoms[count][2])+'\t'+str(atoms[count][3])+'\n')
			count += 1
			continue
		fl.write(line)
	f.close()
	fl.close()
	os.system('mv '+f_cif+' '+str(cyclenumber)+'.cif')
	os.system('mv temp_cif.cif ' +f_cif)
	return

def get_adsorbate_connectivity(adsorbate = None, path = None):
	if path is None:
		print('Using default RASPA path')
		f = open(RASPA_DIR + 'raspa/molecules/TraPPE-UA/' + adsorbate + '.def')
	else:
		f = open(path + adsorbate + '.def')
	start_stretch = False
	start_bend = False
	start_torsion = False
	stretch, bend, torsion = [],[],[]
	count = 0
	for line in f:
		count += 1
		if not line.strip():
			continue
		row = line.split()
		if count == 6:
			number_of_atoms = row[0]
		if len(row) <= 2:
			continue
		if row[0] == '#' and start_stretch == True or row[0] == '#' and start_bend == True or row[0] == '#' and start_torsion == True:
		        start_stretch = False
		        start_bend = False
		        start_torsion = False
		if row[2] == 'stretch:':
			start_stretch = True
			continue
                elif row[2] == 'bending:':
                        start_bend = True
                        continue
                elif row[1] == 'Torsion':
                        start_torsion = True
                        continue
		
		if start_stretch == True:
			stretch.append([row[0], row[1]])
		elif start_bend == True:
			bend.append([row[0], row[1], row[2]])
		elif start_torsion == True:
			torsion.append([row[0], row[1], row[2], row[3]])
	bonding = [stretch, bend, torsion]
	return number_of_atoms, bonding
def write_raspa_input(temp=None, pressure=None, cyclenumber=None, cycles=None,framework=None, adsorbate=None, rosenbluth=None):
	if temp is None:
		print('you must specify a temperature')
		return
	if pressure is None:
		print('you must specify a pressure')
		return
	if cyclenumber is None:
		print('you must specify a cycle number')
		return
	if cycles is None:
		print('you must specify the number of cycles to run RASPA for')
		return
	if adsorbate is None:
		print('you must specify the name of the adsorbate specieis')
		return
	if framework is None:
		print('you must specify the name of the framework material')
		return
	if rosenbluth is None:
		print('defaulting to a rosenbluth weight of 1.000')
		rosenbluth = '1.000'
	f = open('simulation.input', 'w')

	f.write('SimulationType                MonteCarlo\n')
	f.write('NumberOfCycles                ' + str(cycles)+'\n')
	f.write('NumberOfInitializationCycles  ' + str(cycles)+'\n')
	f.write('PrintEvery                    100\n')

	if cyclenumber!=0:
		f.write('RestartFile                   yes\n\n')
	else:
		f.write('\n')

	f.write('Forcefield                    polymer\n')
	f.write('UseChargesFromCIFFile         yes\n\n')
	f.write('Framework                     0\n')
	f.write('FrameworkName                 '+framework+'\n')
	f.write('UnitCells                     1 1 1\n')
	f.write('HeliumVoidFraction            0.3\n')
	f.write('ExternalTemperature           '+str(temp)+'\n')
	f.write('ExternalPressure              '+str(pressure)+'\n\n')
	f.write('Movies                        no\n\n')

	f.write('Component 0 MoleculeName             '+adsorbate+'\n')
	f.write('MoleculeDefinition       TraPPE-UA\n')
	f.write('IdealGasRosenbluthWeight ' + str(rosenbluth) + '\n')
	f.write('TranslationProbability   1.0\n')
	f.write('RotationProbability      1.0\n')
	f.write('ReinsertionProbability   1.0\n')
	f.write('SwapProbability          1.0\n')
	f.close()
	return

def write_lammps_input(temp=None, pressure=None, md_time=None, cyclenumber=None, mode = None):
        if temp is None:
                print('you must specify a temperature')
                return
        if pressure is None:
                print('you must specify a pressure')
                return
        if cyclenumber is None:
                print('you must specify a cycle number')
                return
	if md_time is None:
		print('you must specify a length of MD simulation')
		return
	if mode is None:
		print('Defaulting to explicit hydrogen timestep')
		mode = 'EH'
	'''
	convert pressure to atm
	'''
	pressure*=9.86923e-6

	f=open('npt.in', 'w')
	f.write('units           real\n')
	f.write('boundary        p p p\n\n')
	f.write('atom_style full\n')
	f.write('pair_style lj/cut/coul/long 14.0 14.0\n')
	f.write('kspace_style pppm 1e-5\n')
	f.write('pair_modify shift yes mix arithmetic\n')
	f.write('bond_style harmonic\n')
	f.write('angle_style harmonic\n')
	f.write('dihedral_style harmonic\n')
	f.write('improper_style harmonic\n')
	f.write('special_bonds dreiding\n')
	f.write('read_data temp_lmps.lmps\n')
	f.write('neighbor 1.0 bin\n')
	f.write('neigh_modify once no every 1 delay 0 check yes\n\n')
	f.write('thermo_style     custom step vol temp press etotal density\n')
	f.write('thermo           1000\n\n')
	f.write('min_style        quickmin\n')
	f.write('minimize         1.0e-6 1.0e-6 100000 100000\n')
	f.write('fix              1 all npt temp '+str(temp)+' '+str(temp)+' 100 aniso '+str(pressure)+' '+str(pressure)+' 1000\n')
	f.write('velocity         all create '+str(temp)+' 577263458 dist gaussian\n')
	f.write('run              0\n')
	f.write('velocity         all scale '+str(temp)+'\n')
	if mode == 'EH':
		f.write('run              '+str(md_time)+'\n')
	elif mode == 'UA':
		f.write('timestep 1.0\n')
		f.write('run              '+str(int(md_time*0.2))+'\n')
		f.write('unfix         1\n')
		f.write('fix              1 all npt temp '+str(temp)+' '+str(temp)+' 200 aniso '+str(pressure)+' '+str(pressure)+' 2000\n')
		f.write('timestep 2.0\n')
		f.write('run              '+str(int(md_time*0.8))+'\n')
	f.write('unfix            1\n\n')
	f.write('write_data       ' +str(cyclenumber)+'.lmps\n')
	f.close()
	return
def run_raspa():
	global RASPA_EXEC
	p = Popen([RASPA_EXEC, 'simulation.input'], stdin=PIPE, stdout=PIPE, stderr=PIPE).communicate()
	del(p)
	return
def run_lammps(prefix=None, np=None):
	if prefix is None:
		print('Defaulting to the srun --mpi=pmix_v2 prefix')
		prefix=['srun', '--mpi=pmix_v2']
	global LAMMPS_EXEC
	f = open('npt.in', 'r').read()
        p = Popen(prefix + [LAMMPS_EXEC, '-e', 'both'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
	stdo, stde = p.communicate(f)
	del(p)
	return
def archive(cwd=None, cyclenumber=None):
	if cwd is None:
		cwd = os.getcwd()
	if os.path.isdir(cwd+'/Archive/'):
		print('archiving the results')
	else:
		os.system('mkdir Archive')
	os.system('mv '+str(cyclenumber)+'.lmps Archive/'+str(cyclenumber)+'_aftermd.lmps')
	lst=os.listdir(cwd+'/Output/System_0')
	for i in lst:
		if i.split('_')[0] == 'output':
			os.system('cp Output/System_0/'+i+' Archive/'+str(cyclenumber)+'.data')
	os.system('mv log.lammps Archive/'+str(cyclenumber)+'_lammps.log')
	os.system('mv '+str(cyclenumber)+'.cif Archive/'+str(cyclenumber)+'.cif')
	os.system('rm -r VTK Movies temp_lmps.lmps')
	return
def get_parameters(f_lmps = None, adsorbate=None):
	if f_lmps is None:
		print('You must specify the lammps file')
		return
	if adsorbate is None:
		print('You must provide the name of the adsorbate')
		return
	f = open(f_lmps, 'r')
	for line in f:
                row = line.split()
                if not line.strip():
                        continue
                if len(row) > 1:
                        if row[1] == 'atoms':
                                framework_atoms = int(row[0])
			if row[1] == 'atom':
				framework_types = int(row[0])
	f.close()
	f = open(RASPA_DIR+'share/raspa/lmps/'+adsorbate+'.lmps', 'r')
	for line in f:
                row = line.split()
                if not line.strip():
                        continue
                if len(row) > 1:
                        if row[1] == 'atom':
                                adsorbate_types = int(row[0])
	return framework_atoms, framework_types, adsorbate_types
def get_cycle():
	cwd = os.getcwd()
	lst = os.listdir(cwd+'/Archive/')
	cycle = 0
	for i in lst:
		if i.split('.')[0].isdigit():
			if int(i.split('.')[0]) > cycle:
				cycle = int(i.split('.')[0])
	return cycle
def averages():
	cwd = os.getcwd()
	lst = os.listdir(cwd+'/Archive/')
	hold = []
	for i in lst:
		if i.split('.')[-1] != 'data':
			hold.append(i)
		else:
			continue
	for i in hold:
		lst.remove(i)
	output = []
	start = False
	count = 0
	for i in range(0, len(lst), 1):
		f = open(cwd+'/Archive/'+str(i)+'.data','r')
		for line in f:
			row = line.split()
		        if not line.strip():
				continue
			if row[0] == 'Average' and row[1] == 'loading' and row[2] == 'absolute' and row[3] == '[mol/kg':
				print('I found one')
				hold_3 = (row[5])
				hold_4 = (row[7])
			if row[0] == 'absolute':
				start = True
				try:
					hold_1 = row[8].split(')')[0]
				except:
					pass
			elif start == True:
				start = False
				try:
					hold_2 = row[7].split(')')[0]
				except:
					pass
		output.append([str(count), hold_1, hold_2, hold_4, hold_3])
		count += 1
		f.close()
	f = open('average_loading.txt', 'w')
	for i in output:
		for t in i:
			f.write(t+'\t')
		f.write('\n')
	f.close()
	return
def r_factor(cyclenumber=None):
	if cyclenumber is None:
		print('A cyclenumber must be provided')
		return
	f = open('r_factor.txt', 'w')
	if cyclenumber < 4:
		f.write('1000')
		r_factor = 1000
	else:
		means = []
		variance = []
		average_file = open('average_loading.txt', 'r')
		for line in average_file:
			row = line.split()
			if float(row[0]) >= cyclenumber-4:
				print(row[0])
				means.append(float(row[4]))
				variance.append(float(row[3])**2)
		r_factor = 1+6.0/5.0*(np.std(means)**2)/(np.mean(variance))
		f.write(str(r_factor))
	        average_file.close()
	f.close()
	return r_factor

def write_mc_md(temperature = None, pressure = None, adsorbate = None, adsorbent = None, rosenbluth=None, mode = 'UA'):
	f = open('MCMD.py', 'w')
	f.write('import taxi\n')
	f.write('import os\n\n')
	f.write('cwd = os.getcwd()\n')
	f.write('temperature = ' + str(temperature) +'\n')
	f.write('pressure = ' + str(pressure) + '\n\n')
	f.write('adsorbate = "' + str(adsorbate) +'"'+'\n')
	f.write('f_prefix = "' + str(adsorbent) +'"'+'\n')
	f.write("f_lmps = f_prefix+'.lmps'\n")
	f.write("f_cif = f_prefix+'.cif'\n")
	f.write('fa, ft, at = taxi.get_parameters(f_lmps=f_lmps, adsorbate=adsorbate)\n')
        f.write("r_factor = 1000\n")
	f.write("if os.path.isdir('Archive'):\n")
	f.write("\tlst=os.listdir(cwd+'/Archive/')\n")
	f.write("\tstart_cycle = taxi.get_cycle()\n")
        f.write("\tlst = os.listdir(cwd+'/Restart/System_0/')\n")
	f.write("\ttry:\n")
	f.write("\t\tlst.remove('temprestart.txt')\n")
	f.write("\texcept:\n")
	f.write("\t\tpass\n")
        f.write("\tstart_cycle -= 1\n")
        f.write("\tos.system('cp ' + cwd+'/Archive/'+str(start_cycle)+'.restart '+cwd+'/Restart/System_0/'+lst[0])\n")
	f.write("\tos.system('cp ' + cwd+'/Archive/'+str(start_cycle)+'.cif '+cwd+'/'+f_cif)\n")
        f.write("\ttaxi.updateRaspaRestart(f_lmps=cwd+'/Archive/'+str(start_cycle-1)+'_aftermd.lmps',num_framework_atoms=fa,num_framework_types=ft,num_adsorbate_types=at,cyclenumber = start_cycle)\n")
	f.write("\ttaxi.updatecif(f_lmps=cwd+'/Archive/'+str(start_cycle-1)+'_aftermd.lmps',f_cif=f_cif, num_framework_atoms=fa,cyclenumber=start_cycle)\n")
	f.write('else:\n')
	f.write("\tos.system('mkdir Archive')\n")
        f.write('\ttaxi.write_pseudo_atoms(f=f_lmps, adsorbate=adsorbate)\n')
        f.write('\ttaxi.write_mixing_rules(f=f_lmps, adsorbate=adsorbate)\n')
	f.write('\ttaxi.write_force_field()\n')
	f.write("\tstart_cycle = 0\n\n")
	f.write('while r_factor > 1.5:\n')
	f.write('\ti = start_cycle\n')
	f.write('\ttaxi.write_raspa_input(temp=temperature,pressure=pressure,cyclenumber=i,cycles=5000,adsorbate=adsorbate,framework=f_prefix,rosenbluth='+str(rosenbluth)+')\n')
	f.write("\tos.system('rm -r Restart')\n")
        f.write('\ttaxi.run_raspa()\n')
	f.write('\ttaxi.combine_lmps_RaspaRestart(f_lmps=f_lmps, adsorbate=adsorbate,cyclenumber=i)\n')
	if mode == 'UA':
        	f.write('\ttaxi.write_lammps_input(temp=temperature, pressure=pressure, md_time=500000, cyclenumber=i, mode="UA")\n')
	else:
                f.write('\ttaxi.write_lammps_input(temp=temperature, pressure=pressure, md_time=1000000, cyclenumber=i, mode="EH")\n')
        f.write('\ttaxi.run_lammps()\n')
        f.write("\ttaxi.updatecif(f_lmps=str(i)+'.lmps',f_cif=f_cif, num_framework_atoms=fa,cyclenumber=i)\n")
        f.write("\ttaxi.updateRaspaRestart(f_lmps=str(i)+'.lmps',num_framework_atoms=fa,num_framework_types=ft,num_adsorbate_types=at,cyclenumber = i)\n")
        f.write('\ttaxi.archive(cyclenumber=i)\n')
        f.write('\tif i == 0:\n')
        f.write("\t\tos.system('mv Restart RestartInitial')\n")
	f.write('\ttaxi.averages()\n')
	f.write('\tr_factor = taxi.r_factor(cyclenumber=i)\n')
	f.write('\tstart_cycle += 1\n')
	f.close()
	return

def screening_setup(temp = None, adsorbent=None, adsorbate=None, number_of_samples=None, mode = None):
	if temp is None:
		print('The user must specify a temperature')
		return
	if adsorbent is None:
		print('The user must specify an adsorbent material')
		return
	if adsorbate is None:
		print('The user must specify the adsorbate species')
		return
	if number_of_samples is None:
		print('The user must specify the number adsorbent sample replicates')
		return
	if mode is None:
		print('Defaulting to explicit hydrogen timetstep')
		mode = 'EH'
	try:
		f=open('pressure.txt', 'r')
	except:
		print('Something went wrong reading the pressure file!')
		return
	cwd = os.getcwd()
	temperature = temp
	pressures = []
	for line in f:
		row = line.split()
		pressures.append(float(row[0]))	
	f.close()
	try:
		f=open(TAXI_DIR+'/rosenbluth/'+adsorbate, 'r')
		for line in f:
			row = line.split()
			rosenbluth = row[0]
		f.close()
	except:
		print('Something is wrong with the adsorbate rosenbluth file')
		return
	if os.path.isdir(TAXI_DIR+'/adsorbent/'+adsorbent):
		pass
	else:
		print('No directory in TAXI_DIR/adsorbent that matches the adsorbent name')
		return	
	for p in pressures:
		os.system('mkdir ' +str(p))
		for nm in range(1, number_of_samples+1, 1):
			os.system('mkdir '+cwd+'/'+str(p)+'/'+str(nm))
			write_mc_md(temperature=temp, pressure=str(p), adsorbate=adsorbate, adsorbent=adsorbent, rosenbluth=rosenbluth, mode=mode)
			os.system('mv MCMD.py ' + cwd+'/'+str(p)+'/'+str(nm))
			os.system('cp '+TAXI_DIR+'/adsorbent/'+adsorbent+'/'+str(temp)+'/'+str(nm)+'.lmps '+cwd+'/'+str(p)+'/'+str(nm)+'/'+adsorbent+'.lmps')
                        os.system('cp '+TAXI_DIR+'/adsorbent/'+adsorbent+'/'+str(temp)+'/'+str(nm)+'.cif '+cwd+'/'+str(p)+'/'+str(nm)+'/'+adsorbent+'.cif')
			os.system('cp '+TAXI_DIR+'/submit/submit.sh '+cwd+'/'+str(p)+'/'+str(nm)+'/')
	return
def extract_framework_info(f = None):
	if f is None:
		print('A framework must be given to extract from')
		return
	framework = open(f, 'r')
	for line in framework:
                if line.strip():
			row = line.split()
			if len(row) == 4:
				if row[2] == 'xlo':
					x_off = float(row[0])
					dx = float(row[1]) - x_off
				elif row[2] == 'ylo':
					y_off = float(row[0])
                                        dy = float(row[1]) - y_off
				elif row[2] == 'zlo':
					z_off = float(row[0])
                                        dz = float(row[1]) - z_off
	return  x_off, y_off, z_off, dx, dy, dz
def extract_adsorbent_info(f = None):
        if f is None:
                print('A framework must be given to extract from')
                return
        framework = open(f, 'r')
        for line in framework:
                if line.strip():
                        row = line.split()
                        if len(row) == 2:
                                if row[1] == 'atoms':
                                        n_atoms = int(row[0])
                                elif row[1] == 'bonds':
                                        n_bonds = int(row[0])
                                elif row[1] == 'angles':
                                        n_angles = int(row[0])
                                elif row[1] == 'dihedrals':
                                        n_dihedrals = int(row[0])
                        if len(row) == 3:
                                if row[2] == 'types':
                                        if row[1] == 'atom':
                                                n_atom_types = int(row[0])
                                        elif row[1] == 'bond':
                                                n_bond_types = int(row[0])
                                        elif row[1] == 'angle':
                                                n_angle_types = int(row[0])
                                        elif row[1] == 'dihedral':
                                                n_dihedral_types = int(row[0])
        return n_atoms, n_bonds, n_angles, n_dihedrals, n_atom_types, n_bond_types, n_angle_types, n_dihedral_types
