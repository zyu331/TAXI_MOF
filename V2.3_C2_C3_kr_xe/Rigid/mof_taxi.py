import os
from subprocess import call, Popen, PIPE
from operator import itemgetter
#from Queue import Queue, Empty
from threading import Thread
import numpy as np
from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove
import subprocess


LAMMPS_EXEC = os.environ.get('LAMMPS_EXEC')
RASPA_EXEC = os.environ.get('RASPA_DIR')+'bin/simulate'
RASPA_DIR = os.environ.get('RASPA_DIR')
TAXI_DIR = os.environ.get('TAXI_DIR')

"""
Created on Mon Feb  1 20:12:21 2021

@author: jace Yu

The folloing part is tested and mostly developed by Jace Yu
"""


def replace(file_path,pattern, subst ):
    #Create temp file
    fh, abs_path = mkstemp()
    
    substLine = pattern.split()[0] + ' ' + subst
    for i in range(1,len(pattern.split())):
        substLine = substLine +' '+ pattern.split()[i]
    substLine=substLine+'\n'
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, substLine))
    #Copy the file permissions from the old file to the new file
    copymode(file_path, abs_path)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)


def pairInfoGenandDelete(framework):
    file1 = open('force_field_mixing_rules.def', 'r')
    Lines = file1.readlines()
    file1.close()
    data=Lines[7:7+int(Lines[5])]
    epsilon=np.array([ x.split()[2] for x in data],dtype=float)*0.001987191686485529
    sigma = np.array([ x.split()[3] for x in data],dtype=float)

    print('please check VDW order is consistent between force_field_mixing_rules.def and .lmps')
    totalAtomType=int(Lines[5])-1
    f = open('data.pair','w')
    for i in range(0,totalAtomType):
        for j in range(i,totalAtomType):
            f.write('pair_coeff '+str(i+1)+' '+str(j+1)+' lj/cut ' + str( np.sqrt(epsilon[i]*epsilon[j])) + ' ' + str( (sigma[i]+sigma[j])/2 )+'\n')
    
    file2 = open('data.'+framework, 'r')
    Lines = file2.readlines()
    file2.close()
    
    for i in range(int(Lines[8].split()[0])+1,totalAtomType+1):
        f.write('pair_coeff '+'* '+str(i)+' coul/long \n')
       
    file3 = open('in.'+framework)
    Lines3 = file3.readlines()
    file3.close()
    
    
    file4=open('data.'+framework,'r') 
    Lines4 = file4.readlines()
    file4.close()
    keywordList=['bond_style','angle_style','dihedral_style','improper_style']
    keywordList_datafile={'bond_style':'Bond Coeffs\n','angle_style':'Angle Coeffs\n','dihedral_style':'Dihedral Coeffs\n','improper_style':'Improper Coeffs\n'}
    
    j=0
    for x in Lines3:
        x = x.split()
        try:
            if x[0] == keywordList[j] and x[1]!= 'hybrid':
                j+=1
                for m,y in enumerate(Lines4):
                    if Lines4[m-1] == keywordList_datafile[x[0]] and y =='\n':
                        newlineCount,newm=0,m
                        while newlineCount < 2:
                            if Lines4[newm]=='\n':
                                newlineCount+=1
                            else:
                                replace('data.'+framework,Lines4[newm],x[1])
                            newm+=1
        except:
            pass 
    
    file4=open('data.'+framework,'r') 
    Lines4 = file4.readlines()
    file4.close()
    file5=open('New_data.'+framework,'w+')  
    skipMarker=False
    for i,x in enumerate(Lines4):
     if x == 'Pair Coeffs\n':
         skipMarker=True
     elif x == 'Atoms\n':
         skipMarker=False     
         
     if skipMarker==False:
         file5.write(x)
         
    os.system('rm '+ 'data.'+framework)
    os.system('mv New_data.'+framework+' data.'+framework)
         
    return

def GenAdsorbateData(framework = None, adsorbate = None):
    
    ############# read file needed ##############
    file_framework = open('data.'+framework, 'r')
    Lines = file_framework.readlines()
    keywordDict={}
    for i in range(0,len(Lines)):
        if Lines[i-1]=='\n' and Lines[i+1] =='\n' and len(Lines[i].split())<4:
            keywordDict[Lines[i]]=i
    keywordDict['end']=len(Lines)+1
        
    file_adsorbate = open('data.'+adsorbate, 'r')
    Lines_ad = file_adsorbate.readlines()
    keywordDict_ad={}
    for i in range(0,len(Lines_ad)):
        if Lines_ad[i-1]=='\n' and Lines_ad[i+1] =='\n' and len(Lines_ad[i].split())<4:
            keywordDict_ad[Lines_ad[i]]=i    
    keywordDict_ad['end']=len(Lines_ad)+1
    
    raspa_file = os.listdir('Restart/System_0/')[0]
    file_restart = open('Restart/System_0/'+raspa_file,'r')
    Lines_restart=file_restart.readlines()
    
    for mm in range(0,len(Lines_restart)):
        line=Lines_restart[mm].split()
        try:
            if line[0]=='Component:' and line[2]=='Adsorbate':
                NUMofad=int(line[3])
        except:
            pass 
    
    file_framework.close()
    file_adsorbate.close()
    file_restart.close()
    
    ############## output file #################
    
    ### incorprate restart into adsrobate.lmps
    def adsorbate_write(NUM_adsorbate,restart_start):
        desiredOrder=['Atoms\n', 'Bonds\n', 'Angles\n', 'Dihedrals\n', 'Impropers\n','end']
        for name in desiredOrder:
            try:
                keywordDict_ad[name]
            except:
                desiredOrder.remove(name)
                
        for i in range(5,len(keyOutput)):
            if desiredOrder[i-5]!=keyOutput[i]:
                return("framwrok.lmps keyword order not correct, please change the order")
            adsorbate_lmp_output.write(Lines_ad[keywordDict_ad[desiredOrder[i-5]]])
            adsorbate_lmp_output.write('\n')
            for j in range(0,NUM_adsorbate):
                atomCount=0
                count_ad = keywordDict_ad[desiredOrder[i-5]]+2
                try:
                    while count_ad < keywordDict_ad[desiredOrder[i-4]]-1:
                        line_list=Lines_ad[count_ad].split()
                        line_list[0]= int( line_list[0]) + j*int(Lines_ad[2+i-5][0])
                        newline=''
                        if i==5:
                            for k in range(0,4):
                                newline=newline+str(line_list[k])+' '
                            for k in range(0,3):
                                newline=newline+Lines_restart[restart_start+j*int(Lines_ad[2][0])+atomCount].split()[k+3]+' '
                            atomCount+=1
                            newline=newline+'\n'
                            adsorbate_lmp_output.write(newline)
                            count_ad+=1
                        else:
                            for m in range(2,len(line_list)):
                                line_list[m]=int(line_list[m])+j*int(Lines_ad[2][0])
                            for x in line_list:
                                newline=newline+str(x)+' '
                            newline=newline+'\n'
                            adsorbate_lmp_output.write(newline)
                            count_ad+=1
                except:
                    pass
                                
            adsorbate_lmp_output.write('\n')

    keyOutput=sorted(keywordDict, key=keywordDict.get)
    keyOutput_ad=sorted(keywordDict_ad, key=keywordDict_ad.get)
    desiredOrder=['Masses\n','Bond Coeffs\n', 'Angle Coeffs\n', 'Dihedral Coeffs\n', 'Improper Coeffs\n']
    count_framework=0
    offset={}
    adsorbate_lmp_output=open('data.adsorbate','w+')
    
    for i in range(0,5): # default to improper!
    
    ### record offset
        if desiredOrder[i]!=keyOutput[i]:
            return "framwrok.lmps keyword order not correct, please change the order"
        
        while count_framework < keywordDict[keyOutput[i+1]]-1:
            try:
                NUM= Lines[count_framework].split()[0]
            except:
                pass
            count_framework+=1
        offset[keyOutput[i]]=NUM
        
    ### output head of adsorbate.lmps
        if i==0:
            count_ad=0
            while count_ad < keywordDict_ad[keyOutput_ad[i+1]]:
                if len(Lines_ad[count_ad].split()) == 2:
                    line=Lines_ad[count_ad].split()
                    newLine=str(int(line[0])*NUMofad)+' '+line[1]+'\n'
                else:
                    newLine=Lines_ad[count_ad]
                adsorbate_lmp_output.write(newLine)
                count_ad+=1
     ### output coeff part of the adsorbate.lmps
        elif keyOutput_ad[i]!='Pair Coeffs\n':
            try:
                count_ad = keywordDict_ad[keyOutput_ad[i]]
                while count_ad < keywordDict_ad[keyOutput_ad[i+1]]-1:
                    try:
                        adsorbate_lmp_output.write(Lines_ad[count_ad])
                        count_ad+=1
                    except:
                        pass     
                adsorbate_lmp_output.write('\n')
            except:
                continue 
        
    ### output pos/angle/dihedral/im part of the adsorbate.lmps        
    for j in range(0,len(Lines_restart)):
        line=Lines_restart[j].split()
        try:
            if line[0]=='Component:' and line[2]=='Adsorbate':
                adsorbate_write(int(line[3]),j+2)
        except:
            continue

    adsorbate_lmp_output.close()    
    
    return NUMofad


def preRASPA(temperature,f_prefix,adsorbate,superCell):
    write_raspa_input_voidHe(temperature,f_prefix,superCell)
    run_raspa()
    voidFrac=float(subprocess.check_output("echo `grep 'Average Widom Rosenbluth-weight' Output/System_0/output_* | awk '{print $5}'`",shell=True))
    os.system('mv Output HeResult')
    write_raspa_input_rosenbluth(temperature,adsorbate)
    run_raspa()
    rosenbluth=float(subprocess.check_output("echo `grep 'Average Widom Rosenbluth-weight' Output/System_0/output_* | awk '{print $5}'`",shell=True))
    os.system('mv Output RosenbluthResult')   
    
    return voidFrac,rosenbluth


def write_raspa_input_voidHe(temp=None,framework=None,superCell=None):
	if temp is None:
		print('you must specify a temperature')
		return
	if framework is None:
		print('you must specify the name of the framework material')
		return
	f = open('simulation.input', 'w')
	f.write('SimulationType                MonteCarlo\n')
	f.write('NumberOfCycles 10000\n')
	f.write('PrintEvery                    1000\n')
	f.write('PrintPropertiesEvery                    1000\n')

	f.write('CutOff                           12.5\n')
	f.write('UseChargesFromCIFFile         yes\n\n')
	
	f.write('Framework                     0\n')
	f.write('FrameworkName                 '+framework+'\n')
	f.write('UnitCells                     '+" ".join(superCell)+'\n')
	f.write('ExternalTemperature           '+str(temp)+'\n')
	f.write('RemoveAtomNumberCodeFromLabel yes \n\n')

	f.write('Component 0 MoleculeName helium\n')
	f.write('MoleculeDefinition TraPPE\n')
	f.write('WidomProbability          1.0\n')
	f.write('CreateNumberOfMolecules    0\n')
	f.close()
	return		

def write_raspa_input_rosenbluth(temp=None,adsorbate=None):
	if temp is None:
		print('you must specify a temperature')
		return
	if adsorbate is None:
		print('you must specify the name of the framework material')
		return
	f = open('simulation.input', 'w')
	f.write('SimulationType                MonteCarlo\n')
	f.write('NumberOfCycles 10000\n')
	f.write('PrintEvery                    1000\n')
	f.write('PrintPropertiesEvery                    1000\n')

	f.write('CutOff                           12.5\n')
	f.write('UseChargesFromCIFFile         yes\n\n')
	
	f.write('Box 0\n')
	f.write('BoxLengths 30 30 30\n')
	f.write('ExternalTemperature           '+str(temp)+'\n')

	f.write('Component 0 MoleculeName ' + adsorbate + '\n')
	f.write('WidomProbability          1.0\n')
	f.write('CreateNumberOfMolecules    0\n')
	f.close()
	return	

def write_raspa_input(temp=None, pressure=None, cyclenumber=None, cycles=None,framework=None, adsorbate=None, rosenbluth=None,HeliumVoidFraction=None,startNUM=None,superCell=None):
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
    f.write('PrintEvery                    1000\n')
    f.write('RestartFile                   no\n\n')

    f.write('CutOff                           12.5\n')
    f.write('EwaldPrecision                   1e-6\n')
    f.write('UseChargesFromCIFFile         yes\n\n')
    f.write('Framework                     0\n')
    f.write('FrameworkName                 '+framework+'\n')
    f.write('UnitCells                     '+" ".join(superCell)+'\n')
    f.write('HeliumVoidFraction            '+str(HeliumVoidFraction)+'\n')
    f.write('ExternalTemperature           '+str(temp)+'\n')
    f.write('ExternalPressure              '+str(pressure)+'\n\n')
    f.write('Movies                        yes\n')
    f.write('WriteMoviesEvery	1000\n\n')
    f.write('RemoveAtomNumberCodeFromLabel yes \n\n')
    
    f.write('Component 0 MoleculeName             '+adsorbate+'\n')
    f.write('IdealGasRosenbluthWeight ' + str(rosenbluth) + '\n')
    f.write('TranslationProbability   1.0\n')
    f.write('RotationProbability      1.0\n')
    f.write('ReinsertionProbability   1.0\n')
    f.write('SwapProbability          1.0\n')
    if cyclenumber ==0:
        f.write('CreateNumberOfMolecules    0\n')
    else:
        f.write('CreateNumberOfMolecules    '+str(startNUM)+'\n')
            
    f.close()
    return

def write_lammps_input(temp=None, pressure=None, md_time=None, cyclenumber=None,framework=None,adsorbate_Types=None,f_Types=None,elementArray=None):
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

    '''
	convert pressure to atm
	'''
    pressure*=9.86923e-6
    f=open('npt.in', 'w')	
    f.write("""units           real
atom_style      full
boundary        p p p

atom_style full\n
pair_style hybrid/overlay lj/cut 12.5 coul/long 12.5
pair_modify tail yes mix arithmetic
kspace_style	ewald 1.0E-06
special_bonds   lj/coul 0.0 0.0 1.0
bond_style hybrid harmonic\nangle_style hybrid cosine/periodic harmonic
dihedral_style  hybrid harmonic opls\nimproper_style  hybrid fourier
dielectric      1.0
box tilt        large
neighbor 1.0 bin
neigh_modify once no every 1 delay 0 check yes

read_data data."""+framework+"""\n
include   	 data.pair\n""")
    
    f_temp=open('in.'+framework)
    for line in f_temp:
        row=line.split()
        try:
            if row[0]=='group':
                f.write(line)
                break
        except:
            pass
        
    f.write('\n')
    f_temp.close()        
    
    f.write("""
############################ Record thermodynamic data #################################################
            thermo 		10000 
            thermo_style custom step etotal evdwl ecoul vol
            compute 	t_temp all temp
            change_box all triclinic""")
    f.write("""            
############################### minimization ###########################################################
            dump opt_frame all xyz 10000 fig/opt_framework_*.xyz
            dump_modify	opt_frame element """+ " ".join(elementArray)+"""
            
            label		loop
            variable	a loop 200    
            print		"------> GEO OPT Run $a"
            min_style	cg
            min_modify	line quadratic dmax 0.001
            minimize	1e-10 1e-10 10000 1000000           
            print		"------> CELL OPT Run $a"
            fix		cell all box/relax tri 0.0 vmax 0.001
            min_style	cg
            min_modify	line quadratic dmax 0.001
            minimize	1e-10 1e-10 10000 1000000
            unfix		cell
            next		a
            jump npt.in loop
            label break
            undump opt_frame""")
    f.write("""
################################# pre-run setting #####################################################
            timestep	1
            run_style	verlet
            fix		TETHER1 all momentum 1 linear 1 1 1 angular
            fix		COM fram recenter INIT INIT INIT
            velocity	all create """+str(temp)+ """ 123456 rot yes mom yes dist gaussian
############################### NVT ###########################################################
            dump nvt_struc all xyz 10000 fig/nvt_framework_*.xyz
            dump_modify nvt_struc element """+" ".join(elementArray)+"""
            reset_timestep 0
            fix NVT_ENSEMBLE all nvt temp """+str(temp)+' '+str(temp)+""" 100
            run              """+str(md_time)+"""
            unfix NVT_ENSEMBLE
            undump nvt_struc\n\n""")

# ############################### NPT ###########################################################
#             dump npt_struc all xyz 100000 fig/npt_framework_*.xyz
#             dump_modify	npt_struc element """+" ".join(elementArray)+"""
#             fix NPT_ENSEMBLE all npt temp """+str(temp)+' '+str(temp)+' 100 tri '+str(pressure)+' '+str(pressure)+' 100'+"""
#             timestep 1.0
#             run              """+str(md_time)+"""
#             unfix         NPT_ENSEMBLE
#             unfix         COM
    f.write("""write_data       """ +str(cyclenumber)+""".lmps""")	
    f.close()
    return

def run_raspa():
    global RASPA_EXEC
    p = Popen([RASPA_EXEC, 'simulation.input'], stdin=PIPE, stdout=PIPE, stderr=PIPE).communicate()
    print(p)
    del(p)
    return
		
def run_lammps(cyclenumber,framework,prefix=None, np=None):
    if prefix is None:
        print('Defaulting to the srun --mpi=pmix_v2 prefix(changed by yzz)')
        prefix=['mpirun', '-np','16']
    global LAMMPS_EXEC
    f = open('npt.in', 'rb').read()
    p = Popen(prefix + [LAMMPS_EXEC,'-var','fname','npt.in','-e','both'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdo, stde = p.communicate(f)
    print(p)
    del(p)
	
    return

"""
@author: jace Yu /Anstine

The folloing part is writen by Anstine, and modified by Jace from python2 to python3 without further testing

"""

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
    
    def _cellParameterInput2_(cell_info):
        xy,xz,yz=cell_info[0,2],cell_info[1,2],cell_info[2,2]          
        lx,ly,lz=cell_info[0,1]-cell_info[0,0],cell_info[1,1]-cell_info[1,0],cell_info[2,1]-cell_info[2,0]        
        a,b,c=lx,np.sqrt(ly**2+xy**2),np.sqrt(lz**2+xz**2+yz**2)
        
        alpha,beta,gamma=(xy*xz+ly*yz)/(b*c),xz/c,xy/b
        box=np.array([a,b,c,np.arccos(alpha)*180/np.pi,np.arccos(beta)*180/np.pi,np.arccos(gamma)*180/np.pi])
        omega=a*b*c*np.sqrt(1-alpha**2-beta**2-gamma**2+2*alpha*beta*gamma)            
        transM=np.array(([a,b*gamma,c*beta],[0,b*np.sin(np.arccos(gamma)), c*(alpha-beta*gamma)/np.sin(np.arccos(gamma))],[0,0,omega/(a*b*np.sin(np.arccos(gamma)))]))
        
        return box,transM
    
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
        elif len(row)==6 and row[3] == 'xy':
            xy,xz,yz=row[0],row[1],row[2]
            cellInfo,transM = _cellParameterInput2_(np.array([[xlo,xhi,xy],[ylo,yhi,xz],[zlo,zhi,yz]],dtype=float))
        elif len(row) > 0 and row[0] == 'Atoms':
                start_parsing = True
        if start_parsing == True and len(row) > 3 and int(row[0]) <= num_framework_atoms:
            pos_old=np.array([float(row[4]),float(row[4]),float(row[4])])
            pos = np.matmul(np.linalg.inv(transM),pos_old)
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
            fl.write('_cell_length_a                         ' + str(cellInfo[0]) + '\n')
            continue
        elif row[0] == '_cell_length_b':
            fl.write('_cell_length_b                         ' + str(cellInfo[1]) + '\n')
            continue
        elif row[0] == '_cell_length_c':
            fl.write('_cell_length_c                         ' + str(cellInfo[2]) + '\n')
            continue
        elif row[0] == '_cell_angle_alpha':
            fl.write('_cell_angle_alpha  '+str(cellInfo[3])+ '\n')
            continue
        elif row[0] == '_cell_angle_beta':
            fl.write('_cell_angle_beta  '+str(cellInfo[4])+ '\n')
            continue
        elif row[0] == '_cell_angle_gamma':
            fl.write('_cell_angle_gamma  '+str(cellInfo[5])+ '\n')
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
    f_atomNUM, f_Types,adsorbate_atomNUM,adsorbate_Types = [0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]
    
    f = open(f_lmps, 'r')
    for line in f:
        row = line.split()
        if not line.strip():
            continue
        if len(row) > 1:
            if row[1] == 'atoms':
                f_atomNUM[0] = int(row[0])
            if row[1] == 'atom':
                f_Types[0] = int(row[0])
            if row[1] == 'bonds':
                f_atomNUM[1] = int(row[0])
            if row[1] == 'bond':
                f_Types[1] = int(row[0])
            if row[1] == 'angles':
                f_atomNUM[2] = int(row[0])
            if row[1] == 'angle':
                f_Types[2] = int(row[0])
            if row[1] == 'dihedrals':
                f_atomNUM[3] = int(row[0])
            if row[1] == 'dihedral':
                f_Types[3] = int(row[0])
            if row[1] == 'impropers':
                f_atomNUM[4] = int(row[0])
            if row[1] == 'improper':
                f_Types[4] = int(row[0])

    f.close()
    f = open('data.'+adsorbate, 'r')
    for line in f:
        row = line.split()
        if not line.strip():
            continue
        if len(row) > 1:
            if row[1] == 'atoms':
                adsorbate_atomNUM[0] = int(row[0])
            if row[1] == 'atom':
                adsorbate_Types[0] = int(row[0])
            if row[1] == 'bonds':
                adsorbate_atomNUM[1] = int(row[0])
            if row[1] == 'bond':
                adsorbate_Types[1] = int(row[0])
            if row[1] == 'angles':
                adsorbate_atomNUM[2] = int(row[0])
            if row[1] == 'angle':
                adsorbate_Types[2] = int(row[0])
            if row[1] == 'dihedrals':
                adsorbate_atomNUM[3] = int(row[0])
            if row[1] == 'dihedral':
                adsorbate_Types[3] = int(row[0])
            if row[1] == 'impropers':
                adsorbate_atomNUM[4] = int(row[0])
            if row[1] == 'improper':
                adsorbate_Types[4] = int(row[0])
                
    return f_atomNUM, f_Types, adsorbate_atomNUM, adsorbate_Types


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
    f.write('\ttaxi.write_force_field()\n')
    f.write("\tstart_cycle = 0\n\n")
    f.write('while r_factor > 1.5:\n')
    f.write('\ti = start_cycle\n')
    f.write('\ttaxi.write_raspa_input(temp=temperature,pressure=pressure,cyclenumber=i,cycles=50,adsorbate=adsorbate,framework=f_prefix,rosenbluth='+str(rosenbluth)+')\n')
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


def updateFrameworkData(cyclenumber,framework):
    lmpsName=str(cyclenumber)+'.lmps'
    f = open(lmpsName)
    count =0
    for line in f:
        if not line.strip():
            count += 1
            continue
        count += 1
        row = line.split()
        if len(row)>3:
            if row[2]=='xlo':
                propertyEnd=count
                continue
        if row[0]=='Masses':
            result_begin=count
        if row[0]=='Atoms':
            result_end=count
    f.close()
    
    
    f = open('data.'+framework)
    count2 =0
    for line in f:
        if not line.strip():
            count2 += 1
            continue
        count2 += 1
        row = line.split()
        if len(row)>3:
            if row[2]=='xlo':
                propertyEnd2=count2
        if row[0]=='Masses':
            origin_begin=count2
        if row[0]=='Atoms':
            origin_end=count2
            break	
    f.close()

    os.system('head -n '+str(origin_end-1)+' data.'+framework+'| tail -n '+str(origin_begin-origin_end-1)+' > output_property')
    os.system('head -n'+str(propertyEnd2-1)+ ' data.'+framework+' > output_head' )
    os.system('head -n '+str(result_begin-2)+' '+ str(cyclenumber)+'.lmps' + ' | tail -n ' +str(result_begin-propertyEnd)+ ' > output_cell')
    os.system('tail -n '+str(count-result_end+1)+' '+str(cyclenumber)+ '.lmps > output_pos')
    os.system('cat output_head outpt_cell > output_temp')
    os.system('cat output_temp output_property > output_temp2')
    os.system('cat output_temp2 output_pos > output_final')    
    os.system('rm '+ 'data.'+framework)
    os.system('mv output_final data.'+framework)
    os.system('rm output_*')

