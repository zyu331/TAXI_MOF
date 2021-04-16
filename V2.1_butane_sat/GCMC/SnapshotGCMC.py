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

"""
Created on Mon Feb  1 20:12:21 2021

@author: jace Yu

The folloing part is tested and mostly developed by Jace Yu
"""




def preRASPA(temperature,f_prefix,adsorbate):
    write_raspa_input_voidHe(temperature,f_prefix)
    run_raspa()
    voidFrac=float(subprocess.check_output("echo `grep 'Average Widom Rosenbluth-weight' Output/System_0/output_* | awk '{print $5}'`",shell=True))
    os.system('mv Output HeResult')
    write_raspa_input_rosenbluth(temperature,adsorbate)
    run_raspa()
    rosenbluth=float(subprocess.check_output("echo `grep 'Average Widom Rosenbluth-weight' Output/System_0/output_* | awk '{print $5}'`",shell=True))
    os.system('mv Output RosenbluthResult')   
    
    return voidFrac,rosenbluth


def write_raspa_input_voidHe(temp=None,framework=None):
	if temp is None:
		print('you must specify a temperature')
		return
	if framework is None:
		print('you must specify the name of the framework material')
		return
	f = open('simulation.input', 'w')
	f.write('SimulationType                MonteCarlo\n')
	f.write('NumberOfCycles 1000\n')
	f.write('PrintEvery                    1000\n')
	f.write('PrintPropertiesEvery                    1000\n')

	f.write('CutOff                           12.5\n')
	f.write('UseChargesFromCIFFile         yes\n\n')
	
	f.write('Framework                     0\n')
	f.write('FrameworkName                 '+framework+'\n')
	f.write('UnitCells                     1 1 1\n')
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
	f.write('NumberOfCycles 1000\n')
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

def write_raspa_input_Henry(temp=None, pressure=None, cyclenumber=None, cycles=None,framework=None, adsorbate=None, rosenbluth=None,HeliumVoidFraction=None,startNUM=None,superCell=None):
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
    f = open('simulation.input_Henry', 'w')

    f.write('SimulationType                MonteCarlo\n')
    f.write('NumberOfCycles                ' + str(cycles)+'\n')
    f.write('NumberOfInitializationCycles  ' + str(cycles)+'\n')
    f.write('PrintEvery                    1000\n')
    f.write('RestartFile                   no\n\n')

    f.write('CutOff                           12.5\n')
    f.write('EwaldPrecision                   1e-6\n')
    f.write('UseChargesFromCIFFile         yes\n\n')
    f.write('Framework                     0\n')
    f.write('FrameworkName                 mya\n')
    f.write('UnitCells                     '+" ".join(superCell)+'\n')
    f.write('HeliumVoidFraction            '+str(HeliumVoidFraction)+'\n')
    f.write('ExternalTemperature           '+str(temp)+'\n')

    f.write('Movies                        yes\n')
    f.write('WriteMoviesEvery	1000\n\n')
    f.write('RemoveAtomNumberCodeFromLabel yes \n\n')
    
    f.write('Component 0 MoleculeName             '+adsorbate+'\n')
    f.write('IdealGasRosenbluthWeight ' + str(rosenbluth) + '\n')
    f.write('MoleculeDefinition        TraPPE\n')
    f.write('WidomProbability   1.0\n')
    if cyclenumber ==0:
        f.write('CreateNumberOfMolecules    0\n')
    else:
        f.write('CreateNumberOfMolecules    '+str(startNUM)+'\n')
            
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
    f.write('FrameworkName                 mya\n')
    f.write('UnitCells                     '+" ".join(superCell)+'\n')
    f.write('HeliumVoidFraction            '+str(HeliumVoidFraction)+'\n')
    f.write('ExternalTemperature           '+str(temp)+'\n')
    f.write('ExternalPressure              '+str(pressure)+'\n\n')
    f.write('Movies                        yes\n')
    f.write('WriteMoviesEvery	1000\n\n')
    f.write('RemoveAtomNumberCodeFromLabel yes \n\n')
    
    f.write('Component 0 MoleculeName             '+adsorbate+'\n')
    f.write('IdealGasRosenbluthWeight ' + str(rosenbluth) + '\n')
    f.write('MoleculeDefinition        TraPPE\n')
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


def run_raspa():
    RASPA_EXEC = os.environ.get('RASPA_DIR')+'bin/simulate'
    p = Popen([RASPA_EXEC, 'simulation.input'], stdin=PIPE, stdout=PIPE, stderr=PIPE).communicate()
    print(p)
    del(p)
    return



cwd = os.getcwd()
temperature = 300
pressure = 1000000.0   
adsorbate = "butane"
f_prefix = "AROFET"
r_factor = 1000
superCell=['1','1','1']
f_cif = f_prefix+'.cif'

os.system('cp ../cifs/0.cif '+f_cif)
voidFrac,rosenbluth = preRASPA(temperature,f_prefix,adsorbate)
write_raspa_input(temp=temperature,pressure=pressure,cyclenumber=0,cycles=100000,adsorbate=adsorbate,framework=f_prefix,rosenbluth=rosenbluth,HeliumVoidFraction=voidFrac,startNUM=0,superCell=superCell)
write_raspa_input_Henry(temp=temperature,pressure=pressure,cyclenumber=0,cycles=100000,adsorbate=adsorbate,framework=f_prefix,rosenbluth=rosenbluth,HeliumVoidFraction=voidFrac,startNUM=0,superCell=superCell)
os.system('./makeJobs.sh')
