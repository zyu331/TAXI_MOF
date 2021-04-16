import mof_taxi
import os
import subprocess
from hMOFcifGen import DataInput

cwd = os.getcwd()
temperature = 300
pressure = 1000000
adsorbate = 'myAd'
f_prefix = "myFramework"
r_factor = 1000
dumpModifyArray=myList  # see data.$framework atomType Section

f_lmps = 'data.'+f_prefix
f_cif = f_prefix+'.cif'
f_atomNUM, f_Types, adsorbate_atomNUM, adsorbate_Types = mof_taxi.get_parameters(f_lmps=f_lmps, adsorbate=adsorbate)

mof_taxi.pairInfoGenandDelete(f_prefix)
mof_taxi.write_lammps_input(temp=temperature, pressure=pressure, md_time=500000, cyclenumber=0, framework=f_prefix,elementArray=dumpModifyArray)
mof_taxi.run_lammps(0,f_prefix)

os.system('rm -r cifs')
os.system('mkdir cifs')
SnapshotList=['nvt_framework_400000','nvt_framework_410000','nvt_framework_420000','nvt_framework_430000','nvt_framework_440000','nvt_framework_450000','nvt_framework_460000','nvt_framework_470000','nvt_framework_480000','nvt_framework_490000']
for j,x in enumerate(SnapshotList):
    a=DataInput()
    a._read_lammps_data_(f_prefix,x)
    a._out2cif_('cifs',dumpModifyArray)
    os.system('mv cifs/'+f_cif+ ' cifs/'+str(j)+'.cif')
    
optList=os.listdir('fig')
maxNow=0
for x in optList:
    line=x.split('_')
    if line[0]=='opt' and int(line[-1]) > maxNow:
        maxNow=int(line[-1])
        strNow=x
        
b=DataInput()
b._read_lammps_data_(f_prefix,x)
b._out2cif_('cifs',dumpModifyArray)
os.system('mv cifs/'+f_cif+ ' cifs/'+'relaxed.cif')
