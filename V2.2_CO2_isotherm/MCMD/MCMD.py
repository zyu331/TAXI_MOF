import mof_taxi
import os
import subprocess

cwd = os.getcwd()
temperature = 300
pressure = myP   
adsorbate = "CO2"
f_prefix = "XXXXXX"
r_factor = 1000
dumpModifyArray=['Zn','Zn','H','C','C','O','O','O','C']

f_lmps = 'data.'+f_prefix
f_cif = f_prefix+'.cif'
f_atomNUM, f_Types, adsorbate_atomNUM, adsorbate_Types = mof_taxi.get_parameters(f_lmps=f_lmps, adsorbate=adsorbate)

if os.path.isdir('Archive'):
    lst=os.listdir(cwd+'/Archive/')
    start_cycle = mof_taxi.get_cycle()
    lst = os.listdir(cwd+'/Restart/System_0/')
    try:
        lst.remove('temprestart.txt')
    except:
        pass
    start_cycle -= 1
    os.system('cp ' + cwd+'/Archive/'+str(start_cycle)+'.restart '+cwd+'/Restart/System_0/'+lst[0])
    os.system('cp ' + cwd+'/Archive/'+str(start_cycle)+'.cif '+cwd+'/'+f_cif)
    mof_taxi.updatecif(f_lmps=cwd+'/Archive/'+str(start_cycle-1)+'_aftermd.lmps',f_cif=f_cif, num_framework_atoms=f_atomNUM[0],cyclenumber=start_cycle)
else:
    os.system('mkdir Archive')
    voidFrac,rosenbluth = mof_taxi.preRASPA(temperature,f_prefix,adsorbate)
    start_cycle = 0
    mof_taxi.pairInfoGenandDelete(f_prefix)

NUMofad = 0
while r_factor > 1.5:
########################Main MC MD loop###########################
    i = start_cycle
    mof_taxi.write_raspa_input(temp=temperature,pressure=pressure,cyclenumber=i,cycles=10000,adsorbate=adsorbate,framework=f_prefix,rosenbluth=rosenbluth,HeliumVoidFraction=voidFrac,startNUM= NUMofad )
    os.system('rm -r Restart Movies Output VTK')
    mof_taxi.run_raspa()
    os.system('rm data.adsorbate')
    NUMofad=mof_taxi.GenAdsorbateData(framework=f_prefix, adsorbate=adsorbate,cyclenumber=i)
    mof_taxi.write_lammps_input(temp=temperature, pressure=pressure, md_time=500000, cyclenumber=i, framework=f_prefix,adsorbate_Types=adsorbate_Types,f_Types=f_Types,elementArray=dumpModifyArray)
    mof_taxi.run_lammps(i,f_prefix)
    mof_taxi.updatecif(f_lmps=str(i)+'.lmps',f_cif=f_cif, num_framework_atoms=f_atomNUM[0],cyclenumber=i)
    mof_taxi.updateFrameworkData(cyclenumber=i,framework=f_prefix)
########################Statistics###########################
    mof_taxi.archive(cyclenumber=i)
    mof_taxi.averages()
    r_factor = mof_taxi.r_factor(cyclenumber=i)
    start_cycle += 1
