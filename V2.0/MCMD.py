import mof_taxi
import os
import subprocess

cwd = os.getcwd()
temperature = 303
pressure = 100.0   

adsorbate = "butane"
f_prefix = "ANUGIA"
f_lmps = 'data.'+f_prefix
f_cif = f_prefix+'.cif'
fa, ft, at = mof_taxi.get_parameters(f_lmps=f_lmps, adsorbate=adsorbate)
r_factor = 1000
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
    mof_taxi.updateRaspaRestart(f_lmps=cwd+'/Archive/'+str(start_cycle-1)+'_aftermd.lmps',num_framework_atoms=fa,num_framework_types=ft,num_adsorbate_types=at,cyclenumber = start_cycle)
    mof_taxi.updatecif(f_lmps=cwd+'/Archive/'+str(start_cycle-1)+'_aftermd.lmps',f_cif=f_cif, num_framework_atoms=fa,cyclenumber=start_cycle)
else:
    os.system('mkdir Archive')
    mof_taxi.write_raspa_input_voidHe(temperature,f_prefix)
    mof_taxi.run_raspa()
    voidFrac=float(subprocess.check_output("echo `grep 'Average Widom Rosenbluth-weight' Output/System_0/output_* | awk '{print $5}'`",shell=True))
    os.system('mv Output HeResult')
    mof_taxi.write_raspa_input_rosenbluth(temperature,adsorbate)
    mof_taxi.run_raspa()
    rosenbluth=float(subprocess.check_output("echo `grep 'Average Widom Rosenbluth-weight' Output/System_0/output_* | awk '{print $5}'`",shell=True))
    os.system('mv Output RosenbluthResult')
    start_cycle = 0
    mof_taxi.pairInfoGenandDelete(f_prefix)
	

while r_factor > 1.5:
    i = start_cycle
    mof_taxi.write_raspa_input(temp=temperature,pressure=pressure,cyclenumber=i,cycles=100,adsorbate=adsorbate,framework=f_prefix,rosenbluth=rosenbluth,HeliumVoidFraction=voidFrac)
    os.system('rm -r Restart')
    mof_taxi.run_raspa()
    mof_taxi.combine_lmps_RaspaRestart(framework=f_prefix, adsorbate=adsorbate,cyclenumber=i)
    mof_taxi.write_lammps_input(temp=temperature, pressure=pressure, md_time=500000, cyclenumber=i, mode="UA",framework=f_prefix)
    mof_taxi.run_lammps(i,f_prefix)
    mof_taxi.updatecif(f_lmps=str(i)+'.lmps',f_cif=f_cif, num_framework_atoms=fa,cyclenumber=i)
    mof_taxi.updateRaspaRestart(f_lmps=str(i)+'.lmps',num_framework_atoms=fa,num_framework_types=ft,num_adsorbate_types=at,cyclenumber = i)
    mof_taxi.archive(cyclenumber=i)
    if i == 0:
        os.system('mv Restart RestartInitial')
    mof_taxi.averages()
    r_factor = mof_taxi.r_factor(cyclenumber=i)
    start_cycle += 1
