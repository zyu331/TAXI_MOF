import mof_taxi
import os
import subprocess

cwd = os.getcwd()
temperature = 300
pressure = 1000000.0   
adsorbate = "butane"
f_prefix = "XXXXXX"
r_factor = 1000
superCell=['X','X','X']

f_lmps = 'data.'+f_prefix
f_cif = f_prefix+'.cif'


voidFrac,rosenbluth = mof_taxi.preRASPA(temperature,f_prefix,adsorbate,superCell)
mof_taxi.write_raspa_input(temp=temperature,pressure=pressure,cyclenumber=0,cycles=100000,adsorbate=adsorbate,framework=f_prefix,rosenbluth=rosenbluth,HeliumVoidFraction=voidFrac,startNUM=0,superCell=superCell)
mof_taxi.run_raspa()

