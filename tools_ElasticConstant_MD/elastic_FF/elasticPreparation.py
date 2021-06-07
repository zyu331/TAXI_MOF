#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 15:38:22 2021

@author: jace
"""

import os 

MOFName = 'relaxed'

f = open('data.'+MOFName)
dataLine=f.readlines()
f.close()

f_pair=open("data.pair", "w")
f_data=open("data.MOF","w")
f_lj=open("lammps_lj","w")

keywords = ['Bond Coeffs\n','Atoms\n']
keyword_head=['Bond Coeffs\n','Angle Coeffs\n','Dihedral Coeffs\n','Improper Coeffs\n']
newKey = {'Bond Coeffs\n':'bond_coeff','Angle Coeffs\n':'angle_coeff','Dihedral Coeffs\n':'dihedral_coeff','Improper Coeffs\n':'improper_coeff','Pair Coeffs\n':'pair_coeff'}
out2pair=False
writeljPair=False
for x in dataLine:
    if any( [x ==element for element in keywords]):
        out2pair = not out2pair
        
    if out2pair:
        if any([x ==element for element in keyword_head]):
            currentHead = newKey[x]
            continue
        elif x == 'Pair Coeffs\n':
            currentHead = newKey[x]
            writeljPair = True
            continue
        elif len(x.split()) > 3:
            if writeljPair==True:
                f_lj.write(x)
                continue
            else:
                xNew = currentHead+' '+x
                f_pair.write(xNew)
    else:
        f_data.write(x)

f_pair.close()
f_data.close()
f_lj.close()

os.system('python pairDataGen.py')


f = open('in.'+MOFName)
f2 = open('init.mod','w')
f3 = open('potential.mod','w')

f2.write("""
# real units, elastic constants in GPa
variable up equal 1.0e-2
variable atomjiggle equal 1.0e-5
units		real
variable cfac equal 1.01325e-4
variable cunits string GPa
atom_style      full

# Define minimization parameters
variable etol equal 1.0e-10 
variable ftol equal 1.0e-10
variable maxiter equal 100000
variable maxeval equal 10000000
variable dmax equal 1.0e-3\n""")
f3.write("pair_style      lj/cut 12.5\n")

keywords = ['bond_style','angle_style','dihedral_style','improper_style']
lines = f.readlines()
for x in lines:
    try:
        if any([x.split()[0]==element for element in keywords]):
            f2.write(x)
            f3.write(x)
    except:
        continue
f2.write('box tilt large\n read_data       data.MOF\n')
f2.write('change_box all triclinic')
f3.write("""\n
special_bonds   lj 0.0 0.0 1.0
pair_modify     tail yes mix arithmetic
include data.pair

# Setup neighbor style
neighbor 1.0 nsq
neigh_modify once no every 1 delay 0 check yes

# Setup minimization style
min_style	     cg
min_modify	     dmax ${dmax} line quadratic

# Setup output
thermo		1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no\n""")

f.close()
f2.close()
f3.close()
os.system('rm lammps_lj')
