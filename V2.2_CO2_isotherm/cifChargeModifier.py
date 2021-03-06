#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 15:59:39 2021

@author: jace
"""
import os 
import sys

MOFName = sys.argv[1]


f= open(MOFName+'.cif','r')
Lines=f.readlines()
f.close()

f_out=open(MOFName+'_new.cif','a')
flag=False

for i,x in enumerate(Lines):
    if x=='_atom_site_charge\n':
        f_out.write(x)
        break
    else:
        f_out.write(x)
     
charge=0
for j in range(i+1,len(Lines)):
    element=Lines[j].split()
    if len(element)>0:
        charge+=float(element[-1])

chargeModify=charge/(len(Lines)-i)
for j in range(i+1,len(Lines)):
    element=Lines[j].split()
    outStr=str(element[0])+' '+str(element[1])+' '+str(element[2])+' '+str(element[3])+' '+str(element[4])+' '+str(float(element[5])-chargeModify)+'\n'
    f_out.write(outStr)
    
f_out.close()
    
os.system('rm '+MOFName+'.cif')
os.system('mv '+MOFName+'_new.cif'+ ' '+MOFName+'.cif')        
os.system('lammps-interface -ff UFF4MOF '+MOFName+'.cif')


