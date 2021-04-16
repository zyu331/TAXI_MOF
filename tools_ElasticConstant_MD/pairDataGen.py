#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 11:01:35 2021

@author: jace
"""

import pandas as pd 
import numpy as np

data=pd.read_table('lammps_lj',delim_whitespace=True,header=None)


f = open("data.pairVDW", "w")


for i in range(0,len(data)):
    for j in range(i,len(data)):
        f.write('pair_coeff '+str(i+1)+' '+str(j+1)+' ' + str( np.sqrt(data.loc[i][1]*data.loc[j][1]) ) + ' ' + str( (data.loc[i][2]+data.loc[j][2])/2 )+'\n')
        
        
f.close()
