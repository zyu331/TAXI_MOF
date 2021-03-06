#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 19:02:00 2020

@author: zyu
"""
import sys
import numpy as np
import pandas as pd
import os

class DataInput(object):
        def __init__(self): #No need to implement
             pass
        
        def _read_lammps_data_(self,dataName,snapName):

            self.dataName=dataName
            self.snapName=snapName
            self._read_snapshots_LAMMPS_()
        # def _fetchAtomName_(self):
        #     with open('data.'+self.dataName, 'r') as the_file:
        #         all_data = [line.strip() for line in the_file.readlines()] 
        #     self.NameModify={}
        #     for i,x in enumerate(all_data):
        #         if x == 'Masses':
        #             readMark=True
        #             j=i
        #             while readMark==True:
        #                 line=all_data[j].split()
        #                 if len(line) == 4:
        #                     self.NameModify[line[0]]=line[3]
        #                 elif all_data[j]=='Bond Coeffs':
        #                     readMark=False
        #                 j+=1
        #             break
        
        def _read_snapshots_LAMMPS_(self):
            File='fig/'+self.snapName
            with open(File, 'r') as the_file:
                all_data = [line.strip() for line in the_file.readlines()]
                temp=np.array([element.split(' ', 6) for element in all_data[9:]],dtype=float)
                data=temp.view(np.ndarray)[np.lexsort((temp[:, 0], ))]
            
                self.cell_info=np.array([line.split() for line in np.array(all_data[5:8])],dtype=float)
                self.box,self.charge,self.position,self.atomType,index=self.cell_info[:,1],data[:,5],data[:,2:5],np.array(data[:,1],dtype=int),np.array(data[:,0],dtype=int)
                self.positionFrac=self.position[:]
                self._cellParameterInput2_()
                
            return (self.position,self.positionFrac,self.charge,self.atomType,self.lattice,index)

        def _cellParameterInput2_(self):
            xy,xz,yz=self.cell_info[0,2],self.cell_info[1,2],self.cell_info[2,2]          
            lx=(self.cell_info[0,1]-max((0,xy,xz,xy+xz)))-(self.cell_info[0,0]-min((0,xy,xz,xy+xz)))
            ly=(self.cell_info[1,1]-max(0,yz))-(self.cell_info[1,0]-min(0,yz))
            lz=self.cell_info[2,1]-self.cell_info[2,0]
            a,b,c=lx,np.sqrt(ly**2+xy**2),np.sqrt(lz**2+xz**2+yz**2)
            self.box=[a,b,c]
            alpha,beta,gamma=(xy*xz+ly*yz)/(b*c),xz/c,xy/b
            omega=a*b*c*np.sqrt(1-alpha**2-beta**2-gamma**2+2*alpha*beta*gamma)            
    # add PBC
            self.box_lammps=np.array([lx,ly,lx,np.arccos(alpha),np.arccos(beta),np.arccos(gamma)])
            transM=np.array(([a,b*gamma,c*beta],[0,b*np.sin(np.arccos(gamma)), c*(alpha-beta*gamma)/np.sin(np.arccos(gamma))],[0,0,omega/(a*b*np.sin(np.arccos(gamma)))]))
            self.lattice=transM

#for cif output
            self.cell_para_4cif=[a,b,c,np.rad2deg(np.arccos(alpha)),np.rad2deg(np.arccos(beta)),np.rad2deg(np.arccos(gamma))]
    
        def _out2cif_(self,fileFolder,dumpModifyArray):
            # mkdir
            label=[dumpModifyArray[x-1] for x in self.atomType]
            cellParameter=self.cell_para_4cif
            position=np.column_stack((self.positionFrac,self.charge))
            #if np.shape(position)[2]==3:
                #filename = fileFolder+'/'+'RigidOrSingleOut'+'.cif'
            #else:
            filename = fileFolder+'/'+self.dataName+'.cif'
                
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            all_info=np.column_stack((label,position))
        #           
            t=list(['_cell_length_a '+str(cellParameter[0])])
            t.append('_cell_length_b '+str(cellParameter[1]))
            t.append('_cell_length_c '+str(cellParameter[2]))
            t.append('_cell_angle_alpha '+str(cellParameter[3]))
            t.append('_cell_angle_beta '+str(cellParameter[4]))
            t.append('_cell_angle_gamma '+str(cellParameter[5]))
            t.append("_symmetry_space_group_name_H-M         'P 1'")
            t.append('_symmetry_Int_Tables_number            1')
            t.append('loop_')
            t.append('_symmetry_equiv_pos_as_xyz')
            t.append("   'x, y, z'")
            t.append('                         ')
            t.append('loop_')
            t.append('_atom_site_label') 
            t.append('_atom_site_fract_x ')
            t.append('_atom_site_fract_y')
            t.append('_atom_site_fract_z')
            t.append('_atom_site_charge')
                    
            position=np.array(position)
            for j in range (0,18):
                print(t[j],file = open(filename,"a"))
                
            for k in range(0,np.shape(position)[0]):
                print (*all_info[k,:], file = open(filename,"a")) 
                k
        
            return