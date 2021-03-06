import pandas as pd
import os


adsorbate=['ethane','ethene','propane','propene','xe','kr']
adsorbate_atom_dict={'ethane':['CH3_sp3'], 'ethene':['CH2_sp2'], 'propane':['CH3_sp3','CH2_sp3'], 'propene':['CH2_sp2','CH_sp2','CH3_sp3'],'xe':['xe'],'kr':['kr']}

framework=os.getcwd().split('/')[-3]
os.system('sed -i "s,myFramework,"'+framework+'",g" MCMD.py' )
os.system('sed -i "s,myFramework,"'+framework+'",g" submit.sh')

#### read atom type from .data ####
atomType=[]
record=False
f=open('data.'+framework,'r')
for x in f:
    if not x.strip():
        continue
    element=x.split()
    if x=='Bond Coeffs\n':
        break
    if x=='Masses\n':
        record=True
        continue
    if record==True:
        atomType.append(element[3])
f.close()

#### read ff ####
record=False
ff={}
f=open('force_field_mixing_rules.def','r')
for x in f:
    if not x.strip():
        continue
    element=x.split()
    if len(element)==3:
        record=True
        continue
    if len(element)==6:
        break
    if record==True:
        ff[element[0]]=[element[2],element[3]]
f.close()
    

### makeJobs ###
for x in adsorbate:
    myList='['
    for ii,z in enumerate(atomType):
        if ii!=0:
            myList=myList+","
        myList=myList+"TT'"+z[0:2]+"TT'"
    
    for z in adsorbate_atom_dict[x]:
        myList=myList+","
        temp=z[0:2].replace('H','')
        myList=myList+"TT'"+temp+"TT'"
        
    myList=myList.replace('_','')
    myList=myList+']'
    
    os.system('sed "s/myList/"'+myList+'"/g" MCMD.py > MCMD.py_T' )
    os.system('sed -i "s/TT/\'/g" MCMD.py_T' )
    
    os.system('mkdir '+x)
    os.system('sed "s,myAd,"'+x+'",g" MCMD.py_T > '+x+'/MCMD.py')
    os.system('rm MCMD.py_T')
    os.system('sed "s,myAd,"'+x+'",g" submit.sh > '+x+'/submit.sh')
    os.system('cp in.* data.* clean.sh *.cif force_field_mixing_rules.def pseudo_atoms.def force_field.def '+x)
    os.system('cp dataFile/'+x+'.def '+' dataFile/data.'+x+' '+x)
    if x=='propane' or x =='propene':
        os.system('cp mof_taxi.py '+x)
    else:
        os.system('cp mof_taxi_noBond.py '+x+'/mof_taxi.py')
      
    f=open('force_field_mixing_rules_4lmps.def','r')
    f_out=open(x+'/force_field_mixing_rules_4lmps.def','a+')
    for i,y in enumerate(f):
        if i==5:
            f_out.write(str(len(atomType)+len(adsorbate_atom_dict[x]))+'\n')
        elif i==7:
            for j in range(0,len(atomType)):
                f_out.write(atomType[j][0:2]+' '+'LENNARD_JONES '+ff[atomType[j][0:2]][0]+' '+ff[atomType[j][0:2]][1]+'\n')
            for m in range(len(adsorbate_atom_dict[x])):
                f_out.write(adsorbate_atom_dict[x][m]+' '+'LENNARD_JONES '+str(ff[adsorbate_atom_dict[x][m]][0])+' '+str(ff[adsorbate_atom_dict[x][m]][1])+'\n')                
        else:
            f_out.write(y)
    f.close()
    f_out.close()
                
            
    
    