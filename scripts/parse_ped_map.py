#!/usr/bin/python


import csv
import os.path
from collections import namedtuple
import sn
import os
import sys,string
import numpy as np
import math
import vcf



#take the file .map
dtype = [('chromosome', 'S10'), ('rs', 'S10'), ('morgans', 'S10'),('position', 'S10')]
data = np.genfromtxt("test_all_files_blue_brown.map", dtype=dtype, delimiter=" ")
print data,data[[1,2,100]]
print  len(data['rs'])


#take the file .ped
ifile_ped = open("test_all_files_blue_brown.ped",'r') #folder_test.ped test_all_files_blue_brown.ped zeros( (3,4) )
ped_array=[]

for i in ifile_ped.readlines():
    ped_array.append(i.split())  
ped_array = np.array(ped_array)  

pheno_array=ped_array[:,:6]
snp_array= ped_array[:,6:]


print pheno_array[0],len(pheno_array[0]),"\n"
print snp_array[0][:10],snp_array[0][-10:],len(snp_array[0]),"\n"
print ped_array,"\n"

print ped_array[0][:10],ped_array[0][-10:],"\n"
idx = ped_array[0][:10]=="A"
print ped_array[0][idx],"\n"
print len(ped_array[0]),"\n"



ACTG_matrix = np.zeros( (4,len(snp_array[0])) )

print ACTG_matrix, len(ACTG_matrix[0])

ACTG_matrix[0,1]=1

print ACTG_matrix,"\n"



for i in snp_array:

    idx_A= i=="A"
    idx_C= i=="C"
    idx_T= i=="T"
    idx_G= i=="G"
    
    ACTG_matrix[0][idx_A]=1
    ACTG_matrix[1][idx_C]=1
    ACTG_matrix[2][idx_T]=1
    ACTG_matrix[3][idx_G]=1



print "ACTG_matrix parsed:\n",ACTG_matrix, len(ACTG_matrix[0]),"\n"



diffallelos=ACTG_matrix.sum(axis=0) 

print "diff allelos parsed:\n",diffallelos, "\n"

print "idx",np.flatnonzero(diffallelos[:]>4),len(np.flatnonzero(diffallelos[:]>4)),"\n"

idx = np.flatnonzero(diffallelos[:]>2)

ee = idx[idx%2 == 0]
oo = idx[idx%2 == 1]

eo = ee+1
oe = oo-1


snp_array[:,ee]="R"
snp_array[:,oo]="R"
snp_array[:,eo]="R"
snp_array[:,oe]="R"

    


print snp_array,snp_array[:,ee],snp_array[0][-10:],len(snp_array[0]),"\n"






#delete

idx_del = snp_array[:,:]=="R"
np.delete(snp_array, np.flatnonzero(idx_del[0]) ,axis=1)
idx_del_map= np.flatnonzero(idx_del[0])
















