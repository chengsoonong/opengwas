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



try:	
	file_map = sys.argv[1];file_ped = sys.argv[2];outfile_map = sys.argv[3];outfile_ped = sys.argv[4]

except:
     print "Usage:",sys.argv[0], "file.map dir_files_phenotype1 dir_files_phenotype2 outfile"; sys.exit(1)





#take the file .ped
ifile_ped = open(file_ped,'r') #folder_test.ped test_all_files_blue_brown.ped zeros( (3,4) )

#create a matrix of the zeros to count how much allelos
ACTG_matrix = np.zeros( (4,len(ifile_ped.readline().split())) )
print "ACTG_matrix.shape ", ACTG_matrix.shape, ACTG_matrix,"\n"
ifile_ped.close()

#take the file .ped again
ifile_ped = open(file_ped,'r')


for i in ifile_ped:
    print "individual",i[:20]
    
    line=i.split()
    idx_alle= np.array(line)=="A"
    ACTG_matrix[0][idx_alle]=1

    idx_alle= np.array(line)=="C"
    ACTG_matrix[1][idx_alle]=1

    idx_alle= np.array(line)=="T"
    ACTG_matrix[2][idx_alle]=1

    idx_alle= np.array(line)=="G"
    ACTG_matrix[3][idx_alle]=1

ifile_ped.close()




ACTG_matrix= ACTG_matrix[:,6:]

print "ACTG_matrix.shape ", ACTG_matrix.shape, ACTG_matrix,"\n"

print "calculating..."

idx_rem=[]
idx_keep=[]
for c in np.flatnonzero(np.array(range(ACTG_matrix.shape[1]))%2==0):   
    
    u=np.sum(ACTG_matrix[:,c:c+2], axis=1)    
    if len(np.delete(u, np.flatnonzero(u==0))) >2:
        idx_rem.append(c)
        idx_rem.append(c+1)
    else:
        idx_keep.append(c)
        idx_keep.append(c+1)




idx_rem=np.array(idx_rem)
idx_keep=np.array(idx_keep)


#1 0 0 1  1  AA
#1 0 0 other  1  AA




print "\nwriting the .ped file ..."
ofile_ped = open(outfile_ped,'w')
ifile_ped = open(file_ped,'r')   

for l in ifile_ped:

    line=np.array(l.split())
    ofile_ped.write(line[0]+" "+line[1]+" "+line[2]+" "+line[3]+" "+line[4]+"  "+line[5]+"  ") 

    lines=line[6:][idx_keep]
    

    for c in np.flatnonzero(np.array(range( len(lines)  ))%2==0) :
        ofile_ped.write(lines[c]+" "+lines[c+1]+"  ")    
    ofile_ped.write("\n")

    

ofile_ped.close()






#take the file .map
dtype = [('chromosome', 'S10'), ('rs', 'S10'), ('morgans', 'S10'),('position', 'S10')]
map_array = np.genfromtxt(file_map, dtype=dtype, delimiter=" ")

print "map_array ", map_array.shape, "\n",map_array,"\n"
#get the  idx the columns to be removed in map_array
idx_del_map= idx_rem
idx_del_map = idx_del_map[idx_del_map%2 == 0]
idx_del_map=idx_del_map/2
map_array = np.delete(map_array,idx_del_map, axis=0)
print "map_array", map_array.shape
print "\nwriting the .map file ..."
np.savetxt(outfile_map,map_array,delimiter=' ',fmt='%s',newline='\n')




