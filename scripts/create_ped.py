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
	file_map = sys.argv[1];dir_files_phenotype1 = sys.argv[2];dir_files_phenotype2 = sys.argv[3];outfilename = sys.argv[4]

except:
     print "Usage:",sys.argv[0], "file.map dir_files_phenotype1 dir_files_phenotype2 outfile"; sys.exit(1)










handle = csv.DictReader(open(file_map, "r"),
        fieldnames=["chromosome", "rs","morgans", "position"],
        delimiter=" ")


all_rs=[]  #take all rs
for i in handle:
    all_rs.append(i["rs"])
all_rs=np.array(all_rs)
print len(all_rs)

ofile = open(outfilename,'w')    # open file for writing


print "search"
print "phe1\n"

Family_ID=1


for i in os.listdir(dir_files_phenotype1):

    all_rs_ind_tmp = np.array(["0 0"]*len(all_rs), dtype='S3') #initialaze every genotype with "0 0"    
    print len(all_rs_ind_tmp),'initialaze array with every genotype with "0 0"'    


    sex=""

    print i
    if "XY" in i:
        sex="1"
        print "m"
    elif "XX" in i:
        sex="2"
        print "f"
    else:
        sex="other"
        print "o"

    snps = sn.parse(dir_files_phenotype1+"/"+i)       #take another files

    all_names=[]
    all_gen=[]

    for cur_snp in snps:
        if len(cur_snp.genotype)==2 and cur_snp.genotype != "--":
            all_names.append(cur_snp.name)
            all_gen.append("%s %s" % (cur_snp.genotype[0], cur_snp.genotype[1]))


    idx = np.flatnonzero(np.in1d(all_rs, np.array(all_names))) #return the index of the elements of "all_names" in "all_rs"
    all_rs_ind_tmp[idx] = np.array(all_gen)

   
    
    print "write"

    ofile.write( str(Family_ID)+" "+ "1 0 0 "+sex+"  1  ")
    for i in all_rs_ind_tmp:
            ofile.write(i+"  ")
    ofile.write("\n")
    Family_ID=1+Family_ID





print "phe2\n"



for i in os.listdir(dir_files_phenotype2):

    all_rs_ind_tmp = np.array(["0 0"]*len(all_rs), dtype='S3') #initialaze every genotype with "0 0"    
    print len(all_rs_ind_tmp),'initialaze array with every genotype with "0 0"'    


    sex=""

    print i
    if "XY" in i:
        sex="1"
        print "m"
    elif "XX" in i:
        sex="2"
        print "f"
    else:
        sex="other"
        print "o"

    snps = sn.parse(dir_files_phenotype2+"/"+i)       #take another files

    all_names=[]
    all_gen=[]

    for cur_snp in snps:
        if len(cur_snp.genotype)==2 and cur_snp.genotype != "--":
            all_names.append(cur_snp.name)
            all_gen.append("%s %s" % (cur_snp.genotype[0], cur_snp.genotype[1]))


    idx = np.flatnonzero(np.in1d(all_rs, np.array(all_names)))
    all_rs_ind_tmp[idx] = np.array(all_gen)

   
    print "write"   

    ofile.write( str(Family_ID)+" "+ "1 0 0 "+sex+"  2  ")
    for i in all_rs_ind_tmp:
            ofile.write(i+"  ")
    ofile.write("\n")
    Family_ID=1+Family_ID



ofile.close()
