#!/usr/bin/python


import sn
import os
import sys,string
import numpy as np
import math


try:	
	outfilename = sys.argv[1]

except:
     print "Usage:",sys.argv[0], "dir_files_input outfile"; sys.exit(1)





def creat_map(output_file="output.map"):
    """Create the .map file"""
    input_file=raw_input("\n-------folders with files: ")
    dir_list=input_file.split(";")
    files=[]
    file_dir={}    
    for i in dir_list:
             files.extend(os.listdir(i)) #all files names in a dir
             for j in os.listdir(i):   
                file_dir[j]=i 
                #print j,i  


    rs_list=[]    
    map_list=[]
    idx=0

    print "\nParse the file ", file_dir[files[0]]+"/"+files[0],  

    try:
        snps = sn.parse(file_dir[files[0]]+"/"+files[0])#take the first file	        	
        for i in snps:              #initialaze rs_list and map_list with the first file
            rs_list.append(i[0])
            map_list.append((i[0],i[2],i[3]))
    except:
             print "   ERRO 1"             

    print "  OK"

    dtype = [('rs', 'S10'), ('chrom', 'S10'), ('position', int)]       
    map_list = np.array( map_list, dtype=dtype)

 
    for j in files[1:]:

        print "Parse the file  ",file_dir[j]+"/"+j, 
        
        try:	
            snps = sn.parse(file_dir[j]+"/"+j)       #take another files
            rs_list_tmp=[]
            map_list_tmp=[]
        except:
             print  "   ERRO 2"
             continue   

        try:	        	
          for i in snps:
            rs_list_tmp.append(i[0])
            map_list_tmp.append((i[0],i[2],i[3]))
        except:
             print "   ERRO 3"
             continue   

        try:	                                                        #erro 3 in this files for exemplo
            rs_list_dif=np.setdiff1d( rs_list_tmp,  rs_list)            #user36_file207_yearofbirth_1986_sex_XY.23andme-exome-vcf.txt
        except:                                                         #user36_file327_yearofbirth_1986_sex_XY.23andme-exome-vcf.txt
             print "   ERRO 4"
             continue   


        idx = np.flatnonzero(np.in1d(np.array(rs_list_tmp),np.array(rs_list_dif))) #return the index of the elements of "all_names" in "all_rs"
        map_list_tmp=np.array(map_list_tmp, dtype=dtype)             
        map_list=np.array(np.concatenate((  map_list, map_list_tmp[idx]  )) , dtype=dtype)

        print "  OK"
 
    print "sort..."
    map_list = np.sort(map_list, order=['chrom', 'position']) 

    chrommosomos=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"]        


    ofile = open(output_file,'w')    # open file for writing
    
    print "write the output file..."
    for i in chrommosomos:
        if i in map_list['chrom']:

            idx = map_list['chrom']==i                                  

            for i in map_list[:][idx]:
                ofile.write(str(i[1])+" "+str(i[0])+" "+str(0)+" "+str(i[2])+"\n")
        
        
    ofile.close()
    




creat_map(outfilename)










