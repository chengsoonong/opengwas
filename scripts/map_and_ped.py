#!/usr/bin/python

import sn
import sys,string
import numpy as np
import math
import csv
import os.path
from collections import namedtuple
import os
import vcf
import fnmatch
from optparse import OptionParser
import time



class MapAndPed:
    """
    This classe allow create and parse .map and .ped files to be used in PLINK.
    """

    def __init__(self, outputmap, outputped,outputmap_parse=None, outputped_parse=None):
        """
        Initialaze the output names of the files and other variables
        """

        #self.outputmap = "Phobia-test5.map"
        #self.outputped = "Phobia-test5.ped"
        self.outputmap = outputmap
        self.outputped = outputped
        self.outputmap_parse = outputmap_parse
        self.outputped_parse = outputped_parse
        self.idx_rem=[]
        self.all_rs=[]
        self.valid_alleles = ['AA', 'AT', 'AG', 'AC',
                              'CA', 'CT', 'CG', 'CC', 
                              'TA', 'TT', 'TG', 'TC',
                              'GA', 'GT', 'GG', 'GC']
        self.Family_ID=1
        self.Individual_ID=1
        self.chrommosomos=np.array(["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"])        



    def create_map(self, folders_case, folders_control):
        """Create the .map file"""
        print "\n\nMAP FILE   (",time.asctime(),")\n"
        dir_list=folders_case[:]
        dir_list.extend(folders_control)
        files=[]
        file_dir={}
        map_list=[]    
        idx=0

        for i in dir_list:
                 files.extend(os.listdir(i)) #it get all files names in a dir
                 for j in os.listdir(i):   
                    file_dir[j]=i            #dictionari with file:dir

        print "Reading the files:\n\n",file_dir[files[0]]+"/"+files[0],  #parse the first file
        try:
            snps = sn.parse(file_dir[files[0]]+"/"+files[0])#take the first file	        	
            for i in snps:              #initialaze rs_list and map_list with the first file
                map_list.append((i[0],i[2],i[3]))
        except:
                 print "   ERRO 1"             
        print ""
    
        dtype = [('rs', 'S10'), ('chrom', 'S10'), ('position', int)]       
        map_list = np.array( map_list, dtype=dtype)        
    
        for j in files[1:]:
            map_list_tmp=[]    
            print file_dir[j]+"/"+j,             
            try:	
                snps = sn.parse(file_dir[j]+"/"+j)       #take another files
            except:
                 print  "   ERRO 2"
                 continue   
            try:	        	
              for i in snps:
                map_list_tmp.append((i[0],i[2],i[3]))
            except:
                 print "   ERRO 3"
                 continue   
            print ""

            map_list_tmp=np.array(map_list_tmp, dtype=dtype)             
            map_list=np.array(np.concatenate((  map_list, map_list_tmp)) , dtype=dtype)    
            u, indices = np.unique( map_list['rs'] , return_index=True )
            map_list = map_list[indices]
        
        array_chrom=np.unique( map_list['chrom'])    #add new elements to end of the self.chrommosomos    
        idx_chr=np.in1d(array_chrom,self.chrommosomos)
        self.chrommosomos=np.concatenate(( self.chrommosomos , array_chrom[idx_chr==False]))

        map_list = np.sort(map_list, order=['chrom', 'position']) 
        ofile = open(self.outputmap,'w')    # open file for writing
        print "there are",len(map_list['rs']),"SNPs.\nwriting the",self.outputmap,"file..."
        for i in self.chrommosomos:
            if i in map_list['chrom']:
                idx = map_list['chrom']==i                                  
                for i in map_list[:][idx]:
                    ofile.write(str(i[1])+" "+str(i[0])+" "+str(0)+" "+str(i[2])+"\n")                
        ofile.close()



    def create_ped( self, folders_case, folders_control):
        """Create the .ped file"""

        print "\n\n\nPED FILE   (",time.asctime(),")\n"
        handle = csv.DictReader(open(self.outputmap, "r"),
                fieldnames=["chromosome", "rs","morgans", "position"],
                delimiter=" ")
        for i in handle:
            self.all_rs.append(i["rs"])
        self.all_rs=np.array(self.all_rs)
        ofile = open(self.outputped,'w')    # open file for writing        
        print "\nReading the file to be cases (affected: 2):\n"    
        for folder in folders_case:    
            self.write_ped( folder, ofile, "2")
        print "\nReading the file to be controls (unaffected: 1):\n"
        for folder in folders_control:
            self.write_ped(folder, ofile, "1")
        ofile.close()



    def write_ped(self, dirfilespheno, outputfile, pheno):
        """
        read the file inside a folder and parse to write a .ped.
        """

        for i in os.listdir(dirfilespheno):
            all_rs_ind_tmp = np.array(["0 0"]*len(self.all_rs), dtype='S3') #initialaze every genotype with "0 0"
            sex=""
            all_names=[]
            all_gen=[]
            print dirfilespheno+"/"+i,
            if "XY" in i:
                sex="1"        
            elif "XX" in i:
                sex="2"       
            else:
                sex="9"        #sex="other"
            try:
                snps = sn.parse(dirfilespheno+"/"+i)       #take another files
            except:
                print  "   ERRO 1"
                continue 
            try:
                for cur_snp in snps:
                    if len(cur_snp.genotype)==2 and cur_snp.genotype in self.valid_alleles:# "--" and cur_snp.genotype != "II" and cur_snp.genotype != "DI":
                        all_names.append(cur_snp.name)
                        all_gen.append("%s %s" % (cur_snp.genotype[0], cur_snp.genotype[1]))
            except:
                print  "   ERRO 2"
                continue    
            try:            
                idx = np.flatnonzero(np.in1d(self.all_rs, np.array(all_names)))
            except:
                print  "   ERRO 3"
                continue    
            print ""                 
            all_rs_ind_tmp[idx] = np.array(all_gen)
            outputfile.write( str(self.Family_ID)+" "+ str(self.Individual_ID)+" "+"0 0 "+sex+"  "+ pheno+"  ")
            for i in all_rs_ind_tmp:
                    outputfile.write(i+"  ")
            outputfile.write("\n")
            self.Family_ID=self.Family_ID+1
            self.Individual_ID=self.Individual_ID+1





    def parse_ped (self):
        """
        Parse the .ped to avoid more than 2 alleles. 
        """
        print "\n\nPARSE PED  (",time.asctime(),")\n"
        print "\nparsing the",self.outputped,"file ..."
        #take the file .ped
        ifile_ped = open(self.outputped,'r') #folder_test.ped test_all_files_blue_brown.ped zeros( (3,4) )
        #create a matrix of the zeros to count how much allelos
        ACTG_matrix = np.zeros( (4,len(ifile_ped.readline().split())) )
        ifile_ped.close()
        #take the file .ped again
        ifile_ped = open(self.outputped,'r')

        for i in ifile_ped:   
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
        self.idx_rem=[]
        idx_keep=[]

        for c in np.flatnonzero(np.array(range(ACTG_matrix.shape[1]))%2==0):       
            u=np.sum(ACTG_matrix[:,c:c+2], axis=1)    
            if len(np.delete(u, np.flatnonzero(u==0))) >2:
                self.idx_rem.append(c)
                self.idx_rem.append(c+1)
            else:
                idx_keep.append(c)
                idx_keep.append(c+1)

        self.idx_rem=np.array(self.idx_rem)
        idx_keep=np.array(idx_keep)
        ofile_ped = open(self.outputped_parse,'w')
        ifile_ped = open(self.outputped,'r')   

        print "writing the",self.outputped_parse,"file ..."
        for l in ifile_ped:
            line=np.array(l.split())
            ofile_ped.write(line[0]+" "+line[1]+" "+line[2]+" "+line[3]+" "+line[4]+"  "+line[5]+"  ")    
            lines=line[6:][idx_keep]
            for c in np.flatnonzero(np.array(range( len(lines)  ))%2==0) :
                ofile_ped.write(lines[c]+" "+lines[c+1]+"  ")    
            ofile_ped.write("\n")
        ifile_ped.close()
        ofile_ped.close()






    def parse_map (self):
        """
        Parse the .ped to avoid more than 2 alleles.
        """
        print "\n\nPARSE MAP  (",time.asctime(),")\n"
        print "\nparsing the",self.outputmap ,"file ..."
        #take the file .map
        dtype = [('chromosome', 'S10'), ('rs', 'S10'), ('morgans', 'S10'),('position', 'S10')]
        map_array = np.genfromtxt(self.outputmap, dtype=dtype, delimiter=" ")
        #get the  idx the columns to be removed in map_array
        idx_del_map= self.idx_rem
        idx_del_map = idx_del_map[idx_del_map%2 == 0]
        idx_del_map=idx_del_map/2
        map_array = np.delete(map_array,idx_del_map, axis=0)
        print "writing the",self.outputmap_parse ,"file ..."
        np.savetxt(self.outputmap_parse,map_array,delimiter=' ',fmt='%s',newline='\n')




    def parse_map_ped(self):
        """
        Parse the .map and .ped file to avoid more than 2 alleles.
        """
        self.parse_ped()
        self.parse_map()

       


if __name__ == '__main__':

    # build option parser:
    class MyParser(OptionParser):
        def format_epilog(self, formatter):
            return self.epilog
    
    usage = "usage: python %prog [options] filename\n"    

    description = """
This program allow us create the .map and .ped files to be used in plink.\n"""

    epilog = """

For example:

python map_and_ped.py -m filename.map -p filename.ped --case "folder1 folder2 folder3" --control "folder4 folder5" 

    INPUT:  
        "folder1 folder2 folder3"
        "folder4 folder5"                
    OUTPUT
        filename.map 
        filename.ped

    If you use the output files filename.map and filename.ped in PLINK. You will get a error similar to below:
        ERROR: Locus rs2055204 has >2 alleles:
           individual 2 1 has genotype [ C C ]
           but we've already seen [ A ] and [ G ]

python map_and_ped.py -m filename.map -p filename.ped --case "folder1 folder2 folder3" --control "folder4 folder5" --omp filename-parsed.map --opp filename-parsed.ped

    INPUT:  
        "folder1 folder2 folder3"
        "folder4 folder5"                
    OUTPUT
        filename.map 
        filename.ped
        filename-parsed.map
        filename-parsed.ped  

    You can use the output files filename-parsed.map and filename-parsed.ped in PLINK.



"""

    parser = MyParser(usage, description=description,epilog=epilog)
    parser.add_option("--case",  dest="case", action="store",
                      help='input  -   folders with the files representing case. Put the folders inside "". for example: --case "folder1 folder2 folder3"')
    parser.add_option("--control", dest="control",   action="store",
                      help='input  -   folders with the files representing control. Put the folders inside "". for example: --control "folder4 folder5 folder6"')
    parser.add_option("-m", "--outfile_map", dest="outfile_map",  action="store",
                      help="output  -  file name of the .map.")
    parser.add_option("-p","--outfile_ped",  dest="outfile_ped",  action="store",
                      help="output  -  file name of the .ped.")    
    parser.add_option("--omp", dest="outfile_map_parse",   action="store",
                      help="output  -  file name of the .map to be parsed to be used in plink.")    
    parser.add_option("--opp", dest="outfile_ped_parse",   action="store",
                      help="output  -  file name of the .ped to be parsed to be used in plink")    

    (options, args) = parser.parse_args()

    if len(sys.argv) != 9 and len(sys.argv) != 13:
        parser.error("incorrect number of arguments. Use -h to help you.")
    
    outfile_map = options.outfile_map
    outfile_ped = options.outfile_ped  
    outfile_map_parse = options.outfile_map_parse
    outfile_ped_parse = options.outfile_ped_parse
    case = options.case.split() 
    control = options.control.split() 
    

    if (outfile_ped_parse == None or outfile_map_parse == None):
 
        mp = MapAndPed(outfile_map, outfile_ped)
        mp.create_map(case, control)
        mp.create_ped(case, control)

    else:

        mp = MapAndPed(outfile_map, outfile_ped, outfile_map_parse, outfile_ped_parse )
        mp.create_map(case, control)
        mp.create_ped(case, control)
        mp.parse_map_ped()









