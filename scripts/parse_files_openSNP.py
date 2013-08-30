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
import fnmatch
from optparse import OptionParser


def get_dataset(inputfile):
    """return """
    handle = csv.DictReader(open(inputfile, "r"), delimiter=";")
    return handle

def get_user(pheno, variation,inputfile):
    """Return list of the user with a specific variation """
    dataset = get_dataset(inputfile)
    user_list = []
    for i in dataset:
        if i[pheno] == variation:
            user_list.append(i["user_id"])
    dataset=[]
    return user_list


def create_dir(user_list,variation,folder):
    """Create a folder from a list of the user"""
    user_list=list(set(user_list))
    print "total of the user", len(user_list), user_list
    files= os.listdir(folder)
    os.system("mkdir "+variation)
    n=0
    for j in user_list: 
        for i in files:
          if fnmatch.fnmatch(i, '*.txt'):   
            u="user"+j+"_"
            if u in i:
                print i
                os.system("cp "+folder+"/"+i +" " +variation+"/")
                n=1+n
    print "total of the files copied", n        


def dir_pheno_variation(file_csv,folder_genotype):
    """ Read the openSNP file with all phenotypes to create
        the folders with phenotype and insede the dir with 
        files insede each variations folder. """

    fieldnames=open(file_csv).readline().split(';')  #read all phenotypes and sort
    fieldnames.sort()

    print "\n\n---------------------------  fieldnames (Phenotypes)\n"

    for i in fieldnames:
        print i

    p=raw_input("\n---------------------------  Choose the Phenotype above: ")
    variations_list=[] 
    for i in get_dataset(file_csv):
        if not i[p] in variations_list: 
            variations_list.append(i[p])
            print i[p]        

    v=raw_input("\n---------------------------  Choose the Variations above (use ';' between them): ")
    v=v.split(";")

    print "\n"

    os.system("mkdir "+"_".join(p.split()))

    for i in v:
        print "Variations: ", i
        l=get_user( p, i,file_csv)
        variation="_".join(i.split())
        create_dir(l,variation,folder_genotype)
        os.system("mv "+ variation+" "+"_".join(p.split()))    
        print "\n"


if __name__ == '__main__':

    # build option parser:
    class MyParser(OptionParser):
        def format_epilog(self, formatter):
            return self.epilog
    
    usage = "usage: python %prog [options] filename and path/folder\n"    

    description = """
This program allow us parse the .csv that has information about genotype and phenotype. 
From this information we can create a folder for a phenotype and inside we can have 
folders with  genotype files for each variations\n"""

    epilog = """
Input: 
    file .csv that contains infomation about phenotypes and genotypes for each individual.
Output:
    phenotype folder with all variations\n



    when the program run it will request your preferences. for exemple:


----------------------------------  fieldnames (Phenotypes)

...
ADHD
Asthma
Ability to Tan
...

---------------------------  Phenotype: Asthma
...
False
Allergy induced, mainly when younger
No
Chronic Asthma
slight
...
---------------------------  Variations: No;Chronic Asthma;slight

"""

    parser = MyParser(usage, description=description,epilog=epilog)
    parser.add_option("-f", "--file_csv", dest="file_csv", 
                      help="file name of the .csv to be parse")
    parser.add_option("-g","--folder_genotyphe",  dest="folder_genotype", 
                      help="folder with all genotype files")    


    (options, args) = parser.parse_args()
    if len(sys.argv) != 5:
        parser.error("incorrect number of arguments. Use -h to help you.")
    file_csv = options.file_csv
    folder_genotype = options.folder_genotype   
    dir_pheno_variation(file_csv, folder_genotype)





