from numpy import genfromtxt
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




def get_data_cc(filename):
    """
    This function read and return the information inside the file .epi.cc"
    """
    data = genfromtxt(filename,skip_header=1,
                     dtype={'names':['CHR1', 'SNP1', 'CHR2','SNP2','OR_INT', 'STAT','P'],
                   'formats':['S16','S16','S16','S16',float,float,float]})
    return data


def get_data_bim(filename):
    """
    This function read and return the information inside the file .bim"
    """
    data = genfromtxt(filename,skip_header=1,
                     dtype={'names':["chromosome", "rs","morgans", "position","al1","al2"],
                   'formats':['S16','S16','S16',int,'S16','S16']})
    return data



def write_file(file_cc,file_bim,outfile):
    """
    This function write a output file (outfile) from information in .epi.cc and .bim.
    The output file has 11 columns and the last two columns are ID (prb1 and prb2) for each snp.

    columns:
    CHR1 | SNP1 | CHR2 | SNP2 | OR_INT | STAT | P | POSITION1 | POSITION2 | prb1 | prb2
    """

    print "Reading ..."
    data_cc = get_data_cc(file_cc)
    data_bim = get_data_bim(file_bim)

    rs =    np.unique(np.append( data_cc['SNP1'],data_cc['SNP2'] ))

    shape_data=np.shape(data_cc)

    ofile = open(outfile,'w')    # open file for writing 
    print "Writing ..."
    for i in range(shape_data[0]):

        ofile.write(str(data_cc[i]['CHR1'])+" "
        +str(data_cc[i]['SNP1'])+" "
        +str(data_cc[i]['CHR2'])+" "
        +str(data_cc[i]['SNP2'])+" "
        +str(data_cc[i]['OR_INT'])+" "
        +str(data_cc[i]['STAT'])+" "
        +str(data_cc[i]['P'])+" "  
        +str(data_bim['position'][data_bim['rs']==data_cc[i]['SNP1']][0])+" " #position1
        +str(data_bim['position'][data_bim['rs']==data_cc[i]['SNP2']][0])+" " #position2  
        +str(np.flatnonzero(rs==data_cc[i]['SNP1'])[0])+" " #id1
        +str(np.flatnonzero(rs==data_cc[i]['SNP2'])[0])+" " #id2  np.flatnonzero(
        +"\n"   )
    
    ofile.close


  

if __name__ == '__main__':

    # build option parser:
    class MyParser(OptionParser):
        def format_epilog(self, formatter):
            return self.epilog
    
    usage = "usage: python %prog [options] filename\n"    

    description = """ 
This script write a output file (outfile) from information in files .epi.cc and .bim.
The output file has 11 columns and the last two columns are ID (prb1 and prb2) for each snp.\n
The 11 columns:\t
CHR1 | SNP1 | CHR2 | SNP2 | OR_INT | STAT | P | POSITION1 | POSITION2 | prb1 | prb2
"""

    epilog = """ """

    parser = MyParser(usage, description=description,epilog=epilog)
    parser.add_option("--cc",  dest="cc", action="store",
                      help='file cc - case and control')
    parser.add_option("--bim", dest="bim",   action="store",
                      help='file .bim')
    parser.add_option("--out", dest="out",   action="store",
                      help='out file')
 

    (options, args) = parser.parse_args()   
    cc = options.cc
    bim = options.bim
    out = options.out
    

    write_file(cc,bim,out)











