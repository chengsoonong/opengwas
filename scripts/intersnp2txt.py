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







def get_data(filename):
    """
    This function read the file bestMarkerCombi2.txt and return it.
    """
    data = genfromtxt(filename,skip_header=1,
    dtype={'names':[ "No","Chr_1","rs_No_1","Pos_No_1","p_Single_marker_1","Chr_2","rs_No_2","Pos_No_2","p_Single_marker_2","p_value","p_value_corr"], 'formats':['S16','S16','S16','S16','S16','S16','S16','S16','S16','S16','S16']})

    return data

def file_with_prb(infile,outfile):
    """
    This function read the file bestMarkerCombi2.txt (that has results of two-marker-analysis)
    and write a file with the name 'outfile'. This outfile has the same information that 
    bestMarkerCombi2.txt, but it has more two columns than bestMarkerCombi2.txt. This two columns
    are ID (prb1 and prb2) for each snp . 
    """

    print "Reading ..."
    data = get_data(infile)
    rs =    np.unique(np.append( data["rs_No_1"],data["rs_No_2"] ))  # create the prb1 and prb2
    shape_data=np.shape(data)
    ofile = open(outfile,'w')    # open file for writing 

    print "Writing ..."
    for i in range(shape_data[0]):
        ofile.write(str(data[i]['Chr_1'])+" "
        +str(data[i]['rs_No_1'])+" "
        +str(data[i]['Chr_2'])+" "
        +str(data[i]['rs_No_2'])+" "
        +str(data[i]['p_value'])+" "
        +str(data[i]['p_value_corr'])+" "
        +str(data[i]['Pos_No_1'])+" "
        +str(data[i]['Pos_No_2'])+" "
        +str(np.flatnonzero(rs==data[i]["rs_No_1"])[0])+" "    # create the prb1
        +str(np.flatnonzero(rs==data[i]["rs_No_2"])[0])+"\n")  # create the prb2


    
    ofile.close


  

if __name__ == '__main__':

    # build option parser:
    class MyParser(OptionParser):
        def format_epilog(self, formatter):
            return self.epilog
    
    usage = "usage: python %prog [options] filename\n"    

    description = """ 
This program read the file bestMarkerCombi2.txt (that has results of two-marker-analysis)
and write a file with the name 'outfile'. This outfile has the same information that 
bestMarkerCombi2.txt, but it has more two columns than bestMarkerCombi2.txt. This two columns
are ID (prb1 and prb2) for each snp.\n 

Columns in the file outfile:\n

Chr_1 | rs_No_1 | Chr_2 | rs_No_2 | p_value | p_value_corr | Pos_No_1 | Pos_No_2 | prb1 | prb2
"""

    epilog = """ """

    parser = MyParser(usage, description=description,epilog=epilog)
    parser.add_option("--ifile",  dest="input_file", action="store",
                      help='file name of output of the INTERSNP (bestMarkerCombi2.txt)')
    parser.add_option("--ofile", dest="output_file",   action="store",
                      help='file name of output.')
 

    (options, args) = parser.parse_args()   
    input_file = options.input_file
    output_file = options.output_file
 
    file_with_prb(input_file, output_file)   
#    file_with_prb(ifile,"outfile.txt")


