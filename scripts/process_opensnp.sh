#!/bin/bash
# A brief step-by-step to create the plink files from openSNP. 
# Please, use -h to see information about the programs. 
# =======================================================

WORK_DIR=/path/to/data/

# 1 - Run the parse_files_openSNP.py  to create the folders for a phenotype.

# This program will show all phenotypes inside the folder that 
# has all user genotype which were provide to download in openSNP web-site.
# Next choose the phenotype, you need choose the Variations.
# After choose the phenotype and variations the program will create one folder
# with the name of the phenotype that was choosen and inside this folder the program will
# create the folders with the name of each variation and the each variation folders
# has the files with genotype of the users.
# Please, use -h to get more information 
# (python parse_files_openSNP.py -h).

# -f -> name of the file .csv of openSNP
# -g -> name of the folder with the genotype file of the users 
GENO_DIR=$WORK_DIR/opensnp_datadump.201303070733
PHENO_FILE=$GENO_DIR/phenotypes_201303070733.csv

python parse_files_openSNP.py -f $PHENO_FILE  -g $GENO_DIR

# 2 - Run the map_and_ped.py. 

# This program will create the files in plink format map and ped 
# from the phenotype folder created by parse_files_openSNP.py.
# The input will be the variations folders that has the genotype file of the users.
# The variantions that are case should have like input, for example, "phenotype/variation1" 
# and control should have like input, for example, "phenotype/variation2".

# example:
#
# python ../map_and_ped.py -m asthma.map -p asthma.ped --case "Asthma,_mostly_during_childhood Had_as_a_child,_not_anymore Allergy_induced,_mainly_when_younger  Asthma,_mostly_during_Childhood  Exercise-induced  Asthmatic, slight Cold-induced,_worse_when_younger  Chronic_Asthma"   --control "No  False  no"   --omp asthma-parsed.map --opp asthma-parsed.ped

# TODO: Write the exact command line instructions for generating
#       the three datasets referred to in the paper: Asthma, Dyslexia and Lactose Intolerance

python ../scripts/map_and_ped.py -m filename-output.map -p filename-output.ped --case "phenotype/variation1" --control "phenotype/variation2" --omp filename-parsed-output.map --opp filename-parsed-output.ped






#6 - Run the gss2graph.py that is inside the chillo/tmp-scripts with the small implementation to create the json file from openSNP
# TODO: Move this to the chillo repository

python  chillo/plink-scripts gss2graph.py --inputtxt output_snp_pairs.txt  --h5name filename-plink --outjson filename-plink.json



