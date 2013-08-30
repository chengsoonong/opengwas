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

# To create the folders to Asthma, Lactose intolerance and Dyslexia use:

# phenotype = Asthma 
# variations = No;Asthma, mostly during childhood;Had as a child, not anymore;False;Allergy induced, mainly when younger;no;Asthma, mostly during Childhood;Exercise-induced;slight;Asthmatic,;Cold-induced, worse when younger;Chronic Asthma

# phenotype = Lactose intolerance
# variations = Partially lactose intolerant - surfaces with a lot of dairy products in one day;Lactose-intolerant;Lactose-tolerant;Lactose tolerant;lactose-tolerant;False;Genetically intolerant but drink raw milk and eat lots of dairy.;lactose-intolerant;AA;lactose tolerant

# phenotype = Dyslexia
# variations = No;Yes


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

# Lactose intolerance
python map_and_ped.py  -m lactose_int.map -p lactose_int.ped --opp lactose_int-parsed.ped --omp lactose_int-parsed.map --case  "Lactose-intolerant  Partially_lactose_intolerant_-_surfaces_with_a_lot_of_dairy_products_in_one_day Genetically_intolerant_but_drink_raw_milk_and_eat_lots_of_dairy."  --control "Lactose-tolerant  Lactose_tolerant  False"

# Asthma
python map_and_ped.py -m asthma.map -p asthma.ped --case "Asthma,_mostly_during_childhood Had_as_a_child,_not_anymore Allergy_induced,_mainly_when_younger  Asthma,_mostly_during_Childhood  Exercise-induced  Asthmatic, slight Cold-induced,_worse_when_younger  Chronic_Asthma"   --control "No  False  no"   --omp asthma-parsed.map --opp asthma-parsed.ped

# Dyslexia
python map_and_ped.py -m dyslexia.map -p dyslexia.ped --case "Yes"  --control "No"  --omp dyslexia-parsed.map --opp dyslexia-parsed.ped



