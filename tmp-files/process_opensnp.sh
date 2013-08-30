# script that converts a file downloaded from openSNP to a PLINK format GWAS dataset 
# Please, use -h to see information about the programs. 
# =======================================================


# 1 - Run the parse_files_openSNP.py  to create the folders for a phenotype.

# This program will show all phenotypes inside the folder that has all user genotype which were provide to download in openSNP web-site.
# Next choose the phenotype, you need choose the Variations. After choose the phenotype and variations the program will create one folder
# with the name of the phenotype that was choosen and inside this folder the program will create the folders with the name of each variation and the each variation folders has the files with genotype of the users.
# Please, use -h to get more information 
# (python parse_files_openSNP.py -h).

# -f -> name of the file .csv of openSNP
# -g -> name of the folder with the genotype file of the users 

python ../scripts/parse_files_openSNP.py -f your/path/pensnp_datadump.201303070733/phenotypes_201303070733.csv  -g your/path/opensnp_datadump.201303070733

# To create the folders to Asthma, Lactose intolerance and Dyslexia use:

# phenotype = Asthma 
# variations = No;Asthma, mostly during childhood;Had as a child, not anymore;False;Allergy induced, mainly when younger;no;Asthma, mostly during Childhood;Exercise-induced;slight;Asthmatic,;Cold-induced, worse when younger;Chronic Asthma

# phenotype = Lactose intolerance
# variations = Partially lactose intolerant - surfaces with a lot of dairy products in one day;Lactose-intolerant;Lactose-tolerant;Lactose tolerant;lactose-tolerant;False;Genetically intolerant but drink raw milk and eat lots of dairy.;lactose-intolerant;AA;lactose tolerant

# phenotype = Dyslexia
# variations = No;Yes


# 2 - Run the map_and_ped.py. 

# This program will create the files in plink formart map and ped 
# from the phenotype folder created by parse_files_openSNP.py.
# The input will be the variations folders that has the genotype file of the users.
# The variantions that are case should have like input, for example, "phenotype/variation1" 
# and control should have like input, for example, "phenotype/variation2".

# example to create the plink files for Asthma, Lactose intolerance and Dyslexia:

# Lactose intolerance
python ../scripts/map_and_ped.py  -m lactose_int.map -p lactose_int.ped --opp lactose_int-parsed.ped --omp lactose_int-parsed.map --case  "Lactose-intolerant  Partially_lactose_intolerant_-_surfaces_with_a_lot_of_dairy_products_in_one_day Genetically_intolerant_but_drink_raw_milk_and_eat_lots_of_dairy."  --control "Lactose-tolerant  Lactose_tolerant  False"

# Asthma
python ../scripts/map_and_ped.py -m asthma.map -p asthma.ped --case "Asthma,_mostly_during_childhood Had_as_a_child,_not_anymore Allergy_induced,_mainly_when_younger  Asthma,_mostly_during_Childhood  Exercise-induced  Asthmatic, slight Cold-induced,_worse_when_younger  Chronic_Asthma"   --control "No  False  no"   --omp asthma-parsed.map --opp asthma-parsed.ped

# Dyslexia
python ../scripts/map_and_ped.py -m dyslexia.map -p dyslexia.ped --case "Yes"  --control "No"  --omp dyslexia-parsed.map --opp dyslexia-parsed.ped



# Bellow produces the JSON file from plink results
# ================================================

# 3 - Run the plink with map and ped files for you to create the files .cc with SNP pairs.

# Lactose intolerance
plink  --noweb --allow-no-sex --file lactose_int-parsed  --geno 0.1      --make-bed      --out  lactose_int-parsed-geno
plink  --noweb --allow-no-sex --bfile lactose_int-parsed-geno --mind 0.1 --make-bed --out lactose_int
plink  --bfile lactose_int  --out lactose_int-plink-epi-0001  --noweb  --allow-no-sex  --epi1 0.001   --epistasis


# Asthma
plink  --noweb --allow-no-sex --file asthma-parsed  --geno 0.1      --make-bed      --out  asthma-parsed-geno
plink  --noweb --allow-no-sex --bfile asthma-parsed-geno --mind 0.1 --make-bed --out asthma
plink  --bfile asthma  --out asthma-plink-epi-0001  --noweb  --allow-no-sex      --epi1 0.001       --epistasis


# Dyslexia
plink  --noweb --allow-no-sex --file dyslexia-parsed  --geno 0.1  --make-bed  --out  dyslexia-parsed-geno
plink  --noweb --allow-no-sex --bfile dyslexia-parsed-geno --mind 0.1 --make-bed --out dyslexia
plink  --bfile dyslexia  --out dyslexia-plink-epi  --noweb --allow-no-sex  --epistasis


# 4 - Run the plinkcc2txt.py with the .epi.cc with information about the SNPs pairs
# This script write a output file (outfile) from information in files .epi.cc and .bim.
# The output file has 11 columns and the last two columns are ID (prb1 and prb2) for each snp.
# The 11 columns: CHR1 | SNP1 | CHR2 | SNP2 | OR_INT | STAT | P | POSITION1 | POSITION2 | prb1 | prb2

# Lactose intolerance
python ../scripts/plinkcc2txt.py --cc  lactose_int-plink-epi-0001.epi.cc --bim lactose_int.bim --out lactose_int-plink-epi-0001.epi.txt
# Asthma
python ../scripts/plinkcc2txt.py --cc  asthma-plink-epi-0001.epi.cc --bim asthma.bim --out asthma-plink-epi-0001.epi.txt
# Dyslexia
python ../scripts/plinkcc2txt.py --cc  dyslexia-plink-epi.epi.cc --bim dyslexia.bim --out dyslexia-plink-epi.epi.txt


#5 - Run plink2hdf5.py to create the h5 file

# Lactose intolerance
python ../scripts/plink2hdf5.py lactose_int
# Asthma
python ../scripts/plink2hdf5.py asthma
# Dyslexia
python ../scripts/plink2hdf5.py dyslexia


#6 - Run the gss2graph.py that is inside the chillo/tmp-scripts with the small implementation to create the json file from openSNP

# Lactose intolerance
python chillo/plink-scripts/scripts/gss2graph.py --inputtxt lactose_int-plink-epi-0001.epi.txt  --h5name lactose_int --outjson lactose_int-plink-epi-0001.epi.json
# Asthma
python chillo/plink-scripts/scripts/gss2graph.py --inputtxt asthma-plink-epi-0001.epi.txt  --h5name asthma --outjson asthma-plink-epi-0001.epi.json
# Dyslexia
python chillo/plink-scripts/scripts/gss2graph.py --inputtxt dyslexia-plink-epi.epi.txt  --h5name dyslexia --outjson dyslexia-plink-epi.epi.json







