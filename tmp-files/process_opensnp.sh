# A brief step-by-step to create the plink files from openSNP. 
# Please, use -h to see information about the programs. 
# =======================================================

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

python ../scripts/parse_files_openSNP.py -f /home/cristovao/Desktop/AUS_project/opengwas/public_datas/opensnp_datadump.201303070733/phenotypes_201303070733.csv  -g /home/cristovao/Desktop/AUS_project/opengwas/public_datas/opensnp_datadump.201303070733

# 2 - Run the map_and_ped.py. 

# This program will create the files in plink formart map and ped 
# from the phenotype folder created by parse_files_openSNP.py.
# The input will be the variations folders that has the genotype file of the users.
# The variantions that are case should have like input, for example, "phenotype/variation1" 
# and control should have like input, for example, "phenotype/variation2".

# example:
#
# python ../map_and_ped.py -m asthma.map -p asthma.ped --case "Asthma,_mostly_during_childhood Had_as_a_child,_not_anymore Allergy_induced,_mainly_when_younger  Asthma,_mostly_during_Childhood  Exercise-induced  Asthmatic, slight Cold-induced,_worse_when_younger  Chronic_Asthma"   --control "No  False  no"   --omp asthma-parsed.map --opp asthma-parsed.ped

python ../scripts/map_and_ped.py -m filename-output.map -p filename-output.ped --case "phenotype/variation1" --control "phenotype/variation2" --omp filename-parsed-output.map --opp filename-parsed-output.ped


# Bellow produces the JSON file from plink results
# ==================================================

# 3 - Run the plink with map and ped files for you to create the files .cc with SNP pairs.

plink  --noweb --allow-no-sex --file filename-parsed-output  --geno 0.1      --make-bed      --out  filename-parsed-output-geno

plink  --noweb --allow-no-sex --bfile filename-parsed-output-geno --mind 0.1 --make-bed --out filename-plink

plink  --bfile filename-plink  --out filename-plink-geno-01-epi  --noweb    --allow-no-sex     --geno 0.1  --epistasis 

# 4 - Run the plinkcc2txt.py with the .epi.cc with information about the SNPs pairs
# This script write a output file (outfile) from information in files .epi.cc and .bim.
# The output file has 11 columns and the last two columns are ID (prb1 and prb2) for each snp.
# The 11 columns: CHR1 | SNP1 | CHR2 | SNP2 | OR_INT | STAT | P | POSITION1 | POSITION2 | prb1 | prb2

python ../scripts/plinkcc2txt.py   --cc  file.cc   --bim  filename-plink.bim  --out output_snp_pairs.txt

#5 - Run plink2hdf5.py to create the h5 file

python ../scripts/plink2hdf5.py filename-plink

#6 - Run the gss2graph.py that is inside the chillo/tmp-scripts with the small implementation to create the json file from openSNP

python  chillo/plink-scripts gss2graph.py --inputtxt output_snp_pairs.txt  --h5name filename-plink --outjson filename-plink.json



