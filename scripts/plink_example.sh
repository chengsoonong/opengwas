#!/bin/bash
# Running plink to filter for missing genotypes (per-SNP and per-person).
# Use plink again to find epistatic interactions and collect the results.
# Also create some files that are useful for downstream processing.

# Run the plink with map and ped files for you to create the files .cc with SNP pairs.

plink --noweb --allow-no-sex --geno 0.1 --make-bed\
    --file filename-parsed-output --out filename-parsed-output-geno

plink --noweb --allow-no-sex --mind 0.1 --make-bed\
    --bfile filename-parsed-output-geno --out filename-plink

# Use plink again to find epistatic interactions and collect the results.
plink --noweb --allow-no-sex --geno 0.1 --epistasis\
    --bfile filename-plink --out filename-plink-geno-01-epi 

# Run the plinkcc2txt.py with the .epi.cc with information about the SNPs pairs
# This script write a output file (outfile) from information in files .epi.cc and .bim.
# The output file has 11 columns and the last two columns are ID (prb1 and prb2) for each snp.
# The 11 columns: CHR1 | SNP1 | CHR2 | SNP2 | OR_INT | STAT | P | POSITION1 | POSITION2 | prb1 | prb2

python plinkcc2txt.py   --cc  file.cc   --bim  filename-plink.bim  --out output_snp_pairs.txt

#Run plink2hdf5.py to create the h5 file
python plink2hdf5.py filename-plink
