#!/bin/bash

# Author: Stefan Sevelda
# 1: parse json file and filter all rs numbers in seperate files for the
# chromsomes
# json file must be stored in the ld_generator directory
# make temporary directory 
mkdir tmp && cd tmp

if [ -z "$1" ]
then
    echo "you have to write the path of the input file"
    cd ..
    rm -r tmp
    exit 1
elif [ ! -f "../$1" ]
then 
    echo "this json file does not exist"
    cd ..
    rm -r tmp
    exit 1
elif [ -z "$2" ]
then
    echo "you forget to write a name for the result file"
    cd ..
    rm -r tmp
    exit 1
else
    echo "Json file parses to snap tabular file"
    python ../generatsnap.py ../$1
fi

# calculate the lds for the rede json file 

echo "send snp id's to snap"
for i in {1..24}
    do 
        curl --form snpListFile=@opensnp$i.txt \
            --form searchPairwise= \
            --form submit=search \
            --form snpFileColumn=1 \
            --form snpFileHeader=1 \
            --form hapMapRelease=onekgpilot \
            --form hapMapPanel=CEU \
            --form RSquaredLimit=0.8 \
            --form distanceLimit=500000 \
            --form downloadType=File \
            --form includeQuerySnp=on \
            --form arrayFilter=query \
            --form columnList[]=DP \
            --form columnList[]=GP \
            --form columnList[]=AM \
            http://www.broadinstitute.org/mpg/snap/ldsearch.php > SNAPResults$i.txt
    done

# remove all the opensnp..txt files
rm opensnp*.txt

# parses all the information back to the json schema
python ../parse_back_to_json.py ../$1 ../$2

# removes all the txt files of snap
rm SNAPR*.txt

# remove temporary directory
cd ..
rm -r tmp
