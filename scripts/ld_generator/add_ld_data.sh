# Author: Stefan Sevelda
# 1: parse json file and filter all rs numbers in seperate files for the
# chromsomes
python generatsnap.py $1

# calculate the lds for the rede json file 

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
python parse_back_to_json.py $1 $2

# removes all the txt files of snap
rm SNAPR*.txt
