import json
import jsonparse 
import sys
from pprint import pprint

def snp_rs_detection(n):
    json_data=open(infile)
    rsnumbers = []
    data = json.load(json_data)
    json_data.close();
    for a in data["nodes"]:
        if a["chrom"] == n:
            rsnumbers.append(a["rs"])
    pprint(rsnumbers)
    file = 'opensnp%d.txt' % (n,)
    handler = open(file, 'w')
    for i in rsnumbers:
        handler.write(i + '\n')
    handler.close()

infile = sys.argv[1]

for i in range(24):
        snp_rs_detection(i+1)



    
    
