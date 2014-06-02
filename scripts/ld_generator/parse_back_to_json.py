import json 
import jsonparse
import sys

infile = sys.argv[1]
outputfile = sys.argv[2]

a = 1
json_data=open(infile)
data = json.load(json_data)
json_data.close()

handle = open(outputfile, 'w')
handle.write('{"name": "asthma", "nodes": ')
handle.write(json.dumps(data["nodes"]))
handle.write(', "links": ')
handle.write(json.dumps(data["links"]))
handle.write(', "cont_table": ')
handle.write(json.dumps(data["cont_table"]))
handle.write(', "ld_plot": [{')
while a < 25:
    if a == 1:
        handle.write('"local_links%d": ' %(a,))
    else:
        handle.write(', "local_links%d": ' %(a,))
    file_object = open('SNAPResults%d.txt' % (a,), "r")
    lines = file_object.readlines();
    file_object.close()
    lines = lines[1:]

    snparray = list()
    snpdict = dict()
    snpdict = jsonparse.snpdict(a, infile) 

    for i in lines:
        thisline = i.split();
        linkage_dis  = thisline[0:4] + thisline[6:] 
        if linkage_dis[1] != "WARNING" and linkage_dis[2] != "0":
            linkage_dis[0] = snpdict.get(linkage_dis[0]) or False 
            linkage_dis[1] = snpdict.get(linkage_dis[1]) or False 
            linkage_dis[2] = int(linkage_dis[2]);
            linkage_dis[3] = float(linkage_dis[3]);
            linkage_dis = { 'source': linkage_dis[0], 'target':
                    linkage_dis[1], 'distance': linkage_dis[2], 'R2':
                    linkage_dis[3], 'chrom': a}
            if linkage_dis['target'] != False:
                snparray.append(linkage_dis)
    handle.write(json.dumps(snparray))
    a = a+1
handle.write('}]}')
handle.close()

