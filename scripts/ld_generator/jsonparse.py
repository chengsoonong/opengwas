import json
from pprint import pprint

def snpdict(n, file):
    json_data=open(file)
    id_json={}
    data = json.load(json_data)
    for a in data["nodes"]:
        if a["chrom"] == n:
            id_json[a["rs"]] = a["id"] 
    json_data.close();
    return id_json
