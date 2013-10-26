"""
Read the results of GWES (http://bioinfo.utu.fi/GWESserver/) from gwes_file
Convert into the graph format stored in a JSON file
Write to rede_file

"""

import sys
import json
import numpy
from numpy import genfromtxt, unique, flatnonzero
from numpy import array, delete, in1d

class ProbeGraph(object):
    """A container for representing epistatic interaction graphs"""
    def __init__(self, name):
        self.name = name
        self.probes = []
        self.pairs = []

    def get_nodes(self, data):
        """Populate the nodes"""
        self.all_probes = unique(data['id1'].tolist()+data['id2'].tolist())
        bad_probes = []
        for idx,elem in enumerate(self.all_probes):
            prb = 1
            row = flatnonzero(elem == data['id1'])
            if len(row) == 0:
                prb = 2
                row = flatnonzero(elem == data['id2'])
            node = {}
            chrom = int(data['chrom%d' % prb][row][0])
            if chrom < 1 or chrom > 24:
                bad_probes.append(data['id%d' % prb][row][0])
                continue
            node['chrom'] = chrom
            node['bp_position'] = int(data['pos%d' % prb][row][0])
            node['prb'] = idx+1
            node['id'] = idx+1
            node['prbCode'] = elem
            node['rs'] = elem
            node['degree'] = len(flatnonzero(elem == data['id1']))+len(flatnonzero(elem == data['id2']))
            node['probe_group'] = 1
            self.probes.append(node)
        self.bad_probes = array(bad_probes)
        self.all_probes = delete(self.all_probes, flatnonzero(in1d(self.all_probes, self.bad_probes)))

    def get_links(self, data):
        """Populate the links"""
        for idx,row in enumerate(data):
            if row['id1'] in self.bad_probes or row['id2'] in self.bad_probes:
                print('Probe unknown')
                print(row)
                continue
            link = {}
            link['source'] = int(flatnonzero(row['id1'] == self.all_probes)[0])
            link['target'] = int(flatnonzero(row['id2'] == self.all_probes)[0])
            link['probe_group'] = int(row['genclass'])
            link['ct_id'] = idx+1
            link['gwes'] = float(row['score'])
            self.pairs.append(link)

    def report_stats(self):
        """Summary of the graph"""
        print('Graph: %s' % self.name)
        print('%d nodes, %d links' % (len(self.probes), len(self.pairs)))
        print('%d unmapped probes' % len(self.bad_probes))
        
    def save_json(self, rede_file):
        """Export to rede_file"""
        handle = open(rede_file, 'w')
        handle.write('{"name": "%s", ' % self.name)
        handle.write('"nodes": ')
        handle.write(json.dumps(self.probes))
        handle.write(', "links": ')
        handle.write(json.dumps(self.pairs))
        handle.write('}')
        handle.close()
            
    def load_gwes(self, gwes_file):
        """
        Read the results of GWES (http://bioinfo.utu.fi/GWESserver/) from gwes_file
        GWES file has columns:
        SNP1_ID SNP2_ID SNP1_pValue     SNP2_pValue     p-value p-value(-log10) 
        GenoClassNo SNP1_Chr        SNP2_Chr        SNP1_Position   SNP2_Position
        """
        data = genfromtxt('SNPpairs.txt', skiprows=1,
                          dtype={'names':['id1','id2','uni1','uni2','pval','score',
                                          'genclass','chrom1','chrom2','pos1','pos2'],
                                 'formats':['S15','S15',float,float,float,float,
                                            int,int,int,int,int]})
        self.get_nodes(data)
        self.get_links(data)
        
        

def gwes2rede(gwes_file, rede_file):
    snp_graph = ProbeGraph(gwes_file)
    snp_graph.load_gwes(gwes_file)
    print('File loaded with ...')
    snp_graph.report_stats()
    print('Writing to %s' % rede_file)
    snp_graph.save_json(rede_file)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python %s gwes_file.txt rede_file.json' % sys.argv[0])
        exit(1)
    gwes2rede(sys.argv[1], sys.argv[2])
