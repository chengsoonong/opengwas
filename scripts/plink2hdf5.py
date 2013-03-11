"""Convert PLINK bed, bim and fam files into HDF5."""

import os
from opengwas.io_hdf5 import GenotypeDataPlink

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print 'Usage: python %s <name of plink file>' % sys.argv[0]
        exit(1)
    data = GenotypeDataPlink(sys.argv[1])
    data.plink2hdf5()
