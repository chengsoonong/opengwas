
import os, subprocess
import tables
from opengwas.io_hdf5 import GenotypeDataBase, GenotypeDataPlink, GenotypeData
from numpy import array, flatnonzero

def create_hdf5():
    """Helper code to create HDF5 file"""
    dataname = 'test/plink_test'
    test_create = GenotypeDataPlink(dataname)
    test_create.plink2hdf5()

def test_convert():
    """Convert PLINK bed, bim and fam files into HDF5,
    and check that we can read slices.
    """
    dataname = 'test/plink_test'
    test_create = GenotypeDataPlink(dataname)
    test_create.plink2hdf5()

    # The direct interface via PyTables
    print('Test read')
    test_read = tables.openFile(dataname+'.h5', 'r')
    print(test_read.root.genotype[2:6, :10])
    print(test_read.root.phenotype[0,:10])
    print(test_read.root.probes[3:7]['ID'])
    print(test_read.root.individuals[:10]['phenotype'])
    test_read.close()

    subprocess.check_call('rm test/plink_test.h5', shell=True)


def test_interface():
    """Check that the interfaces work"""
    create_hdf5()

    # The recommended way to access the data
    print('Test read interface in GenotypeDataBase')
    test_face = GenotypeData('test/plink_test')
    test_face.open_file()
    print(test_face.genotype[2:6, :10])
    print(test_face.phenotype[0,:10])
    print(test_face.probes[3:7]['ID'])
    print(test_face.individuals[:10]['phenotype'])

    print('Test interface in GenotypeData')
    print(test_face.get_idx_case())
    print(test_face.get_idx_control())

    print('Example contingency table calculation')
    print('Number of individuals = %d' % test_face.num_individuals)
    idx_snp = 42
    geno_vec = test_face.get_genotype(idx_snp)
    print(geno_vec)
    idx_cases = test_face.get_idx_case()
    idx_controls = test_face.get_idx_control()
    for genotype in [0,1,2]:
        for group in [idx_cases, idx_controls]:
            count = len(flatnonzero(geno_vec[group] == genotype))
            print(count)
        count = len(flatnonzero(geno_vec == genotype))
        print('Total=%d' % count)
    test_face.close_file()

    subprocess.check_call('rm test/plink_test.h5', shell=True)


if __name__ == '__main__':
    test_interface()
