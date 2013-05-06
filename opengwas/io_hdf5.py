"""Interface with HDF5 file containing PLINK data.

The convention used in the container
Each probe is a row in the genotype array.
Each individual is a column in the genotype array.

import tables
raw_data = tables.openFile('test.h5', 'r')
raw_data.root.probes[1:10]['ID']
raw_data.root.individuals[3:6]
raw_data.root.genotype[:10,29:35]
raw_data.root.phenotype[:10]


To access in matlab, note that matrices are transposed:
probes = h5read('test.h5','/probes');
probes.ID(:,2:10)'
indiv = h5read('test.h5','/individuals');
indiv.individual(:,4:6)'
data = h5read('test.h5','/genotype')';
data(1:10,30:35)
labels = h5read('test.h5','/phenotype');
labels(1:10)
"""


import os
import tables
from numpy import genfromtxt
from numpy import array, flatnonzero

class GenotypeData(object):
    """Genotype data for genome wide association studies.
    Provides methods for reading, writing and slicing data.
    
    Container for genotype data in HDF5 format using PyTables.
    """
    def __init__(self, filename):
        self.file_name = filename
        self.data_name = os.path.basename(filename)
        self.h5_name = '%s.h5' % filename

        # The file pointer to the raw data
        self.h5_file = None

        # probes and individuals are similar to PLINK BIM and FAM files resp.
        self.probes = None
        self.individuals = None

        # phenotype and genotype points to a HDF5 carray
        self.genotype = None
        self.phenotype = None
        
        # Faster when using pytables
        #self.pytable_filters = tables.Filters(complevel=5, complib='blosc')
        # for portability and matlab compatibility 
        self.pytable_filters = tables.Filters(complevel=5, complib='zlib')

    @property
    def num_probes(self):
        """The number of loci tested, the number of genotypes per invididual"""
        return len(self.probes)

    @property
    def num_individuals(self):
        """The number of examples, the number of individuals in study"""
        return len(self.individuals)

    def open_file(self):
        """Open the HDF5 file and initialize variables"""
        self.h5_file = tables.openFile(self.h5_name, 'r')
        self.probes = self.h5_file.root.probes[:]
        self.individuals = self.h5_file.root.individuals[:]
        self.genotype = self.h5_file.root.genotype
        self.phenotype = self.h5_file.root.phenotype

    def init_file(self, overwrite=True):
        """Open a new HDF5 file for writing"""
        if (not overwrite) and os.path.isfile(self.h5_name):
            print('File %s found. Please remove manually.' % self.h5_name)
            return
        else:
            self.h5_file = tables.openFile(self.h5_name, 'w', title=self.data_name)

    def init_individuals(self, individuals):
        """Create the array on file to store the individuals.
        Also store the corresponding phenotype.
        TODO: typecheck input
        """
        self.h5_file.createTable(self.h5_file.root, 'individuals', individuals,
                                 title='Individuals', filters=self.pytable_filters)
        self.h5_file.root.individuals[:] = individuals
        self.individuals = self.h5_file.root.individuals[:]

    def init_phenotypes(self):
        """Extract the phenotype from the individuals"""
        #phenotypes = numpy.loadtxt(self.fam_name, usecols=[5], dtype=int)
        phenotypes = self.individuals['phenotype']
        phenotypes[phenotypes==1] = -1
        phenotypes[phenotypes==2] = 1
        phenotypes.shape = (len(phenotypes),1)
        #check for undefined phenotype
        undefined_phenotype = flatnonzero(phenotypes==0)
        if len(undefined_phenotype) > 0:
            print('Some phenotypes were undefined')
            print(undefined_phenotype)

        atom = tables.Int8Atom()
        self.h5_file.createCArray(self.h5_file.root, 'phenotype',
                                  atom, (1,self.num_individuals),
            title='Phenotype', filters=self.pytable_filters)
        self.h5_file.root.phenotype[:] = phenotypes.flatten()
        self.phenotypes = self.h5_file.root.phenotype

    def init_probes(self, probes):
        """Create the array on file to store the probe information.
        TODO: typecheck input
        """
        self.h5_file.createTable(self.h5_file.root, 'probes', probes,
                                 title='Probes', filters=self.pytable_filters)
        self.h5_file.root.probes[:] = probes
        self.probes = self.h5_file.root.probes[:]
        
    def init_genotypes(self):
        """Create the array on file (which is empty) to store genotypes."""
        atom = tables.Int8Atom()
        self.genotype = self.h5_file.createCArray(self.h5_file.root, 'genotype', atom,
                                                   (self.num_probes, self.num_individuals),
            title='Genotype', filters=self.pytable_filters)

    def close_file(self):
        """Close the HDF5 file."""
        self.h5_file.close()

    def get_genotype(self, idx_snp):
        """Return the array of genotypes for snp idx_snp
        for all individuals"""
        return array(self.genotype[idx_snp,:])
    
    def get_idx_case(self):
        """Return the index of individuals who are cases"""
        return flatnonzero(array(self.phenotype) == 1)

    def get_idx_control(self):
        """Return the index of individuals who are controls"""
        return flatnonzero(array(self.phenotype) == -1)

class GenotypeDataPlink(GenotypeData):
    """
    Class for converting data in PLINK bed, bim, fam files
    """
    def __init__(self, filename):
        """Just check that the files exists"""
        super(GenotypeDataPlink, self).__init__(filename)
        if not os.path.isfile('%s.bed' % filename):
            print('%s.bed not found' % filename)
            return
        if not os.path.isfile('%s.bim' % filename):
            print('%s.bim not found' % filename)
            return
        if not os.path.isfile('%s.fam' % filename):
            print('%s.fam not found' % filename)
            return
        self.bed_name = '%s.bed' % filename
        self.bim_name = '%s.bim' % filename
        self.fam_name = '%s.fam' % filename

    def plink2hdf5(self):
        """
        Load the genotypes from the plink bed file,
        and the phenotypes from the plink fam file.
        Load probe information from bim file,
        load individual information from fam file.
        """
        print('Creating %s' % self.h5_name)
        probes = self.load_probes()
        individuals = self.load_individuals()
        self.init_HDF5(individuals, probes)
        print('Number of Probes: %d' % self.num_probes)
        print('Number of Individuals: %d' % self.num_individuals)
        print('Reading features from %s' % self.bed_name)
        self.load_genotypes()
        self.h5_file.close()
        print('Compressed HDF5 file created at %s' % self.h5_name)

    def init_HDF5(self, individuals, probes):
        """
        Create an HDF5 file to contain the genotypes, phenotypes and probe information.
        Assumes that phenotype and probe information already loaded.
        """
        self.init_file()
        self.init_individuals(individuals)
        self.init_phenotypes()
        self.init_probes(probes)
        self.init_genotypes()

    def load_genotypes(self):
        """
        Load the plink BED format genotype data file.
        Assumes samples in columns and SNP loci in rows.

        Needs plinkio.
        https://github.com/fadern/libplinkio
        """
        from plinkio import plinkfile

        bed_file = plinkfile.open(self.file_name)
        for counter, row in enumerate(bed_file):
            self.genotype[counter,:] = list(row)
            if counter % 100000 == 99999:
                print(counter+1)
        bed_file.close()
        
    def load_probes(self):
        """Read the BIM file to get the probe locations
        chromosome (1-22, X, Y or 0 if unplaced)
        rs# or snp identifier
        Genetic distance (morgans)
        Base-pair position (bp units)
        """
        probes = genfromtxt(open(self.bim_name, 'r'), delimiter='\t',
        dtype={'names': ['chrom','ID','distance','bp_position','allele1','allele2'],
               'formats':[int, 'S16', int, int, 'S1', 'S1']})
        return probes

    def load_individuals(self):
        """Read the FAM file to get information about individuals
        Family ID
        Individual ID
        Paternal ID
        Maternal ID
        Sex (1=male; 2=female; other=unknown)
        Phenotype
        """
        individuals = genfromtxt(open(self.fam_name, 'r'), delimiter=' ',
        dtype={'names':['family','individual','paternal','maternal','sex','phenotype'],
               'formats':['S10','S16',int,int,int,int]})
        return individuals



def create_subset_probe(target, source, probe_idx):
    """Read the pytables HDF5 file source
    create a pytables HDF5 file target
    with the features with indices given by probe_idx.
    """
    orig_data = GenotypeData(source)
    orig_data.open_file()
    
    new_data = GenotypeData(target)
    new_data.init_file()
    new_data.h5_file.copyNode(orig_data.h5_file.root.individuals, newparent=new_data.h5_file.root)
    new_data.individuals = new_data.h5_file.root.individuals[:]
    new_data.init_phenotypes()
    probes = orig_data.probes[probe_idx]
    new_data.init_probes(probes)
    new_data.init_genotypes()
    new_data.h5_file.root.genotype[:] = orig_data.genotype[probe_idx,:].copy()

    orig_data.close_file()
    new_data.close_file()
    
def summary(file_name, peek_ind=5, peek_prb=10):
    """Report some statistics about the file"""
    data = GenotypeData(file_name)
    data.open_file()
    print('Number of Probes: %d' % data.num_probes)
    print('Number of Individuals: %d' % data.num_individuals)
    print('The genotype matrix is of size %d by %d' % data.genotype.shape)
    print('First %d individuals' % peek_ind)
    print(data.individuals[:peek_ind])
    print('First %d probes' % peek_prb)
    print(data.probes[:peek_prb])
    print('The corresponding genotypes')
    print(data.genotype[:peek_prb,:peek_ind])
    data.close_file()
    
