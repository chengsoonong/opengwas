========
opengwas
========

Converting between Genome Wide Association Study (GWAS) file formats

Currently supported formats
==================

openSNP
-------

* Reading using https://github.com/superbobry/snpy


PLINK
-----

http://pngu.mgh.harvard.edu/~purcell/plink/

* Writing PED and MAP files
* Reading PED, BIM, FAM files using https://github.com/fadern/libplinkio


HDF5
----

Convenient for python and matlab

* Reading and Writing using http://www.pytables.org



Create the files in format plink from openSNP data
==================================================


Using scripts in folder scripts/

Generate BED, BIM, FAM file from openSNP data for phenotypes of interest
./process_opensnp.sh

Use plink to find epistatic interaction pairs.
Also create some files that are useful for downstream processing
./plink_example.sh

