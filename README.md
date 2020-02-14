usable_vcf - Run multiple programs to check if a VCF is usable or not

Requirements
============

The following programs need to be in your PATH

- bcftools
- vcftools
- tabix

The following libraries need to be installed in your python environment

- pip install pysam
- pip install pyvcf

How to run
==========

usage: usable_vcf [-h] [-v] [-m MAX_ENTRIES] vcf

Runs commands and VCF parsing libraries to check if a vcf is usable or not

positional arguments:
  vcf                   vcf to parse

  optional arguments:
    -h, --help            show this help message and exit
    -v, --verbose         Output errors returned by commands/libraries
    -m MAX_ENTRIES, --max-entries MAX_ENTRIES
    			  Maximum number of entries to parse (10000)

Philosophy
==========
There are many different programs/libraries that intake VCFs. However, many tools that output VCFs that
aren't usable by a large number of tools. This simple script gives a report of how many common tools/use-cases
a VCF would work with. Passing all of these checks requires strict adherence to the VCF specification.

This project is in beta. Please add tickets or pull-requests if you'd like to see other tools and use-cases added
to the suite of those currently checked. 

If you believe a VCF is incorrecly failing a test, please add a ticket. Our checks may be invalid or we may have
an interesting discussion on how to interpret the VCF specification.
