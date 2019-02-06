#!/usr/bin/env python3.6

from setuptools import setup, find_packages
import sys, os, json

print("Trying to download BUSCO datasets from https://busco.ezlab.org")
print("Checking for eukaroyte dataset:")
if os.path.exists("GUESSmyLT/eukaryota_odb9"):
    print("Eukaryote dataset seems to exist")
else:
    print("No eukaryote dataset found, downloading...")
    os.system("wget -qO- https://busco.ezlab.org/datasets/eukaryota_odb9.tar.gz | tar xvz -C ./GUESSmyLT")
print("Checking for prokaryote dataset:")
if os.path.exists("GUESSmyLT/bacteria_odb9"):
    print("Prokaryote dataset seems to exist")
else:
    print("No prokaryote dataset found, downloading...")
    os.system("wget -qO- https://busco.ezlab.org/datasets/bacteria_odb9.tar.gz | tar xvz -C ./GUESSmyLT")

script_dir=os.getcwd() + "/"
config_path = script_dir+"GUESSmyLT/config.json"
with open(config_path,"r+") as configfile:
    data=json.load(configfile)
    data["eukaryote_db"]=script_dir+"GUESSmyLT/eukaryota_odb9"
    data["prokaryote_db"]=script_dir+"GUESSmyLT/bacteria_odb9"
    configfile.seek(0)
    json.dump(data,configfile,indent=4)
    configfile.truncate()

setup(
    name='GUESSmyLT',
    version='0.1',

    description='An efficient way to guess the library type of RNA-reads',

    url='https://github.com/eobw/BioProject/',
    author='Hampus Olin, Erik Berner-Wik, Caitlin Haughey',

    license='GPL-3.0',
    packages=find_packages(),

    install_requires=['biopython==1.67', 'bcbio-gff==0.6.4', 'pysam==0.15.1', 'snakemake==5.4.0', 'busco==3.0.2'], #Can't get bowtie2 dependency to work
    dependency_links=['https://github.com/bioconda/bioconda-recipes/tree/master/recipes//busco/meta.yaml', 'https://github.com/bioconda/bioconda-recipes/tree/master/recipes//bowtie2/meta.yaml'],
    include_package_data=True,

    entry_points={
        'console_scripts': ['GUESSmyLT = GUESSmyLT.GUESSmyLT:main',
        ],
    }
)
