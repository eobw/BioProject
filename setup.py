#!/usr/bin/env python3.6

from setuptools import setup, find_packages
import sys, os

setup(
    name='GUESSmyLT',
    version='0.1',

    description='An efficient way to guess the library type of RNA-reads',

    url='https://github.com/eobw/BioProject/',
    author='Hampus Olin, Erik Berner-Wik, Caitlin Haughey',

    license='GPL-3.0',
    packages=find_packages(),

    install_requires=['biopython>=1.67', 'bcbio-gff>=0.6.4'],
    include_package_data=False,

    entry_points={
        'console_scripts': ['GUESSmyLT = GUESSmyLT.GUESSmyLT:main',
        ],
    }
)

print("Trying to download BUSCO datasets from https://busco.ezlab.org")
print("Checking for eukaroyte dataset:")
if os.path.exists("data/eukaryota_odb9"):
    print("Eukatyote dataset seems to exist")
else:
    print("No eukaryote dataset found, downloading...")
    os.system("wget -qO- https://busco.ezlab.org/datasets/eukaryota_odb9.tar.gz | tar xvz -C ./data")
print("Checking for prokaryote dataset:")
if os.path.exists("data/bacteria_odb9"):
    print("Prokaryote dataset seems to exist")
else:
    print("No prokaryote dataset found, downloading...")
    os.system("wget -qO- https://busco.ezlab.org/datasets/bacteria_odb9.tar.gz | tar xvz -C ./data")