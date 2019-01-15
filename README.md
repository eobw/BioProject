# BioProject

## Prerequisites
Python 3, biopython, bcbio-gff, pysam, Trinity, BUSCO, Bowtie2, Snakemake

## Installation
Installation is not necessary, it adds the ability to run the tool using the command 'GUESSmyLT' globally. Tool can be used by running the python script 'GUESSmyLT.py' from within the tool folder.

#### Installation with git:

Clone the repository:

```bash
git clone https://github.com/eobw/BioProject.git
```

Move to the folder and install:

```bash
cd BioProject/
python setup.py install
```

#### Check installation

Executing:
```bash
GUESSmyLT
```

or

```bash
GUESSmyLT -h
```

to display help.

## Usage

### File formats


### Supported header formats

### Example commands

#### Paired end reads and reference
```bash
GUESSmyLT --reads /home/.../read_1.fastq /home/.../read_2.fastq --organism pro --reference /home/.../ref.fa
```

#### Only paired end reads eukaryotic
```bash
GUESSmyLT --reads /home/.../read_1.fastq /home/.../read_2.fastq --organism euk
```

#### Single end reads and reference
```bash
GUESSmyLT --reads /home/.../reads.fastq --organism pro --reference /home/.../ref.fa
```

## Parameters

### Mandatory

| Parameter | Input | Description |
| --- | --- | --- |
| --reads | .fastq file(s) | Full path(s) to RNA-Seq read file(s). For paired end /1 first and /2 second. Order is important. If one file given single end is presumed. |
| --organism | euk or pro | Eukaryote or prokaryote (euk/pro) is an option needed for the BUSCO anotation. |



### Optional
| Parameter | Input | Description |
| --- | --- | --- |
| --reference | .fa file | Full paths to reference genome/transcriptome for mapping reads to. |
| --annotation | .gff file | Full path to annotation file for skipping BUSCO step. |
| --threads | Integer | Number of threads to use. |
| --memory | Number of GB ex: 10G | Maximum memory that can can be used in GB. |

## Known issues
Spaces in headers
