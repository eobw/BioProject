# GUESSmyLT
Software for Linux to guess the RNA-Seq library type of paired and single end read files using mapping and gene annotation.  

## Table of contents

* [GUESSmyLT](#GUESSmyLT)
* [Background](#background)
* [Prerequisites](#Prerequisites)
* [Installation](#Installation)
  * [Installation with git](#installation-with-git)
  * [Check installation](#check-installation)
* [Usage](#usage)
  * [File formats](#file-formats)
  * [Supported header formats](#supported-header-formats)
  * [Supported interleaved formats](#supported-interleaved-formats)
  * [Example commands](#example-commands)
  * [Output](#output)
  * [Parameters](#parameters)
* [Known issues](#known-issues)
* [TO DO](#to-do)

## Background
The choice of RNA-Seq library type defines the read orientation of the sequencing and the order in which the strands of cDNA are sequenced, which means that RNA-Seq reads from different library types can differ significantly. The information regarding library type can be very useful for reads to be assembled into a transcriptome or mapped to a reference assembly. This is because the library type can help to discern where in the transcriptome shorter ambiguous reads belong by using the read’s relative orientation and from which strand it was sequenced. Unfortunately, this information regarding the library type used is not included in sequencing output files and is usually lost before the assembly of the data. Even when working with RNA-Seq data from public repositories there is no guarantee that the library type information is correct or that it exists at all. This is what GUESSmyLT aims to fix by looking at how reads map to a reference and together with gene annotation guess which library was used to generate the data.

## Prerequisites
Only developed for Linux systems.
Python 3:
  - biopython (1.67)
  - bcbio-gff (0.6.4)
  - pysam (0.15.1)

Other programs:
  - Trinity (2.8.4) - Reference assembly
  - BUSCO (3.0.2) - Gene annotation
  - Bowtie2 (2.3.4.3) - Mapping
  - Snakemake (5.4.0) - Workflow management

## Installation
Installation using setup.py is not necessary, it adds the ability to run the tool using the command 'GUESSmyLT' globally. Tool can be used by running the python script 'GUESSmyLT.py' from within the tool folder.

#### Installation with git:

Clone the repository to your home directory:

```bash
cd ~
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
Read files: .fastq
Mapping:    .bam
Reference:  .fa


### Supported header formats
Tested for Old/New Illumina headers and downloads from SRA.
Should work, but not tested for all fastq header formats at: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/

### Supported interleaved formats
If headers are in Old/New Illumina or if reads are alternating.  
Add example..

### Example commands

#### Paired end reads and reference with specified subsampled reads
```bash
GUESSmyLT --reads /home/.../read_1.fastq /home/.../read_2.fastq --organism pro --reference /home/.../ref.fa --subsample 100000
```

#### Only paired end reads eukaryotic
```bash
GUESSmyLT --reads /home/.../read_1.fastq /home/.../read_2.fastq --organism euk
```

#### Single end reads and reference
```bash
GUESSmyLT --reads /home/.../reads.fastq --organism pro --reference /home/.../ref.fa
```

#### Without installed single end reads and reference
```bash
cd GUESSmyLT/
python3 GUESSmyLT.py --reads /home/.../reads.fastq --organism pro --reference /home/.../ref.fa
```
### Output
GUESSmyLT will print the result in the command line as well as write it to a file:
```bash
GUESSmyLT/data/output/[read_name]_result.txt
```
Results from intermediate steps, such as the mapping from Bowtie2 or annotation from BUSCO are saved in
```bash
GUESSmyLT/data/intermediate/
```

## Parameters

### Mandatory

| Parameter | Input | Description |
| --- | --- | --- |
| --reads | .fastq file(s) | Full path(s) to RNA-Seq read file(s). Can be compressed or uncompressed. Order is not important. Can handle two paired end read files, one interleaved read file and single end read file. |
| --organism | euk or pro | Eukaryote or prokaryote (euk/pro) is an option needed for the BUSCO annotation. |



### Optional
| Parameter | Input | Description |
| --- | --- | --- |
| --subsample | Even integer | Number of reads that will be used for subsampling. |
| --reference | .fa file | Full paths to reference genome/transcriptome for mapping reads to. |
| --threads | Integer | Number of threads to use. |
| --memory | Number of GB ex: 10G | Maximum memory that can be used in GB. |
| --annotation | .gff file | Full path to annotation file for skipping BUSCO step. NOT DEVELOPED YET. |
| --mapped | Sorted .bam file | Full path to mapped read file for skipping Bowtie2 step. NOT DEVELOPED YET. |
| --output | File path | Full path to result file. NOT DEVELOPED YET. |

## Known issues
1) Complains about gzip broken pipe when subsampling with compressed files (but works anyway).  
2) BUSCO sometimes looses the config path. Fix manually in terminal:
```bash
export AUGUSTUS_CONFIG_PATH=~/miniconda3/pkgs/augustus-3.2.3-boost1.60_0/config
```
3) BUSCO might not find any core genes. Fix by using more reads or by providing reference.
4) Mapping, annotation, assembly or the entire pipeline is skipped.  This is most likely due to the fact that Snakemake checks which output files need to be generated and from there only performs the necessary steps of the pipeline. The result of this is that is you already have a .bam file, BUSCO/Trinity output folder or a result .txt file for the reads Snakemake will skip steps

## TO DO
1) Handling of --output  
2) Handling of --mapped
3) Handling of --annotation  
4) Subsampling only takes top n reads in .fastq files. Can be improved by selecting reads randomly.  
