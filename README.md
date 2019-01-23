# GUESSmyLT
Software for Linux to guess the RNA-Seq library type of paired and single end read files using mapping and gene annotation.  

## Table of contents

* [GUESSmyLT](#GUESSmyLT)
* [Background](#background)
* [Overview of pipeline](#Overview of pipeline)
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

## Overview of pipeline
![alt text](https://github.com/eobw/BioProject/blob/master/Overview_of_GUESSmyLT.png "Pipeline of GUESSmyLT")  
GUESSmyLT uses Snakemake to build the pipeline it needs in order to predict the library type. Required arguments are organism (euk/pro) and reads (read file(s) in fastq format). Reference (genome or transcriptome in .fasta format) is optional, and if it is not provided, Trinity will first be executed to create a De novo assembly of the reads. Next, BUSCO is used for annotation. This is also a QC step because BUSCO looks for core genes, so called BUSCOs, in the reference. If they cannot be found, it indicates that the reference has bad quality and therefore the pipeline will terminate. If BUSCOs are found, the process continues with mapping the reads to the reference using Bowtie2. The mapping is done with unstranded option so that the reads can be mapped on both the strands and in both directions. Finally, the mapping and annotation is used for inference, which is done with a python script and the library type is returned. 
On top of Snakemake, we have a python script, GUESSmyLT.py. Its purpose is to handle user arguments by:
1.	Checking that arguments are correct, files exists and are in correct format.
2.	Subsamples reads into new read files that are used in the analysis. This makes GUESSmyLT faster and protects the original files from being modified.
3.	Modifying files:
a.	Changes read files that are in wrong format. Trinity and Pysam can only handle old Illumina format: @read_ID/pair#, where pair# is 1 or 2. They do not work with whitespaces, punctutations nor undescrores. Therefore, the script makes sure that the headers are converted into the correct format.
b.	Deinterleaves paired end read files if they are interleaved.
4.	Telling Snakemake what files exist by updating the config file.
5.	Executing snakemake.

## Result
The results are printed as stdout and to a result file. One example of a result would be:
```
Results of library inferring fruit_fly: 
Library type 	 Reads 	 Percent 
fr_first 	     4019 	 47.2%
fr_second 	    4454 	 52.3%
rf_first 	     21 	   0.2%
rf_second 	    19 	   0.2%
ff_first 	     5 	    0.1%
ff_second 	    2 	    0.0%
undecided 	    1 	    0.0%
```
Based on the orientations of the reads we would assume that the library type is fr-unstranded as there is roughly a 50-50 split between fr-first and fr-second.

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

Old Illumina: @HWUSI-EAS100R:6:73:941:1973#0/1  
New Illumina: @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG  
Alternating:  
&nbsp; &nbsp; &nbsp; @read1 (first mate)  
&nbsp; &nbsp; &nbsp; ..  
&nbsp; &nbsp; &nbsp; @read1 (second mate)  
&nbsp; &nbsp; &nbsp; ..  
&nbsp; &nbsp; &nbsp; @read2 (first mate)  
&nbsp; &nbsp; &nbsp; ..   
&nbsp; &nbsp; &nbsp; @read2 (second mate)  
&nbsp; &nbsp; &nbsp; ..  

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
5) For Linux subsystem on windows, the correct path to GUESSmyLT have to be changed manually. This is done by changing the path of the variable 'script_dir' to the full path of GUESSmyLT. For example: script_dir="/mnt/c/Users/erik_/Documents/GitHub/BioProject/GUESSmyLT"  

## TO DO
1) Handling of --output  
2) Handling of --mapped  
3) Handling of --annotation  
4) Subsampling only takes top n reads in .fastq files. Can be improved by selecting reads randomly.  
5) Make BIOCONDA package for easy access.  
6) Add test data that check that updates of GUESSmyLT are correct.  
7) Write Wiki for developers that want to help with this project.  
8) Optimize BUSCO for transcriptome reference. It is now implemented so that BUSCO always analyses reference with mode 'genome', because otherwise BUSCO genes will not be returned.  
9) Look more into why some reads get undecided orientation. This is when a read's mate cannot be found and is probably due to a read is at the end of a gene and its mate is outside of the selected region.  
10) Logging of last step, Inference.  

## Citing 
If you use GUESSmyLT in your work, please cite us:  
Add DOI: Berner Wik, E., Olin, H., Vigetun Haughey C.
