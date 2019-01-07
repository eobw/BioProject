# Top script for pipeline.
# Execute with: python execute.py --reads read1 (read2) --organism pro/euk
# For more options type: python execute.py --help





# ADD CHECKER FOR IF FILES EXIST!!!!!!!!!!!!!!!!!!!!!!!

import argparse, sys, json, os, re, subprocess
from BCBio import GFF
import pprint

parser=argparse.ArgumentParser(description="Predict library type.")
parser._action_groups.pop()
required=parser.add_argument_group("required arguments")
required.add_argument("--organism",type=str)
optional=parser.add_argument_group("optional arguments")
required.add_argument("--reads",nargs="+",type=str,help="Reads in (un)zipped .fastq format.")
optional.add_argument("--reference",type=str,help="Reference genome/transcriptome in .fasta format.")
optional.add_argument("--annotation",type=str,help="Annotation in .gff format")
optional.add_argument("--mapped",type=str,help="") #Maybe add this?
# maybe add group..
#group=parser.add_mutually_exclusive_group(required=True)
#group.add_argument("--annotation",type=str,help="Annotation in .gff format")
#group.add_argument("--organism",type=str)
optional.add_argument("--output",type=str,help="...")
# These are not added to config.json yet.. 
optional.add_argument("--threads", type=int, default=10,help="")
optional.add_argument("--memory",type=str, default="8G",help="Maximum memory that can be used in GB. Ex. '10G'.")


args=parser.parse_args()

# Maybe add exception class

def check_organism():
    if args.organism.lower() in "prokaryote":
        args.organism="prokaryote"
    
    elif args.organism.lower() in "eukaryote":
        args.organism="eukaryote"
    else: 
        print("Error. Unrecognized organism --organism "+args.organism+". Only eukaryote/prokaryote are valid organism.")
        sys.exit()

def check_annotation(annotation_file):
    in_handle=open(annotation_file)
    for record in GFF.parse(in_handle, limit_info=dict(gff_type=["gene"])):
        if record:
            return 
    # If no genes are found, gff file cannot be used in analysis.
    print("Error. No genes could be found in --annotation "+annotation_file+". Please, submit a .gff file containing genes or no .gff file.")
    sys.exit()

def check_single_reads(read):
    if not any(x in read for x in [".fq",".fastq"]):
        print("Error. Cannot find .fq or fastq extension in read file --read "+read+".")
        sys.exit()
    checked_reads=subprocess.check_output(
    "./check_single_read_files.sh "+
    read,
    shell=True,
    encoding="utf8"
    )
    if len(checked_reads.split("\n"))==1:
        # Do nothing, data file is single end.
        return [read]
    else:
        # Interleaved paired end reads. 
        # They have been deinterleaved with check_single_read_files.sh
        # Store new, deinterleaved readfiles in args.reads:
        return checked_reads.split("\n")[1:2]

            
def check_paired_reads(read1,read2):
    if not any(x in read1 for x in [".fq",".fastq"]) and not any(x in read2 for x in [".fq",".fastq"]):
        print("Error. Cannot find .fq or .fastq extension in read files --read "+read1+" "+read2+".")
        sys.exit()
    ordered_reads=subprocess.check_output(
    "./check_paired_read_files.sh "+
    read1+" "+
    read2,
    shell=True,
    encoding="utf8")
    if len(ordered_reads.split("\n"))==2:
        print("Error. Unrecognized header format for paired end reads.")
        sys.exit()
    else:
        return ordered_reads.split("\n")[1:3]

def check_mapped():
    print("Checker for mapper has not been developed yet.")
    sys.exit()
    
def check_reference(ref):
    if ".gz" in ref:
        zip_command="zcat"
    else:
        zip_command="cat"
    num_headers=subprocess.check_output(
    zip_command+" "+ref+" | grep '^>' | wc -l",
    shell=True,
    encoding="utf8")
    if int(num_headers.split("\n")[0]) > 1:
        # Multiple headers --> dealing with transcriptome.
        # So far we return genome anyway because we have not optimized busco for transcriptome.
        # Need to optimize busco script!
        # This is however a later improvement since we can lie to busco that we are dealing with genome.
        #return "transcriptome"
        return "genome"
    elif int(num_headers.split("\n")[0])==1:
        # One header --> dealing with genome.
        return "genome"
    else:
        print("Error. Reference file is not in fasta format. Missing '>' in beginning of fasta header.")
        sys.exit()

def check_memory():
    if "G" == args.memory[-1] and "0" != args.memory[0] and args.memory[:-1].isdigit():
        return True
    else:
        print("Invalid --memory argument: '"+args.memory+"'. Memory should be given in GIGABYTES. For example --memory '4G'")
        sys.exit()

# ---Main From Here ---

# ---Check reads---
# At least one readfile must be provided. (single end or interleaved paired end reads)
# At most two readiles can be provided. (paired end reads)
#
# Check if readfile is provided.
if not args.reads:
    print("Error. No read files provided.")
    sys.exit()
# Check if one read file is provided.
elif len(args.reads)==1:
    # Check that single reads are OK
    # or deinterleave paired end reads. 
    args.reads=check_single_reads(args.reads[0])
# Check if two read files are provided.    
elif len(args.reads)==2:
    # Check that paired end reads are OK.
    args.reads[0],args.reads[1]=check_paired_reads(args.reads[0],args.reads[1])    
# Too many readfiles were provided.    
else:
    print("Error. Too many read files.")
    sys.exit()

# After reads have been checked, they are sorted so that read[1/left/forward] is first and read[2/right/reverse] is second.
# Base readname on first read: 
# Get readname, i.e: readname from path/to/reads/[readname].fq.gz
readname=re.split(".fq|.fastq",os.path.basename(args.reads[0]))[0]

    
# ---Check threads---
#if args.threads:
#    check_threads()
#else:
#    # Default threads
#    args.threads=6


# ---Check mapping file---
# If a mapping file (.bam/.sam) is provided, the reference used for mapping 
# must also be provided.
if args.mapped and not args.reference:
    print("Error. If a mapping file is provided, the reference used for mapping must also be provided.")
    sys.exit()

if args.mapped:
    check_mapped()
else:
    args.mapped="../data/intermediate/"+readname+"_mapped.bam"
 
if args.reference:
    busco_reference_mode=check_reference(args.reference)
else:
    # Add Trinity assembly to config file.
    args.reference="../data/intermediate/"+readname+".fasta"

busco_reference_mode="genome"
# ---Check organism and/or annotation---
# Organism or annotation file must be provided.
# If none is provided, we cannot look for core genes.
#if args.organism:
#    check_organism()
#elif args.annotation:
#    check_annotation()
#else:
#    print("Error. Provide organism or annotation file.")
#    sys.exit()

if args.annotation:
    check_annotation(args.annotation)
else:
    args.annotation="../data/intermediate/run_"+readname


check_memory()


# add checker for fasta files

# add memory and threads options

# Change config file

with open("config.json","r+") as configfile:
    data=json.load(configfile)
    data["input"]["reference"]=args.reference
    data["input"]["reads"]=args.reads
    data["input"]["readname"]=readname
    data["input"]["annotation"]=args.annotation
    data["input"]["organism"]=args.organism
    data["input"]["mapped-reads"]=args.mapped
    data["busco"]["lineage"]=("../data/bacteria_odb9" if args.organism == "prokaryote" else "../data/eukaryota_odb9")
    data["busco"]["mode"]=busco_reference_mode
    data["input"]["threads"]=args.threads
    data["input"]["memory"]=args.memory
    configfile.seek(0)
    json.dump(data,configfile,indent=4)
    configfile.truncate()
    
# Add execution for snakefile...
#os.system("snakemake ../data/intermediate/trinity --dag | dot -Tsvg > dag.svg -forceall")
#os.system("snakemake ../data/output/result_4.txt")
#os.system("snakemake -s bowtie2 ../data/intermediate/"+readname+"_sorted.bam --forceall")
#os.system("snakemake -s trinity --cores "+str(args.threads))
os.system("snakemake -s busco ../data/intermediate/run_"+readname+" --forceall")