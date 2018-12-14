# Top script for pipeline.
# Execute with: python execute.py --reads read1 (read2) --organism pro/euk
# For more options type: python execute.py --help




import argparse, sys, json, os, re, subprocess

parser=argparse.ArgumentParser(description="Predict library type.")
parser.add_argument("--reads",nargs="+",help="Reads in (un)zipped .fastq format.")
parser.add_argument("--reference",help="Reference genome/transcriptome in .fasta format.")
parser.add_argument("--annotation",help="Annotation in .gff format")
parser.add_argument("--mapped",help="") #Maybe add this?
parser.add_argument("--output",help="...")
parser.add_argument("--organism")
# These are not added to config.json yet.. 
parser.add_argument("--cpu", help="")
parser.add_argument("--memory",help="")


args=parser.parse_args()

#read_basename=os.path.basename(args.reads[0])
#read_basename_no_extension=re.split(".fq|.fastq",read_basename)[0]
#if len(re.split(".fq|.fastq",read_basename))!=2:
#    print("Error. Reads are not in (zipped) fastq format.")
#    sys.exit()

def check_organism():
    if args.organism.lower() in "prokaryote":
        args.organism="prokaryote"
    
    elif args. organism.lower() in "eukaryote":
        args.organism="eukaryote"
    else: 
        print("Error. Unrecognized organism.")
        sys.exit()

def check_annotation():
    print("Checker for annotation file is developed yet.")
    print("Exiting...")
    sys.exit()
    
def check_single_reads():
        checked_reads=subprocess.check_output(
        "./check_single_read_files.sh "+
        args.reads[0],
        shell=True,
        encoding="utf8"
        )
        if len(checked_reads.split("\n"))==2:
            # Do nothing, data file is single end.
            pass 
        else:
            # Interleaved paired end reads. 
            # They have been deinterleaved with check_single_read_files.sh
            # Store new, deinterleaved readfiles in args.reads:
            args.format,args.reads[0]=checked_reads.split("\n")[:2]
            args.reads.append(checked_reads.split("\n")[2])

            
def check_paired_reads():
        ordered_reads=subprocess.check_output(
        "./check_paired_read_files.sh "+
        args.reads[0]+" "+
        args.reads[1],
        shell=True,
        encoding="utf8")
        if len(ordered_reads.split("\n"))==2:
            print("Error. Unrecognized header format for paired end reads.")
            sys.exit()
        else:
            args.format,args.reads[0],args.reads[1]=ordered_reads.split("\n")[:3]

def check_mapped():
    print("Checker for mapper has not been developed yet.")
    sys.exit()
    
def check_reference():
    print("Checker for reference has bot been developed yet.")
    #sys.exit()
    
def check_memory():
    print("Checker for memory has bot been developed yet.")
    sys.exit()
    
def check_cpu():
    print("Checker for CPU has not been developed yet.")
    sys.exit()
    


# ---Check CPU---
if args.cpu:
    check_cpu()
else:
    # Default CPU
    args.cpu=6
    
# ---Check memory---
if args.memory:
    check_memory()
else:
    # Default memory
    args.memory=4

# ---Check mapping file---
# If a mapping file (.bam/.sam) is provided, the reference used for mapping 
# must also be provided.
if args.mapped and not args.reference:
    print("Error. If a mapping file is provided, the reference used for mapping must also be provided.")
    sys.exit()

if args.mapped:
    check_mapped()
 
if args.reference:
    check_reference()
else:
    # Add Trinity assembly to config file.
    args.reference="../data/intermediate/trinity/Trinity.fasta"

# ---Check organism and/or annotation---
# Organism or annotation file must be provided.
# If none is provided, we cannot look for core genes.
if args.organism:
    check_organism()
elif args.annotation:
    check_annotation()
else:
    print("Error. Provide organism or annotation file.")
    sys.exit()


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
    check_single_reads()
# Check if two read files are provided.    
elif len(args.reads)==2:
    # Check that paired end reads are OK.
    check_paired_reads()    
# Too many readfiles were provided.    
else:
    print("Error. Too many read files.")
    sys.exit()


# add checker for fasta files

# add memory and CPU options

# Change config file

with open("config.json","r+") as configfile:
    data=json.load(configfile)
    data["input"]["reference"]=args.reference
    data["input"]["reads"]=args.reads
    data["input"]["annotation"]=args.annotation
    data["input"]["organism"]=args.organism
    data["input"]["mapped-reads"]=args.mapped
    configfile.seek(0)
    json.dump(data,configfile,indent=4)
    configfile.truncate()
    
# Add execution for snakefile...
#os.system("snakemake ../data/intermediate/trinity --dag | dot -Tsvg > dag.svg -forceall")
os.system("snakemake ../data/output/result_4.txt")