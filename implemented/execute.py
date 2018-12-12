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

args=parser.parse_args()

read_basename=os.path.basename(args.reads[0])
read_basename_no_extension=re.split(".fq|.fastq",read_basename)[0]
if len(re.split(".fq|.fastq",read_basename))!=2:
    print("Error. Reads are not in (zipped) fastq format.")
    sys.exit()


# Check argument criteria

# 1. If a mapped file is provided, we also need the reference file the reads were mapped to.
if args.mapped and not args.reference:
    print("Error. If map file is provided, the reference it was mapped to must also be provided.")
    sys.exit()

# 2.a An organism or annotation file must be provided.
if not args.organism and not args.annotation:
    print("Error. Provide organism or annotation file.")
    sys.exit()

# 2.b Must recognize organism
elif args.organism.lower() in "prokaryote":
    args.organism="prokaryote"

elif args. organism.lower() in "eukaryote":
    args.organism="eukaryote"
else: 
    print("Error. Unrecognized organism.")
    sys.exit()

# 3. Check that one or two reads are provided.
if not args.reads:
    print("Error. No read files provided.")
    sys.exit()
elif len(args.reads)>2:
    print("Error. Too many read files.")
    sys.exit()
    # 4. If two read files are provided, we need to know which is forward and reverse by checking names.
    # Other solution might be to check headers in read files...
elif len(args.reads)==2:
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
    # 5. If only one read file was provided. Check if it paired end or single.
elif len(args.reads)==1:
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
    


# add checker for fastq files

# add checker for fasta files

# add memory and CPU options

# Change config file

with open("config2.json","r+") as configfile:
    data=json.load(configfile)
    data["input"]["reference"]=args.reference
    data["input"]["reads"]=args.reads
    data["input"]["annotation"]=args.annotation
    data["input"]["organism"]=args.organism
    configfile.seek(0)
    json.dump(data,configfile,indent=4)
    configfile.truncate()
    
# Add execution for snakefile...
# os.system("snakemake + target files")
