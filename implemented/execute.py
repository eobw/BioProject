# Top script for pipeline.
# Execute with: python execute.py --reads read1 (read2) --organism pro/euk
# For more options type: python execute.py --help





# ADD CHECKER FOR IF FILES EXIST!!!!!!!!!!!!!!!!!!!!!!!

import argparse, sys, json, os, re, subprocess
from BCBio import GFF
import pprint
import gzip

parser=argparse.ArgumentParser(description="Predict library type.")
parser._action_groups.pop()
required=parser.add_argument_group("required arguments")
required.add_argument("--organism",type=str)
optional=parser.add_argument_group("optional arguments")
#required.add_argument("--reads",nargs="+",type=lambda x: is_valid_file(parser,x),help="Reads in (un)zipped .fastq format.")
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

def open_zipped():
    pass
# Maybe add exception class
def check_existence(file):
    try:
        with open(file,"r") as f:
            pass
    except IOError:
        print("Could not find file '"+file+"'. Make sure it exists.")
        sys.exit()

def check_readfile_extension(read):
    if not any(x in read for x in [".fq",".fastq"]):
        print("Error. Cannot find .fq or .fastq extension in read file --read "+read+".")
        sys.exit()

def get_first_line(read):
    if read.endswith(".gz"):
        with gzip.open(read,"rt") as f:
            fline=f.readline()
    else:
        with open(read,"r") as f:
            fline=f.readline()
    return fline
            

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

def whitespace2underscore(read):
    if read.endswith(".gz"):
        os.system("""
        f="""+read+"""
        cp "$f" "$f~" &&
        gzip -cd "$f~" | sed '/^@/s/ /_/g' | gzip > "$f"
        rm "$f~"
        """)
    else:
        os.system("sed -i -e '/^@/s/ /_/g' "+read)

def check_single_reads(read):
            
    def is_interleaved(read,pattern1,pattern2):
        
        if read.endswith(".gz"):
            f=gzip.open(read,"rt")
        else:
            f=open(read,"r")    
        
        num_read1=0
        num_read2=0
        counter=0
        headers=[]
        
        for line in f:
            counter+=1
            if re.match(pattern1,line):
                num_read1+=1
            elif re.match(pattern2,line):
                num_read2+=1
    
            if line.startswith("@"):
                ID=line.split(" ")[0]
                if ID in headers:
                    return True
                else:
                    headers.append(ID) 
                    
            if counter==10000:
                break        
        f.close()
        
        if num_read1>0 and num_read2>0:
            return True
        else:
            return False
    
    
    def deinterleave(read,line):
        print("Found interleaved reads in read file: "+read+".")
        print("Deinterleaving data.")
        readname=re.split(".fq|.fastq",os.path.basename(read))[0]
        deinterleaved="../data/input/deinterleaved."+readname
        deinterleaved1=deinterleaved+".r1.fastq.gz"
        deinterleaved2=deinterleaved+".r2.fastq.gz"
        
        zip_command="cat"
        if read.endswith(".gz"):
            zip_command="zcat"
            
        if re.match("^@.+ 1:", line) and re.match("^@.+ 2:",line):
            print("Extracting left reads of '"+read+"' to '"+deinterleaved1+"'.")
            os.system(zip_command+" "+read+" | awk '/^@.*[ ]1:/{c=4} c&&c--' | gzip --stdout > "+deinterleaved+".r1.fastq.gz")
            print("Extracting right reads of '"+read+"' to '"+deinterleaved1+"'.")
            os.system(zip_command+" "+read+" | awk '/^@.*[ ]2:/{c=4} c&&c--' | gzip --stdout > "+deinterleaved+".r2.fastq.gz")
        elif re.match("^@.+/1",line) and re.match("^@.+/2",line):
            print("Extracting left reads of '"+read+"' to '"+deinterleaved1+"'.")
            os.system(zip_command+" "+read+" | awk '/^@.*[/]1/{c=4} c&&c--' | gzip --stdout > "+deinterleaved+".r1.fastq.gz")
            print("Extracting right reads of '"+read+"' to '"+deinterleaved2+"'.")
            os.system(zip_command+" "+read+" | awk '/^@.*[/]2/{c=4} c&&c--' | gzip --stdout > "+deinterleaved+".r2.fastq.gz")
        else:
            # Assuming that interleaved reads in a file are alternating.
            # This extraction can be written with a one line command:
            #   os.system(paste - - - - - - - - <"""+read+""" | tee >(cut -f 1-4 | tr "\t" "\n" > """+deinterleaved1+""") | cut -f 5-8 | tr "\t" "\n" > """+deinterleaved2)
            # The command works if you execute it from the command line, but not using os.system.
            # Right now, two lines are used, and this part can be optimized.
            print("Extracting left reads of '"+read+"' to '"+deinterleaved1+"'.")
            #os.system("""paste - - - - - - - - < """+read+""" | cut -f 1-4 | tr "\t" "\n" | gzip > """+deinterleaved1)
            os.system(zip_command+""" """+read+""" | paste - - - - - - - - | cut -f 1-4 | tr "\t" "\n" | gzip > """+deinterleaved1)
            print("Extracting right reads of '"+read+"' to '"+deinterleaved2+"'.")
            os.system(zip_command+""" """+read+""" | paste - - - - - - - - | cut -f 5-8 | tr "\t" "\n" | gzip > """+deinterleaved2)
            #os.system("""paste - - - - - - - - < """+read+""" | cut -f 5-8 | tr "\t" "\n" | gzip > """+deinterleaved2)
        
        if " " in get_first_line(read):
            print("Detected whitespaces in read file headers.")
            print("Substituting whitespaces in deinterleaved read file: "+deinterleaved1)
            whitespace2underscore(deinterleaved1)
            print("Substituting whitespaces in deinterleaved read file: "+deinterleaved1)
            whitespace2underscore(deinterleaved2)    
        
        return deinterleaved1,deinterleaved2
    
    
    check_existence(read)
        
    fline=get_first_line(read)
    
    if is_interleaved(read,"^@.+/1","^@.+/2") or is_interleaved(read,"^@.+ 1:","^@.+ 2:"):
        read1,read2=deinterleave(read,fline)
        if re.search(" ",fline) and args.reference:
            whitespace2underscore(read1)
            whitespace2underscore(read2)
            
        return read1,read2
    
    else:
        if re.search(" ",fline) and not args.reference:
            # Trinity cannot handle whitespaces in fastq-headers.
            # Therefore, if Trinity is to be executed (when no reference is given),
            # we need to replace the whitespaces with e.g. "_"
            print("WARNING. Noticed whitespaces in headers of file: "+read+".")
            print("Trinity cannot handle this format.")
            print("Would you like to change the header format of your read files?")
            print("Whitespaces (' ') in every header will be substituded with underscores ('_').")
            while True:
                decision=input("Proceed? (y/n): ")
                if decision.lower() in "yes" or decision.lower() in "no":
                    break
                else:
                    continue
            if decision.lower() in "yes":
                print("Substituting whitespaces with underscores in headers of read file: "+read)
                whitespace2underscore(read)

            else: # decision = 'no'
                print("No substitution will be done. Therefore Trinity cannot be executed and GUESSmyLT cannot continue.")
                print("To continue with the analysis you have 3 options:")
                print(" 1. Provide an assmebly/reference so that Trinity can be skipped.")
                print(" 2. Manually substitute the whitespaces in the headers of "+read+".")
                print(" 3. Re-execute GUESSmyLT and choose ('y') to substitute the whitespaces with underscores.")            
                print("Exiting GUESSmyLT.")
                sys.exit()    
            
        return read

            
def check_paired_reads(read1,read2):
    check_existence(read1)
    check_existence(read2)
    
    check_readfile_extension(read1)
    check_readfile_extension(read2)
 
    
    if read1.endswith(".gz"):
        with gzip.open(read1,"rt") as f:
            fline1=f.readline()
    else:
        with open(read1,"r") as f:
            fline1=f.readline()
            
    if read2.endswith(".gz"):
        with gzip.open(read2,"rt") as f:
            fline2=f.readline()
    else:
        with open(read2,"r") as f:
            fline2=f.readline()
            
    # Check format.
    # Either Illumina1.8+: @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<xpos>:<y-pos> <read>:<is filtered>:<control number>:<index>
    # Or Illumina1.8-: @<machine_id>:<lane>:<tile>:<x_coord>:<y_coord>#<index>/<read>
    format="Unknown"
    
    if re.match("^@.+/1",fline1) and re.match("^@.+/2",fline2):
        format="Old Illumina (Illumina1.8-)"
        pass
    elif re.match("^@.+/2",fline1) and re.match("^@.+/1",fline2):
        read1,read2=read2,read1
        format="Old Illumina (Illumina1.8-)"
    
    elif re.match("^@.*\ 1:",fline1) and re.match("^@.*\ 2:",fline2):
        format="New Illumina (Illumina1.8+)"
        pass
    elif re.match("^@.*\ 2:",fline1) and re.match("^@.*\ 1:",fline2):
        format="New Illumina (Illumina1.8+)"
        read1,read2=read2,read1
        
    print("Based on read headers we are working with format: "+format)
    
    if format=="Unknown":
        print("WARNING. Unknown read header format.")
        print("Assuming left read file: "+ read1+ " and right read file: "+read2)
        
        
    
    if re.search(" ",fline1) and re.search(" ",fline2) and not args.reference:
        # Trinity cannot handle whitespaces in fastq-headers.
        # Therefore, if Trinity is to be executed (when no reference is given),
        # we need to replace the whitespaces with e.g. "_"
        print("WARNING. Noticed whitespaces in headers of files: "+read1+", "+read2+".")
        print("Trinity cannot handle this format.")
        print("Would you like to change the header format of your read files?")
        print("Whitespaces (' ') in every header will be substituded with underscores ('_').")
        while True:
            decision=input("Proceed? (y/n): ")
            if decision.lower() in "yes" or decision.lower() in "no":
                break
            else:
                continue

        if decision.lower() in "yes":
            print("Substituting whitespaces with underscores in headers of read file: "+read1)
            whitespace2underscore(read1)
            print("Substituting whitespaces with underscores in headers of read file: "+read2)    
            whitespace2underscore(read2)
        
        else: # decision = 'no'
            print("No substitution will be done. Therefore Trinity cannot be executed and GUESSmyLT cannot continue.")
            print("To continue with the analysis you have 3 options:")
            print(" 1. Provide an assmebly/reference so that Trinity can be skipped.")
            print(" 2. Manually substitute the whitespaces in the headers of "+read1+" and "+read2+".")
            print(" 3. Re-execute GUESSmyLT and choose ('y') to substitute the whitespaces with underscores.")            
            print("Exiting GUESSmyLT.")
            sys.exit()    
 
 
    return read1,read2

        
def check_mapped():
    print("Checker for mapper has not been developed yet.")
    print("Right now you cannot provide a map-file (.bam).")
    print("Therefore, skip the map-file and let GUESSmyLT do the mapping for you.")
    sys.exit()
    
def check_reference(ref):
    try:
        open(ref)
    except:
        print("Error. Cannot open --reference '"+ref+"'. Make sure it exists.")
        sys.exit()
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

 
if args.reference:
    busco_reference_mode=check_reference(args.reference)
else:
    # Add Trinity assembly to config file.
    args.reference="../data/intermediate/"+readname+".fasta"


if args.mapped:
    check_mapped()
else:
    #args.mapped="../data/intermediate/"+readname+"_sorted.bam"
    args.mapped="../data/intermediate/"+readname+"_on_ref_"+re.split('/|\.',args.reference)[-2]+"_sorted.bam"

busco_reference_mode="genome"
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

if args.annotation:
    check_annotation(args.annotation)
else:
    args.annotation="../data/intermediate/run_"+re.split('/|\.',args.reference)[-2]


check_memory()


# add checker for fasta files

# add memory and threads options

# Change config file

with open("config.json","r+") as configfile:
    data=json.load(configfile)
    data["trinity"]["reference"]=args.reference
    data["input"]["reads"]=args.reads
    data["input"]["readname"]=readname
    data["busco"]["annotation"]=args.annotation
    data["input"]["organism"]=args.organism
    data["bowtie2"]["mapped-reads"]=args.mapped
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
#os.system("snakemake -s bowtie2 "+args.mapped+" --forceall")
#os.system("snakemake -s trinity --cores "+str(args.threads))
#os.system("snakemake ../data/output/result_"+readname+".txt --cores "+str(args.threads))