#!/usr/bin/env python3.6
# Top script for pipeline of GUESSmyLT.
# The script handles user arguments, validates inputs, modifies the read files
# (if needed) and executes the pipeline by calling Snakemake.
#
# Example runs:
#   With minimum inputs: python GUESSmyLT.py --reads read1 read2 --organism euk
#   With some optional inputs: python GUESSmyLT.py --reads read1 read2 --organism euk --threads 10 --memory 5G --reference reference.fasta
#
# For more options type: python GUESSmyLT.py --help


import argparse, sys, json, os, re, subprocess
from BCBio import GFF
import pprint
import gzip

# The path to tool directory. Helps with the relative paths.
script_dir = os.path.dirname(os.path.abspath(__file__)) + "/"
working_dir = os.getcwd() + "/"
# IF GUESSmyLT cannot find the correct path, change this manually. e.g.
# script_dir="/mnt/c/Users/erik_/Documents/GitHub/BioProject/GUESSmyLT"

# ---Function for checking user arguments---

def check_existence(file):
    """
    Checks that file exists and that it can be read.
    """
    try:
        with open(file,"r") as f:
            pass
    except IOError:
        print("Could not find file '"+file+"'. Make sure it exists.")
        sys.exit()


def check_readfile_extension(read):
    """
    Checks that a read file has either '.fq' or '.fastq' as file extension.
    """
    if not any(x in read for x in [".fq",".fastq"]):
        print("Error. Cannot find .fq or .fastq extension in read file --read "+read+".")
        sys.exit()


def get_first_line(read):
    """
    Function for getting the first line of a compressed or uncompressed file.
    """
    print(read)
    if read.endswith(".gz"):
        with gzip.open(read,"rt") as f:
            fline=f.readline()
    else:
        with open(read,"r") as f:
            fline=f.readline()
    return fline


def check_organism():
    """
    Checks that organism argument is correct.
    Can be any string within prokaryote or eukaryote.
    For example, euk, euka, eu are all valid arguments.
    """
    if args.organism.lower() in "prokaryote" and args.organism.lower() not in "eukaryote" :
        args.organism="prokaryote"

    elif args.organism.lower() in "eukaryote" and args.organism.lower() not in "prokaryote":
        args.organism="eukaryote"
    else:
        print("Error. Unrecognized organism --organism "+args.organism+". Only eukaryote/prokaryote are valid organism.")
        sys.exit()

def check_annotation(annotation_file):
    """
    Checks that an optional annotation file is correct by containing genes.
    """
    in_handle=open(annotation_file)
    for record in GFF.parse(in_handle, limit_info=dict(gff_type=["gene"])):
        if record:
            return True
    # If no genes are found, gff file cannot be used in analysis.
    print("Error. No genes could be found in --annotation "+annotation_file+". Please, submit a .gff file containing genes or no .gff file.")
    sys.exit()


def change_header_format(read,num):
    """
    Function that reformats the headers of a fastq file into usable format:
        Wrong header -> readID/num  , where num=1or2
    """
    num=str(num)
    tmp=output_dir+"data/intermediate/tmp.gz"
    fixed_header=read
    zip_command="cat"
    if read.endswith(".gz"):
        zip_command="gzip -cd"
    else:
        fixed_header=read+".gz"

    fline=get_first_line(read)
    if re.match("^@.+\.\d+\.\d ",fline):
        os.system(zip_command+" "+read+" | sed '/^@/s/\./:/g; /^@/s/:[0-9] .*/\/"+num+"/g' | gzip > "+tmp)

    elif re.match("^@.+ ",fline) and "." in fline:
        os.system(zip_command+" "+read+" | sed '/^@/s/\./:/g; /^@/s/ .*/\/"+num+"/g' | gzip > "+tmp)

    else:
        os.system(zip_command+" "+read+" | sed '/^@/s/ .*/\/"+num+"/g' | gzip > "+tmp)

    os.system("mv "+tmp+" "+fixed_header)
    return fixed_header

def check_subsample(num):
    """
    Checks that number of reads used for subsampling is even.
    Otherwise, it will not work if we work with paired end data.
    """
    if num%2==0:
        # Number is even --> Valid
        pass
    else:
        # Number is odd --> Not valid
        print("Error. Number of reads used for subsampling (--subsample) must be even.")
        sys.exit()


def subsample(read,num):
    """
    Selects the num first reads of a fastq files (subsamples) and copies them
    to a subsampled file.
    The subsampled file is used throughout the analysis instead of original file.
    This makes the run faster and also protects the original read file from being
    modified.
    """
    header_flag = True
    readname=re.split(".fq|.fastq",os.path.basename(read))[0]
    subsampled=output_dir+"data/intermediate/"+readname+".sub."+str(num)+".fastq.gz"
    if os.path.exists(subsampled):
        print("Already found file with "+str(num)+" subsampled reads.")
        print("Subsampled file "+subsampled+" will be used for analysis.")
        header_flag = False
    else:
        print("Subsampling "+str(num)+" reads of "+read+" to new read file: "+subsampled)
        # Multiply with 4 because a read consists of 4 lines in a fastq file.
        num=str(4*num)
        zip_command="cat"
        if read.endswith(".gz"):
            zip_command="zcat"

            # gzip gives broken pipe for compressed files, but works anyway.
        os.system(zip_command+" "+read+" | head -"+num+" | gzip > "+subsampled)
    return header_flag, subsampled


def check_single_reads(read,args):

    # --- Inner functions ---

    def is_interleaved(read,pattern1,pattern2):
        """
        Function for checking if a single read file is interleaved.
        If both pattern1 and pattern2 can be found in the first 10000 lines,
        True is returned. Else False.
        If two read headers are identical, True is returned. Else false.
        """
        if read.endswith(".gz"):
            f=gzip.open(read,"rt")
        else:
            f=open(read,"r")

        num_read1=0
        num_read2=0
        counter=0
        headers=[]

        for line in f:
            # Only look at read headers. They come every fourth line.
            if counter%4==0:
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
            counter+=1
            if counter==10000:
                break
        f.close()

        if num_read1>0 and num_read2>0:
            return True
        else:
            return False


    def deinterleave(read,line):
        """
        Function for deinterleaving reads into two separate files.
        It first checks for old and new Illumina fastq headers in
        order to split files. Otherwise, it assumes that the reads
        are alternating and therefore split all odd reads into one
        file and all even reads into another file.
        """
        print("Found interleaved reads in read file: "+read+".")

        readname=re.split(".fq|.fastq",os.path.basename(read))[0]
        deinterleaved=output_dir+"data/intermediate/deinterleaved."+readname+"."
        deinterleaved1=deinterleaved+"left.fastq.gz"
        deinterleaved2=deinterleaved+"right.fastq.gz"

        if os.path.exists(deinterleaved1) and os.path.exists(deinterleaved2):
            print("Already found existing deinterleaved files for read file: "+read)
            print("Continuing run with '"+deinterleaved1+"' and '"+deinterleaved2+"'.")
            return deinterleaved1,deinterleaved2


        print("Deinterleaving data.")
        zip_command="cat"
        if read.endswith(".gz"):
            zip_command="gzip -cd"

        # Check for Illumina(1.8+) or Illumina(1.8-) fastq format.
        if re.match("^@.+ 1:", line) and re.match("^@.+ 2:",line):
            print("Extracting left reads of '"+read+"' to '"+deinterleaved1+"'.")
            os.system(zip_command+" "+read+" | awk '/^@.*[ ]1:/{c=4} c&&c--' | gzip --stdout > "+deinterleaved1)
            print("Extracting right reads of '"+read+"' to '"+deinterleaved1+"'.")
            os.system(zip_command+" "+read+" | awk '/^@.*[ ]2:/{c=4} c&&c--' | gzip --stdout > "+deinterleaved2)
        elif re.match("^@.+/1",line) and re.match("^@.+/2",line):
            print("Extracting left reads of '"+read+"' to '"+deinterleaved1+"'.")
            os.system(zip_command+" "+read+" | awk '/^@.*[/]1/{c=4} c&&c--' | gzip --stdout > "+deinterleaved1)
            print("Extracting right reads of '"+read+"' to '"+deinterleaved2+"'.")
            os.system(zip_command+" "+read+" | awk '/^@.*[/]2/{c=4} c&&c--' | gzip --stdout > "+deinterleaved2)
        else:
            # Assuming that interleaved reads in a file are alternating.
            # This extraction can be written with a one line command:
            #   os.system(paste - - - - - - - - <"""+read+""" | tee >(cut -f 1-4 | tr "\t" "\n" > """+deinterleaved1+""") | cut -f 5-8 | tr "\t" "\n" > """+deinterleaved2)
            # The command works if you execute it from the command line, but not using os.system.
            # Right now, two lines are used, and this part can be optimized.

            print("Extracting left reads of '"+read+"' to '"+deinterleaved1+"'.")
            os.system(zip_command+' '+read+' | paste - - - - - - - - | cut -f 1-4 | tr "\t" "\n" | gzip > '+deinterleaved1)
            print("Extracting right reads of '"+read+"' to '"+deinterleaved2+"'.")
            os.system(zip_command+' '+read+' | paste - - - - - - - - | cut -f 5-8 | tr "\t" "\n" | gzip > '+deinterleaved2)


        return deinterleaved1,deinterleaved2

    # --- Outer function ---

    # Check that read files exists and can be read
    check_existence(read)

    # subsample reads from read file and continue analysis on new, subsampled file.
    read=subsample(read,args.subsample)[1]
    # Get first first line of read file
    fline=get_first_line(read)

    # Check if read file is interleaved.
    # If so we need to split reads (deinterleave) into two files.
    if is_interleaved(read,"^@.+/1","^@.+/2") or is_interleaved(read,"^@.+ 1:","^@.+ 2:"):
        # Deinterleaved reads
        read1,read2=deinterleave(read,fline)
        if any(symbol in get_first_line(read1) for symbol in [" ","_","."]) or any(symbol in get_first_line(read2) for symbol in [" ","_","."]):
            # Reformat headers so that no whitespaces underscores or punctutiations are present.
            # NEEDED since pysam and Trinity cannot handle that format.
            # Convertion is: 'Wrong header format' -> 'readID/pair#'        (pair# = 1 or 2 depending on the mate)
            print("Reformatting headers of: "+read1)
            read1=change_header_format(read1,1)
            print("Reformatting headers of: "+read2)
            read2=change_header_format(read2,2)
        return [read1,read2]

    else:
        # Single read file is not interleaved.
        if any(symbol in fline for symbol in [" ","_","."]):
            # Reformat headers so that no whitespaces underscores or punctutiations are present.
            # NEEDED since pysam and Trinity cannot handle that format.
            # Convertion is: 'Wrong header format' -> 'readID/pair#'        (pair# = 1 or 2 depending on the mate)
            print("Reformating read headers in fastq file: "+read)
            read=change_header_format(read,1)

        return [read]


def check_paired_reads(read1,read2,args):
    """
    Checks that read files exist and are in correct format.
    Updates format of read headers if whitespaces, underscores or punctuations are present.
    """

    check_existence(read1)
    check_existence(read2)

    check_readfile_extension(read1)
    check_readfile_extension(read2)

    # Subsample reads from read file and continue analysis on new, subsampled files.
    # Need to divide with 2 since subsample takes the total amount of reads that will be used.
    args.subsample=int(args.subsample/2)
    header_flag1, read1=subsample(read1,args.subsample)
    header_flag2, read2=subsample(read2,args.subsample)

    fline1=get_first_line(read1)
    fline2=get_first_line(read2)

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
        print("Unknown read header format.")
        print("First line in :"+read1+"\n\t"+fline1)
        print("First line in :"+read2+"\n\t"+fline2)
        print("Assuming left read file: "+ read1+ " and right read file: "+read2)


    if header_flag1 and header_flag2:
        if " " in fline1 or " " in fline2 or "." in fline1 or "." in fline2 or "_" in fline1 or "_" in fline2:
            # Reformat headers so that no whitespaces underscores or punctutiations are present.
            # NEEDED since pysam and Trinity cannot handle that format.
            # Convertion is: 'Wrong header format' -> 'readID/pair#'        (pair# = 1 or 2 depending on the mate)
            print("Reformatting headers of: "+read1)
            read1=change_header_format(read1,1)
            print("Reformatting headers of: "+read2)
            read2=change_header_format(read2,2)
    else:
        print("No header reformatting needed, continuing")

    return read1,read2


def check_mapped():
    # Not developed yet. Added in to do list.
    print("Checker for mapper has not been developed yet.")
    print("Right now you cannot provide a map-file (.bam).")
    print("Therefore, skip the map-file and let GUESSmyLT do the mapping for you.")
    sys.exit()


def check_reference(ref):
    """
    Checks that provided reference is valid according to:
        1. Can be opened.
        2. Has at least one line that begins with '>'.

    Also checks if we are dealing with genome or transcriptome by looking at
    the number of lines that begin with '>'.
    If one line begins with '>' we are dealing with genome.
    If more than one line begin with '>' we are dealing with transcriptome.
    """
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


def check_memory(args):
    """
    Checks that the maximum memory used are written correctly in gigabytes, e.g. 8G.
    """
    if "G" == args.memory[-1] and "0" != args.memory[0] and args.memory[:-1].isdigit():
        return True
    else:
        print("Invalid --memory argument: '"+args.memory+"'. Memory should be given in GIGABYTES. For example --memory '4G'")
        sys.exit()


# ---Main from here ---

def main():

    # Use argparse for handling user arguments.
    # Mandatory arguments are:
    # --organism and --reads
    # Optional arguments are:
    # --reference, --annotation, --mapped (For providing ref, ann, map files).
    # --threads, --memory (How many threads and maximum memory that GUESSmyLT can use.)

    parser=argparse.ArgumentParser(description="GUESSmyLT, GUESS my Library Type. Can predict the library type used for RNA-Seq. The prediction is based on the orientaion of your read files in .fastq format. Knowing the library type helps you with downstream analyses since it greatly improves the assembly.")
    parser._action_groups.pop()
    required=parser.add_argument_group("REQUIRED ARGUMENTS")
    required.add_argument("--organism",type=str, help="What organism are you dealing with? prokaryote or eukaryote.")
    optional=parser.add_argument_group("OPTIONAL ARGUMENTS")
    required.add_argument("--reads",nargs="+",type=str,help="One or two read files in .fastq format. Files can be compressed or uncrompressed. Handles interleaved read files and any known .fastq header format. ")
    optional.add_argument("--subsample",type=int,default=100000,help="Number of subsampled reads that will be used for analysis. Must be an even number. Default value is 100,000 reads.")
    optional.add_argument("--reference",type=str,help="Reference file in .fasta format. Reference can be either transcriptome or genome.")
    optional.add_argument("--annotation",type=str,help="Annotation file in .gff format. Needs to contain genes.")
    optional.add_argument("--mapped",type=str,help="Mapped file in sorted .bam format. Reference that reads have been mapped to has to be provided. Checker for this has not been developed yet. Therefore, do not provide a mapping file.") #Maybe add this?
    optional.add_argument("--threads", type=int, default=10,help="The number of threads that can be used by GUESSmyLT. Needs to be an integer. Defualt value is 10.")
    optional.add_argument("--memory",type=str, default="8G",help="Maximum memory that can be used by GUESSmyLT in GB. E.g. '10G'. Default value is 8G.")
    optional.add_argument("--output",type=str,default=working_dir,help="Full path to output directory")
    args=parser.parse_args()

    global output_dir
    output_dir = args.output
    if not os.path.exists(args.output+"data/intermediate"):
        os.system("mkdir "+args.output+"data")
        os.system("mkdir "+args.output+"data/intermediate")

    if not os.path.exists(args.output+"config.json"):
        os.system("cp "+script_dir+"config.json "+args.output+"config.json")
    # ---Check subsampling---
    check_subsample(args.subsample)

    # ---Check reads---
    # At least one readfile must be provided. (single end or interleaved paired end reads)
    # At most two readiles can be provided. (paired end reads)
    #
    # Check if readfile is provided.
    if not args.reads:
        print("Error. No read files provided. At least one read file should be provided.")
        sys.exit()
    # Check if one read file is provided.
    elif len(args.reads)==1:
        # Check that single reads are OK
        # or deinterleave paired end reads.
        args.reads=check_single_reads(args.reads[0],args)
    # Check if two read files are provided.
    elif len(args.reads)==2:
        # Check that paired end reads are OK.
        args.reads[0],args.reads[1]=check_paired_reads(args.reads[0],args.reads[1],args)

    # Too many readfiles were provided.
    else:
        print("Error. Too many read files. Only one or two read files can be provided.")
        sys.exit()

    # After reads have been checked, they are sorted so that read[1/left/forward] is first and read[2/right/reverse] is second.
    # Base readname on first read:
    # Get readname, i.e: readname from path/to/reads/[readname].fq.gz
    # Readname is used to make each run unique.
    readname=re.split(".fq|.fastq",os.path.basename(args.reads[0]))[0]

    # ---Check mapping file---
    # If a mapping file (.bam/.sam) is provided, the reference used for mapping
    # must also be provided.
    if args.mapped and not args.reference:
        print("Error. If a mapping file is provided, the reference used for mapping must also be provided.")
        sys.exit()

    # Tells busco what type of reference: genome or transcriptome.
    # Right now we always tell busco that it's working with a genome,
    # because it gives less result files if it is a transcriptome.
    # This could be optimized.
    busco_reference_mode="genome"
    if args.reference:
        busco_reference_mode=check_reference(args.reference)
    else:
        # Add Trinity assembly to config file.
        args.reference="data/intermediate/"+readname+".fasta"


    if args.mapped:
        check_mapped()
    else:
        args.mapped="data/intermediate/"+readname+"_on_ref_"+re.split('/|\.',args.reference)[-2]+"_sorted.bam"


    if args.annotation:
        check_annotation(args.annotation)
    else:
        args.annotation="data/intermediate/run_"+re.split('/|\.',args.reference)[-2]


    check_memory(args)


    # Update config file
    #script_dir="."
    config_path = args.output+"config.json"
    with open(config_path,"r+") as configfile:
        data=json.load(configfile)
        data["trinity"]["reference"]=args.reference
        data["input"]["reads"]=args.reads
        data["input"]["readname"]=readname
        data["busco"]["annotation"]=args.annotation
        data["input"]["organism"]=args.organism
        data["bowtie2"]["mapped-reads"]=args.mapped
        data["busco"]["lineage"]=data["prokaryote_db"] if args.organism == "prokaryote" else data["eukaryote_db"]
        data["busco"]["mode"]=busco_reference_mode
        data["input"]["threads"]=args.threads
        data["input"]["memory"]=args.memory
        data["output"]=args.output
        data["script_dir"]=script_dir
        configfile.seek(0)
        json.dump(data,configfile,indent=4)
        configfile.truncate()

    # Execute Snakemake
    os.system("snakemake -s "+script_dir+"Snakefile -d "+args.output+" data/output/result_"+readname+".txt --cores "+str(args.threads))

if __name__ == "__main__":
    main()
