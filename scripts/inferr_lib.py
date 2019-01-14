"""
Script for inferring library type of mapped reads
$ inferr_lib.py [run_name] [single OR paired]
"""

import sys
import re
from difflib import SequenceMatcher
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pysam

# ---------- FUNCTIONS ----------

def extract_genes(run_name):
    '''
    Function for extracting genes corresponding to BUSCO hits (genes).
    Returns a SeqRecord object with one feature per BUSCO hit.
    '''

    file_tsv = open("../data/intermediate/run_"+run_name+"/full_table_"+run_name+".tsv", 'r')

    # Extract BUSCO IDs, start and end from table of hits into SeqRecord object, each BUSCO as a SeqFeature
    busco_record = SeqRecord(seq='', id='hits')
    for line in file_tsv.readlines():
        hit = (re.search(r'(\S*)\s(Complete)\s(\S*)\s(\S*)\s(\S*)\s\S*', line))
        if hit:
            busco_record.features.append(SeqFeature(FeatureLocation(int(hit.group(4)), int(hit.group(5))), id=hit.group(1), type='gene', qualifiers={'contig': hit.group(3)}))

    file_tsv.close()

    # Match the BUSCOs to augustus predicted genes in gff file
    gff_records = []
    correct_genes = SeqRecord(seq='', id='correct_genes')
    limit_infos = dict(
            gff_type = ["gene"]) # Only want genes

    for busco in busco_record.features:
        filename = busco.id # gff filenames are [busco_id].out.[1-999]
        i = 1
        while True:
            try:
                file_gff = open("../data/intermediate/run_"+run_name+"/augustus_output/predicted_genes/"+filename+".out."+str(i))
                for record in GFF.parse(file_gff, limit_info=limit_infos):
                    gff_records.append(record)
                i += 1
            except: # Finished reading all files for that BUSCO
                break

    # Find augustus predicted genes from the gff that match BUSCOs
    for rec in gff_records:
        for hit in busco_record.features:
            for feature in rec.features:
                if hit.location.start-1 == feature.location.start and hit.location.end == feature.location.end: # For some reason start has 1 nt diff...
                    feature.id = rec.id
                    correct_genes.features.append(feature)
                    break

    file_gff.close()
    print("Number of genes extracted: %d\n" % (len(correct_genes.features)))
    return correct_genes

def infer_paired_region(genes):
    '''
    Function for inferring paired library-type by looking at a regions corresponding to genes
    '''
        # Counters for the different lib-types
    libs = {
        'fr_first': 0,
        'fr_second': 0,
        'rf_first': 0,
        'rf_second': 0,
        'ff_first': 0,
        'ff_second': 0,
        'undecided': 0
    }

    for gene in genes.features:
        contig = gene.id
        start = int(gene.location.start)
        stop = int(gene.location.end)
        strand = gene.strand
        reads = []
        # Get reads mapped to a specific contig and in a sequence range
        # TODO: Look into optimizing this step, only take a subset (1000ish) reads?
        # samfile.mate is not made for high throughput
        for read in samfile.fetch(contig, start, stop):
            if not read.mate_is_unmapped and read.is_read1:
                reads.append([read, samfile.mate(read)])



        # Check lib-type of reads
        for read in reads:
            first = read[0]
            second = read[1]
            try:
                lib = ''
                if not first.is_reverse:
                    lib += 'f'
                else:
                    lib += 'r'
                if not second.is_reverse:
                    lib += 'f'
                else:
                    lib += 'r'
                # Gene on sense strand
                if strand == 1:
                    if first.reference_start > second.reference_start:
                        # Flip order of reads
                        lib = lib[::-1]
                        lib += '_first'
                    elif first.reference_start < second.reference_start:
                        lib += '_second'
                    else:
                        lib = 'undecided'
                # Gene on antisense
                elif strand == -1:
                    if first.reference_start > second.reference_start:
                        # Flip order of reads
                        lib = lib[::-1]
                        lib += '_second'
                    elif first.reference_start < second.reference_start:
                        lib += '_first'
                    else:
                        lib = 'undecided'
                libs[lib] += 1
            except: libs['undecided'] +=1 #Some reads missing start or end-values

    return libs


def infer_single_region(genes):
    """
    Function for inferring library type of single-ended library types
    """

    libs = {
        'f_first': 0,
        'f_second': 0,
        'r_first': 0,
        'r_second': 0,
        'undecided': 0
    }

    og_reads = {}
    for record in SeqIO.parse(run_path, "fastq"):
        og_reads[record.id]= str(record.seq)

    for gene in genes.features:
        contig = gene.id
        start = int(gene.location.start)
        stop = int(gene.location.end)
        strand = gene.strand
        reads = []
        # Get reads mapped to a specific contig and in a sequence range
        # Slow, only take 1000 reads?
        for read in samfile.fetch(contig, start, stop):
            if not read.is_unmapped:
                reads.append(read)

        # Check lib-type of reads
        for read in reads:
            try:
                flag = SequenceMatcher(None, og_reads[read.query_name], read.query_sequence).ratio() >= 0.8
                lib = ''
                if strand == 1:
                    if flag and not read.is_reverse:
                        lib += 'f_second'
                    elif not flag and read.is_reverse:
                        lib += 'f_first'
                    elif flag and read.is_reverse:
                        lib += 'r_second'
                    elif not flag and not read.is_reverse:
                        lib += 'r_first'
                    else:
                        lib = 'undecided'
                elif strand == -1:
                    if not flag and read.is_reverse:
                        lib += 'f_second'
                    elif flag and  not read.is_reverse:
                        lib += 'f_first'
                    elif flag and read.is_reverse:
                        lib += 'r_first'
                    elif not flag and  not read.is_reverse:
                        lib += 'r_second'
                    else:
                        lib = 'undecided'
                else:
                    print(read)
                    lib = 'undecided'
                libs[lib] += 1
            except: libs['undecided'] +=1 # Some reads missing start or end-values
    return libs

def write_result(lib_dict):
    """
    Function for calculating precentages and writing results of inferring to resultfile.
    """
    output = open('../data/output/result_'+run_name+'.txt', 'w+')
    output.write("Results of library inferring "+run_name+": \nLibrary type \t Reads \t Percent \n")
    print("Results of library inferring: \nLibrary type \t Reads \t Percent \n")

    total_reads = 0
    for i in lib_dict:
        total_reads += lib_dict[i]
    if total_reads > 0:
        for i in lib_dict:
            percent = '{:.1%}'.format(lib_dict[i]/total_reads)
            output.write("%s \t %d \t %s\n" % (i, lib_dict[i], percent))
            print("%s \t %d \t %s\n" % (i, lib_dict[i], percent))
    else:
        for i in lib_dict:
            output.write("%s \t %d \t %s\n" % (i, 0, 0))
            print("%s \t %d \t %s\n" % (i, 0, 0))

# ---------- RUNNING ----------

reference = sys.argv[1]
mapped_reads = sys.argv[2]
run_name = sys.argv[3]
run_path = sys.argv[4]
state = sys.argv[5]




samfile = pysam.AlignmentFile(mapped_reads, "rb")

print("Extracting genes...")
genes = extract_genes(reference)

if state == 'single':
    print("Running single end inferring")
    result = infer_single_region(genes)
elif state == 'paired':
    print("Running paired end inferring")
    result = infer_paired_region(genes)


print("Prediction finished:\n")
write_result(result)
print("Results also written to file data/output/result_"+run_name+".txt")
