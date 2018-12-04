import re
import sys
"""
Function for fixing the header of SRA-downloaded files
$ python fix_header.py in.fastq > out.fastq
"""

try:
    fastaopen = open(sys.argv[1], 'rU')
except:
    sys.stderr.write("The program expects the name of the fasta file. For example:\n   $ python fix_header.py in.fastq > out.fastq\n")

for line in fastaopen:
    #print(line)
    conten_new = re.sub('^(\S+)\.(\d+)\.(\d+)\s.+', r'\1.\2/\3', line)
    sys.stdout.write(conten_new)
