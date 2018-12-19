#!/bin/bash

# Readfiles
read=$1

inputdatapath="../data/input/"

# maybe add so that output files changes according to read files name.
# out_basename=$...

# Check if file is compressed and adjust command after that.
# 
if [[ $read =~ \.gz$ ]]; then
  fline=$(zcat $read | head -1)
  zip_command=zcat
else 
  zip_command=cat
  fline=$(cat $read | head -1)
fi



# Check format.
# Either Illumina1.8+: @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<xpos>:<y-pos> <read>:<is filtered>:<control number>:<index>
# Or Illumina1.8-: @<machine_id>:<lane>:<tile>:<x_coord>:<y_coord>#<index>/<read>

# Illumina1.8-
if [[ $fline =~ ^@.*/1 ]] || [[ $fline =~ ^@.*/2 ]]; then
  #hej=$(awk '/^@.*[/]1/{c=4} c&&c--' $fline)
  num_read1=$($zip_command $read | head -10000 | grep '^@.*/1' | wc -l)
  num_read2=$($zip_command $read | head -10000 | grep '^@.*/2' | wc -l)
  
  if [[ $num_read1 > 0 ]] && [[ $num_read2 > 0 ]]; then
    $zip_command $read | awk '/^@.*[/]1/{c=4} c&&c--' | gzip --stdout > ${inputdatapath}deinterleaved.left.fastq.gz
    $zip_command $read | awk '/^@.*[/]2/{c=4} c&&c--' | gzip --stdout > ${inputdatapath}deinterleaved.right.fastq.gz
    echo "Illumina1.8-"
    echo "deinterleaved.left.fastq.gz"
    echo "deinterleaved.right.fastq.gz"
  fi

# Illumina1.8+
elif [[ $fline =~ ^@.*\ 1: ]] || [[ $fline =~ ^@.*\ 2: ]]; then
  num_read1=$($zip_command $read | head -10000 | grep '^@.* 1:' | wc -l)
  num_read2=$($zip_command $read | head -10000 | grep '^@.* 2:' | wc -l)

  if [[ $num_read1 > 0 ]] && [[ $num_read2 > 0 ]]; then
    $zip_command $read | awk '/^@.*[ ]1:/{c=4} c&&c--' | gzip --stdout > ${inputdatapath}deinterleaved.left.fastq.gz
    $zip_command $read | awk '/^@.*[ ]2:/{c=4} c&&c--' | gzip --stdout > ${inputdatapath}deinterleaved.right.fastq.gz
    echo "Illumina1.8+"
    echo "deinterleaved.left.fastq.gz"
    echo "deinterleaved.right.fastq.gz"
  fi

else
  # Single end data. Don't do anything.
  echo 
fi