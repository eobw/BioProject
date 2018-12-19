#!/bin/bash

# Readfiles
read1=$1
read2=$2


# Check if file is compressed and get first line of files.
#
# Get first line of readfile 1 
if [[ $read1 =~ \.gz$ ]]; then
  fline1=$(zcat $read1 | head -1)
else 
  fline1=$(cat $read1 | head -1)
fi
#
# Get first line of readfile 2
if [[ $read2 =~ \.gz$ ]]; then
  fline2=$(zcat $read2 | head -1)
else 
  fline2=$(cat $read2 | head -1)
fi


# Check format.
# Either Illumina1.8+: @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<xpos>:<y-pos> <read>:<is filtered>:<control number>:<index>
# Or Illumina1.8-: @<machine_id>:<lane>:<tile>:<x_coord>:<y_coord>#<index>/<read>

# Illumina1.8-
if [[ $fline1 =~ ^@.*/1 ]] && [[ $fline2 =~ ^@.*/2 ]]; then
  echo "Illumina1.8-"
  echo $read1
  echo $read2  
elif [[ $fline1 =~ ^@.*/2 ]] && [[ $fline2 =~ ^@.*/1 ]]; then
  echo "Illumina1.8-"
  echo $read2 
  echo $read1

# Illumina1.8+
elif [[ $fline1 =~ ^@.*\ 1: ]] && [[ $fline2 =~ ^@.*\ 2: ]]; then
  echo "Illumina1.8+"
  echo $read1
  echo $read2 
elif [[ $fline1 =~ ^@.*\ 2: ]] && [[ $fline2 =~ ^@.*\ 1: ]]; then
  echo "Illumina1.8+"
  echo $read2
  echo $read1
else
  # Unknown format... Cannot work with this data...
  echo 
fi