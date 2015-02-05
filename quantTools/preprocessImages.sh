#!/bin/bash
#
# Image pre-processing before quantification
#   1. Rotate images by 180 degrees
#   2. Rename images to start with "tile" if not already
#   3. Assign images into appropriate time series image directories
#
# Usage: 
#   imagePreprocessing.sh image_dir image_files_suffix output_dir output_subdir_prefix output_file_prefix num_timepoints "tiles"
#
# Anthony Ho, ahho@stanford.edu, 1/14/2015
# Last update 1/14/2015


# Define usage function
usage() {
    cat <<EOF

$0 image_dir image_files_suffix output_dir output_subdir_prefix output_file_prefix num_timepoints "tiles" 

This script preprocesses images before quantification by: 
  1. Rotate images by 180 degrees
  2. Rename images to start with "tile" if not already
  3. Assign images into appropriate time series image directories

EOF
}


# Check number of arguments
num_arguments=7

if [ $# -gt $num_arguments ]
then
    echo "Too many arguments!"
    usage
    exit 1
elif [ $# -lt $num_arguments ]
then
    echo "Not enough arguments!"
    usage
    exit 1
fi


# Parse arguments
image_dir=$1
image_files_suffix=$2
output_dir=$3
output_subdir_prefix=$4
output_file_prefix=$5
num_timepoints=$6
tiles=$7


# Make empty directories for each time point if not exist already
for (( timepoint=1; timepoint<=$num_timepoints; timepoint++ ))
do
    mkdir -p ${output_dir}/${output_subdir_prefix}${timepoint}
done


# Declare a function to convert for each tile
convert_each_tile() {

    # Parse argument
    image_dir=$1
    image_files_suffix=$2
    output_dir=$3
    output_subdir_prefix=$4
    output_file_prefix=$5
    tile=$6

    # Initialize timepoint variable
    timepoint=0
    
    # Go through all timepoints of the same tile
    shopt -s nullglob
    for F in ${image_dir}/{,${output_file_prefix}}${tile}${image_files_suffix}
    do 
	((timepoint++))
	Fbase=$(basename $F)
	# Rename images to start with "tile" if not already
	if [[ $Fbase != ${output_file_prefix}* ]] 
	then
	    Fbase=${output_file_prefix}${Fbase}
	fi
	convert $F -rotate 180 ${output_dir}/${output_subdir_prefix}${timepoint}/${Fbase}
    done
}
export -f convert_each_tile


# Convert in parallel for each tile
parallel \
    convert_each_tile $image_dir $image_files_suffix $output_dir $output_subdir_prefix $output_file_prefix {} \
    ::: $tiles
