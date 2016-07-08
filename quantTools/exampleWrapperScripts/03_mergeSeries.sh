#!/bin/bash
#
# Merge time series data of all tile into single files
#
# Anthony Ho, ahho@stanford.edu, 7/8/2016


## Defining input and output parameters
base_dir=$HOME/"analysis/imaging/20160609_CNv3_3_ANUUN_01"
exp_dir=$1
tiles=$2
num_timepoints=$3

CPfluor_dir=${base_dir}/${exp_dir}/"images/t##/CPfluor"
CPfluor_file_prefix="ANUUN_ALL_tile"
CPfluor_file_suffix="_Bottom_filtered*CPfluor"

output_dir=${base_dir}/${exp_dir}/"mergedSeries"
output_file_prefix="ANUUN_ALL_tile"
output_file_suffix="_Bottom_filtered"


## Make empty output directory if not exist already
mkdir -p $output_dir


## Merge in parallel for each tile
parallel \
    mergeSeries.py $num_timepoints \
    ${CPfluor_dir}/${CPfluor_file_prefix}{}${CPfluor_file_suffix} \
    ${output_dir}/${output_file_prefix}{}${output_file_suffix} \
    ::: $tiles
