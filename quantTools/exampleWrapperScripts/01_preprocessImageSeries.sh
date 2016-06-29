#!/bin/bash
#
# Image series pre-processing before feeding into quantification pipeline
#   1. Assign relative symbolic links to images into the appropriate time series image directories
#   2. Rename image file prefix (e.g. "tile") if necessary
#
# Anthony Ho, ahho@stanford.edu, 6/28/2015


## Define input and output parameters
image_dir="$HOME/data/imaging/20160609_CNv3_3_ANUUN_01/05_sm_kd_1"
image_files_suffix="_green_2*_500ms_*.tif"

output_dir="images"
output_subdir_prefix="s"
output_file_prefix="tile"

num_timepoints=19
tiles=$(echo {1..15})


## Run
$HOME/scripts/arrayAnalysisTools/quantTools/preprocessImageSeries.sh $image_dir $image_files_suffix $output_dir $output_subdir_prefix $output_file_prefix $num_timepoints "$tiles" 
