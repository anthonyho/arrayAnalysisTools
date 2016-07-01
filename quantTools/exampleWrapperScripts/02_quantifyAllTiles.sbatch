#!/bin/bash
#
# Bash script to call the quantifyAllTiles.py script using SLURM scheduler on Sherlock
# 
# Usage:
#   02_quantifyAllTiles.sh exp_dir "timepoints"
# 
# Anthony Ho, ahho@stanford.edu, 6/28/2016


## Define SLURM options
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH --time=40:00:00
#################
#quality of service; think of it as job priority
#SBATCH --partition=owners
#SBATCH --qos=normal
##################
#number of nodes you are requesting
#SBATCH --nodes=1
##################
#task to run per node;
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15
#################


## Load matlab
module load matlab/R2012b

## Define paths, input files and options
base_dir=$HOME/"analysis/imaging/20160609_CNv3_3_ANUUN_01"
exp_dir=$1
read_dir=$HOME/"data/sequencing/20160607_CNv3_3_L16v2a3t2WTt2cs_9_ANUUN"

# Define timepoints to be quantified with curly brace expansion
timepoints=$2

image_dir=$base_dir/$exp_dir/"images"
image_subdir_prefix="t"

filter=$base_dir/"CNv3_L16v2a3t2c.CPfilter"
tile_dir=$base_dir/"tiles"
filtered_tile_dir=$base_dir/"filtered_tiles"

gv_dir=$base_dir

# Enable if using the same registration offset map 
#roff_dir=$base_dir/$exp_dir/$image_dir/"roff"
roff_dir="roff"

num_cores="15"
reg_subset1="FMv1"
reg_subset2="CNv3"
data_scaling="MiSeq_to_TIRFStation1"
flowcell_side="bottom"

script=$(basename $0)
script_name=${script%.*}
log_dir=$base_dir/$exp_dir/$script_name"Logs"
log_file_prefix="t"
log_file_suffix=".log"
err_file_suffix=".err"

mkdir -p $log_dir


## Quantification using SLURM scheduler
echo "Submitting jobs via SLURM..."
for i in $timepoints
do

    image_subdir=$image_dir/$image_subdir_prefix$i
    log_file=$log_dir/$log_file_prefix$i$log_file_suffix
    err_file=$log_dir/$log_file_prefix$i$err_file_suffix

    start_time=$SECONDS
    echo "Starting quantification for $exp_dir at timepoint $i..."

    srun python2.7 $HOME/scripts/array_tools/CPscripts/quantifyAllTiles.py \
	-rd $read_dir -td $tile_dir -id $image_subdir -f $filter -ftd $filtered_tile_dir \
	-rod $image_subdir/$roff_dir -fd $image_subdir/"CPfluor" \
	-n $num_cores -rs $reg_subset1 -rs $reg_subset2 -sf $data_scaling -sd $flowcell_side -gv $gv_dir -s \
	1> $log_file 2> $err_file 
    
    duration=$(( SECONDS - start_time )) 
    echo "Done with quantification for $exp_dir at timepoint $i. Duration: $duration" | tee -a $log_file

done
