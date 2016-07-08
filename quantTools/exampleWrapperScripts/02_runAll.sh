#!/bin/bash
#
# Quantify images on Sherlock
#
# Anthony Ho, ahho@stanford.edu, 6/30/2016


exp_dir_prefix="05_sm_kd_"
for i in {1..18}
do
    sbatch $HOME/analysis/imaging/20160609_CNv3_3_ANUUN_01/scripts/02_quantifyAllTiles.sbatch $exp_dir_prefix$i
done


exp_dir_prefix="05_cm_"
for i in {1..18}
do
    sbatch $HOME/analysis/imaging/20160609_CNv3_3_ANUUN_01/scripts/02_quantifyAllTiles.sbatch $exp_dir_prefix$i
done

