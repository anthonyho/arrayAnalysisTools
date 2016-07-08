#!/bin/bash
#
# Preprocess images on Greendragon
#
# Anthony Ho, ahho@stanford.edu, 6/30/2016


input_dir="20160609_CNv3_3_ANUUN_01"
exp_dir_prefix="05_sm_kd_"
for i in {1..10}
do
    $HOME/analysis/imaging/20160609_CNv3_3_ANUUN_01/scripts/01_preprocessImageSeries.sh $input_dir $exp_dir_prefix$i 
done


input_dir="20160614_CNv3_3_ANUUN_02"
exp_dir_prefix="05_sm_kd_"
for i in {11..18}
do
    $HOME/analysis/imaging/20160609_CNv3_3_ANUUN_01/scripts/01_preprocessImageSeries.sh $input_dir $exp_dir_prefix$i 
done


input_dir="20160609_CNv3_3_ANUUN_01"
exp_dir_prefix="05_cm_"
for i in {1..2}
do
    $HOME/analysis/imaging/20160609_CNv3_3_ANUUN_01/scripts/01_preprocessImageSeries.sh $input_dir $exp_dir_prefix$i 
done


input_dir="20160614_CNv3_3_ANUUN_02"
exp_dir_prefix="05_cm_"
for i in {3..18}
do
    $HOME/analysis/imaging/20160609_CNv3_3_ANUUN_01/scripts/01_preprocessImageSeries.sh $input_dir $exp_dir_prefix$i 
done

