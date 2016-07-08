#!/bin/bash
#
# Merge series on Greendragon
#
# Anthony Ho, ahho@stanford.edu, 7/8/2016


exp_dir_prefix="05_sm_kd_"
tiles=$(echo {001,002,004,006,008,009,010,011,014})
num_timepoints=18
for i in {{1..5},13}
do
    $HOME/analysis/imaging/20160609_CNv3_3_ANUUN_01/scripts/03_mergeSeries.sh $exp_dir_prefix$i "$tiles" $num_timepoints
done

exp_dir_prefix="05_sm_kd_"
tiles=$(echo {001,002,004,006,008,009,010,011,014})
num_timepoints=20
for i in {{6..10},12,14,18}
do
    $HOME/analysis/imaging/20160609_CNv3_3_ANUUN_01/scripts/03_mergeSeries.sh $exp_dir_prefix$i "$tiles" $num_timepoints
done

exp_dir_prefix="05_sm_kd_"
tiles=$(echo {001,002,004,006,008,010,011,014})
num_timepoints=18
i=11
$HOME/analysis/imaging/20160609_CNv3_3_ANUUN_01/scripts/03_mergeSeries.sh $exp_dir_prefix$i "$tiles" $num_timepoints

exp_dir_prefix="05_sm_kd_"
tiles=$(echo {001,002,004,008,009,010,011,014})
num_timepoints=20
i=15
$HOME/analysis/imaging/20160609_CNv3_3_ANUUN_01/scripts/03_mergeSeries.sh $exp_dir_prefix$i "$tiles" $num_timepoints

exp_dir_prefix="05_sm_kd_"
tiles=$(echo {001,002,004,006,009,010,011,014})
num_timepoints=18
i=16
$HOME/analysis/imaging/20160609_CNv3_3_ANUUN_01/scripts/03_mergeSeries.sh $exp_dir_prefix$i "$tiles" $num_timepoints

exp_dir_prefix="05_sm_kd_"
tiles=$(echo {001,002,004,006,008,009,010,011,014})
num_timepoints=14
i=17
$HOME/analysis/imaging/20160609_CNv3_3_ANUUN_01/scripts/03_mergeSeries.sh $exp_dir_prefix$i "$tiles" $num_timepoints

