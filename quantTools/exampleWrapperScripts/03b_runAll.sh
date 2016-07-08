#!/bin/bash
#
# Merge series on Greendragon
#
# Anthony Ho, ahho@stanford.edu, 7/8/2016


exp_dir_prefix="05_cm_"
tiles=$(echo {001,002,004,006,008,009,010,011,014})
num_timepoints=4
for i in {1..18}
do
    $HOME/analysis/imaging/20160609_CNv3_3_ANUUN_01/scripts/03_mergeSeries.sh $exp_dir_prefix$i "$tiles" $num_timepoints
done
