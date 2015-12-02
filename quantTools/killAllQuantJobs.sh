#!/bin/bash
#
# Bash script to kill all the quantification jobs running under my name all at once
#
# Anthony Ho, ahho@stanford.edu, 1/15/2015
# Last update 12/2/20015

user=$(whoami)

read -p "Are you sure you want to delete all quantification jobs under your account? (y/n)? " answer
case ${answer:0:1} in
    y|Y )
        # Kill the parent quantifyAllTiles.sh script
	bash_jobs=$(ps aux | grep -E "$user.*[b]ash.*quantifyAllTiles" | awk '{print $2}')
	echo "Found bash jobs" $bash_jobs
	if [ ! -z "$bash_jobs" ]
	then
	    echo "Killing bash jobs" $bash_jobs "now..."
	    kill $bash_jobs
	fi
	# Kill the quantifyAllTiles.py script
	python_jobs=$(ps aux | grep -E "$user.*[p]ython.*quantifyAllTiles" | awk '{print $2}')
	echo "Found python jobs" $python_jobs
	if [ ! -z "$python_jobs" ]
	then
	    echo "Killing python jobs" $python_jobs "now..."
	    kill $python_jobs
	fi
	# Kill the MATLAB instances
	matlab_regJobs=$(ps aux | grep -E "$user.*[m]atlab.*GenerateRegistrationOffsetMap" | awk '{print $2}')
	echo "Found matlab registration jobs" $matlab_regJobs
	if [ ! -z "$matlab_regJobs" ]
	then
	    echo "Killing matlab registration jobs" $matlab_regJobs "now..."
	    kill $matlab_regJobs
	fi
	matlab_quantJobs=$(ps aux | grep -E "$user.*[m]atlab.*AnalyseImage" | awk '{print $2}')
	echo "Found matlab quantification jobs" $matlab_quantJobs
	if [ ! -z "$matlab_quantJobs" ]
	then
	    echo "Killing matlab quantification jobs" $matlab_quantJobs "now..."
	    kill $matlab_quantJobs
	fi
    ;;
    * )
    ;;
esac
