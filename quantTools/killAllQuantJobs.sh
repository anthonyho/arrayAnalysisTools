#!/bin/bash
#
# Bash script to kill all the quantification jobs running under my name all at once
#
# Anthony Ho, ahho@stanford.edu, 1/15/2015
# Last update 1/15/20015

user=$(whoami)

read -p "Are you sure you want to delete all quantification jobs under your account? (y/n)? " answer
case ${answer:0:1} in
    y|Y )
        # Kill the parent quantifyAllTiles.sh script
	kill $(ps aux | grep -E "$user.*[b]ash.*quantifyAllTiles" | awk '{print $2}')
	# Kill the quantifyAllTiles.py script
	kill $(ps aux | grep -E "$user.*[p]ython.*quantifyAllTiles" | awk '{print $2}')
	# Kill the MATLAB instances
	kill $(ps aux | grep -E "$user.*[m]atlab.*AnalyseImage" | awk '{print $2}')
    ;;
    * )
    ;;
esac
