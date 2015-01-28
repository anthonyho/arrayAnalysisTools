#!/usr/bin/env python
#
# Alignemnt and annotation tools
#
# Anthony Ho, ahho@stanford.edu, 1/28/2015
# Last update 1/28/2015


## Import libraries


## Hamming distance calculator, assuming equal length
def hammingDistance(s1,s2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

## Stupid way to find the first (and only) deletion
def findIndel(q,refU):
    for i in range(0,len(refU)):
        if q[i] != refU[i]:
            bail = hammingDistance(q[i+1:i+9],refU[i+2:i+10]) <= 2
            if bail:
                return i+1
    return -1000

