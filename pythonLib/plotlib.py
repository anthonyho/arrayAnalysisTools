#!/usr/bin/env python
#
# Python module containing some handy plotting tools
#
# Anthony Ho, ahho@stanford.edu, 1/15/2015
# Last update 1/15/2015

import matplotlib.pyplot as plt

def makepretty(ax):
    
    # Change figure face color to white
    fig = ax.figure
    fig.patch.set_facecolor('white')

    # Eliminate axis ticks on the right and top
    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left()  

    # Remove the top and right plot frame lines
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    
    # Set label and title font size
    ax.xaxis.label.set_fontsize(16)
    ax.yaxis.label.set_fontsize(16)
    ax.title.set_fontsize(16)

    return
