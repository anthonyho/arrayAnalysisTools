#!/usr/bin/env python
#
# Make a summary of a CPfluor file
#
# Anthony Ho, ahho@stanford.edu, 1/15/2015
# Last update 1/15/2015


# Import libraries
import os, sys
import argparse
from itertools import izip
from numpy import pi
import numpy as np
import matplotlib.pyplot as plt
import plotlib


# Print analysis summaries
def printSummary(text, number, totalNumber):
    print text, str(number), "("+str(round(float(number)/totalNumber*100, 2))+"%)"
    return


# Show histograms
def showHistograms(plotTitle, data1, xlabel1, data2, xlabel2, data3, xlabel3):

    fig = plt.figure(figsize=(24, 9))
    numBins = 50

    data = [data1, data2, data3]
    xLabels = [xlabel1, xlabel2, xlabel3]
    ax = []

    for i in range(0, 3):
        ax.append(fig.add_subplot(1, 3, i+1))
        ax[i].hist(data[i], numBins, color='#3F5D7D', alpha=0.8)

        ax[i].set_xlabel(xLabels[i])
        plotlib.makepretty(ax[i])

    ax[0].set_ylabel("Count")
    ax[1].set_title(plotTitle, y=1.05)

    return


def main():

    # Get options and arguments from command line
    parser = argparse.ArgumentParser(description="summarize a CPfluor file")
    parser.add_argument('-p', '--plot', action='store_true', default=False, help="option to plot figures")
    parser.add_argument('CPfluorPath', help="CPfluor file to be analyzed")
    parser.add_argument('filterName', help="name of the filter to be used")
    parser.add_argument('filteredFilePath', help="location of filtered tile file")
    parser.add_argument('outputDir', help="output directory to store the analyzed data and figure")
    args = parser.parse_args()

    # Initialization
    numSeqs = 0
    numFittedSeqs = 0
    twoPi = 2 * pi

    allIntFluor = []
    allAmp = []
    allSigma = []

    CPfluorBaseName = os.path.basename(args.CPfluorPath)

    print "Analyzing "+CPfluorBaseName+"..."

    # Output files
    outputFilePath = os.path.join(args.outputDir, CPfluorBaseName+".sum")
    figurePath = os.path.join(args.outputDir, CPfluorBaseName+".sum.eps")

    # Go through the CPfluor file
    with open(args.CPfluorPath, 'r') as r1, open(args.filteredFilePath, 'r') as r2, open(outputFilePath, 'w') as w:
        for line1, line2 in izip(r1, r2):

            # Read the CPfluor file
            CPfluorLine = line1.rstrip().split(':')
            fitted = CPfluorLine[7]

            # Read the filter file
            filterTags = line2.rstrip().split('\t')[1]
            lastTag = filterTags.rstrip().split(':')[-1]

            # If the last element matches with the desired filter
            if lastTag == args.filterName:

                numSeqs += 1

                if fitted == '1':

                    numFittedSeqs += 1

                    # Parse data
                    amp = CPfluorLine[8]
                    sigma = CPfluorLine[9]

                    # Compute integrated fluorescence intensity
                    intFluor = str(twoPi*float(amp)*float(sigma)**2)

                    # Write to output file
                    w.write(intFluor+'\t'+amp+'\t'+sigma+'\n')

                    # Append to lists
                    allIntFluor.append(intFluor)
                    allAmp.append(amp)
                    allSigma.append(sigma)

    # Print summaries
    printSummary("  Number of clusters with the correct filter:", numSeqs, numSeqs)
    printSummary("  Number of clusters fitted:", numFittedSeqs, numSeqs)

    # Plot and save histograms
    showHistograms(CPfluorBaseName+" (fitted = " + str(numFittedSeqs) + "/" + str(numSeqs) + " = "
                   + str(round(float(numFittedSeqs)/numSeqs*100, 2)) + "%)",
                   np.array(map(float, allIntFluor)), "Integrated fluorescence",
                   np.array(map(float, allAmp)), "Amplitude",
                   np.array(map(float, allSigma)), "Sigma")
    plt.savefig(figurePath, bbox_inches='tight')

    # Show the histograms if requested
    if args.plot:
        plt.show()

    r1.close()
    r2.close()
    w.close()

    return 1

if __name__ == "__main__":
    main()
