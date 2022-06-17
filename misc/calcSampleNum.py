"""
Calculate the number of MC samples
It should equal to the number of reocnstructed lljj after 137/fb
Input files are produced by extractRootFiles.sh
"""

import ROOT
import numpy as np
from root_numpy import tree2array

ROOT.ROOT.EnableImplicitMT

def main():
    # Names within ttree
    treeFolder = "analysis/allEvents/"
    treeName = "massData"

    weightBranch =  "weight"


    #LOADING THE TTREE
    fileNames = [ "DY100to200.root",
                 "DY200to400.root",
                 "DY400to600.root",
                 "DY600to800.root",
                 "DY800to1200.root",
                 "DY1200to2500.root",
                 "DY2500toInf.root",
                 "TTTo2L2Nu.root"]

    # weight related
    crossSections = [147.4,  40.99, 5.678, 2.198, 0.6304, 0.1514, 0.003565, 88.29]
    counts2 = [2751187, 962195, 1070454, 8292957, 2673066, 596079, 399492, 79140880]
    lumi = 137 #/fb

    sampleNum = 0
    for fileName, xSec, count2 in zip(fileNames, crossSections, counts2):
    	rootfile= ROOT.TFile.Open(fileName, "read")

    	massTree = rootfile.Get(treeFolder+treeName)

    	weightArray = tree2array(massTree, branches=weightBranch)

        allInteractionNum = xSec*lumi
        recoRatio = weightArray.sum()/count2
        sampleNum += recoRatio*allInteractionNum

    sampleNum = int(sampleNum)

    print("You should generate",sampleNum,"samples")

if(__name__=="__main__"):
    main()
