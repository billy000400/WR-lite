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

    mumujjNtupleName = "bg_mumujjRecoMass"
    mumujjEventWeightBranch = "eventWeight"

    eejjNtupleName = "bg_eejjRecoMass"
    eejjEventWeightBranch = "eventWeight"

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
    lumi = 137000 #/pb

    sampleNum_mumujj = 0.0
    sampleNum_eejj = 0.0

    for fileName, xSec, count2 in zip(fileNames, crossSections, counts2):
    	rootfile= ROOT.TFile.Open(fileName, "read")

    	mumujjNtuple = rootfile.Get(analysisFolder+mumujjNtupleName)
        eejjNtuple = rootfile.Get(analysisFolder+eejjNtupleName)

        mumujjEventWeightArray = tree2array(mumujjNtuple, branches=mumujjEventWeightBranch)

        eejjEventWeightArray = tree2array(eejjNtuple, branches=eejjEventWeightBranch)

        allInteractionNum = xSec*lumi
        mumujjRecoRatio = mumujjEventWeightArray.sum()/count2
        eejjRecoRatio = eejjEventWeightArray.sum()/count2

        sampleNum_mumujj += allInteractionNum*mumujjRecoRatio
        sampleNum_eejj += allInteractionNum*eejjRecoRatio

    sampleNum_mumujj = int(sampleNum_mumujj)
    sampleNum_eejj = int(sampleNum_eejj)

    print "You should generate",sampleNum_mumujj,"samples for eejj"
    print "You should generate",sampleNum_eejj,"samples for mumujj"

if(__name__=="__main__"):
    main()
