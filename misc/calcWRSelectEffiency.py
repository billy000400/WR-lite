# @Author: Billy Li <billyli>
# @Date:   05-14-2022
# @Email:  li000400@umn.edu
# @Last modified by:   billyli
# @Last modified time: 06-09-2022



"""
Creates new root files containing just the TTrees with output data
"""

import ROOT

import numpy as np
from root_numpy import tree2array

ROOT.ROOT.EnableImplicitMT()

def main():

    # reco mass ntuples
    analysisFolder = "analysis/"
    treeFolder = "analysis/allEvents/"

    mumujjNtupleName = "WR_RecoMass_mumu"
    mumujjMassBranch = "WR_RecoMass_mumu"
    mumujjEventWeightBranch = "eventWeight"

    eejjNtupleName = "WR_RecoMass_ee"
    eejjMassBranch = "WR_RecoMass_ee"
    eejjEventWeightBranch = "eventWeight"

    #TTree file names
    fileNames = ["WR800N200", "WR800N400", "WR800N600", "WR800N700",
                       "WR1000N200", "WR1000N400", "WR1000N600", "WR1000N800",
                       "WR1000N900", "WR1200N400", "WR1200N600", "WR1200N800",
                       "WR1200N1000", "WR1200N1100", "WR1400N400", "WR1400N600",
                       "WR1400N800", "WR1400N1000", "WR1400N1200", "WR1400N1300",
                       "WR1600N400", "WR1600N600", "WR1600N800", "WR1600N1000",
                       "WR1600N1200", "WR1600N1400", "WR1600N1500", "WR1800N400",
                       "WR1800N600", "WR1800N800", "WR1800N1000", "WR1800N1200",
                       "WR1800N1400", "WR1800N1600", "WR1800N1700", "WR2000N600",
                       "WR2000N800", "WR2000N1000", "WR2000N1200", "WR2000N1400",
                       "WR2000N1600", "WR2000N1800", "WR2000N1900"]
    #Cross sections for different files
    crossSections = [14.46, 10.41, 4.351, 1.473,
    	6.083, 5.023, 3.323, 1.256, 0.4125,
    	2.511, 1.947, 1.207, .4296, .1406,
    	1.321, 1.11, .825, .4887, .166, .05455,
    	.7256, .6375,.5188, .3731, .2129, .07031,.02318,
    	.4125, .3718, .3192, .2528, .1771, .09803,.03186, .01069,
    	.2217, .1962, .1649, .1281, .08747, .04720, .01517, .005147]


    for (files, xSec) in zip(fileNames, crossSections):

    	# create 1 dimensional float arrays as fill variables, in this way the float
    	# array serves as a pointer which can be passed to the branch



    	rootfile= ROOT.TFile.Open(files+".root", "read")

    	countHisto = rootfile.Get(treeFolder+"countHisto")
    	counts = countHisto.GetBinContent(1)

        # fill lljjRecoMass Ntuple
        mumujjNtuple = rootfile.Get(analysisFolder+mumujjNtupleName)
        eejjNtuple = rootfile.Get(analysisFolder+eejjNtupleName)


        mumujjEventWeightArray = tree2array(mumujjNtuple, branches=mumujjEventWeightBranch)
        eejjEventWeightArray = tree2array(eejjNtuple, branches=eejjEventWeightBranch)

        branch_ratio_mumu = 0.5
        branch_ratio_ee = 1-branch_ratio_mumu

        ratio_mumu = mumujjEventWeightArray.sum()/(counts*branch_ratio_mumu)
        ratio_ee = eejjEventWeightArray.sum()/(counts*branch_ratio_ee)

        print files, ratio_mumu, ratio_ee



if(__name__ == "__main__"):
    main()
