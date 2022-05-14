# @Author: Billy Li <billyli>
# @Date:   04-23-2022
# @Email:  li000400@umn.edu
# @Last modified by:   billyli
# @Last modified time: 05-12-2022


"""
Combines all dy ttrees into a single root file
"""
import ROOT
import numpy as np
from root_numpy import tree2array
#import tdrstyle

ROOT.ROOT.EnableImplicitMT()


def main():
    # Names within ttree
    treeFolder = "analysis/allEvents/"
    treeName = "massData"
    wrMassBranch = "WRMass"
    SRmassBranch = "superResolvedNNMass"
    RmassBranch = "resolvedNNMass"
    correctMassBranch = "correctNMass"
    incorrectMassBranch = "incorrectNMass"
    leadMassBranch = "leadNMass"
    subleadMassBranch = "subNMass"
    weightBranch = "weight"

    # reco mass ntuples
    analysisFolder = "analysis/"

    mumujjNtupleName = "bg_mumujjRecoMass"
    mumujjMassBranch = "bg_mumujjRecoMass"
    mumujjEventWeightBranch = "eventWeight"

    eejjNtupleName = "bg_eejjRecoMass"
    eejjMassBranch = "bg_eejjRecoMass"
    eejjEventWeightBranch = "eventWeight"


    # LOADING THE TTREE
    fileNames = ["DY100to200.root",
                 "DY200to400.root",
                 "DY400to600.root",
                 "DY600to800.root",
                 "DY800to1200.root",
                 "DY1200to2500.root",
                 "DY2500toInf.root"]

    crossSections = [147.4,  40.99, 5.678, 2.198, 0.6304, 0.1514, 0.003565]
    counts2 = [2751187, 962195, 1070454, 8292957, 2673066, 596079, 399492]

    # make new root file with new tree and new ntuples
    file = ROOT.TFile("fullDY.root", 'recreate')
    tree = ROOT.TTree("fullDY", "fullDY")
    mumujjNtuple_new = ROOT.TNtuple("invm_mumujj", "invm reco from DY mumujj", "invm_mumujj:rowWeight")
    eejjNtuple_new = ROOT.TNtuple("invm_eejj", "invm reco from DY eejj", "invm_eejj:rowWeight")
    mumujjLogMassNtuple_new = ROOT.TNtuple("log(invm)_mumujj", "log invm reco from DY mumujj", "log_invm_mumujj:rowWeight")
    eejjLogMassNtuple_new = ROOT.TNtuple("log(invm)_eejj", "log invm reco from DY eejj", "log_invm_eejj:rowWeight")

    # create 1 dimensional float arrays as fill variables, in this way the float
    # array serves as a pointer which can be passed to the branch

    # create some random numbers, assign them into the fill variables and call Fill()

    WRMass = np.zeros(1, dtype=float)
    resolvedNNMass = np.zeros(1, dtype=float)
    superResolvedNNMass = np.zeros(1, dtype=float)
    correctNMass = np.zeros(1, dtype=float)
    treeWeight = np.zeros(1, dtype=float)
    treeWeight2 = np.zeros(1, dtype=float)
    incorrectNMass = np.zeros(1, dtype=float)
    leadNMass = np.zeros(1, dtype=float)
    subNMass = np.zeros(1, dtype=float)

    tree.Branch("WRMass", WRMass, "WRMass/D")
    tree.Branch("resolvedNNMass", resolvedNNMass, "resolvedNNMass/D")
    tree.Branch("superResolvedNNMass", superResolvedNNMass,
                "superResolvedNNMass/D")
    tree.Branch("correctNMass", correctNMass, "correctNMass/D")
    tree.Branch("incorrectNMass", incorrectNMass, "incorrectNMass/D")
    tree.Branch("leadNMass", leadNMass, "leadNMass/D")
    tree.Branch("subNMass", subNMass, "subNMass/D")
    tree.Branch("weight", treeWeight, "weight/D")
    tree.Branch("weight2", treeWeight2, "weight2/D")

    for fileName, xSec, count2 in zip(fileNames, crossSections, counts2):
        print(fileName)
        rootfile = ROOT.TFile.Open(fileName, "read")

        massTree = rootfile.Get(treeFolder+treeName)

        countHisto = rootfile.Get(treeFolder+"countHisto")
        counts = countHisto.GetBinContent(1)
        print(counts)

        WRmassArray = tree2array(massTree, branches=wrMassBranch)
        SRmassArray = tree2array(massTree, branches=SRmassBranch)
        RmassArray = tree2array(massTree, branches=RmassBranch)
        correctMassArray = tree2array(massTree, branches=correctMassBranch)
        incorrectMassArray = tree2array(massTree, branches=incorrectMassBranch)
        leadMassArray = tree2array(massTree, branches=leadMassBranch)
        subleadMassArray = tree2array(massTree, branches=subleadMassBranch)
        weightArray = tree2array(massTree, branches=weightBranch)
        weightArray2 = weightArray*xSec/count2
        weightArray = weightArray*xSec/counts

        print(WRmassArray.shape)
        print(WRmassArray.shape[0])
        for i in range(WRmassArray.shape[0]):
            WRMass[0] = WRmassArray[i]
            resolvedNNMass[0] = RmassArray[i]
            superResolvedNNMass[0] = SRmassArray[i]
            correctNMass[0] = correctMassArray[i]
            treeWeight[0] = weightArray[i]
            treeWeight2[0] = weightArray2[i]
            incorrectNMass[0] = incorrectMassArray[i]
            leadNMass[0] = leadMassArray[i]
            subNMass[0] = subleadMassArray[i]
            tree.Fill()

        # fill bgRecoMass Ntuple
        mumujjNtuple = rootfile.Get(analysisFolder+mumujjNtupleName)
        eejjNtuple = rootfile.Get(analysisFolder+eejjNtupleName)

        mumujjMassArray = tree2array(mumujjNtuple, branches=mumujjMassBranch)
        mumujjEventWeightArray = tree2array(mumujjNtuple, branches=mumujjEventWeightBranch)
        mumujjRowWeightArray = mumujjEventWeightArray*xSec/count2

        eejjMassArray = tree2array(eejjNtuple, branches=eejjMassBranch)
        eejjEventWeightArray = tree2array(eejjNtuple, branches=eejjEventWeightBranch)
        eejjRowWeightArray = eejjEventWeightArray*xSec/count2

        for i in range(mumujjMassArray.shape[0]):
            mumujjNtuple_new.Fill(mumujjMassArray[i], mumujjRowWeightArray[i])
            mumujjLogMassNtuple_new.Fill(np.log(mumujjMassArray[i]), mumujjRowWeightArray[i])
        for i in range(eejjMassArray.shape[0]):
            eejjNtuple_new.Fill(eejjMassArray[i], eejjRowWeightArray[i])
            eejjLogMassNtuple_new.Fill(np.log(eejjMassArray[i]), eejjRowWeightArray[i])

    mumujjMassArray_new = tree2array(mumujjNtuple_new, branches="invm_mumujj")
    mumujjRowWeightArray_new = tree2array(mumujjNtuple_new, branches="rowWeight")
    np.save("../../python_analysis/DY_mumujjMassArray.npy", mumujjMassArray_new)
    np.save("../../python_analysis/DY_mumujjRowWeightArray.npy", mumujjRowWeightArray_new)

    # write the tree into the output file and close the file
    file.Write()
    file.Close()


if(__name__ == "__main__"):
    main()
