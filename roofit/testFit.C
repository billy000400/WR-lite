#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
using namespace RooFit;

void testFit()
{
  // importing ntuples
  RooRealVar WR_RecoMass("WR_RecoMass", "WR_RecoMass", 0, 3000);
  RooDataSet ds("ds", "ds", RooArgSet(WR_RecoMass), ImportFromFile("test.root","WR_mass"));
}
