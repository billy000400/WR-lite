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
  RooRealVar N_RecoMass("N_RecoMass", "N_RecoMass", 0, 2000);
  RooDataSet ds("ds", "ds",
                RooArgSet(WR_RecoMass, N_RecoMass),
                ImportFromFile("test.root","analysis/WR_N_Mass_1"));
}
