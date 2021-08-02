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
  // importing ntuples into RooDataSet
  RooRealVar WR_RecoMass("WR_mass", "WR_RecoMass", 0, 1500);
  RooRealVar N_RecoMass("N_mass", "N_RecoMass", 0, 1200);
  RooDataSet ds("ds", "ds",
                RooArgSet(WR_RecoMass, N_RecoMass),
                ImportFromFile("/data/cmszfs1/user/li000400/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/roofit/test.root","analysis/WR_N_Mass_1"));

  // preparing the signal distribution
  RooRealVar m0("m0","m0",1000, 900, 1100);
  RooRealVar sigma("sigma","sigma", 200, 0, 500);
  RooRealVar alpha("alpha", "alpha", 1, 0, 500);
  RooRealVar n("n","n", 1, 0, 500);
  RooCBShape cb("signal", "cb signal",
                WR_RecoMass,
                m0, sigma, alpha, n);

  // fit distribution to data
  cb.fitTo(ds);


  // Draw ntuples
  RooPlot *frame1 = WR_RecoMass.frame(Title("WR Reco Mass"));
  ds.plotOn(frame1, Binning(32));
  cb.plotOn(frame1);

  RooPlot *frame2 = N_RecoMass.frame(Title("N Reco Mass"));
  ds.plotOn(frame2, Binning(32));

  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 1000, 800);
  c->Divide(1,2);
  c->cd(1);
  frame1->Draw();
  c->cd(2);
  frame2->Draw();
}
