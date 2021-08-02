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
  RooRealVar WR_RecoMass("WR_mass", "WR_RecoMass", 0, 3000);
  RooRealVar N_RecoMass("N_mass", "N_RecoMass", 0, 1500);
  RooDataSet ds1("ds1", "ds1",
                RooArgSet(WR_RecoMass, N_RecoMass),
                ImportFromFile("/data/cmszfs1/user/li000400/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/roofit/test.root","analysis/WR_N_Mass_1"));

  RooPlot *frame1 = WR_RecoMass.frame(Title("WR Reco Mass"));
  ds1.plotOn(frame1, Binning(32));

  // preparing the signal distribution
  RooRealVar m0("m0","m0",1000, 1000, 1200);
  RooRealVar sigma("sigma","sigma", 80, 30, 230);
  RooRealVar tail("tail", "tail", 10, -300, 300);

  RooNovosibirsk novo("novo", "novo pdf", WR_RecoMass, m0, sigma, tail);
  novo.plotOn(frame1);

  // fit distribution to data
  novo.fitTo(ds1);
  novo.plotOn(frame1, LineColor(kRed));
  // Draw ntuples


  RooPlot *frame2 = N_RecoMass.frame(Title("N Reco Mass"));
  ds1.plotOn(frame2, Binning(32));

  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 1000, 800);
  c->Divide(1,2);
  c->cd(1);
  frame1->Draw();
  c->cd(2);
  frame2->Draw();
}
