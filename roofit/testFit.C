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
  RooRealVar WR_RecoMass("WR_RecoMass", "WR_RecoMass", 0, 1500);
  RooRealVar N_RecoMass("N_RecoMass", "N_RecoMass", 0, 1200);
  RooDataSet ds("ds", "ds",
                RooArgSet(WR_RecoMass, N_RecoMass),
                ImportFromFile("test.root","analysis/WR_N_Mass_1"));
  ds.get(20)->Print("V");

  // Draw ntuples
  RooPlot *frame1 = WR_RecoMass.frame(Title("WR Reco Mass:Number:Reco Mass (GeV)"));
  ds.plotOn(frame1, Binning(32));

  RooPlot *frame2 = N_RecoMass.frame(Title("N Reco Mass:Number:Reco Mass (GeV)"));
  ds.plotOn(frame2, Binning(32));

  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 1000, 800);
  c->Divide(1,2);
  c->cd(1);
  frame1->Draw();
  c->cd(2);
  frame2->Draw();
}
