/**
 * @Author: Billy Li <billyli>
 * @Date:   08-05-2021
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 08-18-2021
 */



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
  RooRealVar* WR_RecoMass = new RooRealVar("WR_RecoMass", "WR_RecoMass", 0, 5000);
  RooRealVar* N_RecoMass_Match = new RooRealVar("N_RecoMass_Match", "N_RecoMass_Match", 0, 3000);
  RooRealVar* N_RecoMass_NN = new RooRealVar("N_RecoMass_NN", "N_RecoMass_NN", 0, 3000);

  RooDataSet ds1("ds1", "ds1",
                RooArgSet(*WR_RecoMass, *N_RecoMass_Match, *N_RecoMass_NN),
                ImportFromFile("../WR2000_N1900/out_WR2000N1900_4.root","analysis/WR_N_RecoMass"));

  RooDataSet ds2("ds2", "ds2",
                RooArgSet(*WR_RecoMass, *N_RecoMass_Match, *N_RecoMass_NN),
                ImportFromFile("../WR2000_N1900/out_WR2000N1900_4.root","analysis/WR_N_RecoMass"));

  RooPlot *frame1 = WR_RecoMass->frame(Title("WR Reco Mass"));
  ds1.plotOn(frame1, Binning(32));

  // preparing the signal distribution
  RooRealVar m0("m0","m0",2000, 1900, 2100);
  RooRealVar sigma("sigma","sigma", 80, 30, 230);
  RooRealVar alpha("alpha", "alpha", 0.1, -0.25, 0.25);
  RooRealVar n("n","n", 0.1, -10, 10);

  RooGaussian gauss("gauss", "gaussian pdf",
                  *WR_RecoMass,
                  m0, sigma);
  RooCBShape cb("signal", "cb signal",
                *WR_RecoMass,
                m0, sigma, alpha, n);
  cb.plotOn(frame1);

  // fit distribution to data
  gauss.fitTo(ds1, Range(900,2000));
  gauss.plotOn(frame1, LineColor(kYellow));

  m0.setConstant(kTRUE);
  sigma.setConstant(kTRUE);
  cb.fitTo(ds2);
  cb.plotOn(frame1, LineColor(kRed));

  // Draw ntuples


  RooPlot *frame2 = N_RecoMass_Match.frame(Title("N Reco Mass"));
  ds1.plotOn(frame2, Binning(32));

  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 1000, 800);
  c->Divide(1,2);
  c->cd(1);
  frame1->Draw();
  c->cd(2);
  frame2->Draw();
}
