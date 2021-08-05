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

RooAddPdf* DoubleCB(RooRealVar* rrv_x);

void testFit_DoubleCB()
{
  // importing ntuples into RooDataSet
  RooRealVar* WR_RecoMass = new RooRealVar("WR_RecoMass", "WR_RecoMass", 0, 3000);
  RooRealVar* N_RecoMass = new RooRealVar("N_RecoMass", "N_RecoMass", 0, 1500);
  RooDataSet ds1("ds1", "ds1",
                RooArgSet(WR_RecoMass, N_RecoMass),
                ImportFromFile("test.root","analysis/WR_N_RecoMass"));

  RooPlot *frame1 = WR_RecoMass->frame(Title("WR Reco Mass"));
  ds1->plotOn(frame1, Binning(128));

  // preparing the signal distribution
  RooAddPdf* DoubleCB(RooRealVar);

  // fit distribution to data
  DoubleCB->fitTo(ds1);

  // Draw ntuples
  DoubleCB->plotOn(frame1);


  RooPlot *frame2 = N_RecoMass->frame(Title("N Reco Mass"));
  ds1.plotOn(frame2, Binning(128));

  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 1000, 800);
  c->Divide(1,2);
  c->cd(1);
  frame1->Draw();
  c->cd(2);
  frame2->Draw();
}

RooAddPdf* DoubleCB(RooRealVar* rrv_x)
{
  RooRealVar* rrv_mean_CB = new RooRealVar("rrv_mean_CB", "rrv_mean_CB", 1000, 900, 1200);
  RooRealVar* rrv_sigma_CB = new RooRealVar("rrv_sigma_CB", "rrv_sigma_CB", 100, 50, 300);
  RooRealVar* rrv_tail_CB_I = new RooRealVar("rrv_tail_CB_I", "rrv_tail_CB_I", 2,0., 40);
  RooRealVar* rrv_tail_CB_II = new RooRealVar("rrv_tail_CB_II", "rrv_tail_CB_II", -2., -40., 0.);

  RooRealVar* rrv_normalization_CB_I = new RooRealVar("rrv_normalization_CB_I", "rrv_normalization_CB_I", 2, 0., 400);
  RooRealVar* rrv_frac_CB = new RooRealVar("rrv_frac_CB", "rrv_frac_CB", 0.5);


  RooCBShape* Crystal_Ball_I = new RooCBShape("CrystalBall_I", "CrystalBall_I",
                                              *rrv_x,
                                              *rrv_mean_CB,*rrv_sigma_CB,*rrv_tail_CB_I,*rrv_normalization_CB_I);

  RooCBShape* Crystal_Ball_II = new RooCBShape("CrystalBall_II", "CrystalBall_II",
                                              *rrv_x,
                                              *rrv_mean_CB,*rrv_sigma_CB,*rrv_tail_CB_II,*rrv_normalization_CB_I);

  RooAddPdf* model_pdf = new RooAddPdf("model_pdf", "model_pdf",
                                      RooArgList(*Crystal_Ball_I,*Crystal_Ball_II),
                                      RooArgList(*rrv_frac_CB));


  return model_pdf;
}
