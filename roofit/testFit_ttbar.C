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
  RooRealVar* lljjRecoMass = new RooRealVar("lljjRecoMass", "lljjRecoMass", 0, 3000);
  RooRealVar* ljjRecoMass_Res = new RooRealVar("ljjRecoMass_Res", "ljjRecoMass_Res", 0, 2500);
  RooRealVar* ljjRecoMass_SpRes = new RooRealVar("ljjRecoMass_SpRes", "ljjRecoMass_SpRes", 0, 2500);

  RooDataSet ds1("ds1", "ds1",
                RooArgSet(*lljjRecoMass, *ljjRecoMass_Res, *ljjRecoMass_SpRes),
                ImportFromFile("ttTest.root","analysis/bgRecoMass"));

  RooPlot *frame1 = WR_RecoMass->frame(Title("ttbar lljj Reco Mass (top 2 pT lepton)"));
  ds1.plotOn(frame1, Binning(128));
  RooPlot *frame2 = N_RecoMass_Match->frame(Title("ttbar ljj Mass, Reco by Resolved NN"));
  ds1.plotOn(frame2, Binning(128));
  RooPlot *frame3 = N_RecoMass_NN->frame(Title("ttbar ljj Mass, Reco by SuperResolved NN"));
  ds1.plotOn(frame3, Binning(128));

  // preparing the signal distribution
  RooAddPdf* WR_pdf = DoubleCB(lljjRecoMass);

  // fit distribution to data
  WR_pdf->fitTo(ds1);

  // Draw ntuples
  WR_pdf->plotOn(frame1);




  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 1000, 800);
  c->Divide(2,2);
  c->cd(1);
  frame1->Draw();
  c->cd(2);
  frame2->Draw();
  c->cd(3);
  frame3->Draw();
}

RooAddPdf* DoubleCB(RooRealVar* rrv_x)
{
  RooRealVar* rrv_mean_CB = new RooRealVar("rrv_mean_CB", "rrv_mean_CB", 500, 200, 600);
  RooRealVar* rrv_sigma_CB = new RooRealVar("rrv_sigma_CB", "rrv_sigma_CB", 200, 100, 400);
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
