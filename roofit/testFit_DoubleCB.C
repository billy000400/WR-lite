/**
 * @Author: Billy Li <billyli>
 * @Date:   08-10-2021
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

RooAddPdf* DoubleCB(RooRealVar* rrv_x);

void testFit_DoubleCB()
{
  // importing ntuples into RooDataSet
  RooRealVar* WR_RecoMass = new RooRealVar("WR_RecoMass", "WR_RecoMass", 0, 5000);
  RooRealVar* N_RecoMass_Match = new RooRealVar("N_RecoMass_Match", "N_RecoMass_Match", 0, 3000);
  RooRealVar* N_RecoMass_NN = new RooRealVar("N_RecoMass_NN", "N_RecoMass_NN", 0, 3000);

  RooDataSet ds1("ds1", "ds1",
                RooArgSet(*WR_RecoMass, *N_RecoMass_Match, *N_RecoMass_NN),
                ImportFromFile("../WR2000_N1900/out_WR2000N1900_4.root","analysis/WR_N_RecoMass"));

  RooPlot *frame1 = WR_RecoMass->frame(Title("1000 GeV WR Mass, Reco by Matching"));
  ds1.plotOn(frame1, Binning(128));
  RooPlot *frame2 = N_RecoMass_Match->frame(Title("400 GeV N Mass, Reco by Matching"));
  ds1.plotOn(frame2, Binning(128));
  RooPlot *frame3 = N_RecoMass_NN->frame(Title("400 GeV N Mass, Reco by NN"));
  ds1.plotOn(frame3, Binning(128));

  // preparing the signal distribution
  RooRealVar* rrv_mean_CB = new RooRealVar("rrv_mean_CB", "rrv_mean_CB", 2000, 1000, 3000);
  RooRealVar* rrv_sigma_CB = new RooRealVar("rrv_sigma_CB", "rrv_sigma_CB", 100, 50, 300);
  RooRealVar* rrv_tail_CB_I = new RooRealVar("rrv_tail_CB_I", "rrv_tail_CB_I", 2,0., 40);
  RooRealVar* rrv_tail_CB_II = new RooRealVar("rrv_tail_CB_II", "rrv_tail_CB_II", -2., -40., 0.);

  RooRealVar* rrv_normalization_CB_I = new RooRealVar("rrv_normalization_CB_I", "rrv_normalization_CB_I", 2, 0., 400);
  RooRealVar* rrv_frac_CB = new RooRealVar("rrv_frac_CB", "rrv_frac_CB", 0.5);


  RooCBShape* Crystal_Ball_I = new RooCBShape("CrystalBall_I", "CrystalBall_I",
                                              *WR_RecoMass,
                                              *rrv_mean_CB,*rrv_sigma_CB,*rrv_tail_CB_I,*rrv_normalization_CB_I);

  RooCBShape* Crystal_Ball_II = new RooCBShape("CrystalBall_II", "CrystalBall_II",
                                              *WR_RecoMass,
                                              *rrv_mean_CB,*rrv_sigma_CB,*rrv_tail_CB_II,*rrv_normalization_CB_I);

  RooAddPdf* WR_pdf = new RooAddPdf("model_pdf", "model_pdf",
                                      RooArgList(*Crystal_Ball_I,*Crystal_Ball_II),
                                      RooArgList(*rrv_frac_CB));

  // fit distribution to data
  RooFitResult *r = WR_pdf->fitTo(ds1, Save());
  r->Print();

  std::cout << "The value of rrv_frac_CB is " << rrv_frac_CB.evaluate() << std::endl;

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
  RooRealVar* rrv_mean_CB = new RooRealVar("rrv_mean_CB", "rrv_mean_CB", 2000, 1000, 3000);
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
