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

void testFit_WR2000()
{
  // importing ntuples into RooDataSet
  RooRealVar* WR_RecoMass_ee = new RooRealVar("WR_RecoMass_ee", "WR_RecoMass_ee", 0, 3500);
  RooRealVar* WR_RecoMass_mumu = new RooRealVar("WR_RecoMass_mumu", "WR_RecoMass_mumu", 0, 3500);

  RooDataSet ds1("ds1", "ds1",
                RooArgSet(*WR_RecoMass_ee),
                ImportFromFile("../WR2000_N1400/out_WR2000N1400_1.root","analysis/WR_RecoMass_ee"));

  RooDataSet ds2("ds2", "ds2",
                RooArgSet(*WR_RecoMass_mumu),
                ImportFromFile("../WR2000_N1400/out_WR2000N1400_1.root","analysis/WR_RecoMass_mumu"));

  RooPlot *frame1 = WR_RecoMass_ee->frame(Title("2000 GeV WR Mass, Reco by Matching ee"));
  ds1.plotOn(frame1, Binning(128));

  RooPlot *frame2 = WR_RecoMass_mumu->frame(Title("2000 GeV N Mass, Reco by Matching mumu"));
  ds2.plotOn(frame2, Binning(128));


  // preparing the signal distribution
  RooAddPdf* WR_ee_pdf = DoubleCB(WR_RecoMass_ee);
  RooAddPdf* WR_mumu_pdf = DoubleCB(WR_RecoMass_mumu);

  // fit distribution to data
  RooFitResult *r1 = WR_ee_pdf->fitTo(ds1, Save(), Range(1100,2700));
  RooFitResult *r2 = WR_mumu_pdf->fitTo(ds2, Save(), Range(1100,2700));

  // Draw ntuples
  WR_ee_pdf->plotOn(frame1);
  WR_mumu_pdf->plotOn(frame2);

  // preparing the single CB distribution
  // ee
  RooRealVar m0_ee("m0_ee","m0 for ee",2000, 1500, 2500);
  RooRealVar sigma_ee("sigma_ee","sigma for ee", 200, 50, 2000);
  RooRealVar alpha_ee("alpha_ee", "alpha for ee", 2.0, 0., 200.0);
  RooRealVar n_ee("n_ee","n for ee", 2.0, 0.0, 400.0);
  RooCBShape cb_ee("signal_ee", "cb signal for ee",
                *WR_RecoMass_ee,
                m0_ee, sigma_ee, alpha_ee, n_ee);
  // mumu
  RooRealVar m0_mumu("m0_mumu","m0 for mumu",2000, 1500, 2500);
  RooRealVar sigma_mumu("sigma_mumu","sigma for mumu", 200, 50, 2000);
  RooRealVar alpha_mumu("alpha_mumu", "alpha for mumu", 2.0, 0., 200.0);
  RooRealVar n_mumu("n_mumu","n for mumu", 2.0, 0.0, 400.0);
  RooCBShape cb_mumu("signal_mumu", "cb signal for mumu",
                *WR_RecoMass_mumu,
                m0_mumu, sigma_mumu, alpha_mumu, n_mumu);

  // fit distribution to data
  RooFitResult *r3 = cb_ee.fitTo(ds1, Save(), Range(1100,2700));
  RooFitResult *r4 = cb_mumu.fitTo(ds2, Save(), Range(1100,2700));

  // Draw ntuples
  cb_ee.plotOn(frame1, LineColor(kRed));
  cb_mumu.plotOn(frame2, LineColor(kRed));


  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 1000, 800);
  c->Divide(1,2);
  c->cd(1);
  frame1->Draw();
  c->cd(2);
  frame2->Draw();
}

RooAddPdf* DoubleCB(RooRealVar* rrv_x)
{
  RooRealVar* rrv_mean_CB = new RooRealVar("rrv_mean_CB", "rrv_mean_CB", 2000, 1500, 2500);
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
