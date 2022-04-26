/**
 * @Author: Billy Li <billyli>
 * @Date:   08-10-2021
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 04-26-2022
 */

// This script is to figure out the best strategy to fit data into a
// double-side (DS) CB distribution. It will be put in testFit.C to compare the result
// of DSCB and single CB

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "RooDSCBShape.h"
using namespace RooFit;

RooDSCBShape* DSCB_init(RooRealVar* rrv_x, double mean, std::string label);
RooDataSet Hist2Pulls(RooHist* pullPlot, std::string label, bool print=false);

void testFit_DSCB(std::string filePath)
{
  //// Extract useful information and set fitting parameters
  // Extract WR and N mean value via the file name
  std::cout << "Openning file " << filePath << std::endl;
  size_t fileNamePos = filePath.find_last_of("/");
  std::string fileName = filePath.substr(fileNamePos+1);
  size_t RPos = fileName.find_last_of("R");
  size_t NPos = fileName.find_last_of("N");
  size_t dotPos = fileName.find_last_of(".");
  double WRGenMean = std::stod(fileName.substr(RPos+1, NPos-RPos));
  double NGenMean = std::stod(fileName.substr(NPos+1, dotPos-NPos));
  std::cout << "Target WR: " << WRGenMean << ", Target N" << NGenMean << std::endl;


  std::string prefix = "../";
  RooRealVar* WR_RecoMass_ee = new RooRealVar("WR_RecoMass_ee", "WR_RecoMass_ee", WRGenMean*0.45, WRGenMean*1.45);
  RooRealVar* WR_RecoMass_mumu = new RooRealVar("WR_RecoMass_mumu", "WR_RecoMass_mumu", WRGenMean*0.45, WRGenMean*1.45);
  RooDataSet ds_WR_RecoMass_ee("ds1", "ds1",
                RooArgSet(*WR_RecoMass_ee),
                ImportFromFile((prefix+filePath).c_str(), "analysis/WR_RecoMass_ee"));
  RooDataSet ds_WR_RecoMass_mumu("ds2", "ds2",
                RooArgSet(*WR_RecoMass_mumu),
                ImportFromFile((prefix+filePath).c_str(), "analysis/WR_RecoMass_mumu"));

  //// Preparing probability distirbution functions for fitting
  // preparing the double CB distributions
  RooDSCBShape* WR_ee_doubleCB = DSCB_init(WR_RecoMass_ee, WRGenMean, "eejj");
  RooDSCBShape* WR_mumu_doubleCB = DSCB_init(WR_RecoMass_mumu, WRGenMean, "mumujj");


  //// fit distribution to data
  RooFitResult *r1 = WR_ee_doubleCB->fitTo(ds_WR_RecoMass_ee, Save(), Range(WRGenMean*0.45,WRGenMean*1.45));
  RooFitResult *r2 = WR_mumu_doubleCB->fitTo(ds_WR_RecoMass_mumu, Save(), Range(WRGenMean*0.45,WRGenMean*1.45));

  //// Prepare frames for plotting
  RooPlot *eeFrame_doubleCB = WR_RecoMass_ee->frame(Title("eejj Double CB"));
  RooPlot *mumuFrame_doubleCB = WR_RecoMass_mumu->frame(Title("mumujj Double CB"));

  //// Plot on frames
  // plot data on frames
  ds_WR_RecoMass_ee.plotOn(eeFrame_doubleCB, Binning(40), DataError(RooAbsData::SumW2));
  ds_WR_RecoMass_mumu.plotOn(mumuFrame_doubleCB, Binning(40), DataError(RooAbsData::SumW2));
  // plot fitted pdfs on frames
  WR_ee_doubleCB->plotOn(eeFrame_doubleCB);
  WR_mumu_doubleCB->plotOn(mumuFrame_doubleCB);

  //// pull related
  // Prepare pulls
  RooRealVar* pullVar = new RooRealVar("pullVar", "pull value", -6, 6);
  std::cout << "Making the pull plots" << std::endl;
  RooHist *eeHist_doubleCBPull = eeFrame_doubleCB->pullHist();
  RooHist *mumuHist_doubleCBPull = mumuFrame_doubleCB->pullHist();
  // Extract pulls from RooHist
  RooDataSet ee2CBPulls = Hist2Pulls(eeHist_doubleCBPull,"eejj", true);
  RooDataSet mumu2CBPulls = Hist2Pulls(mumuHist_doubleCBPull, "mumujj", true);
  // Prepare frame for the pull histograms
  RooPlot* ee2CBPullFrame = pullVar->frame(Title("ee Double CB Pull Hist"));
  RooPlot* mumu2CBPullFrame = pullVar->frame(Title("mumu Double CB pull Hist"));
  // plot pull histograms on frames
  std::cout << "Making the pull histograms" << std::endl;
  ee2CBPulls.plotOn(ee2CBPullFrame, Binning(15));
  mumu2CBPulls.plotOn(mumu2CBPullFrame, Binning(15));

  //// plot the errors of the fitted functions
  //// This needs to be done after the pulls was calculated
  //// otherwise it will interfer the pull calculations
  // WR_ee_doubleCB->plotOn(eeFrame_doubleCB, VisualizeError(*r1, 1, kFALSE));
  // WR_mumu_doubleCB->plotOn(mumuFrame_doubleCB, VisualizeError(*r2, 1, kFALSE));

  //// Draw Frames on TCanvas
  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 600, 600);
  c->Divide(2,2);
  c->cd(1);
  eeFrame_doubleCB->Draw();
  c->cd(2);
  mumuFrame_doubleCB->Draw();
  c->cd(3);
  ee2CBPullFrame->Draw();
  c->cd(4);
  mumu2CBPullFrame->Draw();

  //// old code
  // importing ntuples into RooDataSet
  // RooRealVar* WR_RecoMass = new RooRealVar("WR_RecoMass", "WR_RecoMass", 0, 5000);
  // RooRealVar* N_RecoMass_Match = new RooRealVar("N_RecoMass_Match", "N_RecoMass_Match", 0, 3000);
  // RooRealVar* N_RecoMass_NN = new RooRealVar("N_RecoMass_NN", "N_RecoMass_NN", 0, 3000);
  //
  // RooDataSet ds1("ds1", "ds1",
  //               RooArgSet(*WR_RecoMass, *N_RecoMass_Match, *N_RecoMass_NN),
  //               ImportFromFile("../WR2000_N1900/out_WR2000N1900_4.root","analysis/WR_N_RecoMass"));
  //
  // RooPlot *frame1 = WR_RecoMass->frame(Title("1000 GeV WR Mass, Reco by Matching"));
  // ds1.plotOn(frame1, Binning(128));
  // RooPlot *frame2 = N_RecoMass_Match->frame(Title("400 GeV N Mass, Reco by Matching"));
  // ds1.plotOn(frame2, Binning(128));
  // RooPlot *frame3 = N_RecoMass_NN->frame(Title("400 GeV N Mass, Reco by NN"));
  // ds1.plotOn(frame3, Binning(128));
  //
  // // preparing the signal distribution
  // RooRealVar* rrv_mean_CB = new RooRealVar("rrv_mean_CB", "rrv_mean_CB", 2000, 1000, 3000);
  // RooRealVar* rrv_sigma_CB = new RooRealVar("rrv_sigma_CB", "rrv_sigma_CB", 100, 50, 300);
  // RooRealVar* rrv_tail_CB_I = new RooRealVar("rrv_tail_CB_I", "rrv_tail_CB_I", 2,0., 40);
  // RooRealVar* rrv_tail_CB_II = new RooRealVar("rrv_tail_CB_II", "rrv_tail_CB_II", -2., -40., 0.);
  //
  // RooRealVar* rrv_normalization_CB_I = new RooRealVar("rrv_normalization_CB_I", "rrv_normalization_CB_I", 2, 0., 400);
  // RooRealVar* rrv_frac_CB = new RooRealVar("rrv_frac_CB", "rrv_frac_CB", 0.5);
  //
  //
  // RooCBShape* Crystal_Ball_I = new RooCBShape("CrystalBall_I", "CrystalBall_I",
  //                                             *WR_RecoMass,
  //                                             *rrv_mean_CB,*rrv_sigma_CB,*rrv_tail_CB_I,*rrv_normalization_CB_I);
  //
  // RooCBShape* Crystal_Ball_II = new RooCBShape("CrystalBall_II", "CrystalBall_II",
  //                                             *WR_RecoMass,
  //                                             *rrv_mean_CB,*rrv_sigma_CB,*rrv_tail_CB_II,*rrv_normalization_CB_I);
  //
  // RooAddPdf* WR_pdf = new RooAddPdf("model_pdf", "model_pdf",
  //                                     RooArgList(*Crystal_Ball_I,*Crystal_Ball_II),
  //                                     RooArgList(*rrv_frac_CB));
  //
  // // fit distribution to data
  // RooFitResult *r = WR_pdf->fitTo(ds1, Save());
  // r->Print();
  //
  // std::cout << "The value of rrv_frac_CB is " << rrv_frac_CB->getVal() << std::endl;
  //
  // // Draw ntuples
  // WR_pdf->plotOn(frame1);
  //
  //
  //
  //
  // TCanvas *c = new TCanvas("Test Fit", "Test Fit", 1000, 800);
  // c->Divide(2,2);
  // c->cd(1);
  // frame1->Draw();
  // c->cd(2);
  // frame2->Draw();
  // c->cd(3);
  // frame3->Draw();
}

RooDSCBShape* DSCB_init(RooRealVar* rrv_x, double mean, std::string label)
{
 RooRealVar* rrv_mean_CB = new RooRealVar((std::string("rrv_mean_CB")+label).c_str(), label.c_str(), mean, 0.8*mean, 1.1*mean);
 RooRealVar* rrv_sigma_CB = new RooRealVar((std::string("rrv_sigma_CB")+label).c_str(), label.c_str(), 260, 50, 2000);
 RooRealVar* rrv_alpha_CB_I = new RooRealVar((std::string("rrv_alpha_CB_I")+label).c_str(), label.c_str(), 1, 0., 20);
 RooRealVar* rrv_alpha_CB_II = new RooRealVar((std::string("rrv_alpha_CB_II")+label).c_str(), label.c_str(), -1., -20., 0.);

 RooRealVar* rrv_n_CB_I = new RooRealVar((std::string("rrv_n_CB_I")+label).c_str(), label.c_str(), 300, 0., 3000.);
 RooRealVar* rrv_n_CB_II = new RooRealVar((std::string("rrv_n_CB_II")+label).c_str(), label.c_str(), 300, 0., 3000.);

 return new RooDSCBShape((std::string("DoubleSideCrystallBall")+label).c_str(), label.c_str(), *rrv_x, *rrv_mean_CB,*rrv_sigma_CB,*rrv_alpha_CB_I,*rrv_alpha_CB_II,*rrv_n_CB_I,*rrv_n_CB_II);
}

RooDataSet Hist2Pulls(RooHist* pullPlot, std::string label, bool print=false)
{
  RooRealVar* pullVar = new RooRealVar("pullVar", label.c_str(), -100.0, 100.0);
  RooDataSet pulls("pulls", "pulls", RooArgSet(*pullVar));
  TH1* hist;

  for (Int_t i=0; i<40; i++){
    Double_t binX;
    Double_t pull;
    pullPlot->GetPoint(i, binX, pull);

    if (print==true)
    {
      std::ofstream csv;
      csv.open("pulls.csv", std::ios_base::app);
      csv << pull << "\n";
    }


    RooRealVar pull_i = RooRealVar("pullVar", label.c_str(), pull);
    pulls.add(RooArgSet(pull_i));
  }

  return pulls;
}
