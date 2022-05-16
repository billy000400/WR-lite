/**
 * @Author: Billy Li <billyli>
 * @Date:   08-10-2021
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 05-15-2022
 */

// This script is to figure out the best strategy to fit data into a
// Exponential (Exp) CB distribution. It will be put in testFit.C to compare the result
// of ExpCB and single CB

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "RooExpCBShape.h"
using namespace RooFit;

RooExpCBShape* ExpCB_init(RooRealVar* rrv_x, double mean, std::string label);
RooDataSet Hist2Pulls(RooHist* pullPlot, std::string label, bool print=false);

void fit_ExpCB(std::string filePath)
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


  std::string prefix = "../analysis/allEvents/";
  RooRealVar* eejjMass_WR = new RooRealVar("invm_eejj",\
                    "invm reco from eejj",WRGenMean*0.5, WRGenMean*1.4);
  RooRealVar* eejjRowWeight_WR = new RooRealVar("rowWeight",\
                    "row weight for WR ntuple eejj rows", -1.5, 1.5);

  RooRealVar* mumujjMass_WR = new RooRealVar("invm_mumujj",\
                    "invm reco from mumujj",WRGenMean*0.5, WRGenMean*1.4);
  RooRealVar* mumujjRowWeight_WR = new RooRealVar("rowWeight",\
                    "row weight for WR ntuple mumujj rows", -1.5, 1.5);

  RooDataSet ds_WR_eejj("ds_WR_eejj", "ds_WR_eejj",
                RooArgSet(*eejjMass_WR, *eejjRowWeight_WR),
                ImportFromFile((prefix+filePath).c_str(), "invm_eejj"),
                WeightVar(*eejjRowWeight_WR));
  RooDataSet ds_WR_mumujj("ds_WR_mumujj", "ds_WR_mumujj",
                RooArgSet(*mumujjMass_WR, *mumujjRowWeight_WR),
                ImportFromFile((prefix+filePath).c_str(), "invm_mumujj"),
                WeightVar(*mumujjRowWeight_WR));

  //// Preparing probability distirbution functions for fitting
  // preparing the double CB distributions
  RooExpCBShape* WR_ee_ExpCB = ExpCB_init(eejjMass_WR, WRGenMean, "eejj");
  RooExpCBShape* WR_mumu_ExpCB = ExpCB_init(mumujjMass_WR, WRGenMean, "mumujj");


  //// fit distribution to data
  RooFitResult *r1 = WR_ee_ExpCB->fitTo(ds_WR_eejj, Save(), SumW2Error(kTRUE), Range(WRGenMean*0.45,WRGenMean*1.45));
  RooFitResult *r2 = WR_mumu_ExpCB->fitTo(ds_WR_mumujj, Save(), SumW2Error(kTRUE), Range(WRGenMean*0.45,WRGenMean*1.45));

  std::cout << "BELOW IS THE RESULT" << std::endl;
  r1->Print();
  r2->Print();
  std::cout << "ABOVE IS THE RESULTS" << std::endl;

  //// Prepare frames for plotting
  RooPlot *eeFrame_ExpCB = eejjMass_WR->frame(Title(("eejj ExpCB "+filePath).c_str()));
  RooPlot *mumuFrame_ExpCB = mumujjMass_WR->frame(Title(("mumujj ExpCB "+filePath).c_str()));

  //// Plot on frames
  // plot data on frames
  ds_WR_eejj.plotOn(eeFrame_ExpCB, Binning(150), DataError(RooAbsData::SumW2));
  ds_WR_mumujj.plotOn(mumuFrame_ExpCB, Binning(150), DataError(RooAbsData::SumW2));
  // plot fitted pdfs on frames
  WR_ee_ExpCB->plotOn(eeFrame_ExpCB);
  WR_mumu_ExpCB->plotOn(mumuFrame_ExpCB);

  //// pull related
  // Prepare pulls
  RooRealVar* pullVar = new RooRealVar("pullVar", "pull value", -6, 6);
  std::cout << "Making the pull plots" << std::endl;
  RooHist *eeHist_ExpCBPull = eeFrame_ExpCB->pullHist();
  RooHist *mumuHist_ExpCBPull = mumuFrame_ExpCB->pullHist();
  // Extract pulls from RooHist
  RooDataSet ee_ExpCBPulls = Hist2Pulls(eeHist_ExpCBPull,"eejj", true);
  RooDataSet mumu_ExpCBPulls = Hist2Pulls(mumuHist_ExpCBPull, "mumujj", true);
  // Prepare frame for the pull histograms
  RooPlot* ee_ExpCBPullFrame = pullVar->frame(Title("ee ExpCB Pull Hist"));
  RooPlot* mumu_ExpCBPullFrame = pullVar->frame(Title("mumu ExpCB pull Hist"));
  // plot pull histograms on frames
  std::cout << "Making the pull histograms" << std::endl;
  ee_ExpCBPulls.plotOn(ee_ExpCBPullFrame, Binning(15));
  mumu_ExpCBPulls.plotOn(mumu_ExpCBPullFrame, Binning(15));


  // chi2
  double chi2_ee_ExpCB = eeFrame_ExpCB->chiSquare(6);
  double chi2_mumu_ExpCB = mumuFrame_ExpCB->chiSquare(6);
  std::cout << "chi2_ee_ExpCB: " << chi2_ee_ExpCB <<std::endl;
  std::cout << "chi2_mumu_ExpCB: " << chi2_mumu_ExpCB <<std::endl;


  //// plot the errors of the fitted functions
  //// This needs to be done after the pulls was calculated
  //// otherwise it will interfer the pull calculations
  // WR_ee_ExpCB->plotOn(eeFrame_ExpCB, VisualizeError(*r1, 1, kFALSE));
  // WR_mumu_ExpCB->plotOn(mumuFrame_ExpCB, VisualizeError(*r2, 1, kFALSE));

  //// Draw Frames on TCanvas
  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 600, 600);
  c->Divide(2,2);
  c->cd(1);
  eeFrame_ExpCB->Draw();
  c->cd(2);
  mumuFrame_ExpCB->Draw();
  c->cd(3);
  ee_ExpCBPullFrame->Draw();
  c->cd(4);
  mumu_ExpCBPullFrame->Draw();

  TCanvas *d = new TCanvas("Test Fit", "Test Fit", 1000, 1000);
  mumuFrame_DSCB->Draw();
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

RooExpCBShape* ExpCB_init(RooRealVar* rrv_x, double mean, std::string label)
{
 RooRealVar* rrv_mean_CB = new RooRealVar((std::string("rrv_mean_ExpCB_")+label).c_str(), label.c_str(), mean, 0.8*mean, 1.1*mean);
 RooRealVar* rrv_sigma_CB = new RooRealVar((std::string("rrv_sigma_ExpCB_")+label).c_str(), label.c_str(), 0.15*mean, 0.01*mean, 0.5*mean);
 RooRealVar* rrv_alpha_CB = new RooRealVar((std::string("rrv_alpha_ExpCB_")+label).c_str(), label.c_str(), 1, 0., 5.0);
 RooRealVar* rrv_n_CB = new RooRealVar((std::string("rrv_n_ExpCB_")+label).c_str(), label.c_str(), 5, 0., 100.);
 RooRealVar* rrv_beta_CB = new RooRealVar((std::string("rrv_beta_ExpCB_")+label).c_str(), label.c_str(), 1., 0., 5.);


 return new RooExpCBShape((std::string("ExponentialCrystallBall_")+label).c_str(), label.c_str(), *rrv_x, *rrv_mean_CB,*rrv_sigma_CB,*rrv_alpha_CB,*rrv_n_CB,*rrv_beta_CB);
}

RooDataSet Hist2Pulls(RooHist* pullPlot, std::string label, bool print=false)
{
  RooRealVar* pullVar = new RooRealVar("pullVar", label.c_str(), -100.0, 100.0);
  RooDataSet pulls("pulls", "pulls", RooArgSet(*pullVar));
  TH1* hist;

  for (Int_t i=0; i<150; i++){
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
