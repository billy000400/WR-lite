/**
 * @Author: Billy Li <billyli>
 * @Date:   08-10-2021
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 05-05-2022
 */

// This script is to figure out the best strategy to fit data into a
// Exponential (Exp) Gauss distribution. It will be put in testFit.C to compare the result
// of ExpmGauss and single Gauss

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "RooExpmGauss.h"
using namespace RooFit;

RooExpmGauss* ExpmGauss_init(RooRealVar* rrv_x, double mean, std::string label);
RooDataSet Hist2Pulls(RooHist* pullPlot, std::string label, bool print=false);

void testFit_ExpmGauss(std::string filePath)
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
  RooRealVar* WR_RecoMass_ee = new RooRealVar("WR_RecoMass_ee", "WR_RecoMass_ee", WRGenMean*0.5, WRGenMean*1.45);
  RooRealVar* WR_RecoMass_mumu = new RooRealVar("WR_RecoMass_mumu", "WR_RecoMass_mumu", WRGenMean*0.5, WRGenMean*1.45);
  RooDataSet ds_WR_RecoMass_ee("ds1", "ds1",
                RooArgSet(*WR_RecoMass_ee),
                ImportFromFile((prefix+filePath).c_str(), "analysis/WR_RecoMass_ee"));
  RooDataSet ds_WR_RecoMass_mumu("ds2", "ds2",
                RooArgSet(*WR_RecoMass_mumu),
                ImportFromFile((prefix+filePath).c_str(), "analysis/WR_RecoMass_mumu"));

  //// Preparing probability distirbution functions for fitting
  // preparing the double Gauss distributions
  RooExpmGauss* WR_ee_doubleGauss = ExpmGauss_init(WR_RecoMass_ee, WRGenMean, "eejj");
  RooExpmGauss* WR_mumu_doubleGauss = ExpmGauss_init(WR_RecoMass_mumu, WRGenMean, "mumujj");


  //// fit distribution to data
  RooFitResult *r1 = WR_ee_doubleGauss->fitTo(ds_WR_RecoMass_ee, Save(), Range(WRGenMean*0.45,WRGenMean*1.45));
  RooFitResult *r2 = WR_mumu_doubleGauss->fitTo(ds_WR_RecoMass_mumu, Save(), Range(WRGenMean*0.45,WRGenMean*1.45));

  //// Prepare frames for plotting
  RooPlot *eeFrame_doubleGauss = WR_RecoMass_ee->frame(Title("eejj ExpmGauss"));
  RooPlot *mumuFrame_doubleGauss = WR_RecoMass_mumu->frame(Title("mumujj ExpmGauss"));

  //// Plot on frames
  // plot data on frames
  ds_WR_RecoMass_ee.plotOn(eeFrame_doubleGauss, Binning(300), DataError(RooAbsData::SumW2));
  ds_WR_RecoMass_mumu.plotOn(mumuFrame_doubleGauss, Binning(300), DataError(RooAbsData::SumW2));
  // plot fitted pdfs on frames
  WR_ee_doubleGauss->plotOn(eeFrame_doubleGauss);
  WR_mumu_doubleGauss->plotOn(mumuFrame_doubleGauss);

  //// pull related
  // Prepare pulls
  RooRealVar* pullVar = new RooRealVar("pullVar", "pull value", -6, 6);
  std::cout << "Making the pull plots" << std::endl;
  RooHist *eeHist_doubleGaussPull = eeFrame_doubleGauss->pullHist();
  RooHist *mumuHist_doubleGaussPull = mumuFrame_doubleGauss->pullHist();
  // Extract pulls from RooHist
  RooDataSet ee2GaussPulls = Hist2Pulls(eeHist_doubleGaussPull,"eejj", true);
  RooDataSet mumu2GaussPulls = Hist2Pulls(mumuHist_doubleGaussPull, "mumujj", true);
  // Prepare frame for the pull histograms
  RooPlot* ee2GaussPullFrame = pullVar->frame(Title("ee ExpmGauss Pull Hist"));
  RooPlot* mumu2GaussPullFrame = pullVar->frame(Title("mumu ExpmGauss pull Hist"));
  // plot pull histograms on frames
  std::cout << "Making the pull histograms" << std::endl;
  ee2GaussPulls.plotOn(ee2GaussPullFrame, Binning(15));
  mumu2GaussPulls.plotOn(mumu2GaussPullFrame, Binning(15));

  //// plot the errors of the fitted functions
  //// This needs to be done after the pulls was calculated
  //// otherwise it will interfer the pull calculations
  WR_ee_doubleGauss->plotOn(eeFrame_doubleGauss, VisualizeError(*r1, 1, kFALSE));
  WR_mumu_doubleGauss->plotOn(mumuFrame_doubleGauss, VisualizeError(*r2, 1, kFALSE));

  //// Draw Frames on TCanvas
  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 600, 600);
  c->Divide(2,2);
  c->cd(1);
  eeFrame_doubleGauss->Draw();
  c->cd(2);
  mumuFrame_doubleGauss->Draw();
  c->cd(3);
  ee2GaussPullFrame->Draw();
  c->cd(4);
  mumu2GaussPullFrame->Draw();

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
  // RooRealVar* rrv_mean_Gauss = new RooRealVar("rrv_mean_Gauss", "rrv_mean_Gauss", 2000, 1000, 3000);
  // RooRealVar* rrv_sigma_Gauss = new RooRealVar("rrv_sigma_Gauss", "rrv_sigma_Gauss", 100, 50, 300);
  // RooRealVar* rrv_tail_Gauss_I = new RooRealVar("rrv_tail_Gauss_I", "rrv_tail_Gauss_I", 2,0., 40);
  // RooRealVar* rrv_tail_Gauss_II = new RooRealVar("rrv_tail_Gauss_II", "rrv_tail_Gauss_II", -2., -40., 0.);
  //
  // RooRealVar* rrv_normalization_Gauss_I = new RooRealVar("rrv_normalization_Gauss_I", "rrv_normalization_Gauss_I", 2, 0., 400);
  // RooRealVar* rrv_frac_Gauss = new RooRealVar("rrv_frac_Gauss", "rrv_frac_Gauss", 0.5);
  //
  //
  // RooGaussShape* Crystal_Ball_I = new RooGaussShape("CrystalBall_I", "CrystalBall_I",
  //                                             *WR_RecoMass,
  //                                             *rrv_mean_Gauss,*rrv_sigma_Gauss,*rrv_tail_Gauss_I,*rrv_normalization_Gauss_I);
  //
  // RooGaussShape* Crystal_Ball_II = new RooGaussShape("CrystalBall_II", "CrystalBall_II",
  //                                             *WR_RecoMass,
  //                                             *rrv_mean_Gauss,*rrv_sigma_Gauss,*rrv_tail_Gauss_II,*rrv_normalization_Gauss_I);
  //
  // RooAddPdf* WR_pdf = new RooAddPdf("model_pdf", "model_pdf",
  //                                     RooArgList(*Crystal_Ball_I,*Crystal_Ball_II),
  //                                     RooArgList(*rrv_frac_Gauss));
  //
  // // fit distribution to data
  // RooFitResult *r = WR_pdf->fitTo(ds1, Save());
  // r->Print();
  //
  // std::cout << "The value of rrv_frac_Gauss is " << rrv_frac_Gauss->getVal() << std::endl;
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

RooExpmGauss* ExpmGauss_init(RooRealVar* rrv_x, double mean, std::string label)
{
 RooRealVar* rrv_mean_Gauss = new RooRealVar((std::string("rrv_mean_ExpmGauss_")+label).c_str(), label.c_str(), mean, 0.8*mean, 1.1*mean);
 RooRealVar* rrv_sigma_Gauss = new RooRealVar((std::string("rrv_sigma_ExpmGauss_")+label).c_str(), label.c_str(), 0.15*mean, 0.01*mean, 0.5*mean);
 RooRealVar* rrv_beta_Gauss = new RooRealVar((std::string("rrv_beta_ExpmGauss_")+label).c_str(), label.c_str(), 1., 1e-4, 5.);
 RooRealVar* rrv_m_Gauss = new RooRealVar((std::string("rrv_m_ExpmGauss_")+label).c_str(), label.c_str(), 1., 1e-4, 5.);


 return new RooExpmGauss((std::string("Exp(-omega*t^m)CrystallBall_")+label).c_str(), label.c_str(), *rrv_x, *rrv_mean_Gauss,*rrv_sigma_Gauss,*rrv_beta_Gauss,*rrv_m_Gauss);
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
