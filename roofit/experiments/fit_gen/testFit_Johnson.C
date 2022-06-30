/**
 * @Author: Billy Li <billyli>
 * @Date:   08-10-2021
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 06-29-2022
 */

// A sample fit using simple initialization and trivial numerical stability
// #include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TApplication.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooFitResult.h"

using namespace RooFit;

RooDataSet Hist2Pulls(RooHist* pullPlot, std::string label, bool print);

void testFit_Johnson(std::string filePath)
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


  std::string prefix = "../../../analysis/allEvents/";
  RooRealVar* WR_GenMass = new RooRealVar("WR_GenMass",\
                    "gen mass extracted from simulation",WRGenMean*0.45, WRGenMean*1.55);

  RooDataSet ds_WR_GenMass("ds_WR_GenMass", "ds_WR_GenMass",
                RooArgSet(*WR_GenMass),
                ImportFromFile((prefix+filePath).c_str(), "WR_GenMass"));

  //// Prepare frames for plotting
  RooPlot *frame= WR_GenMass->frame(Title(("RooJohnson "+filePath).c_str()));

  // Plot on frames
  // plot data on frames
  ds_WR_GenMass.plotOn(frame, Binning(50));

  //// Preparing probability distirbution functions for fitting
  // preparing the Johnson distributions
  RooRealVar *mu = new RooRealVar("mu", "mu", WRGenMean, 0.8*WRGenMean, 1.2*WRGenMean);
  RooRealVar *lm = new RooRealVar("lambda", "lambda", 50, 10, 500);
  RooRealVar *gm = new RooRealVar("gamma", "gamma", 5, 0.1, 20);
  RooRealVar *dt = new RooRealVar("delta", "delta", 2, 0.01, 20);
  double massThreshold = NGenMean;
  RooJohnson* model = new RooJohnson("RooJohnson", "RooJohnson", *WR_GenMass,\
                *mu, *lm, *gm, *dt, massThreshold);


  //// fit distribution to data
  RooFitResult *r = model->fitTo(ds_WR_GenMass, Save(kTRUE),\
                                  Range(WRGenMean*0.45,WRGenMean*1.55), Offset(kTRUE),\
                                  Strategy(2), Strategy(2));

  model->plotOn(frame);

  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 1000, 1000);
  c->Divide(1,1);
  c->cd(1);
  frame->Draw();

  //// compare which initialization gives the lowest chi2
  // calculate chi2
  // RooPlot *eeFrame = eejjMass_WR->frame(Title(("eejj ExpmCB "+filePath).c_str()));
  // RooPlot *mumuFrame = mumujjMass_WR->frame(Title(("mumujj ExpmCB "+filePath).c_str()));
  //
  // ds_WR_eejj.plotOn(eeFrame, Binning(150), DataError(RooAbsData::SumW2));
  // ds_WR_mumujj.plotOn(mumuFrame, Binning(150), DataError(RooAbsData::SumW2));
  //
  // WR_ee_ExpmCB->plotOn(eeFrame);
  // WR_mumu_ExpmCB->plotOn(mumuFrame);
  //
  // double chi2_ee_ExpmCB = eeFrame->chiSquare(6);
  // double chi2_mumu_ExpmCB = mumuFrame->chiSquare(6);
  //
  // // save parameters
  // std::cout << "BELOW IS THE RESULT" << std::endl;
  // r->Print();
  // std::cout << "ABOVE IS THE RESULTS" << std::endl;
  //
  // std::cout << "BELOW IS THE CORRELATION MATRIX" << std::endl;
  // r->correlationMatrix().Print();
  // std::cout << "ABOVE IS THE CORRELATION MATRIX" << std::endl;
  //
  // std::ofstream result1("results_simple_ee/"+filePath+".txt");
  // std::ofstream result2("results_simple_mumu/"+filePath+".txt");
  // result1 << "WR:" << WRGenMean << std::endl;
  // result1 << "N:" << NGenMean << std::endl;
  // result2 << "WR:" << WRGenMean << std::endl;
  // result2 << "N:" << NGenMean << std::endl;
  // r_ee->printMultiline(result1, 1, kTRUE, "");
  // r_mumu->printMultiline(result2, 1, kTRUE, "");

  // // save chi2
  // std::cout << "chi2_ee_ExpmCB: " << chi2_ee_ExpmCB <<std::endl;
  // std::cout << "chi2_mumu_ExpmCB: " << chi2_mumu_ExpmCB <<std::endl;
  // result1 << "chi2_ee_ExpmCB: " << chi2_ee_ExpmCB <<std::endl;
  // result2 << "chi2_mumu_ExpmCB: " << chi2_mumu_ExpmCB <<std::endl;
  //
  // // plot picked pdfs on frames
  // WR_ee_ExpmCB->plotOn(eeFrame_ExpmCB);
  // WR_mumu_ExpmCB->plotOn(mumuFrame_ExpmCB);
  //
  // //// pull related
  // // Prepare pulls
  // RooRealVar* pullVar = new RooRealVar("pullVar", "pull value", -6, 6);
  // std::cout << "Making the pull plots" << std::endl;
  // RooHist *eeHist_ExpmCBPull = eeFrame_ExpmCB->pullHist();
  // RooHist *mumuHist_ExpmCBPull = mumuFrame_ExpmCB->pullHist();
  // // Extract pulls from RooHist
  // RooDataSet ee_ExpmCBPulls = Hist2Pulls(eeHist_ExpmCBPull,"eejj", false);
  // RooDataSet mumu_ExpmCBPulls = Hist2Pulls(mumuHist_ExpmCBPull, "mumujj", false);
  // // Prepare frame for the pull histograms
  // RooPlot* ee_ExpmCBPullFrame = pullVar->frame(Title("ee ExpmCB Pull Hist"));
  // RooPlot* mumu_ExpmCBPullFrame = pullVar->frame(Title("mumu ExpmCB pull Hist"));
  // // plot pull histograms on frames
  // std::cout << "Making the pull histograms" << std::endl;
  // ee_ExpmCBPulls.plotOn(ee_ExpmCBPullFrame, Binning(15));
  // mumu_ExpmCBPulls.plotOn(mumu_ExpmCBPullFrame, Binning(15));

  //// plot the errors of the fitted functions
  // This needs to be done after the pulls was calculated
  // otherwise it will interfer the pull calculations
  // WR_ee_ExpmCB->plotOn(eeFrame_ExpmCB, VisualizeError(*r1, 1, kFALSE));
  // WR_mumu_ExpmCB->plotOn(mumuFrame_ExpmCB, VisualizeError(*r2, 1, kFALSE));

  //// Draw Frames on TCanvas
  // TCanvas *c = new TCanvas("Test Fit", "Test Fit", 1800, 1600);
  // c->Divide(2,2);
  // c->cd(1);
  // eeFrame_ExpmCB->Draw();
  // c->cd(2);
  // mumuFrame_ExpmCB->Draw();
  // c->cd(3);
  // ee_ExpmCBPullFrame->Draw();
  // c->cd(4);
  // mumu_ExpmCBPullFrame->Draw();
  // std::string plot_file_prefix = "plots_simple/";
  // std::string plotPath = plot_file_prefix+filePath;
  // plotPath.erase(plotPath.length()-5); // remove .root
  // plotPath = plotPath+ ".png";
  // c->SaveAs((plotPath).c_str());
  // c->Close();
}

// RooExpmCB* ExpmCB_init(RooRealVar* rrv_x, double mean, std::string label)
// {
//  RooRealVar* rrv_mean_CB = new RooRealVar((std::string("rrv_mean_ExpmCB_")+label).c_str(), label.c_str(), mean, 0.8*mean, 1.1*mean);
//  RooRealVar* rrv_sigma_CB = new RooRealVar((std::string("rrv_sigma_ExpmCB_")+label).c_str(), label.c_str(), 0.15*mean, 0.01*mean, 0.5*mean);
//  RooRealVar* rrv_alpha_CB = new RooRealVar((std::string("rrv_alpha_ExpmCB_")+label).c_str(), label.c_str(), 0.5, 1e-1, 10.0);
//  RooRealVar* rrv_n_CB = new RooRealVar((std::string("rrv_n_ExpmCB_")+label).c_str(), label.c_str(), 1, 0., 5.);
//  RooRealVar* rrv_beta_CB = new RooRealVar((std::string("rrv_beta_ExpmCB_")+label).c_str(), label.c_str(), 0.2, 1e-1, 3.);
//  RooRealVar* rrv_m_CB = new RooRealVar((std::string("rrv_m_ExpmCB_")+label).c_str(), label.c_str(), 0.35, 1e-2, 2.);
//
//
//  return new RooExpmCB((std::string("Exp(-omega*t^m)CrystallBall_")+label).c_str(), label.c_str(), *rrv_x, *rrv_mean_CB,*rrv_sigma_CB,*rrv_beta_CB,*rrv_m_CB,*rrv_alpha_CB,*rrv_n_CB);
// }


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
