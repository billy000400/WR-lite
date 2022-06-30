/**
 * @Author: Billy Li <billyli>
 * @Date:   08-10-2021
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 05-16-2022
 */

// This script is to figure out the best strategy to fit data into a
// Exponential (Exp) CB distribution. It will be put in testFit.C to compare the result
// of ExpmCB and single CB

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

RooDataSet Hist2Pulls(RooHist* pullPlot, std::string label, bool print=false);

// void testFit_ExpmCB(std::string filePath)
void testFit_JxCB()
{
  // // Extract useful information and set fitting parameters
  // Extract WR and N mean value via the file name
  // std::cout << "Openning file " << filePath << std::endl;
  // size_t fileNamePos = filePath.find_last_of("/");
  // std::string fileName = filePath.substr(fileNamePos+1);
  // size_t RPos = fileName.find_last_of("R");
  // size_t NPos = fileName.find_last_of("N");
  // size_t dotPos = fileName.find_last_of(".");
  // double WRGenMean = std::stod(fileName.substr(RPos+1, NPos-RPos));
  // double NGenMean = std::stod(fileName.substr(NPos+1, dotPos-NPos));
  // std::cout << "Target WR: " << WRGenMean << ", Target N" << NGenMean << std::endl;

  // Preparing RooRealVars
  RooRealVar* mumujjMass_WR = new RooRealVar("invm_mumujj",
                            "invm reco from WR mumujj", 500, 3000);
  RooRealVar* mumujjRowWeight_WR = new RooRealVar("rowWeight",
                            "row weight for WR ntuple mumujj rows", -1.5, 1.5);


  RooRealVar* eejjMass_WR = new RooRealVar("invm_eejj",
                              "invm reco from WR eejj", 500, 3000);
  RooRealVar* eejjRowWeight_WR = new RooRealVar("rowWeight",
                              "row weight for WR ntuple eejj rows", -1.5, 1.5);

  //// load dataset
  std::string prefix = "../../../analysis/allEvents/";

  RooDataSet ds_WR_mumujj("ds_WR_mumujj", "ds_WR_mumujj",
                RooArgSet(*mumujjMass_WR, *mumujjRowWeight_WR),
                ImportFromFile((prefix+"fullWR2000N1900.root").c_str(), "invm_mumujj"),
                WeightVar(*mumujjRowWeight_WR));

  RooDataSet ds_WR_eejj("ds_WR_eejj", "ds_WR_eejj",
                RooArgSet(*eejjMass_WR, *eejjRowWeight_WR),
                ImportFromFile((prefix+"fullWR2000N1900.root").c_str(), "invm_eejj"),
                WeightVar(*eejjRowWeight_WR));


  //// Preparing probability distirbution functions for fitting
  // preparing the gen distributions
  RooRealVar *mu = new RooRealVar("mu", "mu", 1.99438e+03);
  RooRealVar *lm = new RooRealVar("lambda", "lambda", 1.79013e+01);
  RooRealVar *gm = new RooRealVar("gamma", "gamma",  -5.54826e-01);
  RooRealVar *dt = new RooRealVar("delta", "delta", 6.43057e-01);
  double massThreshold = 1900;
  RooJohnson* gen_mm = new RooJohnson("gen_mm", "RooJohnson mm", *mumujjMass_WR,\
                *mu, *lm, *gm, *dt, massThreshold);
  RooJohnson* gen_ee = new RooJohnson("gen_ee", "RooJohnson ee", *eejjMass_WR,\
                *mu, *lm, *gm, *dt, massThreshold);
  // preparing the resolution model
  RooRealVar *m0 = new RooRealVar("m0", "m0 for CB res", 0.0);
  RooRealVar *sigma_mm = new RooRealVar("sigma_mm", "sigma mumu", 100.0, 5.0, 400.0);
  RooRealVar *sigma_ee = new RooRealVar("sigma_ee", "sigma ee", 100.0, 5.0, 400.0);
  RooRealVar *alpha_mm = new RooRealVar();
  RooRealVar *n_mm = new RooRealVar();
  RooGaussian* res_mm = new RooGaussian("res_mm", "Gaussian mm", *mumujjMass_WR, *mu_g, *sigma_mm);
  RooGaussian* res_ee = new RooGaussian("res_ee", "Gaussian ee", *eejjMass_WR, *mu_g, *sigma_ee);
  // conv model
  mumujjMass_WR->setBins(10000,"fft");
  RooFFTConvPdf conv_mm("conv_mm", "conv mm", *mumujjMass_WR, *gen_mm, *res_mm);
  eejjMass_WR->setBins(10000,"fft");
  RooFFTConvPdf conv_ee("conv_ee", "conv ee", *eejjMass_WR, *gen_ee, *res_ee);
  //// fit distribution to data
  RooFitResult *r_mm = conv_mm.fitTo(ds_WR_mumujj, Save(), SumW2Error(kTRUE), Strategy(2), Range(1000,3000));
  RooFitResult *r_ee = conv_ee.fitTo(ds_WR_eejj, Save(), SumW2Error(kTRUE), Strategy(2), Range(1000,3000));
  //// test plot
  RooPlot *frame_mm = mumujjMass_WR->frame("Johnson x Gaussian mumu");
  RooPlot *frame_ee = eejjMass_WR->frame("Johnson x Gaussian ee");
  ds_WR_mumujj.plotOn(frame_mm, Binning(100));
  ds_WR_eejj.plotOn(frame_ee, Binning(100));
  conv_mm.plotOn(frame_mm, Binning(100));
  conv_ee.plotOn(frame_ee, Binning(100));

  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 2000, 1000);
  c->Divide(2,1);
  c->cd(1);
  frame_mm->Draw();
  c->cd(2);
  frame_ee->Draw();
  //
  // std::cout << "BELOW IS THE RESULT" << std::endl;
  // r1->Print();
  // r2->Print();
  // std::cout << "ABOVE IS THE RESULTS" << std::endl;
  //
  // // Prepare frames for plotting
  // RooPlot *eeFrame = eejjMass_all->frame(Title("eejj"));
  // RooPlot *mumuFrame  = mumujjMass_all->frame(Title("mumujj WR+bg"));
  //
  // //// Plot on frames
  // // plot data on frames
  // ds_all_eejj.plotOn(eeFrame, Binning(100), DataError(RooAbsData::SumW2));
  // ds_all_mumujj.plotOn(mumuFrame, Binning(100), DataError(RooAbsData::SumW2));
  // // plot fitted pdfs on frames
  // model_ee->plotOn(eeFrame);
  // model_mumu->plotOn(mumuFrame);

  // //// pull related
  // // Prepare pulls
  // RooRealVar* pullVar = new RooRealVar("pullVar", "pull value", -6, 6);
  // std::cout << "Making the pull plots" << std::endl;
  // RooHist *eeHist_ExpmCBPull = eeFrame_ExpmCB->pullHist();
  // RooHist *mumuHist_ExpmCBPull = mumuFrame_ExpmCB->pullHist();
  // // Extract pulls from RooHist
  // RooDataSet ee_ExpmCBPulls = Hist2Pulls(eeHist_ExpmCBPull,"eejj", true);
  // RooDataSet mumu_ExpmCBPulls = Hist2Pulls(mumuHist_ExpmCBPull, "mumujj", true);
  // // Prepare frame for the pull histograms
  // RooPlot* ee_ExpmCBPullFrame = pullVar->frame(Title("ee ExpmCB Pull Hist"));
  // RooPlot* mumu_ExpmCBPullFrame = pullVar->frame(Title("mumu ExpmCB pull Hist"));
  // // plot pull histograms on frames
  // std::cout << "Making the pull histograms" << std::endl;
  // ee_ExpmCBPulls.plotOn(ee_ExpmCBPullFrame, Binning(15));
  // mumu_ExpmCBPulls.plotOn(mumu_ExpmCBPullFrame, Binning(15));
  //
  //
  // // chi2
  // double chi2_ee_ExpmCB = eeFrame_ExpmCB->chiSquare(6);
  // double chi2_mumu_ExpmCB = mumuFrame_ExpmCB->chiSquare(6);
  // std::cout << "chi2_ee_ExpmCB: " << chi2_ee_ExpmCB <<std::endl;
  // std::cout << "chi2_mumu_ExpmCB: " << chi2_mumu_ExpmCB <<std::endl;
  //
  //
  // //// plot the errors of the fitted functions
  // //// This needs to be done after the pulls was calculated
  // //// otherwise it will interfer the pull calculations
  // WR_ee_ExpmCB->plotOn(eeFrame_ExpmCB, VisualizeError(*r1, 1, kFALSE));
  // WR_mumu_ExpmCB->plotOn(mumuFrame_ExpmCB, VisualizeError(*r2, 1, kFALSE));
  //
  //// Draw Frames on TCanvas
  // TCanvas *c = new TCanvas("Test Fit", "Test Fit", 600, 300);
  // c->Divide(2,1);
  // c->cd(1);
  // eeFrame->Draw();
  // c->cd(2);
  // mumuFrame->Draw();
  //
  // TCanvas *d = new TCanvas("composite", "composite", 1300, 1000);
  // mumuFrame->Draw();
  // // c->cd(3);
  // ee_ExpmCBPullFrame->Draw();
  // c->cd(4);
  // mumu_ExpmCBPullFrame->Draw();
}

RooDataSet Hist2Pulls(RooHist* pullPlot, std::string label, bool print=false)
{
  RooRealVar* pullVar = new RooRealVar("pullVar", label.c_str(), -100.0, 100.0);
  RooDataSet pulls("pulls", "pulls", RooArgSet(*pullVar));
  TH1* hist;

  for (Int_t i=0; i<100; i++){
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
