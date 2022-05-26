/**
 * @Author: Billy Li <billyli>
 * @Date:   08-10-2021
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 05-26-2022
 */

// This script is to figure out the best strategy to fit data into a
// Exponential (Exp) CB distribution. It will be put in testFit.C to compare the result
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

void fit_DSCB(std::string filePath)
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
                    "mumujj invm (GeV)",WRGenMean*0.5, WRGenMean*1.4);
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
  RooDSCBShape* WR_ee_DSCB = DSCB_init(eejjMass_WR, WRGenMean, "eejj");
  RooDSCBShape* WR_mumu_DSCB = DSCB_init(mumujjMass_WR, WRGenMean, "mumujj");


  //// fit distribution to data
  RooFitResult *r1 = WR_ee_DSCB->fitTo(ds_WR_eejj, Save(), SumW2Error(kTRUE), Range(WRGenMean*0.45,WRGenMean*1.45));
  RooFitResult *r2 = WR_mumu_DSCB->fitTo(ds_WR_mumujj, Save(), SumW2Error(kTRUE), Range(WRGenMean*0.45,WRGenMean*1.45));

  std::cout << "BELOW IS THE RESULT" << std::endl;
  r1->Print();
  r2->Print();
  std::cout << "ABOVE IS THE RESULTS" << std::endl;

  std::ofstream result1("results_DSCB_ee/"+filePath+".txt");
  std::ofstream result2("results_DSCB_mumu/"+filePath+".txt");
  result1 << "WR:" << WRGenMean << std::endl;
  result1 << "N:" << NGenMean << std::endl;
  result2 << "WR:" << WRGenMean << std::endl;
  result2 << "N:" << NGenMean << std::endl;
  r1->printMultiline(result1, 1, kTRUE, "");
  r2->printMultiline(result2, 1, kTRUE, "");

  //// Prepare frames for plotting
  RooPlot *eeFrame_DSCB = eejjMass_WR->frame(Title(("eejj DSCB "+filePath).c_str()));
  RooPlot *mumuFrame_DSCB = mumujjMass_WR->frame(Title(("mumujj DSCB "+filePath).c_str()));

  //// Plot on frames
  // plot data on frames
  ds_WR_eejj.plotOn(eeFrame_DSCB, Binning(150), DataError(RooAbsData::SumW2));
  ds_WR_mumujj.plotOn(mumuFrame_DSCB, Binning(150), DataError(RooAbsData::SumW2));
  // plot fitted pdfs on frames
  WR_ee_DSCB->plotOn(eeFrame_DSCB);
  WR_mumu_DSCB->plotOn(mumuFrame_DSCB);

  //// pull related
  // Prepare pulls
  RooRealVar* pullVar = new RooRealVar("pullVar", "pull value", -6, 6);
  std::cout << "Making the pull plots" << std::endl;
  RooHist *eeHist_DSCBPull = eeFrame_DSCB->pullHist();
  RooHist *mumuHist_DSCBPull = mumuFrame_DSCB->pullHist();
  // Extract pulls from RooHist
  RooDataSet ee_DSCBPulls = Hist2Pulls(eeHist_DSCBPull,"eejj", true);
  RooDataSet mumu_DSCBPulls = Hist2Pulls(mumuHist_DSCBPull, "mumujj", true);
  // Prepare frame for the pull histograms
  RooPlot* ee_DSCBPullFrame = pullVar->frame(Title("ee DSCB Pull Hist"));
  RooPlot* mumu_DSCBPullFrame = pullVar->frame(Title("mumu DSCB pull Hist"));
  // plot pull histograms on frames
  std::cout << "Making the pull histograms" << std::endl;
  ee_DSCBPulls.plotOn(ee_DSCBPullFrame, Binning(15));
  mumu_DSCBPulls.plotOn(mumu_DSCBPullFrame, Binning(15));


  // chi2
  double chi2_ee_DSCB = eeFrame_DSCB->chiSquare(6);
  double chi2_mumu_DSCB = mumuFrame_DSCB->chiSquare(6);
  std::cout << "chi2_ee_DSCB: " << chi2_ee_DSCB <<std::endl;
  std::cout << "chi2_mumu_DSCB: " << chi2_mumu_DSCB <<std::endl;
  result1 << "chi2_ee_DSCB: " << chi2_ee_DSCB <<std::endl;
  result2 << "chi2_mumu_DSCB: " << chi2_mumu_DSCB <<std::endl;


  //// plot the errors of the fitted functions
  //// This needs to be done after the pulls was calculated
  //// otherwise it will interfer the pull calculations
  // WR_ee_DSCB->plotOn(eeFrame_DSCB, VisualizeError(*r1, 1, kFALSE));
  // WR_mumu_DSCB->plotOn(mumuFrame_DSCB, VisualizeError(*r2, 1, kFALSE));

  //// Draw Frames on TCanvas
  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 1800, 1600);
  c->Divide(2,2);
  c->cd(1);
  eeFrame_DSCB->Draw();
  c->cd(2);
  mumuFrame_DSCB->Draw();
  c->cd(3);
  ee_DSCBPullFrame->Draw();
  c->cd(4);
  mumu_DSCBPullFrame->Draw();
  std::string plot_file_prefix = "plots_DSCB/";
  std::string plotPath = plot_file_prefix+filePath;
  plotPath.erase(plotPath.length()-5); // remove .root
  plotPath = plotPath+ ".png";
  c->SaveAs((plotPath).c_str());
  c->Close();
}

RooDSCBShape* DSCB_init(RooRealVar* rrv_x, double mean, std::string label)
{
 RooRealVar* rrv_mean_CB = new RooRealVar((std::string("rrv_mean_CB_")+label).c_str(), label.c_str(), mean, 0.8*mean, 1.1*mean);
 RooRealVar* rrv_sigma_CB = new RooRealVar((std::string("rrv_sigma_CB_")+label).c_str(), label.c_str(), 260, 50, 2000);
 RooRealVar* rrv_alpha_CB_I = new RooRealVar((std::string("rrv_alpha_CB_I_")+label).c_str(), label.c_str(), 1, 0., 20);
 RooRealVar* rrv_alpha_CB_II = new RooRealVar((std::string("rrv_alpha_CB_II_")+label).c_str(), label.c_str(), -1., -20., 0.);

 RooRealVar* rrv_n_CB_I = new RooRealVar((std::string("rrv_n_CB_I_")+label).c_str(), label.c_str(), 5, 0., 3000.);
 RooRealVar* rrv_n_CB_II = new RooRealVar((std::string("rrv_n_CB_II_")+label).c_str(), label.c_str(), 5, 0., 3000.);

 return new RooDSCBShape((std::string("DoubleSideCrystallBall_")+label).c_str(), label.c_str(), *rrv_x, *rrv_mean_CB,*rrv_sigma_CB,*rrv_alpha_CB_I,*rrv_alpha_CB_II,*rrv_n_CB_I,*rrv_n_CB_II);
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
