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
#include "RooExpmCB.h"
using namespace RooFit;

RooExpmCB* ExpmCB_init(RooRealVar* rrv_x, std::string label);
RooExpm* Expm_init(RooRealVar* x, std::string label);
RooDataSet Hist2Pulls(RooHist* pullPlot, std::string label, bool print=false);

// void testFit_ExpmCB(std::string filePath)
void testFit_WR_bg()
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

  RooRealVar* mumujjMass_DY = new RooRealVar("invm_mumujj",
                            "invm reco from DY mumujj", 500, 3000);
  RooRealVar* mumujjRowWeight_DY = new RooRealVar("rowWeight",
                            "row weight for DY ntuple mumujj rows", -1.5, 1.5);

  RooRealVar* mumujjMass_ttbar = new RooRealVar("invm_mumujj",
                            "invm reco from ttbar mumujj", 500, 3000);
  RooRealVar* mumujjRowWeight_ttbar = new RooRealVar("rowWeight",
                            "row weight for ttbar ntuple mumujj rows", -1.5, 1.5);

  RooRealVar* mumujjMass_all = new RooRealVar("invm_mumujj",
                            "invm reco from all mumujj", 500, 3000);
  RooRealVar* mumujjRowWeight_all = new RooRealVar("rowWeight",
                            "row weight for all ntuple mumujj rows", -1.5, 1.5);

  RooRealVar* eejjMass_WR = new RooRealVar("invm_eejj",
                              "invm reco from WR eejj", 500, 3000);
  RooRealVar* eejjRowWeight_WR = new RooRealVar("rowWeight",
                              "row weight for WR ntuple eejj rows", -1.5, 1.5);

  RooRealVar* eejjMass_DY = new RooRealVar("invm_eejj",
                              "invm reco from DY eejj", 500, 3000);
  RooRealVar* eejjRowWeight_DY = new RooRealVar("rowWeight",
                              "row weight for DY ntuple eejj rows", -1.5, 1.5);

  RooRealVar* eejjMass_ttbar = new RooRealVar("invm_eejj",
                              "invm reco from ttbar eejj", 500, 3000);
  RooRealVar* eejjRowWeight_ttbar = new RooRealVar("rowWeight",
                              "row weight for ttbar ntuple eejj rows", -1.5, 1.5);

  RooRealVar* eejjMass_all = new RooRealVar("invm_eejj",
                              "invm reco from all eejj", 500, 3000);
  RooRealVar* eejjRowWeight_all = new RooRealVar("rowWeight",
                              "row weight for all ntuple eejj rows", -1.5, 1.5);

  //// load dataset
  std::string prefix = "../analysis/allEvents/";

  RooDataSet ds_WR_mumujj("ds_WR_mumujj", "ds_WR_mumujj",
                RooArgSet(*mumujjMass_WR, *mumujjRowWeight_WR),
                ImportFromFile((prefix+"fullWR1600N00.root").c_str(), "invm_mumujj"),
                WeightVar(*mumujjRowWeight_WR));

  RooDataSet ds_DY_mumujj("ds_DY_mumujj", "ds_DY_mumujj",
                RooArgSet(*mumujjMass_DY, *mumujjRowWeight_DY),
                ImportFromFile((prefix+"fullDY.root").c_str(), "invm_mumujj"),
                WeightVar(*mumujjRowWeight_DY));

  RooDataSet ds_ttbar_mumujj("ds_ttbar_mumujj", "ds_ttbar_mumujj",
                RooArgSet(*mumujjMass_ttbar, *mumujjRowWeight_ttbar),
                ImportFromFile((prefix+"fullttbar.root").c_str(), "invm_mumujj"),
                WeightVar(*mumujjRowWeight_ttbar));

  RooDataSet ds_all_mumujj("ds_all_mumujj", "ds_all_mumujj",
                RooArgSet(*mumujjMass_all, *mumujjRowWeight_all),
                WeightVar(*mumujjRowWeight_all));

  ds_all_mumujj.append(ds_WR_mumujj);
  ds_all_mumujj.append(ds_DY_mumujj);
  ds_all_mumujj.append(ds_ttbar_mumujj);

  RooDataSet ds_WR_eejj("ds_WR_eejj", "ds_WR_eejj",
                RooArgSet(*eejjMass_WR, *eejjRowWeight_WR),
                ImportFromFile((prefix+"fullWR1600N1000.root").c_str(), "invm_eejj"),
                WeightVar(*eejjRowWeight_WR));

  RooDataSet ds_DY_eejj("ds_DY_eejj", "ds_DY_eejj",
                RooArgSet(*eejjMass_DY, *eejjRowWeight_DY),
                ImportFromFile((prefix+"fullDY.root").c_str(), "invm_eejj"),
                WeightVar(*eejjRowWeight_DY));

  RooDataSet ds_ttbar_eejj("ds_ttbar_eejj", "ds_ttbar_eejj",
                RooArgSet(*eejjMass_ttbar, *eejjRowWeight_ttbar),
                ImportFromFile((prefix+"fullttbar.root").c_str(), "invm_eejj"),
                WeightVar(*eejjRowWeight_ttbar));

  RooDataSet ds_all_eejj("ds_all_eejj", "ds_all_eejj",
                RooArgSet(*eejjMass_all, *eejjRowWeight_all),
                WeightVar(*eejjRowWeight_all));

  ds_all_eejj.append(ds_WR_eejj);
  ds_all_eejj.append(ds_DY_eejj);
  ds_all_eejj.append(ds_ttbar_eejj);


  //// Preparing probability distirbution functions for fitting
  // preparing the signal distributions
  RooExpmCB* WR_ee_ExpmCB = ExpmCB_init(eejjMass_all, "eejj");
  RooExpmCB* WR_mumu_ExpmCB = ExpmCB_init(mumujjMass_all, "mumujj");
  // preparing the background model
  RooExpm* bg_ee = Expm_init(eejjMass_all, "eejj");
  RooExpm* bg_mumu = Expm_init(mumujjMass_all, "mumujj");
  // add model
  RooRealVar *fsig_ee = new RooRealVar("fsig_ee", "signal fraction eejj", 0.5, 0., 1.);
  RooRealVar *fsig_mumu = new RooRealVar("fsig_mumu", "signal fraction mumujj", 0.5, 0., 1.);
  RooAddPdf *model_ee = new RooAddPdf("model ee", "model ee", RooArgList(*WR_ee_ExpmCB, *bg_ee), *fsig_ee);
  RooAddPdf *model_mumu = new RooAddPdf("model mumu", "model mumu", RooArgList(*WR_mumu_ExpmCB, *bg_mumu), *fsig_mumu);


  //// fit distribution to data
  RooFitResult *r1 = model_ee->fitTo(ds_all_eejj, Save(), SumW2Error(kTRUE), Range(1000,2500));
  RooFitResult *r2 = model_mumu->fitTo(ds_all_mumujj, Save(), SumW2Error(kTRUE), Range(1000,2500));

  std::cout << "BELOW IS THE RESULT" << std::endl;
  r1->Print();
  r2->Print();
  std::cout << "ABOVE IS THE RESULTS" << std::endl;

  // Prepare frames for plotting
  RooPlot *eeFrame = eejjMass_all->frame(Title("eejj"));
  RooPlot *mumuFrame  = mumujjMass_all->frame(Title("mumujj"));

  //// Plot on frames
  // plot data on frames
  ds_all_eejj.plotOn(eeFrame, Binning(100), DataError(RooAbsData::SumW2), Range(1000,2500));
  ds_all_mumujj.plotOn(mumuFrame, Binning(100), DataError(RooAbsData::SumW2), Range(1000,2500));
  // plot fitted pdfs on frames
  model_ee->plotOn(eeFrame);
  model_mumu->plotOn(mumuFrame);

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
  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 600, 300);
  c->Divide(2,1);
  c->cd(1);
  eeFrame->Draw();
  c->cd(2);
  mumuFrame->Draw();
  // c->cd(3);
  // ee_ExpmCBPullFrame->Draw();
  // c->cd(4);
  // mumu_ExpmCBPullFrame->Draw();
}

RooExpmCB* ExpmCB_init(RooRealVar* rrv_x, std::string label)
{
 RooRealVar* rrv_mean_CB = new RooRealVar((std::string("rrv_mean_ExpmCB_")+label).c_str(), label.c_str(), 1400., 600., 2200.);
 RooRealVar* rrv_sigma_CB = new RooRealVar((std::string("rrv_sigma_ExpmCB_")+label).c_str(), label.c_str(), 85);
 RooRealVar* rrv_alpha_CB = new RooRealVar((std::string("rrv_alpha_ExpmCB_")+label).c_str(), label.c_str(), 1.448);
 RooRealVar* rrv_beta_CB = new RooRealVar((std::string("rrv_beta_ExpmCB_")+label).c_str(), label.c_str(), 3.7e-1);
 RooRealVar* rrv_m_CB = new RooRealVar((std::string("rrv_m_ExpmCB_")+label).c_str(), label.c_str(), 1.18);
 RooRealVar* rrv_n_CB = new RooRealVar((std::string("rrv_n_ExpmCB_")+label).c_str(), label.c_str(), 1.9);

 return new RooExpmCB((std::string("Exp(-omega*t^m)CrystallBall_")+label).c_str(), label.c_str(), *rrv_x, *rrv_mean_CB,*rrv_sigma_CB,*rrv_beta_CB,*rrv_m_CB,*rrv_alpha_CB,*rrv_n_CB);
}

RooExpm* Expm_init(RooRealVar* x, std::string label)
{
    RooRealVar* a = new RooRealVar((std::string("a_")+label).c_str(), label.c_str(), -9.4550e-2);
    RooRealVar* b = new RooRealVar((std::string("b_")+label).c_str(), label.c_str(), 6.0258e-1);

    return new RooExpm((std::string("Expm_")+label).c_str(), label.c_str(), *x, *a, *b);
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
