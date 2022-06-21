/**
 * @Author: Billy Li <billyli>
 * @Date:   08-10-2021
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 06-20-2022
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
void fitRatio_WR1600N800()
{
  //// message service
  // RooFit::RooMsgService::instance().getStream(1).removeTopic(NumericIntegration) ;

  //// set sample number
  int sampleNum = 2;
  // int mumujjEventNum = 515750;
  // int eejjEventNum = 727838;

  //// set data dir
  char prefix[64] = "../../../data/ratio_1e-3/WR1600N800/"; // will concatenate with sample name

  TFile fResult_file("test_1e-3.root","RECREATE");
  TTree tree("fit_result","WR1600 N800 ratio");

  double fsig_mumu_val = -1;
  double fsig_ee_val = -1;
  tree.Branch("fsig_mumu", &fsig_mumu_val);
  tree.Branch("fsig_ee", &fsig_ee_val);
  for (int i=0; i<sampleNum; i++){
    std::cout << "Fitting sample: " << i+1 << "/" << sampleNum << std::endl;

    //// init WR distribution
    // mumujj WR
    double alpha_mm_val = 1.4021e+00;
    double alpha_mm_err = 0;
    double beta_mm_val = 4.0786e-01;
    double beta_mm_err = 0;
    double m_mm_val = 1.1626e+00;
    double m_mm_err = 0;
    double mu_mm_val = 1.5965e+03;
    double mu_mm_err = 0;
    double n_mm_val = 1.6055e+00;
    double n_mm_err = 0;
    double sigma_mm_val = 8.0404e+01;
    double sigma_mm_err = 0;


    RooRealVar* mumujjMass = new RooRealVar("invm_mumujj", "invm reco from mumujj", 400, 3000);
    RooRealVar* mu_mm= new RooRealVar("mu", "mu mumujj", mu_mm_val, mu_mm_val-mu_mm_err, mu_mm_val+mu_mm_err);
    RooRealVar* sigma_mm = new RooRealVar("sigma", "sigma mumujj", sigma_mm_val, sigma_mm_val-sigma_mm_err, sigma_mm_val+sigma_mm_err);
    RooRealVar* alpha_mm = new RooRealVar("alpha", "alpha mumujj", alpha_mm_val, alpha_mm_val-alpha_mm_err, alpha_mm_val+alpha_mm_err);
    RooRealVar* n_mm = new RooRealVar("n", "n mumujj", n_mm_val, n_mm_val-n_mm_err, n_mm_val+n_mm_err);
    RooRealVar* beta_mm = new RooRealVar("beta", "beta", beta_mm_val, beta_mm_val-beta_mm_err, beta_mm_val+beta_mm_err);
    RooRealVar* m_mm = new RooRealVar("m", "m mumujj", m_mm_val, m_mm_val-m_mm_err, m_mm_val+m_mm_err);

    RooExpmCB* WR_mumujj = new RooExpmCB("WR mumujj", "WR mumujj",\
            *mumujjMass, *mu_mm,*sigma_mm,*beta_mm,*m_mm,*alpha_mm,*n_mm);

    // eejj WR
    double alpha_ee_val = 1.3094;
    double alpha_ee_err = 0;
    double beta_ee_val = 2.6121e-01;
    double beta_ee_err = 0;
    double m_ee_val = 1.2494;
    double m_ee_err = 0;
    double mu_ee_val = 1.6053e+03;
    double mu_ee_err = 0;
    double n_ee_val = 1.6276e+00;
    double n_ee_err = 0;
    double sigma_ee_val = 6.7117e+01;
    double sigma_ee_err = 0;

    RooRealVar* eejjMass = new RooRealVar("invm_eejj", "invm reco from eejj", 400, 3000);
    RooRealVar* mu_ee= new RooRealVar("mu", "mu eejj", mu_ee_val, mu_ee_val-mu_ee_err, mu_ee_val+mu_ee_err);
    RooRealVar* sigma_ee = new RooRealVar("sigma", "sigma eejj", sigma_ee_val, sigma_ee_val-sigma_ee_err, sigma_ee_val+sigma_ee_err);
    RooRealVar* alpha_ee = new RooRealVar("alpha", "alpha eejj", alpha_ee_val, alpha_ee_val-alpha_ee_err, alpha_ee_val+alpha_ee_err);
    RooRealVar* n_ee = new RooRealVar("n", "n eejj", n_ee_val, n_ee_val-n_ee_err, n_ee_val+n_ee_err);
    RooRealVar* beta_ee = new RooRealVar("beta", "beta", beta_ee_val, beta_ee_val-beta_ee_err, beta_ee_val+beta_ee_err);
    RooRealVar* m_ee = new RooRealVar("m", "m eejj", m_ee_val, m_ee_val-m_ee_err, m_ee_val+m_ee_err);

    RooExpmCB* WR_eejj = new RooExpmCB("WR eejj", "WR eejj",\
            *eejjMass, *mu_ee,*sigma_ee,*beta_ee,*m_ee,*alpha_ee,*n_ee);

    //// init bg distribution
    // mumu
    double a_mm_val = -1.7249e-01;
    double a_mm_err = 0;
    double b_mm_val =  5.3266e-01;
    double b_mm_err = 0;

    RooRealVar* a_mm = new RooRealVar("a_mm", "a_mm", a_mm_val, a_mm_val-a_mm_err, a_mm_val+a_mm_err);
    RooRealVar* b_mm = new RooRealVar("b_mm", "b_mm", b_mm_val, b_mm_val-b_mm_err, b_mm_val+b_mm_err);

    RooExpm* bg_mumujj = new RooExpm("bg mumujj", "Expm Bg", *mumujjMass, *a_mm, *b_mm);
    // ee
    double a_ee_val = -1.7834e-01;
    double a_ee_err = 0;
    double b_ee_val =  5.3223e-01;
    double b_ee_err = 0;

    RooRealVar* a_ee = new RooRealVar("a_ee", "a_ee", a_ee_val, a_ee_val-a_ee_err, a_ee_val+a_ee_err);
    RooRealVar* b_ee = new RooRealVar("b_ee", "b_ee", b_ee_val, b_ee_val-b_ee_err, b_ee_val+b_ee_err);

    RooExpm* bg_eejj = new RooExpm("bg eejj", "Expm Bg", *eejjMass, *a_ee, *b_ee);


    // add distribution
    RooRealVar *fsig_mumu = new RooRealVar("fsig_mumu", "signal fraction mumujj", 7e-4, 5e-4, 15e-4);
    RooRealVar *fsig_ee = new RooRealVar("fsig_ee", "signal fraction eejj", 7e-4, 5e-4, 15e-4);

    RooAddPdf *model_ee = new RooAddPdf("composite_ee", "model ee", RooArgList(*WR_eejj, *bg_eejj), *fsig_ee);
    RooAddPdf *model_mumu = new RooAddPdf("composite_mumu", "model mumu", RooArgList(*WR_mumujj, *bg_mumujj), *fsig_mumu);

    // prepare TFile
    char sample_file_name[32] = "RooFitMC_WR1600N800_";
    char sample_index_str[32];
    sprintf(sample_index_str, "%d", i+1);
    strcat(sample_file_name, sample_index_str);
    char sample_file_path[32];
    strcpy(sample_file_path, prefix);
    strcat(sample_file_path, sample_file_name);
    strcat(sample_file_path, ".root");

    RooDataSet ds_mumujj("ds_mumujj", "ds_mumujj",\
                  RooArgSet(*mumujjMass),\
                  ImportFromFile(sample_file_path, "composite_mumuData"));

    RooDataSet ds_eejj("ds_eejj", "ds_eejj",\
                  RooArgSet(*eejjMass),\
                  ImportFromFile(sample_file_path, "composite_eeData"));

    RooFitResult *r_mumu = model_mumu->fitTo(ds_mumujj, Save(), SumW2Error(kTRUE), Range(800,2000));
    RooFitResult *r_ee= model_ee->fitTo(ds_eejj, Save(), SumW2Error(kTRUE), Range(800,2000));

    fsig_mumu_val = fsig_mumu->getVal();
    fsig_ee_val = fsig_ee->getVal();
    // f.cd();
    // r_mumu->Write("mumu", TObject::kSingleKey);
    // r_ee->Write("ee", TObject::kSingleKey);
    tree.Fill();
  }
  fResult_file.Write();
  fResult_file.Close();

  // Prepare frames for plotting
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
  // c->cd(3);
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
