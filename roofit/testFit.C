/**
 * @Author: Billy Li <billyli>
 * @Date:   08-10-2021
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 09-20-2021
 */



#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
using namespace RooFit;

RooAddPdf* DoubleCB(RooRealVar* rrv_x, double mean);
double Nll2L(double& Nll);
double geoAvg(double& product, double& dFree);
double Nll2LAvg(double& Nll, double& dFree);
double NEvtInRange(RooDataSet& ds, std::string name, double min, double max);
RooDataSet Hist2Pulls(RooHist* pullPlot);

void testFit(std::string filePath)
{
  //// Extract useful information and set fitting parameters
  // set bin number for pullHist and chi2
  double bin_size = 256;
  // Extract WR and N mean value via the file name
  std::cout << "Openning file " << filePath << std::endl;
  size_t fileNamePos = filePath.find_last_of("/");
  std::string fileName = filePath.substr(fileNamePos+1);
  size_t RPos = fileName.find_last_of("R");
  size_t NPos = fileName.find_last_of("N");
  size_t dotPos = fileName.find_last_of(".");
  double WRGenMean = std::stod(fileName.substr(RPos+1, NPos-RPos));
  double NGenMean = std::stod(fileName.substr(NPos+1, dotPos-NPos));
  std::cout << "Target WR: " WRGenMean << "; Target N" << NGenMean << std::endl;
  // calculate bin number, bin_lo and bin_hi for each bin
  int binNum = (int)WRGenMean*0.6/bin_size;
  std::vector<double> bin_left, bin_right;
  for (int i=0;i<binNum;i++){
    bin_left.push_back(WRGenMean*0.65+i*bin_size);
    bin_right.push_back(WRGenMean*0.65+(i+1)*bin_size);
  }

  //// importing TNtuples from root files
  std::string prefix = "../";
  RooRealVar* WR_RecoMass_ee = new RooRealVar("WR_RecoMass_ee", "WR_RecoMass_ee", 0, 2*WRGenMean);
  RooRealVar* WR_RecoMass_mumu = new RooRealVar("WR_RecoMass_mumu", "WR_RecoMass_mumu", 0, 2*WRGenMean);
  RooDataSet ds_WR_RecoMass_ee("ds1", "ds1",
                RooArgSet(*WR_RecoMass_ee),
                ImportFromFile((prefix+filePath).c_str(), "analysis/WR_RecoMass_ee"));
  RooDataSet ds_WR_RecoMass_mumu("ds2", "ds2",
                RooArgSet(*WR_RecoMass_mumu),
                ImportFromFile((prefix+filePath).c_str(), "analysis/WR_RecoMass_mumu"));

  //// Preparing probability distirbution functions for fitting
  // preparing the double CB distributions
  RooAddPdf* WR_ee_doubleCB = DoubleCB(WR_RecoMass_ee, WRGenMean);
  RooAddPdf* WR_mumu_doubleCB = DoubleCB(WR_RecoMass_mumu, WRGenMean);
  // preparing the single CB distributions
  // ee
  RooRealVar m0_ee("m0_ee","m0 for ee", WRGenMean, 0.8*WRGenMean, 1.1*WRGenMean);
  RooRealVar sigma_ee("sigma_ee","sigma for ee", 200, 50, 2000);
  RooRealVar alpha_ee("alpha_ee", "alpha for ee", 2.0, 0., 200.0);
  RooRealVar n_ee("n_ee","n for ee", 2.0, 0.0, 400.0);
  RooCBShape cb_ee("signal_ee", "cb signal for ee",
                *WR_RecoMass_ee,
                m0_ee, sigma_ee, alpha_ee, n_ee);
  // mumu
  RooRealVar m0_mumu("m0_mumu","m0 for mumu", WRGenMean, 0.8*WRGenMean, 1.1*WRGenMean);
  RooRealVar sigma_mumu("sigma_mumu","sigma for mumu", 200, 50, 2000);
  RooRealVar alpha_mumu("alpha_mumu", "alpha for mumu", 2.0, 0., 200.0);
  RooRealVar n_mumu("n_mumu","n for mumu", 2.0, 0.0, 400.0);
  RooCBShape cb_mumu("signal_mumu", "cb signal for mumu",
                *WR_RecoMass_mumu,
                m0_mumu, sigma_mumu, alpha_mumu, n_mumu);

  //// fit distribution to data
  RooFitResult *r1 = WR_ee_doubleCB->fitTo(ds_WR_RecoMass_ee, Save(), Range(WRGenMean*0.65,WRGenMean*1.25));
  RooFitResult *r2 = WR_mumu_doubleCB->fitTo(ds_WR_RecoMass_mumu, Save(), Range(WRGenMean*0.65,WRGenMean*1.25));
  RooFitResult *r3 = cb_ee.fitTo(ds_WR_RecoMass_ee, Save(), Range(WRGenMean*0.65,WRGenMean*1.25));
  RooFitResult *r4 = cb_mumu.fitTo(ds_WR_RecoMass_mumu, Save(), Range(WRGenMean*0.65,WRGenMean*1.25));

  //// Prepare frames for plotting
  RooPlot *eeFrame_doubleCB = WR_RecoMass_ee->frame(Title("eejj Double CB"));
  RooPlot *eeFrame_CB = WR_RecoMass_ee->frame(Title("eejj CB"));
  RooPlot *mumuFrame_doubleCB = WR_RecoMass_mumu->frame(Title("mumujj Double CB"));
  RooPlot* mumuFrame_CB = WR_RecoMass_mumu->frame(Title("mumujj CB"));


  //// Plot on frames
  // plot data on frames
  ds_WR_RecoMass_ee.plotOn(eeFrame_doubleCB, Binning(256), DataError(RooAbsData::SumW2));
  ds_WR_RecoMass_ee.plotOn(eeFrame_CB, Binning(256), DataError(RooAbsData::SumW2));
  ds_WR_RecoMass_mumu.plotOn(mumuFrame_doubleCB, Binning(256), DataError(RooAbsData::SumW2));
  ds_WR_RecoMass_mumu.plotOn(mumuFrame_CB, Binning(256), DataError(RooAbsData::SumW2));
  // plot fitted pdfs on frames
  WR_ee_doubleCB->plotOn(eeFrame_doubleCB);
  WR_mumu_doubleCB->plotOn(mumuFrame_doubleCB);
  cb_ee.plotOn(eeFrame_CB);
  cb_mumu.plotOn(mumuFrame_CB);

  //// pull related
  // Prepare pulls
  RooRealVar* pullVar = new RooRealVar("pullVar", "pull value", -2, 2);
  RooHist *eeHist_doubleCBPull = eeFrame_doubleCB->pullHist();
  RooHist *mumuHist_doubleCBPull = mumuFrame_doubleCB->pullHist();
  RooHist *eeHist_CBPull = eeFrame_CB->pullHist();
  RooHist *mumuHist_CBPull = mumuFrame_CB->pullHist();
  // Extract pulls from RooHist
  RooDataSet ee2CBPulls = Hist2Pulls(eeHist_doubleCBPull);
  RooDataSet mumu2CBPulls = Hist2Pulls(mumuHist_doubleCBPull);
  RooDataSet eeCBPulls = Hist2Pulls(eeHist_CBPull);
  RooDataSet mumuCBPulls = Hist2Pulls(mumuHist_CBPull);
  // Prepare frame for the pull histograms
  RooPlot* ee2CBPullFrame = pullVar->frame(Title("ee Double CB Pull Hist"));
  RooPlot* mumu2CBPullFrame = pullVar->frame(Title("mumu Double CB pull Hist"));
  RooPlot* eeCBPullFrame = pullVar->frame(Title("ee CB Pull Hist"));
  RooPlot* mumuCBPullFrame = pullVar->frame(Title("mumu CB Pull Hist"));
  // plot pull histograms on frames
  ee2CBPulls.plotOn(ee2CBPullFrame, Binning(32));
  mumu2CBPulls.plotOn(mumu2CBPullFrame, Binning(32));
  eeCBPulls.plotOn(eeCBPullFrame, Binning(32));
  mumuCBPulls.plotOn(mumuCBPullFrame, Binning(32));

  //// calculate and print fit parameters
  // minimum NLL
  double minNll_eeDoubleCB = r1->minNll();
  double minNll_eeCB = r2->minNll();
  double minNll_mumuDoubleCB = r3->minNll();
  double minNll_mumuCB = r4->minNll();
  // number of fitted points and dof
  double dFree_ee2CB = NEvtInRange(ds_WR_RecoMass_ee, "WR_RecoMass_ee", WRGenMean*0.65, WRGenMean*1.25)-5.0;
  double dFree_eeCB = NEvtInRange(ds_WR_RecoMass_ee, "WR_RecoMass_ee", WRGenMean*0.65, WRGenMean*1.25)-4.0;
  double dFree_mumu2CB = NEvtInRange(ds_WR_RecoMass_mumu, "WR_RecoMass_mumu", WRGenMean*0.65, WRGenMean*1.25)-5.0;
  double dFree_mumuCB = NEvtInRange(ds_WR_RecoMass_mumu, "WR_RecoMass_mumu", WRGenMean*0.65, WRGenMean*1.25)-4.0;
  std::cout << dFree_ee2CB << std::endl;
  std::cout << dFree_eeCB << std::endl;
  std::cout << dFree_mumu2CB << std::endl;
  std::cout << dFree_mumuCB << std::endl;
  // average likelihood
  double LAvg_eeDoubleCB = Nll2LAvg(minNll_eeDoubleCB, dFree_ee2CB);
  double LAvg_eeCB = Nll2LAvg(minNll_eeCB, dFree_eeCB);
  double LAvg_mumuDoubleCB = Nll2LAvg(minNll_mumuDoubleCB, dFree_mumu2CB);
  double LAvg_mumuCB = Nll2LAvg(minNll_mumuCB, dFree_mumuCB);
  std::cout << "Average L for eeDoubleCB: "<< LAvg_eeDoubleCB << "\n";
  std::cout << "Average L for eeCB: " << LAvg_eeCB << "\n";
  std::cout << "Average L for mumuDoubleCB: " << LAvg_mumuDoubleCB << "\n";
  std::cout << "Average L for mumuCB: " << LAvg_mumuCB << "\n";
  // chi2
  double chi2_ee2CB = eeFrame_doubleCB->chiSquare();
  double chi2_eeCB = eeFrame_CB->chiSquare();
  double chi2_mumu2CB = mumuFrame_doubleCB->chiSquare();
  double chi2_mumuCB = mumuFrame_CB->chiSquare();
  std::cout << "chi2_ee2CB: " << chi2_ee2CB <<std::endl;
  std::cout << "chi2_eeCB: " << chi2_eeCB <<std::endl;
  std::cout << "chi2_mumu2CB: " << chi2_mumu2CB <<std::endl;
  std::cout << "chi2_mumuCB: " << chi2_mumuCB <<std::endl;

  //// Draw Frames on TCanvas
  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 1000, 800);
  c->Divide(4,2);
  c->cd(1);
  eeFrame_doubleCB->Draw();
  c->cd(2);
  eeFrame_CB->Draw();
  c->cd(3);
  mumuFrame_doubleCB->Draw();
  c->cd(4);
  mumuFrame_CB->Draw();
  c->cd(5);
  ee2CBPullFrame->Draw();
  c->cd(6);
  eeCBPullFrame->Draw();
  c->cd(7);
  mumu2CBPullFrame->Draw();
  c->cd(8);
  mumuCBPullFrame->Draw();
}

RooAddPdf* DoubleCB(RooRealVar* rrv_x, double mean)
{
  RooRealVar* rrv_mean_CB = new RooRealVar("rrv_mean_CB", "rrv_mean_CB", mean, 0.8*mean, 1.1*mean);
  RooRealVar* rrv_sigma_CB = new RooRealVar("rrv_sigma_CB", "rrv_sigma_CB", 200, 50, 2000);
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

// Negative log likelihood to likelihood
double Nll2L(double& Nll)
{
  return std::exp(-Nll);
}

// calculate geometric average of data
double geoAvg(double& product, double& dFree)
{
  return std::pow(product, 1.0/dFree);
}

double Nll2LAvg(double& Nll, double& dFree)
{
  return std::exp(-Nll/dFree);
}

double NEvtInRange(RooDataSet& ds, std::string name, double min, double max)
{
  double num=0;
  Int_t numEntries=ds.numEntries();
  for (Int_t i=0; i<numEntries; i++){
    auto data = ds.get(i)->getRealValue(name.c_str());
    // std::cout << ds.get(i)->contentsString() << " " << data << std::endl;
    if ((data>min)&&(data<max)) num++;
  }
  return num;
}

RooDataSet Hist2Pulls(RooHist* pullPlot)
{
  RooRealVar* pullVar = new RooRealVar("pullVar", "pull variable", -100.0, 100.0);
  RooDataSet pulls("pulls", "pulls", RooArgSet(*pullVar));
  TH1* hist;

  for (Int_t i=0; i<256; i++){
    Double_t binX;
    Double_t pull;
    pullPlot->GetPoint(i, binX, pull);
    RooRealVar pull_i = RooRealVar("pullVar", "pull variable", pull);
    pulls.add(RooArgSet(pull_i));
  }

  return pulls;
}
