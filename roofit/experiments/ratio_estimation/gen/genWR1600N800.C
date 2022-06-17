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

void genWR1600N800()
{
  // set data dir
  std::string prefix = "../../../data/WR1600N800/";

  double mean = 1600.0;

  // init WR distribution
  RooRealVar* rrv_mean_CB = new RooRealVar((std::string("mu")).c_str(), "mu", mean, 0.8*mean, 1.1*mean);
  RooRealVar* rrv_sigma_CB = new RooRealVar((std::string("sigma")).c_str(), "sigma", 0.05*mean, 0.01*mean, 0.1*mean);
  RooRealVar* rrv_alpha_CB = new RooRealVar((std::string("alpha")).c_str(), "alpha", 2, 0.1, 10.0);
  RooRealVar* rrv_n_CB = new RooRealVar((std::string("n")).c_str(), "n", 1, 0.5, 2);
  RooRealVar* rrv_beta_CB = new RooRealVar((std::string("beta")).c_str(), "beta", 0.5, 0.01, 3.);
  RooRealVar* rrv_m_CB = new RooRealVar((std::string("m")).c_str(), "m", 1.5, 1e-2, 2.);

  RooExpmCB* WR = new RooExpmCB((std::string("Expm_CB")).c_str(), "WR", *rrv_x, *rrv_mean_CB,*rrv_sigma_CB,*rrv_beta_CB,*rrv_m_CB,*rrv_alpha_CB,*rrv_n_CB);

  // init bg distribution

  // add distribution

  // generate sample
  return
}
