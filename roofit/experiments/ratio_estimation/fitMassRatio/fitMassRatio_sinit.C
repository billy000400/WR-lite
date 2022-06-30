/**
 * @Author: Billy Li <billyli>
 * @Date:   08-10-2021
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 06-29-2022
 */

// This script is to figure out the best strategy to fit data into a
// Exponential (Exp) CB distribution. It will be put in testFit.C to compare the result
// of ExpmCB and single CB

#include <sys/stat.h>

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

// void testFit_ExpmCB(std::string filePath)
void fitMassRatio_sinit(std::string sample_file_name, std::string fsig_tr)
{
  //// message service
  // RooFit::RooMsgService::instance().getStream(1).removeTopic(NumericIntegration) ;

  //// parse file name
  std::string MassPairIdxRoot = sample_file_name.substr(sample_file_name.find("_")+1);
  std::string MassPair = MassPairIdxRoot.substr(0,MassPairIdxRoot.find("_"));
  std::string Idx = MassPairIdxRoot.substr(MassPairIdxRoot.find("_")+1, MassPairIdxRoot.find(".")-MassPairIdxRoot.find("_")-1);
  int Idx_val = std::stoi(Idx);

  //// set sample dir
  char data_dir[64] = "../../../data/RooFitMC/ratio_";
  strcat(data_dir, fsig_tr.c_str());
  strcat(data_dir, "/");
  strcat(data_dir, MassPair.c_str());
  strcat(data_dir, "/");

  //// set result dir
  char MCFitResult_dir[128] = "../../../data/MCFitResult/";
  char fitMassRatio_dir[128] = "../../../data/MCFitResult/fitMassRatio_sinit/";
  char result_dir[128] = "../../../data/MCFitResult/fitMassRatio_sinit/";

  strcat(result_dir, MassPair.c_str());
  strcat(result_dir, "/");
  char mass_dir[128];
  strcpy(mass_dir, result_dir);

  strcat(result_dir, "ratio_");
  strcat(result_dir, fsig_tr.c_str());
  strcat(result_dir, "/");

  //// mkdir if not exist
  mkdir(MCFitResult_dir, 0700);
  mkdir(fitMassRatio_dir, 0700);
  mkdir(mass_dir, 0700);
  mkdir(result_dir, 0700);

  //// prepare random generator for fit parameter initializtaion
  TRandom2 *fsig_mumu_gen = new TRandom2(0);
  TRandom2 *fsig_ee_gen = new TRandom2(1);
  TRandom2 *mu_mumu_gen = new TRandom2(2);
  TRandom2 *mu_ee_gen = new TRandom2(3);

  Double_t mu_low = 700;
  Double_t mu_high = 2500;


  //// init WR distributio
  // mumujj WR
  double alpha_mm_val = 1.4021e+00;
  double alpha_mm_err = 0;
  double beta_mm_val = 4.0786e-01;
  double beta_mm_err = 0;
  double m_mm_val = 1.1626e+00;
  double m_mm_err = 0;
  double n_mm_val = 1.6055e+00;
  double n_mm_err = 0;
  double sigma_mm_val = 8.0404e+01;
  double sigma_mm_err = 0;

  Double_t mu_mumu_init = mu_mumu_gen->Rndm()*(mu_high-mu_low)+mu_low;

  RooRealVar* mumujjMass = new RooRealVar("invm_mumujj", "invm reco from mumujj", 700, 2500);
  RooRealVar* mu_mm= new RooRealVar("mu", "mu mumujj", mu_mumu_init, mu_low, mu_high);
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
  double n_ee_val = 1.6276e+00;
  double n_ee_err = 0;
  double sigma_ee_val = 6.7117e+01;
  double sigma_ee_err = 0;

  Double_t mu_ee_init = mu_ee_gen->Rndm()*(mu_high-mu_low)+mu_low;

  RooRealVar* eejjMass = new RooRealVar("invm_eejj", "invm reco from eejj", 700, 2500);
  RooRealVar* mu_ee= new RooRealVar("mu", "mu eejj", mu_ee_init, mu_low, mu_high);
  RooRealVar* sigma_ee = new RooRealVar("sigma", "sigma eejj", sigma_ee_val, sigma_ee_val-sigma_ee_err, sigma_ee_val+sigma_ee_err);
  RooRealVar* alpha_ee = new RooRealVar("alpha", "alpha eejj", alpha_ee_val, alpha_ee_val-alpha_ee_err, alpha_ee_val+alpha_ee_err);
  RooRealVar* n_ee = new RooRealVar("n", "n eejj", n_ee_val, n_ee_val-n_ee_err, n_ee_val+n_ee_err);
  RooRealVar* beta_ee = new RooRealVar("beta", "beta", beta_ee_val, beta_ee_val-beta_ee_err, beta_ee_val+beta_ee_err);
  RooRealVar* m_ee = new RooRealVar("m", "m eejj", m_ee_val, m_ee_val-m_ee_err, m_ee_val+m_ee_err);

  RooExpmCB* WR_eejj = new RooExpmCB("WR eejj", "WR eejj",\
          *eejjMass, *mu_ee,*sigma_ee,*beta_ee,*m_ee,*alpha_ee,*n_ee);

  //// init bg distribution
  // mumu
  double a_mm_val =  -3.5652e-03;
  double a_mm_err = 0;
  double b_mm_val =  1;
  double b_mm_err = 0;

  RooRealVar* a_mm = new RooRealVar("a_mm", "a_mm", a_mm_val, a_mm_val-a_mm_err, a_mm_val+a_mm_err);
  RooRealVar* b_mm = new RooRealVar("b_mm", "b_mm", b_mm_val, b_mm_val-b_mm_err, b_mm_val+b_mm_err);

  RooExpm* bg_mumujj = new RooExpm("bg mumujj", "Expm Bg", *mumujjMass, *a_mm, *b_mm);
  // ee
  double a_ee_val = -3.5191e-03;
  double a_ee_err = 0;
  double b_ee_val =  1;
  double b_ee_err = 0;

  RooRealVar* a_ee = new RooRealVar("a_ee", "a_ee", a_ee_val, a_ee_val-a_ee_err, a_ee_val+a_ee_err);
  RooRealVar* b_ee = new RooRealVar("b_ee", "b_ee", b_ee_val, b_ee_val-b_ee_err, b_ee_val+b_ee_err);

  RooExpm* bg_eejj = new RooExpm("bg eejj", "Expm Bg", *eejjMass, *a_ee, *b_ee);


  // add distribution
  Double_t fsig_low = 0;
  Double_t fsig_high = 1e-1;
  Double_t fsig_mumu_init = fsig_mumu_gen->Rndm()*(fsig_high-fsig_low)+fsig_low;
  Double_t fsig_ee_init = fsig_ee_gen->Rndm()*(fsig_high-fsig_low)+fsig_low;
  RooRealVar *fsig_mumu = new RooRealVar("fsig_mumu", "signal fraction mumujj", fsig_mumu_init, fsig_low, fsig_high);
  RooRealVar *fsig_ee = new RooRealVar("fsig_ee", "signal fraction eejj", fsig_ee_init, fsig_low , fsig_high);

  RooAddPdf *model_ee = new RooAddPdf("composite_ee", "model ee", RooArgList(*WR_eejj, *bg_eejj), *fsig_ee);
  RooAddPdf *model_mumu = new RooAddPdf("composite_mumu", "model mumu", RooArgList(*WR_mumujj, *bg_mumujj), *fsig_mumu);

  // prepare TFile
  char sample_file_path[128];
  strcpy(sample_file_path, data_dir);
  strcat(sample_file_path, sample_file_name.c_str());

  // file = TFile::Open(fname.Data()); if(!file||file->IsZombie()){delete file; continue;}
  TFile sample_file = TFile(sample_file_path, "READ");
  if (!sample_file.IsOpen() || sample_file.IsZombie()){
      std::cout << "Zombie file " << sample_file_path << std::endl;
      return;
  }

  RooDataSet ds_mumujj("ds_mumujj", "ds_mumujj",\
                RooArgSet(*mumujjMass),\
                ImportFromFile(sample_file_path, "composite_mumuData"));

  RooDataSet ds_eejj("ds_eejj", "ds_eejj",\
                RooArgSet(*eejjMass),\
                ImportFromFile(sample_file_path, "composite_eeData"));

  char result_file_path[128];
  strcpy(result_file_path, result_dir);
  char result_file_name[64] = "MCFitResult_";
  strcat(result_file_name, MassPairIdxRoot.c_str());
  strcat(result_file_path, result_file_name);
  TFile fResult_file(result_file_path,"RECREATE");
  TTree tree("fit_result","mass ratio and their uncertainties");

  RooFitResult *r_mumu = model_mumu->fitTo(ds_mumujj, Save(), SumW2Error(kTRUE), Range(700,2500), Strategy(2), Minos(kTRUE));
  RooFitResult *r_ee= model_ee->fitTo(ds_eejj, Save(), SumW2Error(kTRUE), Range(700,2500), Strategy(2), Minos(kTRUE));

  double fsig_mumu_val = fsig_mumu->getVal();
  double fsig_mumu_hi = fsig_mumu->getAsymErrorHi();
  double fsig_mumu_lo = fsig_mumu->getAsymErrorLo();

  double fsig_ee_val = fsig_ee->getVal();
  double fsig_ee_hi = fsig_ee->getAsymErrorHi();
  double fsig_ee_lo = fsig_ee->getAsymErrorLo();

  double mu_mumu_val = mu_mm->getVal();
  double mu_mumu_hi = mu_mm->getAsymErrorHi();
  double mu_mumu_lo = mu_mm->getAsymErrorLo();

  double mu_ee_val = mu_ee->getVal();
  double mu_ee_hi = mu_ee->getAsymErrorHi();
  double mu_ee_lo = mu_ee->getAsymErrorLo();

  tree.Branch("fsig_mumu_val", &fsig_mumu_val);
  tree.Branch("fsig_mumu_hi", &fsig_mumu_hi);
  tree.Branch("fsig_mumu_lo", &fsig_mumu_lo);

  tree.Branch("fsig_ee_val", &fsig_ee_val);
  tree.Branch("fsig_ee_hi", &fsig_ee_hi);
  tree.Branch("fsig_ee_lo", &fsig_ee_lo);

  tree.Branch("mu_mumu_val", &mu_mumu_val);
  tree.Branch("mu_mumu_hi", &mu_mumu_hi);
  tree.Branch("mu_mumu_lo", &mu_mumu_lo);

  tree.Branch("mu_ee_val", &mu_ee_val);
  tree.Branch("mu_ee_hi", &mu_ee_hi);
  tree.Branch("mu_ee_lo", &mu_ee_lo);

  tree.Fill();

  fResult_file.Write();
  fResult_file.Close();
}
