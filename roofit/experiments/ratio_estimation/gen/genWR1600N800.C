#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include <sys/stat.h>
using namespace RooFit;

void genWR1600N800()
{
  //// set data dir
  char sample_file_path[32] = "../../../data/WR1600N800/"; // will concatenate with sample name

  //// init WR distribution
  // mumujj WR
  double alpha_mm_val = 1.4021e+00;
  double alpha_mm_err = 5.36e-02;
  double beta_mm_val = 4.0786e-01;
  double beta_mm_err = 4.31e-02;
  double m_mm_val = 1.1626e+00;
  double m_mm_err = 3.15e-02;
  double mu_mm_val = 1.5965e+03;
  double mu_mm_err = 2.42e+00;
  double n_mm_val = 1.6055e+00;
  double n_mm_err = 1.14e-01;
  double sigma_mm_val = 8.0404e+01;
  double sigma_mm_err = 2.71e+00;


  RooRealVar* mumujjMass = new RooRealVar("invm_mumujj", "invm reco from mumujj", 800, 2000);
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
  double alpha_ee_err =3.56e-02;
  double beta_ee_val = 2.6121e-01;
  double beta_ee_err = 3.08e-02;
  double m_ee_val = 1.2494;
  double m_ee_err = 2.88e-02;
  double mu_ee_val = 1.6053e+03;
  double mu_ee_err = 2.49e+00;
  double n_ee_val = 1.6276e+00;
  double n_ee_err = 8.19e-02;
  double sigma_ee_val = 6.7117e+01;
  double sigma_ee_err = 2.85e+00;

  RooRealVar* eejjMass = new RooRealVar("invm_eejj", "invm reco from eejj", 800, 2000);
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
  double a_mm_val = -1.4000e-01;
  double a_mm_err = 1.76e-02;
  double b_mm_val =  5.5761e-01;
  double b_mm_err = 3.22e-01;

  RooRealVar* a_mm = new RooRealVar("a_mm", "a_mm", a_mm_val, a_mm_val-a_mm_err, a_mm_val+a_mm_err);
  RooRealVar* b_mm = new RooRealVar("b_mm", "b_mm", b_mm_val, b_mm_val-b_mm_err, b_mm_val+b_mm_err);

  RooExpm* bg_mumujj = new RooExpm("bg mumujj", "Expm Bg", *mumujjMass, *a_mm, *b_mm);
  // ee
  double a_ee_val = -1.3995e-01;
  double a_ee_err = 6.32e-02;
  double b_ee_val =  5.6057e-01;
  double b_ee_err = 3.65e-01;

  RooRealVar* a_ee = new RooRealVar("a_ee", "a_ee", a_ee_val, a_ee_val-a_ee_err, a_ee_val+a_ee_err);
  RooRealVar* b_ee = new RooRealVar("b_ee", "b_ee", b_ee_val, b_ee_val-b_ee_err, b_ee_val+b_ee_err);

  RooExpm* bg_eejj = new RooExpm("bg eejj", "Expm Bg", *eejjMass, *a_ee, *b_ee);


  // add distribution
  RooRealVar *fsig_mumu = new RooRealVar("fsig_mumu", "signal fraction mumujj", 5e-1);
  RooRealVar *fsig_ee = new RooRealVar("fsig_ee", "signal fraction eejj", 5e-1);

  RooAddPdf *model_ee = new RooAddPdf("model ee", "model ee", RooArgList(*WR_eejj, *bg_eejj), *fsig_ee);
  RooAddPdf *model_mumu = new RooAddPdf("model mumu", "model mumu", RooArgList(*WR_mumujj, *bg_mumujj), *fsig_mumu);

  for (int i=0; i<3; i++){
    char sample_file_name[32] = "RooFitMC_WR1600N800_";
    char sample_index_str[32];
    sprintf(sample_index_str, "%d", i);
    strcat(sample_file_name, sample_index_str);
    RooDataSet ds_new("RooFitMC", "RooFit MC WR1600 N800", RooArgSet(*eejjMass, *mumujjMass));
    RooDataSet* ds_new_ee = model_ee->generate(RooArgSet(*eejjMass), Name("ee_tmp"), NumEvents(100));
    RooDataSet* ds_new_mumu = model_mumu->generate(RooArgSet(*mumujjMass), Name("mumu_tmp"), NumEvents(100));
    ds_new_ee->merge(ds_new_mumu);
    ds_new.append(*ds_new_ee);

    strcat(sample_file_path, sample_file_name);
    std::cout << sample_file_path << std::endl;
    TFile outputFile(sample_file_path, "RECREATE");
    ds_new.convertToTreeStore();
  }

  // generate sample

  return;
}
