/**
 * @Author: Billy Li <billyli>
 * @Date:   05-03-2022
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 06-23-2022
 */



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


void testFit_bg_Expm()
{
  //// Preparing RooRealVars
  // mumujj
  RooRealVar* mumujjMass_DY = new RooRealVar("invm_mumujj",\
                            "invm reco from DY mumujj", 400, 3000);
  RooRealVar* mumujjRowWeight_DY = new RooRealVar("rowWeight",\
                            "row weight for DY ntuple mumujj rows", -1.5, 1.5);

  RooRealVar* mumujjMass_ttbar = new RooRealVar("invm_mumujj",\
                            "invm reco from ttbar mumujj", 400, 3000);
  RooRealVar* mumujjRowWeight_ttbar = new RooRealVar("rowWeight",\
                            "row weight for ttbar ntuple mumujj rows", -1.5, 1.5);

  RooRealVar* mumujjMass_bg = new RooRealVar("invm_mumujj",\
                            "invm reco from bg mumujj (GeV)", 400, 3000);
  RooRealVar* mumujjRowWeight_bg = new RooRealVar("rowWeight",\
                            "row weight for bg ntuple mumujj rows", -1.5, 1.5);

  // eejj
  RooRealVar* eejjMass_DY = new RooRealVar("invm_eejj",\
                            "invm reco from DY eejj", 400, 3000);
  RooRealVar* eejjRowWeight_DY = new RooRealVar("rowWeight",\
                            "row weight for DY ntuple eejj rows", -1.5, 1.5);

  RooRealVar* eejjMass_ttbar = new RooRealVar("invm_eejj",\
                            "invm reco from ttbar eejj", 400, 3000);
  RooRealVar* eejjRowWeight_ttbar = new RooRealVar("rowWeight",\
                            "row weight for ttbar ntuple eejj rows", -1.5, 1.5);

  RooRealVar* eejjMass_bg = new RooRealVar("invm_eejj",\
                            "invm reco from bg eejj (GeV)", 400, 3000);
  RooRealVar* eejjRowWeight_bg = new RooRealVar("rowWeight",\
                            "row weight for bg ntuple eejj rows", -1.5, 1.5);

  // RooRealVar* eejjMass_DY = new RooRealVar();
  // RooRealVar* eejjRowWeight_DY = new RooRealVar("rowWeight", "rowWeight", -1.5, 1.5);

  //// importing ntuples into RooDataSet
  std::string prefix = "../analysis/allEvents/";
  // mumujj
  RooDataSet ds_DY_mumujj("ds_DY_mumujj", "ds_DY_mumujj",
                RooArgSet(*mumujjMass_DY, *mumujjRowWeight_DY),
                ImportFromFile((prefix+"fullDY.root").c_str(), "invm_mumujj"),
                WeightVar(*mumujjRowWeight_DY));

  RooDataSet ds_ttbar_mumujj("ds_ttbar_mumujj", "ds_ttbar_mumujj",
                RooArgSet(*mumujjMass_ttbar, *mumujjRowWeight_ttbar),
                ImportFromFile((prefix+"fullttbar.root").c_str(), "invm_mumujj"),
                WeightVar(*mumujjRowWeight_ttbar));

  RooDataSet ds_bg_mumujj("ds_bg_mumujj", "ds_bg_mumujj",
                RooArgSet(*mumujjMass_bg, *mumujjRowWeight_bg),
                WeightVar(*mumujjRowWeight_bg));

  ds_bg_mumujj.append(ds_DY_mumujj);
  ds_bg_mumujj.append(ds_ttbar_mumujj);
  // eejj
  RooDataSet ds_DY_eejj("ds_DY_eejj", "ds_DY_eejj",
                RooArgSet(*eejjMass_DY, *eejjRowWeight_DY),
                ImportFromFile((prefix+"fullDY.root").c_str(), "invm_eejj"),
                WeightVar(*eejjRowWeight_DY));

  RooDataSet ds_ttbar_eejj("ds_ttbar_eejj", "ds_ttbar_eejj",
                RooArgSet(*eejjMass_ttbar, *eejjRowWeight_ttbar),
                ImportFromFile((prefix+"fullttbar.root").c_str(), "invm_eejj"),
                WeightVar(*eejjRowWeight_ttbar));

  RooDataSet ds_bg_eejj("ds_bg_eejj", "ds_bg_eejj",
                RooArgSet(*eejjMass_bg, *eejjRowWeight_bg),
                WeightVar(*eejjRowWeight_bg));

  ds_bg_eejj.append(ds_DY_eejj);
  ds_bg_eejj.append(ds_ttbar_eejj);

  //// declare model
  // mumujj
  RooRealVar *a_mm_dy = new RooRealVar("a_mm_dy", "a_mm_dy", -2, -5, -1);
  RooRealVar *b_mm_dy = new RooRealVar("b_mm_dy", "b_mm_dy", 5e-1, 2e-1, 15e-1);
  RooExpm *model_mm_dy = new RooExpm("mumujj bg from DY", "exponential bg", *mumujjMass_bg, *a_mm_dy, *b_mm_dy);

  RooRealVar *a_mm_tt = new RooRealVar("a_mm_tt", "a_mm_tt", -2, -5, -1);
  RooRealVar *b_mm_tt = new RooRealVar("b_mm_tt", "b_mm_tt", 5e-01, 2e-2, 15e-1);
  RooExpm *model_mm_tt = new RooExpm("mumujj bg from TTbar", "exponential bg", *mumujjMass_bg, *a_mm_tt, *b_mm_tt);

  RooRealVar *a_mm_tot = new RooRealVar("a_mm_tot", "a_mm_tot", -2, -5, -1);
  RooRealVar *b_mm_tot = new RooRealVar("b_mm_tot", "b_mm_tot", 5e-1, 2e-1, 15e-1);
  RooExpm *model_mm_tot = new RooExpm("mumujj bg from DY+TTbar", "exponential bg", *mumujjMass_bg, *a_mm_tot, *b_mm_tot);

  // eejj
  RooRealVar *a_ee_dy = new RooRealVar("a_ee_dy", "a_ee_dy", -2, -5, -1);
  RooRealVar *b_ee_dy = new RooRealVar("b_ee_dy", "b_ee_dy", 5e-1, 2e-1, 15e-1);
  RooExpm *model_ee_dy = new RooExpm("mumujj bg from DY", "exponential bg", *eejjMass_bg, *a_ee_dy, *b_ee_dy);

  RooRealVar *a_ee_tt = new RooRealVar("a_ee_tt", "a_ee_tt", -2, -5, -1);
  RooRealVar *b_ee_tt = new RooRealVar("b_ee_tt", "b_ee_tt", 5e-01, 2e-2, 15e-1);
  RooExpm *model_ee_tt = new RooExpm("mumujj bg from TTbar", "exponential bg", *eejjMass_bg, *a_ee_tt, *b_ee_tt);

  RooRealVar *a_ee_tot = new RooRealVar("a_ee_tot", "a_ee_tot", -2, -5, -1);
  RooRealVar *b_ee_tot = new RooRealVar("b_ee_tot", "b_ee_tot", 5e-1, 2e-1, 15e-1);
  RooExpm *model_ee_tot = new RooExpm("mumujj bg from DY+TTbar", "exponential bg", *eejjMass_bg, *a_ee_tot, *b_ee_tot);

  //// fit model
  RooFitResult *r_mm_dy = model_mm_dy->fitTo(ds_DY_mumujj, Save(), SumW2Error(kTRUE), Range(700, 2500), Strategy(2));
  RooFitResult *r_mm_tt = model_mm_tt->fitTo(ds_ttbar_mumujj, Save(), SumW2Error(kTRUE), Range(700, 2500), Strategy(2));
  RooFitResult *r_mm_tot = model_mm_tot->fitTo(ds_bg_mumujj, Save(), SumW2Error(kTRUE), Range(700, 2500), Strategy(2));

  RooFitResult *r_ee_dy = model_ee_dy->fitTo(ds_DY_eejj, Save(), SumW2Error(kTRUE), Range(700, 2500), Strategy(2));
  RooFitResult *r_ee_tt = model_ee_tt->fitTo(ds_ttbar_eejj, Save(), SumW2Error(kTRUE), Range(700, 2500), Strategy(2));
  RooFitResult *r_ee_tot = model_ee_tot->fitTo(ds_bg_eejj, Save(), SumW2Error(kTRUE), Range(700, 2500), Strategy(2));

  //// prepare frames for plotting
  RooPlot *frame_mm_dy = mumujjMass_bg->frame(Title("DY mumujj Reco Mass"));
  RooPlot *frame_mm_tt = mumujjMass_bg->frame(Title("TTbar mumujj Reco Mass"));
  RooPlot *frame_mm_tot = mumujjMass_bg->frame(Title("DY+TTbar mumujj Reco Mass"));

  RooPlot *frame_ee_dy = eejjMass_bg->frame(Title("DY eejj Reco Mass"));
  RooPlot *frame_ee_tt = eejjMass_bg->frame(Title("TTbar eejj Reco Mass"));
  RooPlot *frame_ee_tot = eejjMass_bg->frame(Title("DY+TTbar eejj Reco Mass"));

  std::cout << "BELOW IS THE RESULT of DY" << std::endl;
  r_mm_dy->Print();
  r_ee_dy->Print();
  std::cout << "ABOVE IS THE RESULTS of DY" << std::endl;

  std::cout << "BELOW IS THE RESULT of TTbar" << std::endl;
  r_mm_tt->Print();
  r_ee_tt->Print();
  std::cout << "ABOVE IS THE RESULTS of TTbar" << std::endl;

  std::cout << "BELOW IS THE RESULT of DY+TTbar" << std::endl;
  r_mm_tot->Print();
  r_ee_tot->Print();
  std::cout << "ABOVE IS THE RESULTS of DY+TTbar" << std::endl;

  //// Plot on frames
  // plot data on frames
  ds_DY_mumujj.plotOn(frame_mm_dy, Binning(100), DataError(RooAbsData::SumW2));
  ds_ttbar_mumujj.plotOn(frame_mm_tt, Binning(100), DataError(RooAbsData::SumW2));
  ds_bg_mumujj.plotOn(frame_mm_tot, Binning(100), DataError(RooAbsData::SumW2));

  ds_DY_eejj.plotOn(frame_ee_dy, Binning(100), DataError(RooAbsData::SumW2));
  ds_ttbar_eejj.plotOn(frame_ee_tt, Binning(100), DataError(RooAbsData::SumW2));
  ds_bg_eejj.plotOn(frame_ee_tot, Binning(100), DataError(RooAbsData::SumW2));

  // plot fitted pdfs on frames
  model_mm_dy->plotOn(frame_mm_dy);
  model_mm_tt->plotOn(frame_mm_tt);
  model_mm_tot->plotOn(frame_mm_tot);

  model_ee_dy->plotOn(frame_ee_dy);
  model_ee_tt->plotOn(frame_ee_tt);
  model_ee_tot->plotOn(frame_ee_tot);

  double chi2_mm_dy = frame_mm_dy->chiSquare(2);
  double chi2_mm_tt = frame_mm_tt->chiSquare(2);
  double chi2_mm_tot = frame_mm_tot->chiSquare(2);

  double chi2_ee_dy = frame_ee_dy->chiSquare(2);
  double chi2_ee_tt = frame_ee_tt->chiSquare(2);
  double chi2_ee_tot = frame_ee_tot->chiSquare(2);

  std::cout << "chi2_mm_dy: " << chi2_mm_dy <<std::endl;
  std::cout << "chi2_mm_tt: " << chi2_mm_tt <<std::endl;
  std::cout << "chi2_mm_tot: " << chi2_mm_tot <<std::endl;

  std::cout << "chi2_ee_dy: " << chi2_ee_dy <<std::endl;
  std::cout << "chi2_ee_tt: " << chi2_ee_tt <<std::endl;
  std::cout << "chi2_ee_tot: " << chi2_ee_tot <<std::endl;

  model_mm_dy->plotOn(frame_mm_dy, VisualizeError(*r_mm_dy, 1, kFALSE));
  model_mm_tt->plotOn(frame_mm_tt, VisualizeError(*r_mm_tt, 1, kFALSE));
  model_mm_tot->plotOn(frame_mm_tot, VisualizeError(*r_mm_tot, 1, kFALSE));

  model_ee_dy->plotOn(frame_ee_dy, VisualizeError(*r_ee_dy, 1, kFALSE));
  model_ee_tt->plotOn(frame_ee_tt, VisualizeError(*r_ee_tt, 1, kFALSE));
  model_ee_tot->plotOn(frame_ee_tot, VisualizeError(*r_ee_tot, 1, kFALSE));


  TCanvas *canvas = new TCanvas("Test Fit", "Test Fit", 1200, 800);
  canvas->Divide(3,2);
  canvas->cd(1);
  frame_mm_dy->Draw();
  canvas->cd(2);
  frame_mm_tt->Draw();
  canvas->cd(3);
  frame_mm_tot->Draw();

  canvas->cd(4);
  frame_ee_dy->Draw();
  canvas->cd(5);
  frame_ee_tt->Draw();
  canvas->cd(6);
  frame_ee_tot->Draw();
}
