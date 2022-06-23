/**
 * @Author: Billy Li <billyli>
 * @Date:   05-03-2022
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 06-17-2022
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
  RooRealVar *a_mm = new RooRealVar("a_mm", "a_mm", -2e-1, -20e-1, -5e-2);
  RooRealVar *b_mm = new RooRealVar("b_mm", "b_mm", 5e-01, 2e-1, 8e-1);
  RooExpm *model_mm = new RooExpm("mumujj bg from DY+ttbar ", "exponential bg", *mumujjMass_bg, *a_mm, *b_mm);
  // eejj
  RooRealVar *a_ee = new RooRealVar("a_ee", "a_ee", -2e-01, -20e-1, -5e-2);
  RooRealVar *b_ee = new RooRealVar("b_ee", "b_ee", 5e-01, 2e-2, 8e-1);
  RooExpm *model_ee = new RooExpm("eejj bg from DY+ttbar ", "exponential bg", *eejjMass_bg, *a_ee, *b_ee);

  //// fit model
  RooFitResult *r_mm = model_mm->fitTo(ds_bg_mumujj, Save(), SumW2Error(kTRUE), Range(650, 2500), Strategy(2));
  RooFitResult *r_ee = model_ee->fitTo(ds_bg_eejj, Save(), SumW2Error(kTRUE), Range(650, 2500), Strategy(2));

  //// prepare frames for plotting
  RooPlot *frame_mm = mumujjMass_bg->frame(Title("DY+TTbar mumujj Reco Mass"));
  RooPlot *frame_ee = eejjMass_bg->frame(Title("DY+TTbar eejj Reco Mass"));

  std::cout << "BELOW IS THE RESULT" << std::endl;
  r_mm->Print();
  r_ee->Print();
  std::cout << "ABOVE IS THE RESULTS" << std::endl;

  //// Plot on frames
  // plot data on frames
  ds_bg_mumujj.plotOn(frame_mm, Binning(100), DataError(RooAbsData::SumW2));
  ds_bg_eejj.plotOn(frame_ee, Binning(100), DataError(RooAbsData::SumW2));

  // plot fitted pdfs on frames
  model_mm->plotOn(frame_mm, VisualizeError(*r_mm, 1, kFALSE));
  model_ee->plotOn(frame_ee, VisualizeError(*r_ee, 1, kFALSE));

  double chi2_mumu = frame_mm->chiSquare(2);
  double chi2_ee = frame_ee->chiSquare(2);
  std::cout << "chi2_mumu_bg: " << chi2_mumu <<std::endl;
  std::cout << "chi2_ee_bg: " << chi2_ee <<std::endl;


  TCanvas *canvas = new TCanvas("Test Fit", "Test Fit", 1000, 500);
  canvas->Divide(2,1);
  canvas->cd(1);
  frame_mm->Draw();
  canvas->cd(2);
  frame_ee->Draw();
}
