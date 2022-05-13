/**
 * @Author: Billy Li <billyli>
 * @Date:   05-03-2022
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 05-12-2022
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


void testFit_bg()
{
  // Preparing RooRealVars
  RooRealVar* mumujjMass_DY = new RooRealVar("invm_mumujj",\
                            "invm reco from DY mumujj", 400, 3000);
  RooRealVar* mumujjRowWeight_DY = new RooRealVar("rowWeight",\
                            "row weight for DY ntuple mumujj rows", -1.5, 1.5);

  RooRealVar* mumujjMass_ttbar = new RooRealVar("invm_mumujj",\
                            "invm reco from ttbar mumujj", 400, 3000);
  RooRealVar* mumujjRowWeight_ttbar = new RooRealVar("rowWeight",\
                            "row weight for ttbar ntuple mumujj rows", -1.5, 1.5);

  RooRealVar* mumujjMass_bg = new RooRealVar("invm_mumujj",\
                            "invm reco from bg mumujj", 400, 3000);
  RooRealVar* mumujjRowWeight_bg = new RooRealVar("rowWeight",\
                            "row weight for bg ntuple mumujj rows", -1.5, 1.5);

  // RooRealVar* eejjMass_DY = new RooRealVar();
  // RooRealVar* eejjRowWeight_DY = new RooRealVar("rowWeight", "rowWeight", -1.5, 1.5);

  // importing ntuples into RooDataSet
  std::string prefix = "../analysis/allEvents/";

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


  // declare model
  RooRealVar *c = new RooRealVar("c", "c", -5e-2, -1e-1, -1e-7);
  RooExponential *model = new RooExponential("exponential bg", "exponential bg", *mumujjMass_bg, *c);

  // fit model
  RooFitResult *r = model->fitTo(ds_bg_mumujj, Save(), SumW2Error(kTRUE), Range(500, 3000));

  // prepare frames for plotting
  RooPlot *frame1 = mumujjMass_bg->frame(Title("DY+TTbar mumujj Reco Mass"));



  //// Plot on frames
  // plot data on frames
  ds_bg_mumujj.plotOn(frame1, Binning(128));

  // plot fitted pdfs on frames
  model->plotOn(frame1);


  TCanvas *canvas = new TCanvas("Test Fit", "Test Fit", 1500, 500);
  canvas->Divide(1,1);
  canvas->cd(1);
  frame1->Draw();
  // canvas->cd(2);
  // frame2->Draw();
  // canvas->cd(3);
  // frame3->Draw();
}
