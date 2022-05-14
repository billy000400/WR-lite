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

RooExpPoly2* ExpPoly2_init(RooRealVar* x, std::string label);

void testFit_bg_ExpPoly2()
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
  RooExpPoly2 *model1 = ExpPoly2_init(mumujjMass_DY, "DY");
  RooExpPoly2 *model2 = ExpPoly2_init(mumujjMass_ttbar, "TTbar");
  RooExpPoly2 *model3 = ExpPoly2_init(mumujjMass_bg, "FullBg");

  // fit model
  RooFitResult *r1 = model1->fitTo(ds_DY_mumujj, Save(), SumW2Error(kTRUE), Range(500, 3000));
  RooFitResult *r2 = model2->fitTo(ds_ttbar_mumujj, Save(), SumW2Error(kTRUE), Range(500, 3000));
  RooFitResult *r3 = model3->fitTo(ds_bg_mumujj, Save(), SumW2Error(kTRUE), Range(500, 3000));

  // prepare frames for plotting
  RooPlot *frame1 = mumujjMass_DY->frame(Title("DY mumujj Reco Mass"));
  RooPlot *frame2 = mumujjMass_ttbar->frame(Title("TTbar mumujj Reco Mass"));
  RooPlot *frame3 = mumujjMass_bg->frame(Title("DY+TTbar mumujj Reco Mass"));

  //// Plot on frames
  // plot data on frames
  ds_DY_mumujj.plotOn(frame1, Binning(128));
  ds_ttbar_mumujj.plotOn(frame2, Binning(128));
  ds_bg_mumujj.plotOn(frame3, Binning(128));

  // plot fitted pdfs on frames
  model1->plotOn(frame1);
  model2->plotOn(frame2);
  model3->plotOn(frame3);


  TCanvas *canvas = new TCanvas("Test Fit", "Test Fit", 1500, 500);
  canvas->Divide(3,1);
  canvas->cd(1);
  frame1->Draw();
  canvas->cd(2);
  frame2->Draw();
  canvas->cd(3);
  frame3->Draw();

  TCanvas *logCanvas = new TCanvas("log", "log", 1500, 500);
  logCanvas->Divide(3,1);
  logCanvas->cd(1);
  logCanvas->GetPad(1)->SetLogy();
  frame1->Draw();
  logCanvas->cd(2);
  logCanvas->GetPad(2)->SetLogy();
  frame2->Draw();
  logCanvas->cd(3);
  logCanvas->GetPad(3)->SetLogy();
  frame3->Draw();
}

RooExpPoly2* ExpPoly2_init(RooRealVar* x, std::string label)
{
    RooRealVar* a = new RooRealVar((std::string("a")+label).c_str(), label.c_str(), 1., 1e-7, 5.);
    RooRealVar* b = new RooRealVar((std::string("b")+label).c_str(), label.c_str(), -6e5., -6e7, 6e7);

    return new RooExpPoly2((std::string("ExpPoly2_")+label).c_str(), label.c_str(), *x, *a, *b);
}
