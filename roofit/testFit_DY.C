/**
 * @Author: Billy Li <billyli>
 * @Date:   05-03-2022
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 05-06-2022
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


void testFit_DY()
{
  // Preparing RooRealVars
  RooRealVar* lljjRecoMass = new RooRealVar("lljjRecoMass", "lljjRecoMass", 400, 3000);
  RooRealVar* ljjRecoMass_Res = new RooRealVar("ljjRecoMass_Res", "ljjRecoMass_Res", 50, 3000);
  RooRealVar* ljjRecoMass_SpRes = new RooRealVar("ljjRecoMass_SpRes", "ljjRecoMass_SpRes", 50, 3000);
  RooRealVar* rowWeight = new RooRealVar("rowWeight", "rowWeight", -1.5, 1.5);

  // importing ntuples into RooDataSet
  std::string prefix = "../analysis/allEvents/";

  RooDataSet ds_lljjRecoMass("ds_lljjRecoMass", "ds_lljjRecoMass",
                RooArgSet(*lljjRecoMass, *rowWeight),
                ImportFromFile((prefix+"fullDY.root").c_str(), "fullBgRecoMass"),
                WeightVar(*rowWeight));

  RooDataSet ds_ljjRecoMass_Res("ds_ljjRecoMass_Res", "ds_ljjRecoMass_Res",
                RooArgSet(*ljjRecoMass_Res, *rowWeight),
                ImportFromFile((prefix+"fullDY.root").c_str(), "fullBgRecoMass"),
                WeightVar(*rowWeight));

  RooDataSet ds_ljjRecoMass_SpRes("ds_ljjRecoMass_SpRes", "ds_ljjRecoMass_SpRes",
                RooArgSet(*ljjRecoMass_SpRes, *rowWeight),
                ImportFromFile((prefix+"fullDY.root").c_str(), "fullBgRecoMass"),
                WeightVar(*rowWeight));


  // declare model
  RooRealVar *c = new RooRealVar("c", "c", -5e-2, -1e-1, -1e-7);
  RooExponential *model = new RooExponential("exponential DY", "exponential DY", *lljjRecoMass, *c);

  // fit model
  RooFitResult *r = model->fitTo(ds_lljjRecoMass, Save(), Range(1000, 3000));

  // prepare frames for plotting
  RooPlot *frame1 = lljjRecoMass->frame(Title("DY lljj Reco Mass (top 2 pT lepton)"));
  RooPlot *frame2 = ljjRecoMass_Res->frame(Title("DY ljj Mass, Reco by Resolved NN"));
  RooPlot *frame3 = ljjRecoMass_SpRes->frame(Title("DY ljj Mass, Reco by SuperResolved NN"));


  //// Plot on frames
  // plot data on frames
  ds_lljjRecoMass.plotOn(frame1, Binning(128));
  ds_ljjRecoMass_Res.plotOn(frame2, Binning(128));
  ds_ljjRecoMass_SpRes.plotOn(frame3, Binning(128));
  // plot fitted pdfs on frames
  model->plotOn(frame1);


  TCanvas *canvas = new TCanvas("Test Fit", "Test Fit", 1500, 500);
  canvas->Divide(3,1);
  canvas->cd(1);
  frame1->Draw();
  canvas->cd(2);
  frame2->Draw();
  canvas->cd(3);
  frame3->Draw();
}
