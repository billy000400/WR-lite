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
  RooRealVar* lljjRecoMass = new RooRealVar("lljjRecoMass", "lljjRecoMass", 0, 3000);
  RooRealVar* ljjRecoMass_Res = new RooRealVar("ljjRecoMass_Res", "ljjRecoMass_Res", 0, 2500);
  RooRealVar* ljjRecoMass_SpRes = new RooRealVar("ljjRecoMass_SpRes", "ljjRecoMass_SpRes", 0, 2500);
  RooRealVar* rowWeight = new RooRealVar("rowWeight", "rowWeight", -1.5, 1.5);

  // importing ntuples into RooDataSet
  std::string prefix = "../analysis/allEvents";
  RooDataSet dsWeight("dsWeight", "dsWeight",
                    RooArgSet(*rowWeight),
                    ImportFromFile((prefix+"fullDY.root").c_str(), "analysis/fullRowWeight"));
  RooDataSet dsDY("dsDY", "dsDY",
                RooArgSet(*lljjRecoMass, *ljjRecoMass_Res, *ljjRecoMass_SpRes),
                ImportFromFile((prefix+"fullDY.root").c_str(), "analysis/bgRecoMass"),
                WeightVar(*rowWeight));


  RooPlot *frame1 = lljjRecoMass->frame(Title("DY lljj Reco Mass (top 2 pT lepton)"));
  dsDY.plotOn(frame1, Binning(128));
  RooPlot *frame2 = ljjRecoMass_Res->frame(Title("DY ljj Mass, Reco by Resolved NN"));
  dsDY.plotOn(frame2, Binning(128));
  RooPlot *frame3 = ljjRecoMass_SpRes->frame(Title("DY ljj Mass, Reco by SuperResolved NN"));
  dsDY.plotOn(frame3, Binning(128));



  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 1500, 500);
  c->Divide(1,3);
  c->cd(1);
  frame1->Draw();
  c->cd(2);
  frame2->Draw();
  c->cd(3);
  frame3->Draw();
}
