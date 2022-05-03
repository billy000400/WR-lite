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

RooAddPdf* DoubleCB(RooRealVar* rrv_x);

void testFit_DY()
{
  // Preparing RooRealVars
  RooRealVar* lljjRecoMass = new RooRealVar("lljjRecoMass", "lljjRecoMass", 0, 3000);
  RooRealVar* ljjRecoMass_Res = new RooRealVar("ljjRecoMass_Res", "ljjRecoMass_Res", 0, 2500);
  RooRealVar* ljjRecoMass_SpRes = new RooRealVar("ljjRecoMass_SpRes", "ljjRecoMass_SpRes", 0, 2500);

  // importing ntuples into RooDataSet
  std::string prefix = "../";
  std::vector<std::string> folders;
  folders.push_back((prefix+"DYJetsToLL_M-50_HT-100to200/").c_str());
  folders.push_back((prefix+"DYJetsToLL_M-50_HT-100to200/").c_str());
  folders.push_back((prefix+"DYJetsToLL_M-50_HT-100to200/").c_str());
  folders.push_back((prefix+"DYJetsToLL_M-50_HT-100to200/").c_str());
  folders.push_back((prefix+"DYJetsToLL_M-50_HT-100to200/").c_str());
  RooDataSet ds1("ds1", "ds1",
                RooArgSet(*lljjRecoMass, *ljjRecoMass_Res, *ljjRecoMass_SpRes),
                ImportFromFile(prefix"dyTest.root","analysis/bgRecoMass"));


  RooPlot *frame1 = lljjRecoMass->frame(Title("DY lljj Reco Mass (top 2 pT lepton)"));
  ds1.plotOn(frame1, Binning(128));
  RooPlot *frame2 = ljjRecoMass_Res->frame(Title("DY ljj Mass, Reco by Resolved NN"));
  ds1.plotOn(frame2, Binning(128));
  RooPlot *frame3 = ljjRecoMass_SpRes->frame(Title("DY ljj Mass, Reco by SuperResolved NN"));
  ds1.plotOn(frame3, Binning(128));

  // preparing the background distribution
  RooAddPdf* WR_pdf = DoubleCB(lljjRecoMass);

  // fit distribution to data
  WR_pdf->fitTo(ds1);

  // Draw ntuples
  WR_pdf->plotOn(frame1);




  TCanvas *c = new TCanvas("Test Fit", "Test Fit", 1000, 800);
  c->Divide(2,2);
  c->cd(1);
  frame1->Draw();
  c->cd(2);
  frame2->Draw();
  c->cd(3);
  frame3->Draw();
}
