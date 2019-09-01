#include <iostream>
#include <string>
#include <vector>
#include <functional>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TColor.h>
#include "TSystem.h"

#include "RooFit.h"
#include "RooAbsArg.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooFFTConvPdf.h"
#include "RooLandau.h"
#include "RooKeysPdf.h"
#include "RooBinning.h"
#include "RooWorkspace.h"

#include "testbeam.h"

#define NUM_CHANNEL 2
#define NUM_PADS 16
#define GEOM 4

using namespace RooFit;
using namespace std;

class Unbinned{

 public:

  Unbinned(){};

  virtual ~Unbinned(){};

  TChain *tch;

  TString chainPath;
  TString tag;
  TString outDir;
  TString inDir;

  Bool_t isDebug = false;

  vector<Float_t> xcuts;
  vector<Float_t> ycuts;
  vector<Float_t> frame_margins;

  Float_t geom_cuts[NUM_PADS+1][GEOM];

  Float_t satCut;
  Int_t nbins = 100;
  Float_t xl = 0.;
  Float_t xh = -999.;
  Float_t prec = 0.05;

  void PrintGeomMap();
  void SavePadInfo();
  void SavePadNoise();
  void SavePadSignal();

 private:

  TString name;

  Int_t pad_id;
  Int_t vecsize = 0;

  TFile* fIn;
  TFile* fOut;

  //Pad cuts
  Float_t xmin = -999.;
  Float_t xmax = -999.;
  Float_t ymin = -999.;
  Float_t ymax = -999.;

  //new branches
  Int_t ntracks;
  Int_t npix;
  Int_t nback;
                                                                                                              
  Float_t amp[NUM_CHANNEL];
  Float_t x_dut;
  Float_t y_dut;

  //Pad histograms
  vector<TH1F*> pad_amp;
  vector<TH1F*> pad_noise;
  vector<TH2F*> pad2D;
  vector<TH2F*> pad2D_all;

  //Roofit vectors
  vector<RooRealVar*> amp_rrvar;
  vector<RooArgSet*> amp_argset;
  vector<RooDataSet*> noise_dataset;
  vector<RooDataSet*> mip_dataset;
  vector<RooRealVar*> ml,sl,mg,sg;
  vector<RooLandau*> landau;
  vector<RooGaussian*> gauss;
  vector<RooFFTConvPdf*> lxg;
  vector<RooAddPdf*> pdf_lxg;
  vector<RooKeysPdf*> rkpdf;
  vector<RooRealVar*> Nn,Ns;
  vector<RooPlot*> rplot;
  vector<RooWorkspace*> ws;
  vector<RooWorkspace*> wsIn;

  void InitBranches();
  void PrepareTree();
  void PrepareCuts();
  void AssignCuts();
  void PrepVectors();
  void ClearVectors();

  void AssignCuts(Int_t pad_id);
  void PadInfoSelection(Int_t pad_id);
  void NoiseSelection(Int_t pad_id);
  void SignalSelection(Int_t pad_id);
  void ProcessEntries(Int_t pad_id, std::function<void(Int_t)> selection);

  void ProcessPadInfo(Int_t pad_id);
  void ProcessPadNoise(Int_t pad_id);
  void FitPadNoise(Int_t pad_id);
  void ProcessPadSignal(Int_t pad_id);
  void FitPadSignal(Int_t pad_id);

  //pad# and index of xmin,xmax,ymin,ymax values
  int geom_map[NUM_PADS+1][GEOM] = {
    {0,0,0,0},
    {3,4,2,3},//pad1
    {3,4,3,4},//pad2
    {2,3,2,3},//pad3
    {2,3,3,4},//pad4
    {1,2,3,4},//pad5
    {1,2,2,3},//pad6
    {0,1,3,4},//pad7
    {0,1,2,3},//pad8
    {0,1,1,2},//pad9
    {0,1,0,1},//pad10
    {1,2,1,2},//pad11
    {1,2,0,1},//pad12
    {2,3,0,1},//pad13
    {2,3,1,2},//pad14
    {3,4,0,1},//pad15
    {3,4,1,2} //pad16
  };
};


