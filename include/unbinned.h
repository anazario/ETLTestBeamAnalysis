#include <iostream>
#include <string>
#include <vector>

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

using namespace RooFit;
using namespace std;

class Unbinned{

 public:

  TChain *tch;

  TString chainPath;
  TString tag;
  TString outDir;
  TString inDir;

  vector<int> *run_start;
  vector<int> *run_end;

  static const int npad = 16;
  static const int nlgads = 3;
  static const int geom = 4;

  TString name;

  vector<TH1F*> pad_amp;
  vector<TH1F*> pad_noise;
  vector<TH2F*> pad2D;
  vector<TH2F*> pad2D_all;

  vector<RooPlot*> noise_plot;
  vector<RooPlot*> mip_plot;

  vector<RooWorkspace*> ws;
  vector<RooRealVar*> amp_rrvar;
  vector<RooArgSet*> amp_argset;
  vector<RooDataSet*> mip_dataset;
  vector<RooRealVar*> ml,sl,mg,sg;
  vector<RooLandau*> landau;
  vector<RooGaussian*> gauss;
  vector<RooFFTConvPdf*> lxg;
  vector<RooAddPdf*> pdf_lxg;

  vector<RooDataSet*> noise_dataset;
  vector<RooKeysPdf*> rkpdf;
  vector<RooRealVar*> Nn,Ns;

  vector<float> xcuts;
  vector<float> ycuts;

  float angle;
  int pad_id;

  Float_t amp[nlgads];
  Int_t ntracks;
  Int_t nback;
  Int_t npix;
  Float_t y_dut[3];
  Float_t x_dut[3];
  Float_t geom_cuts[npad+1][geom];


  Int_t run;
  Float_t satCut;
  Float_t prec;
  vector<int> *pads = new vector<int>();

  void PrintGeomMap();
  void InitBranches();
  void PrepareTree();
  void PrepareCuts();
  void SavePadNoise(Float_t xmin, Float_t xmax, int nbins);
  void SavePadSignalFit(Float_t xmin, Float_t xmax, int nbins);
  void Save2DPads(Float_t xl, Float_t xh, Float_t yl, Float_t yh);

 private:
  //pad# and index of xmin,xmax,ymin,ymax values
  int geom_map[npad+1][geom] = {
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


