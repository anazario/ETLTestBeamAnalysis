#ifndef SIGNAL_FIT_H
#define SIGNAL_FIT_H

#include <iostream>
#include <string>
#include <vector>

#include <TROOT.h>

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

#include "Plot.h"

using namespace RooFit;
using namespace std;

class SignalFit{

 public:

  SignalFit(){
    ml = nullptr;
    sl = nullptr;
    mg = nullptr;
    sg = nullptr;
    Nn = nullptr;
    Ns = nullptr;
    landau = nullptr;
    gauss = nullptr;
    lxg = nullptr;
    pdf_lxg = nullptr;
};

  virtual ~SignalFit(){
    if(ml != nullptr)
      delete ml;
    if(sl != nullptr)
      delete sl;
    if(mg != nullptr)
      delete mg;
    if(sg != nullptr)
      delete sg;
    if(Nn != nullptr)
      delete Nn;
    if(Ns != nullptr)
      delete Ns;
    if(landau != nullptr)
      delete landau;
    if(gauss != nullptr)
      delete gauss;
    if(lxg != nullptr)
      delete lxg;
    if(pdf_lxg != nullptr)
      delete pdf_lxg;
};

  Float_t GetMPV(){return MPV;};
  Float_t GetMPVerr(){return MPV_err;};
  Float_t GetSigmaL(){return sigmaL;};
  Float_t GetSigmaG(){return sigmaG;};
  RooPlot* GetPlot(){return rplot;};

  void SetFitNorm(Int_t norm_noise, Int_t norm_sig);
  void SetFitRange(Float_t xmin, Float_t xmax, TString label);
  void AmplitudeFit(RooDataSet* data, RooRealVar* amp, RooKeysPdf* error_fit, Float_t mean, TString label);
  void NoiseFit(RooDataSet* data, RooRealVar* amp, TString label);
  void PlotFit(Float_t xl, Float_t xh, Int_t nbins, TString label);  
  void PlotSignal(Float_t xl, Float_t xh, Int_t nbins, TString label);
  void PlotNoise(Float_t xl, Float_t xh, Int_t nbins, TString label);
  void PrintSignalPlot(Float_t xl, Float_t xh, Int_t nbins, TString label);
  void PrintNoisePlot(Float_t xl, Float_t xh, Int_t nbins, TString label);
  
 private:

  Int_t noise_n = -999;
  Int_t sig_n = -999;

  Float_t MPV = -999.;
  Float_t MPV_err = -999.;
  Float_t sigmaL = -999.;
  Float_t sigmaG = -999.;
  Float_t xl = -999.;
  Float_t xh = -999.;

  Bool_t set_range = false;

  TString name;

  RooRealVar *amp_rrvar;
  RooRealVar *ml;
  RooRealVar *sl;
  RooRealVar *mg;
  RooRealVar *sg;
  RooRealVar *Nn;
  RooRealVar *Ns;
  RooDataSet *dataset;
  RooKeysPdf *rkpdf;
  RooLandau *landau;
  RooGaussian *gauss;
  RooFFTConvPdf *lxg;
  RooAddPdf *pdf_lxg;
  RooPlot *rplot;

  void PrintParameters(TString label);
  void PrintPlot(Float_t xl, Float_t xh, Int_t nbins, TString label,std::function<void(Float_t,Float_t,Int_t,TString)> PlotFit);
  void CleanUp();
};

#endif
