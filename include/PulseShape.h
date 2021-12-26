#ifndef PULSESHAPE_H
#define PULSESHAPE_H

#define NUM_CHANNEL 4
#define NUM_SAMPLES 1600

#include <iostream>
#include <map>
#include <vector>
#include <math.h>

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TAxis.h"
#include "TH1.h"
#include "TTree.h"
#include "TRandom.h"
#include "TEventList.h"
#include "TColor.h"
//#include "Plot.h"
#include "TStyle.h"
#include "TLatex.h"
#include "PulseTools.h"


class PulseShape{

 public:

  PulseShape() {};

  virtual ~PulseShape() {};

  void PlotHists();
  void GetMeanErrors(vector<TH1F> pulse_hist, float *&mean, float *&high_err, float *&low_err, int CI);
  void GetErrors(float *&low_err, float *&high_err, float CI);
  void GetMaxT0andAmp0(float* sample, float& t0, float& amp0);
  void GetMaxT0Amp0Vec(vector<float>& t0_vec, vector<float>& amp0_vec);
  float* GetMeanPulse(vector<float> t0, vector<float> max_amp, float step_size, std::function<float*(float*,float,float,float)> PulseType);
  static float* MakeInterpolation(float* sample, float t0, float amp_peak, float step_size);
  static float* GetPulseCDF(float* sample, float t0, float max_amp, float step_size);

  void PlotAllPulses(vector<float> t0, vector<float> max_amp);
  void PlotAllCDFs(vector<float> t0, vector<float> max_amp);
  void PlotPulseMean(vector<float> t0, vector<float> max_amp);
  void PlotCDFMean(vector<float> t0, vector<float> max_amp);
  void PlotPulseMeanErr(vector<float> t0, vector<float> max_amp);
  void PlotCDFMeanErr(vector<float> t0, vector<float> max_amp);
  void PlotRangeCDF(float xmin, float xmax, int total_slices);
  void PlotRangePulse(float xmin, float xmax, int total_slices);

 private:

  vector<float*> sample; 
  vector<float> t_peak;
  vector<float> max_amp;

  vector<TH1F> hist_vec;

  void PlotAll(vector<float> t0, vector<float> max_amp, float step_size, TString name, TString y_label, 
	       std::function<float*(float*,float,float,float)> PulseType);
  void PlotMean(vector<float> t0, vector<float> max_amp, float step_size, TString name, TString y_label, 
		std::function<float*(float*,float,float,float)> PulseType);
  void PlotMeanErr(vector<float> t0, vector<float> max_amp, float step_size, TString name, TString y_label, TString opt, 
		   std::function<float*(float*,float,float,float)> PulseType);
  void PlotRange(float xmin, float xmax, int total_slices, TString name, TString y_label, TString opt, 
		 std::function<float*(float*,float,float,float)> PulseType);
};

#endif
