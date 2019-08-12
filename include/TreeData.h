#include<iostream>
#include<vector>

#include<TChain.h>
#include<TString.h>
#include<TFile.h>
#include<TBranch.h>

#define NUM_CHANNEL 4
#define NUM_SAMPLES 1600
#define NUM_PADS 16

using namespace std;

class TreeData{

 public:

  TreeData(){};

  virtual ~TreeData(){};

  TChain* tch;
  TString chainPath;
  TString outName;
  TString tag;

  Int_t pad_ntracks[NUM_PADS+1];
  Int_t pad_npix[NUM_PADS+1];
  Int_t pad_nback[NUM_PADS+1];
  Int_t pad_run[NUM_PADS+1];

  Float_t pad_amp[NUM_PADS+1];
  Float_t pad_y_dut[NUM_PADS+1];
  Float_t pad_x_dut[NUM_PADS+1];
  Float_t *pad_channel[NUM_PADS+1];// = new Float_t[NUM_SAMPLES];
  Float_t *pad_time[NUM_PADS+1];// = new Float_t[NUM_SAMPLES];
  
  void FillVariables(vector<int> *run_start, vector<int> *run_end);

 private:

  Int_t ntracks;
  Int_t npix;
  Int_t nback;
  Int_t run[NUM_CHANNEL];
  
  Float_t amp[NUM_CHANNEL];
  Float_t y_dut[3];
  Float_t x_dut[3];
  Float_t channel[NUM_CHANNEL][NUM_SAMPLES];
  Float_t time[1][NUM_SAMPLES];

  vector<int> *pads = new vector<int>();

  void InitBranches();
  void PrepareTree(vector<int> *run_start, vector<int> *run_end);

};
 
