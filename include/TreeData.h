#include<iostream>
#include<vector>

#include<TSystem.h>
#include<TChain.h>
#include<TString.h>
#include<TFile.h>
#include<TBranch.h>
#include<TMath.h>

#define NUM_CHANNEL 4
#define NUM_SAMPLES 1600
#define NUM_PADS 16
#define TREE_CH 2

using namespace std;

class TreeData{

 public:

  TreeData(){};

  virtual ~TreeData(){};

  TChain* tch;
  TString chainPath;
  TString outDir;
  TString outName;
  TString tag;

  vector<int> *run_start;
  vector<int> *run_end;

  float angle = 0;

  bool isReduced;

  void FillTrees();
  Float_t* RotateXY(Float_t x, Float_t y);

 private:

  Float_t* rp = new Float_t[2];

  //original branches
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

  //new branches
  Int_t pad_ntracks[NUM_PADS];
  Int_t pad_npix[NUM_PADS];
  Int_t pad_nback[NUM_PADS];

  Float_t zeroes[1600] = {};
  Int_t i_evt[NUM_PADS] = {};

  Float_t pad_amp[NUM_PADS][TREE_CH];
  Float_t pad_channel[NUM_PADS][NUM_SAMPLES];
  Float_t mcp_channel[NUM_PADS][NUM_SAMPLES];
  Float_t pad_time[NUM_PADS][NUM_SAMPLES];
  Float_t pad_y_dut[NUM_PADS];
  Float_t pad_x_dut[NUM_PADS];

  void InitBranches();
  void PrepareTree();

};
 
