//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug  6 17:19:27 2019 by ROOT version 6.10/09
// from TTree pulse/Digitized waveforms
// found on file: root://cmseos.fnal.gov//store/group/cmstestbeam/2019_04_April_CMSTiming/KeySightScope/RecoData/TimingDAQRECO/RecoWithTracks/v1/confInfo/run_scope10473_info.root
//////////////////////////////////////////////////////////

#ifndef testbeam_h
#define testbeam_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
using namespace std;

class testbeam {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          i_evt;
   Float_t         channel[4][1600];
   Float_t         time[1][1600];
   Float_t         baseline[4];
   Float_t         baseline_RMS[4];
   Float_t         noise[4];
   Float_t         amp[4];
   Float_t         t_peak[4];
   Float_t         integral[4];
   Float_t         intfull[4];
   Float_t         risetime[4];
   Float_t         decaytime[4];
   Float_t         LP2_5[4];
   Float_t         LP2_10[4];
   Float_t         LP2_15[4];
   Float_t         LP2_20[4];
   Float_t         LP2_25[4];
   Float_t         LP2_30[4];
   Float_t         LP2_35[4];
   Float_t         LP2_40[4];
   Float_t         LP2_45[4];
   Float_t         LP2_50[4];
   Float_t         LP2_55[4];
   Float_t         LP2_60[4];
   Float_t         LP2_65[4];
   Float_t         LP2_30mV[4];
   Float_t         gaus_mean[4];
   Float_t         gaus_sigma[4];
   Float_t         gaus_chi2[4];
   Float_t         xIntercept;
   Float_t         yIntercept;
   Float_t         xSlope;
   Float_t         ySlope;
   Float_t         x_dut[3];
   Float_t         y_dut[3];
   Float_t         chi2;
   Int_t           ntracks;
   Int_t           nplanes;
   Int_t           npix;
   Int_t           nback;
   Int_t           run;
   Int_t           gconf;
   vector<string>  *sensors;
   vector<int>     *pads;
   vector<string>  *sensorsHV;
   vector<int>     *HVs;

   // List of branches
   TBranch        *b_i_evt;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_time;   //!
   TBranch        *b_baseline;   //!
   TBranch        *b_baseline_RMS;   //!
   TBranch        *b_noise;   //!
   TBranch        *b_amp;   //!
   TBranch        *b_t_peak;   //!
   TBranch        *b_integral;   //!
   TBranch        *b_intfull;   //!
   TBranch        *b_risetime;   //!
   TBranch        *b_decaytime;   //!
   TBranch        *b_LP2_5;   //!
   TBranch        *b_LP2_10;   //!
   TBranch        *b_LP2_15;   //!
   TBranch        *b_LP2_20;   //!
   TBranch        *b_LP2_25;   //!
   TBranch        *b_LP2_30;   //!
   TBranch        *b_LP2_35;   //!
   TBranch        *b_LP2_40;   //!
   TBranch        *b_LP2_45;   //!
   TBranch        *b_LP2_50;   //!
   TBranch        *b_LP2_55;   //!
   TBranch        *b_LP2_60;   //!
   TBranch        *b_LP2_65;   //!
   TBranch        *b_LP2_30mV;   //!
   TBranch        *b_gaus_mean;   //!
   TBranch        *b_gaus_sigma;   //!
   TBranch        *b_gaus_chi2;   //!
   TBranch        *b_xIntercept;   //!
   TBranch        *b_yIntercept;   //!
   TBranch        *b_xSlope;   //!
   TBranch        *b_ySlope;   //!
   TBranch        *b_x_dut;   //!
   TBranch        *b_y_dut;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_ntracks;   //!
   TBranch        *b_nplanes;   //!
   TBranch        *b_npix;   //!
   TBranch        *b_nback;   //!
   TBranch        *b_run;   //!
   TBranch        *b_gconf;   //!
   TBranch        *b_sensors;   //!
   TBranch        *b_pads;   //!
   TBranch        *b_sensorsHV;   //!
   TBranch        *b_HVs;   //!

   testbeam(TTree *tree=0);
   virtual ~testbeam();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef testbeam_cxx
testbeam::testbeam(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://cmseos.fnal.gov//store/group/cmstestbeam/2019_04_April_CMSTiming/KeySightScope/RecoData/TimingDAQRECO/RecoWithTracks/v1/confInfo/run_scope10473_info.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://cmseos.fnal.gov//store/group/cmstestbeam/2019_04_April_CMSTiming/KeySightScope/RecoData/TimingDAQRECO/RecoWithTracks/v1/confInfo/run_scope10473_info.root");
      }
      f->GetObject("pulse",tree);

   }
   Init(tree);
}

testbeam::~testbeam()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t testbeam::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t testbeam::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void testbeam::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   sensors = 0;
   pads = 0;
   sensorsHV = 0;
   HVs = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("i_evt", &i_evt, &b_i_evt);
   fChain->SetBranchAddress("channel", channel, &b_channel);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("baseline", baseline, &b_baseline);
   fChain->SetBranchAddress("baseline_RMS", baseline_RMS, &b_baseline_RMS);
   fChain->SetBranchAddress("noise", noise, &b_noise);
   fChain->SetBranchAddress("amp", amp, &b_amp);
   fChain->SetBranchAddress("t_peak", t_peak, &b_t_peak);
   fChain->SetBranchAddress("integral", integral, &b_integral);
   fChain->SetBranchAddress("intfull", intfull, &b_intfull);
   fChain->SetBranchAddress("risetime", risetime, &b_risetime);
   fChain->SetBranchAddress("decaytime", decaytime, &b_decaytime);
   fChain->SetBranchAddress("LP2_5", LP2_5, &b_LP2_5);
   fChain->SetBranchAddress("LP2_10", LP2_10, &b_LP2_10);
   fChain->SetBranchAddress("LP2_15", LP2_15, &b_LP2_15);
   fChain->SetBranchAddress("LP2_20", LP2_20, &b_LP2_20);
   fChain->SetBranchAddress("LP2_25", LP2_25, &b_LP2_25);
   fChain->SetBranchAddress("LP2_30", LP2_30, &b_LP2_30);
   fChain->SetBranchAddress("LP2_35", LP2_35, &b_LP2_35);
   fChain->SetBranchAddress("LP2_40", LP2_40, &b_LP2_40);
   fChain->SetBranchAddress("LP2_45", LP2_45, &b_LP2_45);
   fChain->SetBranchAddress("LP2_50", LP2_50, &b_LP2_50);
   fChain->SetBranchAddress("LP2_55", LP2_55, &b_LP2_55);
   fChain->SetBranchAddress("LP2_60", LP2_60, &b_LP2_60);
   fChain->SetBranchAddress("LP2_65", LP2_65, &b_LP2_65);
   fChain->SetBranchAddress("LP2_30mV", LP2_30mV, &b_LP2_30mV);
   fChain->SetBranchAddress("gaus_mean", gaus_mean, &b_gaus_mean);
   fChain->SetBranchAddress("gaus_sigma", gaus_sigma, &b_gaus_sigma);
   fChain->SetBranchAddress("gaus_chi2", gaus_chi2, &b_gaus_chi2);
   fChain->SetBranchAddress("xIntercept", &xIntercept, &b_xIntercept);
   fChain->SetBranchAddress("yIntercept", &yIntercept, &b_yIntercept);
   fChain->SetBranchAddress("xSlope", &xSlope, &b_xSlope);
   fChain->SetBranchAddress("ySlope", &ySlope, &b_ySlope);
   fChain->SetBranchAddress("x_dut", x_dut, &b_x_dut);
   fChain->SetBranchAddress("y_dut", y_dut, &b_y_dut);
   fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
   fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
   fChain->SetBranchAddress("nplanes", &nplanes, &b_nplanes);
   fChain->SetBranchAddress("npix", &npix, &b_npix);
   fChain->SetBranchAddress("nback", &nback, &b_nback);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("gconf", &gconf, &b_gconf);
   fChain->SetBranchAddress("sensors", &sensors, &b_sensors);
   fChain->SetBranchAddress("pads", &pads, &b_pads);
   fChain->SetBranchAddress("sensorsHV", &sensorsHV, &b_sensorsHV);
   fChain->SetBranchAddress("HVs", &HVs, &b_HVs);
   Notify();

   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("ntracks",1);
   fChain->SetBranchStatus("npix",1);
   fChain->SetBranchStatus("nback",1);
   fChain->SetBranchStatus("amp",1);
   fChain->SetBranchStatus("x_dut",1);
   fChain->SetBranchStatus("y_dut",1);
   fChain->SetBranchStatus("LP2_20",1);
   fChain->SetBranchStatus("gaus_mean",1);
   fChain->SetBranchStatus("pads",1);
}

Bool_t testbeam::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void testbeam::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t testbeam::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef testbeam_cxx
