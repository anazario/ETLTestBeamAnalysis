#include "../include/TreeData.h"


inline void TreeData::FillVariables(vector<int> *run_start, vector<int> *run_end){

  int pad_id;

  PrepareTree(run_start,run_end);
  InitBranches();

  outName = Form("pad_tree_%s",tag.Data());
  TFile* file = new TFile(outName, "RECREATE");
  TTree* tree = new TTree("pulse","Tree data per pad");

  tree->Branch("channel",pad_channel,"channel[17][1600]/S");
  tree->Branch("time", pad_time, "time[17][1600]/F");
  tree->Branch("amp", &pad_amp, "amp[17]/F");
  tree->Branch("ntracks", &pad_ntracks, "ntracks[17]/I");
  tree->Branch("npix", &pad_npix, "npix[17]/I");
  tree->Branch("nback", &pad_nback, "nback[17]/I");
  tree->Branch("x_dut", &pad_x_dut, "x_dut[17]/F");
  tree->Branch("y_dut", &pad_y_dut, "y_dut[17]/F");

  int Nentry = tch->GetEntries();

  for(int e = 0; e < Nentry; e++){
    if (e % 1000 == 0) {
      fprintf(stdout, "\r  Processed events: %8d of %8d ", e, Nentry);
    }
    fflush(stdout);
    tch->GetEntry(e);

    for(int ichan = 0; ichan < NUM_CHANNEL; ichan++){
      pad_id = pads->at(ichan);

      //cout << pad_id << endl;
      //pad_channel[pad_id] = channel[ichan];
      //pad_time[pad_id] = time[0];
      pad_amp[pad_id] = amp[ichan];
      pad_ntracks[pad_id] = ntracks;
      pad_npix[pad_id] = npix;
      pad_nback[pad_id] = nback;
      pad_x_dut[pad_id] = x_dut[2];
      pad_y_dut[pad_id] = y_dut[2];
      //tree->Fill();
      //cout << "here" << endl;
      if(ichan == NUM_CHANNEL-1){
	//pad_channel[0] = channel[ichan];
	//pad_time[0] = time[0];
	pad_amp[0] = amp[ichan];
	pad_ntracks[0] = ntracks;
	pad_npix[0] = npix;
	pad_nback[0] = nback;
	pad_x_dut[0] = x_dut[2];
	pad_y_dut[0] = y_dut[2];
      }     
      
    }

    for(int i = 0; i < NUM_PADS+1; i++){
      //cout << "channel: " << pad_channel[i][100] << endl;
      //cout << "time: " << pad_time[i][100] << endl;
      cout << "amp: " << pad_amp[i] << endl;
      cout << "ntracks: " << pad_ntracks[i] << endl;
      cout << "npix: " << pad_npix[i] << endl;
      cout << "nback: " << pad_nback[i] << endl;
      cout << "x: " << pad_x_dut[i] << endl;
      cout << "y: " << pad_y_dut[i] << endl;
    }
    

    cout << "here" << endl; 
    tree->Fill();
  }
  file->Write();
  file->Close();

}

inline void TreeData::PrepareTree(vector<int> *run_start, vector<int> *run_end){

  tch = new TChain("pulse");

  int run_lsts = run_start->size();

  for(int i_runrange = 0; i_runrange < run_lsts; i_runrange++){

    cout << "Preparing runs: " << run_start->at(i_runrange) << " through " << run_end->at(i_runrange) << "." << endl;

    for(int irun = run_start->at(i_runrange); irun <= run_end->at(i_runrange);irun++)
      tch->Add(Form("%s/run_scope%i_info.root",chainPath.Data(),irun));
  }
}

inline void TreeData::InitBranches(){

  //turn off all branches
  tch->SetBranchStatus("*", 0);

  //analog info                                                                                                                                  
  tch->SetBranchStatus("amp", 1); tch->SetBranchAddress("amp", &amp);
  tch->SetBranchStatus("channel", 1); tch->SetBranchAddress("channel", &channel);
  tch->SetBranchStatus("time", 1); tch->SetBranchAddress("time", &time);

  //Track info                                                                                                                                   
  tch->SetBranchStatus("ntracks", 1); tch->SetBranchAddress("ntracks", &ntracks);
  tch->SetBranchStatus("nback", 1); tch->SetBranchAddress("nback", &nback);
  tch->SetBranchStatus("npix", 1); tch->SetBranchAddress("npix", &npix);
  tch->SetBranchStatus("x_dut", 1); tch->SetBranchAddress("x_dut", &x_dut);
  tch->SetBranchStatus("y_dut", 1); tch->SetBranchAddress("y_dut", &y_dut);

  //Run conf info                                                                                                                                
  tch->SetBranchStatus("run", 1); tch->SetBranchAddress("run", &run);
  tch->SetBranchStatus("pads", 1); tch->SetBranchAddress("pads", &pads);
}
