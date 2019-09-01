#include "../include/TreeData.h"

inline void TreeData::FillTrees(){

  int pad_id = -999;
  Float_t* r;

  PrepareTree();
  InitBranches();

  TFile* file[NUM_PADS+1];
  TTree* tree[NUM_PADS+1];

  outDir = Form("../ROOT/%s_trees/",tag.Data());

  if(gSystem->OpenDirectory(outDir) == 0)
    gSystem->mkdir(outDir);
  
  for(int ipad = 1; ipad < NUM_PADS+1; ipad++){

    if (!isReduced)
      outName = Form(outDir+"%s_tree_pad%i.root",tag.Data(),ipad);
    else
      outName = Form(outDir+"%s_tree_pad%i_reduced.root",tag.Data(),ipad);

    file[ipad] = new TFile(outName, "RECREATE");
    tree[ipad] = new TTree("pulse",Form("Pad %i Tree",ipad));
   
    if (!isReduced){ 
      tree[ipad]->Branch("i_evt",&i_evt[ipad],"i_evt/I");
      tree[ipad]->Branch("channel",&pad_channel[ipad],"channel[1600]/F");
      tree[ipad]->Branch("photek",&mcp_channel[ipad],"photek[1600]/F");
      tree[ipad]->Branch("time", &pad_time[ipad], "time[1600]/F");
    }

    tree[ipad]->Branch("amp", &pad_amp[ipad], "amp[2]/F");
    tree[ipad]->Branch("ntracks", &pad_ntracks[ipad], "ntracks/I");
    tree[ipad]->Branch("npix", &pad_npix[ipad], "npix/I");
    tree[ipad]->Branch("nback", &pad_nback[ipad], "nback/I");
    tree[ipad]->Branch("x_dut", &pad_x_dut[ipad], "x_dut/F");
    tree[ipad]->Branch("y_dut", &pad_y_dut[ipad], "y_dut/F");
  }

  int Nentry = tch->GetEntries();
  int Nentry2 = 10000;

  for(int e = 0; e < Nentry; e++){
    if (e % 1000 == 0) {
      fprintf(stdout, "\r  Processed events: %8d of %8d ", e, Nentry);
    }
    fflush(stdout);
    tch->GetEntry(e);

    r = RotateXY(x_dut[2],y_dut[2]);

    for(int ichan = 0; ichan < NUM_CHANNEL-1; ichan++){
      pad_id = pads->at(ichan);

      pad_ntracks[pad_id] = ntracks;
      pad_npix[pad_id] = npix;
      pad_nback[pad_id] = nback;
      pad_x_dut[pad_id] = r[0];
      pad_y_dut[pad_id] = r[1];
      pad_amp[pad_id][0] = amp[ichan];
      pad_amp[pad_id][1] = amp[NUM_CHANNEL-1];

      if (!isReduced){
	for(int samples = 0; samples < NUM_SAMPLES; samples++){
	  pad_time[pad_id][samples] = time[0][samples];
	  mcp_channel[pad_id][samples] = channel[NUM_CHANNEL-1][samples];
	  pad_channel[pad_id][samples] = channel[ichan][samples];
	}
      }

      tree[pad_id]->Fill();
      i_evt[pad_id]++;
    }
  }
  cout << endl;

  for(int ipad = 1; ipad < NUM_PADS+1; ipad++){
    file[ipad]->Write();
    file[ipad]->Close();
  }
}

inline void TreeData::PrepareTree(){

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

  if (!isReduced){                                                                                          
    tch->SetBranchStatus("channel", 1); tch->SetBranchAddress("channel", &channel);
    tch->SetBranchStatus("time", 1); tch->SetBranchAddress("time", &time);
  }

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

inline Float_t* TreeData::RotateXY(Float_t x, Float_t y){

  Float_t costheta = TMath::Cos(angle);
  Float_t sintheta = TMath::Sin(angle);

  rp[0] = x*costheta-y*sintheta;
  rp[1] = y*costheta+x*sintheta;

  return rp;
}
