#include "../include/Unbinned.h"

inline void Unbinned::PadInfoSelection(Int_t pad_id){

  Float_t prec_x = prec*2.;
  Float_t prec_y = prec;

  AssignCuts(pad_id);

  if(x_dut > (xmin-prec_x) && x_dut < (xmax+prec_x) && y_dut > (ymin-prec_y) && y_dut < (ymax+prec_y) && ntracks == 1
     && npix > 0 && nback > 0 && amp[0] > 0. && amp[0] < satCut){
    pad2D[pad_id]->Fill(x_dut,y_dut);
  }
  if(ntracks == 1 && npix > 0 && nback > 0 && amp[0] > 40. && amp[0] < satCut)
    pad2D_all[pad_id]->Fill(x_dut,y_dut);

  if(x_dut > (xmin+prec_x) && x_dut < (xmax-prec_x) && y_dut > (ymin+prec_y) && y_dut < (ymax-prec_y) && ntracks == 1
     && npix > 0 && nback > 0 && amp[0] > 0. && amp[0] < satCut){
    pad_amp[pad_id]->Fill(amp[0]);
  }

  if((x_dut < (xmin-prec_x) || x_dut > (xmax+prec_x)) && (y_dut < (ymin-prec_y) || y_dut > (ymax+prec_y)) && ntracks == 1
     && npix > 0 && nback > 0 && amp[0] > 0. && amp[0] < satCut){
    pad_noise[pad_id]->Fill(amp[0]);
  }
}

inline void Unbinned::NoiseSelection(Int_t pad_id){

  AssignCuts(pad_id);

  if(ntracks==1 && npix>0 && nback>0 && (x_dut<(xmin-prec) || x_dut>(xmax+prec))
     && (y_dut<(ymin-prec) || y_dut>(ymax+prec)) && amp[0]>0. && amp[0]<satCut){
    amp_rrvar[pad_id]->setVal(amp[0]);
    noise_dataset[pad_id]->add(*amp_argset[pad_id]);
  }
}

inline void Unbinned::SignalSelection(Int_t pad_id){

  AssignCuts(pad_id);

  if(ntracks==1 && npix>0 && nback>0 && (x_dut>(xmin+prec) && x_dut<(xmax-prec))
     && (y_dut>(ymin+prec) && y_dut<(ymax-prec)) && amp[0]>0. && amp[0]<satCut){

    amp_rrvar[pad_id]->setVal(amp[0]);
    mip_dataset[pad_id]->add(*amp_argset[pad_id]);
    pad_amp[pad_id]->Fill(amp[0]);
  }
}

inline void Unbinned::ProcessEntries(Int_t pad_id, std::function<void(Int_t)> selection){

  tch = new TChain("pulse");
  tch->Add(Form("%s/%s_tree_pad%i_reduced.root",chainPath.Data(),tag.Data(),pad_id));
  InitBranches();

  Int_t Nentry = tch->GetEntries();

  cout << "Processing " << Nentry << " entries in pad " << pad_id << "." << endl;

  if(isDebug)
    Nentry = 1000;

  for(int e = 0; e < Nentry; e++){
    if (e % 1000 == 0) {
      fprintf(stdout, "\r  Processed events: %8d of %8d ", e, Nentry);
    }
    fflush(stdout);
    tch->GetEntry(e);

    selection(pad_id);
  }
  cout << endl;

  delete tch;
}

inline void Unbinned::ProcessPadInfo(Int_t pad_id){

  PrepareCuts();

  if(vecsize == 0)
    PrepVectors();

  Float_t x1 = frame_margins[0];
  Float_t x2 = frame_margins[1];
  Float_t y1 = frame_margins[2];
  Float_t y2 = frame_margins[3];

  name = Form("xy_pad%i",pad_id);
  pad2D[pad_id] = new TH2F(name,name,180,x1,x2,60,y1,y2);

  name = Form("xy_all_pad%i",pad_id);
  pad2D_all[pad_id] = new TH2F(name,name,180,x1,x2,60,y1,y2);

  name = Form("amp_pad%i",pad_id);
  pad_amp[pad_id] = new TH1F(name,name,100,0.,satCut);

  name = Form("noise_pad%i",pad_id);
  pad_noise[pad_id] = new TH1F(name,name,100,0.,satCut);

  ProcessEntries(pad_id,[this](Int_t pad_id){PadInfoSelection(pad_id);});

}

inline void Unbinned::SavePadInfo(){

  PrepVectors();

  outDir = Form("../output/%s/",tag.Data());
  if(gSystem->OpenDirectory(outDir) == 0)
    gSystem->mkdir(outDir);

  name = Form("pad_summary_%s",tag.Data());
  TString outFileName = Form("%s%s.root",outDir.Data(),name.Data());
  TFile fOut(outFileName,"recreate");

  for(int ipad = 1; ipad < NUM_PADS+1; ipad++){
    ProcessPadInfo(ipad);
    pad2D[ipad]->SetDrawOption("COLZ");
    pad2D[ipad]->Write();
    pad2D_all[ipad]->Write();
    pad_amp[ipad]->Write();
    pad_noise[ipad]->Write();
  }
  fOut.Close();
}

inline void Unbinned::ProcessPadNoise(Int_t pad_id){

  PrepareCuts();
  if(vecsize == 0)
    PrepVectors();

  name = Form("rrvar_pad%i",pad_id);
  amp_rrvar[pad_id] = new RooRealVar(name,"Signal Amplitude [mV]",0.,satCut);

  name = Form("fitNoise_pad%i",pad_id);
  amp_rrvar[pad_id]->setRange(name, xl,xh);

  amp_argset[pad_id] = new RooArgSet(*amp_rrvar[pad_id]);

  name = Form("noise_ds_pad%i",pad_id);
  noise_dataset[pad_id] = new RooDataSet(name,name,*amp_argset[pad_id]);

  ProcessEntries(pad_id,[this](Int_t pad_id){NoiseSelection(pad_id);});
  FitPadNoise(pad_id);
}

inline void Unbinned::FitPadNoise(Int_t pad_id){

  cout << endl;
  cout << "Fitting noise distribution on pad " << pad_id <<  "." << endl;

  name = Form("pad%i_noise",pad_id);
  ws[pad_id] = new RooWorkspace(name);

  name = Form("rkPDF_pad%i",pad_id);
  rkpdf[pad_id] = new RooKeysPdf(name,name,*amp_rrvar[pad_id],*noise_dataset[pad_id],RooKeysPdf::NoMirror);

  name = Form("Nn_pad%i",pad_id);
  Nn[pad_id] = new RooRealVar(name, "N_{n}", 8000, "events");
  Nn[pad_id]->setConstant(kFALSE);

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  rplot[pad_id] = amp_rrvar[pad_id]->frame(xl,xh,nbins);

  name = Form("fitNoise_pad%i",pad_id);
  rkpdf[pad_id]->fitTo(*noise_dataset[pad_id],Range(name));
  noise_dataset[pad_id]->plotOn(rplot[pad_id]);
  rkpdf[pad_id]->plotOn(rplot[pad_id]);

}

inline void Unbinned::SavePadNoise(){

  PrepVectors();

  outDir = "../ROOT/";
  name = Form("noise_fits_%s",tag.Data());
  TString outFileName = Form("%s%s.root",outDir.Data(),name.Data());
  TFile fOut(outFileName,"recreate");

  for(int ipad = 1; ipad < NUM_PADS+1; ipad++){
    ProcessPadNoise(ipad);
    ws[ipad]->import(*rkpdf[ipad]);
    ws[ipad]->Write();
    rplot[ipad]->Write();
  }
}

inline void Unbinned::ProcessPadSignal(Int_t pad_id){

  PrepareCuts();
  if(vecsize == 0)
    PrepVectors();

  name = Form("rrvar_pad%i",pad_id);
  amp_rrvar[pad_id] = new RooRealVar(name,"Signal Amplitude [mV]",0.,satCut);

  name = Form("fitSignal_pad%i",pad_id);
  amp_rrvar[pad_id]->setRange(name, xl,xh);

  amp_argset[pad_id] = new RooArgSet(*amp_rrvar[pad_id]);

  name = Form("signal_pad%i",pad_id);
  mip_dataset[pad_id] = new RooDataSet(name,name,*amp_argset[pad_id]);

  name = Form("amp_pad%i",pad_id);
  pad_amp[pad_id] = new TH1F(name,name,100,0.,satCut);

  ProcessEntries(pad_id,[this](Int_t pad_id){SignalSelection(pad_id);});
  FitPadSignal(pad_id);
}

inline void Unbinned::FitPadSignal(Int_t pad_id){

  Float_t max;
  
  name = Form("pad%i_noise", pad_id);
  wsIn[pad_id] = (RooWorkspace*)fIn->Get(name);
  name = Form("rkPDF_pad%i",pad_id);
  rkpdf[pad_id] = (RooKeysPdf*)wsIn[pad_id]->pdf(name);

  name = Form("Nn_pad%i",pad_id);
  Nn[pad_id] = new RooRealVar(name, "N_{n}", 8000, "events");
  Nn[pad_id]->setConstant(kFALSE);

  //Get approximate mean                                                                                                         
  max = pad_amp[pad_id]->GetBinCenter(pad_amp[pad_id]->GetMaximumBin());

  //Landau fit                                                                                                                   
  name = Form("ml_pad%i",pad_id);
  ml[pad_id] = new RooRealVar(name,"MPV", max, -1000000., 1000000.);
  name = Form("sl_pad%i",pad_id);
  sl[pad_id] = new RooRealVar(name,"Landau Width",25., 0., 10000.);
  name = Form("landau_pad%i",pad_id);
  landau[pad_id] = new RooLandau(name,name,*amp_rrvar[pad_id],*ml[pad_id],*sl[pad_id]);

  //Gaussian fit                                                                                                                 
  name = Form("mg_pad%i",pad_id);
  mg[pad_id] = new RooRealVar(name,"Gauss Mean",0.0);
  name = Form("sg_pad%i",pad_id);
  sg[pad_id] = new RooRealVar(name,"Gauss Width",25., 0., 10000.);
  name = Form("gauss_pad%i",pad_id);
  gauss[pad_id] = new RooGaussian(name,name,*amp_rrvar[pad_id],*mg[pad_id],*sg[pad_id]);
  mg[pad_id]->setConstant( kTRUE );
  sg[pad_id]->setConstant( kFALSE );

  //LandauXGaussian convolution                                                                                                  
  name = Form("lxg_pad%i",pad_id);
  lxg[pad_id] = new RooFFTConvPdf(name, name, *amp_rrvar[pad_id], *landau[pad_id], *gauss[pad_id]);

  name = Form("Ns_pad%i",pad_id);
  Ns[pad_id] = new RooRealVar( name, "N_{s}", 8000, "events");
  Ns[pad_id]->setConstant(kFALSE);

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  
  name = Form("conv_gauss_landau_pad%i",pad_id);
  pdf_lxg[pad_id] = new RooAddPdf(name, name, RooArgList(*lxg[pad_id], *rkpdf[pad_id]), RooArgList(*Ns[pad_id],*Nn[pad_id]));

  name = Form("fitSignal_pad%i",pad_id);
  pdf_lxg[pad_id]->fitTo(*mip_dataset[pad_id], RooFit::PrintLevel(-1), RooFit::Strategy(0),RooFit::Extended( kTRUE ), RooFit::Range(name) );

  cout << endl;
  cout << "Pad " << pad_id << " fit parameters: " << endl;
  cout << "  MPV: " << ml[pad_id]->getValV() << " +/- " << ml[pad_id]->getError() << endl;
  cout << "  Landau Width: " << sl[pad_id]->getValV() << endl;
  cout << "  Gauss Width: " << sg[pad_id]->getValV() << endl;
  cout << endl;

  rplot[pad_id] = amp_rrvar[pad_id]->frame(xl,xh,nbins);
  mip_dataset[pad_id]->plotOn(rplot[pad_id]);
  pdf_lxg[pad_id]->plotOn(rplot[pad_id]);
  pdf_lxg[pad_id]->plotOn(rplot[pad_id],Components(*rkpdf[pad_id]),LineColor(kRed));
  pdf_lxg[pad_id]->plotOn(rplot[pad_id],Components(*lxg[pad_id]),LineColor(kGreen));
  pdf_lxg[pad_id]->paramOn(rplot[pad_id],Format("NE", AutoPrecision(2)), Layout(0.65,0.89,0.89));
  rplot[pad_id]->getAttText()->SetTextSize(0.02);
  name = Form("fitSignal_pad%i",pad_id);
  rplot[pad_id]->SetTitle(name);
}

inline void Unbinned::SavePadSignal(){

  PrepVectors();

  inDir ="../ROOT/";
  name = Form("noise_fits_%s.root",tag.Data());
  fIn = new TFile(inDir+name,"READ");

  for(int ipad = 1; ipad < NUM_PADS+1; ipad++)
    ProcessPadSignal(ipad);

  outDir = Form("../output/%s/",tag.Data());
  if(gSystem->OpenDirectory(outDir) == 0)
    gSystem->mkdir(outDir);
  name = Form("signal_fits_%s",tag.Data());
  TString outFileName = Form("%s%s.root",outDir.Data(),name.Data());
  fOut = new TFile(outFileName,"recreate");

  for(int ipad = 1; ipad < NUM_PADS+1; ipad++){
    name = Form("fitSignal_pad%i",ipad);
    rplot[ipad]->Write(name);    
  }

  cout << endl;

  fIn->Close();
  fOut->Close();

  delete fIn;
  delete fOut;
}

inline void Unbinned::PrepVectors(){

  if(xh < 0)
    xh = satCut;

  if(vecsize != 0)
    ClearVectors();

  for(int ipad = 0; ipad < NUM_PADS+1; ipad++){

    //Hist object vectors
    pad_amp.push_back(new TH1F());
    pad_noise.push_back(new TH1F());
    pad2D.push_back(new TH2F());
    pad2D_all.push_back(new TH2F());

    //Roofit object vectors
    amp_rrvar.push_back(new RooRealVar());
    amp_argset.push_back(new RooArgSet());
    noise_dataset.push_back(new RooDataSet());
    mip_dataset.push_back(new RooDataSet());
    ml.push_back(new RooRealVar());
    sl.push_back(new RooRealVar());
    mg.push_back(new RooRealVar());
    sg.push_back(new RooRealVar());
    landau.push_back(new RooLandau());
    gauss.push_back(new RooGaussian());
    lxg.push_back(new RooFFTConvPdf());
    pdf_lxg.push_back(new RooAddPdf());
    rkpdf.push_back(new RooKeysPdf());
    Nn.push_back(new RooRealVar());
    Ns.push_back(new RooRealVar());
    rplot.push_back(new RooPlot());
    ws.push_back(new RooWorkspace());
    wsIn.push_back(new RooWorkspace());

    vecsize++;
  }
}

inline void Unbinned::ClearVectors(){

  //histograms
  pad_amp.clear();
  pad_noise.clear();
  pad2D.clear();
  pad2D_all.clear();

  //RooFit
  amp_rrvar.clear();
  amp_argset.clear();
  noise_dataset.clear();
  mip_dataset.clear();
  ml.clear();
  sl.clear();
  mg.clear();
  sg.clear();
  landau.clear();
  gauss.clear();
  lxg.clear();
  pdf_lxg.clear();
  rkpdf.clear();
  Nn.clear();
  Ns.clear();
  rplot.clear();
  ws.clear();

  vecsize = 0;
}

inline void Unbinned::PrepareCuts(){

  for(int ipad = 0; ipad < NUM_PADS+1; ipad++){
    geom_cuts[ipad][0] = xcuts[geom_map[ipad][0]];
    geom_cuts[ipad][1] = xcuts[geom_map[ipad][1]];
    geom_cuts[ipad][2] = ycuts[geom_map[ipad][2]];
    geom_cuts[ipad][3] = ycuts[geom_map[ipad][3]];
  }
}

inline void Unbinned::AssignCuts(Int_t pad_id){

  xmin = geom_cuts[pad_id][0];
  xmax = geom_cuts[pad_id][1];
  ymin = geom_cuts[pad_id][2];
  ymax = geom_cuts[pad_id][3];

}

inline void Unbinned::PrintGeomMap(){

  for(int ipad = 1; ipad < NUM_PADS+1; ipad++){
    cout << "Indices in pad " << ipad << ": " << endl;
    cout << geom_map[ipad][0] << " <  x < " << geom_map[ipad][1] << endl;
    cout << geom_map[ipad][2] << " <  y < " << geom_map[ipad][3] << endl;
  }

  cout << endl;
  cout << "            HPK type 3.1 Pad Layout           " << endl;
  cout << "4-|---------|----------|----------|----------|" << endl;
  cout << "  |    7    |     5    |     4    |     2    |" << endl;
  cout << "3-|---------|----------|----------|----------|" << endl;
  cout << "  |    8    |     6    |     3    |     1    |" << endl;
  cout << "2-|---------|----------|----------|----------|" << endl;
  cout << "  |    9    |    11    |    14    |    16    |" << endl;
  cout << "1-|---------|----------|----------|----------|" << endl;
  cout << "  |   10    |    12    |    13    |    15    |" << endl;
  cout << "0-|---------|----------|----------|----------|" << endl;
  cout << "  0         1          2          3          4" << endl;
  cout << endl;
}

inline void Unbinned::InitBranches(){

  tch->SetBranchStatus("*", 0);
  //analog info                                                                                                                  
  tch->SetBranchStatus("amp", 1); tch->SetBranchAddress("amp", &amp);
  //Track info                                                                                                                   
  tch->SetBranchStatus("ntracks", 1); tch->SetBranchAddress("ntracks", &ntracks);
  tch->SetBranchStatus("nback", 1); tch->SetBranchAddress("nback", &nback);
  tch->SetBranchStatus("npix", 1); tch->SetBranchAddress("npix", &npix);
  tch->SetBranchStatus("x_dut", 1); tch->SetBranchAddress("x_dut", &x_dut);
  tch->SetBranchStatus("y_dut", 1); tch->SetBranchAddress("y_dut", &y_dut);
}
