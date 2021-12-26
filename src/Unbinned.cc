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

    pad_amp[pad_id]->Fill(amp[0]);

    amp_rrvar[pad_id]->setVal(amp[0]);
    x_rrvar[pad_id]->setVal(x_dut);
    y_rrvar[pad_id]->setVal(y_dut);
    mip_dataset[pad_id]->add(*amp_argset[pad_id]);

    h3[pad_id]->Fill(x_dut,y_dut,amp[0]);
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

  amp_argset[pad_id] = new RooArgSet(*amp_rrvar[pad_id]);

  name = Form("noise_ds_pad%i",pad_id);
  noise_dataset[pad_id] = new RooDataSet(name,name,*amp_argset[pad_id]);

  ProcessEntries(pad_id,[this](Int_t pad_id){NoiseSelection(pad_id);});
}

inline void Unbinned::FitPadNoise(Int_t pad_id){

  name = Form("pad%i_noise",pad_id);
  ws[pad_id] = new RooWorkspace(name);

  name = Form("pad%i",pad_id);
  sf[pad_id]->SetFitRange(xl,xh,name);
  sf[pad_id]->NoiseFit(noise_dataset[pad_id], amp_rrvar[pad_id], name);
  sf[pad_id]->PlotNoise(xl, xh, nbins, name);
  rplot[pad_id] = sf[pad_id]->GetPlot();
}

inline void Unbinned::SavePadNoise(){

  PrepVectors();

  outDir = "../ROOT/";
  name = Form("noise_fits_%s",tag.Data());
  TString outFileName = Form("%s%s.root",outDir.Data(),name.Data());
  TFile fOut(outFileName,"recreate");

  for(int ipad = 1; ipad < NUM_PADS+1; ipad++){
    ProcessPadNoise(ipad);
    FitPadNoise(ipad);
    ws[ipad]->import(*rkpdf[ipad]);
    ws[ipad]->Write();
    rplot[ipad]->Write();
  }
}

inline void Unbinned::ProcessPadSignal(Int_t pad_id){

  Float_t xlo,xhi,ylo,yhi;

  if(nbinsX < 0 || nbinsY < 0){
    cout << "Using default binning." << endl;
    nbinsX = 120;
    nbinsY = 60;
    nbinsZ = nbins;
  }

  if(frame_margins.size() == 0 || frame_margins.size() > 4){
    cout << "Enter valid margin values for the histogram!" << endl;
    exit(0);
  }
  else{
    xlo = frame_margins[0];
    xhi = frame_margins[1];
    ylo = frame_margins[2];
    yhi = frame_margins[3];
  }

  PrepareCuts();
  if(vecsize == 0)
    PrepVectors();

  name = Form("rrvar_pad%i",pad_id);
  amp_rrvar[pad_id] = new RooRealVar(name,"Signal Amplitude [mV]",0.,satCut);
  name = Form("x_rrvar_pad%i",pad_id);
  x_rrvar[pad_id] = new RooRealVar(name,"x [mm]",xlo,xhi);
  name = Form("y_rrvar_pad%i",pad_id);
  y_rrvar[pad_id] = new RooRealVar(name,"y [mm]",ylo,yhi);
  
  amp_argset[pad_id] = new RooArgSet(*amp_rrvar[pad_id],*x_rrvar[pad_id],*y_rrvar[pad_id]);

  name = Form("signal_pad%i",pad_id);
  mip_dataset[pad_id] = new RooDataSet(name,name,*amp_argset[pad_id]);

  name = Form("amp_pad%i",pad_id);
  pad_amp[pad_id] = new TH1F(name,name,100,0.,satCut);

  name = Form("h3_pad%i",pad_id);
  h3[pad_id] = new TH3F(name,name,nbinsX,xlo,xhi,nbinsY,ylo,yhi,nbinsZ,xl,xh);

  ProcessEntries(pad_id,[this](Int_t pad_id){SignalSelection(pad_id);});
}

inline void Unbinned::FitPadSignal(Int_t pad_id){

  //Get approximate mean
  Float_t max = pad_amp[pad_id]->GetBinCenter(pad_amp[pad_id]->GetMaximumBin());

  name = Form("pad%i_noise", pad_id);
  wsIn[pad_id] = (RooWorkspace*)fIn->Get(name);
  name = Form("rkPDF_pad%i",pad_id);
  rkpdf[pad_id] = (RooKeysPdf*)wsIn[pad_id]->pdf(name);

  name = Form("pad%i",pad_id);
  sf[pad_id]->SetFitRange(xl,xh,name);
  sf[pad_id]->AmplitudeFit(mip_dataset[pad_id], amp_rrvar[pad_id], rkpdf[pad_id], max, name);
  sf[pad_id]->PlotSignal(xl, xh, nbins, name);
  rplot[pad_id] = sf[pad_id]->GetPlot();
}

inline void Unbinned::SavePadSignal(){

  PrepVectors();
  
  inDir ="../ROOT/";
  name = Form("noise_fits_%s.root",tag.Data());
  fIn = new TFile(inDir+name,"READ");
  
  for(int ipad = 1; ipad < NUM_PADS+1; ipad++){
    ProcessPadSignal(ipad);
    FitPadSignal(ipad);
  }

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

  fIn->Close();
  fOut->Close();

  delete fIn;
  delete fOut;
}

inline std::pair<RooDataSet*,RooRealVar*> Unbinned::GetDataSet(Int_t pad_id, Float_t &mean){
  
  //ProcessPadNoise(pad_id); 
  ProcessPadSignal(pad_id); 
  mean = pad_amp[pad_id]->GetBinCenter(pad_amp[pad_id]->GetMaximumBin());
  return std::make_pair(mip_dataset[pad_id],amp_rrvar[pad_id]); 
  //return std::make_pair(noise_dataset[pad_id],amp_rrvar[pad_id]);
}

inline void Unbinned::PrepVectors(){

  if(xh < 0)
    xh = satCut;

  if(vecsize != 0)
    ClearVectors();

  for(int ipad = 0; ipad < NUM_PADS+1; ipad++){

    h3.push_back(new TH3F());

    //Hist object vectors
    pad_amp.push_back(new TH1F());
    pad_noise.push_back(new TH1F());
    pad2D.push_back(new TH2F());
    pad2D_all.push_back(new TH2F());

    //Roofit object vectors
    amp_rrvar.push_back(new RooRealVar());
    x_rrvar.push_back(new RooRealVar());
    y_rrvar.push_back(new RooRealVar());
    amp_argset.push_back(new RooArgSet());
    noise_dataset.push_back(new RooDataSet());
    mip_dataset.push_back(new RooDataSet());
    rkpdf.push_back(new RooKeysPdf());
    rplot.push_back(new RooPlot());
    ws.push_back(new RooWorkspace());
    wsIn.push_back(new RooWorkspace());

    sf.push_back(new SignalFit());

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
  rkpdf.clear();
  rplot.clear();
  ws.clear();
  wsIn.clear();
  sf.clear();

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

inline void Unbinned::MakeMap(){
  //Int_t pad_id = 2;

  Int_t tmp_entries;

  Float_t x1, x2, y1, y2;
  Float_t MPV = -999.;
  Float_t MPVerr = -999.;

  SignalFit *sf;

  name = "h2";
  TH2F *h2 = new TH2F(name,name,nbinsX,4.,17.,nbinsY,30.6,35.6);

  PrepVectors();
  inDir ="../ROOT/";
  name = Form("noise_fits_%s.root",tag.Data());
  fIn = new TFile(inDir+name,"READ");

  outDir = Form("../output/%s/",tag.Data());
  if(gSystem->OpenDirectory(outDir) == 0)
    gSystem->mkdir(outDir);
  name = Form("MPV_fits_%s",tag.Data());
  TString outFileName = Form("%s%s.root",outDir.Data(),name.Data());

  TFile *outfile = new TFile(outFileName,"recreate");

  TDirectory *folder;

  for(int ipad = 1; ipad < NUM_PADS+1; ipad++){
    ProcessPadSignal(ipad);

    folder = outfile->mkdir(Form("pad%i",ipad));
    folder->cd();

    //sf = new SignalFit();

    Float_t max = pad_amp[ipad]->GetBinCenter(pad_amp[ipad]->GetMaximumBin());

    name = Form("pad%i_noise", ipad);
    wsIn[ipad] = (RooWorkspace*)fIn->Get(name);
    name = Form("rkPDF_pad%i",ipad);
    rkpdf[ipad] = (RooKeysPdf*)wsIn[ipad]->pdf(name);
    
    RooAbsData* tmp_data;
    
    for(int xi = 1; xi < nbinsX+1; xi++){
      x1 = h2->GetXaxis()->GetBinLowEdge(xi);
      x2 = h2->GetXaxis()->GetBinUpEdge(xi);
      for(int yi = 1; yi < nbinsY+1; yi++){
	y1 = h2->GetYaxis()->GetBinLowEdge(yi);
	y2 = h2->GetYaxis()->GetBinUpEdge(yi);

	sf = new SignalFit();
	tmp_data = ReduceDataSet(mip_dataset[ipad],x1,x2,y1,y2,ipad);
	tmp_entries = tmp_data->numEntries();
	if(tmp_entries > 20 && (x1 >= xmin && x2 <= xmax && y1 >= ymin && y2 <= ymax)){
	  sf->SetFitRange(xl,xh,Form("pad%i",ipad));
	  sf->SetFitNorm(tmp_entries,tmp_entries);
	  sf->AmplitudeFit((RooDataSet*)tmp_data,amp_rrvar[ipad],rkpdf[ipad],max,Form("pad%i",ipad));
	  sf->PlotSignal(xl,xh,nbins,Form("pad%i_binx%f_biny%f_ent%i",ipad,x1,y2,tmp_entries));
	  MPV = sf->GetMPV();
	  MPVerr = sf->GetMPVerr();
	  sf->GetPlot()->Write(Form("pad%i_binx%f_biny%f_ent%i",ipad,x1,y2,tmp_entries));

	  if (MPV > 0. && MPV < 2*max && MPVerr < 3.)
	    h2->SetBinContent(xi,yi,MPV);
	}

	delete sf;	
      }
    }
  }
  gROOT->SetBatch(kTRUE);
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = (TCanvas*) new TCanvas("can","can",1100,500);

  //can->GetXaxis()->SetTitle

  h2->GetZaxis()->SetRangeUser(13.,25.);
  h2->GetXaxis()->SetTitle("x [mm]");
  h2->GetYaxis()->SetTitle("y [mm]");
  h2->Draw("COLZ");
  CMSmark("HPK 3.1 8e14");
  can->SaveAs("map_test.pdf");

  fIn->Close();
}

inline RooAbsData* Unbinned::ReduceDataSet(RooDataSet* data, Float_t x1, Float_t x2, Float_t y1, Float_t y2, Int_t pad_id){

  TString cutstring = Form("x_rrvar_pad%i > %f && x_rrvar_pad%i < %f && y_rrvar_pad%i > %f && y_rrvar_pad%i < %f",pad_id,x1,pad_id,x2,pad_id,y1,pad_id,y2);
  RooAbsData* reducedData = (RooAbsData*)mip_dataset[pad_id]->reduce(*amp_argset[pad_id],cutstring);

  return reducedData;
}

inline void Unbinned::ReduceTest(Float_t x1, Float_t x2, Float_t y1, Float_t y2){

  Int_t pad_id = 1;

  PrepVectors();

  ProcessPadSignal(pad_id);
  Float_t max = pad_amp[pad_id]->GetBinCenter(pad_amp[pad_id]->GetMaximumBin());

  RooPlot* newplot;
  /*
  Float_t x1 = 13.53;
  Float_t x2 = 14.03;
  Float_t y1 = 33.02;
  Float_t y2 = 33.7;
  */
  name = Form("(x_rrvar_pad1 > %f && x_rrvar_pad1 < %f) && (y_rrvar_pad1 > %f && y_rrvar_pad1 < %f)",x1,x2,y1,y2);
  RooAbsData* reducedData = (RooAbsData*)mip_dataset[pad_id]->reduce(*amp_argset[pad_id],name);

  cout << "entries: " << reducedData->numEntries() << endl;
  SignalFit* sigfit = new SignalFit();

  inDir ="../ROOT/";
  name = Form("noise_fits_%s.root",tag.Data());
  fIn = new TFile(inDir+name,"READ");
  name = Form("pad%i_noise", pad_id);
  wsIn[pad_id] = (RooWorkspace*)fIn->Get(name);
  name = Form("rkPDF_pad%i",pad_id);
  rkpdf[pad_id] = (RooKeysPdf*)wsIn[pad_id]->pdf(name);

  sigfit->SetFitRange(xl,xh,"pad1");
  sigfit->AmplitudeFit((RooDataSet*)reducedData,amp_rrvar[pad_id],rkpdf[pad_id],max,"pad1");
  sigfit->PlotSignal(xl,xh,nbins,"pad1");
  newplot = sigfit->GetPlot();

  gROOT->SetBatch(kTRUE);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = (TCanvas*) new TCanvas("can","can",600,800);
  can->SetLogy();

  newplot->Draw();
  can->SaveAs("test.pdf");
  can->Close();

  delete can;
}

void Unbinned::CMSmark(TString plotTitle){
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();
  l.SetTextSize(0.035);
  l.SetTextFont(42);
  // l.DrawLatex(0.17,0.855,g_PlotTitle.c_str());
  l.DrawLatex(0.51,0.91,plotTitle);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(0.11,0.91,"#bf{CMS} #it{Preliminary}");
}
