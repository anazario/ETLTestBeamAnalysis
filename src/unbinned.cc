#include "../include/unbinned.h"
using namespace std;

inline void Unbinned::PrintGeomMap(){

  for(int i = 1; i < npad+1; i++){
    cout << "Indices in pad " << i << ": " << endl;
    cout << "   xmin: " << geom_map[i][0] << endl;
    cout << "   xmax: " << geom_map[i][1] << endl;
    cout << "   ymin: " << geom_map[i][2] << endl;
    cout << "   ymax: " << geom_map[i][3] << endl;
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

inline void Unbinned::SavePadNoise(Float_t xl, Float_t xh, int nbins){

  float xmin = -999.;
  float xmax = -999.;
  float ymin = -999.;
  float ymax = -999.;

  float xt,yt,xp,yp;

  double costheta = TMath::Cos(angle);
  double sintheta = TMath::Sin(angle);

  //Create roodatasets for each pad (pads start at 1)
  for(int ipad = 0; ipad < npad+1; ipad++){

    name = Form("rrvar_pad%i",ipad);
    amp_rrvar.push_back(new RooRealVar(name,"Signal Amplitude [mV]",0.,satCut));
    name = Form("fitNoise_pad%i",ipad);
    amp_rrvar[ipad]->setRange(name, xl,xh);

    noise_plot.push_back(new RooPlot()); 

    amp_argset.push_back(new RooArgSet(*amp_rrvar[ipad]));

    name = Form("noise_ds_pad%i",ipad);
    noise_dataset.push_back(new RooDataSet(name,name,*amp_argset[ipad]));

  }
  
  //Fill datasets from trees
  PrepareTree();
  InitBranches();
  
  int Nentry = tch->GetEntries(); 
  int Nentry2 = 10000;
  
  for(int e = 0; e < Nentry; e++){
    if (e % 1000 == 0) {
      fprintf(stdout, "\r  Processed events: %8d of %8d ", e, Nentry);
      }
    fflush(stdout);
    tch->GetEntry(e);
    
    xt = x_dut[2];
    yt = y_dut[2];
    
    xp = xt*costheta-yt*sintheta;
    yp = yt*costheta+xt*sintheta;
    
    for(int ilgad = 0; ilgad < nlgads; ilgad++){
      
      pad_id = pads->at(ilgad);
      xmin = geom_cuts[pad_id][0];
      xmax = geom_cuts[pad_id][1];
      ymin = geom_cuts[pad_id][2];
      ymax = geom_cuts[pad_id][3];
      
      if(ntracks==1 && npix>0 && nback>0 && (xp<(xmin-prec) || xp>(xmax+prec))
	 && (yp<(ymin-prec) || yp>(ymax+prec)) && amp[ilgad]>0. && amp[ilgad]<satCut){
	amp_rrvar[pad_id]->setVal(amp[ilgad]);
	noise_dataset[pad_id]->add(*amp_argset[pad_id]);
      }
    }
  }  

  outDir = "../ROOT/";
  name = Form("noise_fits_%s",tag.Data());
  TString outFileName = Form("%s%s.root",outDir.Data(),name.Data());
  TFile fOut(outFileName,"recreate");

  for (int ipad = 0; ipad < npad+1; ipad++){

    name = Form("pad%i_noise",ipad);
    ws.push_back(new RooWorkspace(name));

    name = Form("rkPDF_pad%i",ipad);
    rkpdf.push_back(new RooKeysPdf(name,name,*amp_rrvar[ipad],*noise_dataset[ipad],RooKeysPdf::NoMirror));

    name = Form("Nn_pad%i",ipad);
    Nn.push_back(new RooRealVar(name, "N_{n}", 8000, "events"));
    Nn[ipad]->setConstant(kFALSE);

    if (ipad > 0){      

      cout << endl;
      cout << "Fitting noise on pad " << ipad <<  "." << endl;

      RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
      
      noise_plot[ipad] = amp_rrvar[ipad]->frame(xl,xh,nbins);
      
      name = Form("fitNoise_pad%i",ipad);
      rkpdf[ipad]->fitTo(*noise_dataset[ipad],Range(name));
      noise_dataset[ipad]->plotOn(noise_plot[ipad]);
      rkpdf[ipad]->plotOn(noise_plot[ipad]);

      ws[ipad]->import(*rkpdf[ipad]);
      ws[ipad]->Write();
      noise_plot[ipad]->Write();
    }
  } 

  fOut.Close();
}

inline void Unbinned::Save2DPads(Float_t xl, Float_t xh, Float_t yl, Float_t yh){

  Float_t xmin = -999.;
  Float_t xmax = -999.;
  Float_t ymin = -999.;
  Float_t ymax = -999.;

  Float_t xt,yt,xp,yp;

  double costheta = TMath::Cos(angle);
  double sintheta = TMath::Sin(angle);

  Float_t prec_x = prec*2.;
  Float_t prec_y = prec;

  for(int ipad = 0; ipad<npad+1; ipad++){

    name = Form("xy_pad%i",ipad);
    pad2D.push_back(new TH2F(name,name,180,xl,xh,60,yl,yh));

    name = Form("xy_all_pad%i",ipad);
    pad2D_all.push_back(new TH2F(name,name,180,xl,xh,60,yl,yh));

    name = Form("amp_pad%i",ipad);
    pad_amp.push_back(new TH1F(name,name,100,0.,satCut));

    name = Form("noise_pad%i",ipad);
    pad_noise.push_back(new TH1F(name,name,100,0.,satCut));
  }

  PrepareTree();
  InitBranches();

  int Nentry = tch->GetEntries();
  int Nentry2 = 10000;

  for(int e = 0; e < Nentry; e++){
    if (e % 1000 == 0) {
      fprintf(stdout, "\r  Processed events: %8d of %8d ", e, Nentry);
    }
    fflush(stdout);
    tch->GetEntry(e);

    xt = x_dut[2];
    yt = y_dut[2];

    xp = xt*costheta-yt*sintheta;
    yp = yt*costheta+xt*sintheta;

    for(int ilgad = 0; ilgad < nlgads; ilgad++){

      pad_id = pads->at(ilgad);

      xmin = geom_cuts[pad_id][0];
      xmax = geom_cuts[pad_id][1];
      ymin = geom_cuts[pad_id][2];
      ymax = geom_cuts[pad_id][3];

      if(xp > (xmin-prec_x) && xp < (xmax+prec_x) && yp > (ymin-prec_y) && yp < (ymax+prec_y) && ntracks == 1 
	 && npix > 0 && nback > 0 && amp[ilgad] > 0. && amp[ilgad] < satCut){
        
        pad2D[pad_id]->Fill(xp,yp);
      }

      if(ntracks == 1 && npix > 0 && nback > 0 && amp[ilgad] > 20. && amp[ilgad] < satCut)
	pad2D_all[pad_id]->Fill(xp,yp);

      if(xp > (xmin+prec_x) && xp < (xmax-prec_x) && yp > (ymin+prec_y) && yp < (ymax-prec_y) && ntracks == 1
         && npix > 0 && nback > 0 && amp[ilgad] > 0. && amp[ilgad] < satCut){
        pad_amp[pad_id]->Fill(amp[ilgad]);
      }

      if((xp < (xmin-prec_x) || xp > (xmax+prec_x)) && (yp < (ymin-prec_y) || yp > (ymax+prec_y)) && ntracks == 1
         && npix > 0 && nback > 0 && amp[ilgad] > 0. && amp[ilgad] < satCut){
        pad_noise[pad_id]->Fill(amp[ilgad]);
      }
    }

  }

  cout << endl;

  outDir = Form("../output/%s/",tag.Data());
  if(gSystem->OpenDirectory(outDir) == 0)
    gSystem->mkdir(outDir);
  name = Form("pad_summary_%s",tag.Data());
  TString outFileName = Form("%s%s.root",outDir.Data(),name.Data());
  TFile fOut(outFileName,"recreate");

  for(int ipad = 1; ipad < npad+1; ipad++){
    pad2D[ipad]->SetDrawOption("COLZ");
    pad2D[ipad]->Write();
    pad2D_all[ipad]->Write();
    pad_amp[ipad]->Write();
    pad_noise[ipad]->Write();
  }
  fOut.Close();

}

inline void Unbinned::SavePadSignalFit(Float_t xl, Float_t xh, int nbins){

  float xmin = -999.;
  float xmax = -999.;
  float ymin = -999.;
  float ymax = -999.;

  float xt,yt,xp,yp;

  double costheta = TMath::Cos(angle);
  double sintheta = TMath::Sin(angle);

  Float_t prec_x = prec*2.;
  Float_t prec_y = prec;

  vector<RooWorkspace*> wsIn;

  for(int ipad = 0; ipad < npad+1; ipad++){

    wsIn.push_back(new RooWorkspace());
    rkpdf.push_back(new RooKeysPdf);

    name = Form("rrvar_pad%i",ipad);
    amp_rrvar.push_back(new RooRealVar(name,"Signal Amplitude [mV]",0.,satCut));

    name = Form("fitSignal_pad%i",ipad);
    amp_rrvar[ipad]->setRange(name, xl,xh);

    mip_plot.push_back(new RooPlot());

    amp_argset.push_back(new RooArgSet(*amp_rrvar[ipad]));

    name = Form("signal_pad%i",ipad);
    mip_dataset.push_back(new RooDataSet(name,name,*amp_argset[ipad]));

    name = Form("amp_pad%i",ipad);
    pad_amp.push_back(new TH1F(name,name,100,0.,satCut));
  }

  PrepareTree();
  InitBranches();

  int Nentry = tch->GetEntries();
  int Nentry2 = 10000;

  for(int e = 0; e < Nentry; e++){
    if (e % 1000 == 0) {
      fprintf(stdout, "\r  Processed events: %8d of %8d ", e, Nentry);
    }
    fflush(stdout);
    tch->GetEntry(e);

    xt = x_dut[2];
    yt = y_dut[2];

    xp = xt*costheta-yt*sintheta;
    yp = yt*costheta+xt*sintheta;

    for(int ilgad = 0; ilgad < nlgads; ilgad++){

      pad_id = pads->at(ilgad);

      xmin = geom_cuts[pad_id][0];
      xmax = geom_cuts[pad_id][1];
      ymin = geom_cuts[pad_id][2];
      ymax = geom_cuts[pad_id][3];

      if(ntracks==1 && npix>0 && nback>0 && (xp>(xmin+prec) && xp<(xmax-prec))
         && (yp>(ymin+prec) && yp<(ymax-prec)) && amp[ilgad]>0. && amp[ilgad]<satCut){

        amp_rrvar[pad_id]->setVal(amp[ilgad]);
        mip_dataset[pad_id]->add(*amp_argset[pad_id]);
	pad_amp[pad_id]->Fill(amp[ilgad]);
      }
    }
  }
  cout << endl;

  inDir ="../ROOT/";
  name = Form("noise_fits_%s.root",tag.Data());
  TFile fIn(inDir+name,"READ");

  outDir = Form("../output/%s/",tag.Data());
  if(gSystem->OpenDirectory(outDir) == 0)
    gSystem->mkdir(outDir);
  name = Form("signal_fits_%s",tag.Data());
  TString outFileName = Form("%s%s.root",outDir.Data(),name.Data());
  TFile fOut(outFileName,"recreate");

  Float_t max;

  for(int ipad = 0; ipad < npad+1; ipad++){
    
    //Get noise fit from root file
    if(ipad > 0){
      name = Form("pad%i_noise",ipad);
      wsIn[ipad] = (RooWorkspace*)fIn.Get(name);
      name = Form("rkPDF_pad%i",ipad);
      rkpdf[ipad] = (RooKeysPdf*)wsIn[ipad]->pdf(name);
    }

    name = Form("Nn_pad%i",ipad);
    Nn.push_back(new RooRealVar( "Nn", "N_{n}", 8000, "events"));
    Nn[ipad]->setConstant(kFALSE);

    //Get approximate mean
    max = pad_amp[ipad]->GetBinCenter(pad_amp[ipad]->GetMaximumBin());

    //Landau fit
    name = Form("ml_pad%i",ipad);
    ml.push_back(new RooRealVar(name,name, max, -1000000., 1000000.));
    name = Form("sl_pad%i",ipad);
    sl.push_back(new RooRealVar(name,name,25., 0., 10000.));
    name = Form("landau_pad%i",ipad);
    landau.push_back(new RooLandau(name,name,*amp_rrvar[ipad],*ml[ipad],*sl[ipad]));

    //Gaussian fit
    name = Form("mg_pad%i",ipad);
    mg.push_back(new RooRealVar(name,name,0.0));
    name = Form("sg_pad%i",ipad);
    sg.push_back(new RooRealVar(name,name,25., 0., 10000.));
    name = Form("gauss_pad%i",ipad);
    gauss.push_back(new RooGaussian(name,name,*amp_rrvar[ipad],*mg[ipad],*sg[ipad]));
    mg[ipad]->setConstant( kTRUE );
    sg[ipad]->setConstant( kFALSE );

    //LandauXGaussian convolution
    name = Form("lxg_pad%i",ipad);
    lxg.push_back(new RooFFTConvPdf(name, name, *amp_rrvar[ipad], *landau[ipad], *gauss[ipad]));

    name = Form("Ns_pad%i",ipad);
    Ns.push_back(new RooRealVar( "Ns", "N_{s}", 8000, "events"));
    Ns[ipad]->setConstant(kFALSE);

    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    name = Form("conv_gauss_landau_pad%i",ipad);
    pdf_lxg.push_back(new RooAddPdf(name, name, RooArgList(*lxg[ipad], *rkpdf[ipad]), RooArgList(*Ns[ipad],*Nn[ipad])));

    //
    if (ipad > 0){
      name = Form("fitSignal_pad%i",ipad);
      pdf_lxg[ipad]->fitTo(*mip_dataset[ipad], RooFit::PrintLevel(-1), RooFit::Strategy(0),RooFit::Extended( kTRUE ), RooFit::Range(name) );

      cout << endl;
      cout << "Pad " << ipad << " fit parameters: " << endl; 
      cout << "  MPV: " << ml[ipad]->getValV() << " +/- " << ml[ipad]->getError() << endl;
      cout << "  Landau Width: " << sl[ipad]->getValV() << endl;
      cout << "  Gauss Width: " << sg[ipad]->getValV() << endl;

      mip_plot[ipad] = amp_rrvar[ipad]->frame(xl,xh,nbins);
      mip_dataset[ipad]->plotOn(mip_plot[ipad]);
      pdf_lxg[ipad]->plotOn(mip_plot[ipad]);
      pdf_lxg[ipad]->plotOn(mip_plot[ipad],Components(*rkpdf[ipad]),LineColor(kRed));
      pdf_lxg[ipad]->plotOn(mip_plot[ipad],Components(*lxg[ipad]),LineColor(kGreen));

      name = Form("fitSignal_pad%i",ipad);
      mip_plot[ipad]->Write(name);
      //fOut.WriteTObject(mip_plot[ipad],name);
    }
  } 
  cout << endl;

  fOut.Close();
  fIn.Close();
}

inline void Unbinned::PrepareCuts(){

  Float_t xmin, xmax, ymin, ymax; 
  
  for(int ipad = 0; ipad < npad+1; ipad++){
    
    xmin = xcuts[geom_map[ipad][0]];
    xmax = xcuts[geom_map[ipad][1]];
    ymin = ycuts[geom_map[ipad][2]];
    ymax = ycuts[geom_map[ipad][3]];

    geom_cuts[ipad][0] = xmin;
    geom_cuts[ipad][1] = xmax;
    geom_cuts[ipad][2] = ymin;
    geom_cuts[ipad][3] = ymax;
  }  
}

inline void Unbinned::PrepareTree(){

  tch = new TChain("pulse");

  int run_lsts = run_start->size();

  for(int i_runrange = 0;i_runrange < run_lsts;i_runrange++){

    cout << "Preparing runs: " << run_start->at(i_runrange) << " through " << run_end->at(i_runrange) << "." << endl; 

    for(int irun = run_start->at(i_runrange); irun <= run_end->at(i_runrange);irun++)
      tch->Add(Form("%s/run_scope%i_info.root",chainPath.Data(),irun));
  }
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
    //Run conf info
    tch->SetBranchStatus("run", 1); tch->SetBranchAddress("run", &run);
    tch->SetBranchStatus("pads", 1); tch->SetBranchAddress("pads", &pads);
}
