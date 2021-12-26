#include "../include/SignalFit.h"

inline void SignalFit::SetFitNorm(Int_t norm_noise, Int_t norm_sig){
  noise_n = norm_noise;
  sig_n = norm_sig;
}

inline void SignalFit::SetFitRange(Float_t xmin, Float_t xmax, TString label){
 
  set_range = true;
  xl = xmin;
  xh = xmax;
}

inline void SignalFit::AmplitudeFit(RooDataSet* data, RooRealVar* amp,  RooKeysPdf* error_fit, Float_t mean, TString label){

  CleanUp();

  if(noise_n < 0 || sig_n < 0)
    SetFitNorm(8000,8000);

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  dataset = data;
  amp_rrvar = amp;
  rkpdf = error_fit;

  if(set_range){
    name = "fitSignal_"+label;
    amp_rrvar->setRange(name, xl,xh);
  }

  name = "Nn_"+label;
  Nn = new RooRealVar(name, "N_{n}", noise_n, "events");
  Nn->setConstant(kFALSE);

  //Landau fit
  name = "ml_"+label;
  ml = new RooRealVar(name,"MPV", mean, -1000000., 1000000.);
  name = "sl_"+label;
  sl = new RooRealVar(name,"Landau Width",25., 0., 10000.);
  name = "landau_"+label;
  landau = new RooLandau(name,name,*amp_rrvar,*ml,*sl);

  //Gaussian fit
  name = "mg_"+label;
  mg = new RooRealVar(name,"Gauss Mean",0.0);
  name = "sg_"+label;
  sg = new RooRealVar(name,"Gauss Width",25., 0., 10000.);
  name = "gauss_"+label;
  gauss = new RooGaussian(name,name,*amp_rrvar,*mg,*sg);
  mg->setConstant(kTRUE);
  sg->setConstant(kFALSE);

  //LandauXGaussian convolution
  name = "lxg_"+label;
  lxg = new RooFFTConvPdf(name, name, *amp_rrvar, *landau, *gauss);

  name = "Ns_"+label;
  Ns = new RooRealVar( name, "N_{s}", sig_n, "events");
  Ns->setConstant(kFALSE);

  name = "conv_gauss_landau_"+label;
  pdf_lxg = new RooAddPdf(name, name, RooArgList(*lxg, *rkpdf), RooArgList(*Ns,*Nn));

  if(set_range){
    name = "fitSignal_"+label;
    pdf_lxg->fitTo(*dataset, RooFit::PrintLevel(-1), RooFit::Strategy(0),RooFit::Extended(kTRUE), RooFit::Range(name));
  }
  else 
    pdf_lxg->fitTo(*dataset, RooFit::PrintLevel(-1), RooFit::Strategy(0),RooFit::Extended(kTRUE));

  MPV =  ml->getValV();
  MPV_err = ml->getError();
  sigmaL = sl->getValV();
  sigmaG = sg->getValV();

  //PrintParameters(label);
}

inline void SignalFit::PlotSignal(Float_t xl, Float_t xh, Int_t nbins, TString label){

  rplot = amp_rrvar->frame(xl,xh,nbins);
  dataset->plotOn(rplot);
  pdf_lxg->plotOn(rplot);
  pdf_lxg->plotOn(rplot,Components(*rkpdf),LineColor(kRed));
  pdf_lxg->plotOn(rplot,Components(*lxg),LineColor(kGreen));
  pdf_lxg->paramOn(rplot,Format("NE", AutoPrecision(2)), Layout(0.65,0.89,0.89));
  rplot->getAttText()->SetTextSize(0.02);
  name = "fitSignal_"+label;
  rplot->SetTitle(name);
}

inline void SignalFit::NoiseFit(RooDataSet* data, RooRealVar* amp, TString label){

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  cout << endl;
  cout << "Fitting noise distribution on " << label <<  "." << endl;

  dataset = data;
  amp_rrvar = amp;

  if(set_range){
    name = "fitNoise_"+label;
    amp_rrvar->setRange(name, xl,xh);
  }

  name = "rkPDF_"+label;
  rkpdf = new RooKeysPdf(name,name,*amp_rrvar,*dataset,RooKeysPdf::NoMirror);

  name = "Nn_"+label;
  Nn = new RooRealVar(name, "N_{n}", 8000, "events");
  Nn->setConstant(kFALSE);

  if(set_range){
    name = "fitNoise_"+label;
    rkpdf->fitTo(*dataset, RooFit::Range(name));
  }
  else
    rkpdf->fitTo(*dataset);
}

inline void SignalFit::PlotNoise(Float_t xl, Float_t xh, Int_t nbins, TString label){

  rplot = amp_rrvar->frame(xl,xh,nbins);
  dataset->plotOn(rplot);
  rkpdf->plotOn(rplot);
  name = "fitNoise_"+label;
  rplot->SetTitle(name);
}

inline void SignalFit::PrintPlot(Float_t xl, Float_t xh, Int_t nbins, TString label,
				 std::function<void(Float_t,Float_t,Int_t,TString)> PlotFit){
  gROOT->SetBatch(kTRUE);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = (TCanvas*) new TCanvas("can","can",600,800);

  can->SetLogy();
  can->SetLeftMargin(0.15);

  PlotFit(xl,xh,nbins,label);

  rplot->Draw();
  CMSmark("");
  can->SaveAs(label+".pdf");        
  can->Close();

  delete can;
}

inline void SignalFit::PrintSignalPlot(Float_t xl, Float_t xh, Int_t nbins, TString label){
  PrintPlot(xl,xh,nbins,"signal_fit_"+label,
	    [this](Float_t xl, Float_t xh, Int_t nbins, TString label){PlotSignal(xl, xh, nbins, label);});
}

inline void SignalFit::PrintNoisePlot(Float_t xl, Float_t xh, Int_t nbins, TString label){
  PrintPlot(xl,xh,nbins,"noise_fit_"+label,
            [this](Float_t xl, Float_t xh, Int_t nbins, TString label){PlotNoise(xl, xh, nbins, label);});
}

inline void SignalFit::PrintParameters(TString label){

  cout << endl;
  cout << label << " fit parameters: " << endl;
  cout << "  MPV: " << MPV << " +/- " << MPV_err << endl;
  cout << "  Landau Width: " << sigmaL << endl;
  cout << "  Gauss Width: " << sigmaG << endl;
  cout << endl;
}

inline void SignalFit::CleanUp(){
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
}
