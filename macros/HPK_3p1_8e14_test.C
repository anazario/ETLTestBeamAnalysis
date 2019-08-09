#include "../src/unbinned.cc"

void HPK_3p1_8e14_test(){
  Unbinned ub;

  ub.chainPath = "root://cmseos.fnal.gov//store/group/cmstestbeam/2019_04_April_CMSTiming/KeySightScope/RecoData/TimingDAQRECO/RecoWithTracks/v1/confInfo/";
  ub.run_start = new vector<int>{10473,11068};
  ub.run_end = new vector<int>{10602,11121};

  ub.tag = "HPK_3p1_8e14";

  float dx = 5.7;
  float dy = 0.34;
  float theta = TMath::ATan(dy/dx);

  ub.angle = theta;
  ub.satCut = 130.;
  ub.prec = 0.05;//window to cut noise (mm)

  ub.xcuts = {4.59,7.54,10.58,13.53,16.49};
  ub.ycuts = {31.2,32.2,33.2,34.2,35.2};
  ub.PrepareCuts();

  //ub.PrintGeomMap();
  //ub.SavePadNoise(0.,130.,100);
  //ub.Save2DPads(3.5,17.5,30.,36.);
  ub.SavePadSignalFit(0.,130.,100);
  
}

