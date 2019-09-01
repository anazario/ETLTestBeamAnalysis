#include "../src/TreeData.cc"

void Reprocess_HPK_3p1_8e14(){

  float dx = 5.7;
  float dy = 0.34;
  float theta = TMath::ATan(dy/dx);
  TreeData td;

  td.chainPath = "root://cmseos.fnal.gov//store/group/cmstestbeam/2019_04_April_CMSTiming/KeySightScope/RecoData/TimingDAQRECO/RecoWithTracks/v1/confInfo/";
  td.run_start = new vector<int>{10473,11068};
  td.run_end = new vector<int>{10602,11121};

  td.angle = theta;

  td.tag = "HPK_3p1_8e14";

  td.isReduced = true;
  td.FillTrees();

}
