#include "../src/Unbinned.cc"

void test_ub(){
  
  Unbinned ub;
  ub.tag = "HPK_3p1_8e14";
  ub.chainPath = Form("../ROOT/%s_trees/",ub.tag.Data());
  //ub.isDebug = true;
  ub.satCut = 130.;

  ub.frame_margins = {3.5,17.5,30.,36.};
  ub.xcuts = {4.59,7.54,10.58,13.53,16.49};
  ub.ycuts = {31.2,32.2,33.2,34.2,35.2};

  //ub.SavePadNoise(1);
  //ub.SavePadNoise(2);
  //ub.SavePadInfo();
  ub.SavePadSignal();
}
