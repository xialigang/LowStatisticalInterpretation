#include "quantity.h"
#include "function.h"
#include "minuit.h"
#include "AtlasStyle.C"

int main(){
/*
   gROOT -> SetStyle("Plain");
   gROOT -> ForceStyle();
   gStyle -> SetOptStat(0);
   gStyle -> SetFillColor(10);
   gStyle -> SetFrameLineWidth(2);
*/
//gROOT -> LoadMacro("AtlasStyle.C")
   SetAtlasStyle();

   load_hists();
   first_fit();
   cal_p0(2000,1);
   //cal_upperlimit(1,1000,1);
   //cal_upperlimit(10,1000,1);
   //cal_upperlimit(30,1000,1);
   //cal_upperlimit(50,1000,1);
   cal_upperlimit(75,1000,1);
   //cal_upperlimit(100,1000,1);
   
   
   return 0;
}
