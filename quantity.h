#ifndef QUANTITY_H
#define QUANTITY_H
#include<iomanip>

const int Nfit = 1;
double xfit[Nfit];
double xerr[Nfit];
double xlow[Nfit];
double xup[Nfit];
double globcc[Nfit];

//ifstream infile_sb("sig_bkg.txt");
//ifstream infile_unc("syst_unc.txt");
//ofstream outfile_results;


TH1F* h_sig;
TH1F* h_bkg;
TH1F* h_data;
TH1F* h_data_tmp;
TH1F* h_syst_1_sig;
TH1F* h_syst_1_bkg;

int nbins=-1;

double qobs = 0;
double q0 = 0;
double qmu = 0;
double t0 = 0;
double tmu = 0;
double tobs0 = 0;
double tobsmu = 0;

#endif
