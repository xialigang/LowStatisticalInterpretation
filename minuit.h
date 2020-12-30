#ifndef MINUIT_H
#define MINUIT_H

void fcn(Int_t &npar, Double_t *gin, Double_t &chi2, Double_t *par, Int_t iflag);

double minuit_fit(bool fix_mu = false, float mu0=1., int printlevel=-1, bool save_fit_results=0)
{
   // Get data points
   //get_input_data();
   // Prepare for fit
   const int npar=Nfit;
   TMinuit *gMinuit = new TMinuit(npar); 
   gMinuit->SetPrintLevel(printlevel);
   gMinuit->SetFCN(fcn);

   Double_t arglist[10];
   Int_t ierflg = 0;


   arglist[0] = 1;
   arglist[1] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

   // Set starting values and step sizes for parameters
   Double_t vstart[npar];// = {5, 0.0, 0.0, 0, 0, -2.0};
   //Double_t vstart[npar] = {10, 0, 0,  80, -20, -2};
   Double_t step[npar];// = {0.01, 0.01, 0.01, 0.1,0.1, 0.01};
   for(int i=0; i<npar; i++){
	vstart[i] = 1.0;
	step[i] = 0.001;
   }

   if(fix_mu){
       gMinuit->mnparm(0, "mu", mu0, 0, -10, 200,ierflg);
   }
   else{
       gMinuit->mnparm(0, "mu", 1, 0.1, -10, 200,ierflg);
   }
   // Now ready for minimization step
   arglist[0] = 1.;
   gMinuit->mnexcm("SIMPLEX", arglist ,0,ierflg);
   gMinuit->mnexcm("MIGRAD", arglist ,0,ierflg);
   gMinuit->mnexcm("HESSE", arglist ,0,ierflg);
   //gMinuit->mnexcm("MINOS", arglist ,0,ierflg);

   // Print results
   Double_t fmin,fedm,errdef;
   Int_t nvpar,nparx,icstat;
   //Double_t covmat[npar][npar];
   //Double_t rho[npar][npar];
   gMinuit->mnstat(fmin,fedm,errdef,nvpar,nparx,icstat);
   if(save_fit_results){
       for(int i=0; i<Nfit; i ++){
           gMinuit->GetParameter(i, xfit[i], xerr[i]);
       }
       for(int i=0; i<npar; i++){
           gMinuit->mnerrs(i,xup[i],xlow[i],xerr[i],globcc[i]);
       }
   }
   //gMinuit->mnemat(&covmat[0][0],npar); 
   //rho[0][1]=covmat[0][1]/(mass_err*eff_err);
   //rho[0][2]=covmat[0][2]/(mass_err*xbg_err);
   //rho[1][2]=covmat[1][2]/(eff_err*xbg_err);

   //cout<<"fmin = "<<fmin<<endl;
   return fmin;
}

void fcn(Int_t &npar, Double_t *gin, Double_t &chi2, Double_t *par, Int_t iflag)
{
    // Calculate log-likelihood


   chi2 = 0;

   double mu = par[0];

   double logL = 0.;
   double bi, si, ni;
   for(int i=0; i<nbins; i++){
       bi = h_bkg -> GetBinContent(i+1);
       si = h_sig -> GetBinContent(i+1);
       ni = h_data_tmp -> GetBinContent(i+1);

       logL += ni*log(bi+mu*si) - (bi+mu*si);
   }
   chi2 = -2*logL;

   //cout<<"lnL = "<<log_l<<endl;
   //cout<<"chi2/ndof = "<<(double)chi2/(Nscan-3)<<endl; 
}


#endif
