#ifndef FUNCTION_H
#define FUNCTION_H
#include "TMath.h"
#include "minuit.h"
#include <iomanip>

void load_hists(){


    TFile* input = new TFile("input.root", "read");
    h_sig = (TH1F*) input -> Get("h_sig");
    h_bkg = (TH1F*) input -> Get("h_bkg");
    h_data = (TH1F*) input -> Get("h_data");

    cout<<"Reading histograms ..."<<endl;
    nbins = h_sig -> GetNbinsX();
    cout<<"nbins = "<<nbins<<endl;
    for(int i=0; i<nbins; i++){
        float s = h_sig -> GetBinContent(i+1);
        float b = h_bkg -> GetBinContent(i+1);
        float n = h_data -> GetBinContent(i+1);
        cout<<"bin "<<i+1<<"    sig: "<<s<<"    bkg: "<<b<<"    data: "<<n<<endl;
    }
}
double get_qobs(double mu_min =0, double mu_max = 10){

    double dmu = 0.01;
    const int Nmu = (mu_max - mu_min)/dmu;
    double qobs = 9999;
    for(int i=0; i<Nmu; i++){
        double mu = mu_min + i*dmu;
        double qtmp = minuit_fit(1,mu);
        if(qtmp < qobs){
            qobs = qtmp;
        }
    }
    return qobs;
}

void first_fit(){
    h_data_tmp = (TH1F*) h_data -> Clone("h_data_tmp");
    q0 = minuit_fit(1,0,-1);
    cout<<"q0 = "<<q0<<endl;
    qobs = minuit_fit(0,0,-1,1);
    cout<<"qobs = "<<qobs<<endl;
    for(int i=0; i<Nfit; i ++){
        cout<<"par "<<i<<"   "<<xfit[i]<<"  +/- "<<xerr[i]<<endl;
    }
    if(xfit[0]<0 || qobs>q0){
        qobs = q0;
        xfit[0] = 0.;
        cout<<"force mu=0 and qobs=q0!"<<endl;
    }
    tobs0 = q0 - qobs;
    cout<<"tobs0 = "<<tobs0<<endl;

}
double get_tobsmu(double mu=1.){

    h_data_tmp = (TH1F*) h_data -> Clone("h_data_tmp");
    //for(int i=0; i<nbins; i++) cout<<"bin "<<i<<" ndata = "<<h_data_tmp -> GetBinContent(i+1)<<endl;
    //cout<< "prob(0|8.877) = "<<TMath::Poisson(0,8.877)<<endl;
    double _qobs = minuit_fit(0,0,-1,1);
    cout<<"lazy qobs = "<<_qobs<<endl;
    double qobs = get_qobs();
    cout<<"qobs = "<<qobs<<endl;
    double _q0 = minuit_fit(1,0,-1);
    cout<<"q0 = "<<_q0<<endl;
    for(int i=0; i<Nfit; i ++){
        cout<<"par "<<i<<"   "<<xfit[i]<<"  +/- "<<xerr[i]<<endl;
    }
    /*
    if(xfit[0]<0 || _qobs>_q0){
        _qobs = _q0;
        xfit[0] = 0.;
        cout<<"force mu=0 and qobs=q0!"<<endl;
    }
    */
    double _qmu = minuit_fit(1,mu,-1);
    cout<<"qmu = "<<_qmu<<endl;
   

    double _tobs0 = _q0 - _qobs;
    double _tobsmu = _qmu - _qobs;

    cout<<"tobs0 = "<<_tobs0<<endl;
    cout<<"tobsmu = "<<_tobsmu<<" for mu = "<<mu<<endl;
    return _tobsmu;

}

void get_pdf_tmu(double &alpha, double &beta, TH1F* &h_muhat, TH1F* &h_tmu, double mu0=0, double mu=0, TString nametag="0", int ntoys = 100, int nmax=100){
//first mu is the "mu" in test statistics t_mu
//second mu is the "mu" in hypothesis data
    tobsmu = get_tobsmu(mu0);
    h_data_tmp = (TH1F*) h_data -> Clone("h_data_tmp");
    TH1F* h_data_mu = (TH1F*) h_bkg -> Clone("h_data_mu");
    double s = h_sig -> Integral();
    double b = h_bkg -> Integral();
    double n_mu = b + mu*s;
    //cout<<"n_mu = "<<n_mu<<endl;
    if(mu>0){
        double bi = 0;
        double si = 0;
        for(int i=0; i<nbins; i++){
            bi = h_bkg -> GetBinContent(i+1);
            si = h_sig -> GetBinContent(i+1);
            h_data_mu -> SetBinContent(i+1, bi+mu*si);
            //cout<<"bin i :"<<bi<<", "<<si<<", "<<h_data_mu -> GetBinContent(i+1)<<endl;
        }
    }
    //TH1F* h_tmu = new TH1F("h_tmu_"+nametag, "", 100, 0, 50);
    h_tmu = new TH1F("h_tmu_"+nametag, "", 50, 0, 50);
    double muhatmax = 100;
    if(mu0>=50) muhatmax = 200;
    h_muhat = new TH1F("h_muhat_"+nametag, "", int(muhatmax), 0, muhatmax);
    double muhat = 0;
    double tmu = -1;
    
    alpha = 0;
    beta = 0;
    double tot_prob = 0;
    //special case n = 0;
    muhat = 0.;
    tmu = 2*mu0*s;
    double prob0 = TMath::Poisson(0,n_mu);
    h_tmu -> Fill(tmu, ntoys*prob0);
    h_muhat -> Fill(muhat, ntoys*prob0);
    if(tmu <= tobsmu) alpha += ntoys*prob0;
    if(tmu >= tobsmu) beta += ntoys*prob0;
    tot_prob = TMath::Poisson(0,n_mu);
    //rest cases n > 0
    for(int in = 1; in<=nmax; in++){
        double prob_in = TMath::Poisson(in,n_mu);
        if(prob_in < 1e-7) {
            cout<<"prob( "<<in<<" | "<<n_mu<<" ) = "<<prob_in<<" < 1e-7"<<endl;
            break;
        }
        tot_prob += prob_in;
        for(int itoy=0; itoy<ntoys; itoy++){
            //initialize toy data 
            for(int i=0; i<nbins; i++){
                h_data_tmp -> SetBinContent(i+1,0);
            }
            //generate toy data
            for(int i=0; i<in; i++){
                //cout<<"before"<<endl;
                double x = h_data_mu -> GetRandom();
                //cout<<"after"<<endl;
                h_data_tmp -> Fill(x);
            }
            //perform fits
            //fit with fixing mu=0
            q0 = minuit_fit(1,0);
            qobs = minuit_fit(0, 0, -1, 1);
            muhat = xfit[0];
            if(xfit[0]<0 || qobs > q0){
                muhat = 0.;
                qobs = q0;
            }
            qmu = minuit_fit(1,mu0);
            tmu = qmu - qobs;
            h_muhat -> Fill(muhat, prob_in);
            h_tmu -> Fill(tmu, prob_in);
            if(tmu <= tobsmu) alpha += prob_in;
            if(tmu >= tobsmu) beta += prob_in;
        }
    }
    cout<<"n_mu = "<<n_mu<<endl;
    cout<<"tobsmu = "<<tobsmu<<endl;
    alpha = alpha/ntoys/tot_prob;
    beta = beta/ntoys/tot_prob;
    cout<<"Prob(t"<<mu<<"<=tobsmu) = "<<alpha<<endl;
    cout<<"Prob(t"<<mu<<">=tobsmu) = "<<beta<<endl;
    cout<<"h_tmu integral = "<<h_tmu -> Integral()<<endl;
    h_tmu -> Scale(1./h_tmu->Integral());
    h_muhat -> Scale(1./h_muhat->Integral());
    //return h_tmu;
}

void cal_p0(int ntoys=100, bool logy=0){
   
    TH1F* h_muhat;
    TH1F* h_t0;
    double alpha;
    double beta;
    //TH1F* h_t0 = get_pdf_tmu(0., 0., "0", ntoys);
    get_pdf_tmu(alpha, beta, h_muhat, h_t0, 0., 0., "0", ntoys);
    cout<<"alpha, beta = "<<alpha<<" , "<<beta<<endl;
    ofstream outfile("results_p0.txt");
    outfile<<"alpha(t0) = "<<alpha<<endl;
    outfile<<"beta(t0) = "<<beta<<endl;
    outfile.close();
    TCanvas* Cs_p0 = new TCanvas("Cs_p0", "", 10, 10, 1000, 500);
    Cs_p0 -> Divide(2,1);
    Cs_p0 -> cd(1);
    if(logy) gPad -> SetLogy();
    h_muhat -> Draw("PE");
    h_muhat -> GetXaxis() -> SetTitle("#hat{#mu} (#mu_{H}=0)");
    h_muhat -> GetYaxis() -> SetTitle("PDF");
    if(logy)
        h_muhat -> GetYaxis() -> SetRangeUser(1e-8, 9.);
    else
        h_muhat -> GetYaxis() -> SetRangeUser(1e-8, 1.);
    Cs_p0 -> cd(2);
    if(logy) gPad -> SetLogy();
    h_t0 -> Draw("PE");
    h_t0 -> GetXaxis() -> SetTitle("t_{0} (#mu_{H}=0)");
    h_t0 -> GetYaxis() -> SetTitle("PDF");
    if(logy)
        h_t0 -> GetYaxis() -> SetRangeUser(1e-8, 9.);
    else
        h_t0 -> GetYaxis() -> SetRangeUser(1e-8, 1.);
    gPad -> Update();
    TArrow* a_tobs0 = new TArrow(tobs0, 1.0, tobs0, 1e-8, 0.017, "|>");
    a_tobs0 -> Draw();
    a_tobs0 -> SetLineWidth(2);
    a_tobs0 -> SetLineColor(kBlue);
    a_tobs0 -> SetFillColor(kBlue);
    Cs_p0 -> SaveAs("Cs_p0.png");
    Cs_p0 -> SaveAs("Cs_p0.pdf");

}
void cal_upperlimit(double mu=1., int ntoys=100, bool logy=0){
   
    TH1F* h_muhat_H0;
    TH1F* h_tmu_H0;
    double alpha_H0;
    double beta_H0;

    TH1F* h_muhat_Hmu;
    TH1F* h_tmu_Hmu;
    double alpha_Hmu;
    double beta_Hmu;
    
    get_pdf_tmu(alpha_H0, beta_H0, h_muhat_H0, h_tmu_H0, mu, 0, "muH0", ntoys);
    get_pdf_tmu(alpha_Hmu, beta_Hmu, h_muhat_Hmu, h_tmu_Hmu, mu, mu, "muHmu", ntoys);
    tobsmu = get_tobsmu(mu);
    
    double CLs = beta_Hmu/beta_H0;
    cout<<"*****************************CLs = "<<CLs<<endl;
    ofstream outfile(Form("results_upperlimit_mu%.0f.txt",mu));
    outfile<<"mu = "<<mu<<endl;
    outfile<<"alpha_H0 = "<<alpha_H0<<endl;
    outfile<<"beta_H0 = "<<beta_H0<<endl;
    outfile<<"alpha_Hmu = "<<alpha_Hmu<<endl;
    outfile<<"beta_Hmu = "<<beta_Hmu<<endl;
    outfile<<"tobsmu = "<<tobsmu<<endl;
    outfile<<"CLs = "<<CLs<<endl;
    outfile.close();


    TCanvas* Cs_upperlimit = new TCanvas("Cs_upperlimit", "", 10, 10, 1000, 500);
    Cs_upperlimit -> Divide(2,1);
    Cs_upperlimit -> cd(1);
    if(logy) gPad -> SetLogy();
    h_muhat_H0 -> Draw("PE");
    h_muhat_Hmu -> Draw("PE,same");
    h_muhat_Hmu -> SetMarkerStyle(24);
    h_muhat_Hmu -> SetMarkerColor(kRed);
    h_muhat_Hmu -> SetLineColor(kRed);
    h_muhat_H0 -> GetXaxis() -> SetTitle(Form("#hat{#mu} (#mu=%.0f)",mu));
    h_muhat_H0 -> GetYaxis() -> SetTitle("PDF");
    if(logy)
        h_muhat_H0 -> GetYaxis() -> SetRangeUser(1e-8, 9.);
    else
        h_muhat_H0 -> GetYaxis() -> SetRangeUser(1e-8, 1.);
    TLegend* leg_muhat = new TLegend(0.5, 0.7, 0.95, 0.90);
    leg_muhat -> AddEntry(h_muhat_H0, "#mu_{H}=0", "PE");
    leg_muhat -> AddEntry(h_muhat_Hmu, "#mu_{H}=#mu", "PE");
    leg_muhat -> SetBorderSize(0);
    leg_muhat -> SetFillStyle(0);
    leg_muhat -> SetFillColor(0);
    leg_muhat -> Draw();
    Cs_upperlimit -> cd(2);
    if(logy) gPad -> SetLogy();
    h_tmu_H0 -> Draw("PE");
    h_tmu_Hmu -> Draw("PE,same");
    h_tmu_Hmu -> SetMarkerStyle(24);
    h_tmu_Hmu -> SetMarkerColor(kRed);
    h_tmu_Hmu -> SetLineColor(kRed);
    h_tmu_H0 -> GetXaxis() -> SetTitle(Form("t_{#mu} (#mu=%.0f)",mu));
    h_tmu_H0 -> GetYaxis() -> SetTitle("PDF");
    if(logy)
        h_tmu_H0 -> GetYaxis() -> SetRangeUser(1e-8, 9.);
    else
        h_tmu_H0 -> GetYaxis() -> SetRangeUser(1e-8, 1.);
    TArrow* a_tobsmu = new TArrow(tobsmu, 1., tobsmu, 1e-8, 0.017, "|>");
    a_tobsmu -> Draw();
    a_tobsmu -> SetLineWidth(2);
    a_tobsmu -> SetLineColor(kBlue);
    a_tobsmu -> SetFillColor(kBlue);
    TLegend* leg_tmu = new TLegend(0.5, 0.6, 0.95, 0.90);
    leg_tmu -> AddEntry(h_tmu_H0, "#mu_{H}=0", "PE");
    leg_tmu -> AddEntry(h_tmu_Hmu, "#mu_{H}=#mu", "PE");
    leg_tmu -> AddEntry(a_tobsmu, "t_{#mu}(obs.)", "L");
    leg_tmu -> SetBorderSize(0);
    leg_tmu -> SetFillStyle(0);
    leg_tmu -> SetFillColor(0);
    leg_tmu -> Draw();
    Cs_upperlimit -> SaveAs(Form("Cs_upperlimit_mu%.0f.png",mu));
    Cs_upperlimit -> SaveAs(Form("Cs_upperlimit_mu%.0f.pdf",mu));
}




#endif
