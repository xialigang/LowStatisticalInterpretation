/////////////////////// trash codes? /////////////////////////////
#if 0
    for(int itoy=0; itoy<ntoys; itoy++){
        int n = gRandom -> Poisson(n_mu);
        if(n==0) {
            muhat = 0;
            tmu = 2*mu*s;
        }
        else{
            //initialize toy data 
            for(int i=0; i<nbins; i++){
                h_data -> SetBinContent(i+1,0);
            }
            //generate toy data
            for(int i=0; i<n; i++){
                double x = h_data_mu -> GetRandom();
                h_data -> Fill(x);
            }
            //perform fits
            //fit with fixing mu=0
            q0 = minuit_fit(1,0);
            qobs = minuit_fit(0);
            if(xfit[0]<0 || qobs > q0){
                muhat = 0.;
                qobs = q0;
            }
            qmu = minuit_fit(1,mu);
            tmu = qmu - qobs;
        }
        h_tmu -> Fill(tmu);
    }
#endif

