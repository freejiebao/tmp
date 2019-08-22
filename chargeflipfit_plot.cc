#include <iostream>
#include "TFile.h"
#include "TH2D.h"
using namespace std;
// negative log-likelihood function
void negative_log_likelihood(Int_t &npar, Double_t *gin, Double_t &nll, Double_t *par,Int_t iflag) {
 // the count_chargeflip__DY_2018.root is 2D histogram related to of two leptons
 TFile* fin=TFile::Open("count_chargeflip__DY_2018.root");
 TH2D *h_ss=(TH2D*)fin->Get("h_ss");
 TH2D *h_os=(TH2D*)fin->Get("h_os");
 Double_t fun;
 Double_t sum=0.;
 Int_t nss,nos;
 for (Int_t i = 0; i < 5; i++) {
   for (Int_t j = 0; j < 5; j++) {
    nss=Long_t(h_ss->GetBinContent(i+1,j+1)+0.5);
    nos=Long_t(h_os->GetBinContent(i+1,j+1)+0.5);
    Double_t Log_Factorial=0;
    // Factorial of same sign events is quite large, so cannot calculate factorial directly
    for(Int_t k=0; k<nss;k++){
    Log_Factorial+=TMath::Log(Double_t(k+1));
    }
    fun=nss*TMath::Log(nos*(par[i]+par[j]))-nos*(par[i]+par[j])-Log_Factorial;
    sum += fun;
   }
 }
 //Double_t fi = f(&data[i], par);
 //sum += TMath::Log(fi);

 nll = -sum;
 //cout<<"nll:"<<nll<<endl;
}
void chargeflipfit_plot(){
// prepare minuit
 Int_t nPar = 5; // number of fit parameters
 TMinuit m(nPar);
 m.SetFCN(negative_log_likelihood);
 m.SetPrintLevel(0); // -1 quiet, 0 normal, 1 verbose
 // 1 for chi2 fit, 0.5 for negative log-likelihood fir
 // see section 1.4.1 in MINUIT manual, e.g., http://hep.fi.infn.it/minuit.pdf
 m.SetErrorDef(1);
 // parameters:
 // parameter no., name, start value, step size, range min., range max.
 // range min = range max = 0 -> no limits
 m.DefineParameter(0, "p_{#eta0}", 0.000023, 0.000000001, 0, 0.005);
 m.DefineParameter(1, "p_{#eta1}", 0.000077, 0.000000001, 0, 0.005);
 m.DefineParameter(2, "p_{#eta2}", 0.000345, 0.000000001, 0, 0.5);
 m.DefineParameter(3, "p_{#eta3}", 0.002251, 0.000000001, 0, 0.5);
 m.DefineParameter(4, "p_{#eta4}", 0.002174, 0.000000001, 0, 0.5);
 // now ready for minimization step
 m.Migrad();
 m.Command("SHOW COV"); // show covariance matrix
 // draw fit
 Double_t p0, p0_err,p1, p1_err,p2, p2_err,p3, p3_err,p4, p4_err;
 m.GetParameter(0, p0, p0_err);
 m.GetParameter(1, p1, p1_err);
 m.GetParameter(2, p2, p2_err);
 m.GetParameter(3, p3, p3_err);
 m.GetParameter(4, p4, p4_err);
 cout<<p0 <<"\t"<<p1 <<"\t"<<p2 <<"\t"<<p3 <<"\t"<<p4 <<"\n";
 cout<<p0_err <<"\t"<<p1_err <<"\t"<<p2_err <<"\t"<<p3_err <<"\t"<<p4_err <<"\n";
}
