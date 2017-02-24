// a macro that simply extracts the number of jpsis from a list of files containing fit results

#include "TFile.h"
#include "TMath.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include <iostream>

using namespace std;
double RError(double,double,double,double);

void simplePrintResults() {
   vector<string> filenames;
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt4065_rap1519_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt65100_rap1519_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt100300_rap1519_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt5065_rap1015_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt65100_rap1015_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt100300_rap1015_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt65100_rap410_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt100300_rap410_cent0200.root");      
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt65100_rap-54_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt100300_rap-54_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt65100_rap-14-5_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt100300_rap-14-5_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt65100_rap-20-14_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt100300_rap-20-14_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt4065_rap-24-20_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt65100_rap-24-20_cent0200.root");
   filenames.push_back("/home/llr/cms/abdullah/pPbPsi2s/DimuonCADIs/HIN-16-004/Fitter/Output/pPbfitsPsi2sAllBins/ctauMass/DATA/result/FIT_CTAUMASS_DATA_Pbp_Bkg_Exponential_BkgNoPR_TripleDecay_CtauRes_TripleGaussianResolution_Psi2S_SingleCrystalBall_Psi2SNoPR_SingleSidedDecay_pt100300_rap-24-20_cent0200.root");

   const char* parname1 = "b_Psi2S_Pbp";
   const char* parname2 = "N_Psi2S_Pbp";
   double bPsi2s[20]= {0.0},  bPsi2sErr[20]= {0.0};
   double nPsi2s[20]= {0.0},  nPsi2sErr[20]= {0.0};
   
   double prPsi2s[20]= {0.0}, prPsi2sErr[20]= {0.0};
   double nprPsi2s[20]= {0.0}, nprPsi2sErr[20]= {0.0};

   int j = 0; 
   
   vector<string>::iterator it = filenames.begin();
   for (it; it<filenames.end(); it++) {
     TFile *f = new TFile(it->c_str());
     if (!f) {
       cout << "Error, " << *it << " not found" << endl;
       continue;
     }
     RooWorkspace *ws = (RooWorkspace*) f->Get("workspace");
     if (!ws) {
       cout << "Error, workspace not found in " << *it << endl;
       continue;
     }
     
     RooRealVar *var1 = ws->var(parname1);
     RooRealVar *var2 = ws->var(parname2);
     if (!ws) {
       cout << "Error, variable " << parname1<< " not found in " << *it << endl;
       cout << "Error, variable " << parname2<< " not found in " << *it << endl;
       continue;
     }

     bPsi2s[j] = var1->getVal(); bPsi2sErr[j] = var1->getError();
     nPsi2s[j] = var2->getVal(); nPsi2sErr[j] = var2->getError();
     nprPsi2s[j] = bPsi2s[j]*nPsi2s[j];
     prPsi2s[j] = (1-bPsi2s[j])*nPsi2s[j];

     nprPsi2sErr[j] = RError(bPsi2s[j],bPsi2sErr[j],nPsi2s[j],nPsi2sErr[j]);
     double errpr = TMath::Sqrt(pow(nprPsi2sErr[j], 2) + pow(nPsi2sErr[j],2));
     prPsi2sErr[j] = errpr;

     //cout << *it << " " << var1->getVal() << " +- " << var1->getError() << endl;
     //cout << *it << " " << var2->getVal() << " +- " << var2->getError() << endl;
     
     //cout << var1->getVal() << " +- " << var1->getError() << endl;
     //cout << var2->getVal() << " +- " << var2->getError() << endl;
     //cout << var1->getVal() << " +- " << var1->getError() << endl;
     
     j = j +1;
   }

   for(int i = 0; i<j; i++)
     {
       if(i==0) cout <<"bPsi2s { ";
       cout<<bPsi2s[i] <<" ,";
       if(i==j-1) cout <<" }"<<endl;
     }

   for(int i = 0; i<j; i++)
     {
       if(i==0) cout <<"bpsi2sErr { ";
       cout<<bPsi2sErr[i] <<" ,";
       if(i==j-1) cout <<" }"<<endl;
     }

   for(int i = 0; i<j; i++)
     {
       if(i==0) cout <<"nPsi2s { ";
       cout<<nPsi2s[i] <<" ,";
       if(i==j-1) cout <<" }"<<endl;
     }

   for(int i = 0; i<j; i++)
     {
       if(i==0) cout <<"nPsi2sErr { ";
       cout<<nPsi2sErr[i] <<" ,";
       if(i==j-1) cout <<" }"<<endl;
     }


   for(int i = 0; i<j; i++)
     {
       if(i==0) cout <<"prPsi2s { ";
       cout<<prPsi2s[i] <<" ,";
       if(i==j-1) cout <<" }"<<endl;
     }

   for(int i = 0; i<j; i++)
     {
       if(i==0) cout <<"prPsi2sErr { ";
       cout<<prPsi2sErr[i] <<" ,";
       if(i==j-1) cout <<" }"<<endl;
     }

   for(int i = 0; i<j; i++)
     {
       if(i==0) cout <<"nprPsi2s { ";
       cout<<nprPsi2s[i] <<" ,";
       if(i==j-1) cout <<" }"<<endl;
     }

   for(int i = 0; i<j; i++)
     {
       if(i==0) cout <<"nprPsi2sErr { ";
       cout<<nprPsi2sErr[i] <<" ,";
       if(i==j-1) cout <<" }"<<endl;
     }

}


double RError(double A, double eA, double B, double eB){
  double eR =0.0;
  double f=0;
  double fA=0;
  double fB=0;
  f=A/B;
  fA=eA/A;
  fB=eB/B;
  eR=  f*sqrt(fA*fA + fB*fB) ;
  return eR;
}

