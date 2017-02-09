#ifndef results2tree_C
#define results2tree_C

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooHist.h"

#include <vector>
#include <cstring>

#include "Macros/Utilities/resultUtils.h"
#include "Macros/Utilities/EVENTUTILS.h"

using namespace std;
using namespace RooFit;

struct poi {
   Char_t name[64];
   float val;
   float err;
};

const int nBins = 46;

void results2tree(
      const char* workDirName, 
      const char* DSTag, //="DATA", // Data Set tag can be: "DATA","MCPSI2SP", "MCJPSIP" ...
      const char* prependPath, //="",
      const char* fitType, // "mass", "ctau"...
      const char* thePoiNames, //="N_Jpsi,b_Jpsi,f_Jpsi,m_Jpsi,sigma1_Jpsi,alpha_Jpsi,n_Jpsi,sigma2_Jpsi,MassRatio,rSigma21_Jpsi,
      // lambda1_Bkg,lambda2_Bkg,lambda3_Bkg,lambda4_Bkg,lambda5_Bkg,N_Bkg,b_Bkg,
      // ctau1_CtauRes,ctau2_CtauRes,f_CtauRes,rSigma21_CtauRes,sigma1_CtauRes,
      // fDFSS_BkgNoPR,fDLIV_BkgNoPR,lambdaDDS_BkgNoPR,lambdaDF_BkgNoPR,lambdaDSS_BkgNoPR,lambdaDSS_JpsiNoPR,
      // eff,effnp,lumi,taa,ncoll,npart,correl_N_Jpsi_vs_b_Jpsi",
      bool wantPureSMC //=false
      ) {
   // workDirName: usual tag where to look for files in Output
   // thePoiNames: comma-separated list of parameters to store ("par1,par2,par3"). Default: all

   TFile *f = new TFile(treeFileName(workDirName,DSTag,prependPath,fitType),"RECREATE");
   TTree *tr = new TTree("fitresults","fit results");

   // bin edges
   float ptmin, ptmax, ymin, ymax, centmin, centmax;
   // model names
   Char_t jpsiName[128], psi2sName[128], bkgName[128];
   // collision system
   Char_t collSystem[8];
   // goodness of fit
   float nll, chi2, chi2_invMass, chi2_ctau, chi2prob, normchi2; int npar, nparbkg, ndof , ndof_invMass, ndof_ctau;
   // parameters to store: make it a vector
   vector<poi> thePois;
   TString thePoiNamesStr(thePoiNames);
   TString t; Int_t from = 0;
   while (thePoiNamesStr.Tokenize(t, from , ",")) {
      poi p; strcpy(p.name, t.Data());
      cout << p.name << endl;
      thePois.push_back(p);
   }

   // create tree branches
   tr->Branch("ptmin",&ptmin,"ptmin/F");
   tr->Branch("ptmax",&ptmax,"ptmax/F");
   tr->Branch("ymin",&ymin,"ymin/F");
   tr->Branch("ymax",&ymax,"ymax/F");
   tr->Branch("centmin",&centmin,"centmin/F");
   tr->Branch("centmax",&centmax,"centmax/F");
   tr->Branch("jpsiName",jpsiName,"jpsiName/C");
   tr->Branch("psi2sName",psi2sName,"psi2sName/C");
   tr->Branch("bkgName",bkgName,"bkgName/C");
   tr->Branch("collSystem",collSystem,"collSystem/C");
   tr->Branch("nll",&nll,"nll/F");
   tr->Branch("chi2",&chi2,"chi2/F");
   tr->Branch("chi2_invMass",&chi2_invMass,"chi2_invMass/F");
   tr->Branch("chi2_ctau",&chi2_ctau,"chi2_ctau/F");
   tr->Branch("chi2prob",&chi2prob,"chi2prob/F");
   tr->Branch("normchi2",&normchi2,"normchi2/F");
   tr->Branch("npar",&npar,"npar/I");
   tr->Branch("nparbkg",&nparbkg,"nparbkg/I");
   tr->Branch("ndof",&ndof,"ndof/I");
   tr->Branch("ndof_invMass",&ndof_invMass,"ndof_invMass/I");
   tr->Branch("ndof_ctau",&ndof_ctau,"ndof_ctau/I");

   for (vector<poi>::iterator it=thePois.begin(); it!=thePois.end(); it++) {
      tr->Branch(Form("%s_val",it->name),&(it->val),Form("%s_val/F",it->name));
      tr->Branch(Form("%s_err",it->name),&(it->err),Form("%s_err/F",it->name));
   }

   // list of files
   vector<TString> theFiles = fileList(workDirName,"",DSTag,"",fitType);
        
        cout << theFiles.size() << endl;

   int cnt=0;
   for (vector<TString>::const_iterator it=theFiles.begin(); it!=theFiles.end(); it++) {
      cout << "Parsing file " << cnt << " / " << theFiles.size() << ": " << *it << endl;

      // parse the file name to get info
      anabin thebin = binFromFile(*it);
      ptmin = thebin.ptbin().low();
      ptmax = thebin.ptbin().high();
      ymin = thebin.rapbin().low();
      ymax = thebin.rapbin().high();
      centmin = thebin.centbin().low();
      centmax = thebin.centbin().high();
      strcpy(collSystem, (it->Contains("Pbp")) ? "Pbp" : "PP");
      bool isPP = !(it->Contains("Pbp"));

      //cout<<"==================ok1========================"<<endl;
      
      // get the model names
      from = 0;
      bool catchjpsi=false, catchbkg=false, catchtype=false, catchpsi2s=false;
      TString modelType;
      while (it->Tokenize(t, from, "_")) {
	//cout <<"t.Data()  "<<t.Data()  << endl;
	if (catchjpsi) {strcpy(jpsiName, t.Data()); catchjpsi=false;}
	if (catchpsi2s) {strcpy(psi2sName, t.Data()); catchpsi2s=false;}
	if (catchbkg) {strcpy(bkgName, t.Data()); catchbkg=false;}
	if (catchtype) {modelType = t; catchtype=false;}
	if (t=="Jpsi") catchjpsi=true;
	if (t=="Psi2S") catchpsi2s=true;
	if (t=="Bkg") catchbkg=true;
	if (t.EndsWith("FIT")) catchtype=true;
      }

      TFile *f = TFile::Open(*it); RooWorkspace *ws = NULL;
      if (!f) {
         cout << "Error, file " << *it << " does not exist." << endl;
      } else {
         ws = (RooWorkspace*) f->Get("workspace");
         if (!ws) {
            cout << "Error, workspace not found in " << *it << "." << endl;
         }
      }

      nll=0; chi2=0; npar=0; ndof=0;
      chi2_invMass=0; ndof_invMass=0;
      chi2_ctau=0; ndof_ctau=0;
      if (f && ws) {
         // get the model for nll and npar
         RooAbsPdf *model = pdfFromWS(ws, Form("_%s",collSystem), Form("pdf%s_Tot",modelType.Data()));

         RooAbsPdf *model_bkg = pdfFromWS(ws, Form("_%s",collSystem), Form("pdf%s_Bkg",modelType.Data()));
         const char* token = (strcmp(DSTag,"DATA") && wantPureSMC) ? Form("_%s_NoBkg",collSystem) : Form("_%s",collSystem);
         RooAbsData *dat = dataFromWS(ws, token, Form("dOS_%s", DSTag));
         if (dat) {
            if (model) {
               RooAbsReal *NLL = model->createNLL(*dat);
               if (NLL) nll = NLL->getVal();
               npar = model->getParameters(dat)->selectByAttrib("Constant",kFALSE)->getSize();

               // compute the chi2 and the ndof
               RooRealVar *chi2var = ws->var("chi2");
               RooRealVar *ndofvar = ws->var("ndof");
               RooRealVar *chi2var_invMass = ws->var("chi2_invMass");
               RooRealVar *ndofvar_invMass = ws->var("ndof_invMass");
               RooRealVar *chi2var_ctau = ws->var("chi2_ctau");
               RooRealVar *ndofvar_ctau = ws->var("ndof_ctau");
               if (chi2var && ndofvar) {
                  chi2 = chi2var->getVal();
                  ndof = ndofvar->getVal();
               } else {
                  RooPlot* frame = ws->var("invMass")->frame(Bins(nBins));
                  dat->plotOn(frame, DataError(RooAbsData::SumW2), XErrorSize(0));
                  model->plotOn(frame, Precision(1e-4), Range("invMass"));
                  TH1 *hdatact = dat->createHistogram("hdatact", *(ws->var("invMass")), Binning(nBins));
                  RooHist *hpull = frame->pullHist(0,0, true);
                  double* ypulls = hpull->GetY();
                  unsigned int nFullBins = 0;
                  for (int i = 0; i < nBins; i++) {
                     if (hdatact->GetBinContent(i+1) > 0.0) {
                        chi2 += ypulls[i]*ypulls[i];
                        nFullBins++;
                     }
                  }
                  ndof = nFullBins - npar;
               }

               normchi2 = chi2/ndof;
               chi2prob = TMath::Prob(chi2,ndof);
            }
            if (model_bkg) {
               nparbkg = model_bkg->getParameters(dat)->selectByAttrib("Constant",kFALSE)->getSize();
            }
         }

         // get the POIs
         for (vector<poi>::iterator itpoi=thePois.begin(); itpoi!=thePois.end(); itpoi++) {
            RooRealVar *thevar = poiFromWS(ws, Form("_%s",collSystem), itpoi->name);
            itpoi->val = thevar ? thevar->getVal() : 0;
            itpoi->err = thevar ? thevar->getError() : 0;
	    
           
            if (TString(itpoi->name).Contains("correl_")) {
               // correlation between two variables
               TString toparse = TString(itpoi->name).ReplaceAll("correl_","");
               TString var1name, var2name;
               from = 0; int i=0;
               while (toparse.Tokenize(t, from, "_vs_")) {
                  if (i==0) var1name = t;
                  if (i==1) var2name = t;
                  i++;
               }
              
               RooFitResult *fr = (RooFitResult*) ws->obj(Form("fitResult_pdf%s_Tot_%s",modelType.Data(),collSystem));
               RooRealVar *thevar1 = poiFromWS(ws, Form("_%s",collSystem), var1name);
               RooRealVar *thevar2 = poiFromWS(ws, Form("_%s",collSystem), var2name);
               if (fr) {
                  itpoi->val = fr->correlation(Form("%s_%s",var1name.Data(),collSystem),Form("%s_%s",var2name.Data(),collSystem));
                  itpoi->err = thevar1 && thevar2 ? itpoi->val * thevar1->getError() * thevar2->getError() : 0;
               }
            } else if (TString(itpoi->name).Contains("eff")) {
               // efficiency
               TFile *feff = TFile::Open(Form("../Efficiency/files/histos_%s_%s.root", 
                        TString(itpoi->name)=="effnp" ? "npjpsi" : "jpsi", 
                        isPP ? "pp" : "pbpb"));
               bool isallcent = (thebin.centbin() == binI(0,200));
               int catmin, catmax;
               bool isallrap = (thebin.rapbin() == binF(0,2.4));
               bool islowpt = (thebin.ptbin() == binF(3,6.5));
               bool ishighpt = ((thebin.ptbin() == binF(6.5,50)) || (thebin.ptbin() == binF(6.5,30)) );
               if (isallcent || (!isallcent && ishighpt)) {
                  catmin = thebin.rapbin().low()*10;
                  catmax = thebin.rapbin().high()*10;
               } else {
                  catmin = thebin.centbin().low()/2;
                  catmax = thebin.centbin().high()/2;
               }

               TString hnumname, hdenname;
               TString tag1, tag2;
               if (isallcent && ishighpt) {
                  tag1 = "rap";
                  hnumname = "hnum_rap";
                  hdenname = "hden_rap";
               } else if (!islowpt) {
                  tag2 = (isallcent || (!isallcent && ishighpt)) ? "rap" : "cent";
                  tag1 = ishighpt ? "cent" : "pt";
                  hnumname = Form("hnum_%s_%s%02i%02i", tag1.Data(), tag2.Data(), catmin, catmax);
                  hdenname = Form("hden_%s_%s%02i%02i", tag1.Data(), tag2.Data(), catmin, catmax);
               } else {
                  tag1 = "cent";
                  hnumname = "hnum_cent_rap1824_pt3065";
                  hdenname = "hden_cent_rap1824_pt3065";
               }
               TH1F *hnum = (TH1F*) feff->Get(hnumname);
               TH1F *hden = (TH1F*) feff->Get(hdenname);
               if (!hnum || !hden) {
                  thebin.print();
                  cout << hnumname << " not found!" << endl;
                  itpoi->val = 0;
                  itpoi->err = 0;
                  continue;
               }

               double numval, numerr, denval, denerr;
               int ibin = hnum->FindBin((thebin.centbin().low()+thebin.centbin().high())/4.);
               if (tag1 == "pt") ibin = hnum->FindBin((thebin.ptbin().low()+thebin.ptbin().high())/2.);
               if (tag1 == "rap") ibin = hnum->FindBin((thebin.rapbin().low()+thebin.rapbin().high())/2.);
               numval = hnum->GetBinContent(ibin);
               numerr = hnum->GetBinError(ibin);
               denval = hden->GetBinContent(ibin);
               denerr = hden->GetBinError(ibin);
               // special case of the all integrated efficiency
               if ((isallcent && isallrap && ishighpt) || (isallcent && islowpt)) {
                  numval = hnum->IntegralAndError(1,hnum->GetNbinsX(),numerr);
                  denval = hden->IntegralAndError(1,hden->GetNbinsX(),denerr);
               }
               double efficiency = (denval>0) ? numval / denval : 0;
               itpoi->val = efficiency;
               itpoi->err = (numval>0 && denval>0) ? efficiency*sqrt(pow(numerr/numval,2)+pow(denerr/denval,2)) : 0;
               delete feff;
            } else if (TString(itpoi->name)=="lumi") {
               // luminosity and Ncoll
               if (isPP) {
                  itpoi->val = lumipp;
                  itpoi->err = 0.023; // from LUM-16-001
               } else {
                  if (thebin.centbin().low()>=60) itpoi->val = lumipbpb_peri;
                  else itpoi->val = lumipbpb_ABCD;
                  itpoi->err = 0; // FIXME
               }
            } else if (TString(itpoi->name)=="ncoll") {
               if (isPP) {
                  itpoi->val=1; 
                  itpoi->err=0;
               } else {
                  itpoi->val = HI::findNcollAverage(thebin.centbin().low(),thebin.centbin().high());
                  itpoi->err = 0; //FIXME
               }
            } else if (TString(itpoi->name)=="npart") {
               if (isPP) {
                  itpoi->val=2; 
                  itpoi->err=0;
               } else {
                  itpoi->val = HI::findNpartAverage(thebin.centbin().low(),thebin.centbin().high());
                  itpoi->err = 0; //FIXME
               }
            }
            if (TString(itpoi->name)=="taa") {
               itpoi->val=HI::findTaaAverage(thebin.centbin().low(),thebin.centbin().high());
               itpoi->err=HI::findTaaAverage_err(thebin.centbin().low(),thebin.centbin().high());
            }

         }

         // delete model;
         // delete model_bkg;
         // delete dat;
         delete ws;
         delete f;
      } else {
         for (vector<poi>::iterator itpoi=thePois.begin(); itpoi!=thePois.end(); itpoi++) {
            itpoi->val = 0;
            itpoi->err = 0;
         }
      }

      // fill the tree
      tr->Fill();
      cnt++;
   } // loop on the files

   f->Write();
   f->Close();
}

#endif // #ifndef results2tree_C
