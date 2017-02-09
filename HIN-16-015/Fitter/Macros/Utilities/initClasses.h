#ifndef initClasses_h
#define initClasses_h

#include "TSystem.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH1.h"
#include "TF1.h"
#include "TProfile.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TPave.h"
#include "TPad.h"
#include "TFrame.h"
#include "TAxis.h"
#include "TLegend.h"

#include "RooBinning.h"
#include "RooNumIntConfig.h"
#include "RooWorkspace.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"
#include "RooGaussModel.h"
#include "RooChi2Var.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooPlotable.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAbsPdf.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooExponential.h"

#include "../CMS/tdrstyle.C"
#include "../CMS/CMS_lumi.C"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>


using namespace std;
using namespace RooFit;

typedef struct StartEnd {
  int Start, End;
} StartEnd;

typedef struct MinMax {
  double Min, Max;
} MinMax;

typedef struct EvtPar {
  StartEnd RunNb;
  int      TriggerBit;
} EvtPar;

typedef struct DiMuonPar {
  MinMax ctau, ctauNRes, ctauRes, ctauErr, ctauTrue, M, Pt, AbsRap;
  string ctauCut;
} DiMuonPar;

typedef struct SiMuonPar {
  MinMax Pt, Eta;
} SiMuonPar;

typedef struct InputOpt {
  int        oniaMode;
  EvtPar     PbPb, Pbp, pp;
} InputOpt;

typedef struct KinCuts {
  StartEnd   Centrality;
  SiMuonPar  sMuon;
  DiMuonPar  dMuon;
} KinCuts;

bool isEqualKinCuts(struct KinCuts cutA, struct KinCuts cutB, bool isPbp) 
{
  bool cond = true;

  if (isPbp) {
    cond = cond && (cutA.Centrality.Start    == cutB.Centrality.Start);
    cond = cond && (cutA.Centrality.End      == cutB.Centrality.End);
  }
  cond = cond && (cutA.sMuon.Pt.Min        == cutB.sMuon.Pt.Min);
  cond = cond && (cutA.sMuon.Pt.Max        == cutB.sMuon.Pt.Max);
  cond = cond && (cutA.sMuon.Eta.Min       == cutB.sMuon.Eta.Min);
  cond = cond && (cutA.sMuon.Eta.Max       == cutB.sMuon.Eta.Max);

  cond = cond && (cutA.dMuon.ctau.Min      == cutB.dMuon.ctau.Min);
  cond = cond && (cutA.dMuon.ctau.Max      == cutB.dMuon.ctau.Max);
  cond = cond && (cutA.dMuon.ctauErr.Min   == cutB.dMuon.ctauErr.Min);
  cond = cond && (cutA.dMuon.ctauErr.Max   == cutB.dMuon.ctauErr.Max);
  cond = cond && (cutA.dMuon.ctauTrue.Min  == cutB.dMuon.ctauTrue.Min);
  cond = cond && (cutA.dMuon.ctauTrue.Max  == cutB.dMuon.ctauTrue.Max);
  cond = cond && (cutA.dMuon.ctauNRes.Min  == cutB.dMuon.ctauNRes.Min);
  cond = cond && (cutA.dMuon.ctauRes.Min   == cutB.dMuon.ctauRes.Min);
  cond = cond && (cutA.dMuon.ctauCut       == cutB.dMuon.ctauCut);
  cond = cond && (cutA.dMuon.M.Min         == cutB.dMuon.M.Min);
  cond = cond && (cutA.dMuon.M.Max         == cutB.dMuon.M.Max);
  cond = cond && (cutA.dMuon.Pt.Min        == cutB.dMuon.Pt.Min);
  cond = cond && (cutA.dMuon.Pt.Max        == cutB.dMuon.Pt.Max);
  cond = cond && (cutA.dMuon.AbsRap.Min    == cutB.dMuon.AbsRap.Min);
  cond = cond && (cutA.dMuon.AbsRap.Max    == cutB.dMuon.AbsRap.Max);

  return cond;
}


struct ParticleMass { double JPsi, Psi2S, Y1S, Y2S, Y3S, Z; };
ParticleMass Mass = {3.096, 3.686, 9.460, 10.023, 10.355, 91.188};


enum class MassModel 
{
    InvalidModel =0,
    SingleGaussian=1, 
    DoubleGaussian=2, 
    SingleCrystalBall=3, 
    DoubleCrystalBall=4, 
    GaussianAndCrystalBall=6, 
    Uniform=7, 
    Chebychev1=8, 
    Chebychev2=9, 
    Chebychev3=10, 
    Chebychev4=11,
    Chebychev5=12,
    Chebychev6=13,
    ExpChebychev1=14,
    ExpChebychev2=15,
    ExpChebychev3=16,
    ExpChebychev4=17,
    ExpChebychev5=18,
    ExpChebychev6=19,
    Exponential=20
};
map< string , MassModel > MassModelDictionary = {
  {"InvalidModel",            MassModel::InvalidModel},
  {"SingleGaussian",          MassModel::SingleGaussian},
  {"DoubleGaussian",          MassModel::DoubleGaussian},
  {"SingleCrystalBall",       MassModel::SingleCrystalBall},
  {"DoubleCrystalBall",       MassModel::DoubleCrystalBall},
  {"GaussianAndCrystalBall",  MassModel::GaussianAndCrystalBall},
  {"Uniform",                 MassModel::Uniform},
  {"Chebychev1",              MassModel::Chebychev1},
  {"Chebychev2",              MassModel::Chebychev2},
  {"Chebychev3",              MassModel::Chebychev3},
  {"Chebychev4",              MassModel::Chebychev4},
  {"Chebychev5",              MassModel::Chebychev5},
  {"Chebychev6",              MassModel::Chebychev6},
  {"ExpChebychev1",           MassModel::ExpChebychev1},
  {"ExpChebychev2",           MassModel::ExpChebychev2},
  {"ExpChebychev3",           MassModel::ExpChebychev3},
  {"ExpChebychev4",           MassModel::ExpChebychev4},
  {"ExpChebychev5",           MassModel::ExpChebychev5},
  {"ExpChebychev6",           MassModel::ExpChebychev6},
  {"Exponential",             MassModel::Exponential}
};


enum class CtauModel 
{     
    InvalidModel=0,
    QuadrupleGaussianResolution=1,
    TripleGaussianResolution=2,
    DoubleGaussianResolution=3, 
    SingleGaussianResolution=4,
    TripleDecay=5,
    DoubleSingleSidedDecay=6,
    SingleSidedDecay=7,
    Delta=8
};
map< string , CtauModel > CtauModelDictionary = {
  {"InvalidModel",             CtauModel::InvalidModel},
  {"QuadrupleGaussianResolution", CtauModel::QuadrupleGaussianResolution},
  {"TripleGaussianResolution", CtauModel::TripleGaussianResolution},
  {"DoubleGaussianResolution", CtauModel::DoubleGaussianResolution},
  {"SingleGaussianResolution", CtauModel::SingleGaussianResolution},
  {"TripleDecay",              CtauModel::TripleDecay},
  {"DoubleSingleSidedDecay",   CtauModel::DoubleSingleSidedDecay},
  {"SingleSidedDecay",         CtauModel::SingleSidedDecay},
  {"Delta",                    CtauModel::Delta}
};

typedef struct CtauPNP {
  CtauModel    Prompt, NonPrompt;
} CtauPNP;

typedef struct CtauMassModel {
  MassModel  Mass;
  CtauPNP    Ctau;
} FitModel;

typedef struct CharmModel {
  CtauMassModel  Jpsi, Psi2S, Bkg;
  CtauModel CtauRes;
  CtauModel CtauTrue, CtauTrueRes;
} CharmModel;

typedef struct OniaModel {
  CharmModel  PbPb, Pbp, PP;
} OniaModel;




#include <boost/algorithm/string/replace.hpp>

void setFixedVarsToContantVars(RooWorkspace& ws)
{
  RooArgSet listVar = ws.allVars();
  TIterator* parIt = listVar.createIterator();
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    if ( it->getMin()==it->getMax() && !it->isConstant() ) {
      cout << "[INFO] Setting " << it->GetName() << " constant!" << endl;
      it->setConstant(kTRUE);
    }
  }
};


bool compareSnapshots(RooArgSet *pars1, const RooArgSet *pars2) {
  TIterator* parIt = pars1->createIterator(); 
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    double val = pars2->getRealValue(it->GetName(),-1e99);
    if ( strcmp(it->GetName(),"ctauErr")==0 || strcmp(it->GetName(),"ctau")==0 || strcmp(it->GetName(),"ctauTrue")==0 || strcmp(it->GetName(),"ctauRes")==0 || strcmp(it->GetName(),"ctauNRes")==0 ) continue;
    if (val==-1e99) return false;          // the parameter was not found!
    if (val != it->getVal()) return false;  // the parameter was found, but with a different value!
    if ( ((RooRealVar&)(*pars2)[it->GetName()]).getMin() != it->getMin() ) return false;  // the parameter has different lower limit
    if ( ((RooRealVar&)(*pars2)[it->GetName()]).getMax() != it->getMax() ) return false;  // the parameter has different upper limit
  }
  return true;
};


bool saveWorkSpace(RooWorkspace& myws, string outputDir, string FileName)
{
  // Save the workspace
  gSystem->mkdir(outputDir.c_str(), kTRUE);
  cout << FileName << endl;
  TFile *file =  new TFile(FileName.c_str(), "RECREATE");
  if (!file) {
    file->Close(); delete file;
    cout << "[ERROR] Output root file with fit results could not be created!" << endl; return false; 
  } else {
    file->cd();    
    myws.Write("workspace"); 
    file->Write(); file->Close(); delete file;
  }
  return true;
};


bool isFitAlreadyFound(RooArgSet *newpars, string FileName, string pdfName) 
{
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[INFO] FileName: " << FileName << " was not found" << endl;
    return false; // File was not found
  }
  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;  
  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    cout << "[INFO] Workspace was not found" << endl;
    file->Close(); delete file;
    return false;
  }
  string snapShotName = Form("%s_parIni", pdfName.c_str());
  const RooArgSet *params = ws->getSnapshot(snapShotName.c_str());
  if (!params) {
    cout << "[INFO] Snapshot parIni was not found" << endl;
    delete ws;
    file->Close(); delete file;
    return false;
  }
  bool result = compareSnapshots(newpars, params);
  delete ws;
  file->Close(); delete file; 

  return result;
};


bool loadPreviousFitResult(RooWorkspace& myws, string FileName, string DSTAG, bool isPbp, bool cutSideBand=false)
{
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[INFO] File " << FileName << " was not found!" << endl;
    return false; // File was not found
  }
  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;

  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    cout << "[INFO] Workspace was not found in: " << FileName << endl;
    file->Close(); delete file;
    return false;
  }

  if (DSTAG.find("_")!=std::string::npos) DSTAG.erase(DSTAG.find("_"));

  cout <<  "[INFO] Loading variables and functions from: " << FileName << endl;
  RooArgSet listVar = ws->allVars();
  TIterator* parIt = listVar.createIterator();
  string print = "[INFO] Variables loaded: ";
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    string name = it->GetName();
    if (myws.var(name.c_str()) &&  myws.var(name.c_str())->isConstant()) continue;
    if ( name=="invMass" || name=="ctau" || name=="ctauErr" || name=="ctauRes" || name=="ctauNRes" ||
         name=="ctauTrue" || name=="pt" || name=="cent" || name=="rap" || name=="One" ||
	 name=="Ntracks" || name=="SumET_HFEta4" || name=="Zvtx") continue;
    if ( (DSTAG.find("MC")!=std::string::npos || cutSideBand) && (name.find("N_")!=std::string::npos) ) continue; 
    if (myws.var(name.c_str())) {
      print = print + Form("  %s: %.5f->%.5f ", name.c_str(), myws.var(name.c_str())->getValV(), ws->var(name.c_str())->getValV()) ;
      myws.var(name.c_str())->setVal  ( ws->var(name.c_str())->getValV()  );
      myws.var(name.c_str())->setError( ws->var(name.c_str())->getError() );
    } else {
      if ( (name==Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))) && myws.var(Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))) ) {
        print = print + Form("  %s: %.5f->%.5f", Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), myws.var(Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))->getValV(), ws->var(name.c_str())->getValV()) ;
        myws.var(Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))->setVal  ( ws->var(name.c_str())->getValV()  );
        myws.var(Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))->setError( ws->var(name.c_str())->getError() ); 
      }
      if ( (name==Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))) && myws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))) {
        print = print + Form("  %s: %.5f->%.5f  ", Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), myws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))->getValV(), ws->var(name.c_str())->getValV()) ;
        myws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))->setVal  ( ws->var(name.c_str())->getValV()  );
        myws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))->setError( ws->var(name.c_str())->getError() ); 
      }
    }
  }
  cout << print << endl;
  RooArgSet listFun = ws->allFunctions();
  TIterator* parFunIt = listFun.createIterator();
  string printFun = "[INFO] Functions loaded: ";
  for (RooRealVar* it = (RooRealVar*)parFunIt->Next(); it!=NULL; it = (RooRealVar*)parFunIt->Next() ) {
    string name = it->GetName();
    if (myws.var(name.c_str()) &&  myws.var(name.c_str())->isConstant()) continue;
    if ( name=="invMass" || name=="ctau" || name=="ctauErr" || name=="ctauRes" || name=="ctauNRes" ||
         name=="ctauTrue" || name=="pt" || name=="cent" || name=="rap" || name=="One" ||
	 name=="Ntracks" || name=="SumET_HFEta4" || name=="Zvtx") continue;
    if ( (DSTAG.find("MC")!=std::string::npos || cutSideBand) && (name.find("N_")!=std::string::npos) ) continue; 
    if (myws.var(name.c_str())) { 
      printFun = printFun + Form("  %s: %.5f->%.5f  ", name.c_str(), myws.var(name.c_str())->getValV(), ws->function(name.c_str())->getValV()) ;
      myws.var(name.c_str())->setVal  ( ws->function(name.c_str())->getValV()  );
      myws.var(name.c_str())->setError( 0.0 );
    } else {
      Double_t MassRatio = (Mass.Psi2S/Mass.JPsi);
      string reName = name.c_str();
      boost::replace_all(reName, "Jpsi", "Psi2S");
      if (myws.var(reName.c_str())) {
        Double_t value = 0.0;
        if ( (reName==Form("sigma1_Psi2S_%s", (isPbp?"Pbp":"PP"))) ) { value = ws->var(name.c_str())->getValV() * MassRatio; }
        else if ( (reName==Form("m_Psi2S_%s", (isPbp?"Pbp":"PP"))) ) { value = ws->var(name.c_str())->getValV() * MassRatio; }
        else { value = ws->var(name.c_str())->getValV(); }
        printFun = printFun + Form("  %s: %.5f->%.5f  ", reName.c_str(), myws.var(reName.c_str())->getValV(), value) ;
        myws.var(reName.c_str())->setVal  ( value );
        myws.var(reName.c_str())->setError( 0.0 );
      }
    }
  }
  setFixedVarsToContantVars(myws);
  cout << printFun << endl;

  delete ws;
  file->Close(); delete file;
  return true;
};


bool loadCtauErrRange(RooWorkspace& myws, string FileName, struct KinCuts& cut)
{
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[INFO] File " << FileName << " was not found!" << endl;
    return false; // File was not found
  }
  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;
  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    cout << "[INFO] Workspace was not found in: " << FileName << endl;
    file->Close(); delete file;
    return false;
  }

  if (ws->var("ctauErr")) {
    cut.dMuon.ctauErr.Min = ws->var("ctauErr")->getMin();
    cut.dMuon.ctauErr.Max = ws->var("ctauErr")->getMax();
  } else { 
    cout << Form("[WARNING] ctauErr was not found!") << endl;
    delete ws;
    file->Close(); delete file;
    return false;
  }

  delete ws;
  file->Close(); delete file;
  return true;
};


int importDataset(RooWorkspace& myws, RooWorkspace& inputWS, struct KinCuts cut, string label, bool cutSideBand=false)
{
  string indMuonMass    = Form("(%.6f < invMass && invMass < %.6f)",       cut.dMuon.M.Min,       cut.dMuon.M.Max);
  if (cutSideBand) {
    //indMuonMass =  indMuonMass + "&&" + "((2.0 < invMass && invMass < 2.8) || (3.3 < invMass && invMass < 3.5) || (3.9 < invMass && invMass < 5.0))";
    indMuonMass = indMuonMass + "&&" +  "((3.35 <= invMass && invMass <= 3.45) || (3.85 < invMass && invMass < 4.2) )";
  }
  //string indMuonRap     = Form("(%.6f <= abs(rap) && abs(rap) < %.6f)",    cut.dMuon.AbsRap.Min,   cut.dMuon.AbsRap.Max);
  string indMuonRap     = Form("(%.6f <= rap && rap < %.6f)",    cut.dMuon.AbsRap.Min,  cut.dMuon.AbsRap.Max);
  string indMuonPt      = Form("(%.6f <= pt && pt < %.6f)",                cut.dMuon.Pt.Min,       cut.dMuon.Pt.Max);
  string indMuonCtau    = Form("(%.6f < ctau && ctau < %.6f)",             cut.dMuon.ctau.Min,     cut.dMuon.ctau.Max); 
  if(cut.dMuon.ctauCut!=""){ indMuonCtau = cut.dMuon.ctauCut; }
  string indMuonCtauErr = Form("(%.12f < ctauErr && ctauErr < %.12f)",     cut.dMuon.ctauErr.Min,  cut.dMuon.ctauErr.Max);
  string inCentrality   = Form("(%d <= cent && cent < %d)",                cut.Centrality.Start,   cut.Centrality.End);
  string indMuonCtauTrue = Form("(%.12f < ctauTrue && ctauTrue < %.12f)",  cut.dMuon.ctauTrue.Min, cut.dMuon.ctauTrue.Max);
  string indMuonCtauRes = Form("(%.12f < ctauRes && ctauRes < %.12f)",     cut.dMuon.ctauRes.Min,  cut.dMuon.ctauRes.Max);
  string indMuonCtauNRes = Form("(%.12f < ctauNRes && ctauNRes < %.12f)",  cut.dMuon.ctauNRes.Min, cut.dMuon.ctauNRes.Max);
  string strCut         = indMuonMass +"&&"+ indMuonRap +"&&"+ indMuonPt +"&&"+ indMuonCtau +"&&"+ indMuonCtauErr;
  //if (label.find("Pbp")!=std::string::npos){ strCut = strCut +"&&"+ inCentrality; }
  if (label.find("MC")!=std::string::npos){ strCut = strCut +"&&"+ indMuonCtauTrue +"&&"+ indMuonCtauNRes +"&&"+ indMuonCtauRes;  }

  // Reduce and import the datasets
  if (!(inputWS.data(Form("dOS_%s", label.c_str())))){ 
    cout << "[ERROR] The dataset " <<  Form("dOS_%s", label.c_str()) << " was not found!" << endl;
    return -1;
  }
  cout << "[INFO] Importing local RooDataSet with cuts: " << strCut << endl;
  RooDataSet* dataOS = (RooDataSet*)inputWS.data(Form("dOS_%s", label.c_str()))->reduce(strCut.c_str());
  if (dataOS->sumEntries()==0){ 
    cout << "[WARNING] No events from dataset " <<  Form("dOS_%s", label.c_str()) << " passed the kinematic cuts!" << endl;
    return 0;
  }  
  myws.import(*dataOS);
  delete dataOS;
  
  if (label.find("NoBkg")==std::string::npos && label.find("AccEff")==std::string::npos && label.find("lJpsiEff")==std::string::npos) // Don't try to find SS dataset if label contais NoBkg or correction
  {
    if (!(inputWS.data(Form("dSS_%s", label.c_str())))){
      cout << "[ERROR] The dataset " <<  Form("dSS_%s", label.c_str()) << " was not found!" << endl;
      return -1;
    }
    RooDataSet* dataSS = (RooDataSet*)inputWS.data(Form("dSS_%s", label.c_str()))->reduce(strCut.c_str());
    if (dataSS->sumEntries()==0){
      cout << "[WARNING] No events from dataset " <<  Form("dSS_%s", label.c_str()) << " passed the kinematic cuts!" << endl;
    }
    myws.import(*dataSS);
    delete dataSS;
  }

  // Set the range of each global parameter in the local roodataset
  const RooArgSet* rowOS = myws.data(Form("dOS_%s", label.c_str()))->get();
  ((RooRealVar*)rowOS->find("invMass"))->setMin(cut.dMuon.M.Min);        
  ((RooRealVar*)rowOS->find("invMass"))->setMax(cut.dMuon.M.Max);
  ((RooRealVar*)rowOS->find("pt"))->setMin(cut.dMuon.Pt.Min);            
  ((RooRealVar*)rowOS->find("pt"))->setMax(cut.dMuon.Pt.Max);
  ((RooRealVar*)rowOS->find("ctau"))->setMin(cut.dMuon.ctau.Min);        
  ((RooRealVar*)rowOS->find("ctau"))->setMax(cut.dMuon.ctau.Max);
  ((RooRealVar*)rowOS->find("ctauErr"))->setMin(cut.dMuon.ctauErr.Min);
  ((RooRealVar*)rowOS->find("ctauErr"))->setMax(cut.dMuon.ctauErr.Max);
  //if (label.find("Pbp")!=std::string::npos){
  //((RooRealVar*)rowOS->find("cent"))->setMin(cut.Centrality.Start);      
  //((RooRealVar*)rowOS->find("cent"))->setMax(cut.Centrality.End);
  //}
  if (label.find("MC")!=std::string::npos){
    ((RooRealVar*)rowOS->find("ctauTrue"))->setMin(cut.dMuon.ctauTrue.Min);      
    ((RooRealVar*)rowOS->find("ctauTrue"))->setMax(cut.dMuon.ctauTrue.Max);
    ((RooRealVar*)rowOS->find("ctauRes"))->setMin(cut.dMuon.ctauRes.Min);      
    ((RooRealVar*)rowOS->find("ctauRes"))->setMax(cut.dMuon.ctauRes.Max);
    ((RooRealVar*)rowOS->find("ctauNRes"))->setMin(cut.dMuon.ctauNRes.Min);      
    ((RooRealVar*)rowOS->find("ctauNRes"))->setMax(cut.dMuon.ctauNRes.Max);
  }
  // Set the range of each global parameter in the local workspace
  myws.var("invMass")->setMin(cut.dMuon.M.Min);        
  myws.var("invMass")->setMax(cut.dMuon.M.Max);
  myws.var("pt")->setMin(cut.dMuon.Pt.Min);            
  myws.var("pt")->setMax(cut.dMuon.Pt.Max);
  myws.var("rap")->setMin(cut.dMuon.AbsRap.Min);       
  myws.var("rap")->setMax(cut.dMuon.AbsRap.Max);
  myws.var("ctau")->setMin(cut.dMuon.ctau.Min);        
  myws.var("ctau")->setMax(cut.dMuon.ctau.Max);
  myws.var("ctauErr")->setMin(cut.dMuon.ctauErr.Min);  
  myws.var("ctauErr")->setMax(cut.dMuon.ctauErr.Max);
  //if (label.find("Pbp")!=std::string::npos){
  //myws.var("cent")->setMin(cut.Centrality.Start);      
  //myws.var("cent")->setMax(cut.Centrality.End);
  //}
  if (label.find("MC")!=std::string::npos){
    myws.var("ctauTrue")->setMin(cut.dMuon.ctauTrue.Min);      
    myws.var("ctauTrue")->setMax(cut.dMuon.ctauTrue.Max);
    myws.var("ctauRes")->setMin(cut.dMuon.ctauRes.Min);      
    myws.var("ctauRes")->setMax(cut.dMuon.ctauRes.Max);
    myws.var("ctauNRes")->setMin(cut.dMuon.ctauNRes.Min);      
    myws.var("ctauNRes")->setMax(cut.dMuon.ctauNRes.Max);
  }
  cout << "[INFO] Analyzing bin: " << Form(
                                           "%.3f < pt < %.3f, %.3f < rap < %.3f", 
                                           cut.dMuon.Pt.Min,
                                           cut.dMuon.Pt.Max,
                                           cut.dMuon.AbsRap.Min,
                                           cut.dMuon.AbsRap.Max
					   ) << endl;

  return 1;
};


void printChi2(RooWorkspace& myws, TPad* Pad, RooPlot* frame, string varLabel, string dataLabel, string pdfLabel, int nBins, bool useDefaultName=true)
{
  double chi2=0; unsigned int ndof=0;
  Pad->cd();
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.1); 
  unsigned int nFitPar = myws.pdf(pdfLabel.c_str())->getParameters(*myws.data(dataLabel.c_str()))->selectByAttrib("Constant",kFALSE)->getSize();
  TH1* hdatact = myws.data(dataLabel.c_str())->createHistogram("hdatact", *myws.var(varLabel.c_str()), Binning(frame->GetNbinsX(),frame->GetXaxis()->GetXmin(),frame->GetXaxis()->GetXmax()));
//  RooHist *hpull = frame->pullHist("hdatact",pdfLabel.c_str(), true);
  RooHist *hpull = frame->pullHist(0,0, true);
  double* ypulls = hpull->GetY();
  unsigned int nFullBins = 0;
  for (int i = 0; i < nBins; i++) {
    if (hdatact->GetBinContent(i+1) > 0.0) {
      chi2 += ypulls[i]*ypulls[i];
      nFullBins++;
    }
  }
  ndof = nFullBins - nFitPar;
  t->DrawLatex(0.7, 0.85, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));
  if (useDefaultName) {
    RooRealVar chi2Var("chi2","chi2",chi2);
    RooRealVar ndofVar("ndof","ndof",ndof);
    myws.import(chi2Var); myws.import(ndofVar);
  } else {
    RooRealVar chi2Var((string("chi2_")+varLabel).c_str(),(string("chi2_")+varLabel).c_str(),chi2);
    RooRealVar ndofVar((string("ndof_")+varLabel).c_str(),(string("ndof_")+varLabel).c_str(),ndof);
    myws.import(chi2Var); myws.import(ndofVar);
  }
  delete hdatact; 
  delete hpull;
};


void getCtauErrRange(TH1* hist, int nMaxBins, vector<double>& rangeErr)
{
  // 1) Find the bin with the maximum Y value
  int binMaximum = hist->GetMaximumBin();
  // 2) Loop backward and find the first bin
  int binWithContent = -1;
  int firstBin = 1;
  for (int i=binMaximum; i>0; i--) {
    if (hist->GetBinContent(i)>0.0) {
      if ( binWithContent>0 && ((binWithContent-i) > nMaxBins) && hist->GetBinContent(i)<5.0 ) { firstBin = binWithContent; break; }
      else { binWithContent = i; }
    }
  }
  // 3) Loop forward and find the last bin
  binWithContent = -1;
  int lastBin = hist->GetNbinsX();
  for (int i=binMaximum; i<hist->GetNbinsX(); i++) {
    if (hist->GetBinContent(i)>0.0) {
      if ( binWithContent>0 && ((i - binWithContent) > nMaxBins) && hist->GetBinContent(i)<5.0 ) { lastBin = binWithContent+1; break; }
      else { binWithContent = i; }
    }
  }
  // 4) Build the set of bins
  int startBin = ( (firstBin>1) ? (firstBin-1) : firstBin );
  const int nNewBins = lastBin - startBin + 1;
  double binning[nNewBins+1];
  binning[0] = hist->GetXaxis()->GetXmin();
  binning[nNewBins] = hist->GetXaxis()->GetXmax();
  for (int i=1; i<nNewBins; i++) {
    int iBin = startBin + i;
    binning[i] = hist->GetBinLowEdge(iBin); 
  }
  rangeErr.push_back(binning[(firstBin>1)?1:0]);
  rangeErr.push_back(binning[nNewBins-1]);

  cout << "[INFO] Ctau error range set to be: [ " <<  rangeErr[0] << ", " << rangeErr[1] << " ]" << endl;

  return;
};


bool setConstant( RooWorkspace& myws, string parName, bool CONST)
{
  if (myws.var(parName.c_str())) { 
    myws.var(parName.c_str())->setConstant(CONST);
    if (CONST) { cout << "[INFO] Setting parameter " << parName << " : " << myws.var(parName.c_str())->getVal() << " to constant value!" << endl; }
  }
  else if (!myws.function(parName.c_str())) { 
    cout << "[ERROR] Parameter " << parName << " was not found!" << endl;
    return false;
  }

  return true;
};


#endif // #ifndef initClasses_h
