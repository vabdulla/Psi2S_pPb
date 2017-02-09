#ifndef fitCharmoniaCtauModel_C
#define fitCharmoniaCtauModel_C

#include "Utilities/initClasses.h"
#include "buildCharmoniaCtauModel.C"
#include "fitCharmoniaMassModel.C"
#include "fitCharmoniaCtauResModel.C"
#include "fitCharmoniaCtauErrModel.C"
#include "drawCtauPlot.C"

bool setCtauModel( struct OniaModel& model, map<string, string>&  parIni, bool isPbp, bool incJpsi, bool incPsi2S, bool incBkg, bool incPrompt, bool incNonPrompt );
void setCtauFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbp, bool fitSideBand);
void setCtauGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label, double binWidth, bool fitCtauRes=false);
void setCtauCutParameters(struct KinCuts& cut, bool incNonPrompt);


bool fitCharmoniaCtauModel( RooWorkspace& myws,             // Local Workspace
                            RooWorkspace& inputWorkspace,   // Workspace with all the input RooDatasets
                            struct KinCuts& cut,            // Variable containing all kinematic cuts
                            map<string, string>&  parIni,   // Variable containing all initial parameters
                            struct InputOpt& opt,           // Variable with run information (kept for legacy purpose)
                            string outputDir,               // Path to output directory
                            // Select the type of datasets to fit
                            string DSTAG,                   // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
                            bool isPbp        = false,     // isPbp = false for pp, true for Pbp
                            bool importDS      = true,      // Select if the dataset is imported in the local workspace
                            // Select the type of object to fit
                            bool incJpsi       = true,      // Includes Jpsi model
                            bool incPsi2S      = false,     // Includes Psi(2S) model
                            bool incBkg        = true,      // Includes Background model
                            bool incPrompt     = true,      // Includes Prompt ctau model
                            bool incNonPrompt  = false,     // Includes NonPrompt ctau model
                            // Select the fitting options
                            bool doFit         = true,      // Flag to indicate if we want to perform the fit
                            bool wantPureSMC   = false,     // Flag to indicate if we want to fit pure signal MC
                            bool loadFitResult = false,     // Load previous fit results
                            map<string, string> inputFitDir={},// User-defined Location of the fit results
                            int  numCores      = 2,         // Number of cores used for fitting
                            // Select the drawing options
                            bool setLogScale   = true,      // Draw plot with log scale
                            bool incSS         = false,     // Include Same Sign data
                            map<string, double> binWidth={} // User-defined Location of the fit results
                            )  
{

  bool usePerEventError = true;
  
  if (DSTAG.find("_")!=std::string::npos) DSTAG.erase(DSTAG.find("_"));

  // Check if input dataset is MC
  bool isMC = false;
  if (DSTAG.find("MC")!=std::string::npos) {
    if (incJpsi && incPsi2S) { 
      cout << "[ERROR] We can only fit one type of signal using MC" << endl; return false; 
    }
    isMC = true;
  }
  if (!isMC) { wantPureSMC=false; }
  bool fitSideBand = (incBkg && (!incPsi2S && !incJpsi));

  string COLL = (isPbp ? "Pbp" : "PP" );
  string fitType = "CTAU";
  if (isMC && incPrompt && !incNonPrompt) { fitType = "CTAURES"; }
  if (!isMC && fitSideBand) { fitType = "CTAUSB";  }

  if (importDS) {
    setCtauCutParameters(cut, incNonPrompt);
    setMassCutParameters(cut, incJpsi, incPsi2S, isMC, true);
    if (usePerEventError) {
      // check if we have already done the ctauErr fits. If yes, load their parameters
      string FileName = "";
      string pdfName = Form("pdfCTAUERR_Tot_%s", COLL.c_str());
      string plotLabel = "";
      //bool incJpsi = true;
      bool incPsi2S = true;
      bool fitSideBand = false;
      if (incJpsi)  { plotLabel = plotLabel + "_Jpsi";     }
      if (incPsi2S) { plotLabel = plotLabel + "_Psi2S";    }
      plotLabel = plotLabel + "_Bkg";
      setCtauErrFileName(FileName, (inputFitDir["CTAUERR"]=="" ? outputDir : inputFitDir["CTAUERR"]), "DATA", plotLabel, cut, isPbp, fitSideBand);
      bool foundFit = false;
      if ( loadCtauErrRange(myws, FileName, cut) ) { foundFit = true; }
      if (foundFit) { cout << "[INFO] The ctauErr fit was found and I'll load the ctau Error range used." << endl; }
    }
    setCtauErrCutParameters(cut);
  }

  // Import the local datasets
  double numEntries = 1000000;
  string label = ((DSTAG.find(COLL.c_str())!=std::string::npos) ? DSTAG.c_str() : Form("%s_%s", DSTAG.c_str(), COLL.c_str()));
  if (wantPureSMC) label = Form("%s_NoBkg", label.c_str());
  string dsName = Form("dOS_%s", label.c_str());
  if (importDS) {
    if ( !(myws.data(dsName.c_str())) ) {
      int importID = importDataset(myws, inputWorkspace, cut, label, fitSideBand);
      if (importID<0) { return false; }
      else if (importID==0) { doFit = false; }
    }
    numEntries = myws.data(dsName.c_str())->sumEntries(); if (numEntries<=0) { doFit = false; }
  }
  else if (doFit && !(myws.data(dsName.c_str()))) { cout << "[ERROR] No local dataset was found to perform the fit!" << endl; return false; }

  if (importDS) { 
    // Set global parameters
    setCtauErrGlobalParameterRange(myws, parIni, cut, "", binWidth["CTAUERR"], true);
    setCtauGlobalParameterRange(myws, parIni, cut, label, binWidth[fitType.c_str()], (fitType=="CTAURES"));
    
    RooDataSet* dataToFit = (RooDataSet*)(myws.data(dsName.c_str())->reduce(parIni["CtauRange_Cut"].c_str()))->Clone((dsName+"_CTAUCUT").c_str());
    myws.import(*dataToFit);
  }
  string dsNameCut = dsName+"_CTAUCUT";

  // Define pdf and plot names
  string pdfName = Form("pdfCTAU_Tot_%s", COLL.c_str());
  
  if (!loadFitResult) {
    // Set models based on initial parameters
    struct OniaModel model;
    if (!setCtauModel(model, parIni, isPbp, incJpsi, incPsi2S, incBkg, incPrompt, incNonPrompt)) { return false; }
    //// LOAD CTAU ERROR PDF
    if (usePerEventError) {
      // Setting extra input information needed by each fitter
      bool loadCtauErrFitResult = true;
      bool doCtauErrFit = true;
      bool importDS = isMC;
      //bool incJpsi = true;
      bool incPsi2S = true;
      string DSTAG = Form("DATA_%s", (isPbp?"Pbp":"PP"));
  
      if ( !fitCharmoniaCtauErrModel( myws, inputWorkspace, cut, parIni, opt, outputDir, 
                                      DSTAG, isPbp, importDS, 
                                      incJpsi, incPsi2S, incBkg, 
                                      doCtauErrFit, wantPureSMC, loadCtauErrFitResult, inputFitDir, numCores, 
                                      setLogScale, incSS, binWidth
                                      ) 
           ) { return false; }
    }
    // Build the Fit Model
    if (!buildCharmoniaCtauModel(myws, (isPbp ? model.Pbp : model.PP), parIni, dsName, isPbp, incBkg, incJpsi, incPsi2S, incPrompt, incNonPrompt,  numEntries))  { return false; }

    //// LOAD CTAU RESOLUTION PDF
    if (fitSideBand) {
      cout<<"fitSideBand :  "<<fitSideBand<<endl;
      
      // check if we have already done the resolution fits. If yes, load their results
      string FileName = "";
      string plotLabel = Form("_CtauRes_%s", parIni[Form("Model_CtauRes_%s", COLL.c_str())].c_str());
      //string DSTAG = Form("MCJPSIPR_%s", (isPbp?"Pbp":"PP"));
      string DSTAG = Form("MCPSI2SPR_%s", (isPbp?"Pbp":"PP"));
      setCtauResFileName(FileName, (inputFitDir["CTAURES"]=="" ? outputDir : inputFitDir["CTAURES"]), DSTAG, plotLabel, cut, isPbp);
      if (wantPureSMC) { plotLabel = plotLabel + "_NoBkg"; }
      bool found = false;
      if (!found && gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAURES"]!="") {
        plotLabel = string(Form("_CtauRes_%s_NoBkg", parIni[Form("Model_CtauRes_%s", COLL.c_str())].c_str()));
        setCtauResFileName(FileName, (inputFitDir["CTAURES"]=="" ? outputDir : inputFitDir["CTAURES"]), DSTAG, plotLabel, cut, isPbp);
      } else if (inputFitDir["CTAURES"]!="") { found = true; }
      if (!found && gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAURES"]!="") {
        plotLabel = Form("_CtauRes_%s", parIni[Form("Model_CtauRes_%s", COLL.c_str())].c_str());
        setCtauResFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbp);
      } else if (inputFitDir["CTAURES"]!="") { found = true; }
      if (!found && gSystem->AccessPathName(FileName.c_str())) {
        plotLabel = string(Form("_CtauRes_%s_NoBkg", parIni[Form("Model_CtauRes_%s", COLL.c_str())].c_str()));
        setCtauResFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbp);
      } else { found = true; }
      if (!found && gSystem->AccessPathName(FileName.c_str())) {
        cout << "[ERROR] User Input File : " << FileName << " was not found!" << endl;
        return false;
      }
      if ( !loadPreviousFitResult(myws, FileName, DSTAG, isPbp) ) {
        cout << "[ERROR] The ctau resolution fit results were not loaded!" << endl;
        return false;
      } else { 
        cout << "[INFO] The ctau resolution fits were found, so I'll load the fit results." << endl; 
      }
      if (myws.pdf(Form("pdfCTAU_BkgPR_%s", (isPbp?"Pbp":"PP")))) {
        cout << "[INFO] Setting Prompt Background parameters to constant!" << endl;
        myws.pdf(Form("pdfCTAU_BkgPR_%s", (isPbp?"Pbp":"PP")))->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE); 
      } else { cout << "[ERROR] Prompt Background PDF was not found!" << endl; return false; }
      if (!setConstant(myws, Form("s1_CtauRes_%s", COLL.c_str()), false)) { return false; }
    }

    // save the initial values of the model we've just created
    RooArgSet* params = (RooArgSet*) myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("invMass"), *myws.var("ctauErr")));
    myws.saveSnapshot((pdfName+"_parIni").c_str(),*params,kTRUE);
    delete params;
  }

  string plotLabel = "";
  map<string, bool> plotLabels = {{"JpsiNoPR", (incJpsi&&incNonPrompt)}, 
                                  {"Psi2SNoPR", (incPsi2S&&incNonPrompt)}, 
                                  {"BkgNoPR", (incBkg&&incNonPrompt)}, 
                                  {"CtauRes", (true)}};
  for (map<string, bool>::iterator iter = plotLabels.begin(); iter!=plotLabels.end(); iter++) {
    string obj = iter->first;
    bool cond = iter->second;
    if (cond && parIni.count(Form("Model_%s_%s", obj.c_str(), COLL.c_str()))>0) { 
      plotLabel = plotLabel + Form("_%s_%s", obj.c_str(), parIni[Form("Model_%s_%s", obj.c_str(), COLL.c_str())].c_str()); 
    }
  }
  if (wantPureSMC) { plotLabel = plotLabel + "_NoBkg"; }

  // check if we have already done this fit. If yes, do nothing and return true.
  string FileName = "";
  setCtauFileName(FileName, (inputFitDir[fitType.c_str()]=="" ? outputDir : inputFitDir[fitType.c_str()]), DSTAG, plotLabel, cut, isPbp, fitSideBand);
  if (gSystem->AccessPathName(FileName.c_str()) && inputFitDir[fitType.c_str()]!="") {
    cout << "[WARNING] User Input File : " << FileName << " was not found!" << endl;
    cout<<"===================okfitCtau0========================="<<endl;
    if (loadFitResult) return false;
    cout<<"===================okfitCtau1========================="<<endl;
    setCtauFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbp, fitSideBand);
  }
  cout<<"===================okfitCtau2========================="<<endl;
  bool found =  true; bool skipFit = !doFit;
  RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr")));
  found = found && isFitAlreadyFound(newpars, FileName, pdfName.c_str());
  if (loadFitResult) {
    if ( loadPreviousFitResult(myws, FileName, DSTAG, isPbp) ) { skipFit = true; } else  { skipFit = false; }
    if (skipFit) { cout << "[INFO] This ctau fit was already done, so I'll load the fit results." << endl; }
    myws.saveSnapshot(Form("%s_parLoad", pdfName.c_str()),*newpars,kTRUE);
  } else if (found) {
    cout << "[INFO] This ctau fit was already done, so I'll just go to the next one." << endl;
    return true;
  }

  cout<<"===================okfitCtau3========================="<<endl;

  // Fit the Datasets
  if (skipFit==false) {
    bool isWeighted = myws.data(dsName.c_str())->isWeighted();
    RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsNameCut.c_str()), Extended(kTRUE), NumCPU(numCores), ConditionalObservables(*myws.var("ctauErr")), SumW2Error(isWeighted), Save());
    fitResult->Print("v");
    myws.import(*fitResult, Form("fitResult_%s", pdfName.c_str()));
    // Draw the mass plot
    drawCtauPlot(myws, outputDir, opt, cut, parIni, plotLabel, DSTAG, isPbp, incJpsi, incPsi2S, incBkg, incPrompt, incNonPrompt, wantPureSMC, setLogScale, incSS, binWidth[fitType.c_str()]);
    // Save the results
    string FileName = ""; setCtauFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbp, fitSideBand);
    myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE);
    saveWorkSpace(myws, Form("%sctau%s/%s/result", outputDir.c_str(), (fitSideBand?"SB":""), DSTAG.c_str()), FileName);
  }
  
  return true;
};


bool setCtauModel( struct OniaModel& model, map<string, string>&  parIni, bool isPbp, bool incJpsi, bool incPsi2S, bool incBkg, bool incPrompt, bool incNonPrompt )
{
  if (isPbp) {
    if (parIni.count("Model_CtauRes_Pbp")>0) {
      model.Pbp.CtauRes = CtauModelDictionary[parIni["Model_CtauRes_Pbp"]];
      if (model.Pbp.CtauRes==CtauModel(0)) {
        cout << "[ERROR] The ctau resolution model: " << parIni["Model_CtauRes_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Ctau Resolution model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  } else {
    if (parIni.count("Model_CtauRes_PP")>0) {
      model.PP.CtauRes = CtauModelDictionary[parIni["Model_CtauRes_PP"]];
      if (model.PP.CtauRes==CtauModel(0)) {
        cout << "[ERROR] The ctau resolution model: " << parIni["Model_CtauRes_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Ctau Resolution model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (isPbp && incBkg && incNonPrompt) {
    if (parIni.count("Model_BkgNoPR_Pbp")>0) {
      model.Pbp.Bkg.Ctau.NonPrompt = CtauModelDictionary[parIni["Model_BkgNoPR_Pbp"]];
      if (model.Pbp.Bkg.Ctau.NonPrompt==CtauModel(0)) {
        cout << "[ERROR] The background non-prompt ctau model: " << parIni["Model_BkgNoPR_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Background non-prompt ctau model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  }
  if (isPbp && incBkg && incPrompt) {
    if (parIni.count("Model_BkgPR_Pbp")>0) {
      model.Pbp.Bkg.Ctau.Prompt = CtauModelDictionary[parIni["Model_BkgPR_Pbp"]];
      if (model.Pbp.Bkg.Ctau.Prompt==CtauModel(0)) {
        cout << "[ERROR] The background prompt ctau model: " << parIni["Model_BkgPR_Pbp"] << " is invalid" << endl; return false;
      }
    } else {
      parIni["Model_BkgPR_Pbp"] = "Delta";
      model.Pbp.Bkg.Ctau.Prompt=CtauModel::Delta;
    }
  }
  if (isPbp && incJpsi && incNonPrompt) {
    if (parIni.count("Model_JpsiNoPR_Pbp")>0) {
      model.Pbp.Jpsi.Ctau.NonPrompt = CtauModelDictionary[parIni["Model_JpsiNoPR_Pbp"]];
      if (model.Pbp.Jpsi.Ctau.NonPrompt==CtauModel(0)) {
        cout << "[ERROR] The Jpsi non-prompt ctau model: " << parIni["Model_JpsiNoPR_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Jpsi non-prompt ctau model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  }
  if (isPbp && incJpsi && incPrompt) {
    if (parIni.count("Model_JpsiPR_Pbp")>0) {
      model.Pbp.Jpsi.Ctau.Prompt = CtauModelDictionary[parIni["Model_JpsiPR_Pbp"]];
      if (model.Pbp.Jpsi.Ctau.Prompt==CtauModel(0)) {
        cout << "[ERROR] The Jpsi prompt ctau model: " << parIni["Model_JpsiPR_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      parIni["Model_JpsiPR_Pbp"] = "Delta";
      model.Pbp.Jpsi.Ctau.Prompt=CtauModel::Delta;
    }
  }
  if (isPbp && incPsi2S && incNonPrompt) {
    if (parIni.count("Model_Psi2SNoPR_Pbp")>0) {
      model.Pbp.Psi2S.Ctau.NonPrompt = CtauModelDictionary[parIni["Model_Psi2SNoPR_Pbp"]];
      if (model.Pbp.Psi2S.Ctau.NonPrompt==CtauModel(0)) {
        cout << "[ERROR] The psi(2S) non-prompt ctau model: " << parIni["Model_Psi2SNoPR_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] psi(2S) non-prompt ctau model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  }
  if (isPbp && incPsi2S && incPrompt) {
    if (parIni.count("Model_Psi2SPR_Pbp")>0) {
      model.Pbp.Psi2S.Ctau.Prompt = CtauModelDictionary[parIni["Model_Psi2SPR_Pbp"]];
      if (model.Pbp.Psi2S.Ctau.Prompt==CtauModel(0)) {
        cout << "[ERROR] The psi(2S) prompt ctau model: " << parIni["Model_Psi2SPR_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      parIni["Model_Psi2SPR_Pbp"] = "Delta";
      model.Pbp.Psi2S.Ctau.Prompt=CtauModel::Delta;
    }
  }

  if (!isPbp && incBkg && incNonPrompt) {
    if (parIni.count("Model_BkgNoPR_PP")>0) {
      model.PP.Bkg.Ctau.NonPrompt = CtauModelDictionary[parIni["Model_BkgNoPR_PP"]];
      if (model.PP.Bkg.Ctau.NonPrompt==CtauModel(0)) {
        cout << "[ERROR] The background non-prompt ctau model: " << parIni["Model_BkgNoPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Background non-prompt ctau model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbp && incBkg && incPrompt) {
    if (parIni.count("Model_BkgPR_PP")>0) {
      model.PP.Bkg.Ctau.Prompt = CtauModelDictionary[parIni["Model_BkgPR_PP"]];
      if (model.PP.Bkg.Ctau.Prompt==CtauModel(0)) {
        cout << "[ERROR] The background prompt ctau model: " << parIni["Model_BkgPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      parIni["Model_BkgPR_PP"] = "Delta";
      model.PP.Bkg.Ctau.Prompt=CtauModel::Delta;
    }
  }
  if (!isPbp && incJpsi && incNonPrompt) {
    if (parIni.count("Model_JpsiNoPR_PP")>0) {
      model.PP.Jpsi.Ctau.NonPrompt = CtauModelDictionary[parIni["Model_JpsiNoPR_PP"]];
      if (model.PP.Jpsi.Ctau.NonPrompt==CtauModel(0)) {
        cout << "[ERROR] The Jpsi non-prompt ctau model: " << parIni["Model_JpsiNoPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Jpsi non-prompt ctau model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbp && incJpsi && incPrompt) {
    if (parIni.count("Model_JpsiPR_PP")>0) {
      model.PP.Jpsi.Ctau.Prompt = CtauModelDictionary[parIni["Model_JpsiPR_PP"]];
      if (model.PP.Jpsi.Ctau.Prompt==CtauModel(0)) {
        cout << "[ERROR] The Jpsi prompt ctau model: " << parIni["Model_JpsiPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      parIni["Model_JpsiPR_PP"] = "Delta";
      model.PP.Jpsi.Ctau.Prompt=CtauModel::Delta;
    }
  }
  if (!isPbp && incPsi2S && incNonPrompt) {
    if (parIni.count("Model_Psi2SNoPR_PP")>0) {
      model.PP.Psi2S.Ctau.NonPrompt = CtauModelDictionary[parIni["Model_Psi2SNoPR_PP"]];
      if (model.PP.Psi2S.Ctau.NonPrompt==CtauModel(0)) {
        cout << "[ERROR] The psi(2S) non-prompt ctau model: " << parIni["Model_Psi2SNoPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] psi(2S) non-prompt ctau model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbp && incPsi2S && incPrompt) {
    if (parIni.count("Model_Psi2SPR_PP")>0) {
      model.PP.Psi2S.Ctau.Prompt = CtauModelDictionary[parIni["Model_Psi2SPR_PP"]];
      if (model.PP.Psi2S.Ctau.Prompt==CtauModel(0)) {
        cout << "[ERROR] The psi(2S) prompt ctau model: " << parIni["Model_Psi2SPR_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      parIni["Model_Psi2SPR_PP"] = "Delta";
      model.PP.Psi2S.Ctau.Prompt=CtauModel::Delta;
    }
  }

  return true;
};


void setCtauGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label, double binWidth, bool fitCtauRes)
{
  Double_t ctauMax; Double_t ctauMin;
  myws.data(Form("dOS_%s", label.c_str()))->getRange(*myws.var("ctau"), ctauMin, ctauMax);
  ctauMin -= 0.00001;  ctauMax += 0.00001;
  int nBins = min(int( round((ctauMax - ctauMin)/binWidth) ), 1000);
  if (fitCtauRes) {
    TH1D* hTot = (TH1D*)myws.data(Form("dOS_%s", label.c_str()))->createHistogram("TMP", *myws.var("ctau"), Binning(nBins, ctauMin, ctauMax));
    vector<double> rangeCtau; 
    getCtauErrRange(hTot, (int)(ceil(2)), rangeCtau);
    hTot->Delete();
    ctauMin = rangeCtau[0];
    ctauMax = rangeCtau[1];
    if (abs(ctauMax)>abs(ctauMin)) { ctauMax = abs(ctauMin); } else { ctauMin = -1.0*abs(ctauMax); }
    if (ctauMin<cut.dMuon.ctau.Min) { ctauMin = cut.dMuon.ctau.Min; }
    if (ctauMax>cut.dMuon.ctau.Max) { ctauMax = cut.dMuon.ctau.Max; }
    if (parIni.count("ctauResCut")>0 && parIni["ctauResCut"]!="") {
      parIni["ctauResCut"].erase(parIni["ctauResCut"].find("["), string("[").length());
      parIni["ctauResCut"].erase(parIni["ctauResCut"].find("]"), string("]").length());
      double ctauCut = 0.0;
      try {
        ctauCut = std::stod(parIni["ctauResCut"].c_str());
      } catch (const std::invalid_argument&) {
        cout << "[WARNING] ctauResCut is not a number, will ignore it!" << endl;
      }
      ctauCut = abs(ctauCut); if(ctauCut>0.0) { ctauMax = abs(ctauCut); ctauMin = -1.0*abs(ctauCut); }
    }
  }
  cout << "[INFO] Range from data: ctauMin: " << ctauMin << "  ctauMax: " << ctauMax << endl;
  myws.var("ctau")->setRange("CtauWindow", ctauMin, ctauMax);
  parIni["CtauRange_Cut"]   = Form("(%.12f <= ctau && ctau < %.12f)", ctauMin, ctauMax);
  cut.dMuon.ctau.Max = ctauMax;
  cut.dMuon.ctau.Min = ctauMin;
  myws.var("ctau")->setRange("CtauFullWindow", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max);
  myws.var("ctau")->setRange("FullWindow", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max);
  myws.var("ctau")->setRange("SideBandTOP_FULL", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max); 
  myws.var("ctau")->setRange("SideBandMID_FULL", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max);
  myws.var("ctau")->setRange("SideBandBOT_FULL", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max); 
  myws.var("ctau")->setRange("SideBandMID_JPSI", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max);
  myws.var("ctau")->setRange("SideBandBOT_JPSI", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max);
  myws.var("ctau")->setRange("SideBandTOP_PSI2S", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max); 
  myws.var("ctau")->setRange("SideBandMID_PSI2S", cut.dMuon.ctau.Min, cut.dMuon.ctau.Max);

  return;
};


void setCtauFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbp, bool fitSideBand)
{
  if (TAG.find("_")!=std::string::npos) TAG.erase(TAG.find("_"));
  FileName = Form("%sctau%s/%s/result/FIT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), (fitSideBand?"SB":""), TAG.c_str(), "CTAU", TAG.c_str(), (isPbp?"Pbp":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);

  return;
};
 

void setCtauCutParameters(struct KinCuts& cut, bool incNonPrompt)
{
  // Define the ctau range
  if (cut.dMuon.ctau.Min==-1000. && cut.dMuon.ctau.Max==1000.) { 
    // Default ctau values, means that the user did not specify a ctau range
    if (incNonPrompt) {
      cut.dMuon.ctau.Min = -50.0;
      cut.dMuon.ctau.Max = 100.0;
    } else {
      cut.dMuon.ctau.Min = -2.0;
      cut.dMuon.ctau.Max = 2.0;
    }
  }
        
  cout << "[INFO] Setting ctau range to min: " << cut.dMuon.ctau.Min << " and max " << cut.dMuon.ctau.Max << endl;
  
  return;
};


#endif // #ifndef fitCharmoniaCtauModel_C
