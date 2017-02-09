#ifndef fitCharmoniaCtauMassModel_C
#define fitCharmoniaCtauMassModel_C

#include "Utilities/initClasses.h"
#include "fitCharmoniaMassModel.C"
#include "fitCharmoniaCtauModel.C"
#include "fitCharmoniaCtauErrModel.C"
#include "fitCharmoniaCtauTrueModel.C"
#include "drawMassFrom2DPlot.C"
#include "drawCtauFrom2DPlot.C"
#include "drawCtauMass2DPlot.C"

void setCtauMassFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbp);


bool fitCharmoniaCtauMassModel( RooWorkspace& myws,             // Local Workspace
                                RooWorkspace& inputWorkspace,   // Workspace with all the input RooDatasets
                                struct KinCuts& cut,            // Variable containing all kinematic cuts
                                map<string, string>&  parIni,   // Variable containing all initial parameters
                                struct InputOpt& opt,           // Variable with run information (kept for legacy purpose)
                                string outputDir,               // Path to output directory
                                // Select the type of datasets to fit
                                string DSTAG,                   // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
                                bool isPbp        = false,     // isPbp = false for pp, true for Pbp
                                // Select the type of object to fit
                                bool incJpsi       = true,      // Includes Jpsi model
                                bool incPsi2S      = false,     // Includes Psi(2S) model
                                // Select the fitting options
                                map<string, string> inputFitDir={},// User-defined Location of the fit results
                                int  numCores      = 2,         // Number of cores used for fitting
                                // Select the drawing options
                                bool setLogScale   = true,      // Draw plot with log scale
                                bool incSS         = false,     // Include Same Sign data
                                map<string, double> binWidth={} // User-defined Location of the fit results
                                )  
{
  
  if (DSTAG.find("_")!=std::string::npos) DSTAG.erase(DSTAG.find("_"));

  // Check if input dataset is MC
  if (DSTAG.find("MC")!=std::string::npos) {
    cout << "[ERROR] We can not make 2D fits on MC, change your input settings!" << endl; return false;
  }
  // Check if user has selected an object
  if (!incJpsi && !incPsi2S) {
    cout << "[ERROR] We need to include signal to perform 2D fits, change your input settings!" << endl; return false;
  }

  string pdfType = "pdfCTAUMASS";
  string COLL = (isPbp ? "Pbp" : "PP" );
  string pdfName = Form("%s_Tot_%s", pdfType.c_str(), COLL.c_str());

  bool isMC          = false;
  bool incBkg        = true;      // Includes Background model
  bool incPrompt     = true;      // Includes Prompt ctau model
  bool incNonPrompt  = true;     // Includes NonPrompt ctau model
  bool fitSideBand   = false;
  bool usePerEventError = true;
  bool wantPureSMC = false;
  bool fitMass = true;
  bool fitCtau = true;

  // Set all the cuts on data
  setMassCutParameters(cut, incJpsi, incPsi2S, isMC, true);
  if (usePerEventError) {
    // check if we have already done the ctauErr fits. If yes, load their parameters
    string FileName = "";
    string pdfName = Form("pdfCTAUERR_Tot_%s", COLL.c_str());
    string plotLabel = "";
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
  setCtauCutParameters(cut, incNonPrompt);

  // Import the local datasets
  string label = ((DSTAG.find(COLL.c_str())!=std::string::npos) ? DSTAG.c_str() : Form("%s_%s", DSTAG.c_str(), COLL.c_str()));
  string dsName = Form("dOS_%s", label.c_str());
  if ( !(myws.data(dsName.c_str())) ) {
    int importID = importDataset(myws, inputWorkspace, cut, label, fitSideBand);
    if (importID<0) { return false; }
    else if (importID==0) { cout << "[ERROR] Local dataset failed to import!" << endl; return false; }
  }
  double numEntries = myws.data(dsName.c_str())->sumEntries(); 
  if ((numEntries<=0) || !(myws.data(dsName.c_str()))) { cout << "[ERROR] No local dataset was found to perform the fit!" << endl; return false; }
 
  // Set global parameters
  setMassGlobalParameterRange(myws, parIni, cut, incJpsi, incPsi2S, incBkg, false);
  setCtauErrGlobalParameterRange(myws, parIni, cut, "", binWidth["CTAUERR"], true);
  setCtauGlobalParameterRange(myws, parIni, cut, label, binWidth["CTAU"], false);

  // Cut the RooDataSet
  RooDataSet* dataToFit = (RooDataSet*)(myws.data(dsName.c_str())->reduce(parIni["CtauRange_Cut"].c_str()))->Clone((dsName+"_CTAUCUT").c_str());
  myws.import(*dataToFit); string dsNameCut = dsName+"_CTAUCUT";

  //// LOAD MASS PDF
  if (fitMass) {
    // Setting extra input information needed by each fitter
    string iMassFitDir = inputFitDir["MASS"];
    double ibWidth = binWidth["MASS"];
    bool loadMassFitResult = true;
    bool doMassFit = true;
    bool importDS = false;
    bool getMeanPT = false;
    bool zoomPsi = false;
    const char* applyCorr = "";
    bool doSimulFit = false;
    bool cutCtau = false;

    if ( !fitCharmoniaMassModel( myws, inputWorkspace, cut, parIni, opt, outputDir, 
                                 DSTAG, isPbp, importDS,
                                 incJpsi, incPsi2S, true, 
                                 doMassFit, cutCtau, doSimulFit, wantPureSMC, applyCorr, loadMassFitResult, iMassFitDir, numCores, 
                                 setLogScale, incSS, zoomPsi, ibWidth, getMeanPT
                                 ) 
         ) { return false; }
    if (myws.pdf(Form("pdfMASS_Tot_%s", (isPbp?"Pbp":"PP")))) {
      cout << "[INFO] Setting mass parameters to constant!" << endl;
      myws.pdf(Form("pdfMASS_Tot_%s", (isPbp?"Pbp":"PP")))->getParameters(RooArgSet(*myws.var("invMass")))->setAttribAll("Constant", kTRUE); 
    } else { cout << "[ERROR] Mass PDF was not found!" << endl; return false; }
    //if (myws.var(Form("N_Bkg_%s", (isPbp?"Pbp":"PP")))) { myws.var(Form("N_Bkg_%s", (isPbp?"Pbp":"PP")))->setConstant(kFALSE); }
    //if (myws.var(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")))) { myws.var(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")))->setConstant(kFALSE); }
    //if (myws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))) { myws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))->setConstant(kFALSE); }
  }

  //// LOAD CTAU ERROR PDF
  if (usePerEventError) {
    // Setting extra input information needed by each fitter
    bool loadCtauErrFitResult = true;
    bool doCtauErrFit = false;
    bool importDS = false;
    bool wantPureSMC = false;
    
    if ( !fitCharmoniaCtauErrModel( myws, inputWorkspace, cut, parIni, opt, outputDir, 
                                    DSTAG, isPbp, importDS, 
                                    incJpsi, incPsi2S, incBkg, 
                                    doCtauErrFit, wantPureSMC, loadCtauErrFitResult, inputFitDir, numCores, 
                                    setLogScale, incSS, binWidth
                                    ) 
         ) { return false; }
  }

  // Set models based on initial parameters
  struct OniaModel model;
  if (!setCtauModel(model, parIni, isPbp, incJpsi, incPsi2S, incBkg, incPrompt, incNonPrompt)) { return false; }

  // Build the Fit Model
  if (!buildCharmoniaCtauModel(myws, (isPbp ? model.Pbp : model.PP), parIni, dsName, isPbp, incBkg, incJpsi, incPsi2S, incPrompt, incNonPrompt,  numEntries))  { return false; }

  //// LOAD CTAU SIDEBAND PDF
  if (fitCtau) {
    // check if we have already done the sideband fits. If yes, load their results
    string FileName = "";
    bool fitSB = true;
    string plotLabel = Form("_BkgNoPR_%s_CtauRes_%s", parIni[Form("Model_BkgNoPR_%s", COLL.c_str())].c_str(), parIni[Form("Model_CtauRes_%s", COLL.c_str())].c_str());
    string DSTAG = Form("DATA_%s", (isPbp?"Pbp":"PP"));
    setCtauFileName(FileName, (inputFitDir["CTAUSB"]=="" ? outputDir : inputFitDir["CTAUSB"]), DSTAG, plotLabel, cut, isPbp, fitSB);
    bool found = false;
    if (!found && gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAUSB"]!="") {
      setCtauFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbp, fitSB);
    } else if (inputFitDir["CTAUSB"]!="") { found = true; }
    if (!found && gSystem->AccessPathName(FileName.c_str())) {
      cout << "[ERROR] User Input File : " << FileName << " was not found!" << endl;
      return false;
    }
    if ( !loadPreviousFitResult(myws, FileName, DSTAG, isPbp) ) {
      cout << "[ERROR] The ctau sideband fit results were not loaded!" << endl;
      return false;
    } else { 
      cout << "[INFO] The ctau sideband fits were found, so I'll load the fit results." << endl; 
    }
    if (myws.pdf(Form("pdfCTAU_BkgPR_%s", (isPbp?"Pbp":"PP")))) {
      cout << "[INFO] Setting Prompt Background parameters to constant!" << endl;
      myws.pdf(Form("pdfCTAU_BkgPR_%s", (isPbp?"Pbp":"PP")))->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE); 
    } else { cout << "[ERROR] Prompt Background PDF was not found!" << endl; return false; }
    if (myws.pdf(Form("pdfCTAU_BkgNoPR_%s", (isPbp?"Pbp":"PP")))) {
      cout << "[INFO] Setting NonPrompt Background parameters to constant!" << endl;
      myws.pdf(Form("pdfCTAU_BkgNoPR_%s", (isPbp?"Pbp":"PP")))->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE); 
    } else { cout << "[ERROR] NonPrompt Background PDF was not found!" << endl; return false; }
    if (!setConstant(myws, Form("s1_CtauRes_%s", COLL.c_str()), false)) { return false; }
  }

  //// LOAD CTAU TRUE PDF RESULTS
  if (fitCtau) {
    // check if we have already done the ctau true fits. If yes, load their results
    string FileName = "";
    //string ModelName = Form("Model_JpsiNoPR_%s", COLL.c_str());
    string ModelName = Form("Model_Psi2SNoPR_%s", COLL.c_str());
    string plotLabel = Form("_CtauTrue_%s", parIni[ModelName.c_str()].c_str());
    //string DSTAG = Form("MCJPSINOPR_%s", (isPbp?"Pbp":"PP"));
    string DSTAG = Form("MCPSI2SNOPR_%s", (isPbp?"Pbp":"PP"));
    setCtauTrueFileName(FileName, (inputFitDir["CTAUTRUE"]=="" ? outputDir : inputFitDir["CTAUTRUE"]), DSTAG, plotLabel, cut, isPbp);
    bool found = false;
    if (!found && gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAUTRUE"]!="") {
      plotLabel = string(Form("_CtauTrue_%s_NoBkg", parIni[ModelName.c_str()].c_str()));
      setCtauTrueFileName(FileName, (inputFitDir["CTAUTRUE"]=="" ? outputDir : inputFitDir["CTAUTRUE"]), DSTAG, plotLabel, cut, isPbp);
    } else if (inputFitDir["CTAUTRUE"]!="") { found = true; }
    if (!found && gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAUTRUE"]!="") {
      plotLabel = string(Form("_CtauTrue_%s", parIni[ModelName.c_str()].c_str()));
      setCtauTrueFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbp);
    } else if (inputFitDir["CTAUTRUE"]!="") { found = true; }
    if (!found && gSystem->AccessPathName(FileName.c_str())) {
      plotLabel = string(Form("_CtauTrue_%s_NoBkg", parIni[ModelName.c_str()].c_str()));
      setCtauTrueFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbp);
    } else { found = true; }
    if (!found && gSystem->AccessPathName(FileName.c_str())) {
      cout << "[ERROR] User Input File : " << FileName << " was not found!" << endl;
      return false;
    }
    if ( !loadPreviousFitResult(myws, FileName, DSTAG, isPbp) ) {
      cout << "[ERROR] The ctau true fit results were not loaded!" << endl;
      return false;
    } else { 
      cout << "[INFO] The ctau true fits were found, so I'll load the fit results." << endl; 
    }
  }

  // save the initial values of the model we've just created
  RooArgSet* params = (RooArgSet*) myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("invMass"), *myws.var("ctauErr")));
  myws.saveSnapshot((pdfName+"_parIni").c_str(),*params,kTRUE);
  delete params;

  string plotLabel = "";
  map<string, bool> plotLabels = {{"Jpsi",      (incJpsi)},
                                  {"Psi2S",     (incPsi2S)},
                                  {"Bkg",       (incBkg)},
                                  {"JpsiNoPR",  (incJpsi&&incNonPrompt)}, 
                                  {"Psi2SNoPR", (incPsi2S&&incNonPrompt)}, 
                                  {"BkgNoPR",   (incBkg&&incNonPrompt)}, 
                                  {"CtauRes",   (true)}};
  for (map<string, bool>::iterator iter = plotLabels.begin(); iter!=plotLabels.end(); iter++) {
    string obj = iter->first;
    bool cond = iter->second;
    if (cond && parIni.count(Form("Model_%s_%s", obj.c_str(), COLL.c_str()))>0) { 
      plotLabel = plotLabel + Form("_%s_%s", obj.c_str(), parIni[Form("Model_%s_%s", obj.c_str(), COLL.c_str())].c_str()); 
    }
  }

  // check if we have already done this fit. If yes, do nothing and return true.
  string FileName = "";
  setCtauMassFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbp);
  bool found =  true;
  RooArgSet *newpars = (RooArgSet*) myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr"), *myws.var("invMass")));
  found = found && isFitAlreadyFound(newpars, FileName, pdfName.c_str());
  if (found) {
    cout << "[INFO] This ctau fit was already done, so I'll just go to the next one." << endl;
    return true;
  }

  // Fit the Datasets
  bool isWeighted = myws.data(dsName.c_str())->isWeighted();
  RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsNameCut.c_str()), Extended(kTRUE), NumCPU(numCores), ConditionalObservables(*myws.var("ctauErr")), SumW2Error(isWeighted), Save());
  fitResult->Print("v");
  myws.import(*fitResult, Form("fitResult_%s", pdfName.c_str()));
  // Draw the mass plot
  drawCtauFrom2DPlot(myws, outputDir, opt, cut, parIni, plotLabel, DSTAG, isPbp, incJpsi, incPsi2S, incBkg, setLogScale, incSS, binWidth["CTAU"]);
  drawMassFrom2DPlot(myws, outputDir, opt, cut, parIni, plotLabel, DSTAG, isPbp, incJpsi, incPsi2S, incBkg, setLogScale, incSS, binWidth["MASS"]);
  drawCtauMass2DPlot(myws, outputDir, cut, plotLabel, DSTAG, isPbp, binWidth);
  // Save the results
  myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE);
  saveWorkSpace(myws, Form("%sctauMass/%s/result", outputDir.c_str(), DSTAG.c_str()), FileName);

  return true;
};


void setCtauMassFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbp)
{
  if (TAG.find("_")!=std::string::npos) TAG.erase(TAG.find("_"));
  FileName = Form("%sctauMass/%s/result/FIT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), "CTAUMASS", TAG.c_str(), (isPbp?"Pbp":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);

  return;
};


#endif // #ifndef fitCharmoniaCtauMassModel_C
