#ifndef fitCharmoniaCtauErrModel_C
#define fitCharmoniaCtauErrModel_C

#include "Utilities/initClasses.h"
#include "buildCharmoniaCtauErrModel.C"
#include "fitCharmoniaMassModel.C"
#include "drawCtauErrorPlot.C"

void setCtauErrCutParameters(struct KinCuts& cut);
bool isCtauErrPdfAlreadyFound(RooWorkspace& myws, string FileName, vector<string> pdfNames, bool loadCtauErrPdf=false);
void setCtauErrFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbp, bool cutSideBand);
void setCtauErrGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label, double binWidth, bool useForCtauFits=false);
void reNormMassVar( RooWorkspace& myws, string obj, bool isPbp);


bool fitCharmoniaCtauErrModel( RooWorkspace& myws,             // Local Workspace
                               RooWorkspace& inputWorkspace,   // Workspace with all the input RooDatasets
                               struct KinCuts& cut,            // Variable containing all kinematic cuts
                               map<string, string>&  parIni,   // Variable containing all initial parameters
                               struct InputOpt& opt,           // Variable with run information (kept for legacy purpose)
                               string outputDir,               // Path to output directory
                               // Select the type of datasets to fit
                               string DSTAG,                   // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
                               bool isPbp         = false,    // isPbp = false for pp, true for Pbp
                               bool importDS       = true,     // Select if the dataset is imported in the local workspace
                               // Select the type of object to fit
                               bool incJpsi        = false,     // Includes Jpsi model
                               bool incPsi2S       = true,     // Includes Psi(2S) model
                               bool incBkg         = true,     // Includes Background model
                               // Select the fitting options
                               bool doCtauErrPdf   = true,     // Flag to indicate if we want to make the ctau Error Pdf
                               bool wantPureSMC    = false,    // Flag to indicate if we want to use pure signal MC
                               bool loadCtauErrPdf = false,    // Load previous ctau Error Pdf
                               map<string, string> inputFitDir={},// User-defined Location of the fit results
                               int  numCores       = 2,         // Number of cores used for fitting
                               // Select the drawing options
                               bool setLogScale    = true,      // Draw plot with log scale
                               bool incSS          = false,     // Include Same Sign data
                               map<string, double> binWidth={} // User-defined Location of the fit results
                               )  
{

  if (DSTAG.find("_")!=std::string::npos) DSTAG.erase(DSTAG.find("_"));

  // Check if input dataset is MC
  bool isMC = false;
  if (DSTAG.find("MC")!=std::string::npos) {
    if (incJpsi && incPsi2S) { 
      cout << "[ERROR] We can only fit one type of signal using MC" << endl; return false; 
    }
    isMC = true;
  }
  wantPureSMC = (isMC && wantPureSMC);
  bool cutSideBand = (incBkg && (!incPsi2S && !incJpsi));

  setMassCutParameters(cut, incJpsi, incPsi2S, isMC);
  setCtauErrCutParameters(cut);
  
  string pdfType = "pdfCTAUERR";
  string COLL = (isPbp ? "Pbp" : "PP" );
  incJpsi = false;	
  incPsi2S = true;   

  // Define pdf and plot names
  vector<string> pdfNames;
  string plotLabel = "";
  pdfNames.push_back(Form("%sTot_Tot_%s", pdfType.c_str(), COLL.c_str()));
  if (incJpsi)  { plotLabel = plotLabel + "_Jpsi";   pdfNames.push_back(Form("%s_Jpsi_%s", pdfType.c_str(), COLL.c_str()));  }
  if (incPsi2S) { plotLabel = plotLabel + "_Psi2S";  pdfNames.push_back(Form("%s_Psi2S_%s", pdfType.c_str(), COLL.c_str())); }
  if (!isMC)    { plotLabel = plotLabel + "_Bkg";    pdfNames.push_back(Form("%s_Bkg_%s", pdfType.c_str(), COLL.c_str()));   }
  if (wantPureSMC) { plotLabel = plotLabel + "_NoBkg"; }
  incPsi2S = (incPsi2S ? true : false);  

  // check if we have already done this fit. If yes, do nothing and return true.
  string FileName = "";
  setCtauErrFileName(FileName, (inputFitDir["CTAUERR"]=="" ? outputDir : inputFitDir["CTAUERR"]), DSTAG, plotLabel, cut, isPbp, cutSideBand);
  if (gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAUERR"]!="") {
    cout << "[WARNING] User Input File : " << FileName << " was not found!" << endl;
    if (loadCtauErrPdf) return false;
    setCtauErrFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbp, cutSideBand);
   }
  bool found =  true; bool skipCtauErrPdf = !doCtauErrPdf;
  found = found && isCtauErrPdfAlreadyFound(myws, FileName, pdfNames, loadCtauErrPdf);
  if (found) {
    if (loadCtauErrPdf) {
      cout << "[INFO] This ctauErr Pdf was already made, so I'll load the pdf." << endl;
    } else {
      cout << "[INFO] This ctauErr Pdf was already made, so I'll just go to the next one." << endl;
    }
    return true;
  }
 
  
  // Import the local datasets
  double numEntries = 1000000;
  string label = ((DSTAG.find(COLL.c_str())!=std::string::npos) ? DSTAG.c_str() : Form("%s_%s", DSTAG.c_str(), COLL.c_str()));
  if (wantPureSMC) label = Form("%s_NoBkg", label.c_str());
  string dsName = Form("dOS_%s", label.c_str());
  if (importDS) {
    if ( !(myws.data(dsName.c_str())) ) {
      int importID = importDataset(myws, inputWorkspace, cut, label, cutSideBand);
      if (importID<0) { return false; }
      else if (importID==0) { doCtauErrPdf = false; }
    }
    numEntries = myws.data(dsName.c_str())->sumEntries(); if (numEntries<=0) { doCtauErrPdf = false; }
  }
  else if (doCtauErrPdf && !(myws.data(dsName.c_str()))) { cout << "[ERROR] No local dataset was found to make the ctau Error Pdf!" << endl; return false; }

  // Set global parameters
  setCtauErrGlobalParameterRange(myws, parIni, cut, label, binWidth["CTAUERRFORCUT"]);

  if (!isMC && (incJpsi || incPsi2S)) { 
    // Setting extra input information needed by each fitter
    string iMassFitDir = inputFitDir["MASS"];
    double ibWidth = binWidth["MASS"];
    bool loadMassFitResult = true;
    bool doMassFit = false;
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
    if (myws.pdf(Form("pdfMASS_Bkg_%s", (isPbp?"Pbp":"PP"))))   { reNormMassVar(myws, "Bkg", isPbp);   }
    if (myws.pdf(Form("pdfMASS_Jpsi_%s", (isPbp?"Pbp":"PP"))))  { reNormMassVar(myws, "Jpsi", isPbp);  }
    if (myws.pdf(Form("pdfMASS_Psi2S_%s", (isPbp?"Pbp":"PP")))) { reNormMassVar(myws, "Psi2S", isPbp); }
  }

  if (skipCtauErrPdf==false) {
    // Create the ctau Error Pdf
    // Build the Ctau Error Template
    if (!buildCharmoniaCtauErrModel(myws, parIni, cut, dsName, incJpsi, incPsi2S, binWidth["CTAUERR"], numEntries))  { return false; }
    string pdfName = Form("%s_Tot_%s", pdfType.c_str(), COLL.c_str());
    bool isWeighted = myws.data(dsName.c_str())->isWeighted();
    //RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), SumW2Error(isWeighted), Range("CtauErrWindow"), NumCPU(numCores), Save());
    //fitResult->Print("v");
    //myws.import(*fitResult, Form("fitResult_%s", pdfName.c_str()));
  
    // Draw the mass plot
    drawCtauErrorPlot(myws, outputDir, opt, cut, parIni, plotLabel, DSTAG, isPbp, incJpsi, incPsi2S, incBkg, wantPureSMC, setLogScale, incSS, binWidth["CTAUERR"]);
   
    // Save the results
    string FileName = ""; setCtauErrFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbp, cutSideBand);
    RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(*(myws.var("ctauErr")));
    myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE);
    if (!saveWorkSpace(myws, Form("%sctauErr%s/%s/result", outputDir.c_str(), (cutSideBand?"SB":""), DSTAG.c_str()), FileName)) { return false; };
  }
  
  return true;
};


void setCtauErrGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label, double binWidth, bool useForCtauFits)
{
  if (!useForCtauFits) {
    Double_t ctauErrMax; Double_t ctauErrMin;
    myws.data(Form("dOS_%s", label.c_str()))->getRange(*myws.var("ctauErr"), ctauErrMin, ctauErrMax);
    ctauErrMin -= 0.0001;  ctauErrMax += 0.0001;
    int nBins = min(int( round((ctauErrMax - ctauErrMin)/binWidth) ), 1000);
    TH1D* hTot = (TH1D*)myws.data(Form("dOS_%s", label.c_str()))->createHistogram("TMP", *myws.var("ctauErr"), Binning(nBins, ctauErrMin, ctauErrMax));
    vector<double> rangeErr; getCtauErrRange(hTot, (int)(ceil(2)), rangeErr); // KEEP THIS NUMBER WITH 2, JUST KEEP IT LIKE THAT :D
    hTot->Delete();
    ctauErrMin = rangeErr[0];
    ctauErrMax = rangeErr[1];
    if (ctauErrMin<cut.dMuon.ctauErr.Min) { ctauErrMin = cut.dMuon.ctauErr.Min; }
    if (ctauErrMax>cut.dMuon.ctauErr.Max) { ctauErrMax = cut.dMuon.ctauErr.Max; }
    cut.dMuon.ctauErr.Max = ctauErrMax;
    cut.dMuon.ctauErr.Min = ctauErrMin;
  }
  cout << "[INFO] Range from data: ctauErrMin: " << cut.dMuon.ctauErr.Min << "  ctauErrMax: " << cut.dMuon.ctauErr.Max << endl;
  myws.var("ctauErr")->setRange("CtauErrWindow", cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max);
  parIni["CtauErrRange_Cut"]   = Form("(%.12f <= ctauErr && ctauErr < %.12f)", cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max);
  myws.var("ctauErr")->setRange("CtauErrFullWindow", cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max);
  myws.var("ctauErr")->setRange("FullWindow", cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max);
  myws.var("ctauErr")->setRange("SideBandTOP_FULL", cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max); 
  myws.var("ctauErr")->setRange("SideBandMID_FULL", cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max);
  myws.var("ctauErr")->setRange("SideBandBOT_FULL", cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max); 
  myws.var("ctauErr")->setRange("SideBandMID_JPSI", cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max);
  myws.var("ctauErr")->setRange("SideBandBOT_JPSI", cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max);
  myws.var("ctauErr")->setRange("SideBandTOP_PSI2S", cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max); 
  myws.var("ctauErr")->setRange("SideBandMID_PSI2S", cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max);

  return;
};


void setCtauErrFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbp, bool cutSideBand)
{
  if (TAG.find("_")!=std::string::npos) TAG.erase(TAG.find("_"));
  FileName = Form("%sctauErr%s/%s/result/FIT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), (cutSideBand?"SB":""), TAG.c_str(), "CTAUERR", TAG.c_str(), (isPbp?"Pbp":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);

  return;
};


bool isCtauErrPdfAlreadyFound(RooWorkspace& myws, string FileName, vector<string> pdfNames, bool loadCtauErrPdf)
{
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[INFO] Results not found for: " << FileName << endl;
    return false; // File was not found
  }
  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;

  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    cout << "[INFO] Workspace not found in: " << FileName << endl;
    file->Close(); delete file;
    return false;
  }

  bool found = true;
  for (unsigned int i=0; i<pdfNames.size(); i++) {
    string pdfName = pdfNames.at(i);
    string dataName = pdfName;
    dataName.replace(dataName.find("pdf"), string("pdf").length(), "dh");
    if ( !(ws->pdf(pdfName.c_str())) || !(ws->data(dataName.c_str())) ) {
      cout << "[INFO] " << pdfName << " was not found in: " << FileName << endl; found = false;
    }
    if (loadCtauErrPdf && found) {
      myws.import(*(ws->pdf(pdfName.c_str())));
      myws.import(*(ws->data(dataName.c_str())));
      if (ws->pdf(pdfName.c_str()))   { cout << "[INFO] Pdf " << pdfName << " succesfully imported!" << endl;       }
      else {  cout << "[ERROR] Pdf " << pdfName << " import failed!" << endl; found = false; }
      if (ws->data(dataName.c_str())) { cout << "[INFO] DataHist " << dataName << " succesfully imported!" << endl; }
      else {  cout << "[ERROR] DataHist " << dataName << " import failed!" << endl; found = false; }
    }
  }
  
  delete ws;
  file->Close(); delete file;

  return found;
};


void setCtauErrCutParameters(struct KinCuts& cut)
{
  // Define the ctau error range
  if (cut.dMuon.ctauErr.Min==-1000.0 && cut.dMuon.ctauErr.Max==1000.0) { 
    // Default ctau error values, means that the user did not specify a ctau error range
    cut.dMuon.ctauErr.Min = 0.0; 
    cut.dMuon.ctauErr.Max = 10.0;
  }
  cout << "[INFO] Setting ctauErr range to min: " << cut.dMuon.ctauErr.Min << " and max " << cut.dMuon.ctauErr.Max << endl;

  return;
};


void reNormMassVar( RooWorkspace& myws, string obj, bool isPbp)
{
  string pdfName = Form("pdfMASS_%s_%s", obj.c_str(), (isPbp?"Pbp":"PP"));
  string varName = Form("N_%s_%s", obj.c_str(), (isPbp?"Pbp":"PP"));
  if (obj=="Bkg") { ((RooChebychev*)myws.pdf(pdfName.c_str()))->selectNormalizationRange("InclusiveMassRegion", kTRUE); }
  double NormFactor = 1.0/(myws.pdf(pdfName.c_str())->createIntegral(*myws.var("invMass"), NormSet(*myws.var("invMass")), Range("InclusiveMassRegion"))->getValV());
  myws.var(varName.c_str())->setVal(myws.var(varName.c_str())->getValV());
  if (myws.var(varName.c_str())) { myws.var(varName.c_str())->setConstant(kFALSE); }
  //cout << "[INFO] NORM FACTOR FOR " << obj << " : " << NormFactor << endl;
  return;
};

#endif // #ifndef fitCharmoniaCtauErrModel_C
