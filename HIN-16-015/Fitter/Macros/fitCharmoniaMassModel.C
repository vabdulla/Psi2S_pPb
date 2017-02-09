#ifndef fitCharmoniaMassModel_C
#define fitCharmoniaMassModel_C

#include "Utilities/initClasses.h"
#include "buildCharmoniaMassModel.C"
#include "drawMassPlot.C"

void setCtauCuts(struct KinCuts& cut, bool isPbp);
bool setMassModel( struct OniaModel& model, map<string, string>&  parIni, bool isPbp, bool incJpsi, bool incPsi2S, bool incBkg );
void setMassFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbp, bool cutSideBand=false, bool doSimulFit=false) ;
void setMassGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, bool incJpsi, bool incPsi2S, bool incBkg, bool wantPureSMC=false);
void setMassCutParameters(struct KinCuts& cut, bool incJpsi, bool incPsi2S, bool isMC=false, bool useForCtauFits=false);


bool fitCharmoniaMassModel( RooWorkspace& myws,            // Local Workspace
                            RooWorkspace& inputWorkspace,  // Workspace with all the input RooDatasets
                            struct KinCuts& cut,           // Variable containing all kinematic cuts
                            map<string, string>&  parIni,  // Variable containing all initial parameters
                            struct InputOpt& opt,          // Variable with run information (kept for legacy purpose)
                            string outputDir,              // Path to output directory
                            // Select the type of datasets to fit
                            string DSTAG,                  // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
                            bool isPbp      = false,      // isPbp = false for pp, true for Pbp
                            bool importDS    = true,       // Select if the dataset is imported in the local workspace
                            // Select the type of object to fit
                            bool incJpsi     = true,       // Includes Jpsi model
                            bool incPsi2S    = true,       // Includes Psi(2S) model
                            bool incBkg      = true,       // Includes Background model
                            // Select the fitting options
                            bool doFit       = true,       // Flag to indicate if we want to perform the fit
                            bool cutCtau     = false,      // Apply prompt ctau cuts
                            bool doSimulFit  = false,      // Do simultaneous fit
                            bool wantPureSMC = false,      // Flag to indicate if we want to fit pure signal MC
                            const char* applyCorr ="",     // Flag to indicate if we want corrected dataset and which correction
                            bool loadFitResult = false,    // Load previous fit results
                            string inputFitDir = "",       // Location of the fit results
                            int  numCores    = 2,          // Number of cores used for fitting
                            // Select the drawing options
                            bool setLogScale = true,       // Draw plot with log scale
                            bool incSS       = false,      // Include Same Sign data
                            bool zoomPsi     = false,      // Zoom Psi(2S) peak on extra pad
                            double  binWidth = 0.05,       // Bin width used for plotting
                            bool getMeanPT   = false       // Compute the mean PT (NEED TO FIX)
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
  bool applyWeight_Corr = ( strcmp(applyCorr,"") );
  
  // Define the mass range
  setMassCutParameters(cut, incJpsi, incPsi2S, isMC);
  parIni["invMassNorm"] = Form("RooFormulaVar::%s('( -1.0 + 2.0*( @0 - @1 )/( @2 - @1) )', {%s, mMin[%.6f], mMax[%.6f]})", "invMassNorm", "invMass", cut.dMuon.M.Min, cut.dMuon.M.Max );
  // Apply the ctau cuts to reject non-prompt charmonia
  if (cutCtau) { setCtauCuts(cut, isPbp); }
  
  string COLL = (isPbp ? "Pbp" : "PP" );
  string plotLabelPbp,  plotLabelPP;

  if (doSimulFit || !isPbp) {
    // Set models based on initial parameters
    struct OniaModel model;
    if (!setMassModel(model, parIni, false, incJpsi, incPsi2S, incBkg)) { return false; }

    // Import the local datasets
    double numEntries = 1000000;
    string label = ((DSTAG.find("PP")!=std::string::npos) ? DSTAG.c_str() : Form("%s_%s", DSTAG.c_str(), "PP"));
    if (wantPureSMC) label += "_NoBkg";
    if (applyWeight_Corr) label += Form("_%s",applyCorr);
    string dsName = Form("dOS_%s", label.c_str());
    if (importDS) {
      if ( !(myws.data(dsName.c_str())) ) {
        int importID = importDataset(myws, inputWorkspace, cut, label, cutSideBand);
        if (importID<0) { return false; }
        else if (importID==0) { doFit = false; }
      }
      numEntries = myws.data(dsName.c_str())->sumEntries(); if (numEntries<=0) { doFit = false; }
    }
    else if (doFit && !(myws.data(dsName.c_str()))) { cout << "[ERROR] No local dataset was found to perform the fit!" << endl; return false; }

    // Set global parameters
    setMassGlobalParameterRange(myws, parIni, cut, incJpsi, incPsi2S, incBkg, wantPureSMC);

    // Build the Fit Model
    if (!buildCharmoniaMassModel(myws, model.PP, parIni, false, doSimulFit, incBkg, incJpsi, incPsi2S, numEntries))  { return false; }

    // Define plot names
    if (incJpsi)  { plotLabelPP += Form("_Jpsi_%s", parIni["Model_Jpsi_PP"].c_str());   } 
    if (incPsi2S) { plotLabelPP += Form("_Psi2S_%s", parIni["Model_Psi2S_PP"].c_str()); }
    if (incBkg)   { plotLabelPP += Form("_Bkg_%s", parIni["Model_Bkg_PP"].c_str());     }
    if (wantPureSMC) plotLabelPP +="_NoBkg";
    if (applyWeight_Corr) plotLabelPP +=Form("_%s",applyCorr);
  }

  if (doSimulFit || isPbp) {
    // Set models based on initial parameters
    struct OniaModel model;
    if (!setMassModel(model, parIni, true, incJpsi, incPsi2S, incBkg)) { return false; }

    // Import the local datasets
    double numEntries = 1000000;
    string label = ((DSTAG.find("Pbp")!=std::string::npos) ? DSTAG.c_str() : Form("%s_%s", DSTAG.c_str(), "Pbp"));
    if (wantPureSMC) label += "_NoBkg";
    if (applyWeight_Corr) label += Form("_%s",applyCorr);
    string dsName = Form("dOS_%s", label.c_str());
    if (importDS) {
      if ( !(myws.data(dsName.c_str())) ) {
        int importID = importDataset(myws, inputWorkspace, cut, label, cutSideBand);
        if (importID<0) { return false; }
        else if (importID==0) { doFit = false; }
      }
      numEntries = myws.data(dsName.c_str())->sumEntries(); if (numEntries<=0) { doFit = false; }
    }
    else if (doFit && !(myws.data(dsName.c_str()))) { cout << "[ERROR] No local dataset was found to perform the fit!" << endl; return false; }
      
    // Set global parameters
    setMassGlobalParameterRange(myws, parIni, cut, incJpsi, incPsi2S, incBkg, wantPureSMC);

    // Build the Fit Model
    if (!buildCharmoniaMassModel(myws, model.Pbp, parIni, true, doSimulFit, incBkg, incJpsi, incPsi2S, numEntries))  { return false; }

    // Define plot names
    if (incJpsi)  { plotLabelPbp += Form("_Jpsi_%s", parIni["Model_Jpsi_Pbp"].c_str());   } 
    if (incPsi2S) { plotLabelPbp += Form("_Psi2S_%s", parIni["Model_Psi2S_Pbp"].c_str()); }
    if (incBkg)   { plotLabelPbp += Form("_Bkg_%s", parIni["Model_Bkg_Pbp"].c_str());     }
    if (wantPureSMC) plotLabelPbp += "_NoBkg";
    if (applyWeight_Corr) plotLabelPbp += Form("_%s",applyCorr);
  }

  if (doSimulFit) {
    // Create the combided datasets
    RooCategory* sample = new RooCategory("sample","sample"); sample->defineType("Pbp"); sample->defineType("PP");
    RooDataSet*  combData = new RooDataSet("combData","combined data", *myws.var("invMass"), Index(*sample),
                                           Import("Pbp", *((RooDataSet*)myws.data("dOS_DATA_Pbp"))),
                                           Import("PP",   *((RooDataSet*)myws.data("dOS_DATA_PP")))
                                           );
    myws.import(*sample);

    // Create the combided models
    RooSimultaneous* simPdf = new RooSimultaneous("simPdf", "simultaneous pdf", *sample);
    simPdf->addPdf(*myws.pdf("pdfMASS_Tot_Pbp"), "Pbp"); simPdf->addPdf(*myws.pdf("pdfMASS_Tot_PP"), "PP");
    myws.import(*simPdf);

    // check if we have already done this fit. If yes, do nothing and return true.
    string FileName = "";
    setMassFileName(FileName, (inputFitDir=="" ? outputDir : inputFitDir), DSTAG, (plotLabelPP + plotLabelPbp), cut, isPbp, cutSideBand, doSimulFit);
    if (gSystem->AccessPathName(FileName.c_str()) && inputFitDir!="") {
      cout << "[WARNING] User Input File : " << FileName << " was not found!" << endl;
      if (loadFitResult) return false;
      setMassFileName(FileName, outputDir, DSTAG, (plotLabelPP + plotLabelPbp), cut, isPbp, cutSideBand, doSimulFit);
    }
    bool found =  true; bool skipFit = !doFit;
    RooArgSet *newpars = myws.pdf("simPdf")->getParameters(*(myws.var("invMass")));
    myws.saveSnapshot("simPdf_parIni", *newpars, kTRUE);
    found = found && isFitAlreadyFound(newpars, FileName, "simPdf");
    if (loadFitResult) {
      if ( loadPreviousFitResult(myws, FileName, DSTAG, false, cutSideBand) ) { skipFit = true; } else { skipFit = false; }
      if ( loadPreviousFitResult(myws, FileName, DSTAG, true, cutSideBand)  ) { skipFit = true; } else { skipFit = false; }
      if (skipFit) { cout << "[INFO] This simultaneous mass fit was already done, so I'll load the fit results." << endl; }
      myws.saveSnapshot("simPdf_parLoad", *newpars, kTRUE);
    } else if (found) {
      cout << "[INFO] This simultaneous mass fit was already done, so I'll just go to the next one." << endl;
      return true;
    }

    // Do the simultaneous fit
    if (skipFit==false) {
      RooFitResult* fitResult = simPdf->fitTo(*combData, Offset(kTRUE), Extended(kTRUE), NumCPU(numCores), Range("MassWindow"), Save()); //, Minimizer("Minuit2","Migrad")
      fitResult->Print("v");
      myws.import(*fitResult, "fitResult_simPdf"); 
      // Create the output files
      int nBins = min(int( round((cut.dMuon.M.Max - cut.dMuon.M.Min)/binWidth) ), 1000);
      drawMassPlot(myws, outputDir, opt, cut, parIni, plotLabelPP, DSTAG, false, incJpsi, incPsi2S, incBkg, cutCtau, doSimulFit, false, setLogScale, incSS, zoomPsi, nBins, getMeanPT);
      drawMassPlot(myws, outputDir, opt, cut, parIni, plotLabelPbp, DSTAG, true, incJpsi, incPsi2S, incBkg, cutCtau, doSimulFit, false, setLogScale, incSS, zoomPsi, nBins, getMeanPT);
      // Save the results
      string FileName = ""; setMassFileName(FileName, outputDir, DSTAG, (plotLabelPP + plotLabelPbp), cut, isPbp, cutSideBand, doSimulFit);
      myws.saveSnapshot("simPdf_parFit", *newpars, kTRUE);
      saveWorkSpace(myws, Form("%smass%s/%s/result", outputDir.c_str(), (cutSideBand?"SB":""), DSTAG.c_str()), FileName);
      // Delete the objects used during the simultaneous fit
      delete sample; delete combData; delete simPdf;
    }
  }
  else {
    // Define pdf and plot names
    string pdfName = Form("pdfMASS_Tot_%s", COLL.c_str());
    string plotLabel = (isPbp ? plotLabelPbp : plotLabelPP);

    // Import the local datasets
    string label = ((DSTAG.find(COLL.c_str())!=std::string::npos) ? DSTAG.c_str() : Form("%s_%s", DSTAG.c_str(), COLL.c_str()));
    if (wantPureSMC) label += "_NoBkg";
    if (applyWeight_Corr) label += Form("_%s",applyCorr);
    string dsName = Form("dOS_%s", label.c_str());
      
    // check if we have already done this fit. If yes, do nothing and return true.
    string FileName = "";
    setMassFileName(FileName, (inputFitDir=="" ? outputDir : inputFitDir), DSTAG, plotLabel, cut, isPbp, cutSideBand);
    if (gSystem->AccessPathName(FileName.c_str()) && inputFitDir!="") {
      cout << "[WARNING] User Input File : " << FileName << " was not found!" << endl;
      if (loadFitResult) return false;
      setMassFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbp, cutSideBand);
    }
    bool found =  true; bool skipFit = !doFit;
    RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(*(myws.var("invMass")));
    found = found && isFitAlreadyFound(newpars, FileName, pdfName.c_str());
    if (loadFitResult) {
        if ( loadPreviousFitResult(myws, FileName, DSTAG, isPbp, cutSideBand) ) { skipFit = true; } else { skipFit = false; } 
        if (skipFit) { cout << "[INFO] This mass fit was already done, so I'll load the fit results." << endl; }
        myws.saveSnapshot(Form("%s_parLoad", pdfName.c_str()), *newpars, kTRUE);
    } else if (found) {
      cout << "[INFO] This mass fit was already done, so I'll just go to the next one." << endl;
      return true;
    }

    // Fit the Datasets
    if (skipFit==false) {
      bool isWeighted = myws.data(dsName.c_str())->isWeighted();
      RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), SumW2Error(isWeighted), Range(cutSideBand ? parIni["BkgMassRange_FULL_Label"].c_str() : "MassWindow"), NumCPU(numCores), Save());
      fitResult->Print("v"); 
      myws.import(*fitResult, Form("fitResult_%s", pdfName.c_str())); 
      // Create the output files
      int nBins = min(int( round((cut.dMuon.M.Max - cut.dMuon.M.Min)/binWidth) ), 1000);
      drawMassPlot(myws, outputDir, opt, cut, parIni, plotLabel, DSTAG, isPbp, incJpsi, incPsi2S, incBkg, cutCtau, doSimulFit, wantPureSMC, setLogScale, incSS, zoomPsi, nBins, getMeanPT);
      // Save the results
      string FileName = ""; setMassFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbp, cutSideBand);
      myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()), *newpars, kTRUE);
      saveWorkSpace(myws, Form("%smass%s/%s/result", outputDir.c_str(), (cutSideBand?"SB":""), DSTAG.c_str()), FileName);
    }
  }

  return true;
};


void setCtauCuts(struct KinCuts& cut, bool isPbp) 
{
  if (cut.dMuon.AbsRap.Max<=1.6) {
    cut.dMuon.ctauCut = "( ctau < (0.012 + (0.23/pt)) )";
  }
  if (cut.dMuon.AbsRap.Min>=1.6) {
    cut.dMuon.ctauCut = "( ctau < (0.014 + (0.28/pt)) )";
  }

  return;
};


bool setMassModel( struct OniaModel& model, map<string, string>&  parIni, bool isPbp, bool incJpsi, bool incPsi2S, bool incBkg)
{
  if (isPbp && incBkg) {
    if (parIni.count("Model_Bkg_Pbp")>0) {
      model.Pbp.Bkg.Mass = MassModelDictionary[parIni["Model_Bkg_Pbp"]];
      if (model.Pbp.Bkg.Mass==MassModel(0)) {
        cout << "[ERROR] The background model: " << parIni["Model_Bkg_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Background mass model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  }
  if (isPbp && incJpsi) {
    if (parIni.count("Model_Jpsi_Pbp")>0) {
      model.Pbp.Jpsi.Mass = MassModelDictionary[parIni["Model_Jpsi_Pbp"]];
      if (model.Pbp.Jpsi.Mass==MassModel(0)) {
        cout << "[ERROR] The Jpsi model: " << parIni["Model_Jpsi_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Jpsi mass model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  }
  if (isPbp && incPsi2S) {
    if (parIni.count("Model_Psi2S_Pbp")>0) {
      model.Pbp.Psi2S.Mass = MassModelDictionary[parIni["Model_Psi2S_Pbp"]];
      if (model.Pbp.Psi2S.Mass==MassModel(0)) {
        cout << "[ERROR] The psi2S model: " << parIni["Model_Psi2S_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] psi(2S) mass model for Pbp was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbp && incBkg) {
    if (parIni.count("Model_Bkg_PP")>0) {
      model.PP.Bkg.Mass = MassModelDictionary[parIni["Model_Bkg_PP"]];
      if (model.PP.Bkg.Mass==MassModel(0)) {
        cout << "[ERROR] The background model: " << parIni["Model_Bkg_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Background mass model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbp && incJpsi) {
    if (parIni.count("Model_Jpsi_PP")>0) {
      model.PP.Jpsi.Mass = MassModelDictionary[parIni["Model_Jpsi_PP"]];
      if (model.PP.Jpsi.Mass==MassModel(0)) {
        cout << "[ERROR] The Jpsi model: " << parIni["Model_Jpsi_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Jpsi mass model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbp && incPsi2S) {
    if (parIni.count("Model_Psi2S_PP")>0) {
      model.PP.Psi2S.Mass = MassModelDictionary[parIni["Model_Psi2S_PP"]];
      if (model.PP.Psi2S.Mass==MassModel(0)) {
        cout << "[ERROR] The psi2S model: " << parIni["Model_Psi2S_Pbp"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] psi(2S) mass model for PP was not found in the initial parameters!" << endl; return false;
    }
  }

  return true;
};


void setMassGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, bool incJpsi, bool incPsi2S, bool incBkg,  bool wantPureSMC)
{
  /*if (wantPureSMC)
    {
      if (incPsi2S)
        {
          if (cut.dMuon.AbsRap.Min >= 1.6) {
            myws.var("invMass")->setRange("MassWindow", cut.dMuon.M.Min, 3.95);
            parIni["MassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", cut.dMuon.M.Min, 3.95);
          }
          else { 
            myws.var("invMass")->setRange("MassWindow", cut.dMuon.M.Min, 3.85);
            parIni["MassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", cut.dMuon.M.Min, 3.85);
          }
        }
      if (incJpsi)
        {
          if (cut.dMuon.AbsRap.Min >= 1.6) {
            myws.var("invMass")->setRange("MassWindow", 2.6, 3.32);
            parIni["MassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", 2.6, 3.32);
          }
          else {
            myws.var("invMass")->setRange("MassWindow", 2.6, 3.26);
            parIni["MassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", 2.6, 3.26);
          }
        }
    }
  */
  myws.var("invMass")->setRange("InclusiveMassRegion", cut.dMuon.M.Min, cut.dMuon.M.Max);
  myws.var("invMass")->setRange("FullWindow", cut.dMuon.M.Min, cut.dMuon.M.Max);
  myws.var("invMass")->setRange("MassWindow", cut.dMuon.M.Min, cut.dMuon.M.Max);
  parIni["MassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", cut.dMuon.M.Min, cut.dMuon.M.Max);
  if (incJpsi) {
    myws.var("invMass")->setRange("JpsiWindow", 2.9, 3.3);
    parIni["JpsiMassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", 2.9, 3.3);
  }
  if (incPsi2S) {
    myws.var("invMass")->setRange("Psi2SWindow", 3.45, 3.85);
    parIni["Psi2SMassRange_Cut"] = Form("(invMass>%.6f && invMass<%.6f)", 3.45, 3.85);
  }
  if (incBkg) {
    myws.var("invMass")->setRange("SideBandMID_FULL",  ((cut.dMuon.M.Min<3.3)?3.3:cut.dMuon.M.Min), ((cut.dMuon.M.Max>3.45)?3.5:cut.dMuon.M.Max));
    myws.var("invMass")->setRange("SideBandMID_JPSI",  ((cut.dMuon.M.Min<3.3)?3.3:cut.dMuon.M.Min), ((cut.dMuon.M.Max>3.45)?3.5:cut.dMuon.M.Max));
    myws.var("invMass")->setRange("SideBandMID_PSI2S", ((cut.dMuon.M.Min<3.35)?3.35:cut.dMuon.M.Min), ((cut.dMuon.M.Max>3.5)?3.5:cut.dMuon.M.Max));
    parIni["BkgMassRange_FULL_Label"]  = "SideBandMID_FULL";
    parIni["BkgMassRange_JPSI_Label"]  = "SideBandMID_JPSI";
    parIni["BkgMassRange_PSI2S_Label"] = "SideBandMID_PSI2S";
    if (cut.dMuon.M.Min < 2.8) {
      myws.var("invMass")->setRange("SideBandBOT_FULL", cut.dMuon.M.Min, 2.8);
      myws.var("invMass")->setRange("SideBandBOT_JPSI", ((cut.dMuon.M.Min<2.0)?2.0:cut.dMuon.M.Min), 2.8);
      parIni["BkgMassRange_FULL_Label"] = parIni["BkgMassRange_FULL_Label"] + "," + "SideBandBOT_FULL";
      parIni["BkgMassRange_JPSI_Label"] = parIni["BkgMassRange_JPSI_Label"] + "," + "SideBandBOT_JPSI";
    }
    if (cut.dMuon.M.Max > 3.9) {
      myws.var("invMass")->setRange("SideBandTOP_FULL", 3.85, cut.dMuon.M.Max);
      myws.var("invMass")->setRange("SideBandTOP_PSI2S", 3.85, ((cut.dMuon.M.Max>4.2)?4.2:cut.dMuon.M.Max));
      parIni["BkgMassRange_FULL_Label"] = parIni["BkgMassRange_FULL_Label"] + "," + "SideBandTOP_FULL";
      parIni["BkgMassRange_PSI2S_Label"] = parIni["BkgMassRange_PSI2S_Label"] + "," + "SideBandTOP_PSI2S";
    }
    parIni["BkgMassRange_FULL_Cut"]  = Form("(%.6f < invMass && invMass < %.6f)",       cut.dMuon.M.Min,       cut.dMuon.M.Max);
    parIni["BkgMassRange_FULL_Cut"]  = parIni["BkgMassRange_FULL_Cut"]  + "&&" + "((2.0 < invMass && invMass < 2.8) || (3.3 < invMass && invMass < 3.45) || (3.85 < invMass && invMass < 5.0))";
    parIni["BkgMassRange_JPSI_Cut"]  = parIni["BkgMassRange_FULL_Cut"]  + "&&" + "(2.0 < invMass && invMass < 3.45)";
    parIni["BkgMassRange_PSI2S_Cut"] = parIni["BkgMassRange_FULL_Cut"] + "&&" + "(3.35 < invMass && invMass < 5.0)";
    parIni["BkgMassRange_FULL_Cut"]  = "("+parIni["BkgMassRange_FULL_Cut"]+")";
    parIni["BkgMassRange_JPSI_Cut"]  = "("+parIni["BkgMassRange_JPSI_Cut"]+")";
    parIni["BkgMassRange_PSI2S_Cut"] = "("+parIni["BkgMassRange_PSI2S_Cut"]+")";
  }
  

  return;
};

void setMassFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbp, bool cutSideBand, bool doSimulFit) 
{
  if (TAG.find("_")!=std::string::npos) TAG.erase(TAG.find("_"));
  FileName = Form("%smass%s/%s/result/FIT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), (cutSideBand?"SB":""), TAG.c_str(), "MASS", TAG.c_str(), (doSimulFit ? "COMB" : (isPbp?"Pbp":"PP")), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
};


void setMassCutParameters(struct KinCuts& cut, bool incJpsi, bool incPsi2S, bool isMC, bool useForCtauFits)
{
  // Define the mass range
  if (cut.dMuon.M.Max==5 && cut.dMuon.M.Min==2) { 
    // Default mass values, means that the user did not specify a mass range
    if ( incJpsi && !incPsi2S) {
      if (isMC && !useForCtauFits){
        cut.dMuon.M.Min = 2.2;
        cut.dMuon.M.Max = 3.5;
      } else {
        cut.dMuon.M.Min = 2.6;
        cut.dMuon.M.Max = 3.5;
      }
    }
    else if ( !incJpsi && incPsi2S) {
      if(isMC && !useForCtauFits) {
        cut.dMuon.M.Min = 3.4;
        cut.dMuon.M.Max = 4.6;
      } else {
        cut.dMuon.M.Min = 3.35;
        cut.dMuon.M.Max = 4.2;
      }
    }
    else if ( incJpsi && incPsi2S) {
      cut.dMuon.M.Min = 2.2;
      cut.dMuon.M.Max = 4.2;
    }
    else {
      cut.dMuon.M.Min = 3.35;
      cut.dMuon.M.Max = 4.2;
    }
  }
  cout << "[INFO] Setting mass range to min: " << cut.dMuon.M.Min << " and max " << cut.dMuon.M.Max << endl;

  return;
};


#endif // #ifndef fitCharmoniaMassModel_C
