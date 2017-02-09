#ifndef fitCharmonia_C
#define fitCharmonia_C

#include "Utilities/initClasses.h"
#include "fitCharmoniaMassModel.C" 
#include "fitCharmoniaCtauModel.C"
#include "fitCharmoniaCtauErrModel.C"
#include "fitCharmoniaCtauTrueModel.C"
#include "fitCharmoniaCtauMassModel.C"
#include "fitCharmoniaCtauResModel.C"

void setOptions(struct InputOpt* opt);

bool fitCharmonia( RooWorkspace&  inputWorkspace,  // Workspace with all the input RooDatasets
		   struct KinCuts cut,             // Variable containing all kinematic cuts
		   map<string, string>  parIni,    // Variable containing all initial parameters
		   string outputDir,               // Path to output directory
                   // Select the type of datasets to fit
		   string DSTAG,                   // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
		   bool isPbp      = false,       // isPbp = false for pp, true for Pbp
                   // Select the type of object to fit
                   bool fitMass      = true,       // Fit mass distribution
                   bool fitCtau      = false,      // Fit ctau distribution
                   bool fitCtauTrue  = false,      // Fit ctau truth MC distribution
                   bool incJpsi      = true,       // Includes Jpsi model
                   bool incPsi2S     = true,       // Includes Psi(2S) model
                   bool incBkg       = true,       // Includes Background model
                   bool incPrompt    = true,       // Includes Prompt ctau model
                   bool incNonPrompt = false,      // Includes NonPrompt ctau model
                   bool doCtauErrPDF = false,      // If yes, it builds the Ctau Error PDFs from data
                   // Select the fitting options
                   bool cutCtau      = false,      // Apply prompt ctau cuts
                   bool doSimulFit   = false,      // Do simultaneous fit
                   bool wantPureSMC  = false,      // Flag to indicate if we want to fit pure signal MC
                   map<string, string> inputFitDir={},// User-defined Location of the fit results
                   const char* applyCorr ="",      // Flag to indicate if we want corrected dataset and which correction
                   int  numCores     = 2,          // Number of cores used for fitting
                   // Select the drawing options
                   bool setLogScale  = true,       // Draw plot with log scale
                   bool incSS        = false,      // Include Same Sign data
                   bool zoomPsi      = false,      // Zoom Psi(2S) peak on extra pad
                   map<string, double> binWidth={},// Bin width used for plotting
                   bool getMeanPT    = false       // Compute the mean PT (NEED TO FIX)
		   )  
{

  RooMsgService::instance().getStream(0).removeTopic(Caching);  
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  RooWorkspace     myws("workspace", "local workspace");

  // Preventing issues in the fitter
  cutCtau = (cutCtau && !fitCtau && !fitCtauTrue);
  getMeanPT = false;
  wantPureSMC = (DSTAG.find("MC")!=std::string::npos && wantPureSMC);
  bool isMC = (DSTAG.find("MC")!=std::string::npos);

  // Setting default user-defined input fit directories ( "" means use current working directory )
  if (inputFitDir.count("MASS")==0)     { inputFitDir["MASS"]     = ""; }
  if (inputFitDir.count("CTAUTRUE")==0) { inputFitDir["CTAUTRUE"] = ""; }
  if (inputFitDir.count("CTAUERR")==0)  { inputFitDir["CTAUERR"]  = ""; }
  if (inputFitDir.count("CTAURES")==0)  { inputFitDir["CTAURES"]  = ""; }
  if (inputFitDir.count("CTAUSB")==0)   { inputFitDir["CTAUSB"]   = ""; }
  // Setting default user-defined bin width
  if (binWidth.count("MASS")==0)     { binWidth["MASS"]     = 0.05; }
  if (binWidth.count("CTAUTRUE")==0) { binWidth["CTAUTRUE"] = 0.05; }
  if (binWidth.count("CTAUERR")==0)  { binWidth["CTAUERR"]  = 0.05; }
  if (binWidth.count("CTAUSB")==0)   { binWidth["CTAUSB"]   = 0.05; }
  if (binWidth.count("CTAURES")==0)  { binWidth["CTAURES"]  = 0.05; }
  if (binWidth.count("CTAU")==0)     { binWidth["CTAU"]     = 0.05; }
  binWidth["CTAUERRFORCUT"]  = 0.0025;

  if (isPbp==false) {
    cut.Centrality.Start = 0;
    cut.Centrality.End = 200;
  }

  // Setting run information
  struct InputOpt opt; setOptions(&opt);

  // Starting the fits for each variable (where the magic starts)
  if (fitMass && !fitCtau && !fitCtauTrue && !doCtauErrPDF) {
    
    // Setting extra input information needed by each fitter
    double ibWidth = binWidth["MASS"];
    string iFitDir = inputFitDir["MASS"];
    bool loadFitResult = false;
    bool doFit = true;
    bool importDS = true;

    if ( !fitCharmoniaMassModel( myws, inputWorkspace, cut, parIni, opt, outputDir, 
                                 DSTAG, isPbp, importDS,
                                 incJpsi, incPsi2S, incBkg, 
                                 doFit, cutCtau, doSimulFit, wantPureSMC, applyCorr, loadFitResult, iFitDir, numCores, 
                                 setLogScale, incSS, zoomPsi, ibWidth, getMeanPT 
                                 ) 
         ) { return false; }
  }

  if (fitCtauTrue && !fitCtau && !fitMass && !doCtauErrPDF) {

    // Setting extra input information needed by each fitter
    double ibWidth = binWidth["CTAUTRUE"];
    string iFitDir = inputFitDir["CTAUTRUE"];
    bool loadFitResult = false;
    bool doFit = true;
    bool importDS = true;
    bool incResol = false;
    
    if ( !fitCharmoniaCtauTrueModel( myws, inputWorkspace, cut, parIni, opt, outputDir, 
                                     DSTAG, isPbp, importDS, 
                                     incJpsi, incPsi2S, incResol, 
                                     doFit, wantPureSMC, loadFitResult, iFitDir, numCores, 
                                     setLogScale, incSS, ibWidth
                                     ) 
         ) { return false; }
  }

  if (doCtauErrPDF && !fitCtau && !fitCtauTrue && !fitMass) {

    // Setting extra input information needed by each fitter
    bool loadFitResult = false;
    bool doFit = true;
    bool importDS = true;

    if ( !fitCharmoniaCtauErrModel( myws, inputWorkspace, cut, parIni, opt, outputDir, 
                                    DSTAG, isPbp, importDS, 
                                    incJpsi, incPsi2S, incBkg, 
                                    doFit, wantPureSMC, loadFitResult, inputFitDir, numCores, 
                                    setLogScale, incSS, binWidth
                                    ) 
         ) { return false; }
  }

  //if (fitCtau && !doCtauErrPDF && !fitCtauTrue && !fitMass && (incJpsi!=incBkg) && !isMC) {
  if (fitCtau && !doCtauErrPDF && !fitCtauTrue && !fitMass && (incPsi2S!=incBkg) && !isMC) {

    // Setting extra input information needed by each fitter
    bool loadFitResult = false;
    bool doFit = true;
    bool importDS = true;

    if ( !fitCharmoniaCtauModel( myws, inputWorkspace, cut, parIni, opt, outputDir, 
                           DSTAG, isPbp, importDS, 
                           incJpsi, incPsi2S, incBkg, incPrompt, incNonPrompt, 
                           doFit, wantPureSMC, loadFitResult, inputFitDir, numCores, 
                           setLogScale, incSS, binWidth
                           ) 
  
         ) { return false; }
  }

  if (fitCtau && !doCtauErrPDF && !fitCtauTrue && !fitMass && !incBkg && isMC) {

    // Setting extra input information needed by each fitter
    bool loadFitResult = false;
    bool doFit = true;
    bool importDS = true;

    if ( !fitCharmoniaCtauResModel( myws, inputWorkspace, cut, parIni, opt, outputDir, 
                                    DSTAG, isPbp, importDS, 
                                    incJpsi, incPsi2S,
                                    doFit, wantPureSMC, loadFitResult, inputFitDir, numCores, 
                                    setLogScale, incSS, binWidth
                                    ) 
         ) {return false; }
  }

  if (fitCtau && fitMass && !doCtauErrPDF && !fitCtauTrue ) {

    // Setting extra input information needed by each fitter
   
    if ( !fitCharmoniaCtauMassModel( myws, inputWorkspace, cut, parIni, opt, outputDir, 
                                     DSTAG, isPbp,
                                     incJpsi, incPsi2S, 
                                     inputFitDir, numCores,
                                     setLogScale, incSS, binWidth
                                     ) 
         ) { return false; }
  }

  return true;
};


void setOptions(struct InputOpt* opt) 
{
  opt->pp.RunNb.Start   = 262157; opt->Pbp.RunNb.Start = 210676;
  opt->pp.RunNb.End     = 262328; opt->Pbp.RunNb.End   = 211631;
  opt->pp.TriggerBit    = (int) PP::HLT_HIL1DoubleMu0_v1; //HLT_PAL1DoubleMuOpen_v1; 
  opt->Pbp.TriggerBit  = (int) HI::HLT_PAL1DoubleMuOpen_v1; //HLT_HIL1DoubleMu0_v1; 
  return;
};


#endif // #ifndef fitCharmonia_C
