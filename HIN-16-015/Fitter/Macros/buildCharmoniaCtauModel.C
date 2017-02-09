#ifndef buildCharmoniaCtauModel_C
#define buildCharmoniaCtauModel_C

#include "Utilities/initClasses.h"

void fixCtauParPsi2StoJpsi(map<string, string>& parIni, bool isPbp);
void setCtauDefaultParameters(map<string, string> &parIni, bool isPbp, double numEntries);
bool defineCtauResolModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp); 
bool addSignalCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp); 
bool addBackgroundCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp);

bool buildCharmoniaCtauModel(RooWorkspace& ws, struct CharmModel model, map<string, string>  parIni, string dsName,
                             bool isPbp,                 // Determine if we are working with Pbp (True) or PP (False)
                             bool incBkg,                 // Include background model
                             bool incJpsi,                // Include Jpsi model
                             bool incPsi2S,               // Include Psi(2S) model
                             bool incPrompt,              // Include Prompt models
                             bool incNonPrompt,           // Include NonPrompt models
                             double  numEntries = 300000. // Number of entries in the dataset
                             )
{
  // If the initial parameters are empty, set defaul parameter values
  setCtauDefaultParameters(parIni, isPbp, numEntries);
  
  // Fix all psi2S parameters to jpsi
  if (incJpsi && incPsi2S) {
    fixCtauParPsi2StoJpsi(parIni, isPbp);
  }

  // C r e a t e   m o d e l 

  bool fitMass = false;
  if ( ws.pdf(Form("pdfMASS_Tot_%s", (isPbp?"Pbp":"PP"))) ) { fitMass = true; } 

  string pdfName     = "pdfCTAU";
  if (fitMass) { pdfName = "pdfCTAUMASS"; }
  bool isMC = (dsName.find("MC")!=std::string::npos);
  bool incCtauErrPDF = true;
  
  if (incJpsi) {
    if (incPrompt) {
      if(!defineCtauResolModel(ws, "JpsiPR", model.CtauRes, parIni, isPbp)) { cout << "[ERROR] Defining the Prompt Ctau Resolution Model failed" << endl; return false; }
      if(!addSignalCtauModel(ws, "JpsiPR", model.Jpsi.Ctau.Prompt, parIni, isPbp)) { cout << "[ERROR] Adding Prompt Jpsi Ctau Model failed" << endl; return false; }
    }
    if (incNonPrompt) {
      if(!defineCtauResolModel(ws, "JpsiNoPR", model.CtauRes, parIni, isPbp)) { cout << "[ERROR] Defining the Non-Prompt Ctau Resolution Model failed" << endl; return false; }
      if(!addSignalCtauModel(ws, "JpsiNoPR", model.Jpsi.Ctau.NonPrompt, parIni, isPbp)) { cout << "[ERROR] Adding NonPrompt Jpsi Ctau Model failed" << endl; return false; }
    }
    if (incPrompt && !incNonPrompt) {
      if (incCtauErrPDF) {
        RooProdPdf pdfJpsi(Form("pdfCTAU_JpsiPR_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Jpsi_%s", (isPbp?"Pbp":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_JpsiPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfJpsi);
      } else {
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_JpsiPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUCOND_JpsiPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      }
      if ( ws.pdf(Form("pdfMASS_Jpsi_%s", (isPbp?"Pbp":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_JpsiPR_%s", (isPbp?"Pbp":"PP")),
                       Form("%s_JpsiPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbp?"Pbp":"PP")),
                       Form("pdfMASS_Jpsi_%s",(isPbp?"Pbp":"PP"))
                       ));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_JpsiPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("%s_JpsiPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))
                      ));
    }
    if (!incPrompt && incNonPrompt) {
      if (incCtauErrPDF) {
        RooProdPdf pdfJpsi(Form("pdfCTAU_JpsiNoPR_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Jpsi_%s", (isPbp?"Pbp":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfJpsi);
      } else {
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_JpsiNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUCOND_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      }
      if ( ws.pdf(Form("pdfMASS_Jpsi_%s", (isPbp?"Pbp":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")),
                         Form("%s_JpsiNoPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbp?"Pbp":"PP")),
                         Form("pdfMASS_Jpsi_%s",(isPbp?"Pbp":"PP"))
                         ));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_JpsiNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("%s_JpsiNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))
                      ));
    }
    if (incPrompt && incNonPrompt) {
      if (incCtauErrPDF) {
        RooProdPdf pdfJpsiPR(Form("pdfCTAU_JpsiPR_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Jpsi_%s", (isPbp?"Pbp":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_JpsiPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau")) )
                       );
        ws.import(pdfJpsiPR);
        RooProdPdf pdfJpsiNoPR(Form("pdfCTAU_JpsiNoPR_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Jpsi_%s", (isPbp?"Pbp":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfJpsiNoPR);
      } else {
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_JpsiPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUCOND_JpsiPR_%s", (isPbp?"Pbp":"PP"))
                        ));
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_JpsiNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUCOND_JpsiPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      }
      ws.factory( parIni[Form("b_Jpsi_%s", (isPbp?"Pbp":"PP"))].c_str() );
      if (incCtauErrPDF) {
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU_Jpsi_%s", (isPbp?"Pbp":"PP")),
                        Form("b_Jpsi_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAU_JpsiNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAU_JpsiPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      } else {
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUCOND_Jpsi_%s", (isPbp?"Pbp":"PP")),
                        Form("b_Jpsi_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUCOND_JpsiNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUCOND_JpsiPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      }
      if ( ws.pdf(Form("pdfMASS_Jpsi_%s", (isPbp?"Pbp":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_JpsiPR_%s", (isPbp?"Pbp":"PP")),
                        Form("%s_JpsiPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbp?"Pbp":"PP")),
                        Form("pdfMASS_Jpsi_%s",(isPbp?"Pbp":"PP"))
                        ));
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("%s_JpsiNoPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbp?"Pbp":"PP")),
                        Form("pdfMASS_Jpsi_%s",(isPbp?"Pbp":"PP"))
                        ));
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUMASS_Jpsi_%s", (isPbp?"Pbp":"PP")),
                        Form("b_Jpsi_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUMASS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUMASS_JpsiPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Jpsi_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("%s_Jpsi_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))
                      ));
    }
  }
  if (incPsi2S) {
    if (incPrompt) {
      if(!defineCtauResolModel(ws, "Psi2SPR", model.CtauRes, parIni, isPbp)) { cout << "[ERROR] Defining the Prompt Ctau Resolution Model failed" << endl; return false; }
      if(!addSignalCtauModel(ws, "Psi2SPR", model.Psi2S.Ctau.Prompt, parIni, isPbp)) { cout << "[ERROR] Adding Prompt Psi2S Ctau Model failed" << endl; return false; }
    }
    if (incNonPrompt) {
      if(!defineCtauResolModel(ws, "Psi2SNoPR", model.CtauRes, parIni, isPbp)) { cout << "[ERROR] Defining the Non-Prompt Ctau Resolution Model failed" << endl; return false; }
      if(!addSignalCtauModel(ws, "Psi2SNoPR", model.Psi2S.Ctau.NonPrompt, parIni, isPbp)) { cout << "[ERROR] Adding NonPrompt Psi2S Ctau Model failed" << endl; return false; }
    }
    if (incPrompt && !incNonPrompt) {
      if (incCtauErrPDF) {
        RooProdPdf pdfPsi2S(Form("pdfCTAU_Psi2SPR_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Psi2S_%s", (isPbp?"Pbp":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_Psi2SPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfPsi2S);
      } else {
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_Psi2SPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUCOND_Psi2SPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      }
      if ( ws.pdf(Form("pdfMASS_Psi2S_%s", (isPbp?"Pbp":"PP"))) ){
         ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_Psi2SPR_%s", (isPbp?"Pbp":"PP")),
                         Form("%s_Psi2SPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbp?"Pbp":"PP")),
                         Form("pdfMASS_Psi2S_%s",(isPbp?"Pbp":"PP"))
                         ));
       }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Psi2SPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("%s_Psi2SPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))
                      ));
    }
    if (!incPrompt && incNonPrompt) {
      if (incCtauErrPDF) {
        RooProdPdf pdfPsi2S(Form("pdfCTAU_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Psi2S_%s", (isPbp?"Pbp":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfPsi2S);
      } else {
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUCOND_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      }
      if ( ws.pdf(Form("pdfMASS_Psi2S_%s", (isPbp?"Pbp":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("%s_Psi2SNoPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbp?"Pbp":"PP")),
                        Form("pdfMASS_Psi2S_%s",(isPbp?"Pbp":"PP"))
                        ));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Psi2SNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("%s_Psi2SNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))
                      ));
    }
    if (incPrompt && incNonPrompt) {
      if (incCtauErrPDF) {
        RooProdPdf pdfPsi2SPR(Form("pdfCTAU_Psi2SPR_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Psi2S_%s", (isPbp?"Pbp":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_Psi2SPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfPsi2SPR);
        RooProdPdf pdfPsi2SNoPR(Form("pdfCTAU_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Psi2S_%s", (isPbp?"Pbp":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfPsi2SNoPR);
      } else {
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_Psi2SPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUCOND_Psi2SPR_%s", (isPbp?"Pbp":"PP"))
                        ));
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUCOND_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      }
      ws.factory( parIni[Form("b_Psi2S_%s", (isPbp?"Pbp":"PP"))].c_str() );
      if (incCtauErrPDF) {
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU_Psi2S_%s", (isPbp?"Pbp":"PP")),
                        Form("b_Psi2S_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAU_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAU_Psi2SPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      } else {
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUCOND_Psi2S_%s", (isPbp?"Pbp":"PP")),
                        Form("b_Psi2S_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUCOND_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUCOND_Psi2SPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      }
      if ( ws.pdf(Form("pdfMASS_Psi2S_%s", (isPbp?"Pbp":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_Psi2SPR_%s", (isPbp?"Pbp":"PP")),
                        Form("%s_Psi2SPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbp?"Pbp":"PP")),
                        Form("pdfMASS_Psi2S_%s",(isPbp?"Pbp":"PP"))
                        ));
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("%s_Psi2SNoPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbp?"Pbp":"PP")),
                        Form("pdfMASS_Psi2S_%s",(isPbp?"Pbp":"PP"))
                        ));
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUMASS_Psi2S_%s", (isPbp?"Pbp":"PP")),
                        Form("b_Psi2S_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUMASS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUMASS_Psi2SPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Psi2S_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("%s_Psi2S_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))
                      ));
    }
  }  
  if (incBkg) {
    if (incPrompt) {
      if(!defineCtauResolModel(ws, "BkgPR", model.CtauRes, parIni, isPbp)) { cout << "[ERROR] Defining the Prompt Ctau Resolution Model failed" << endl; return false; }
      if(!addBackgroundCtauModel(ws, "BkgPR", model.Bkg.Ctau.Prompt, parIni, isPbp)) { cout << "[ERROR] Adding Prompt Bkg Ctau Model failed" << endl; return false; }
    }
    if (incNonPrompt) {
      if(!defineCtauResolModel(ws, "BkgNoPR", model.CtauRes, parIni, isPbp)) { cout << "[ERROR] Defining the Non-Prompt Ctau Resolution Model failed" << endl; return false; }
      if(!addBackgroundCtauModel(ws, "BkgNoPR", model.Bkg.Ctau.NonPrompt, parIni, isPbp)) { cout << "[ERROR] Adding NonPrompt Bkg Ctau Model failed" << endl; return false; }
    }
    if (incPrompt && !incNonPrompt) {
      if (incCtauErrPDF) {
        RooProdPdf pdf(Form("pdfCTAU_BkgPR_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Bkg_%s", (isPbp?"Pbp":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_BkgPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdf);
      }
      if ( ws.pdf(Form("pdfMASS_Bkg_%s", (isPbp?"Pbp":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_BkgPR_%s", (isPbp?"Pbp":"PP")),
                        Form("%s_BkgPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbp?"Pbp":"PP")),
                        Form("pdfMASS_Bkg_%s",(isPbp?"Pbp":"PP"))
                        ));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_BkgPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("%s_BkgPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("N_Bkg_%s", (isPbp?"Pbp":"PP"))
                      ));
    }
    if (!incPrompt && incNonPrompt) {
      if (incCtauErrPDF) {
        RooProdPdf pdf(Form("pdfCTAU_BkgNoPR_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Bkg_%s", (isPbp?"Pbp":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_BkgNoPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdf);
      }
      if ( ws.pdf(Form("pdfMASS_Bkg_%s", (isPbp?"Pbp":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_BkgNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("%s_BkgNoPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbp?"Pbp":"PP")),
                        Form("pdfMASS_Bkg_%s",(isPbp?"Pbp":"PP"))
                        ));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_BkgNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("%s_BkgNoPR_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("N_Bkg_%s", (isPbp?"Pbp":"PP"))
                      ));
    }
    if (incPrompt && incNonPrompt) {
      if (incCtauErrPDF) {
        if ( !ws.pdf(Form("pdfCTAUERR_Bkg_%s", (isPbp?"Pbp":"PP"))) ) { return false; }
        RooProdPdf pdfPR(Form("pdfCTAU_BkgPR_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Bkg_%s", (isPbp?"Pbp":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_BkgPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfPR);
        RooProdPdf pdfNoPR(Form("pdfCTAU_BkgNoPR_%s", (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_Bkg_%s", (isPbp?"Pbp":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_BkgNoPR_%s", (isPbp?"Pbp":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfNoPR);
      }
      ws.factory( parIni[Form("b_Bkg_%s", (isPbp?"Pbp":"PP"))].c_str() );
      if (incCtauErrPDF) {
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU_Bkg_%s", (isPbp?"Pbp":"PP")),
                        Form("b_Bkg_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAU_BkgNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAU_BkgPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      } else {
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUCOND_Bkg_%s", (isPbp?"Pbp":"PP")),
                        Form("b_Bkg_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUCOND_BkgNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUCOND_BkgPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      }
      if ( ws.pdf(Form("pdfMASS_Bkg_%s", (isPbp?"Pbp":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_BkgPR_%s", (isPbp?"Pbp":"PP")),
                        Form("%s_BkgPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbp?"Pbp":"PP")),
                        Form("pdfMASS_Bkg_%s",(isPbp?"Pbp":"PP"))
                        ));
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_BkgNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("%s_BkgNoPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbp?"Pbp":"PP")),
                        Form("pdfMASS_Bkg_%s",(isPbp?"Pbp":"PP"))
                        ));
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUMASS_Bkg_%s", (isPbp?"Pbp":"PP")),
                        Form("b_Bkg_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUMASS_BkgNoPR_%s", (isPbp?"Pbp":"PP")),
                        Form("pdfCTAUMASS_BkgPR_%s", (isPbp?"Pbp":"PP"))
                        ));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Bkg_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("%s_Bkg_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")),
                      Form("N_Bkg_%s", (isPbp?"Pbp":"PP"))
                      ));
    }
  }
  

  // Total PDF
  RooAbsPdf *themodel = NULL;
  string tag = "";
  if (incPrompt  && !incNonPrompt) { tag = "PR"; }
  if (!incPrompt && incNonPrompt)  { tag = "NoPR";   }

  if (incJpsi && incPsi2S && incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Jpsi%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP"))), 
                                        *ws.pdf(Form("%s_Psi2S%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP"))), 
                                        *ws.pdf(Form("%s_Bkg%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP"))) 
                                         ),
                             RooArgList( 
                                        *ws.var(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))), 
                                        (ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))?*ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))):*ws.function(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))),
                                        *ws.var(Form("N_Bkg_%s", (isPbp?"Pbp":"PP"))) 
                                         )
                             );
  }
  if (incJpsi && incPsi2S && !incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Jpsi%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP"))), 
                                        *ws.pdf(Form("%s_Psi2S%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP")))
                                         ),
                             RooArgList( 
                                        *ws.var(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))),
                                        (ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))?*ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))):*ws.function(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))))
                                         )
                             );
  }
  if (incJpsi && !incPsi2S && incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Jpsi%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP"))), 
                                        *ws.pdf(Form("%s_Bkg%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP")))
                                         ),
                             RooArgList( 
                                        *ws.var(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))), 
                                        *ws.var(Form("N_Bkg_%s", (isPbp?"Pbp":"PP")))
                                         )
                             );
  }
  if (!incJpsi && incPsi2S && incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Psi2S%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP"))), 
                                        *ws.pdf(Form("%s_Bkg%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP")))
                                         ),
                             RooArgList( 
                                        (ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))?*ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))):*ws.function(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))),
                                        *ws.var(Form("N_Bkg_%s", (isPbp?"Pbp":"PP")))
                                         )
                             );
  }
  if (incJpsi && !incPsi2S && !incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Jpsi%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP")))
                                         ),
                             RooArgList( 
                                        *ws.var(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")))
                                         )
                             );
  }
  if (!incJpsi && incPsi2S && !incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Psi2S%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP")))
                                         ),
                             RooArgList( 
                                        (ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))?*ws.var(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))):*ws.function(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))))
                                         )
                             );
  }
  if (!incJpsi && !incPsi2S && incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Bkg%s_%s", pdfName.c_str(), tag.c_str(), (isPbp?"Pbp":"PP")))
                                         ),
                             RooArgList( 
                                        *ws.var(Form("N_Bkg_%s", (isPbp?"Pbp":"PP")))
                                         )
                             );
  }
  if (!incJpsi && !incPsi2S && !incBkg) {
    cout << "[ERROR] User did not include any model, please fix your input settings!" << endl; return false;
  }
  ws.import(*themodel);
  //ws.pdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbp?"Pbp":"PP")))->setNormRange("CtauWindow");
  
  return true;
};


void fixCtauParPsi2StoJpsi(map<string, string>& parIni, bool isPbp)
{
  cout << "[INFO] Constraining Psi(2S) parameters to Jpsi" << endl;
  parIni[Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0',{%s})", Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP") ));
  parIni[Form("lambdaDSS_Psi2SNoPR_Pbp")] = Form("RooFormulaVar::%s('@0',{%s})", Form("lambdaDSS_Psi2SNoPR_Pbp"), Form("lambdaDSS_JpsiNoPR_Pbp"));
  //parIni[Form("b_Psi2S_Pbp")] = Form("RooFormulaVar::%s('@0',{%s})", Form("b_Psi2S_Pbp"), Form("b_Jpsi_Pbp"));
};


bool defineCtauResolModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp)
{ 
  if (ws.pdf(Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))) { 
    cout << Form("[ERROR] The %s Ctau Resolution Model has already been implemented!", object.c_str()) << endl;
    return false; 
  }

  string pdfType = "pdfCTAURES";
  string varName = "ctau";
  bool usePerEventError = true;

  cout << Form("[INFO] Implementing %s Ctau Resolution Model", object.c_str()) << endl;
  
  switch(model) 
    {
    case (CtauModel::SingleGaussianResolution):  
      if (!( 
            parIni.count(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP"))) 
             )) {
 	cout << Form("[ERROR] Initial parameters where not found for Single Gaussian Ctau Resolution Model in %s", (isPbp?"Pbp":"PP")) << endl; return false; 
      }

      // create the variables for this model  
      if (!ws.var("One")) { ws.factory("One[1.0]"); }
      if (!ws.var(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))))  { ws.factory( parIni[Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() );  }
      if (!ws.var(Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP")))) { ws.factory( parIni[Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() ); }

      // create the PDF
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), varName.c_str(),
 		      Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP")),
 		      Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));
      
      cout << Form("[INFO] %s Single Gaussian Ctau Resolution PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;

    case (CtauModel::DoubleGaussianResolution):  

      if (!( 
            parIni.count(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))) 
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for Double Gaussian Ctau Resolution Model in %s", (isPbp?"Pbp":"PP")) << endl; return false; 
      }
      
      // create the variables for this model  
      if (!ws.var("One")) { ws.factory("One[1.0]"); }
      if (!ws.var(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))))   { ws.factory( parIni[Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str()  ); }
      if (!ws.var(Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP"))))  { ws.factory( parIni[Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() ); }
      if (!ws.var(Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))))   { ws.factory( parIni[Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str()  ); }
      if (!ws.var(Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP"))))  { ws.factory( parIni[Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() ); }
      if (!ws.var(Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))))       { ws.factory( parIni[Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() );      }
      
      // create the two PDFs
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s1_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), varName.c_str(),
 		      Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP")),
 		      Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s2_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), varName.c_str(),
 		      Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP")), 
 		      Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));

      // combine the two PDFs
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%s_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("%s1_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("%s2_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")),  
 		      Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))
 		      ));

      cout << Form("[INFO] %s Double Gaussian Ctau Resolution PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;


    case (CtauModel::TripleGaussianResolution):  

      if (!( 
            parIni.count(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("f2_CtauRes_%s", (isPbp?"Pbp":"PP")))
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for Triple Gaussian Ctau Resolution Model in %s", (isPbp?"Pbp":"PP")) << endl; return false; 
      }
      
      // create the variables for this model  
      if (!ws.var("One")) { ws.factory("One[1.0]"); }
      if (!ws.var(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))))   { ws.factory( parIni[Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str()  ); }
      if (!ws.var(Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP"))))  { ws.factory( parIni[Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() ); }
      if (!ws.var(Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))))   { ws.factory( parIni[Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str()  ); }
      if (!ws.var(Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP"))))  { ws.factory( parIni[Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() ); }
      if (!ws.var(Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP"))))   { ws.factory( parIni[Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str()  ); }
      if (!ws.var(Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP"))))  { ws.factory( parIni[Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() ); }
      if (!ws.var(Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))))       { ws.factory( parIni[Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() );      }
      if (!ws.var(Form("f2_CtauRes_%s", (isPbp?"Pbp":"PP"))))      { ws.factory( parIni[Form("f2_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() );     }
      
      // create the three PDFs
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s1_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), varName.c_str(),
 		      Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP")),
 		      Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s2_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), varName.c_str(),
 		      Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP")), 
 		      Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s3_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), varName.c_str(),
 		      Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP")), 
 		      Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));

      // combine the two PDFs
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%s23_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("%s2_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("%s3_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")),  
 		      Form("f2_CtauRes_%s", (isPbp?"Pbp":"PP"))
 		      ));
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%s_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("%s1_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("%s23_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")),  
 		      Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))
 		      ));

      cout << Form("[INFO] %s Triple Gaussian Ctau Resolution PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;

    case (CtauModel::QuadrupleGaussianResolution):  

      if (!( 
            parIni.count(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("ctau4_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("s4_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("f2_CtauRes_%s", (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("f3_CtauRes_%s", (isPbp?"Pbp":"PP")))
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for Triple Gaussian Ctau Resolution Model in %s", (isPbp?"Pbp":"PP")) << endl; return false; 
      }
      
      // create the variables for this model  
      if (!ws.var("One")) { ws.factory("One[1.0]"); }
      if (!ws.var(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))))   { ws.factory( parIni[Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str()  ); }
      if (!ws.var(Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP"))))  { ws.factory( parIni[Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() ); }
      if (!ws.var(Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))))   { ws.factory( parIni[Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str()  ); }
      if (!ws.var(Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP"))))  { ws.factory( parIni[Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() ); }
      if (!ws.var(Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP"))))   { ws.factory( parIni[Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str()  ); }
      if (!ws.var(Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP"))))  { ws.factory( parIni[Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() ); }
      if (!ws.var(Form("ctau4_CtauRes_%s", (isPbp?"Pbp":"PP"))))   { ws.factory( parIni[Form("ctau4_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str()  ); }
      if (!ws.var(Form("s4_CtauRes_%s", (isPbp?"Pbp":"PP"))))  { ws.factory( parIni[Form("s4_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() ); }
      if (!ws.var(Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))))       { ws.factory( parIni[Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() );      }
      if (!ws.var(Form("f2_CtauRes_%s", (isPbp?"Pbp":"PP"))))      { ws.factory( parIni[Form("f2_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() );     }
      if (!ws.var(Form("f3_CtauRes_%s", (isPbp?"Pbp":"PP"))))      { ws.factory( parIni[Form("f3_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str() );     }
      
      // create the four PDFs
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s1_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), varName.c_str(),
 		      Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP")),
 		      Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s2_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), varName.c_str(),
 		      Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP")), 
 		      Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s3_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), varName.c_str(),
 		      Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP")), 
 		      Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s4_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), varName.c_str(),
 		      Form("ctau4_CtauRes_%s", (isPbp?"Pbp":"PP")), 
 		      Form("s4_CtauRes_%s", (isPbp?"Pbp":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));

      // combine the two PDFs
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%s34_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("%s3_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("%s4_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")),  
 		      Form("f3_CtauRes_%s", (isPbp?"Pbp":"PP"))
 		      ));
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%s234_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("%s2_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("%s34_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")),  
 		      Form("f2_CtauRes_%s", (isPbp?"Pbp":"PP"))
 		      ));
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%s_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("%s1_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("%s234_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")),  
 		      Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))
 		      ));

      cout << Form("[INFO] %s Quadruple Gaussian Ctau Resolution PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;
  
    default :

      cout<< "[ERROR] Selected Ctau Resolution Model has not been implemented"<< endl; return false;

    }
 
  return true;

};


bool addBackgroundCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp) 
{
  if (ws.pdf(Form("pdfCTAUTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))) { 
    cout << Form("[ERROR] The %s Background Ctau Model has already been implemented!", object.c_str()) << endl;
    return false; 
  }

  cout << Form("[INFO] Implementing %s Background Ctau Model", object.c_str()) << endl;
   
  string objectInc = object; 
  if (objectInc.find("NoPR")!=std::string::npos) { objectInc.erase(objectInc.find("NoPR"), objectInc.length()); }
  if (objectInc.find("PR")!=std::string::npos)   { objectInc.erase(objectInc.find("PR"), objectInc.length());   }

  switch(model) 
    {  
    case (CtauModel::TripleDecay): 
      if (!( 
            parIni.count(Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambdaDF_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambdaDDS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("fDFSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("fDLIV_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for %s Background Triple Decay Ctau Model in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false; 
      }
      
      // create the variables for this model 
      if ( !ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))) ){ 
        ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
      } 
      ws.factory( parIni[Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
      ws.factory( parIni[Form("lambdaDF_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()  );
      ws.factory( parIni[Form("lambdaDDS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
      ws.factory( parIni[Form("fDFSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()     ); 
      ws.factory( parIni[Form("fDLIV_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()     );

      // create the three PDFs
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", Form("pdfCTAUDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctau", 
 		      Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", Form("pdfCTAUDF_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctau", 
 		      Form("lambdaDF_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", Form("pdfCTAUDDS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctau", 
 		      Form("lambdaDDS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
      
      // combine the three PDFs
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("fDFSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAUDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAUDF_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("fDLIV_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAU1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAUDDS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
 
      cout << Form("[INFO] %s Background Triple Decay Ctau PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
 
    case (CtauModel::Delta): 
      
      if ( !ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))) ){ 
        ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
      }

      // create the PDF
      ws.factory(Form("SUM::%s(%s)", Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));

      cout << Form("[INFO] %s Background Delta Ctau PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
      
    default :

      cout<< "[ERROR] Selected Background Ctau Model has not been implemented"<< endl; return false;

    }
   
  return true;
};


bool addSignalCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp) 
{
  if (ws.pdf(Form("pdfCTAUTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))) { 
    cout << Form("[ERROR] The %s Signal Ctau Model has already been implemented!", object.c_str()) << endl;
    return false; 
  }

  cout << Form("[INFO] Implementing %s Signal Ctau Model", object.c_str()) << endl;
   
  string objectInc = object; 
  if (objectInc.find("NoPR")!=std::string::npos) { objectInc.erase(objectInc.find("NoPR"), objectInc.length()); }
  if (objectInc.find("PR")!=std::string::npos)   { objectInc.erase(objectInc.find("PR"), objectInc.length());   } 

  switch(model) 
    {  
    case (CtauModel::SingleSidedDecay): 

      if (!( 
            parIni.count(Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for %s Signal Single Sided Decay Ctau Model in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false; 
      }
      
      // create the variables for this model 
      if ( !ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))) ){ 
        ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
      }
      ws.factory( parIni[Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
            
      // create the PDF
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctau", 
 		      Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
 
      cout << Form("[INFO] %s Signal Single Sided Decay Ctau PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
 
    case (CtauModel::Delta):

      if ( !ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))) ){ 
        ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
      }

      // create the three PDFs 
      ws.factory(Form("SUM::%s(%s)", Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAURES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
 		      
      cout << Form("[INFO] %s Signal Delta Ctau PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
 
    default :

      cout<< "[ERROR] Selected Signal Ctau Model has not been implemented"<< endl; return false;

    }
   
  return true;
};
  

void setCtauDefaultParameters(map<string, string> &parIni, bool isPbp, double numEntries)
{

  cout << "[INFO] Setting user undefined initial parameters to their default values" << endl;

  // DEFAULT RANGE OF NUMBER OF EVENTS
  if (parIni.count(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count(Form("N_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("N_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("N_Bkg_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Bkg_%s", (isPbp?"Pbp":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count(Form("b_Jpsi_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("b_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("b_Jpsi_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("b_Jpsi_%s", (isPbp?"Pbp":"PP")), 0.2, 0.0, 1.0);
  }
  if (parIni.count(Form("b_Psi2S_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("b_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("b_Psi2S_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("b_Psi2S_%s", (isPbp?"Pbp":"PP")), 0.2, 0.0, 1.0);
  }
  if (parIni.count(Form("b_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("b_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("b_Bkg_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("b_Bkg_%s", (isPbp?"Pbp":"PP")), 0.2, 0.0, 1.0);
  }

 // CTAU FIT PARAMETERS

 // Resolution Ctau Model
  if (parIni.count(Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.0, -1.0, 1.0);
    
  }
  if (parIni.count(Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP")))==0) {
    parIni[Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0',{%s})", Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP")), Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP") ));
  } else if (parIni[Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.0, -1.0, 1.0);
  }
  if (parIni.count(Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP")))==0) {
    parIni[Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0',{%s})", Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP")), Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP") ));
  } else if (parIni[Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.0, -1.0, 1.0);
  }
  if (parIni.count(Form("ctau4_CtauRes_%s", (isPbp?"Pbp":"PP")))==0) {
    parIni[Form("ctau4_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0',{%s})", Form("ctau4_CtauRes_%s", (isPbp?"Pbp":"PP")), Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP") ));
  } else if (parIni[Form("ctau4_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("ctau4_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("ctau4_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.0, -1.0, 1.0);
  }
  if (parIni.count(Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP")), 1.0, 0.01, 2.0);
  }
  if (parIni.count(Form("rS21_CtauRes_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("rS21_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") {
      parIni[Form("rS21_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("rS21_CtauRes_%s", (isPbp?"Pbp":"PP")), 1.2, 1.0, 3.0);
  }
  if (parIni.count(Form("rS32_CtauRes_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("rS32_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") {
      parIni[Form("rS32_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("rS32_CtauRes_%s", (isPbp?"Pbp":"PP")), 1.2, 1.0, 3.0);
  }
  if (parIni.count(Form("rS43_CtauRes_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("rS43_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") {
      parIni[Form("rS43_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("rS43_CtauRes_%s", (isPbp?"Pbp":"PP")), 1.2, 1.0, 3.0);
  }
  if (parIni.count(Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP")))==0) {
    parIni[Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP")),
                                                                   parIni[Form("rS21_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str(), Form("s1_CtauRes_%s", (isPbp?"Pbp":"PP") ));
  } else if ( parIni[Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP")), 1.2, 0.001, 6.0);
  }
  if (parIni.count(Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP")))==0) {
    parIni[Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP")),
                                                                   parIni[Form("rS32_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str(), Form("s2_CtauRes_%s", (isPbp?"Pbp":"PP") ));
  } else if ( parIni[Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP")), 1.2, 0.001, 6.0);
  }
  if (parIni.count(Form("s4_CtauRes_%s", (isPbp?"Pbp":"PP")))==0) {
    parIni[Form("s4_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("s4_CtauRes_%s", (isPbp?"Pbp":"PP")),
                                                                   parIni[Form("rS43_CtauRes_%s", (isPbp?"Pbp":"PP"))].c_str(), Form("s3_CtauRes_%s", (isPbp?"Pbp":"PP") ));
  } else if ( parIni[Form("s4_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("s4_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("s4_CtauRes_%s", (isPbp?"Pbp":"PP")), 1.2, 0.001, 6.0);
  }
  if (parIni.count(Form("f_CtauRes_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("f_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("f_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.8, 0.4, 1.0);
  }
  if (parIni.count(Form("f2_CtauRes_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("f2_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("f2_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("f2_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.8, 0.4, 1.0);
  }
  if (parIni.count(Form("f3_CtauRes_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("f3_CtauRes_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("f3_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("f3_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.8, 0.4, 1.0);
  }

  // Signal Ctau Model
  if (parIni.count(Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")), 0.8, 0.01, 2.0);
  }
  if (parIni.count(Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), 0.8, 0.01, 2.0);
  }

  // Background Ctau Model
  if (parIni.count(Form("f_BkgNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("f_BkgNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("f_BkgNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("f_BkgNoPR_%s", (isPbp?"Pbp":"PP")), 0.3, 0., 1.);
  }
  if (parIni.count(Form("fDFSS_BkgNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("fDFSS_BkgNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("fDFSS_BkgNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDFSS_BkgNoPR_%s", (isPbp?"Pbp":"PP")), 0.8, 0., 1.);
  }
  if (parIni.count(Form("fDLIV_BkgNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("fDLIV_BkgNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("fDLIV_BkgNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDLIV_BkgNoPR_%s", (isPbp?"Pbp":"PP")), 0.5, 0., 1.);
  }
  if (parIni.count(Form("lambdaDSS_BkgNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDSS_BkgNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDSS_BkgNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDSS_BkgNoPR_%s", (isPbp?"Pbp":"PP")), 0.45, 0.0001, 5.0);
  }
  if (parIni.count(Form("lambdaDF_BkgNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDF_BkgNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDF_BkgNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDF_BkgNoPR_%s", (isPbp?"Pbp":"PP")), 0.30, 0.0001, 5.0);
  }
  if (parIni.count(Form("lambdaDDS_BkgNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDDS_BkgNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDDS_BkgNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDDS_BkgNoPR_%s", (isPbp?"Pbp":"PP")), 0.06, 0.0001, 5.0);
  }
  if (parIni.count(Form("f_BkgPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("f_BkgPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("f_BkgPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("f_BkgPR_%s", (isPbp?"Pbp":"PP")), 0.3, 0., 1.);
  }
  if (parIni.count(Form("fDFSS_BkgPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("fDFSS_BkgPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("fDFSS_BkgPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDFSS_BkgPR_%s", (isPbp?"Pbp":"PP")), 0.9, 0., 1.);
  }
  if (parIni.count(Form("fDLIV_BkgPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("fDLIV_BkgPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("fDLIV_BkgPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDLIV_BkgPR_%s", (isPbp?"Pbp":"PP")), 0.9, 0., 1.);
  }
  if (parIni.count(Form("lambdaDSS_BkgPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDSS_BkgPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDSS_BkgPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDSS_BkgPR_%s", (isPbp?"Pbp":"PP")), 0.42, 0.05, 1.5);
  }
  if (parIni.count(Form("lambdaDF_BkgPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDF_BkgPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDF_BkgPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDF_BkgPR_%s", (isPbp?"Pbp":"PP")), 0.79, 0.001, 1.5);
  }
  if (parIni.count(Form("lambdaDDS_BkgPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDDS_BkgPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDDS_BkgPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDDS_BkgPR_%s", (isPbp?"Pbp":"PP")), 0.69, 0.02, 5.0);
  }

};


#endif // #ifndef buildCharmoniaCtauModel_C
