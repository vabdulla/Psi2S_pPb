#ifndef buildCharmoniaCtauTrueModel_C
#define buildCharmoniaCtauTrueModel_C

#include "Utilities/initClasses.h"

void setCtauTrueDefaultParameters(map<string, string> &parIni, bool isPbp, double numEntries);
bool defineCtauTrueResolModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp); 
bool addSignalCtauTrueModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp);


bool buildCharmoniaCtauTrueModel(RooWorkspace& ws, struct CharmModel model, map<string, string>  parIni, 
                                 bool isPbp,                 // Determine if we are working with Pbp (True) or PP (False)
                                 bool incJpsi,                // Include Jpsi model
                                 bool incPsi2S,                // Include Psi(2S) model
                                 double  numEntries = 300000. // Number of entries in the dataset
                                 )
{

  // If the initial parameters are empty, set defaul parameter values
  setCtauTrueDefaultParameters(parIni, isPbp, numEntries);

  // C r e a t e   m o d e l 
     
  if (incJpsi) {
    if(!defineCtauTrueResolModel(ws, "JpsiNoPR", model.CtauTrueRes, parIni, isPbp)) { cout << "[ERROR] Defining the Ctau Truth Resolution Model failed" << endl; return false; }
    if(!addSignalCtauTrueModel(ws, "JpsiNoPR", model.CtauTrue, parIni, isPbp)) { cout << "[ERROR] Adding Non-Prompt Jpsi Ctau Truth Model failed" << endl; return false; }
  }
  if (incPsi2S) {
    if(!defineCtauTrueResolModel(ws, "Psi2SNoPR", model.CtauTrueRes, parIni, isPbp)) { cout << "[ERROR] Defining the Ctau Truth Resolution Model failed" << endl; return false; }
    if(!addSignalCtauTrueModel(ws, "Psi2SNoPR", model.CtauTrue, parIni, isPbp)) { cout << "[ERROR] Adding Non-Prompt Psi2S Ctau Truth Model failed" << endl; return false; }
  }

  // Total PDF
  string pdfType = "pdfCTAUTRUE";
  string pdfName = Form("%s_Tot_%s", pdfType.c_str(), (isPbp?"Pbp":"PP"));

  RooArgList pdfList;
  if (incJpsi) { pdfList.add( *ws.pdf(Form("%sTot_JpsiNoPR_%s", pdfType.c_str(), (isPbp?"Pbp":"PP"))) );  }
  if (incPsi2S){ pdfList.add( *ws.pdf(Form("%sTot_Psi2SNoPR_%s", pdfType.c_str(), (isPbp?"Pbp":"PP"))) ); }
  if (incJpsi && incPsi2S) {
    cout << "[ERROR] User included Jpsi and Psi2S together, but can only fit them separately, please fix your input settings!" << endl; return false;
  }
  if (!incJpsi && !incPsi2S) {
    cout << "[ERROR] User did not include any model, please fix your input settings!" << endl; return false;
  }
  RooAbsPdf *themodel = new RooAddPdf(pdfName.c_str(), pdfName.c_str(), pdfList);
  ws.import(*themodel);
  ws.pdf(pdfName.c_str())->setNormRange("CtauTrueWindow");
  
  setFixedVarsToContantVars(ws);

  // save the initial values of the model we've just created
  RooArgSet* params = (RooArgSet*) themodel->getParameters(RooArgSet(*ws.var("ctauTrue")));
  ws.saveSnapshot((pdfName+"_parIni").c_str(),*params,kTRUE);
  
  return true;
};


bool defineCtauTrueResolModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp) 
{ 

  if (ws.pdf(Form("pdfCTAUTRUERES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))) { 
    cout << Form("[ERROR] The %s Ctau Truth Resolution Model has already been implemented!", object.c_str()) << endl;
    return false; 
  }

  cout << Form("[INFO] Implementing %s Ctau Truth Resolution Model", object.c_str()) << endl;
  
  switch(model) 
    {  
    case (CtauModel::SingleGaussianResolution):
      if (!( 
            parIni.count(Form("ctauMC_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("sigmaMC_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for Single Gaussian Ctau Truth Resolution Model in %s", (isPbp?"Pbp":"PP")) << endl; return false; 
      }

      // create the variables for this model  
      if (!ws.var("One")) { ws.factory("One[1.0]"); }
      if (!ws.var(Form("ctauMC_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))))  { ws.factory( parIni[Form("ctauMC_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );  }
      if (!ws.var(Form("sigmaMC_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))) { ws.factory( parIni[Form("sigmaMC_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); }

      // create the PDFs
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, One)", Form("pdfCTAUTRUERES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctauTrue", 
 		      Form("ctauMC_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
 		      Form("sigmaMC_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
      
      cout << Form("[INFO] %s Single Gaussian Ctau Truth Resolution PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;
    
    case (CtauModel::Delta):

      // create the Delta PDF
      ws.factory(Form("TruthModel::%s(%s)", Form("pdfCTAUTRUERES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctauTrue"));

      cout << Form("[INFO] %s Delta Ctau Truth Resolution PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;
 
    default :

      cout<< "[ERROR] Selected Ctau Truth Resolution Model has not been implemented"<< endl; return false;

    }
 
  return true;

};


bool addSignalCtauTrueModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbp) 
{
  if (ws.pdf(Form("pdfCTAUTRUE_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))) { 
    cout << Form("[ERROR] The %s Signal Ctau Truth Model has already been implemented!", object.c_str()) << endl;
    return false; 
  }

  cout << Form("[INFO] Implementing %s Signal Ctau Truth Model", object.c_str()) << endl;
   
  string objectInc = object; 
  if (objectInc.find("NoPR")!=std::string::npos) { objectInc.erase(objectInc.find("NoPR"), objectInc.length()); }
  if (objectInc.find("PR")!=std::string::npos)   { objectInc.erase(objectInc.find("PR"), objectInc.length());   } 

  switch(model) 
    {
    case (CtauModel::DoubleSingleSidedDecay): 
      if (!( 
            parIni.count(Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))  &&
	    parIni.count(Form("lambdaDSS2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("fDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))
             )) {
	cout << Form("[ERROR] Initial parameters where not found for %s Signal Double Single Sided Decay Ctau Truth Model in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false; 
      }
      
      // create the variables for this model 
      if ( !ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))) ){ ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); }
      ws.factory( parIni[Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()  );
      ws.factory( parIni[Form("lambdaDSS2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
      ws.factory( parIni[Form("fDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()       ); 

      // create the two PDFs
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", Form("pdfCTAUTRUEDSS1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctauTrue",
		      Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAUTRUERES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));

      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", Form("pdfCTAUTRUEDSS2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctauTrue",
 		      Form("lambdaDSS2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAUTRUERES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
      
      // combine the two PDFs
      ws.factory(Form("AddModel::%s({%s , %s}, %s)", Form("pdfCTAUTRUE_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
		      Form("pdfCTAUTRUEDSS1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAUTRUEDSS2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("fDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));
      
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfCTAUTRUETot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
		      Form("pdfCTAUTRUE_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
		      Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      cout << Form("[INFO] %s Signal Double Single Sided Decay Ctau Truth PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 

    case (CtauModel::SingleSidedDecay): 

      if (!( 
            parIni.count(Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for %s Signal Single Sided Decay Ctau Truth Model in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false; 
      }
      
      // create the variables for this model 
      if (!ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP")))){ ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); } 
      ws.factory( parIni[Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
            
      // create the PDF
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", Form("pdfCTAUTRUE_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "ctauTrue", 
 		      Form("lambdaDSS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAUTRUERES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
      
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfCTAUTRUETot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfCTAUTRUE_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))
                      ));
 
      cout << Form("[INFO] %s Signal Single Sided Decay Ctau Truth PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
 
    case (CtauModel::Delta):

      if (!ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP")))){ ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); } 

      // create the three PDFs 
      ws.factory(Form("AddModel::%s(%s)", Form("pdfCTAUTRUE_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
 		      Form("pdfCTAUTRUERES_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
 		      ));
     
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfCTAUTRUETot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfCTAUTRUE_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("N_%s_%s", objectInc.c_str(), (isPbp?"Pbp":"PP"))
                      ));
 		      
      cout << Form("[INFO] %s Signal Delta Ctau Truth PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
 
    default :

      cout<< "[ERROR] Selected Signal Ctau Truth Model has not been implemented"<< endl; return false;

    }
   
  return true;
};


void setCtauTrueDefaultParameters(map<string, string> &parIni, bool isPbp, double numEntries)
{
  
  cout << "[INFO] Setting user undefined initial parameters to their default values" << endl;
  
  // DEFAULT RANGE OF NUMBER OF EVENTS
  if (parIni.count(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  
  // CTAU TRUE FIT PARAMETERS
  
 // Resolution Ctau Truth Model
  if (parIni.count(Form("ctauMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("ctauMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("ctauMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("ctauMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP")), 0.0, 0.0, 0.0);
  }
  if (parIni.count(Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigmaMC_JpsiNoPR_%s", (isPbp?"Pbp":"PP")), 0.001, 0.000000001, 1.0);
  }
  if (parIni.count(Form("ctauMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("ctauMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("ctauMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("ctauMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), 0.0, 0.0, 0.0);
  }
  if (parIni.count(Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigmaMC_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), 0.001, 0.000000001, 1.0);
  }

  // Signal Ctau Truth Model
  if (parIni.count(Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")), 0.8, 0.001, 2.0);
  }
  if (parIni.count(Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), 0.8, 0.001, 2.0);
  }
  if (parIni.count(Form("lambdaDSS2_JpsiNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDSS2_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDSS2_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDSS2_JpsiNoPR_%s", (isPbp?"Pbp":"PP")), 0.8, 0.001, 2.0);
  }
  if (parIni.count(Form("lambdaDSS2_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambdaDSS2_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("lambdaDSS2_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDSS2_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), 0.8, 0.001, 2.0);
  }
  if (parIni.count(Form("fDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("fDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("fDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDSS_JpsiNoPR_%s", (isPbp?"Pbp":"PP")), 0.8, 0., 1.);
  }
  if (parIni.count(Form("fDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("fDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("fDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDSS_Psi2SNoPR_%s", (isPbp?"Pbp":"PP")), 0.8, 0., 1.);
  }
  
};


#endif // #ifndef buildCharmoniaCtauTrueModel_C
