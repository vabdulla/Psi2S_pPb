#ifndef buildCharmoniaMassModel_C
#define buildCharmoniaMassModel_C

#include "Utilities/initClasses.h"

void fixMassParPsi2StoJpsi(map<string, string>& parIni, bool isPbp);
void fixPbptoPP(map<string, string>& parIni);
void setMassDefaultParameters(map<string, string> &parIni, bool isPbp, double numEntries);
bool addSignalMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni, bool isPbp); 
bool addBackgroundMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni, bool isPbp);


bool buildCharmoniaMassModel(RooWorkspace& ws, struct CharmModel model, map<string, string>  parIni, 
                             bool isPbp,                 // Determine if we are working with Pbp (True) or PP (False)
                             bool doSimulFit,             // Do simultaneous fit
                             bool incBkg,                 // Include background model
                             bool incJpsi,                // Include Jpsi model
                             bool incPsi2S,               // Include Psi(2S) model
                             double  numEntries = 300000. // Number of entries in the dataset
                             )
{

  // If the initial parameters are empty, set defaul parameter values
  setMassDefaultParameters(parIni, isPbp, numEntries);

  // Fix all psi2S parameters to jpsi
  if (incJpsi && incPsi2S) {
    fixMassParPsi2StoJpsi(parIni, isPbp);
  }

  // Let's define the single and double ratio variables
  if (incPsi2S && incJpsi) {
    if (doSimulFit && isPbp) {

      // Fix mean, alpha and n parameters in Pbp to PP values 
      // fixPbptoPP(parIni);

      // create parameters related to the double ratio
      ws.factory( parIni["RFrac2Svs1S_PbpvsPP"].c_str() );                     // Double Ratio
      ws.factory( parIni[Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))].c_str() );    // Number of Jpsi in Pbp

      // Asign NPsi2S(Pbp) = DoubleRatio(Pbp/PP) * SingleRatioPP(Psi2S/Jpsi) * NJpsi(Pbp) 
      ws.factory(Form("RooFormulaVar::%s('@0*@1*@2',{%s,%s,%s})", "N_Psi2S_Pbp", 
		      "RFrac2Svs1S_PbpvsPP", 
		      "RFrac2Svs1S_PP", 
		      Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))
		      )); 
      ws.var("RFrac2Svs1S_PbpvsPP")->setConstant(kFALSE);

    } else {

      // create parameters related to the double ratio
      ws.factory( parIni[Form("RFrac2Svs1S_%s", (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))].c_str() );

      // Asign N(Psi2S) = SingleRatio(Psi2S/Jpsi) * N(Jpsi)
      ws.factory(Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")), 
		      Form("RFrac2Svs1S_%s", (isPbp?"Pbp":"PP")),
		      Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")) 
		      )); 
      ws.var(Form("RFrac2Svs1S_%s", (isPbp?"Pbp":"PP")))->setConstant(kFALSE);

    }
    
    // As we have already declared the Number of Psi2S and Jpsi before, just use their names in their PDFs
    parIni[Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))] = Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"));
    parIni[Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")); 
  }

  // C r e a t e   m o d e l  

  if (incJpsi) {
    if(!addSignalMassModel(ws, "Jpsi", model.Jpsi.Mass, parIni, isPbp)) { cout << "[ERROR] Adding Jpsi Mass Model failed" << endl; return false; }
  }
  if (incPsi2S) { 
    if (!addSignalMassModel(ws, "Psi2S", model.Psi2S.Mass, parIni, isPbp)) { cout << "[ERROR] Adding Psi(2S) Mass Model failed" << endl; return false; }  
  }
  if (incBkg) {
    if(!addBackgroundMassModel(ws, "Bkg", model.Bkg.Mass, parIni, isPbp)) { cout << "[ERROR] Adding Background Mass Model failed" << endl; return false; }
  }

  // Total PDF
  string pdfType = "pdfMASS";
  string pdfName = Form("%s_Tot_%s", pdfType.c_str(), (isPbp?"Pbp":"PP"));

  RooArgList pdfList;
  if (incJpsi) { pdfList.add( *ws.pdf(Form("%sTot_Jpsi_%s", pdfType.c_str(), (isPbp?"Pbp":"PP"))) );  }
  if (incPsi2S){ pdfList.add( *ws.pdf(Form("%sTot_Psi2S_%s", pdfType.c_str(), (isPbp?"Pbp":"PP"))) ); }
  if (incBkg)  { pdfList.add( *ws.pdf(Form("%sTot_Bkg_%s", pdfType.c_str(), (isPbp?"Pbp":"PP"))) );   }
  if (!incJpsi && !incPsi2S && !incBkg) { cout << "[ERROR] User did not include any model, please fix your input settings!" << endl; return false; }
  RooAbsPdf *themodel = new RooAddPdf(pdfName.c_str(), pdfName.c_str(), pdfList );
  ws.import(*themodel);
  ws.pdf(pdfName.c_str())->setNormRange("MassWindow");

  setFixedVarsToContantVars(ws);

  // save the initial values of the model we've just created
  RooArgSet* params = (RooArgSet*) themodel->getParameters(RooArgSet(*ws.var("invMass")));
  ws.saveSnapshot((pdfName+"_parIni").c_str(),*params,kTRUE);
  
  return true;
};

bool addBackgroundMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni, bool isPbp) 
{
  cout << Form("[INFO] Implementing %s Background Mass Model", object.c_str()) << endl;
  
  switch(model) 
    {  
    case (MassModel::Uniform): 

      // create the PDF           
      ws.factory(Form("Uniform::%s(%s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass"));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Background Uniform PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;
 
    case (MassModel::Chebychev1): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background First Order Chebychev in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      } 

      // create the variables for this model
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF           
      ws.factory(Form("Chebychev::%s(%s, {%s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Background 1st Order Chebychev PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;
 
    case (MassModel::Chebychev2): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Second Order Chebychev in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      } 

      // create the variables for this model
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF           
      ws.factory(Form("Chebychev::%s(%s, {%s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Background 2nd Order Chebychev PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 

    case (MassModel::Chebychev3): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Third Order Chebychev in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      }

      // create the variables for this model 
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF                 
      ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Background 3rd Order Chebychev PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 

    case (MassModel::Chebychev4): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Fourth Order Chebychev in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF                 
      ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Background 4th Order Chebychev PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 

    case (MassModel::Chebychev5): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda5_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Fifth Order Chebychev in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda5_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF                 
      ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda5_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Background 5th Order Chebychev PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 

    case (MassModel::Chebychev6): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda5_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("lambda6_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Sixth Order Chebychev in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda5_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda6_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF                 
      ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s, %s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda5_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda6_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Background 6th Order Chebychev PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 

    case (MassModel::ExpChebychev1): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background First Order Exponential Chebychev in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni["invMassNorm"].c_str() );
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF                 
      ws.factory(Form("RooFormulaVar::%s('@1*(@0) + 1.0', {%s, %s})", 
                      Form("pdfMASSPol_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMassNorm", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));               

      ws.factory(Form("Exponential::%s(%s, One[1.0])", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASSPol_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
                      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Background 1st Order Exponential Chebychev PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
   
    case (MassModel::ExpChebychev2): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))  
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Second Order Exponential Chebychev in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      }

      // create the variables for this model    
      ws.factory( parIni["invMassNorm"].c_str() );    
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF            
      ws.factory(Form("RooFormulaVar::%s('@2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s})", 
                      Form("pdfMASSPol_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMassNorm",
		      Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));             

      ws.factory(Form("Exponential::%s(%s, One[1.0])", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASSPol_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
                      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Background 2nd Order Exponential Chebychev PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
   
    case (MassModel::ExpChebychev3): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))  
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Third Order Exponential Chebychev in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      }

      // create the variables for this model  
      ws.factory( parIni["invMassNorm"].c_str() );      
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF
      ws.factory(Form("RooFormulaVar::%s('@3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s})", 
                      Form("pdfMASSPol_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMassNorm",
		      Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));        

      ws.factory(Form("Exponential::%s(%s, One[1.0])", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASSPol_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
                      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Background 3rd Order Exponential Chebychev PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
   
    case (MassModel::ExpChebychev4): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))    
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Fourth Order Exponential Chebychev in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      }

      // create the variables for this model   
      ws.factory( parIni["invMassNorm"].c_str() );     
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF                 
      ws.factory(Form("RooFormulaVar::%s('@4*(8.0*@0*@0*@0*@0 - 8.0*@0*@0 + 1.0) + @3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s, %s})", 
                      Form("pdfMASSPol_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMassNorm",
		      Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("Exponential::%s(%s, One[1.0])", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASSPol_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
                      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Background 4th Order Exponential Chebychev PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
   
    case (MassModel::ExpChebychev5): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda5_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))    
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Fifth Order Exponential Chebychev in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      }

      // create the variables for this model     
      ws.factory( parIni["invMassNorm"].c_str() );   
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda5_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF       
      ws.factory(Form("RooFormulaVar::%s('@5*(16.0*@0*@0*@0*@0*@0 - 20.0*@0*@0*@0 + 5.0*@0) + @4*(8.0*@0*@0*@0*@0 - 8.0*@0*@0 + 1.0) + @3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s, %s, %s})", 
                      Form("pdfMASSPol_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMassNorm", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda5_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("Exponential::%s(%s, One[1.0])", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASSPol_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
                      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Background 5th Order Exponential Chebychev PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 
   
    case (MassModel::ExpChebychev6): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda5_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))),
            parIni.count(Form("lambda6_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))     
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Sixth Order Exponential Chebychev in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      }

      // create the variables for this model    
      ws.factory( parIni["invMassNorm"].c_str() );    
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda5_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda6_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF
      ws.factory(Form("RooFormulaVar::%s('@6*(32.0*@0*@0*@0*@0*@0*@0 - 48.0*@0*@0*@0*@0 + 18.0*@0*@0 - 1.0) + @5*(16.0*@0*@0*@0*@0*@0 - 20.0*@0*@0*@0 + 5.0*@0) + @4*(8.0*@0*@0*@0*@0 - 8.0*@0*@0 + 1.0) + @3*(4.0*@0*@0*@0 - 3.0*@0) + @2*(2.0*@0*@0 - 1.0) + @1*(@0) + 1.0', {%s, %s, %s, %s, %s, %s, %s})", 
                      Form("pdfMASSPol_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMassNorm",
		      Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda3_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda4_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda5_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("lambda6_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("Exponential::%s(%s, One[1.0])", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASSPol_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
                      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Background 6th Order Exponential Chebychev PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 

    case (MassModel::Exponential): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Exponential in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF                 
      ws.factory(Form("Exponential::%s(%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Background Exponential PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;
   
    default :
      
      cout << "[ERROR] Selected Background Mass Model: " << parIni[Form("Model_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))] << " has not been implemented" << endl; return false;
    
    }
  
  return true;
};


bool addSignalMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni, bool isPbp) 
{
  cout << Form("[INFO] Implementing %s Mass Model", object.c_str()) << endl;
  
  switch(model) 
    {    
    case (MassModel::SingleGaussian): 
      
      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Single Gaussian Model in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF                       
      ws.factory(Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass", 
		      Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));
    
      cout << Form("[INFO] %s Single Gaussian PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;  
      
    case (MassModel::DoubleGaussian): 

      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("sigma2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("f_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Double Gaussian Model in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false; 
      }

      // create the variables for this model              
      ws.factory( parIni[Form("f_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("sigma2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the two PDFs             
      ws.factory(Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASS1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass", 
                      Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
                      Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
                      ));
      ws.factory(Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASS2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass", 
                      Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
                      Form("sigma2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
                      ));

      // Sum the PDFs to get the signal PDF
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
		      Form("f_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
		      Form("pdfMASS1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
		      Form("pdfMASS2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Double Gaussian PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break; 

    case (MassModel::SingleCrystalBall):  

      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("alpha_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("n_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))
             )) {
	cout << Form("[ERROR] Initial parameters where not found for %s Single Crystal Ball Model in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false; 
      }

      // create the variables for this model             
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("alpha_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("n_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the PDF              
      ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass", 
		      Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
		      Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
		      Form("alpha_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
		      Form("n_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Single Crystal Ball PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;
      
    case (MassModel::DoubleCrystalBall): 
      
      // check that all input parameters are defined
        if (!(
              parIni.count(Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) &&
              parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) &&
              parIni.count(Form("sigma2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) &&
              parIni.count(Form("alpha_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) &&
              parIni.count(Form("alpha2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) &&
              parIni.count(Form("n_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) &&
              parIni.count(Form("n2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) &&
              parIni.count(Form("f_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))
              )) {
          cout << Form("[ERROR] Initial parameters where not found for %s Double Crystal Ball Model in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
        }
        
      // create the variables for this model             
      ws.factory( parIni[Form("f_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("alpha_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("alpha2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("n_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("n2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("sigma2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );

      // create the two PDFs
      ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASS1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass", 
                      Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
                      Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("alpha_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("n_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
                      ));
      ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASS2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass", 
                      Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
                      Form("sigma2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("alpha2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("n2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
                      ));

      // Sum the PDFs to get the signal PDF
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
		      Form("f_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
		      Form("pdfMASS1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
		      Form("pdfMASS2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Double Crystal Ball PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;
      
    case (MassModel::GaussianAndCrystalBall):

      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("sigma2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) && 
            parIni.count(Form("alpha_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("n_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))) &&
            parIni.count(Form("f_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")))  
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Gaussian and Crystal Ball Model in %s", object.c_str(), (isPbp?"Pbp":"PP")) << endl; return false;
      } 

      // create the variables for this model             
      ws.factory( parIni[Form("f_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("sigma2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("alpha_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      ws.factory( parIni[Form("n_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str() );
      
      // create the two PDFs
      ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASS1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass",
                        Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                        Form("sigma1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                        Form("alpha_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                        Form("n_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
                        ));
        
      ws.factory(Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASS2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), "invMass",
                      Form("m_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")), 
                      Form("sigma2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
                      ));
 
      // Sum the PDFs to get the signal PDF 
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
		      Form("f_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
		      Form("pdfMASS1_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
		      Form("pdfMASS2_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))
		      ));

      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMASSTot_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      Form("pdfMASS_%s_%s", object.c_str(), (isPbp?"Pbp":"PP")),
                      parIni[Form("N_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))].c_str()
                      ));

      cout << Form("[INFO] %s Gaussian and Crystal Ball PDF in %s included", object.c_str(), (isPbp?"Pbp":"PP")) << endl; break;
	
    default :

      cout << "[ERROR] Selected Signal Mass Model: " << parIni[Form("Model_%s_%s", object.c_str(), (isPbp?"Pbp":"PP"))] << " has not been implemented" << endl; return false;

    }
  
  return true;
};

void fixPbptoPP(map<string, string>& parIni)
{
  //parIni["m_Jpsi_Pbp"]  = Form("RooFormulaVar::%s('@0',{%s})", "m_Jpsi_Pbp", "m_Jpsi_PP");
  //parIni["sigma1_Jpsi_Pbp"]  = Form("RooFormulaVar::%s('@0',{%s})", "sigma1_Jpsi_Pbp", "sigma1_JpsiPP");
  //parIni["sigma1_Psi2S_Pbp"] = Form("RooFormulaVar::%s('@0',{%s})", "sigma1_Psi2S_Pbp", "sigma1_Psi2SPP");
  // if (parIni.count("rSigma21_Jpsi_Pbp")!=0 && parIni.count("rSigma21_Jpsi_PP")!=0) {
  //   parIni["sigma2_Jpsi_Pbp"] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", "sigma2_Jpsi_Pbp", "sigma1_Jpsi_Pbp", "rSigma21_Jpsi_PP");
  // }
  // parIni["alpha_Jpsi_Pbp"]  = Form("RooFormulaVar::%s('@0',{%s})", "alpha_Jpsi_Pbp", "alpha_Jpsi_PP");
  // parIni["n_Jpsi_Pbp"] = Form("RooFormulaVar::%s('@0',{%s})", "n_Jpsi_Pbp", "n_Jpsi_PP");
  //parIni["f_Jpsi_Pbp"] = Form("RooFormulaVar::%s('@0',{%s})", "f_Jpsi_Pbp", "f_Jpsi_PP");
  //parIni["f_Psi2S_Pbp"] = Form("RooFormulaVar::%s('@0',{%s})", "f_Psi2S_Pbp", "f_Jpsi_PP");  
};

void fixMassParPsi2StoJpsi(map<string, string>& parIni, bool isPbp)
{
  cout << "[INFO] Constraining Psi(2S) parameters to Jpsi using PDF Mass Ratio" << endl;
  Double_t MassRatio = (Mass.Psi2S/Mass.JPsi);
  parIni[Form("m_Psi2S_%s", (isPbp?"Pbp":"PP"))]      = Form("RooFormulaVar::%s('@0*@1',{MassRatio[%.6f],%s})", Form("m_Psi2S_%s", (isPbp?"Pbp":"PP")), MassRatio, Form("m_Jpsi_%s", (isPbp?"Pbp":"PP") ));
  parIni[Form("sigma1_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{MassRatio,%s})", Form("sigma1_Psi2S_%s", (isPbp?"Pbp":"PP")), Form("sigma1_Jpsi_%s", (isPbp?"Pbp":"PP") ));
  parIni[Form("sigma2_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{MassRatio,%s})", Form("sigma2_Psi2S_%s", (isPbp?"Pbp":"PP")), Form("sigma2_Jpsi_%s", (isPbp?"Pbp":"PP") ));
  parIni[Form("alpha_Psi2S_%s", (isPbp?"Pbp":"PP"))]  = Form("RooFormulaVar::%s('@0',{%s})", Form("alpha_Psi2S_%s", (isPbp?"Pbp":"PP")), Form("alpha_Jpsi_%s", (isPbp?"Pbp":"PP")));
  parIni[Form("alpha2_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0',{%s})", Form("alpha2_Psi2S_%s", (isPbp?"Pbp":"PP")), Form("alpha2_Jpsi_%s", (isPbp?"Pbp":"PP")));
  parIni[Form("n_Psi2S_%s", (isPbp?"Pbp":"PP"))]      = Form("RooFormulaVar::%s('@0',{%s})", Form("n_Psi2S_%s", (isPbp?"Pbp":"PP")), Form("n_Jpsi_%s", (isPbp?"Pbp":"PP")));
  parIni[Form("n2_Psi2S_%s", (isPbp?"Pbp":"PP"))]     = Form("RooFormulaVar::%s('@0',{%s})", Form("n2_Psi2S_%s", (isPbp?"Pbp":"PP")), Form("n2_Jpsi_%s", (isPbp?"Pbp":"PP")));
  parIni[Form("f_Psi2S_%s", (isPbp?"Pbp":"PP"))]      = Form("RooFormulaVar::%s('@0',{%s})", Form("f_Psi2S_%s", (isPbp?"Pbp":"PP")), Form("f_Jpsi_%s", (isPbp?"Pbp":"PP")));
};

void setMassDefaultParameters(map<string, string> &parIni, bool isPbp, double numEntries)
{

  cout << "[INFO] Setting user undefined initial parameters to their default values" << endl;

  // DEFAULT SINGLE AND DOUBLE RATIO PARAMETERS
  if (parIni.count("RFrac2Svs1S_PbpvsPP")==0 || parIni["RFrac2Svs1S_PbpvsPP"]=="") { 
    parIni["RFrac2Svs1S_PbpvsPP"] = Form("%s[%.4f,%.4f,%.4f]", "RFrac2Svs1S_PbpvsPP", 0.26, -3.0, 3.0);
  }
  if (parIni.count(Form("RFrac2Svs1S_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("RFrac2Svs1S_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("RFrac2Svs1S_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("RFrac2Svs1S_%s", (isPbp?"Pbp":"PP")), 0.26, -2.0, 2.0);
  }

  // DEFAULT RANGE OF NUMBER OF EVENTS
  if (parIni.count(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.10f,%.10f,%.10f]", Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")), numEntries, -2.0*numEntries, 2.0*numEntries);
  }
  if (parIni.count(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.10f,%.10f,%.10f]", Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")), numEntries, -2.0*numEntries, 2.0*numEntries);
  }
  if (parIni.count(Form("N_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("N_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("N_Bkg_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.10f,%.10f,%.10f]", Form("N_Bkg_%s", (isPbp?"Pbp":"PP")), numEntries, -2.0*numEntries, 2.0*numEntries);
  }

  // DEFAULT SIGNAL MASS MODEL PARAMETERS 
  if (parIni.count(Form("m_Jpsi_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("m_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("m_Jpsi_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.6f,%.6f,%.6f]", Form("m_Jpsi_%s", (isPbp?"Pbp":"PP")), Mass.JPsi, Mass.JPsi-0.1, Mass.JPsi+0.1);
  }
  if (parIni.count(Form("m_Psi2S_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("m_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("m_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.6f,%.6f,%.6f]", Form("m_Psi2S_%s", (isPbp?"Pbp":"PP")), Mass.Psi2S, Mass.Psi2S-0.1, Mass.Psi2S+0.1);
  }
  if (parIni.count(Form("sigma1_Jpsi_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("sigma1_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("sigma1_Jpsi_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma1_Jpsi_%s", (isPbp?"Pbp":"PP")), 0.02, 0.005, 0.07);
  }
  if (parIni.count(Form("rSigma21_Jpsi_%s", (isPbp?"Pbp":"PP")))==0) {
    if (parIni.count(Form("sigma2_Jpsi_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("sigma2_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="") {
      parIni[Form("sigma2_Jpsi_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma2_Jpsi_%s", (isPbp?"Pbp":"PP")), 0.04, 0.01, 0.10);
    }
  } else {
    if (parIni[Form("rSigma21_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="") {
      parIni[Form("rSigma21_Jpsi_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("rSigma21_Jpsi_%s", (isPbp?"Pbp":"PP")), 2.0, 1.0, 4.0);
    }
    parIni[Form("sigma2_Jpsi_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("sigma2_Jpsi_%s", (isPbp?"Pbp":"PP")),
                                                                parIni[Form("rSigma21_Jpsi_%s", (isPbp?"Pbp":"PP"))].c_str(), Form("sigma1_Jpsi_%s", (isPbp?"Pbp":"PP") ));
  }      
  if (parIni.count(Form("sigma1_Psi2S_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("sigma1_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("sigma1_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma1_Psi2S_%s", (isPbp?"Pbp":"PP")), 0.02, 0.005, 0.07);
  }
  if (parIni.count(Form("rSigma21_Psi2S_%s", (isPbp?"Pbp":"PP")))==0) {
    if (parIni.count(Form("sigma2_Psi2S_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("sigma2_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="") {
      parIni[Form("sigma2_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma2_Psi2S_%s", (isPbp?"Pbp":"PP")), 0.04, 0.01, 0.10);
    }
  }
  else {
    if (parIni[Form("rSigma21_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="") {
      parIni[Form("rSigma21_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("rSigma21_Psi2S_%s", (isPbp?"Pbp":"PP")), 2.0, 1.0, 4.0);
    }
    parIni[Form("sigma2_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("sigma2_Psi2S_%s", (isPbp?"Pbp":"PP")), 
                                                                 parIni[Form("rSigma21_Psi2S_%s", (isPbp?"Pbp":"PP"))].c_str(), Form("sigma1_Psi2S_%s", (isPbp?"Pbp":"PP") ));
  }
  if (parIni.count(Form("alpha_Jpsi_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("alpha_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("alpha_Jpsi_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("alpha_Jpsi_%s", (isPbp?"Pbp":"PP")), 2.0, 0.5, 30.0);
  }
  if (parIni.count(Form("alpha2_Jpsi_%s", (isPbp?"Pbp":"PP")))==0) {
    parIni[Form("alpha2_Jpsi_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0',{%s})", Form("alpha2_Jpsi_%s", (isPbp?"Pbp":"PP")), Form("alpha_Jpsi_%s", (isPbp?"Pbp":"PP")));
  }
  else if (parIni[Form("alpha2_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="")
  {
    parIni[Form("alpha2_Jpsi_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("alpha2_Jpsi_%s", (isPbp?"Pbp":"PP")), 2.0, 0.5, 30.0);
  }
  if (parIni.count(Form("alpha_Psi2S_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("alpha_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("alpha_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("alpha_Psi2S_%s", (isPbp?"Pbp":"PP")), 2.0, 0.5, 30.0);
  }
  if (parIni.count(Form("alpha2_Psi2S_%s", (isPbp?"Pbp":"PP")))==0) {
    parIni[Form("alpha2_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0',{%s})", Form("alpha2_Psi2S_%s", (isPbp?"Pbp":"PP")), Form("alpha_Psi2S_%s", (isPbp?"Pbp":"PP")));
  }
  else if (parIni[Form("alpha2_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="")
  {
    parIni[Form("alpha2_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("alpha2_Psi2S_%s", (isPbp?"Pbp":"PP")), 2.0, 0.5, 30.0);
  }
  if (parIni.count(Form("n_Jpsi_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("n_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("n_Jpsi_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("n_Jpsi_%s", (isPbp?"Pbp":"PP")), 1.8, 0.5, 10.0);
  }
  if (parIni.count(Form("n2_Jpsi_%s", (isPbp?"Pbp":"PP")))==0) {
    parIni[Form("n2_Jpsi_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0',{%s})", Form("n2_Jpsi_%s", (isPbp?"Pbp":"PP")), Form("n_Jpsi_%s", (isPbp?"Pbp":"PP")));
  }
  else if (parIni[Form("n2_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="")
  {
    parIni[Form("n2_Jpsi_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("n2_Jpsi_%s", (isPbp?"Pbp":"PP")), 1.8, 0.5, 10.0);
  }
  if (parIni.count(Form("n_Psi2S_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("n_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("n_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("n_Psi2S_%s", (isPbp?"Pbp":"PP")), 1.8, 0.5, 10.0);
  }
  if (parIni.count(Form("n2_Psi2S_%s", (isPbp?"Pbp":"PP")))==0) {
    parIni[Form("n2_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("RooFormulaVar::%s('@0',{%s})", Form("n2_Psi2S_%s", (isPbp?"Pbp":"PP")), Form("n_Psi2S_%s", (isPbp?"Pbp":"PP")));
  }
  else if (parIni[Form("n2_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="")
  {
        parIni[Form("n2_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("n2_Psi2S_%s", (isPbp?"Pbp":"PP")), 1.8, 0.5, 10.0);
  }
  if (parIni.count(Form("f_Jpsi_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("f_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("f_Jpsi_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("f_Jpsi_%s", (isPbp?"Pbp":"PP")), 0.5, 0.0, 1.0);
  }
  if (parIni.count(Form("f_Psi2S_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("f_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="") {
    parIni[Form("f_Psi2S_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("f_Psi2S_%s", (isPbp?"Pbp":"PP")), 0.5, 0.0, 1.0);
  }

  // DEFAULT BACKGROUND MASS MODEL PARAMETERS
  if (parIni[Form("Model_Bkg_%s",(isPbp?"Pbp":"PP"))].find("ExpChebychev")!=std::string::npos) {
    if (parIni.count(Form("lambda1_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambda1_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
      parIni[Form("lambda1_Bkg_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda1_Bkg_%s", (isPbp?"Pbp":"PP")), 0.0, -10.0, 10.0);
    }
    if (parIni.count(Form("lambda2_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambda2_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
      parIni[Form("lambda2_Bkg_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda2_Bkg_%s", (isPbp?"Pbp":"PP")), 0.0, -10.0, 10.0);
    }
    if (parIni.count(Form("lambda3_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambda3_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
      parIni[Form("lambda3_Bkg_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda3_Bkg_%s", (isPbp?"Pbp":"PP")), 0.0, -10.0, 10.0);
    }
    if (parIni.count(Form("lambda4_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambda4_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
      parIni[Form("lambda4_Bkg_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda4_Bkg_%s", (isPbp?"Pbp":"PP")), 0.0, -10.0, 10.0);
    }
    if (parIni.count(Form("lambda5_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambda5_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
      parIni[Form("lambda5_Bkg_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda5_Bkg_%s", (isPbp?"Pbp":"PP")), 0.0, -10.0, 10.0);
    }
    if (parIni.count(Form("lambda6_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambda6_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
      parIni[Form("lambda6_Bkg_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda6_Bkg_%s", (isPbp?"Pbp":"PP")), 0.0, -10.0, 10.0);
    }
  }
  else if (parIni[Form("Model_Bkg_%s",(isPbp?"Pbp":"PP"))].find("Chebychev")!=std::string::npos) {
    if (parIni.count(Form("lambda1_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambda1_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
      parIni[Form("lambda1_Bkg_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda1_Bkg_%s", (isPbp?"Pbp":"PP")), 0.0, -2.0, 2.0);
    }
    if (parIni.count(Form("lambda2_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambda2_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
      parIni[Form("lambda2_Bkg_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda2_Bkg_%s", (isPbp?"Pbp":"PP")), 0.0, -2.0, 2.0);
    }
    if (parIni.count(Form("lambda3_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambda3_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
      parIni[Form("lambda3_Bkg_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda3_Bkg_%s", (isPbp?"Pbp":"PP")), 0.0, -2.0, 2.0);
    }
    if (parIni.count(Form("lambda4_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambda4_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
      parIni[Form("lambda4_Bkg_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda4_Bkg_%s", (isPbp?"Pbp":"PP")), 0.0, -2.0, 2.0);
    }
    if (parIni.count(Form("lambda5_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambda5_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
      parIni[Form("lambda5_Bkg_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda5_Bkg_%s", (isPbp?"Pbp":"PP")), 0.0, -2.0, 2.0);
    }
    if (parIni.count(Form("lambda6_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambda6_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
      parIni[Form("lambda6_Bkg_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda6_Bkg_%s", (isPbp?"Pbp":"PP")), 0.0, -2.0, 2.0);
    }
  } 
  else if (parIni[Form("Model_Bkg_%s",(isPbp?"Pbp":"PP"))].find("Exponential")!=std::string::npos) {
    if (parIni.count(Form("lambda1_Bkg_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("lambda1_Bkg_%s", (isPbp?"Pbp":"PP"))]=="") { 
      parIni[Form("lambda1_Bkg_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda1_Bkg_%s", (isPbp?"Pbp":"PP")), 0.0, -100.0, 100.0);
    }
  }
 
};


#endif // #ifndef buildCharmoniaMassModel_C
