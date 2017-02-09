#ifndef buildCharmoniaCtauResModel_C
#define buildCharmoniaCtauResModel_C

#include "Utilities/initClasses.h"

void setCtauResDefaultParameters(map<string, string> &parIni, bool isPbp, double numEntries);
void setCtauResMeanToZero(map<string, string> &parIni, bool isPbp);
bool defineCtauResModel(RooWorkspace& ws, string pdfType, string object, string varName, CtauModel model, map<string,string> parIni, bool isPbp, bool usePerEventError=true);


bool buildCharmoniaCtauResModel(RooWorkspace& ws, struct CharmModel model, map<string, string>  parIni, string dsName, 
                                string varName,              // Variable name
                                string pdfType,              // Pdf Type Name
                                bool isPbp,                 // Determine if we are working with Pbp (True) or PP (False)
                                bool usePerEventError,       // Determine if we use the perEventError technique
                                double  numEntries = 300000. // Number of entries in the dataset
                                )
{

  bool isMC = (dsName.find("MC")!=std::string::npos);
  bool incNonPrompt = (dsName.find("NOPR")!=std::string::npos);
  bool incJpsi = (dsName.find("JPSI")!=std::string::npos);
  if (isMC==false) { cout << "[ERROR] The ctau resolution fit can only be performed in MC!" << endl; return false; }

  string objInc = (incJpsi?"Jpsi":"Psi2S");
  string obj    = objInc + (incNonPrompt?"NoPR":"PR");

  // If the initial parameters are empty, set defaul parameter values
  setCtauResDefaultParameters(parIni, isPbp, numEntries);
  setCtauResMeanToZero(parIni, isPbp);

  // C r e a t e   m o d e l 

  if(!defineCtauResModel(ws, pdfType, obj, varName, model.CtauRes, parIni, isPbp, usePerEventError)) { 
    cout << "[ERROR] Defining the " << objInc  << (incNonPrompt?" Non-":" ") << "Prompt MC Ctau Resolution Model failed" << endl; return false; 
  }
  if (usePerEventError) {
    RooProdPdf pdf(Form("%s_%s_%s", pdfType.c_str(), obj.c_str(), (isPbp?"Pbp":"PP")), "", *ws.pdf(Form("pdfCTAUERR_%s_%s", objInc.c_str(), (isPbp?"Pbp":"PP"))),
                   Conditional( *ws.pdf(Form("%sCOND_%s_%s", pdfType.c_str(), obj.c_str(), (isPbp?"Pbp":"PP"))), RooArgList(*ws.var(varName.c_str())) )
                   ); 
    ws.import(pdf);
  } else {
    ws.factory(Form("SUM::%s(%s)", Form("%s_%s_%s", pdfType.c_str(), obj.c_str(), (isPbp?"Pbp":"PP")),
                    Form("%sCOND_%s_%s", pdfType.c_str(), obj.c_str(), (isPbp?"Pbp":"PP"))
                    ));
  }
  if ( !ws.var(Form("N_%s_%s", objInc.c_str(), (isPbp?"Pbp":"PP"))) ){ 
    ws.factory( parIni[Form("N_%s_%s", objInc.c_str(), (isPbp?"Pbp":"PP"))].c_str() ); 
  }
  ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_%s_%s", pdfType.c_str(), obj.c_str(), (isPbp?"Pbp":"PP")),
                  Form("%s_%s_%s", pdfType.c_str(), obj.c_str(), (isPbp?"Pbp":"PP")),
                  Form("N_%s_%s", objInc.c_str(), (isPbp?"Pbp":"PP"))
                  ));
  
  // Total PDF
  RooAbsPdf *themodel = themodel = new RooAddPdf(Form("%s_Tot_%s", pdfType.c_str(), (isPbp?"Pbp":"PP")), Form("%s_Tot_%s", pdfType.c_str(), (isPbp?"Pbp":"PP")),
                                                 *ws.pdf(Form("%s_%s_%s", pdfType.c_str(), obj.c_str(), (isPbp?"Pbp":"PP"))),
                                                 *ws.var(Form("N_%s_%s", objInc.c_str(), (isPbp?"Pbp":"PP")))
                                                 );
  ws.import(*themodel);

  setFixedVarsToContantVars(ws);
  
  return true;
};


bool defineCtauResModel(RooWorkspace& ws, string pdfType, string object, string varName, CtauModel model, map<string,string> parIni, bool isPbp, bool usePerEventError)
{

  if (ws.pdf(Form("%s_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")))) { 
    cout << Form("[ERROR] The %s Ctau Resolution Model has already been implemented!", object.c_str()) << endl;
    return false;
  }

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
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%sCOND_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), varName.c_str(),
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
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%sCOND_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), 
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
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%sCOND_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")), 
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
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%sCOND_%s_%s", pdfType.c_str(), object.c_str(), (isPbp?"Pbp":"PP")),
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


void setCtauResDefaultParameters(map<string, string> &parIni, bool isPbp, double numEntries)
{
  // DEFAULT RANGE OF NUMBER OF EVENTS
  if (parIni.count(Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("N_Jpsi_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Jpsi_%s", (isPbp?"Pbp":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count(Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")))==0 || parIni[Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))]=="") { 
    parIni[Form("N_Psi2S_%s", (isPbp?"Pbp":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Psi2S_%s", (isPbp?"Pbp":"PP")), numEntries, 0.0, numEntries*2.0);
  }

 // CTAU FIT PARAMETERS

 // Resolution Ctau Model

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

};


void setCtauResMeanToZero(map<string, string> &parIni, bool isPbp)
{
  parIni[Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("ctau1_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.0, 0.0, 0.0);
  parIni[Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("ctau2_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.0, 0.0, 0.0);
  parIni[Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("ctau3_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.0, 0.0, 0.0);
  parIni[Form("ctau4_CtauRes_%s", (isPbp?"Pbp":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("ctau4_CtauRes_%s", (isPbp?"Pbp":"PP")), 0.0, 0.0, 0.0);

  return;
};


#endif // #ifndef buildCharmoniaCtauResModel_C
