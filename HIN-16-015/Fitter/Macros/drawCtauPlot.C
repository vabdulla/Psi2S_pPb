#ifndef drawCtauPlot_C
#define drawCtauPlot_C

#include "Utilities/initClasses.h"

void setCtauRange(RooWorkspace& myws, RooPlot* frame, string dsName, bool setLogScale, vector<double> rangeErr, double excEvts=0.0);
void printCtauParameters(RooWorkspace myws, TPad* Pad, bool isPbp, string pdfName, bool isWeighted);

void drawCtauPlot(RooWorkspace& myws,   // Local workspace
                  string outputDir,     // Output directory
		  struct InputOpt opt,  // Variable with run information (kept for legacy purpose)
                  struct KinCuts cut,   // Variable with current kinematic cuts
                  map<string, string>  parIni,   // Variable containing all initial parameters
                  string plotLabel,     // The label used to define the output file name
                  // Select the type of datasets to fit
                  string DSTAG,         // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
                  bool isPbp,          // Define if it is Pbp (True) or PP (False)
                  // Select the type of object to fit
                  bool incJpsi,         // Includes Jpsi model
                  bool incPsi2S,        // Includes Psi(2S) model
                  bool incBkg,          // Includes Background model        
                  bool incPrompt,       // Includes Prompt ctau model       
                  bool incNonPrompt,    // Includes Non-Prompt ctau model
                  // Select the fitting options
                  bool plotPureSMC,     // Flag to indicate if we want to fit pure signal MC
                  // Select the drawing options
                  bool setLogScale,     // Draw plot with log scale
                  bool incSS,           // Include Same Sign data
                  double binWidth       // Bin width
                  ) 
{

  RooMsgService::instance().getStream(0).removeTopic(Caching);  
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  if (DSTAG.find("_")!=std::string::npos) DSTAG.erase(DSTAG.find("_"));

  string pdfTotName  = Form("pdfCTAU_Tot_%s", (isPbp?"Pbp":"PP"));
  string dsOSName = Form("dOS_%s_%s", DSTAG.c_str(), (isPbp?"Pbp":"PP"));
  if (plotPureSMC) dsOSName = Form("dOS_%s_%s_NoBkg", DSTAG.c_str(), (isPbp?"Pbp":"PP"));
  string dsOSNameCut = dsOSName+"_CTAUCUT";
  string hOSName = Form("dhCTAUERRTot_Tot_%s", (isPbp?"Pbp":"PP"));
  string hOSNameBkg  = Form("dhCTAUERR_Bkg_%s", (isPbp?"Pbp":"PP"));
  string hOSNameJpsi = Form("dhCTAUERR_Jpsi_%s", (isPbp?"Pbp":"PP"));
  string hOSNamePsi2S = Form("dhCTAUERR_Psi2S_%s", (isPbp?"Pbp":"PP"));
  string dsSSName = Form("dSS_%s_%s", DSTAG.c_str(), (isPbp?"Pbp":"PP"));

  bool isWeighted = myws.data(dsOSName.c_str())->isWeighted();
  bool isMC = (DSTAG.find("MC")!=std::string::npos);
  vector<double> range; range.push_back(cut.dMuon.ctau.Min); range.push_back(cut.dMuon.ctau.Max);

  double minRange = (double)(floor(range[0]*10.)/10.)+0.1;
  double maxRange = (double)(ceil(range[1]*10.)/10.)+0.1;
  if (!incNonPrompt) {
    if (abs(maxRange)>abs(minRange)) { minRange = -1.0*abs(maxRange); } else { maxRange = abs(minRange); }
  } else {
    minRange = -3.0;
    maxRange = 5.0;
  }
  Double_t outTot = myws.data(dsOSName.c_str())->sumEntries();
  Double_t outErr = myws.data(dsOSName.c_str())->reduce(Form("(ctau>%.6f || ctau<%.6f)", range[1], range[0]))->sumEntries();
  int nBins = min(int( round((maxRange - minRange)/binWidth) ), 1000);

  double normDSTot   = 1.0;  if (myws.data(dsOSNameCut.c_str()))  { normDSTot   = myws.data(dsOSName.c_str())->sumEntries()/myws.data(dsOSNameCut.c_str())->sumEntries();  }
  double normJpsi  = 1.0;  if (myws.data(hOSNameJpsi.c_str()))  { normJpsi  = myws.data(dsOSName.c_str())->sumEntries()*normDSTot/myws.data(hOSNameJpsi.c_str())->sumEntries();  }
  double normPsi2S = 1.0;  if (myws.data(hOSNamePsi2S.c_str())) { normPsi2S = myws.data(dsOSName.c_str())->sumEntries()*normDSTot/myws.data(hOSNamePsi2S.c_str())->sumEntries(); }
  double normBkg   = 1.0;  if (myws.data(hOSNameBkg.c_str()))   { normBkg   = myws.data(dsOSName.c_str())->sumEntries()*normDSTot/myws.data(hOSNameBkg.c_str())->sumEntries();   }
  double normTot   = 1.0;  if (myws.data(hOSName.c_str()))  { normTot   = myws.data(dsOSName.c_str())->sumEntries()*normDSTot/myws.data(hOSName.c_str())->sumEntries();  }

  // Create the main plot of the fit
  RooPlot*   frame     = myws.var("ctau")->frame(Bins(nBins), Range(minRange, maxRange));
  frame->updateNormVars(RooArgSet(*myws.var("invMass"), *myws.var("ctau"), *myws.var("ctauErr"))) ;
  myws.data(dsOSName.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  
  if (incBkg && (!incJpsi && !incPsi2S)) {
    myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKG"),  Normalization(normBkg, RooAbsReal::NumEvent), NumCPU(32),
                                         ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(hOSNameBkg.c_str()), kTRUE),
                                         FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), Precision(1e-4)
                                         );
    myws.data(dsOSName.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
    myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKGPR"),Components(RooArgSet(*myws.pdf(Form("pdfCTAU_BkgPR_%s", (isPbp?"Pbp":"PP"))))),
                                         ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(hOSNameBkg.c_str()), kTRUE),
                                         Normalization(normBkg, RooAbsReal::NumEvent),
                                         LineColor(kRed+2), Precision(1e-4), NumCPU(32)
                                         );
    myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKGNoPR"),Components(RooArgSet(*myws.pdf(Form("pdfCTAU_BkgNoPR_%s", (isPbp?"Pbp":"PP"))))),
                                         ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(hOSNameBkg.c_str()), kTRUE),
                                         Normalization(normBkg, RooAbsReal::NumEvent),
                                         LineColor(kOrange+10), Precision(1e-4), NumCPU(32)
                                         );
    myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PDF"),  Normalization(normBkg, RooAbsReal::NumEvent), NumCPU(32),
                                         ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(hOSNameBkg.c_str()), kTRUE),
                                         LineColor(kBlack), Precision(1e-4)
                                         );
  } 
  if (!incBkg || (incJpsi || incPsi2S)) {
    if (incPrompt) {
      if (incJpsi) {
        int COLOR[] = { kGreen+5, kRed+2, kBlue+3, kViolet-5};
        for (int i=1; i<5; i++) {
          if (myws.pdf(Form("pdfCTAURES%d_JpsiPR_%s", i,(isPbp?"Pbp":"PP")))){
            myws.pdf(pdfTotName.c_str())->plotOn(frame,Name(Form("PDFJPSI%d",i)),Components(RooArgSet(*myws.pdf(Form("pdfCTAURES%d_JpsiPR_%s", i, (isPbp?"Pbp":"PP"))))),
                                                 ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(dsOSNameCut.c_str()), kTRUE),
                                                 Normalization(normDSTot, RooAbsReal::NumEvent),
                                                 LineColor(COLOR[i-1]), Precision(1e-5), NumCPU(32)
                                                 );
          }
        }
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PDF"),Components(RooArgSet(*myws.pdf(Form("pdfCTAU_JpsiPR_%s", (isPbp?"Pbp":"PP"))))),
                                             ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(dsOSNameCut.c_str()), kTRUE),
                                             Normalization(normDSTot, RooAbsReal::NumEvent),
                                             LineColor(kBlack), Precision(1e-5), NumCPU(32)
                                             );
      }
      if (incPsi2S) {
        int COLOR[] = { kGreen+5, kRed+2, kBlue+3, kViolet-5};
        for (int i=1; i<5; i++) {
          if (myws.pdf(Form("pdfCTAURES%d_Psi2SPR_%s", i,(isPbp?"Pbp":"PP")))){
            myws.pdf(pdfTotName.c_str())->plotOn(frame,Name(Form("PDFPSI2S%d",i)),Components(RooArgSet(*myws.pdf(Form("pdfCTAURES%d_Psi2SPR_%s", i, (isPbp?"Pbp":"PP"))))),
                                                 ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(dsOSNameCut.c_str()), kTRUE),
                                                 Normalization(normDSTot, RooAbsReal::NumEvent),
                                                 LineColor(COLOR[i-1]), Precision(1e-5), NumCPU(32)
                                                 );
          }
        }
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PDF"),Components(RooArgSet(*myws.pdf(Form("pdfCTAU_Psi2SPR_%s", (isPbp?"Pbp":"PP"))))),
                                             ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(dsOSNameCut.c_str()), kTRUE),
                                             Normalization(normDSTot, RooAbsReal::NumEvent),
                                             LineColor(kBlack), Precision(1e-5), NumCPU(32)
                                             );
      }
    }
    if (incNonPrompt) {
      if (incJpsi) {
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PDF"),Components(RooArgSet(*myws.pdf(Form("pdfCTAU_JpsiNoPR_%s", (isPbp?"Pbp":"PP"))))),
                                             ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(hOSNameJpsi.c_str()), kTRUE),
                                             Normalization(normJpsi, RooAbsReal::NumEvent),
                                             LineColor(kBlack), Precision(1e-5), NumCPU(32)
                                             );
      }
      if (incPsi2S) {
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PDF"),Components(RooArgSet(*myws.pdf(Form("pdfCTAU_Psi2SNoPR_%s", (isPbp?"Pbp":"PP"))))),
                                             ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(hOSNamePsi2S.c_str()), kTRUE),
                                             ProjWData(*myws.data(hOSNamePsi2S.c_str())), 
                                             Normalization(normPsi2S, RooAbsReal::NumEvent),
                                             LineColor(kBlack), Precision(1e-5), NumCPU(32)
                                             );
      }
    } 
  }
  if (incSS) { 
    myws.data(dsSSName.c_str())->plotOn(frame, Name("dSS"), MarkerColor(kRed), LineColor(kRed), MarkerSize(1.2)); 
  }
  myws.data(dsOSName.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  
  
  // set the CMS style
  setTDRStyle();

  // Create the pull distribution of the fit 
  RooHist *hpull = frame->pullHist(0, "PDF", true);
  hpull->SetName("hpull");
  RooPlot* frame2 = myws.var("ctau")->frame(Title("Pull Distribution"), Bins(nBins), Range(minRange, maxRange));
  frame2->addPlotable(hpull, "PX"); 
  
  // Create the main canvas
  TCanvas *cFig  = new TCanvas(Form("cCtauFig_%s", (isPbp?"Pbp":"PP")), "cCtauFig",800,800);
  TPad    *pad1  = new TPad(Form("pad1_%s", (isPbp?"Pbp":"PP")),"",0,0.23,1,1);
  TPad    *pad2  = new TPad(Form("pad2_%s", (isPbp?"Pbp":"PP")),"",0,0,1,.228);
  TLine   *pline = new TLine(minRange, 0.0, maxRange, 0.0);
  
  TPad *pad4 = new TPad("pad4","This is pad4",0.55,0.46,0.97,0.87);
  pad4->SetFillStyle(0);
  pad4->SetLeftMargin(0.28);
  pad4->SetRightMargin(0.10);
  pad4->SetBottomMargin(0.21);
  pad4->SetTopMargin(0.072);

  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("");
  frame->GetXaxis()->CenterTitle(kTRUE);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetTitleFont(42);
  frame->GetXaxis()->SetTitleOffset(3);
  frame->GetXaxis()->SetLabelOffset(3);
  frame->GetYaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetTitleOffset(1.7);
  frame->GetYaxis()->SetTitleFont(42);
  setCtauRange(myws, frame, dsOSNameCut, setLogScale, range, outErr);
 
  cFig->cd();
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.4);
  pad2->SetFillStyle(4000); 
  pad2->SetFrameFillStyle(4000); 
  pad1->SetBottomMargin(0.015); 
  //plot fit
  pad1->Draw();
  pad1->cd(); 
  frame->Draw();

  printCtauParameters(myws, pad1, isPbp, pdfTotName, isWeighted);
  pad1->SetLogy(setLogScale);

  // Drawing the text in the plot
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.032);
  float dy = 0; 
  
  t->SetTextSize(0.03);
  t->DrawLatex(0.21, 0.86-dy, "2015 HI Soft Muon ID"); dy+=0.045;
  if (isPbp) {
    t->DrawLatex(0.21, 0.86-dy, "HLT_HIL1DoubleMu0_v1"); dy+=0.045;
  } else {
    t->DrawLatex(0.21, 0.86-dy, "HLT_HIL1DoubleMu0_v1"); dy+=0.045;
  } 
  t->DrawLatex(0.21, 0.86-dy, Form("%.1f #leq p_{T}^{#mu#mu} < %.1f GeV/c",cut.dMuon.Pt.Min,cut.dMuon.Pt.Max)); dy+=0.045;
  t->DrawLatex(0.21, 0.86-dy, Form("%.1f #leq |y^{#mu#mu}| < %.1f",cut.dMuon.AbsRap.Min,cut.dMuon.AbsRap.Max)); dy+=0.045;
  if (isPbp) {t->DrawLatex(0.21, 0.86-dy, Form("Cent. %d-%d%%", (int)(cut.Centrality.Start/2), (int)(cut.Centrality.End/2))); dy+=0.045;}
  if (outErr>0.0) {
    if (isPbp && isMC) {
      t->DrawLatex(0.21, 0.86-dy, Form("Excl: (%.4f%%) %.4f evts", (outErr*100.0/outTot), outErr)); dy+=1.5*0.045;
    } else {
      t->DrawLatex(0.21, 0.86-dy, Form("Excl: (%.4f%%) %.0f evts", (outErr*100.0/outTot), outErr)); dy+=1.5*0.045;
    }
  }

  // Drawing the Legend
  double ymin = 0.7802;
  if (incPsi2S && incJpsi && incSS)  { ymin = 0.7202; } 
  if (incPsi2S && incJpsi && !incSS) { ymin = 0.7452; }
  TLegend* leg = new TLegend(0.5175, ymin, 0.7180, 0.8809); leg->SetTextSize(0.03);
  leg->AddEntry(frame->findObject("dOS"), (incSS?"Opposite Charge":"Data"),"pe");
  if (incSS) { leg->AddEntry(frame->findObject("dSS"),"Same Charge","pe"); }
  if(frame->findObject("PDF")) { leg->AddEntry(frame->findObject("PDF"),"Total fit","l"); }
  if((incBkg && (incJpsi || incPsi2S)) && frame->findObject("BKG")) { leg->AddEntry(frame->findObject("BKG"),"Background","fl");  }
  if(incBkg && incJpsi && frame->findObject("JPSI")) { leg->AddEntry(frame->findObject("JPSI"),"J/#psi PDF","l"); }
  if(incBkg && incPsi2S && frame->findObject("PSI2S")) { leg->AddEntry(frame->findObject("PSI2S"),"#psi(2S) PDF","l"); }
  for (int i=1; i<5; i++) {
    if(incJpsi && frame->findObject(Form("PDFJPSI%d",i))) { leg->AddEntry(frame->findObject(Form("PDFJPSI%d",i)),Form("PR Gauss %d", i),"l"); }
  }
  leg->Draw("same");

  //Drawing the title
  TString label;
  if (isPbp) {
    if (opt.Pbp.RunNb.Start==opt.Pbp.RunNb.End){
      label = Form("Pbp Run %d", opt.Pbp.RunNb.Start);
    } else {
      label = Form("%s [%s %d-%d]", "Pbp", "HIOniaL1DoubleMu0", opt.Pbp.RunNb.Start, opt.Pbp.RunNb.End);
    }
  } else {
    if (opt.pp.RunNb.Start==opt.pp.RunNb.End){
      label = Form("PP Run %d", opt.pp.RunNb.Start);
    } else {
      label = Form("%s [%s %d-%d]", "PP", "DoubleMu0", opt.pp.RunNb.Start, opt.pp.RunNb.End);
    }
  }
  
  //CMS_lumi(pad1, isPbp ? 105 : 104, 33, label);
  CMS_lumi(pad1, isPbp ? 108 : 107, 33, "");
  gStyle->SetTitleFontSize(0.05);
  
  pad1->Update();
  cFig->cd(); 

  //---plot pull
  pad2->Draw();
  pad2->cd();
    
  frame2->SetTitle("");
  frame2->GetYaxis()->CenterTitle(kTRUE);
  frame2->GetYaxis()->SetTitleOffset(0.4);
  frame2->GetYaxis()->SetTitleSize(0.1);
  frame2->GetYaxis()->SetLabelSize(0.1);
  frame2->GetYaxis()->SetTitle("Pull");
  frame2->GetXaxis()->CenterTitle(kTRUE);
  frame2->GetXaxis()->SetTitleOffset(1);
  frame2->GetXaxis()->SetTitleSize(0.12);
  frame2->GetXaxis()->SetLabelSize(0.1);
  frame2->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  frame2->GetYaxis()->SetRangeUser(-7.0, 7.0);

  frame2->Draw(); 
  
  // *** Print chi2/ndof 
  printChi2(myws, pad2, frame, "ctau", dsOSName.c_str(), pdfTotName.c_str(), nBins);
  
  pline->Draw("same");
  pad2->Update();
  
  bool SB = (incBkg&&(!incPsi2S&&!incJpsi));
  // Save the plot in different formats
  gSystem->mkdir(Form("%sctau%s/%s/plot/root/", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str()), kTRUE); 
  cFig->SaveAs(Form("%sctau%s/%s/plot/root/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str(), "CTAU", DSTAG.c_str(), (isPbp?"Pbp":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  gSystem->mkdir(Form("%sctau%s/%s/plot/png/", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str()), kTRUE);
  cFig->SaveAs(Form("%sctau%s/%s/plot/png/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.png", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str(), "CTAU", DSTAG.c_str(), (isPbp?"Pbp":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  gSystem->mkdir(Form("%sctau%s/%s/plot/pdf/", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str()), kTRUE);
  cFig->SaveAs(Form("%sctau%s/%s/plot/pdf/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.pdf", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str(), "CTAU", DSTAG.c_str(), (isPbp?"Pbp":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  
  cFig->Clear();
  cFig->Close();

}

void setCtauRange(RooWorkspace& myws, RooPlot* frame, string dsName, bool setLogScale, vector<double> rangeErr, double excEvts)
{ 
  // Find maximum and minimum points of Plot to rescale Y axis
  TH1* h = myws.data(dsName.c_str())->createHistogram("hist", *myws.var("ctau"), Binning(frame->GetNbinsX(),frame->GetXaxis()->GetXmin(),frame->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));

  Double_t Yup(0.),Ydown(0.);
  
  bool isMC = false;
  if (dsName.find("MC")!=std::string::npos) isMC = true;

  if(setLogScale)
  {
    if (isMC) {
      Yup = YMax*TMath::Power((YMax/YMin), (0.5/(1.0-0.5-0.2)));
      Ydown = YMin/(TMath::Power((YMax/YMin), (0.2/(1.0-0.5-0.2))));
    } 
    else {
      Yup = YMax*TMath::Power((YMax/0.1), 0.5);
      Ydown = 0.1;
    }
  }
  else
  {
    Yup = YMax+(YMax-0.0)*0.5;
    Ydown = 0.0;
  }
  frame->GetYaxis()->SetRangeUser(Ydown,Yup);
  delete h;

  if (excEvts>0.0) {
    TLine   *minline = new TLine(rangeErr[0], 0.0, rangeErr[0], (setLogScale?(Ydown*TMath::Power((Yup/Ydown),0.4)):(Ydown + (Yup-Ydown)*0.4)));
    minline->SetLineStyle(2);
    minline->SetLineColor(1);
    minline->SetLineWidth(3);
    frame->addObject(minline);
    TLine   *maxline = new TLine(rangeErr[1], 0.0, rangeErr[1], (setLogScale?(Ydown*TMath::Power((Yup/Ydown),0.4)):(Ydown + (Yup-Ydown)*0.4)));
    maxline->SetLineStyle(2);
    maxline->SetLineColor(1);
    maxline->SetLineWidth(3);
    frame->addObject(maxline);
  }
};


void printCtauParameters(RooWorkspace myws, TPad* Pad, bool isPbp, string pdfName, bool isWeighted)
{
  Pad->cd();
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.026); float dy = 0.045; 
  RooArgSet* Parameters = (RooArgSet*)myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("invMass"), *myws.var("ctau"), *myws.var("ctauErr")))->selectByAttrib("Constant",kFALSE);
  TIterator* parIt = Parameters->createIterator(); 
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    stringstream ss(it->GetName()); string s1, s2, s3, label; 
    getline(ss, s1, '_'); getline(ss, s2, '_'); getline(ss, s3, '_');
    // Parse the parameter's labels
    if(s1=="invMass"){continue;} else if(s1=="ctau"){continue;} else if(s1=="MassRatio"){continue;}   
    else if(s1=="One"){continue;} else if(s1=="mMin"){continue;} else if(s1=="mMax"){continue;}
    else if(s1.find("sigma")!=std::string::npos || s1.find("lambda")!=std::string::npos || s1.find("alpha")!=std::string::npos){
      s1=Form("#%s",s1.c_str());
    }

    if(s2=="CtauRes")  { s2="Res";   } 
    else if(s2=="JpsiNoPR")  { s2="J/#psi[NoPR]";   } 
    else if(s2=="JpsiPR")  { s2="J/#psi[PR]";   } 
    else if(s2=="Jpsi" && (s1=="N" || s1=="b"))  { s2="J/#psi";   } 
    else if(s2=="Psi2SNoPR") { s2="#psi(2S)[NoPR]"; }
    else if(s2=="Psi2SPR") { s2="#psi(2S)[PR]"; } 
    else if(s2=="Psi2S" && (s1=="N" || s1=="b"))  { s2="#psi(2S)";   } 
    else if(s2=="BkgNoPR")   { s2="bkg[NoPR]";      } 
    else if(s2=="BkgPR")   { s2="bkg[PR]";      } 
    else if(s2=="Bkg" && (s1=="N" || s1=="b"))   { s2="bkg";      }
    else {continue;}
    if(s3!=""){
      label=Form("%s_{%s}^{%s}", s1.c_str(), s2.c_str(), s3.c_str());
    }
    // Print the parameter's results
    if(s1=="N"){ 
      t->DrawLatex(0.69, 0.75-dy, Form((isWeighted?"%s = %.6f#pm%.6f ":"%s = %.0f#pm%.0f "), label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else if(s1=="b"){ 
      t->DrawLatex(0.69, 0.75-dy, Form("%s = %.4f#pm%.4f ", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else if(s1.find("sigma")!=std::string::npos) { 
      if (s2.find("Res")!=std::string::npos) {
        t->DrawLatex(0.69, 0.75-dy, Form("%s = %.2f#pm%.2f", label.c_str(), it->getValV(), it->getError())); dy+=0.045;
      } else {
        t->DrawLatex(0.69, 0.75-dy, Form("%s = %.2f#pm%.2f mm", label.c_str(), it->getValV(), it->getError())); dy+=0.045;
      }
    }
    else if(s1.find("lambda")!=std::string::npos){ 
      t->DrawLatex(0.69, 0.75-dy, Form("%s = %.4f#pm%.4f", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else if(s1.find("ctau")!=std::string::npos){ 
      t->DrawLatex(0.69, 0.75-dy, Form("%s = %.4f#pm%.4f mm", (label.insert(1, string("#"))).c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else { 
      t->DrawLatex(0.69, 0.75-dy, Form("%s = %.4f#pm%.4f", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
  }
};


#endif // #ifndef drawCtauPlot_C
