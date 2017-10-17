#include "fitterShape.h"

#include "TROOT.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TF1.h"
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TIterator.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#ifndef __CINT__
#include "RooGlobalFunc.h"
#include "RooMsgService.h"
#endif
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
using namespace RooFit;
using namespace fitterShape;
const int fitMassLo=60;
const int fitMassHi=120;
RooFitResult *fitBin(
  string histName="pass_ptBin0_etaBin0",
  string dataFileName="",
  string templateFileName="",
  shapeType signalModel=kTemplateConvDas,
  shapeType bkgModel=kBkgErfcExp, 
  TString plotTitle="",
  string outputDir=""
  //metadata file name with kinematics info, initial params, info for plot 
) {
  //if(outputDir!="" && outputDir[outputDir.size()-1]!='/') outputDir=outputDir+"/";
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  
  // Open output file for the fit results
  ofstream fitResultFile;
  char txtfname[1000];    
  sprintf(txtfname,"%s/plots/fitres_%s.txt",outputDir.c_str(),histName.c_str());
  fitResultFile.open(txtfname);
  assert(fitResultFile.is_open());
 
  // Open the root files for reading
  TFile *dataFile=TFile::Open(dataFileName.c_str(),"READ"); assert(dataFile && dataFile->IsOpen());
  TH1D *dataHist=(TH1D*)dataFile->Get(histName.c_str()); assert(dataHist); dataHist->SetDirectory(0);
  TFile *templateFile; TH1D *templateHist;
  bool useMCTemplate=(templateFileName!="");
  if(useMCTemplate) {
    templateFile=TFile::Open(templateFileName.c_str(),"READ"); assert(templateFile);
    templateHist=(TH1D*)templateFile->Get(histName.c_str());
  }  

  //todo: declare models somewhere else
  RooRealVar m("m","mass",40.,140.);
  m.setBins(15000);

  // Signal model: template conv. CB for now
  //RooRealVar *mean  = new RooRealVar("sigCB_mean" ,"sigCB_mean" ,-2,-10,10);
  //RooRealVar *sigma = new RooRealVar("sigCB_sigma","sigCB_sigma",0.1,0.01,5);
  //RooRealVar *alpha = new RooRealVar("sigCB_alpha","sigCB_alpha",5,0,20);
  //RooRealVar *n     = new RooRealVar("sigCB_n"    ,"sigCB_n"    ,1,0,10);
  //RooCBShape *cb = new RooCBShape("cb","cb",m,*mean,*sigma,*alpha,*n);
  //RooDataHist *templateRDH = new RooDataHist("templateRDH","templateRDH",RooArgSet(m),templateHist);
  //RooHistPdf  *templateRHP = new RooHistPdf("templateRHP","templateRHP",m,*templateRDH,2);
  //RooAbsPdf *signalModelPdf=new RooFFTConvPdf("signalModel","signalModel",m,*templateRHP,*cb);
  
  // Signal model : das distribution
  fitterShape::fitterShapeBase *sigFitterShape, *bkgFitterShape;
 
  // if statement
  switch(signalModel) {
    case kTemplateConvDas    : sigFitterShape=new fitterShape::templateConvDas(m,templateHist); break;
    case kTemplateConvGaus   : sigFitterShape=new fitterShape::templateConvGaus(m,templateHist); break;
    case kTemplateBWConvGaus : sigFitterShape=new fitterShape::templateBWConvGaus(m,templateHist); break;
    case kBreitWignerConvDas : sigFitterShape=new fitterShape::breitWignerConvDas(m); break;
    default                  : printf("Unsupported signal model\n"); assert(0); break;
  }
  switch(bkgModel) {
    case kBkgErfcExp         : bkgFitterShape=new fitterShape::bkgErfcExp(m); break;
    case kBkgErfcExpPlusExp  : bkgFitterShape=new fitterShape::bkgErfcExpPlusExp(m); break;
    case kBkgDasPlusExp      : bkgFitterShape=new fitterShape::bkgDasPlusExp(m); break;
    case kBkgDas             : bkgFitterShape=new fitterShape::bkgDas(m); break;
    case kBkgDoubleExp       : bkgFitterShape=new fitterShape::bkgDoubleExp(m); break;
    default                  : printf("Unsupported background model\n"); assert(0); break;
  }
  //sigFitterShape->initializeParams(map);
  //bkgFitterShape->initializeParams(map)
  

  // Signal+background model
  double Ndata = dataHist->Integral( dataHist->FindBin((double)fitMassLo), dataHist->FindBin((double)fitMassHi-0.01) );
  RooRealVar Ntotal = RooRealVar("Ntotal","Ntotal",Ndata,0,Ndata); Ntotal.setConstant(kTRUE);
  RooRealVar Nbkg = RooRealVar("Nbkg","Nbkg",0.2*Ndata, 0, .98*Ndata);
  RooFormulaVar Nsig = RooFormulaVar("Nsig","Nsig","Ntotal - Nbkg", RooArgList(Ntotal,Nbkg));
  RooAbsPdf *totalPdf = new RooAddPdf("totalPdf","totalPdf", RooArgList(*(sigFitterShape->theShape),*(bkgFitterShape->theShape)), RooArgList(Nsig,Nbkg));
  // Data histogram in RooFit
  RooDataHist *dataRDH = new RooDataHist("dataRDH","dataRDH",RooArgSet(m), dataHist);
  m.setBins(15000);
  m.setRange(fitMassLo,fitMassHi);
  
  // Initial parameter estimation
  //(bkgFitterShape->theShape)->fitTo(*dataRDH,RooFit::Strategy(0),RooFit::NumCPU(1),Range(fitMassLo,80)); 
  //(bkgFitterShape->theShape)->fitTo(*dataRDH,RooFit::Strategy(0),RooFit::NumCPU(1),Range(100,fitMassHi)); 
  //(sigFitterShape->theShape)->fitTo(*dataRDH,RooFit::Strategy(0),RooFit::NumCPU(1),Range(86,96)); 
  
  //TF1 *bkgEstimator = new TF1("bkgEstimator","expo",fitMassLo,fitMassHi);
  //dataHist->Fit(bkgEstimator, "N0Q", "goff", fitMassLo,80);
  //Nbkg.setVal(TMath::Max(0.07*Ndata,TMath::Min(.95*Ndata,bkgEstimator->Integral(fitMassLo,fitMassHi))));
  //if(bkgModel==kBkgErfcExp) {
  //  //((fitterShape::bkgErfcExp*)bkgFitterShape)->mpv->setVal(dataHist->GetMean());
  //  ((fitterShape::bkgErfcExp*)bkgFitterShape)->gamma->setVal(bkgEstimator->GetParameter(1));
  //} else if(bkgModel==kBkgErfcExpPlusExp) {
  //  ((fitterShape::bkgErfcExpPlusExp*)bkgFitterShape)->gamma->setVal(bkgEstimator->GetParameter(1));
  //  //((fitterShape::bkgErfcExpPlusExp*)bkgFitterShape)->mpv->setVal(dataHist->GetMean());
  //}
  TF1 *peakEstimator = new TF1("peakEstimator","gaus",fitMassLo,fitMassHi);
  if(signalModel==kTemplateConvDas||signalModel==kTemplateConvGaus||signalModel==kTemplateBWConvGaus||signalModel==kBreitWignerConvDas) {
   //TODO: estimate this with a cheap gaussian or pol2 fit in small zmass window
   //double peakPos=dataHist->GetBinCenter(dataHist->GetMaximumBin());
   //double mcPeakPos=templateHist->GetBinCenter(dataHist->GetMaximumBin());
   dataHist->Fit(peakEstimator, "N0Q", "goff", 86.2, 96.2);
   double peakPos=peakEstimator->GetParameter(1);
   double peakSigma=peakEstimator->GetParameter(2);
   double mcPeakPos, mcPeakSigma;
   if(signalModel==kTemplateConvDas||signalModel==kTemplateConvGaus||signalModel==kTemplateBWConvGaus) {
    templateHist->Fit(peakEstimator, "N0Q", "goff", 86.2, 96.2);
    mcPeakPos=peakEstimator->GetParameter(1);
    mcPeakSigma=peakEstimator->GetParameter(2);
   } else if(signalModel==kBreitWignerConvDas) {
    mcPeakPos=91.1876;
    mcPeakSigma=2.495/2.355;
   }
   double meanEst=peakPos-mcPeakPos;
   double sigmaEst=TMath::Min(6.,TMath::Max(0.01, peakSigma-mcPeakSigma));
   if(peakPos<100. && peakPos>80.) {
    if(signalModel==kTemplateConvDas) {
     ((fitterShape::templateConvDas*)sigFitterShape)->mean->setRange(meanEst-2.,meanEst+2.);
     ((fitterShape::templateConvDas*)sigFitterShape)->mean->setVal(meanEst);
     //((fitterShape::templateConvDas*)sigFitterShape)->sigma->setRange(TMath::Max(sigmaEst-.5,0.001),sigmaEst+.5);
     ((fitterShape::templateConvDas*)sigFitterShape)->sigma->setRange(.001,sigmaEst+.5);
     ((fitterShape::templateConvDas*)sigFitterShape)->sigma->setVal(sigmaEst);
    } else if(signalModel==kTemplateConvGaus) {
     ((fitterShape::templateConvGaus*)sigFitterShape)->mean->setRange(meanEst-2.,meanEst+2.);
     ((fitterShape::templateConvGaus*)sigFitterShape)->mean->setVal(peakPos-mcPeakPos);
     ((fitterShape::templateConvGaus*)sigFitterShape)->sigma->setRange(.001,sigmaEst+.5);
     ((fitterShape::templateConvGaus*)sigFitterShape)->sigma->setVal(sigmaEst);
    } else if(signalModel==kTemplateBWConvGaus) {
     ((fitterShape::templateBWConvGaus*)sigFitterShape)->mean->setRange(meanEst-2.,meanEst+2.);
     ((fitterShape::templateBWConvGaus*)sigFitterShape)->mean->setVal(meanEst);
     ((fitterShape::templateBWConvGaus*)sigFitterShape)->sigma->setRange(0.001,sigmaEst+.5);
     ((fitterShape::templateBWConvGaus*)sigFitterShape)->sigma->setVal(sigmaEst);
    } else if(signalModel==kBreitWignerConvDas) {
     ((fitterShape::breitWignerConvDas*)sigFitterShape)->mean->setRange(meanEst-2.,meanEst+2.);
     ((fitterShape::breitWignerConvDas*)sigFitterShape)->mean->setVal(meanEst);
     ((fitterShape::breitWignerConvDas*)sigFitterShape)->sigma->setRange(0.001,sigmaEst+1.);
     ((fitterShape::breitWignerConvDas*)sigFitterShape)->sigma->setVal(sigmaEst);
    }
  }}
  //delete bkgEstimator;
  delete peakEstimator;
  
  RooFitResult *fitResult=0; //prefit
  m.setRange("Zpeak",86,96);
  //totalPdf->fitTo(*dataRDH,RooFit::Extended(),RooFit::Strategy(1),RooFit::NumCPU(1),Range(fitMassLo,80));
  m.setRange("lowMass",fitMassLo,80);
  //fitResult=totalPdf->fitTo(*dataRDH,RooFit::Extended(),RooFit::Strategy(1),RooFit::NumCPU(1),Range("lowMass"));

  int bombingRuns;
  int nInitFloatingPars=((RooArgSet*)totalPdf->getParameters(m)->selectByAttrib("Constant",kFALSE))->getSize();
  RooArgSet *qFloatingPars;
  TIterator *parIter;
  RooRealVar *qPar, *worstPar;  //questionable parameter
  double relParError=0, worstRelParError=0;
  vector<RooRealVar*> parsFrozen;
  // First perform the fit to low mass only
  bombingRuns=1; while(bombingRuns<=nInitFloatingPars-1) {
   fitResult=totalPdf->fitTo(*dataRDH,RooFit::Extended(),RooFit::Strategy(1),RooFit::NumCPU(1),RooFit::Save(),Range("lowMass"));
   bool isGoodFit = fitResult->covQual()>1 && fitResult->status()==0 && fitResult->edm()<1e-3;
   if(isGoodFit) { 
    printf("\nFit result is OK (edm %f, cov qual %d, HESSE status %d, MIGRAD status %d)\n", fitResult->edm(), fitResult->covQual(), fitResult->status()/100, fitResult->status()%10);
    break; // ok with covariance matrix being accurate (3) or made pos-def (2)
   }
   // We know something is wrong with covariance matrix, need to find the problem
   qPar=0; worstPar=0;relParError=0, worstRelParError=0;
   printf("\nBkg est. run %d: edm %f, cov qual %d, HESSE status %d, MIGRAD status %d, need to freeze a parameter\n", bombingRuns, fitResult->edm(), fitResult->covQual(), fitResult->status()/100, fitResult->status()%10);
   qFloatingPars = (RooArgSet*)totalPdf->getParameters(m)->selectByAttrib("Constant",kFALSE);
   // Start looking for questionable parameters
   parIter=qFloatingPars->createIterator(); qPar=(RooRealVar*)parIter->Next(); unsigned iPar=0;
   while(qPar) {
    relParError = qPar->getError() / qPar->getVal();
    if(
     strcmp("Nbkg",qPar->GetName())!=0 && // Cannot treat the background normalization as a questionable parameter
     relParError >= 0.5 && // If the parameter's error is 100% or more, it's questionable
     relParError>worstRelParError
    ) { worstRelParError=relParError; worstPar=qPar;}
    qPar=(RooRealVar*)parIter->Next(); iPar++;
   }
   if(worstPar) {
    printf("The worst parameter is %s , relative parameter error %f. Freezing it and retrying the fit...\n", worstPar->GetName(), worstRelParError);
    worstPar->setConstant();
    parsFrozen.push_back(worstPar);
   } else {
    printf("Couldn't find a questionable parameter, break\n"); break;
   }
   bombingRuns++;
  
  }
  // Unfreeze those parameters that we froze
  for(unsigned long i=0;i<parsFrozen.size();i++) parsFrozen[i]->setConstant(false);
  parsFrozen.clear();

  // Perform the actual fit
  bombingRuns=1; while(bombingRuns<=nInitFloatingPars-1) {
   TString rangeName=Form("fitRange_%d",bombingRuns);
   m.setRange(rangeName,fitMassLo,fitMassHi);
   fitResult = totalPdf->fitTo(*dataRDH,RooFit::Extended(),RooFit::Strategy(2),RooFit::NumCPU(1),RooFit::Save(),Range(rangeName));
   bool isGoodFit = fitResult->covQual()>1 && fitResult->status()==0 && fitResult->edm()<1e-3;
   if(isGoodFit) { 
    printf("\nFit result is OK (edm %f, cov qual %d, HESSE status %d, MIGRAD status %d)\n", fitResult->edm(), fitResult->covQual(), fitResult->status()/100, fitResult->status()%10);
    break; // ok with covariance matrix being accurate (3) or made pos-def (2)
   }
   // We know something is wrong with covariance matrix, need to find the problem
   qPar=0; worstPar=0;relParError=0, worstRelParError=0;
   printf("\nFit run %d: edm %f, cov qual %d, HESSE status %d, MIGRAD status %d, need to freeze a parameter\n", bombingRuns, fitResult->edm(), fitResult->covQual(), fitResult->status()/100, fitResult->status()%10);
   qFloatingPars = (RooArgSet*)totalPdf->getParameters(m)->selectByAttrib("Constant",kFALSE);
   // Start looking for questionable parameters
   parIter=qFloatingPars->createIterator(); qPar=(RooRealVar*)parIter->Next(); unsigned iPar=0;
   while(qPar) {
    relParError = qPar->getError() / qPar->getVal();
    if(
     strcmp("Nbkg",qPar->GetName())!=0 && // Cannot treat the background normalization as a questionable parameter
     relParError >= 0.5 && // If the parameter's error is 100% or more, it's questionable
     relParError>worstRelParError
    ) { worstRelParError=relParError; worstPar=qPar; }
    qPar=(RooRealVar*)parIter->Next(); iPar++;
   }
   if(worstPar) {
    printf("The worst parameter is %s , relative parameter error %f. Freezing it and retrying the fit...\n", worstPar->GetName(), worstRelParError);
    worstPar->setConstant();
   } else {
    printf("Couldn't find a questionable parameter, retry blindly\n");
   }
   bombingRuns++;
  }

   // Print fit results
  fitResult->printStream(fitResultFile,RooPrintable::kValue,RooPrintable::kVerbose);
  fitResultFile << endl;
  ios_base::fmtflags flags = fitResultFile.flags();
 
  // Combined errors and goodness-of-fit test
  TH1D *dataHistWithCombErrors=(TH1D*)dataHist->Clone("dataHistWithCombErrors"); dataHistWithCombErrors->SetDirectory(0);
  TH1D *ratioDataPdf=(TH1D*)dataHist->Clone("ratioDataPdf"); ratioDataPdf->SetDirectory(0); ratioDataPdf->Scale(0); ratioDataPdf->Clear();
  TH1D *pdfErrorBand=(TH1D*)dataHist->Clone("pdfErrorBand"); pdfErrorBand->SetDirectory(0); pdfErrorBand->Scale(0); pdfErrorBand->Clear();
  TH1D *ratioErrorBand=(TH1D*)dataHist->Clone("ratioErrorBand"); ratioErrorBand->SetDirectory(0); ratioErrorBand->Scale(0); ratioErrorBand->Clear();

  double chiSquaredSum=0;
  fitResultFile << "###########################################\n";
  fitResultFile << "# Computing chi square goodness-of-fit    #\n";
  fitResultFile << "###########################################\n";
  // Should we propagate the smearing PDF to the MC errors? not a big deal for now
  for(int nb=1; nb<=dataHistWithCombErrors->GetNbinsX(); nb++) {
    // Integrate the pdf
    double xmin=dataHistWithCombErrors->GetBinLowEdge(nb);
    double xmax=dataHistWithCombErrors->GetBinLowEdge(nb+1);
    if(xmin>=fitMassHi||xmax<=fitMassLo) continue;
    m.setVal((xmin+xmax)/2.);
    double modelVal = totalPdf->getVal(RooArgSet(m)) * Ntotal.getVal() * (xmax-xmin);
    double sigVal = sigFitterShape->theShape->getVal(RooArgSet(m)) * Nsig.getVal() * (xmax-xmin);
    //m.setRange("theIntegral",xmin,xmax);
    //double modelVal = totalPdf->createIntegral(m,m,"theIntegral")->getVal() * Ntotal.getVal();
    //double sigVal = sigFitterShape->theShape->createIntegral(m,m,"theIntegral")->getVal() * Nsig.getVal();
    double theError,theError2;
    double dataError = (dataHistWithCombErrors->GetBinContent(nb)>0)? dataHistWithCombErrors->GetBinError(nb):1;
    double dataPdfRatio = (dataHistWithCombErrors->GetBinContent(nb)>0 && modelVal>0)? dataHistWithCombErrors->GetBinContent(nb) / modelVal : 1;
    if(0!=sigFitterShape->templateHist) {
      int mcBin=templateHist->FindBin(dataHistWithCombErrors->GetBinCenter(nb));
      double absTemplateError = (templateHist->GetBinContent(mcBin)>0)? templateHist->GetBinError(mcBin)/templateHist->GetBinContent(mcBin) * sigVal : sigVal;
      double relTemplateError = absTemplateError / modelVal;
      theError2=
        pow(dataError,2)+
        pow(absTemplateError,2);
      theError=sqrt(theError2);
      ratioErrorBand->SetBinError(nb, dataPdfRatio*relTemplateError);
      pdfErrorBand->SetBinContent(nb, modelVal);
      pdfErrorBand->SetBinError(nb,absTemplateError);
    } else { theError = dataError; theError2 = theError*theError; }
    //dataHistWithCombErrors->SetBinError(nb,theError);
    double chiSquaredContribution=pow(dataHistWithCombErrors->GetBinContent(nb)-modelVal ,2)/theError2;
    if(chiSquaredContribution==chiSquaredContribution) chiSquaredSum+=chiSquaredContribution;
    fitResultFile << Form("\tbin %3d: pdf=%9.2f, data=%9d, error=%9.2f, contributing %8.4f to chi2\n",nb,modelVal,(int)dataHistWithCombErrors->GetBinContent(nb),theError,chiSquaredContribution);

    ratioDataPdf->SetBinContent(nb, dataPdfRatio);
    double ratioDataError = (dataHistWithCombErrors->GetBinContent(nb)>0)? dataHist->GetBinError(nb) / modelVal : 1;
    ratioDataPdf->SetBinError(nb,ratioDataError);
    ratioErrorBand->SetBinContent(nb, 1);
  }
  int ndf=dataHistWithCombErrors->FindBin(fitMassHi)-dataHistWithCombErrors->FindBin(fitMassLo)+1 -((RooArgSet*)totalPdf->getParameters(m)->selectByAttrib("Constant",kFALSE))->getSize();
  double fitProb=TMath::Prob(chiSquaredSum,ndf);
  printf("fit probability: %4.2f%%\n",100.*fitProb);
  printf("signal events  : %9.2f\n",Nsig.getVal());
  printf("bkg events     : %9.2f\n",Nbkg.getVal());

  // Draw fit
  TString rangeName=Form("fitRange_plotting");
  m.setRange(rangeName,fitMassLo,fitMassHi);
  gStyle->SetOptStat(0);
  TCanvas *canvas=new TCanvas(Form("c_%s",histName.c_str()), Form("c_%s",histName.c_str()),600,480);
  canvas->SetTopMargin(0.0); 
  canvas->SetBottomMargin(0); 
  canvas->SetRightMargin(0.02);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1.,1);
  pad1->SetTopMargin(0.1);
  pad1->SetLeftMargin(0.15);
  pad1->SetRightMargin(0.04);
  pad1->SetBottomMargin(0.03); 
  pad1->SetGridx();         
  pad1->Draw();             
  pad1->cd();              
  RooPlot *frame = m.frame(Bins(30));
  plotTitle.ReplaceAll("~"," ");
  frame->SetTitle(plotTitle!=""?plotTitle:histName);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetTitle("");
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetTitle(Form("Events / %.1f GeV",dataHist->GetBinWidth(1)));
  frame->GetYaxis()->SetLabelSize(0.05);
  //totalPdf->plotOn(frame, Name("sigPdfLine"), LineStyle(kSolid), LineColor(kViolet-4), DrawOption("l"), Normalization(Ntotal.getVal(),RooAbsReal::NumEvent) /*, NormRange("fitRange")*/);
  //totalPdf->plotOn(frame, Name("bkgPdfLine"), Components("bkgModel"),LineColor(kCyan-6),LineStyle(kDashed),DrawOption("l"),Normalization(Ntotal.getVal(),RooAbsReal::NumEvent) /*, NormRange("fitRange")*/);
  //totalPdf->plotOn(frame, Components("bkgModel"),VisualizeError(*fitResult,1,kTRUE),FillStyle(3004),FillColor(kCyan-6),Normalization(Ntotal.getVal(),RooAbsReal::NumEvent) /*, NormRange("fitRange")*/);
  totalPdf->plotOn(frame, Name("sigPdfLine"), LineStyle(kSolid), LineColor(kViolet-4), DrawOption("l"), Normalization(Ntotal.getVal(),RooAbsReal::NumEvent), NormRange(rangeName));
  totalPdf->plotOn(frame, Name("bkgPdfLine"), Components("bkgModel"),LineColor(kCyan-6),LineStyle(kDashed),DrawOption("l"),Normalization(Ntotal.getVal(),RooAbsReal::NumEvent), NormRange(rangeName));
  totalPdf->plotOn(frame, Components("bkgModel"),VisualizeError(*fitResult,1,kTRUE),FillStyle(3004),FillColor(kCyan-6),Normalization(Ntotal.getVal(),RooAbsReal::NumEvent), NormRange(rangeName));
  frame->GetYaxis()->SetRangeUser(0, 1.4*dataHist->GetBinContent(dataHist->GetMaximumBin()));
  frame->Draw();
  pdfErrorBand->SetMarkerColor(kViolet-4);
  pdfErrorBand->SetFillStyle(3254);
  pdfErrorBand->SetFillColor(kViolet-4);
  pdfErrorBand->SetLineColor(kBlack);
  pdfErrorBand->Draw("e2 same");
  dataHistWithCombErrors->SetMarkerStyle(20);
  dataHistWithCombErrors->SetMarkerSize(0.8);
  dataHistWithCombErrors->SetMarkerColor(kBlack);
  dataHistWithCombErrors->SetLineColor(kBlack);
  dataHistWithCombErrors->SetLineWidth(2);
  dataHistWithCombErrors->Draw("p e0 same");
  TLegend *legend=new TLegend(.68,.62,.95,.88);
  legend->SetFillColor(0);
  legend->AddEntry(dataHistWithCombErrors, "Data","lp");
  legend->AddEntry("sigPdfLine", sigFitterShape->plotLabel,"l");
  if(0!=sigFitterShape->templateHist) legend->AddEntry(pdfErrorBand, "MC Stat. Unc.","f");
  legend->AddEntry("bkgPdfLine", bkgFitterShape->plotLabel,"l");
  legend->Draw("same");
  //TPaveText *fitStatsPave = new TPaveText(0.68,0.3,0.94,.58, "NDC");
  TPaveText *fitStatsPave = new TPaveText(0.18,0.6,0.44,.88, "NDC");
  fitStatsPave->SetFillColor(0);
  fitStatsPave->SetBorderSize(1);
  fitStatsPave->SetTextFont(42);
  fitStatsPave->AddText(Form("%d data events",(int)Ndata));
  fitStatsPave->AddText(Form("N_{sig} %d #pm %d",(int)round(Nsig.getVal()),(int)round(sqrt(pow(Nbkg.getError(),2)+Ntotal.getVal()))));
  fitStatsPave->AddText(Form("N_{bkg} %d #pm %d",(int)round(Nbkg.getVal()),(int)round(Nbkg.getError())));
  fitStatsPave->AddText(Form("Fit prob. %3.1f%%", fitProb*100.));
  fitStatsPave->AddText(Form("EDM %3.1e", fitResult->edm()));
  fitStatsPave->Draw("SAME");
  canvas->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetLeftMargin(0.15);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.4);
  pad2->SetRightMargin(0.04);
  pad2->SetGridx(); 
  pad2->Draw();
  pad2->cd();
  ratioErrorBand->SetTitle("");
  ratioErrorBand->GetXaxis()->SetTitleSize(0.18);
  ratioErrorBand->GetXaxis()->SetTitleOffset(1.1);
  ratioErrorBand->GetYaxis()->SetNdivisions(503);
  ratioErrorBand->GetXaxis()->SetLabelSize(0.2);
  ratioErrorBand->GetXaxis()->SetTitle("Invariant mass [GeV]");
  ratioErrorBand->GetXaxis()->SetRangeUser(fitMassLo,fitMassHi-0.01);
  ratioErrorBand->GetYaxis()->SetTitle("Data/PDF");
  ratioErrorBand->GetYaxis()->CenterTitle();
  ratioErrorBand->GetYaxis()->SetTitleOffset(0.37);
  ratioErrorBand->GetYaxis()->SetTitleSize(0.12);
  ratioErrorBand->GetYaxis()->SetLabelSize(0.15);
  ratioErrorBand->SetFillStyle(3254);
  ratioErrorBand->SetFillColor(kBlack);
  ratioErrorBand->SetMinimum(.3);
  ratioErrorBand->SetMaximum(1.7);
  ratioErrorBand->Draw("E2");
  ratioDataPdf->SetLineColor(kBlack);
  ratioDataPdf->SetMarkerStyle(20);
  ratioDataPdf->SetMarkerSize(0.8);
  ratioDataPdf->Draw("P E0 x0 SAME");
  TLine *baseline = new TLine(fitMassLo,1,fitMassHi,1);
  baseline->SetLineStyle(kSolid); baseline->Draw("SAME");
  canvas->Print(Form("%s/plots/%s.png",outputDir.c_str(),histName.c_str()));
  canvas->Print(Form("%s/plots/%s.pdf",outputDir.c_str(),histName.c_str()));


  RooArgList parlist = fitResult->floatParsFinal();
  fitResultFile << "\n"; 
  fitResultFile << "###########################################\n";
  fitResultFile << "# Correlation Matrix                      #\n";
  fitResultFile << "###########################################\n";
  for(int i=0; i<parlist.getSize(); i++) {
    for(int j=0; j<parlist.getSize(); j++) 
      fitResultFile << "  " << setw(7) << setprecision(4) << fixed << fitResult->correlationMatrix()(i,j);    
    fitResultFile << endl;
  }
  fitResultFile.flags(flags);
  fitResultFile << "###########################################\n";
  fitResultFile << "# Fit Summary                             #\n";
  fitResultFile << "###########################################\n";
  fitResultFile << Form("summary_Nsig %.2f +/- %.2f\nsummary_EDM %.2e\nsummary_fitProb %.4f\n",Nsig.getVal(),sqrt(pow(Nbkg.getError(),2)+Ntotal.getVal()),fitResult->edm(),fitProb);
  fitResultFile.close();
  return fitResult;
}

