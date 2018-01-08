
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
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "fitterShape.h"
using namespace RooFit;
using namespace fitterShape;
const int fitMassLo=60;
const int fitMassHi=120;
RooFitResult *fitBin(
  string histName="pass_ptBin0_etaBin0",
  string dataFileName="",
  string sigTemplateFileName="",
  string bkgTemplateFileName="",
  string bkgTemplateFileName2="",
  shapeType signalModel=kTemplateConvDas,
  shapeType bkgModel=kBkgErfcExp, 
  TString plotTitle="",
  string outputDir="",
  TString signalLabel="",
  TString bkgLabel=""
  //metadata file name with kinematics info, initial params, info for plot 
) {
  //if(outputDir!="" && outputDir[outputDir.size()-1]!='/') outputDir=outputDir+"/";
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  system(Form("mkdir -p %s/plots",outputDir.c_str()));
  
  // Open output file for the fit results
  ofstream fitResultFile;
  char txtfname[1000];    
  sprintf(txtfname,"%s/plots/fitres_%s.txt",outputDir.c_str(),histName.c_str());
  fitResultFile.open(txtfname);
  assert(fitResultFile.is_open());
  // Open the file with data histograms and load the histogram we are fitting
  TFile *dataFile=TFile::Open(dataFileName.c_str(),"READ"); assert(dataFile && dataFile->IsOpen());
  TH1D *dataHist=(TH1D*)dataFile->Get(histName.c_str()); assert(dataHist); dataHist->SetDirectory(0); dataFile->Close();
  
  // If called for, open the files with signal and background templates and load the histograms to use as a template
  TFile *sigTemplateFile=0,*bkgTemplateFile=0, *bkgTemplateFile2=0;
  TH1D *sigTemplateHist=0,*bkgTemplateHist=0,*bkgTemplateHist2=0;
  bool useSigMCTemplate=(sigTemplateFileName!=""); if(useSigMCTemplate) {
    sigTemplateFile=TFile::Open(sigTemplateFileName.c_str(),"READ"); assert(sigTemplateFile && sigTemplateFile->IsOpen()); 
    sigTemplateHist=(TH1D*)sigTemplateFile->Get(histName.c_str()); assert(sigTemplateHist); sigTemplateHist->SetDirectory(0); sigTemplateFile->Close();
    for(unsigned ibin=0; ibin<=(unsigned)sigTemplateHist->GetNbinsX(); ibin++) if(sigTemplateHist->GetBinContent(ibin)<=0) sigTemplateHist->SetBinContent(ibin,0.0001);
  }  
  bool useBkgMCTemplate=(bkgTemplateFileName!=""); if(useBkgMCTemplate) {
    bkgTemplateFile=TFile::Open(bkgTemplateFileName.c_str(),"READ"); assert(bkgTemplateFile && bkgTemplateFile->IsOpen());
    bkgTemplateHist=(TH1D*)bkgTemplateFile->Get(histName.c_str()); assert(bkgTemplateHist); bkgTemplateHist->SetDirectory(0); bkgTemplateFile->Close();
    for(unsigned ibin=0; ibin<=(unsigned)bkgTemplateHist->GetNbinsX(); ibin++) if(bkgTemplateHist->GetBinContent(ibin)<=0) bkgTemplateHist->SetBinContent(ibin,0.0001);
    if(bkgTemplateFileName2!="") {
      bkgTemplateFile2 = TFile::Open(bkgTemplateFileName2.c_str(),"READ"); assert(bkgTemplateFile2 && bkgTemplateFile2->IsOpen());
      bkgTemplateHist2=(TH1D*)bkgTemplateFile2->Get(histName.c_str()); assert(bkgTemplateHist2); bkgTemplateHist2->SetDirectory(0); bkgTemplateFile2->Close();
      for(unsigned ibin=0; ibin<=(unsigned)bkgTemplateHist2->GetNbinsX(); ibin++) if(bkgTemplateHist2->GetBinContent(ibin)<=0) bkgTemplateHist2->SetBinContent(ibin,0.0001);
    }
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
  //RooDataHist *templateRDH = new RooDataHist("templateRDH","templateRDH",RooArgSet(m),sigTemplateHist);
  //RooHistPdf  *templateRHP = new RooHistPdf("templateRHP","templateRHP",m,*templateRDH,2);
  //RooAbsPdf *signalModelPdf=new RooFFTConvPdf("signalModel","signalModel",m,*templateRHP,*cb);
  
  // Signal model : das distribution
  fitterShape::fitterShapeBase *sigFitterShape=0, *bkgFitterShape=0;
 
  signalLabel.ReplaceAll("~"," ");
  bkgLabel.ReplaceAll("~"," ");
  // if statement
  switch(signalModel) {
    case kTemplateConvDas    : sigFitterShape=new fitterShape::templateConvDas(m,sigTemplateHist); break;
    case kTemplateConvGaus   : sigFitterShape=new fitterShape::templateConvGaus(m,sigTemplateHist); break;
    case kTemplateConvDoubleGaus   : sigFitterShape=new fitterShape::templateConvDoubleGaus(m,sigTemplateHist); break;
    case kTemplateBWConvGaus : sigFitterShape=new fitterShape::templateBWConvGaus(m,sigTemplateHist); break;
    case kBreitWignerConvDas : sigFitterShape=new fitterShape::breitWignerConvDas(m); break;
    case kTemplateConvLandau   : sigFitterShape=new fitterShape::templateConvLandau(m,sigTemplateHist); break;
    case kTemplateConvCrystalBall   : sigFitterShape=new fitterShape::templateConvCrystalBall(m,sigTemplateHist); break;
    case kTemplateConvBeta   : sigFitterShape=new fitterShape::templateConvBeta(m,sigTemplateHist); break;

    default                  : printf("Unsupported signal model\n"); assert(0); break;
  }
  switch(bkgModel) {
    case kBkgErfcExp         : bkgFitterShape=new fitterShape::bkgErfcExp(m); break;
    case kBkgErfcExpPlusExp  : bkgFitterShape=new fitterShape::bkgErfcExpPlusExp(m); break;
    case kBkgDasPlusExp      : bkgFitterShape=new fitterShape::bkgDasPlusExp(m); break;
    case kBkgDas             : bkgFitterShape=new fitterShape::bkgDas(m); break;
    case kBkgDoubleExp       : bkgFitterShape=new fitterShape::bkgDoubleExp(m); break;
    case kBkgTemplate        : bkgFitterShape=new fitterShape::bkgTemplate(m,bkgTemplateHist,bkgLabel); break;
    case kBkgTemplatePlusExp : bkgFitterShape=new fitterShape::bkgTemplatePlusExp(m,bkgTemplateHist,bkgLabel); break;
    case kBkgTwoTemplates    : bkgFitterShape=new fitterShape::bkgTwoTemplates(m,bkgTemplateHist,bkgTemplateHist2,bkgLabel); break;
    default                  : printf("Unsupported background model\n"); assert(0); break;
  }
  //sigFitterShape->initializeParams(map);
  //bkgFitterShape->initializeParams(map)

  // Signal+background model
  double Ndata = dataHist->Integral( dataHist->FindBin((double)fitMassLo), dataHist->FindBin((double)fitMassHi-0.01) );
  RooRealVar Ntotal = RooRealVar("Ntotal","Ntotal",Ndata,0,Ndata); Ntotal.setConstant(kTRUE);
  RooRealVar Nbkg = RooRealVar("Nbkg","Nbkg",0.2*Ndata, 0, .9*Ndata);
  
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
  if(true) { //if(!(bkgFitterShape->templateHist!=0)) {
   printf("###########################################\n");
   printf("# Performing peak estimation procedure    #\n");
   printf("###########################################\n");
   TF1 *peakEstimator = new TF1("peakEstimator","gausn",fitMassLo,fitMassHi);
   peakEstimator->SetParameter(0,.9*Ndata); peakEstimator->SetParLimits(0, 1,Ndata);
   peakEstimator->SetParameter(1,91); peakEstimator->SetParLimits(1, 86.1876, 96.1876);
   peakEstimator->SetParameter(0, 1); peakEstimator->SetParLimits(2, 0.001,10);
   if(signalModel==kTemplateConvDas||signalModel==kTemplateConvGaus||signalModel==kTemplateConvDoubleGaus||signalModel==kTemplateBWConvGaus||signalModel==kBreitWignerConvDas||signalModel==kTemplateConvLandau||signalModel==kTemplateConvCrystalBall||signalModel==kTemplateConvBeta) {
    //TODO: estimate this with a cheap gaussian or pol2 fit in small zmass window
    //double peakPos=dataHist->GetBinCenter(dataHist->GetMaximumBin());
    //double mcPeakPos=sigTemplateHist->GetBinCenter(dataHist->GetMaximumBin());
    printf("\nPerforming data peak estimation\n\n");
    dataHist->Fit(peakEstimator, "MN0", "goff", 86.2, 96.2);
    double peakPos=peakEstimator->GetParameter(1);
    double peakSigma=peakEstimator->GetParameter(2);
    double mcPeakPos, mcPeakSigma;
    if(signalModel==kTemplateConvDas||signalModel==kTemplateConvGaus||signalModel==kTemplateConvDoubleGaus||signalModel==kTemplateBWConvGaus||signalModel==kTemplateConvLandau||signalModel==kTemplateConvCrystalBall||signalModel==kTemplateConvBeta) {
     printf("\nPerforming MC peak estimation\n\n");
     sigTemplateHist->Fit(peakEstimator, "MN0", "goff", 86.2, 96.2);
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
      ((fitterShape::templateConvDas*)sigFitterShape)->sigma->setRange(0.5*sigmaEst,1.5*sigmaEst);
      ((fitterShape::templateConvDas*)sigFitterShape)->sigma->setVal(sigmaEst);
     } else if(signalModel==kTemplateConvGaus) {
      ((fitterShape::templateConvGaus*)sigFitterShape)->mean->setRange(meanEst-2.,meanEst+2.);
      ((fitterShape::templateConvGaus*)sigFitterShape)->mean->setVal(peakPos-mcPeakPos);
      ((fitterShape::templateConvGaus*)sigFitterShape)->sigma->setRange(0.5*sigmaEst,1.5*sigmaEst);
      ((fitterShape::templateConvGaus*)sigFitterShape)->sigma->setVal(sigmaEst);
     } 
     else if(signalModel==kTemplateConvDoubleGaus) {
       ((fitterShape::templateConvDoubleGaus*)sigFitterShape)->mean1->setRange(meanEst-2.,meanEst+2.);
       ((fitterShape::templateConvDoubleGaus*)sigFitterShape)->mean1->setVal(peakPos-mcPeakPos);
       ((fitterShape::templateConvDoubleGaus*)sigFitterShape)->sigma1->setRange(0.5*sigmaEst,1.5*sigmaEst);
       ((fitterShape::templateConvDoubleGaus*)sigFitterShape)->sigma1->setVal(sigmaEst);
       ((fitterShape::templateConvDoubleGaus*)sigFitterShape)->mean2->setRange(meanEst-3,meanEst+3);
       ((fitterShape::templateConvDoubleGaus*)sigFitterShape)->mean2->setVal(peakPos-mcPeakPos-1);
       ((fitterShape::templateConvDoubleGaus*)sigFitterShape)->sigma2->setRange(0.5*sigmaEst,1.5*sigmaEst);
       ((fitterShape::templateConvDoubleGaus*)sigFitterShape)->sigma2->setVal(sigmaEst);

     }else if(signalModel==kTemplateBWConvGaus) {
      ((fitterShape::templateBWConvGaus*)sigFitterShape)->mean->setRange(meanEst-2.,meanEst+2.);
      ((fitterShape::templateBWConvGaus*)sigFitterShape)->mean->setVal(meanEst);
      ((fitterShape::templateBWConvGaus*)sigFitterShape)->sigma->setRange(0.5*sigmaEst,1.5*sigmaEst);
      ((fitterShape::templateBWConvGaus*)sigFitterShape)->sigma->setVal(sigmaEst);
     } else if(signalModel==kBreitWignerConvDas) {
      ((fitterShape::breitWignerConvDas*)sigFitterShape)->mean->setRange(meanEst-2.,meanEst+2.);
      ((fitterShape::breitWignerConvDas*)sigFitterShape)->mean->setVal(meanEst);
      ((fitterShape::breitWignerConvDas*)sigFitterShape)->sigma->setRange(0.5*sigmaEst,1.5*sigmaEst);
      ((fitterShape::breitWignerConvDas*)sigFitterShape)->sigma->setVal(sigmaEst);
     } else if(signalModel==kTemplateConvLandau) {
       ((fitterShape::templateConvLandau*)sigFitterShape)->mean->setRange(meanEst-2.,meanEst+2.);
       ((fitterShape::templateConvLandau*)sigFitterShape)->mean->setVal(peakPos-mcPeakPos);
       ((fitterShape::templateConvLandau*)sigFitterShape)->sigma->setRange(0.5*sigmaEst,1.5*sigmaEst);
       ((fitterShape::templateConvLandau*)sigFitterShape)->sigma->setVal(sigmaEst);
     }else if(signalModel==kTemplateConvCrystalBall) {
       ((fitterShape::templateConvCrystalBall*)sigFitterShape)->mean->setRange(meanEst-2.,meanEst+2.);
       ((fitterShape::templateConvCrystalBall*)sigFitterShape)->mean->setVal(peakPos-mcPeakPos);
       ((fitterShape::templateConvCrystalBall*)sigFitterShape)->sigma->setRange(0.5*sigmaEst,1.5*sigmaEst);
       ((fitterShape::templateConvCrystalBall*)sigFitterShape)->sigma->setVal(sigmaEst);
       ((fitterShape::templateConvCrystalBall*)sigFitterShape)->a->setRange(meanEst-2.,meanEst+2.);
       ((fitterShape::templateConvCrystalBall*)sigFitterShape)->a->setVal(peakPos-mcPeakPos);
       ((fitterShape::templateConvCrystalBall*)sigFitterShape)->n->setRange(0.5*sigmaEst,1.5*sigmaEst);
       ((fitterShape::templateConvCrystalBall*)sigFitterShape)->n->setVal(sigmaEst);
     }else if(signalModel==kTemplateConvBeta) {
       ((fitterShape::templateConvBeta*)sigFitterShape)->a->setRange(meanEst-2.,meanEst+2.);
       ((fitterShape::templateConvBeta*)sigFitterShape)->a->setVal(peakPos-mcPeakPos);
       ((fitterShape::templateConvBeta*)sigFitterShape)->b->setRange(0.5*sigmaEst,1.5*sigmaEst);
       ((fitterShape::templateConvBeta*)sigFitterShape)->b->setVal(sigmaEst);
       }
   }}
   //delete bkgEstimator;
   delete peakEstimator;
  }  
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
  if(bkgFitterShape->templateHist==0) {
   printf("###########################################\n");
   printf("# Low mass fit with dynamic freezing      #\n");
   printf("###########################################\n");
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
  }

  // Perform the actual fit
  printf("###########################################\n");
  printf("# Real fit with dynamic freezing          #\n");
  printf("###########################################\n");
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
  TH1D *totalPdfErrorBand=(TH1D*)dataHist->Clone("totalPdfErrorBand"); totalPdfErrorBand->SetDirectory(0); totalPdfErrorBand->Scale(0); totalPdfErrorBand->Clear();
  TH1D *sigPdfErrorBand=(TH1D*)dataHist->Clone("sigPdfErrorBand"); sigPdfErrorBand->SetDirectory(0); sigPdfErrorBand->Scale(0); sigPdfErrorBand->Clear();
  TH1D *bkgPdfErrorBand=(TH1D*)dataHist->Clone("bkgPdfErrorBand"); bkgPdfErrorBand->SetDirectory(0); bkgPdfErrorBand->Scale(0); bkgPdfErrorBand->Clear();
  TH1D *ratioErrorBand=(TH1D*)dataHist->Clone("ratioErrorBand"); ratioErrorBand->SetDirectory(0); ratioErrorBand->Scale(0); ratioErrorBand->Clear();

  double chiSquaredSum=0;
  fitResultFile << "###########################################\n";
  fitResultFile << "# Computing chi square goodness-of-fit    #\n";
  fitResultFile << "###########################################\n";
  
  // Create the resolution function (x) error - maybe this should be in fitterShape ?
  // Optimize this at some point, to do list
  TH1D *sigTemplateError=0, *bkgTemplateError=0, *bkgTemplateError2=0;
  m.setRange(-100,100);
  if(sigFitterShape->templateHist!=0 && sigFitterShape->resolutionFunction!=0) {
    sigTemplateError=(TH1D*)sigTemplateHist->Clone("sigTemplateError"); sigTemplateError->SetDirectory(0);
    sigTemplateError->Reset(); sigTemplateError->Clear();
    // calculate weighted quadrature sum of the errors from all the bins, affecting bin #nb
    for(int nb=1; nb<=sigTemplateHist->GetNbinsX(); nb++) {
      double sumWeightsSquared=0; double sumWeightedErrorsSquared=0;
      for(int mb=1; mb<=sigTemplateHist->GetNbinsX(); mb++) {
        m.setRange("int",sigTemplateHist->GetBinCenter(nb)-sigTemplateHist->GetBinLowEdge(mb+1),sigTemplateHist->GetBinCenter(nb)-sigTemplateHist->GetBinLowEdge(mb));
        double resFrac=sigFitterShape->resolutionFunction->createIntegral(m,NormSet(m),Range("int"))->getVal();
        sumWeightsSquared+=resFrac*resFrac; sumWeightedErrorsSquared+=pow(resFrac*sigTemplateHist->GetBinError(mb),2);
      }
      sigTemplateError->SetBinContent(nb, sqrt(sumWeightedErrorsSquared/sumWeightsSquared));
    }
    m.removeRange("int");
  }
  if(bkgFitterShape->templateHist!=0 && bkgFitterShape->resolutionFunction!=0) {
    bkgTemplateError=(TH1D*)bkgTemplateHist->Clone("bkgTemplateError"); bkgTemplateError->SetDirectory(0);
    bkgTemplateError->Reset(); bkgTemplateError->Clear();
    // calculate weighted quadrature sum of the errors from all the bins, affecting bin #nb
    for(int nb=1; nb<=bkgTemplateHist->GetNbinsX(); nb++) {
      double sumWeightsSquared=0; double sumWeightedErrorsSquared=0;
      for(int mb=1; mb<=bkgTemplateHist->GetNbinsX(); mb++) {
        m.setRange("int",bkgTemplateHist->GetBinCenter(nb)-bkgTemplateHist->GetBinLowEdge(mb+1),bkgTemplateHist->GetBinCenter(nb)-bkgTemplateHist->GetBinLowEdge(mb));
        double resFrac=bkgFitterShape->resolutionFunction->createIntegral(m,NormSet(m),Range("int"))->getVal();
        sumWeightsSquared+=resFrac*resFrac; sumWeightedErrorsSquared+=pow(resFrac*bkgTemplateHist->GetBinError(mb),2);
        //printf("\tfrom bin %d resFrac %f error %f\n",mb,resFrac,bkgTemplateHist->GetBinError(mb));
      }
      if(bkgModel==kBkgTwoTemplates)
        bkgTemplateError->SetBinContent(nb, sqrt(sumWeightedErrorsSquared/sumWeightsSquared)*((fitterShape::bkgTwoTemplates*)bkgFitterShape)->frac->getVal());
      else bkgTemplateError->SetBinContent(nb, sqrt(sumWeightedErrorsSquared/sumWeightsSquared));
    }
    m.removeRange("int");
  }
  if(bkgFitterShape->templateHist!=0 && bkgFitterShape->resolutionFunction!=0) {
    bkgTemplateError2=(TH1D*)bkgTemplateHist2->Clone("bkgTemplateError2"); bkgTemplateError2->SetDirectory(0);
    bkgTemplateError2->Reset(); bkgTemplateError2->Clear();
    // calculate weighted quadrature sum of the errors from all the bins, affecting bin #nb
    for(int nb=1; nb<=bkgTemplateHist2->GetNbinsX(); nb++) {
      double sumWeightsSquared=0; double sumWeightedErrorsSquared=0;
      for(int mb=1; mb<=bkgTemplateHist2->GetNbinsX(); mb++) {
        m.setRange("int",bkgTemplateHist2->GetBinCenter(nb)-bkgTemplateHist2->GetBinLowEdge(mb+1),bkgTemplateHist2->GetBinCenter(nb)-bkgTemplateHist2->GetBinLowEdge(mb));
        double resFrac=bkgFitterShape->resolutionFunction->createIntegral(m,NormSet(m),Range("int"))->getVal();
        sumWeightsSquared+=resFrac*resFrac; sumWeightedErrorsSquared+=pow(resFrac*bkgTemplateHist2->GetBinError(mb),2);
        //printf("\tfrom bin %d resFrac %f error %f\n",mb,resFrac,bkgTemplateHist2->GetBinError(mb));
      }
      if(bkgModel==kBkgTwoTemplates)
        bkgTemplateError2->SetBinContent(nb, sqrt(sumWeightedErrorsSquared/sumWeightsSquared)*(1.-((fitterShape::bkgTwoTemplates*)bkgFitterShape)->frac->getVal()));
      else bkgTemplateError2->SetBinContent(nb, sqrt(sumWeightedErrorsSquared/sumWeightsSquared));
    }
    m.removeRange("int");
  }
    
  m.setRange(fitMassLo,fitMassHi);
  for(int nb=1; nb<=dataHistWithCombErrors->GetNbinsX(); nb++) {
    double xmin=dataHistWithCombErrors->GetBinLowEdge(nb);
    double xmax=dataHistWithCombErrors->GetBinLowEdge(nb+1);
    if(xmin>=fitMassHi||xmax<=fitMassLo) continue;
    m.setVal((xmin+xmax)/2.);
    double modelVal = totalPdf->getVal(RooArgSet(m)) * Ntotal.getVal() * (xmax-xmin); // Determine the combined pdf value at the bin center
    double bkgVal = bkgFitterShape->theShape->getVal(RooArgSet(m)) * Nbkg.getVal() * (xmax-xmin); // Determine the background pdf value at the bin center
    double sigVal = sigFitterShape->theShape->getVal(RooArgSet(m)) * Nsig.getVal() * (xmax-xmin); // Determine the signal pdf value at the bin center
    double theError2=0; // sum of squares of statistical errors from data, and signal and bkg MC
    
    // compute the data error and add it to the sum of squares of errors
    double dataError = (dataHistWithCombErrors->GetBinContent(nb)>0)? dataHistWithCombErrors->GetBinError(nb):1;
    theError2+=dataError*dataError;    

    // Now compute the absolute and relative stat errors on the templates, adding to sum of squares of errors
    double relTemplateError=0;
    if(0!=sigFitterShape->templateHist) { // compute stat uncertainty from signal template
      int mcBin=sigTemplateHist->FindBin(dataHistWithCombErrors->GetBinCenter(nb));
      double absTemplateError;
      if(sigTemplateError) {
        absTemplateError=sigTemplateError->GetBinContent(mcBin)/sigTemplateHist->GetBinContent(mcBin) * sigVal;
        //absTemplateError = sigTemplateErrorConvRes->getVal(RooArgSet(m)) * sigTemplateError->Integral() * sigVal/sigTemplateHist->GetBinContent(mcBin);
      } else absTemplateError = (sigTemplateHist->GetBinContent(mcBin)>0)? sigTemplateHist->GetBinError(mcBin)/sigTemplateHist->GetBinContent(mcBin) * sigVal : sigVal;
      theError2+=absTemplateError*absTemplateError;
      relTemplateError += pow(absTemplateError,2);
      sigPdfErrorBand->SetBinContent(nb, sigVal);
      sigPdfErrorBand->SetBinError(nb,absTemplateError);
      totalPdfErrorBand->SetBinContent(nb, modelVal);
      totalPdfErrorBand->SetBinError(nb,absTemplateError);
    }
    if(0!=bkgFitterShape->templateHist) { // compute stat uncertainty from bkgnal template
      int mcBin=bkgTemplateHist->FindBin(dataHistWithCombErrors->GetBinCenter(nb));
      double absTemplateError;
      if(bkgTemplateError) {
        absTemplateError=bkgTemplateError->GetBinContent(mcBin)/bkgTemplateHist->GetBinContent(mcBin) * bkgVal;
        //absTemplateError = bkgTemplateErrorConvRes->getVal(RooArgSet(m)) * bkgTemplateError->Integral() * bkgVal/bkgTemplateHist->GetBinContent(mcBin);
      } else absTemplateError = (bkgTemplateHist->GetBinContent(mcBin)>0)? bkgTemplateHist->GetBinError(mcBin)/bkgTemplateHist->GetBinContent(mcBin) * bkgVal : bkgVal;
      theError2+=absTemplateError*absTemplateError;
      relTemplateError += pow(absTemplateError,2);
      bkgPdfErrorBand->SetBinContent(nb, bkgVal);
      bkgPdfErrorBand->SetBinError(nb,absTemplateError);
      totalPdfErrorBand->SetBinContent(nb, modelVal);
      totalPdfErrorBand->SetBinError(nb, sqrt(pow(totalPdfErrorBand->GetBinError(nb),2)+absTemplateError*absTemplateError));
    }
    if(0!=bkgFitterShape->templateHist2) { // compute stat uncertainty from bkgnal template
      int mcBin=bkgTemplateHist2->FindBin(dataHistWithCombErrors->GetBinCenter(nb));
      double absTemplateError;
      if(bkgTemplateError2) {
        absTemplateError=bkgTemplateError2->GetBinContent(mcBin)/bkgTemplateHist2->GetBinContent(mcBin) * bkgVal;
      } else absTemplateError = (bkgTemplateHist2->GetBinContent(mcBin)>0)? bkgTemplateHist2->GetBinError(mcBin)/bkgTemplateHist2->GetBinContent(mcBin) * bkgVal : bkgVal;
      theError2+=absTemplateError*absTemplateError;
      relTemplateError += pow(absTemplateError,2);
      bkgPdfErrorBand->SetBinContent(nb, bkgVal);
      bkgPdfErrorBand->SetBinError(nb,sqrt(pow(bkgPdfErrorBand->GetBinError(nb),2)+absTemplateError*absTemplateError));
      totalPdfErrorBand->SetBinContent(nb, modelVal);
      totalPdfErrorBand->SetBinError(nb, sqrt(pow(totalPdfErrorBand->GetBinError(nb),2)+absTemplateError*absTemplateError));
    }
    double theError = sqrt(theError2);
    relTemplateError = sqrt(relTemplateError) / modelVal; // Computation of relTemplateError complete
    double dataPdfRatio = (dataHistWithCombErrors->GetBinContent(nb)>0 && modelVal>0)? dataHistWithCombErrors->GetBinContent(nb) / modelVal : 1;
    ratioErrorBand->SetBinError(nb, dataPdfRatio*relTemplateError);
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
  dataRDH->plotOn(frame, Name("Data"), MarkerStyle(kFullCircle),MarkerSize(0.8),LineWidth(2),LineColor(kBlack),DrawOption("ZP"));
  totalPdf->plotOn(frame, Name("totalPdfLine"), LineStyle(kSolid), LineColor(kViolet-4), DrawOption("l"),Normalization(Ntotal.getVal(),RooAbsReal::NumEvent));
  //if(!(sigFitterShape->templateHist!=0)) totalPdf->plotOn(frame, Components("signalModel"),VisualizeError(*fitResult,1,kTRUE),FillStyle(3254),FillColor(kPink-5),Normalization(Ntotal.getVal(),RooAbsReal::NumEvent));
  //totalPdf->plotOn(frame, Name("sigPdfLine"), Components("signalModel"),LineColor(kPink-5),LineStyle(9),DrawOption("l"),Normalization(Ntotal.getVal(),RooAbsReal::NumEvent));
  totalPdf->plotOn(frame, Name("bkgPdfLine"), Components("bkgModel"),LineColor(kCyan-6),LineStyle(kDashed),DrawOption("l"),Normalization(Ntotal.getVal(),RooAbsReal::NumEvent));
  if(!(bkgFitterShape->templateHist!=0)) totalPdf->plotOn(frame, Components("bkgModel"),VisualizeError(*fitResult,1,kTRUE),FillStyle(3245),FillColor(kCyan-6),Normalization(Ntotal.getVal(),RooAbsReal::NumEvent));
  frame->SetMinimum(0);
  frame->SetMaximum(dataHist->GetBinContent(dataHist->GetMaximumBin())*1.5);
  frame->Draw();
  //if(0!=sigFitterShape->templateHist) {
  //  sigPdfErrorBand->SetMarkerColor(kPink-5);
  //  sigPdfErrorBand->SetFillStyle(3253);
  //  sigPdfErrorBand->SetFillColor(kPink-5);
  //  sigPdfErrorBand->SetLineColor(kBlack);
  //  sigPdfErrorBand->Draw("e2 same");
  //}
  if(0!=bkgFitterShape->templateHist) {
    bkgPdfErrorBand->SetMarkerColor(kCyan-6);
    bkgPdfErrorBand->SetFillStyle(3245);
    bkgPdfErrorBand->SetFillColor(kCyan-6);
    bkgPdfErrorBand->SetLineColor(kBlack);
    bkgPdfErrorBand->Draw("e2 same");
  }
  if(sigFitterShape->templateHist!=0 || bkgFitterShape->templateHist!=0) {
    totalPdfErrorBand->SetMarkerColor(kViolet-4);
    totalPdfErrorBand->SetFillStyle(3254);
    totalPdfErrorBand->SetFillColor(kViolet-4);
    totalPdfErrorBand->SetLineColor(kBlack);
    totalPdfErrorBand->Draw("e2 same");
  }
  TLegend *legend=new TLegend(.68,.62,.95,.88);
  legend->SetFillColor(0);
  legend->AddEntry("Data", "Data","lp");
  //legend->AddEntry("totalPdfLine", "Total PDF","l");
  //if(sigFitterShape->templateHist!=0 || bkgFitterShape->templateHist!=0) legend->AddEntry(totalPdfErrorBand, "PDF Stat. Unc.","f");
  //legend->AddEntry("sigPdfLine", sigFitterShape->plotLabel,"l");
  //if(0!=sigFitterShape->templateHist) legend->AddEntry(sigPdfErrorBand, "Sig. Stat. Unc.","f");
  legend->AddEntry("totalPdfLine", signalLabel!=""?signalLabel:sigFitterShape->plotLabel,"l");
  if(0!=sigFitterShape->templateHist || 0!=sigFitterShape->templateHist) legend->AddEntry(totalPdfErrorBand, "PDF Stat. Unc.","f");
  legend->AddEntry("bkgPdfLine", bkgFitterShape->plotLabel,"l");
  if(0!=sigFitterShape->templateHist) legend->AddEntry(bkgPdfErrorBand, "Bkg. Stat. Unc.","f");
  legend->Draw("same");
  //TPaveText *fitStatsPave = new TPaveText(0.68,0.3,0.94,.58, "NDC");
  TPaveText *fitStatsPave = new TPaveText(0.20,0.6,0.46,.88, "NDC");
  fitStatsPave->SetFillColor(0);
  fitStatsPave->SetBorderSize(1);
  fitStatsPave->SetTextFont(42);
  fitStatsPave->AddText(Form("%d data events",(int)Ndata));
  fitStatsPave->AddText(Form("N_{sig} %d #pm %d",(int)round(Nsig.getVal()),(int)round(sqrt(pow(Nbkg.getError(),2)+Ntotal.getVal()))));
  fitStatsPave->AddText(Form("N_{bkg} %d #pm %d",(int)round(Nbkg.getVal()),(int)round(Nbkg.getError())));
  fitStatsPave->AddText(Form("Fit prob. %3.1f%%", fitProb*100.));
  fitStatsPave->AddText(Form("EDM %3.1e", fitResult->edm()));
  fitStatsPave->SetTextAlign(12);
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
  canvas->Print(Form("%s/plots/%s.pdf",outputDir.c_str(),histName.c_str()));
  system(Form("gs -sDEVICE=png16m -dTextAlphaBits=4 -g1800x1440 -dUseCropBox -dFIXEDMEDIA -dPDFFitPage -o %s/plots/%s.png %s/plots/%s.pdf >/dev/null 2>&1",outputDir.c_str(),histName.c_str(),outputDir.c_str(),histName.c_str()));
  //canvas->Print(Form("%s/plots/%s.pdf",outputDir.c_str(),histName.c_str()));


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

