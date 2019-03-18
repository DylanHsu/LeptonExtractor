#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include "TMath.h"
#include "TH1D.h"
#include "RooFit.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooCMSShape.cc"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooFormulaVar.h"
#include "RooGaussDoubleSidedExp.cc"
#include "RooGaussian.h"
#include "RooGlobalFunc.h"
#include "RooLandau.h"
#include "RooTFnBinding.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
//#include "RooKeysPdf.h"
#include "RooRealVar.h"
#include "RooVoigtian.h"
#include "RooVoigtianShape.h"
#include "RooPlot.h"
#include "TF1.h"
#include "Math/DistFunc.h"
#include "RooTFnBinding.h" 
#include "RooCFunction1Binding.h"
#include "RooCFunction3Binding.h"

using namespace RooFit;
namespace fitterShape {
 // Parent and child class declarations
 enum shapeType {
  kTemplateConvDas,
  kTemplateConvGaus,
  kTemplateConvDoubleGaus,
  kTemplateBWConvGaus,
  kTemplateConvLandau,
  kTemplateConvBeta,
  kTemplateConvCrystalBall,
  kBreitWignerConvDas,
  kBkgErfcExp,
  kBkgErfcExpPlusExp,
  kBkgDasPlusExp,
  kBkgDas,
  kBkgDoubleExp,
  kBkgTemplate,
  kBkgTemplatePlusExp,
  kBkgTwoTemplates,
  nShapeTypes
 };
 class fitterShapeBase {
  public:
   fitterShapeBase():theShape(0){}
   virtual ~fitterShapeBase() { delete theShape; }
   RooAbsPdf *theShape=0;
   RooAbsPdf *resolutionFunction=0;
   TString plotLabel;
   TH1D *templateHist=0;
   TH1D *templateHist2=0;// hack, make it better later
   void initializeParams(std::map<std::string, double> &initialization);
 };
 class templateConvDas: public fitterShapeBase {
  public:
   templateConvDas(RooRealVar &m, TH1D *templateHist_); ~templateConvDas();
   RooRealVar *mean,*sigma,*kLo,*kHi;
   RooGaussDoubleSidedExp *dd;
   RooDataHist *templateRDH; RooHistPdf *templateRHP;
 };
 class templateConvGaus: public fitterShapeBase {
  public:
   templateConvGaus(RooRealVar &m, TH1D *templateHist_); ~templateConvGaus();
   RooRealVar *mean,*sigma;
   RooGaussian *gaus;
   RooDataHist *templateRDH; RooHistPdf *templateRHP;
 };
 class templateConvDoubleGaus: public fitterShapeBase {
 public:
   templateConvDoubleGaus(RooRealVar &m, TH1D *templateHist_); ~templateConvDoubleGaus();
   RooRealVar *mean1,*sigma1, *mean2, *sigma2, *frac;
   RooGaussian *gaus1, *gaus2;
   RooAbsPdf *sum, *theShape1, *theShape2;
   RooDataHist *templateRDH; RooHistPdf *templateRHP;
 };

 class templateBWConvGaus: public fitterShapeBase {
  public:
   templateBWConvGaus(RooRealVar &m, TH1D *templateHist_); ~templateBWConvGaus();
   RooRealVar *mean,*sigma,*mass,*width,*frac;
   RooGaussian *gaus;
   RooDataHist *templateRDH; RooHistPdf *templateRHP;
   RooBreitWigner *bw;
   RooAbsPdf *templatePlusBw;
 };
 class templateConvLandau: public fitterShapeBase {
 public:
   templateConvLandau(RooRealVar &m, TH1D *templateHist_); ~templateConvLandau();
   RooRealVar *mean,*sigma;
   RooFormulaVar *x;
   RooLandau *landau;
   RooDataHist *templateRDH; RooHistPdf *templateRHP;
 };
 class templateConvCrystalBall: public fitterShapeBase {
 public:
   templateConvCrystalBall(RooRealVar &m, TH1D *templateHist_); ~templateConvCrystalBall();
   RooRealVar *mean, *sigma, *a, *n;
   RooCBShape *CB;
   RooDataHist *templateRDH; RooHistPdf *templateRHP;
 };

 class templateConvBeta: public fitterShapeBase {
 public:
   templateConvBeta(RooRealVar &m, TH1D *templateHist_); ~templateConvBeta();
   RooRealVar *a,*b;
   RooAbsPdf *beta;
   RooDataHist *templateRDH; RooHistPdf *templateRHP;
 };

 class breitWignerConvDas: public fitterShapeBase {
  public:
   breitWignerConvDas(RooRealVar &m, TH1D *breitWignerHist=0); ~breitWignerConvDas();
   RooRealVar *mean,*sigma,*kLo,*kHi,*mass,*width;
   RooGaussDoubleSidedExp *dd;
   RooBreitWigner *bw;
 };
class bkgErfcExp: public fitterShapeBase {
  public:
   bkgErfcExp(RooRealVar &m, TH1D *templateHist_=0); ~bkgErfcExp();
   RooRealVar *mpv,*beta,*gamma,*zmass;
 };
class bkgErfcExpPlusExp: public fitterShapeBase {
  public:
   bkgErfcExpPlusExp(RooRealVar &m, TH1D *templateHist_=0); ~bkgErfcExpPlusExp();
   RooRealVar *mpv,*beta,*gamma,*zmass, *t1, *frac;
   RooErfcExpoPdf *erfcExp;
   RooExponential *exp1;
 };
 class bkgDasPlusExp: public fitterShapeBase {
  public:
   bkgDasPlusExp(RooRealVar &m, TH1D *breitWignerHist=0); ~bkgDasPlusExp();
   RooRealVar *mean,*sigma,*kLo,*kHi,*t1, *frac;
   RooGaussDoubleSidedExp *dd;
   RooExponential *exp1;
 };
 class bkgDas: public fitterShapeBase {
  public:
   bkgDas(RooRealVar &m, TH1D *breitWignerHist=0); ~bkgDas();
   RooRealVar *mean,*sigma,*kLo,*kHi;
 };
 class bkgDoubleExp: public fitterShapeBase {
  public:
   bkgDoubleExp(RooRealVar &m, TH1D *templateHist_=0); ~bkgDoubleExp();
   RooRealVar *t1,*t2,*frac;
   RooExponential *exp1, *exp2;
 };
 class bkgTemplate: public fitterShapeBase {
  public:
   bkgTemplate(RooRealVar &m, TH1D *templateHist_, TString plotLabel_);
   ~bkgTemplate();
   RooDataHist *templateRDH;
 };
 class bkgTemplatePlusExp: public fitterShapeBase {
  public:
   bkgTemplatePlusExp(RooRealVar &m, TH1D *templateHist_, TString plotLabel_);
   ~bkgTemplatePlusExp();
   RooRealVar *t1,*frac;
   RooExponential *exp1;
   RooDataHist *templateRDH; RooHistPdf *templateRHP;
 };
 class bkgTwoTemplates: public fitterShapeBase {
  public:
   bkgTwoTemplates(RooRealVar &m, TH1D *templateHist_, TH1D *templateHist2_, TString plotLabel_);
   ~bkgTwoTemplates();
   RooRealVar *frac;
   RooDataHist *templateRDH, *templateRDH2; RooHistPdf *templateRHP, *templateRHP2;
 };
 // Parent class definitions
 ////////////////////////////////////////////////////
 // Child class definitions
 ////////////////////////////////////////////////////
 templateConvDas::templateConvDas(RooRealVar &m, TH1D *templateHist_) {
   plotLabel="Sig.: MC #otimes Das";
   templateHist=templateHist_; assert(templateHist);
   mean   = new RooRealVar("sig_mean" , "sig_mean" ,    0,    -5,  2);
   sigma  = new RooRealVar("sig_sigma", "sig_sigma", 0.01, 0.001,  2);
   kLo    = new RooRealVar("sig_kLo"  , "sig_kLo"  ,  1.5,    .3, 10);
   kHi    = new RooRealVar("sig_kHi"  , "sig_kHi"  ,  1.5,    .3, 10);
   dd = new RooGaussDoubleSidedExp("dd","dd",m,*mean,*sigma,*kLo,*kHi);
   templateRDH = new RooDataHist("templateRDH","templateRDH",RooArgSet(m),templateHist_);
   templateRHP = new RooHistPdf("templateRHP","templateRHP",m,*templateRDH,2);
   theShape=new RooFFTConvPdf("signalModel","signalModel",m,*templateRHP,*dd);
   resolutionFunction=dd;
 } 
 templateConvDas::~templateConvDas() {
  delete mean; delete sigma; delete kLo; delete kHi;
  delete dd; delete templateRDH; delete templateRHP;
 }
 ////////////////////////////////////////////////////
 templateConvGaus::templateConvGaus(RooRealVar &m, TH1D *templateHist_) {
   plotLabel="Sig.: MC #otimes Gaus";
   templateHist=templateHist_; assert(templateHist);
   mean   = new RooRealVar("sig_mean" , "sig_mean" ,    0,    -5,  2);
   sigma  = new RooRealVar("sig_sigma", "sig_sigma", 0.1, 0.001,  2);
   gaus   = new RooGaussian("gaus","gaus",m,*mean,*sigma);
   templateRDH = new RooDataHist("templateRDH","templateRDH",RooArgSet(m),templateHist_);
   templateRHP = new RooHistPdf("templateRHP","templateRHP",m,*templateRDH,2);
   theShape=new RooFFTConvPdf("signalModel","signalModel",m,*templateRHP,*gaus);
   resolutionFunction=gaus;
 } 
 templateConvGaus::~templateConvGaus() {
  delete templateRHP; delete templateRDH;
  delete gaus; delete mean; delete sigma;
 }

 ////////////////////////////////////////////////////
 templateConvDoubleGaus::templateConvDoubleGaus(RooRealVar &m, TH1D *templateHist_) {
   /*plotLabel="Sig.: MC #otimes Double Gaus";
   templateHist=templateHist_; assert(templateHist);
   mean1   = new RooRealVar("sig_mean1" , "sig_mean1" ,    0,    -5,  2);
   sigma1  = new RooRealVar("sig_sigma1", "sig_sigma1", 0.01, 0.001,  5);
   gaus1   = new RooGaussian("gaus1","gaus1",m,*mean1,*sigma1);
   mean2   = new RooRealVar("sig_mean2" , "sig_mean2" ,    0,    -5,  2);
   sigma2  = new RooRealVar("sig_sigma2", "sig_sigma2", 0.01, 0.001,  5);
   gaus2   = new RooGaussian("gaus2","gaus2",m,*mean2,*sigma2);
   frac = new RooRealVar("sig_frac","sig_frac", 0.05, 0.,1.);
   templateRDH = new RooDataHist("templateRDH","templateRDH",RooArgSet(m),templateHist_);
   templateRHP = new RooHistPdf("templateRHP","templateRHP",m,*templateRDH,2);
   sum = new RooAddPdf("sum","sum",RooArgList(*gaus1,*gaus2),RooArgList(*frac));
   theShape = new RooFFTConvPdf("signalModel","signalModel",m,*templateRHP,*sum);
   resolutionFunction=sum;*/
 
   plotLabel="Sig.: MC #otimes Double Gaus";
   templateHist=templateHist_; assert(templateHist);
   mean1   = new RooRealVar("sig_mean1" , "sig_mean1" ,    0,    -5,  2);
   sigma1  = new RooRealVar("sig_sigma1", "sig_sigma1", 0.01, 0.001,  5);
   gaus1   = new RooGaussian("gaus1","gaus1",m,*mean1,*sigma1);
   mean2   = new RooRealVar("sig_mean2" , "sig_mean2" ,    0,    -5,  2);
   sigma2  = new RooRealVar("sig_sigma2", "sig_sigma2", 0.01, 0.001,  5);
   gaus2   = new RooGaussian("gaus2","gaus2",m,*mean2,*sigma2);
   frac = new RooRealVar("sig_frac","sig_frac", 0.05, 0.,1.);
   templateRDH = new RooDataHist("templateRDH","templateRDH",RooArgSet(m),templateHist_);
   templateRHP = new RooHistPdf("templateRHP","templateRHP",m,*templateRDH,2);
   theShape1 = new RooFFTConvPdf("signalModel1","signalModel1",m,*templateRHP,*gaus1);
   theShape2 = new RooFFTConvPdf("signalModel2","signalModel2",m,*templateRHP,*gaus2);
   theShape = new RooAddPdf("signalModel","signalModel",RooArgList(*theShape1,*theShape2),RooArgList(*frac));              
   sum = new RooAddPdf("sum","sum",RooArgList(*gaus1,*gaus2),RooArgList(*frac));                                                             
   resolutionFunction=sum;


}
 templateConvDoubleGaus::~templateConvDoubleGaus() {
   delete templateRHP; delete templateRDH;
   delete gaus1; delete mean1; delete sigma1;
   delete gaus2; delete mean2; delete sigma2; 
   delete frac; delete sum;
   delete theShape1; delete theShape2; delete theShape;
 }
 ////////////////////////////////////////////
 templateConvLandau::templateConvLandau(RooRealVar &m, TH1D *templateHist_) {
   plotLabel="Sig.: MC #otimes Landau";
   templateHist=templateHist_; assert(templateHist);
   mean   = new RooRealVar("sig_mean" , "sig_mean" ,    0,    -5,  2);
   sigma  = new RooRealVar("sig_sigma", "sig_sigma", 0.01, 0.001,  5);
   x = new RooFormulaVar("arg","arg","-@0",RooArgList(m));
   landau = new RooLandau ("landau", "landau", *x, *mean, *sigma);
   templateRDH = new RooDataHist("templateRDH","templateRDH",RooArgSet(m),templateHist_);
   templateRHP = new RooHistPdf("templateRHP","templateRHP",m,*templateRDH,2);
   theShape=new RooFFTConvPdf("signalModel","signalModel",m,*templateRHP,*landau);
   resolutionFunction=landau;

 }
 templateConvLandau::~templateConvLandau() {
   delete templateRHP; delete templateRDH;
   delete landau; delete mean; delete sigma;
   delete x;
 }

 ////////////////////////////////////////////                                                                                                 
 templateConvBeta::templateConvBeta(RooRealVar &m, TH1D *templateHist_) {
   plotLabel="Sig.: MC #otimes Beta";
   templateHist=templateHist_; assert(templateHist);
   a = new RooRealVar("alpha" , "alpha" ,    0,    0.001,  5);
   b = new RooRealVar("beta", "beta", 0.01, 0.001,  5);
   beta = bindPdf("beta",ROOT::Math::beta_pdf,m,*a,*b) ;   
   templateRDH = new RooDataHist("templateRDH","templateRDH",RooArgSet(m),templateHist_);
   templateRHP = new RooHistPdf("templateRHP","templateRHP",m,*templateRDH,2);
   theShape=new RooFFTConvPdf("signalModel","signalModel",m,*templateRHP,*beta);
   resolutionFunction=beta;
 }
 templateConvBeta::~templateConvBeta() {
   delete templateRHP; delete templateRDH;
   delete beta; delete b; delete a;
   }
 

 ////////////////////////////////////////////////////                                                                                         
 templateConvCrystalBall::templateConvCrystalBall(RooRealVar &m, TH1D *templateHist_) {
   plotLabel="Sig.: MC #otimes CB";
   templateHist=templateHist_; assert(templateHist);
   mean   = new RooRealVar("sig_mean" , "sig_mean" ,    0,    -5,  2);
   sigma  = new RooRealVar("sig_sigma", "sig_sigma", 0.01, 0.001,  5);
   a   = new RooRealVar("sig_a" , "sig_a" ,    0,  -5  ,  2);
   n  = new RooRealVar("sig_n", "sig_n", 0.01, 0.001,  5);
   CB =  new RooCBShape("crystalball","crystalball",m,*mean,*sigma,*a,*n);
   templateRDH = new RooDataHist("templateRDH","templateRDH",RooArgSet(m),templateHist_);
   templateRHP = new RooHistPdf("templateRHP","templateRHP",m,*templateRDH,2);
   theShape=new RooFFTConvPdf("signalModel","signalModel",m,*templateRHP,*CB);
   resolutionFunction=CB;
 }
 templateConvCrystalBall::~templateConvCrystalBall() {
   delete templateRHP; delete templateRDH;
   delete CB; delete mean; delete sigma; delete a; delete n;
}
 ////////////////////////////////////////////////////
 templateBWConvGaus::templateBWConvGaus(RooRealVar &m, TH1D *templateHist_) {
   plotLabel="Sig.: (MC+BW)#otimesGaus";
   templateHist=templateHist_; assert(templateHist);
   mean   = new RooRealVar("sig_mean" , "sig_mean" ,    0,    -5,  2);
   sigma  = new RooRealVar("sig_sigma", "sig_sigma", 0.01, 0.001,  5);
   gaus   = new RooGaussian("gaus","gaus",m,*mean,*sigma);
   templateRDH = new RooDataHist("templateRDH","templateRDH",RooArgSet(m),templateHist_);
   templateRHP = new RooHistPdf("templateRHP","templateRHP",m,*templateRDH,2);
   resolutionFunction=gaus;
   
   mass   = new RooRealVar("sig_mass" , "sig_mass" ,91.1876,    80, 99); mass->setConstant(kTRUE);
   width  = new RooRealVar("sig_width", "sig_width", 2.4952,    .1, 10); width->setConstant(kTRUE);
   bw = new RooBreitWigner("bw","bw",m,*mass,*width);
   frac = new RooRealVar("sig_frac","sig_frac", 0.05, 0.,1.);
   templatePlusBw = new RooAddPdf("templatePlusBw","templatePlusBw", RooArgList(*templateRHP,*bw),RooArgList(*frac));
   theShape=new RooFFTConvPdf("signalModel","signalModel",m,*templatePlusBw,*gaus);
 } 
 templateBWConvGaus::~templateBWConvGaus() {
  delete templateRHP; delete templateRDH;
  delete gaus; delete mean; delete sigma; delete mass; delete width; delete frac;
  delete bw; delete templatePlusBw;
 }
 ////////////////////////////////////////////////////
 breitWignerConvDas::breitWignerConvDas(RooRealVar &m, TH1D *templateHist_) {
   plotLabel="Sig.: BW #otimes Das";
   mean   = new RooRealVar("sig_mean" , "sig_mean"  ,     0,    -5,  2);
   sigma  = new RooRealVar("sig_sigma", "sig_sigma",      1,   0.1,  6);
   kLo    = new RooRealVar("sig_kLo"  , "sig_kLo"  ,    1.5,    .3, 10);
   kHi    = new RooRealVar("sig_kHi"  , "sig_kHi"  ,      5,    .3, 10);
   mass   = new RooRealVar("sig_mass" , "sig_mass" ,91.1876,    80, 99); mass->setConstant(kTRUE);
   width  = new RooRealVar("sig_width", "sig_width", 2.4952,    .1, 10); width->setConstant(kTRUE);
   bw = new RooBreitWigner("bw","bw",m,*mass,*width);
   dd = new RooGaussDoubleSidedExp("dd","dd",m,*mean,*sigma,*kLo,*kHi);
   resolutionFunction=dd;
   theShape=new RooFFTConvPdf("signalModel","signalModel",m,*bw,*dd);
 } 
 breitWignerConvDas::~breitWignerConvDas() {
  delete dd; delete bw;
  delete mean; delete sigma; delete kLo; delete kHi; delete mass; delete width;
 }
 ////////////////////////////////////////////////////
 bkgErfcExp::bkgErfcExp(RooRealVar &m, TH1D *templateHist_) {
  plotLabel="Bkg.: Erfc #times exp"; 
  mpv   = new RooRealVar("bkg_mpv"  ,"bkgEE_mpv"  ,      60,    40,  150);
  beta  = new RooRealVar("bkg_beta" ,"bkgEE_beta" ,    0.05, 0.03,   .2); 
  gamma = new RooRealVar("bkg_gamma","bkgEE_gamma",    0.03, 0.01,    1); 
  zmass = new RooRealVar("bkg_mZ"   ,"bkgEE_mZ"   , 91.1876,   85,   97); zmass->setConstant(kTRUE);
  theShape = new RooErfcExpoPdf("bkgModel","bkgModel",m,*mpv,*beta,*gamma,*zmass);
 }
 bkgErfcExp::~bkgErfcExp() {
  delete mpv; delete beta; delete gamma; delete zmass;
 }
 ////////////////////////////////////////////////////
 bkgErfcExpPlusExp::bkgErfcExpPlusExp(RooRealVar &m, TH1D *templateHist_) {
  plotLabel="Bkg.: Erfc #times exp + exp"; 
  mpv   = new RooRealVar("bkg_mpv"  ,"bkgEE_mpv"  ,      70,    0,  250);
  beta  = new RooRealVar("bkg_beta" ,"bkgEE_beta" ,    0.05, 0.03,   .2); 
  gamma = new RooRealVar("bkg_gamma","bkgEE_gamma",    0.03, 0.01,    1); 
  zmass = new RooRealVar("bkg_mZ"   ,"bkgEE_mZ"   , 91.1876,   85,   97); zmass->setConstant(kTRUE);
  erfcExp = new RooErfcExpoPdf("erfcExp","erfcExp",m,*mpv,*beta,*gamma,*zmass);
  t1   = new RooRealVar("bkg_t1"  ,"bkg_t1"  ,-0.20,-1.,0.);
  frac = new RooRealVar("bkg_frac","bkg_frac", 0.9, 0.,1.);
  exp1 = new RooExponential("bkg_exp1","bkg_exp1",m,*t1);
  theShape = new RooAddPdf("bkgModel","bkgModel",RooArgList(*erfcExp,*exp1),RooArgList(*frac));
 }
 bkgErfcExpPlusExp::~bkgErfcExpPlusExp() {
  delete mpv; delete beta; delete gamma; delete zmass; delete t1; delete frac;
  delete exp1; delete erfcExp;
 }
 ////////////////////////////////////////////////////
 bkgDasPlusExp::bkgDasPlusExp(RooRealVar &m, TH1D *templateHist_) {
  plotLabel="Bkg.: Das + exp"; 
  //mean   = new RooRealVar("bkg_mean" , "bkg_mean" ,   60,    30,200);
  //sigma  = new RooRealVar("bkg_sigma", "bkg_sigma",   12,    10, 60);
  //kLo    = new RooRealVar("bkg_kLo"  , "bkg_kLo"  ,  1.5,   .02, 10);
  //kHi    = new RooRealVar("bkg_kHi"  , "bkg_kHi"  ,  1.5,   .02, 10);
  mean   = new RooRealVar("bkg_mean" , "bkg_mean" ,   90,    30,200);
  sigma  = new RooRealVar("bkg_sigma", "bkg_sigma",   12,    10, 60);
  kLo    = new RooRealVar("bkg_kLo"  , "bkg_kLo"  ,  1.5,   .02, 10);
  kHi    = new RooRealVar("bkg_kHi"  , "bkg_kHi"  ,  1.5,   .02, 10);
  dd = new RooGaussDoubleSidedExp("bkgDas","bkgDas",m,*mean,*sigma,*kLo,*kHi);
  t1   = new RooRealVar("bkg_t1"  ,"bkg_t1"  ,-0.20,-.4,.4);
  frac = new RooRealVar("bkg_frac","bkg_frac", 0.05, 0.,1.);
  exp1 = new RooExponential("bkg_exp1","bkg_exp1",m,*t1);
  theShape = new RooAddPdf("bkgModel","bkgModel",RooArgList(*dd,*exp1),RooArgList(*frac));
 }
 bkgDasPlusExp::~bkgDasPlusExp() {
  delete mean; delete sigma; delete kLo; delete kHi; delete t1; delete frac;
  delete exp1; delete dd;
 }
 ////////////////////////////////////////////////////
 bkgDas::bkgDas(RooRealVar &m, TH1D *templateHist_) {
  plotLabel="Bkg.: Wide Das"; 
  mean   = new RooRealVar("bkg_mean" , "bkg_mean" ,   60,    30,200);
  sigma  = new RooRealVar("bkg_sigma", "bkg_sigma",   12,    10,60);//30
  kLo    = new RooRealVar("bkg_kLo"  , "bkg_kLo"  ,  1.5,   .02, 10); 
  kHi    = new RooRealVar("bkg_kHi"  , "bkg_kHi"  ,  1.5,   .02, 10);
  theShape = new RooGaussDoubleSidedExp("bkgModel","bkgModel",m,*mean,*sigma,*kLo,*kHi);
 }
 bkgDas::~bkgDas() {
  delete mean; delete sigma; delete kLo; delete kHi;
 }
 ////////////////////////////////////////////////////
 bkgDoubleExp::bkgDoubleExp(RooRealVar &m, TH1D *templateHist_) {
  plotLabel="Bkg.: Exp + exp"; 
  t1   = new RooRealVar("bkg_t1"  ,"bkg_t1"  ,-0.20,-1.,0.);
  t2   = new RooRealVar("bkg_t2"  ,"bkg_t2"  ,-0.05,-1.,0.);
  frac = new RooRealVar("bkg_frac","bkg_frac", 0.50, 0.,1.);
  exp1 = new RooExponential("bkg_exp1","bkg_exp1",m,*t1);
  exp2 = new RooExponential("bkg_exp2","bkg_exp2",m,*t2);
  theShape = new RooAddPdf("bkgModel","bkgModel",RooArgList(*exp1,*exp2),RooArgList(*frac));
 }
 bkgDoubleExp::~bkgDoubleExp() {
  delete exp1; delete exp2;
  delete t1; delete t2; delete frac;
 }
 ////////////////////////////////////////////////////
 bkgTemplate::bkgTemplate(RooRealVar &m, TH1D *templateHist_, TString plotLabel_) {
   if(plotLabel_!="")plotLabel=plotLabel_; else plotLabel="Bkg.: Data driven";
   templateHist=templateHist_; assert(templateHist);
   templateRDH = new RooDataHist("bkgTemplateRDH","bkgTemplateRDH",RooArgSet(m),templateHist_);
   theShape=new RooHistPdf("bkgModel","bkgModel",m,*templateRDH,2);
 } 
 bkgTemplate::~bkgTemplate() {
  delete templateRDH;
 }
 ////////////////////////////////////////////////////
 bkgTemplatePlusExp::bkgTemplatePlusExp(RooRealVar &m, TH1D *templateHist_, TString plotLabel_) {
   if(plotLabel_!="")plotLabel=plotLabel_; else plotLabel="Bkg.: Data driven";
   templateHist=templateHist_; assert(templateHist);
   templateRDH = new RooDataHist("bkgTemplatePlusExpRDH","bkgTemplatePlusExpRDH",RooArgSet(m),templateHist_);
   templateRHP = new RooHistPdf("bkgTemplatePlusExpRHP","bkgTemplatePlusExpRHP",m,*templateRDH,2);
   t1   = new RooRealVar("bkg_t1"  ,"bkg_t1"  ,-0.20,-1.,0.);
   frac = new RooRealVar("bkg_frac","bkg_frac", 0.90, 0.,1.);
   exp1 = new RooExponential("bkg_exp1","bkg_exp1",m,*t1);
   theShape=new RooAddPdf("bkgModel","bkgModel",RooArgList(*templateRHP,*exp1),RooArgList(*frac));
 } 
 bkgTemplatePlusExp::~bkgTemplatePlusExp() {
  delete templateRHP; delete templateRDH;
  delete t1; delete frac; delete exp1;
 }
 ////////////////////////////////////////////////////
 bkgTwoTemplates::bkgTwoTemplates(RooRealVar &m, TH1D *templateHist_, TH1D *templateHist2_, TString plotLabel_) {
   if(plotLabel_!="")plotLabel=plotLabel_; else plotLabel="Bkg.: Data driven";
   templateHist =templateHist_ ; assert(templateHist );
   templateHist2=templateHist2_; assert(templateHist2);
   templateRDH = new RooDataHist("bkgTemplateRDH","bkgTemplateRDH",RooArgSet(m),templateHist_);
   templateRHP = new RooHistPdf("bkgTemplateRHP","bkgTemplateRHP",m,*templateRDH,2);
   templateRDH2 = new RooDataHist("bkgTemplateRDH2","bkgTemplateRDH2",RooArgSet(m),templateHist2_);
   templateRHP2 = new RooHistPdf("bkgTemplateRHP2","bkgTemplateRHP2",m,*templateRDH2,2);
   frac = new RooRealVar("bkg_frac","bkg_frac", 0.5, 0.,1.);
   theShape=new RooAddPdf("bkgModel","bkgModel",RooArgList(*templateRHP,*templateRHP2),RooArgList(*frac));
 } 
 bkgTwoTemplates::~bkgTwoTemplates() {
  delete templateRHP; delete templateRDH;
  delete templateRHP2; delete templateRDH2;
  delete frac;
 }
 ////////////////////////////////////////////////////
} 
/*
}void fitBin(
  string histName="pass_ptBin0_etaBin0",
  string dataFileName="",
  string templateFileName="",
  string signalModel="template",
  string bkgModel="erfcexp",
  TString plotTitle=""
  //metadata file name with kinematics info, initial params, info for plot 
) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  TFile *dataFile=TFile::Open(dataFileName.c_str(),"READ"); assert(dataFile && dataFile->IsOpen());
  TH1D *dataHist=(TH1D*)dataFile->Get(histName.c_str()); assert(dataHist); dataHist->SetDirectory(0);
  TFile *templateFile; TH1D *templateHist_;
  if(templateFileName!="") {
    templateFile=TFile::Open(templateFileName.c_str(),"READ"); assert(templateFile);
    templateHist_=(TH1D*)templateFile->Get(histName.c_str());
  }  
  assert(signalModel!="template" || templateFile);

  //todo: declare models somewhere else
  RooRealVar m("m","mass",60.,120.);
  m.setBins(10000);

  // Signal model: template conv. CB for now
  //RooRealVar *mean  = new RooRealVar("sigCB_mean" ,"sigCB_mean" ,-2,-10,10);
  //RooRealVar *sigma = new RooRealVar("sigCB_sigma","sigCB_sigma",0.1,0.01,5);
  //RooRealVar *alpha = new RooRealVar("sigCB_alpha","sigCB_alpha",5,0,20);
  //RooRealVar *n     = new RooRealVar("sigCB_n"    ,"sigCB_n"    ,1,0,10);
  //RooCBShape *cb = new RooCBShape("cb","cb",m,*mean,*sigma,*alpha,*n);
  //RooDataHist *templateRDH = new RooDataHist("templateRDH","templateRDH",RooArgSet(m),templateHist_);
  //RooHistPdf  *templateRHP = new RooHistPdf("templateRHP","templateRHP",m,*templateRDH,2);
  //RooAbsPdf *signalModelPdf=new RooFFTConvPdf("signalModel","signalModel",m,*templateRHP,*cb);
  
  // Signal model : das distribution
  RooRealVar *mean   = new RooRealVar("sigDD_mean","sigDD_mean",0,-5,2);
  RooRealVar *sigma  = new RooRealVar("sigDD_sigma","sigDD_sigma",0.01,0.001,3);
  RooRealVar *kLo    = new RooRealVar("sigDD_kLo","sigDD_kLo",1.5,.3,10);
  RooRealVar *kHi    = new RooRealVar("sigDD_kHi","sigDD_kHi",1.5,.3,10);
  RooGaussDoubleSidedExp *dd = new RooGaussDoubleSidedExp("dd","dd",m,*mean,*sigma,*kLo,*kHi);
  
  RooDataHist *templateRDH = new RooDataHist("templateRDH","templateRDH",RooArgSet(m),templateHist_);
  RooHistPdf  *templateRHP = new RooHistPdf("templateRHP","templateRHP",m,*templateRDH,2);
  RooAbsPdf *signalModelPdf=new RooFFTConvPdf("signalModel","signalModel",m,*templateRHP,*dd);
  //RooRealVar *vMean  = new RooRealVar("sigV_mean","sigV_mean",91.1876, 86.1876, 96.1876); vMean->setConstant(kTRUE);
  //RooRealVar *vWidth = new RooRealVar("sigV_width","sigV_width",2.495,2.49,2.5); vWidth->setConstant(kTRUE);
  //RooRealVar *vSigma = new RooRealVar("sigV_sigma","sigV_sigma",1, 0.01, 2);
  //RooRealVar *vFrac  = new RooRealVar("sigV_frac","sigV_frac",0.005, 0, .5);
  //char formula[100];
  //sprintf(formula, "1 - sigV_frac");
  //RooFormulaVar *oneMinusVFrac = new RooFormulaVar("oneMinusVFrac","oneMinusVFrac",formula, *vFrac);
  //RooVoigtian *voigt = new RooVoigtian("voigt","voigt", m, *vMean, *vWidth, *vSigma);
  //RooAbsPdf *templateConvDD=new RooFFTConvPdf("templateConvDD","templateConvDD",m,*templateRHP,*dd);
  //RooAbsPdf *signalModelPdf = new RooAddPdf("signalModel","signalModel", RooArgList(*voigt, *templateConvDD), RooArgList(*vFrac, *oneMinusVFrac));
  
  // Background model: erfc * expo for now
  RooRealVar *mpv   = new RooRealVar("bkgEE_mpv"  ,"bkgEE_mpv"  ,70,0,250    );
  RooRealVar *beta  = new RooRealVar("bkgEE_beta" ,"bkgEE_beta" ,0.05,0.03,.2); 
  RooRealVar *gamma = new RooRealVar("bkgEE_gamma","bkgEE_gamma",0.03,0,1    ); 
  RooRealVar *zmass = new RooRealVar("bkgEE_mZ"   ,"bkgEE_mZ"   ,91.1876,85,97); zmass->setConstant(kTRUE);
  RooAbsPdf *bkgModelPdf = new RooErfcExpoPdf("bkgModel","bkgModel",m,*mpv,*beta,*gamma,*zmass);
  //if(dataHist->Integral(10,20)/dataHist->Integral(1,30) > 0.8) {
  //  beta->setConstant(kTRUE);
  //  gamma->setConstant(kTRUE);
  //}
  
  // double exponential
  //RooRealVar *k1 = new RooRealVar("bkg2E_k1","bkg2E_k1",-0.03,-.2,0.);
  //RooRealVar *k2 = new RooRealVar("bkg2E_k2","bkg2E_k2",-0.01,-.2,0.);
  //RooRealVar *frac = new RooRealVar("bkg2E_frac", "bkg2E_frac", 0.50, 0.,1.);
  //RooExponential *exp1 = new RooExponential("exp1","exp1",m,*k1);
  //RooExponential *exp2 = new RooExponential("exp2","exp2",m,*k2);
  //RooAbsPdf *bkgModelPdf = new RooAddPdf("bkgModel","bkgModel",RooArgList(*exp1,*exp2),RooArgList(*frac));
  //RooRealVar *bkg_mpv    = new RooRealVar("bkg_mpv","bkg_mpv",40,0,250);
  //RooRealVar *bkg_sigma  = new RooRealVar("bkg_sigma","bkg_sigma",8,6,50);
  //RooRealVar *bkg_kLo    = new RooRealVar("bkg_kLo","bkg_kHi",1.50,.5,10);
  //RooRealVar *bkg_kHi    = new RooRealVar("bkg_kHi","bkg_kHi",1.5,.5,10);
  //RooAbsPdf *bkgModelPdf = new RooGaussDoubleSidedExp("bkgModel","bkgModel",m,*bkg_mpv,*bkg_sigma,*bkg_kLo,*bkg_kHi);
 
  // Signal+background model
  RooRealVar Ntotal = RooRealVar("Ntotal","Ntotal",dataHist->Integral(),0,dataHist->Integral()); Ntotal.setConstant(kTRUE);
  RooRealVar Nbkg = RooRealVar("Nbkg","Nbkg",0.2*dataHist->Integral(), 0, dataHist->Integral());
  RooFormulaVar Nsig = RooFormulaVar("Nsig","Nsig","Ntotal - Nbkg", RooArgList(Ntotal,Nbkg));
  RooAbsPdf *totalPdf = new RooAddPdf("totalPdf","totalPdf", RooArgList(*signalModelPdf,*bkgModelPdf), RooArgList(Nsig,Nbkg));
  // Data histogram in RooFit
  RooDataHist *dataRDH = new RooDataHist("dataRDH","dataRDH",RooArgSet(m), dataHist);
  //Perform the fit
  TF1 *bkgEstimator = new TF1("bkgEstimator","expo",60,120);
  dataHist->Fit(bkgEstimator, "N0Q", "goff", 60,80);
  Nbkg.setVal(TMath::Min(dataHist->Integral(),bkgEstimator->Integral(60,120)));
  delete bkgEstimator;
  sigma->setVal(dataHist->GetBinCenter(dataHist->GetMaximumBin())-templateHist_->GetBinCenter(dataHist->GetMaximumBin()));
  RooFitResult *fitResult=0;
  fitResult = totalPdf->fitTo(*dataRDH,RooFit::Extended(),RooFit::Strategy(2),RooFit::NumCPU(1),RooFit::Save());
  
  // Combined errors and goodness-of-fit test
  TH1D *dataHistWithCombErrors=(TH1D*)dataHist->Clone("dataHistWithCombErrors"); dataHistWithCombErrors->SetDirectory(0);
  TH1D *ratioDataPdf=(TH1D*)dataHist->Clone("ratioDataPdf"); ratioDataPdf->SetDirectory(0); ratioDataPdf->Scale(0); ratioDataPdf->Clear();
  TH1D *pdfErrorBand=(TH1D*)dataHist->Clone("pdfErrorBand"); pdfErrorBand->SetDirectory(0); pdfErrorBand->Scale(0); pdfErrorBand->Clear();
  TH1D *ratioErrorBand=(TH1D*)dataHist->Clone("ratioErrorBand"); ratioErrorBand->SetDirectory(0); ratioErrorBand->Scale(0); ratioErrorBand->Clear();

  double chiSquaredSum=0;
  printf("###########################################\n");
  printf("# computing chi square goodness-of-fit :D #\n");
  printf("###########################################\n");
  // Should we propagate the smearing PDF to the MC errors? not a big deal for now
  RooArgSet *floatingPars = (RooArgSet*)totalPdf->getParameters(m)->selectByAttrib("Constant",kFALSE);
  for(int nb=1; nb<=dataHistWithCombErrors->GetNbinsX(); nb++) {
    // Integrate the pdf
    double xmin=dataHistWithCombErrors->GetBinLowEdge(nb);
    double xmax=dataHistWithCombErrors->GetBinLowEdge(nb+1);
    if(xmin>=120.||xmax<=60.) continue;
    m.setRange("theIntegral",xmin,xmax);
    double modelVal = totalPdf->createIntegral(m,m,"theIntegral")->getVal() * Ntotal.getVal();

    double theError,theError2;
    if(templateFileName!="") {
      int mcBin=templateHist_->FindBin(dataHistWithCombErrors->GetBinCenter(nb));
      theError2=(
      pow(dataHistWithCombErrors->GetBinError(nb)/dataHistWithCombErrors->GetBinContent(nb),2)+
      pow(templateHist_->GetBinError(mcBin)/templateHist_->GetBinContent(mcBin),2)
      )*pow(dataHistWithCombErrors->GetBinContent(nb),2);
      theError=sqrt(theError2);
      ratioErrorBand->SetBinError(nb, modelVal*templateHist_->GetBinError(mcBin)/templateHist_->GetBinContent(mcBin)/dataHist->GetBinContent(nb));
      //m.setVal((xmin+xmax)/2.);
      //pdfErrorBand->SetBinContent(nb, 2.*totalPdf->getVal(RooArgSet(m))*Ntotal.getVal());
      //pdfErrorBand->SetBinError(nb,templateHist_->GetBinError(mcBin)/templateHist_->GetBinContent(mcBin)*2.*totalPdf->getVal(RooArgSet(m))*Ntotal.getVal()); //2 = bin width
      pdfErrorBand->SetBinContent(nb, modelVal);
      pdfErrorBand->SetBinError(nb,templateHist_->GetBinError(mcBin)/templateHist_->GetBinContent(mcBin)*modelVal);
    } else { theError = dataHistWithCombErrors->GetBinError(nb); theError2 = theError*theError; }
    //dataHistWithCombErrors->SetBinError(nb,theError);
    double chiSquaredContribution=pow(dataHistWithCombErrors->GetBinContent(nb)-modelVal ,2)/theError2;
    chiSquaredSum+=chiSquaredContribution;
    printf("\tbin %3d: pdf=%9.2f, data=%9d, error=%9.2f, contributing %8.4f to chi2\n",nb,modelVal,(int)dataHistWithCombErrors->GetBinContent(nb),theError,chiSquaredContribution);

    ratioDataPdf->SetBinContent(nb,dataHistWithCombErrors->GetBinContent(nb) / modelVal);
    ratioDataPdf->SetBinError(nb,dataHist->GetBinError(nb) * modelVal / pow(dataHistWithCombErrors->GetBinContent(nb),2));
    ratioErrorBand->SetBinContent(nb, 1);
  }
  int ndf=dataHistWithCombErrors->FindBin(120)-dataHistWithCombErrors->FindBin(60)+1 - floatingPars->getSize();
  double fitProb=TMath::Prob(chiSquaredSum,ndf);
  printf("fit probability: %4.2f%%\n",100.*fitProb);
  printf("signal events  : %9.2f\n",Nsig.getVal());
  printf("bkg events     : %9.2f\n",Nbkg.getVal());

*/
