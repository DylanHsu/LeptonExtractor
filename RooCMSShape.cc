/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: PhysicsTools/TagAndProbe/RooCMSShape
 *
 *
 * Authors:
 *   Nadia Adam, Princeton - neadam@princeton.edu
 *   Adam Hunt, Princeton  - ahunt@princeton.edu
 *   Kalanand Mishra, Fermilab - kalanand@fnal.gov
 *
 * Description:
 *   Defines a probability density function which has exponential decay 
 *   distribution at high mass beyond the pole position (say, Z peak)  
 *   but turns over (i.e., error function) at low mass due to threshold 
 *   effect. We use this to model the background shape in Z->ll invariant 
 *   mass.
 * History:
 *   
 *
 * Copyright (C) 2008 FNAL 
 *****************************************************************************/

#include "RooCMSShape.h"

//ClassImp(RooCMSShape) 

RooCMSShape::RooCMSShape(const char *name, const char *title, 
		       RooAbsReal& _x,
		       RooAbsReal& _alpha,
		       RooAbsReal& _beta,
		       RooAbsReal& _gamma,
		       RooAbsReal& _peak) :
  RooAbsPdf(name,title), 
  x("x","x",this,_x),
  alpha("alpha","alpha",this,_alpha),
  beta("beta","beta",this,_beta),
  gamma("gamma","gamma",this,_gamma),
  peak("peak","peak",this,_peak)
{ } 


RooCMSShape::RooCMSShape(const RooCMSShape& other, const char* name):
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  alpha("alpha",this,other.alpha),
  beta("beta",this,other.beta),
  gamma("gamma",this,other.gamma),
  peak("peak",this,other.peak)
{ } 



Double_t RooCMSShape::evaluate() const 
{ 
 // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 

 //Double_t erf = TMath::Erfc((alpha - x) * beta);
 Double_t erf = RooMath::erfc((alpha - x) * beta);
 Double_t u = (x - peak)*gamma;

 if(u < -70) u = 1e20;
 else if( u>70 ) u = 0;
 else u = exp(-u);   //exponential decay
 return erf*u;
}

RooErfcExpoPdf::RooErfcExpoPdf(const char *name, const char *title, 
		       RooAbsReal& _x,
		       RooAbsReal& _mpv,
		       RooAbsReal& _beta,
		       RooAbsReal& _gamma,
		       RooAbsReal& _zmass) :
  RooAbsPdf(name,title), 
  x("x","x",this,_x),
  mpv("mpv","mpv",this,_mpv),
  beta("beta","beta",this,_beta),
  gamma("gamma","gamma",this,_gamma),
  zmass("zmass","zmass",this,_zmass)
{ } 


RooErfcExpoPdf::RooErfcExpoPdf(const RooErfcExpoPdf& other, const char* name):
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  mpv("mpv",this,other.mpv),
  beta("beta",this,other.beta),
  gamma("gamma",this,other.gamma),
  zmass("zmass",this,other.zmass)
{ } 



Double_t RooErfcExpoPdf::evaluate() const 
{ 
 // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 

 //Double_t erf = TMath::Erfc((mpv - x) * beta);
 //Double_t erf = RooMath::erfc((mpv - x) * beta + gamma/2./beta);
 Double_t arg1 = (mpv - x) * beta + gamma/2./beta;
 Double_t arg2 = gamma/2./beta;
 Double_t erf = TMath::Erfc(arg1)/TMath::Erfc(arg2);
 
 Double_t u = TMath::Exp(-gamma*(x - zmass));
 //Double_t u = (x - zmass)*gamma;
 //if(u < -70) u = 1e20;
 //else if( u>70 ) u = 0;
 //else u = exp(-u);   //exponential decay
 return erf*u;
} 
