/**
 * @Author: Billy Li <billyli>
 * @Date:   04-25-2022
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 04-26-2022
 */



/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

// Your description goes here...

#include "Riostream.h"

#include "RooCBShape_MN.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"

ClassImp(RooCBShape_MN);

Double_t RooCBShape_MN::ApproxErf(Double_t arg) const
{
  static const double erflim = 5.0;
  if( arg > erflim )
    return 1.0;
  if( arg < -erflim )
    return -1.0;

  return RooMath::erf(arg);
}

RooCBShape_MN::RooCBShape_MN(const char *name, const char *title,
                      RooAbsReal& _m,
                      RooAbsReal& _m0,
                      RooAbsReal& _sigma,
                      RooAbsReal& _alpha,
                      RooAbsReal& _n) :
 RooAbsPdf(name,title),
 m("m","m",this,_m),
 m0("m0","m0",this,_m0),
 sigma("sigma","sigma",this,_sigma),
 alpha("alpha","alpha",this,_alpha),
 n("n","n",this,_n)
{
}


RooCBShape_MN::RooCBShape_MN(const RooCBShape_MN& other, const char* name) :
  RooAbsPdf(other,name),
  m("m",this,other.m),
  m0("m0",this,other.m0),
  sigma("sigma",this,other.sigma),
  alpha("alpha",this,other.alpha),
  n("n",this,other.n)
{
}


Double_t RooCBShape_MN::evaluate() const
{
  Double_t t = (m-m0)/sigma;
  if (alpha < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)alpha);
  if (t >= -absAlpha) {
    return exp(-0.5*t*t);
  }
  else {
    // Double_t a =  TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    // Double_t b= n/absAlpha - absAlpha;

    // return a/TMath::Power(b - t, n);

    Double_t exp_part = exp(-0.5*absAlpha*absAlpha);
    Double_t D = (n-absAlpha*absAlpha-absAlpha*t);
    Double_t arg = n/D;

    return TMath::Power(arg, n)*exp_part;
  }
}

Int_t RooCBShape_MN::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if( matchArgs(allVars,analVars,m) )
    return 0 ; // change this to 1 if you want to noramlize by numerical integral

  return 0;
}

Double_t RooCBShape_MN::analyticalIntegral(Int_t code, const char* rangeName) const
{
 static const double sqrtPiOver2 = 1.2533141373;
 static const double sqrt2 = 1.4142135624;

 R__ASSERT(code==1);
 double result = 0.0;
 bool useLog = false;

 if( fabs(n-1.0) < 1.0e-05 )
   useLog = true;

 double sig = fabs((Double_t)sigma);

 double tmin = (m.min(rangeName)-m0)/sig;
 double tmax = (m.max(rangeName)-m0)/sig;

 if(alpha < 0) {
   double tmp = tmin;
   tmin = -tmax;
   tmax = -tmp;
 }

 double absAlpha = fabs((Double_t)alpha);

 if( tmin >= -absAlpha ) {
   result += sig*sqrtPiOver2*(   ApproxErf(tmax/sqrt2)
                               - ApproxErf(tmin/sqrt2) );
 }
 else if( tmax <= -absAlpha ) {
   double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
   double b = n/absAlpha - absAlpha;

   if(useLog) {
     result += a*sig*( log(b-tmin) - log(b-tmax) );
   }
   else {
     result += a*sig/(1.0-n)*(   1.0/(TMath::Power(b-tmin,n-1.0))
                               - 1.0/(TMath::Power(b-tmax,n-1.0)) );
   }
 }
 else {
   double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
   double b = n/absAlpha - absAlpha;

   double term1 = 0.0;
   if(useLog) {
     term1 = a*sig*(  log(b-tmin) - log(n/absAlpha));
   }
   else {
     term1 = a*sig/(1.0-n)*(   1.0/(TMath::Power(b-tmin,n-1.0))
                             - 1.0/(TMath::Power(n/absAlpha,n-1.0)) );
   }

   double term2 = sig*sqrtPiOver2*(   ApproxErf(tmax/sqrt2)
                                    - ApproxErf(-absAlpha/sqrt2) );


   result += term1 + term2;
 }

 return result != 0. ? result : 1.E-300;
}
