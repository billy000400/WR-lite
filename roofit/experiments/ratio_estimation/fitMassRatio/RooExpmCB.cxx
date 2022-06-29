/**
 * @Author: Billy Li <billyli>
 * @Date:   05-05-2022
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 06-03-2022
 */



/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

// Your description goes here...

#include "Riostream.h"

#include "RooExpmCB.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"

ClassImp(RooExpmCB);

 RooExpmCB::RooExpmCB(const char *name, const char *title,
                        RooAbsReal& _x,
                        RooAbsReal& _mu,
                        RooAbsReal& _sigma,
                        RooAbsReal& _beta,
                        RooAbsReal& _m,
                        RooAbsReal& _alpha,
                        RooAbsReal& _n) :
   RooAbsPdf(name,title),
   x("x","x",this,_x),
   mu("mu","mu",this,_mu),
   sigma("sigma","sigma",this,_sigma),
   beta("beta","beta",this,_beta),
   m("m","m",this,_m),
   alpha("alpha","alpha",this,_alpha),
   n("n","n",this,_n)
 {
 }


 RooExpmCB::RooExpmCB(const RooExpmCB& other, const char* name) :
   RooAbsPdf(other,name),
   x("x",this,other.x),
   mu("mu",this,other.mu),
   sigma("sigma",this,other.sigma),
   beta("beta",this,other.beta),
   m("m",this,other.m),
   alpha("alpha",this,other.alpha),
   n("n",this,other.n)
 {
 }



 Double_t RooExpmCB::evaluate() const
 {
    Double_t t = (x-mu)/sigma;

    Double_t absAlpha = fabs((Double_t)alpha);
    Double_t absBeta = fabs((Double_t)beta);

    Double_t result =0;

    if (t < -absBeta){
      Double_t m_inv = 1/m;
      Double_t A = exp(absBeta*absBeta*(m_inv-0.5));
      Double_t omega = m_inv*TMath::Power(absBeta, 2-m);
      Double_t abs_t = fabs(t);
      result = A*exp(-omega*TMath::Power(abs_t,m));

      // handle nan (assume exp gives a near 0 number)
      // if (result != result){
      //     result = 0.;
      // }

    } else if (t < absAlpha) {

      result = exp(-0.5*t*t);

    }
    else {
      // Double_t a =  TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
      // Double_t b= n/absAlpha - absAlpha;

      // return a/TMath::Power(b - t, n);

      Double_t exp_part = exp(-0.5*absAlpha*absAlpha);
      Double_t D = (n-absAlpha*absAlpha+absAlpha*t);
      Double_t arg = n/D;

      result = TMath::Power(arg, n)*exp_part;
   }
   result += 1e-36;
   return result;
 }
