/**
 * @Author: Billy Li <billyli>
 * @Date:   05-13-2022
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 06-23-2022
 */



/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

// Your description goes here...

#include "Riostream.h"

#include "RooExpm.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"

ClassImp(RooExpm);

 RooExpm::RooExpm(const char *name, const char *title,
                        RooAbsReal& _x,
                        RooAbsReal& _a,
                        RooAbsReal& _b) :
   RooAbsPdf(name,title),
   x("x","x",this,_x),
   a("a","a",this,_a),
   b("b","b",this,_b)
 {
 }


 RooExpm::RooExpm(const RooExpm& other, const char* name) :
   RooAbsPdf(other,name),
   x("x",this,other.x),
   a("a",this,other.a),
   b("b",this,other.b)
 {
 }



 Double_t RooExpm::evaluate() const
 {
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE
   return exp(a*TMath::Power(x,b))+1e-28 ; 
 }
