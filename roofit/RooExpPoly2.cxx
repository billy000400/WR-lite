/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooExpPoly2.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(RooExpPoly2); 

 RooExpPoly2::RooExpPoly2(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _a,
                        RooAbsReal& _b) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   a("a","a",this,_a),
   b("b","b",this,_b)
 { 
 } 


 RooExpPoly2::RooExpPoly2(const RooExpPoly2& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   a("a",this,other.a),
   b("b",this,other.b)
 { 
 } 



 Double_t RooExpPoly2::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   return exp(a*x*x+b*x) ; 
 } 



