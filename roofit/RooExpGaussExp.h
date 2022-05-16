/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOEXPGAUSSEXP
#define ROOEXPGAUSSEXP

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooExpGaussExp : public RooAbsPdf {
public:
  RooExpGaussExp() {} ; 
  RooExpGaussExp(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _mu,
	      RooAbsReal& _sigma,
	      RooAbsReal& _m,
	      RooAbsReal& _n);
  RooExpGaussExp(const RooExpGaussExp& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooExpGaussExp(*this,newname); }
  inline virtual ~RooExpGaussExp() { }

protected:

  RooRealProxy x ;
  RooRealProxy mu ;
  RooRealProxy sigma ;
  RooRealProxy m ;
  RooRealProxy n ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooExpGaussExp,1) // Your description goes here...
};
 
#endif
