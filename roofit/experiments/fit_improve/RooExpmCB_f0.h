/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

#ifndef ROOEXPMCB
#define ROOEXPMCB

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooExpmCB_f0 : public RooAbsPdf {
public:
  RooExpmCB_f0() {} ;
  RooExpmCB_f0(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _mu,
	      RooAbsReal& _sigma,
	      RooAbsReal& _beta,
	      RooAbsReal& _m,
	      RooAbsReal& _alpha,
	      RooAbsReal& _n);
  RooExpmCB_f0(const RooExpmCB_f0& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooExpmCB_f0(*this,newname); }
  inline virtual ~RooExpmCB_f0() { }

protected:

  RooRealProxy x ;
  RooRealProxy mu ;
  RooRealProxy sigma ;
  RooRealProxy beta ;
  RooRealProxy m ;
  RooRealProxy alpha ;
  RooRealProxy n ;

  Double_t evaluate() const ;

private:

  ClassDef(RooExpmCB_f0,1) // Your description goes here...
};

#endif
