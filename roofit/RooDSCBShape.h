/**
 * @Author: Billy Li <billyli>
 * @Date:   04-26-2022
 * @Email:  li000400@umn.edu
 * @Last modified by:   billyli
 * @Last modified time: 04-26-2022
 */



/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

#ifndef ROODSCBSHAPE
#define ROODSCBSHAPE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooDSCBShape : public RooAbsPdf {
public:
  RooDSCBShape() {} ;
  RooDSCBShape(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _mu,
	      RooAbsReal& _sigma,
	      RooAbsReal& _alphaL,
	      RooAbsReal& _alphaR,
	      RooAbsReal& _nL,
	      RooAbsReal& _nR);
  RooDSCBShape(const RooDSCBShape& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooDSCBShape(*this,newname); }
  inline virtual ~RooDSCBShape() { }

protected:

  RooRealProxy x ;
  RooRealProxy mu ;
  RooRealProxy sigma ;
  RooRealProxy alphaL ;
  RooRealProxy alphaR ;
  RooRealProxy nL ;
  RooRealProxy nR ;

  Double_t evaluate() const ;

private:

  ClassDef(RooDSCBShape,1) // Your description goes here...
};

#endif
