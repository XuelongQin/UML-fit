/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef PARBOUND
#define PARBOUND

#include <math.h>
#include "Math/SpecFunc.h"
#include "TMath.h"

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooFit.h"
#include "Riostream.h"
#include "RooObjCacheManager.h"
 
class ParBound : public RooAbsPdf {
 protected:

  RooRealProxy P1 ;
  RooRealProxy P2 ;
  RooRealProxy P3 ;
  RooRealProxy P4p ;
  RooRealProxy P5p ;
  RooRealProxy P6p ;
  RooRealProxy P8p ;

  Double_t evaluate() const ;

 public:
  bool verbose = false;
  int  Q2Bin  = 0;
  ParBound() {} ; 
  ParBound(const char *name, const char *title,
	   RooAbsReal& _P1,
	   RooAbsReal& _P2,
	   RooAbsReal& _P3,
	   RooAbsReal& _P4p,
	   RooAbsReal& _P5p,
	   RooAbsReal& _P6p,
	   RooAbsReal& _P8p);
  ParBound(const ParBound& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new ParBound(*this,newname); }
  inline virtual ~ParBound() { }
  
  ClassDef(ParBound,1) // PDF constraint for physical boundary on angular parameters
    };

#endif
