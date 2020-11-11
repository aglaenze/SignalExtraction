
#ifndef EXTENDEDCRYSTALBALL
#define EXTENDEDCRYSTALBALL

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class ExtendedCrystalBall : public RooAbsPdf {
 public:
	ExtendedCrystalBall(){};
	ExtendedCrystalBall(const char* name, const char* title, RooAbsReal& _m,
	                    RooAbsReal& _m0, RooAbsReal& _sigma, RooAbsReal& _alphaL,
	                    RooAbsReal& _nL, RooAbsReal& _alphaR, RooAbsReal& _nR);

	ExtendedCrystalBall(const ExtendedCrystalBall& other, const char* name = 0);
	virtual TObject* clone(const char* newname) const {
		return new ExtendedCrystalBall(*this, newname);
	}

	inline virtual ~ExtendedCrystalBall() {}

 protected:
	RooRealProxy m;
	RooRealProxy m0;
	RooRealProxy sigma;
	RooRealProxy alphaL;
	RooRealProxy nL;
	RooRealProxy alphaR;
	RooRealProxy nR;

	Double_t evaluate() const;

 private:
	ClassDef(ExtendedCrystalBall, 1)
};

#endif
