#include "RooFit.h"

#include "Riostream.h"
#include <math.h>

#include "ExtendedCrystalBall.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "TMath.h"

ExtendedCrystalBall::ExtendedCrystalBall(const char* name, const char* title,
                                         RooAbsReal& _m, RooAbsReal& _m0,
                                         RooAbsReal& _sigma, RooAbsReal& _alphaL,
                                         RooAbsReal& _nL, RooAbsReal& _alphaR,
                                         RooAbsReal& _nR) :
  RooAbsPdf(name, title),
  m("m", "mass", this, _m),
  m0("m0", "pole", this, _m0),
  sigma("sigma", "Sigma", this, _sigma),
  alphaL("alphaL", "Alpha left tail", this, _alphaL),
  nL("nL", "Order left tail", this, _nL),
  alphaR("alphaR", "Alpha right tail", this, _alphaR),
  nR("nR", "Order right tail", this, _nR) {}

ExtendedCrystalBall::ExtendedCrystalBall(const ExtendedCrystalBall& other,
                                         const char* name) :
  RooAbsPdf(other, name),
  m("m", this, other.m),
  m0("m0", this, other.m0),
  sigma("sigma", this, other.sigma),
  alphaL("alphaL", this, other.alphaL),
  nL("nL", this, other.nL),
  alphaR("alphaR", this, other.alphaR),
  nR("nR", this, other.nR) {}

Double_t ExtendedCrystalBall::evaluate() const {
	Double_t t = (m - m0) / sigma;
	if (alphaL < 0 || alphaR < 0)
		t = -t;

	Double_t absAlphaL = fabs((Double_t)alphaL);
	Double_t absAlphaR = fabs((Double_t)alphaR);

	if (t > -absAlphaL && t < absAlphaR) {
		return exp(-0.5 * t * t);
	} else if (t <= -absAlphaL) {
		Double_t a =
		  TMath::Power(nL / absAlphaL, nL) * exp(-0.5 * absAlphaL * absAlphaL);
		Double_t b = nL / absAlphaL - absAlphaL;

		return a / TMath::Power(b - t, nL);
	} else if (t >= absAlphaR) {
		Double_t c =
		  TMath::Power(nR / absAlphaR, nR) * exp(-0.5 * absAlphaR * absAlphaR);
		Double_t d = nR / absAlphaR - absAlphaR;

		return c / TMath::Power(d + t, nR);
	} else
		return 0; // to avoid warnings
}
