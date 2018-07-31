#include <stdio.h>
#include <iostream>

#include "Math/Minimizer.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include "wave.h"
#include "event.h"

wave::wave(const wave& o)
{
  l = o.l;
  m = o.m;
  name = o.name;
  idx = o.idx;
  histIntensity = o.histIntensity;
  mHistPhase = o.mHistPhase;
  phaseLocked = o.phaseLocked;
}


void
wave::buildHists(int nBins, double lower, double upper)
{
  char histName[999];
  char title[999];

  snprintf(histName, 999, "hIntensity%s", name.c_str());
  snprintf(title, 999, "Intensity of %s wave", name.c_str());
  histIntensity = new TH1D(histName, title, nBins, lower, upper);
}


TH1*
wave::getHistPhase(const wave& other)
{
  map<string, TH1*>::const_iterator it = mHistPhase.find(other.name);
  if (it == mHistPhase.end())
    return 0;
  else
    return it->second;
}


void
wave::fillHistIntensity(int iBin, const ROOT::Math::Minimizer* minuit)
{
  complex<double> a(minuit->X()[idx], minuit->X()[idx+1]);

  histIntensity->SetBinContent(iBin+1,norm(a));
  double error;
  if (minuit->IsFixedVariable(idx+1))
    error = 2*abs(a)*minuit->Errors()[idx];
  else
    {
      double cov = minuit->CovMatrix(idx, idx + 1);

      error = 2*(sqrt(pow(real(a) * minuit->Errors()[idx], 2)
                      + pow(imag(a) * minuit->Errors()[idx+1], 2)
      		      + (2*real(a)*imag(a)*cov)));
    }
  histIntensity->SetBinError(iBin+1, error);
}  


void
wave::fillHistPhase(int iBin, const wave& other, const ROOT::Math::Minimizer* minuit)
{
  TH1* h = 0;
  if (mHistPhase.find(other.name) == mHistPhase.end())
    {
      char histName[999];
      char title[999];

      snprintf(histName, 999, "hPhase%s%s", name.c_str(), other.name.c_str());
      snprintf(title, 999, "arg(%s / %s)", name.c_str(), other.name.c_str());

      int nBins = histIntensity->GetNbinsX();
      double lower = histIntensity->GetXaxis()->GetXmin();
      double upper = histIntensity->GetXaxis()->GetXmax();
      mHistPhase[other.name] = new TH1D(histName, title, nBins, lower, upper);
    }

  h = mHistPhase[other.name];

  complex<double> a1(minuit->X()[idx], minuit->X()[idx + 1]);
  complex<double> a2(minuit->X()[other.getIndex()], minuit->X()[other.getIndex()+1]);

  double phi = arg(a1 / a2);
  double oldPhase = 0;
  if (iBin > 0)
    {
      oldPhase = h->GetBinContent(iBin);
      // but sometimes a number of previous fits failed
      int fail=1;
      while (oldPhase == 0 && iBin-fail>0)
	{
	  oldPhase = h->GetBinContent(iBin-fail);
	  fail++;
	}
    }
  while ((oldPhase - phi) > M_PI)
    phi += 2*M_PI;
  while ((phi - oldPhase) > M_PI)
    phi -= 2*M_PI;
  
  h->SetBinContent(iBin+1, phi);

  // Put together the covariance matrix.  There are several cases:
  // 1. Im(a1) fixed.
  // 2. Im(a2) fixed.
  // 3. Neither fixed.

  if (minuit->IsFixedVariable(idx + 1) || minuit->IsFixedVariable(other.getIndex() + 1))
    {
      if (minuit->IsFixedVariable(idx + 1) && minuit->IsFixedVariable(other.getIndex() + 1))
	{
	  // Can't happen, ignore.
	  h->SetBinError(iBin+1, 0.2);
	}
      else
	{
	  size_t idxFix = idx;
	  size_t idxFree = other.getIndex();
	  if (minuit->IsFixedVariable(other.getIndex() + 1))
	    std::swap(idxFix, idxFree);
	    
	  TMatrixDSym cov(3);

	  cov(0,0) = minuit->CovMatrix(idxFix, idxFix);
	  cov(1,1) = minuit->CovMatrix(idxFree, idxFree);
	  cov(2,2) = minuit->CovMatrix(idxFree + 1, idxFree + 1);

	  cov(0,1) = cov(1,0) = minuit->CovMatrix(idxFix, idxFree);
	  cov(0,2) = cov(2,0) = minuit->CovMatrix(idxFix, idxFree + 1);
	  cov(1,2) = cov(2,1) = minuit->CovMatrix(idxFree, idxFree + 1);

	  double values[3] = { minuit->X()[idxFix],
			       minuit->X()[idxFree],
			       minuit->X()[idxFree + 1] };
	  TVectorD gradient(3);
	  // The algorithm follows the one in TF1::Derivative() :
	  //   df(x) = (4 D(h/2) - D(h)) / 3
	  // with D(h) = (f(x + h) - f(x - h)) / (2 h).
	  for (size_t i = 0; i < 3; i++)
	    {
	      double save = values[i];
	      const double step = 1e-5;

	      values[i] = save + step/2;
	      complex<double> a1(values[0], 0);
	      complex<double> a2(values[1], values[2]);

	      double rightShort = arg(a1 / a2);

	      values[i] = save - step/2;
	      a1 = complex<double>(values[0], 0);
	      a2 = complex<double>(values[1], values[2]);

	      double leftShort = arg(a1 / a2);
	      while(rightShort - leftShort > M_PI)
		leftShort += 2*M_PI;
	      while(leftShort - rightShort > M_PI)
		rightShort += 2*M_PI;

	      values[i] = save + step;
	      a1 = complex<double>(values[0], 0);
	      a2 = complex<double>(values[1], values[2]);

	      double rightLong = arg(a1 / a2);

	      values[i] = save - step;
	      a1 = complex<double>(values[0], 0);
	      a2 = complex<double>(values[1], values[2]);

	      double leftLong = arg(a1 / a2);
	      while(rightLong - leftLong > M_PI)
		leftLong += 2*M_PI;
	      while(leftLong - rightLong > M_PI)
		rightLong += 2*M_PI;

	      double derivFull = (rightLong - leftLong) / 2 / step;
	      double derivShort = (rightShort - leftShort) / step;

	      gradient[i] =  1./3.*(4*derivShort - derivFull);
	    }

	  h->SetBinError(iBin+1, std::min(M_PI,sqrt(cov.Similarity(gradient))));
	}
    }
  else
    {   
      TMatrixDSym cov(4);
      size_t idxOther = other.getIndex();
      cov(0,0) = minuit->CovMatrix(idx, idx);
      cov(1,1) = minuit->CovMatrix(idx + 1, idx + 1);
      cov(2,2) = minuit->CovMatrix(idxOther, idxOther);
      cov(3,3) = minuit->CovMatrix(idxOther + 1, idxOther + 1);

      cov(0,1) = cov(1,0) = minuit->CovMatrix(idx, idx + 1);
      cov(0,2) = cov(2,0) = minuit->CovMatrix(idx, idxOther);
      cov(0,3) = cov(3,0) = minuit->CovMatrix(idx, idxOther + 1);
      cov(1,2) = cov(2,1) = minuit->CovMatrix(idx + 1, idxOther);
      cov(1,3) = cov(3,1) = minuit->CovMatrix(idx + 1, idxOther + 1);
      cov(2,3) = cov(3,2) = minuit->CovMatrix(idxOther, idxOther + 1);

      double values[4] = { minuit->X()[idx], minuit->X()[idx + 1],
			   minuit->X()[other.getIndex()], minuit->X()[other.getIndex()+1] };
      TVectorD gradient(4);
      // The algorithm follows the one in TF1::Derivative() :
      //   df(x) = (4 D(h/2) - D(h)) / 3
      // with D(h) = (f(x + h) - f(x - h)) / (2 h).
      for (size_t i = 0; i < 4; i++)
	{
	  double save = values[i];
	  const double step = 1e-5;

	  values[i] = save + step/2;
	  complex<double> a1(values[0], values[1]);
	  complex<double> a2(values[2], values[3]);

	  double rightShort = arg(a1 / a2);

	  values[i] = save - step/2;
	  a1 = complex<double>(values[0], values[1]);
	  a2 = complex<double>(values[2], values[3]);

	  double leftShort = arg(a1 / a2);
	  while(rightShort - leftShort > M_PI)
	    leftShort += 2*M_PI;
	  while(leftShort - rightShort > M_PI)
	    rightShort += 2*M_PI;

	  values[i] = save + step;
	  a1 = complex<double>(values[0], values[1]);
	  a2 = complex<double>(values[2], values[3]);

	  double rightLong = arg(a1 / a2);

	  values[i] = save - step;
	  a1 = complex<double>(values[0], values[1]);
	  a2 = complex<double>(values[2], values[3]);

	  double leftLong = arg(a1 / a2);
	  while(rightLong - leftLong > M_PI)
	    leftLong += 2*M_PI;
	  while(leftLong - rightLong > M_PI)
	    rightLong += 2*M_PI;

	  double derivFull = (rightLong - leftLong) / 2 / step;
	  double derivShort = (rightShort - leftShort) / step;

	  gradient[i] =  1./3.*(4*derivShort - derivFull);
	}

      h->SetBinError(iBin+1, std::min(M_PI,sqrt(cov.Similarity(gradient))));
    }

}


complex<double>
coherent_waves::sum(const vector<double>& x, const event& e) const
{
  complex<double> result = 0;
  vector<wave>::const_iterator it;
  for (it = this->waves.begin(); it != this->waves.end(); it++)
    {
      complex<double> a(x[it->getIndex()], x[it->getIndex()+1]);
      result += a * e.decayAmplitude(this->reflectivity,*it);
    }

  return result;
}


double
coherent_waves::getEventWeight(const vector<double>& x, const event& e) const
{
  double weight = 0;
  vector<wave>::const_iterator wave1, wave2;
  for (wave1 = this->waves.begin();
       wave1 != this->waves.end();
       wave1++)
    {
      complex<double> a1(x[wave1->getIndex()], x[wave1->getIndex()+1]);
      for (wave2 = this->waves.begin();
	   wave2 != this->waves.end();
	   wave2++)
	{
	  complex<double> conj_a2(x[wave2->getIndex()], -x[wave2->getIndex()+1]);
	  weight += real(a1 * conj_a2
			 * e.decayAmplitude(this->reflectivity, *wave1)
			 * e.decayAmplitude(this->reflectivity, *wave2));
	}
    }
  return 4 * M_PI * weight;
}

waveset::waveset() { return; }

ClassImp(wave)
ClassImp(coherent_waves)
ClassImp(waveset)
