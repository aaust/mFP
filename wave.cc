#include <stdio.h>
#include <iostream>

#include "TFitterMinuit.h"
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


// The covariance matrix doesn't contain columns corresponding to the
// fixed parameters, and there seems to be no way foreseen to obtain
// the correct indices into the covariance matrix.  Therefore this
// contortion ... which carries the assumption also made above that we
// never fix real parts.
size_t
wave::idxInCovariance(const TFitterMinuit* minuit) const
{
  size_t countFixedBelow = 0;
  for (size_t i = 0; i < idx; i++)
    countFixedBelow += minuit->IsFixed(i);
  return idx - countFixedBelow;
}

  
void
wave::fillHistIntensity(int iBin, const TFitterMinuit* minuit)
{
  complex<double> a(minuit->GetParameter(idx), minuit->GetParameter(idx+1));

  histIntensity->SetBinContent(iBin+1,norm(a));
  double error;
  if (minuit->IsFixed(idx+1))
    error = 2*abs(a)*minuit->GetParError(idx);
  else
    {
      size_t idxCov = idxInCovariance(minuit);
      double cov = minuit->GetCovarianceMatrixElement(idxCov, idxCov + 1);

      error = 2*(sqrt(pow(real(a) * minuit->GetParError(idx), 2)
		      + pow(imag(a) * minuit->GetParError(idx+1), 2)
		      + (2*real(a)*imag(a)*cov)));
    }
  histIntensity->SetBinError(iBin+1, error);
}  


void
wave::fillHistPhase(int iBin, const wave& other, const TFitterMinuit* minuit)
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

  complex<double> a1(minuit->GetParameter(idx), minuit->GetParameter(idx + 1));
  complex<double> a2(minuit->GetParameter(other.getIndex()), minuit->GetParameter(other.getIndex()+1));

  double phi = arg(a1 / a2);
  if (iBin > 0)
    {
      double oldPhase = h->GetBinContent(iBin);
      while (oldPhase - phi > M_PI)
	phi += 2*M_PI;
      while (phi - oldPhase > M_PI)
	phi -= 2*M_PI;
    }
  h->SetBinContent(iBin+1, phi);

  // Put together the covariance matrix.  There are several cases:
  // 1. Im(a1) fixed.
  // 2. Im(a2) fixed.
  // 3. Neither fixed.

  if (minuit->IsFixed(idx + 1) || minuit->IsFixed(other.getIndex() + 1))
    {
      if (minuit->IsFixed(idx + 1) && minuit->IsFixed(other.getIndex() + 1))
	{
	  // Can't happen, ignore.
	  h->SetBinError(iBin+1, 0.2);
	}
      else
	{
	  bool swapped = false;
	  size_t idxFix = idx;
	  size_t idxFree = other.getIndex();
	  if (minuit->IsFixed(other.getIndex() + 1))
	    {
	      swapped = true;
	      std::swap(idxFix, idxFree);
	    }
	  TMatrixDSym cov(3);
	  size_t idxCovFix, idxCovFree;
	  if (swapped)
	    {
	      idxCovFix = other.idxInCovariance(minuit);
	      idxCovFree = idxInCovariance(minuit);
	    }
	  else
	    {
	      idxCovFix = idxInCovariance(minuit);
	      idxCovFree = other.idxInCovariance(minuit);
	    }	  

	  cov(0,0) = minuit->GetCovarianceMatrixElement(idxCovFix, idxCovFix);
	  cov(1,1) = minuit->GetCovarianceMatrixElement(idxCovFree, idxCovFree);
	  cov(2,2) = minuit->GetCovarianceMatrixElement(idxCovFree + 1, idxCovFree + 1);

	  cov(0,1) = cov(1,0) = minuit->GetCovarianceMatrixElement(idxCovFix, idxCovFree);
	  cov(0,2) = cov(2,0) = minuit->GetCovarianceMatrixElement(idxCovFix, idxCovFree + 1);
	  cov(1,2) = cov(2,1) = minuit->GetCovarianceMatrixElement(idxCovFree, idxCovFree + 1);

	  double values[3] = { minuit->GetParameter(idxFix),
			       minuit->GetParameter(idxFree),
			       minuit->GetParameter(idxFree + 1) };
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

	  h->SetBinError(iBin+1, std::min(3.142,sqrt(cov.Similarity(gradient))));
	}
    }
  else
    {   
      // First, for the a1, a2:
      TMatrixDSym cov(4);
      size_t idxCov, idxCovOther;
      idxCov = idxInCovariance(minuit);
      idxCovOther = other.idxInCovariance(minuit);
      cov(0,0) = minuit->GetCovarianceMatrixElement(idxCov, idxCov);
      cov(1,1) = minuit->GetCovarianceMatrixElement(idxCov + 1, idxCov + 1);
      cov(2,2) = minuit->GetCovarianceMatrixElement(idxCovOther, idxCovOther);
      cov(3,3) = minuit->GetCovarianceMatrixElement(idxCovOther + 1, idxCovOther + 1);

      cov(0,1) = cov(1,0) = minuit->GetCovarianceMatrixElement(idxCov, idxCov + 1);
      cov(0,2) = cov(2,0) = minuit->GetCovarianceMatrixElement(idxCov, idxCovOther);
      cov(0,3) = cov(3,0) = minuit->GetCovarianceMatrixElement(idxCov, idxCovOther + 1);
      cov(1,2) = cov(2,1) = minuit->GetCovarianceMatrixElement(idxCov + 1, idxCovOther);
      cov(1,3) = cov(3,1) = minuit->GetCovarianceMatrixElement(idxCov + 1, idxCovOther + 1);
      cov(2,3) = cov(3,2) = minuit->GetCovarianceMatrixElement(idxCovOther, idxCovOther + 1);

      double values[4] = { minuit->GetParameter(idx), minuit->GetParameter(idx + 1),
			   minuit->GetParameter(other.getIndex()), minuit->GetParameter(other.getIndex()+1) };
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

      h->SetBinError(iBin+1, std::min(3.142,sqrt(cov.Similarity(gradient))));
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
