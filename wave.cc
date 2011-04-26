#include <stdio.h>
#include <iostream>

#include "TFitterMinuit.h"

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
      // The covariance matrix doesn't contain columns corresponding
      // to the fixed parameters, and there seems to be no way
      // foreseen to obtain the correct indices into the covariance
      // matrix.  Therefore this contortion ... which carries the
      // assumption also made above that we never fix real parts.
      int countFixedBelow = 0;
      for (size_t i = 0; i < idx; i++)
	countFixedBelow += minuit->IsFixed(i);

      double cov = minuit->GetCovarianceMatrixElement(idx - countFixedBelow,
						      idx + 1 - countFixedBelow);

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
  h->SetBinError(iBin+1, 0.2);
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
