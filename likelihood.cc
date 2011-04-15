#include <complex>
#include <vector>
#include <iostream>

#include "TH1.h"
#include "TStopwatch.h"

#include "control.h"
#include "likelihood.h"

using namespace std;


likelihood::likelihood(waveset ws_,
		       vector<event>& RDevents_,
		       vector<event>& MCevents_,
		       vector<event>& MCallEvents_,
		       size_t nBins_, double threshold_, double binWidth_)
  : ws(ws_),
    RDevents(RDevents_),
    MCevents(MCevents_),
    MCallEvents(MCallEvents_),
    nBins(nBins_),
    threshold(threshold_),
    binWidth(binWidth_),
    currentBin(0)
{
  // Bin the data, once and for all.
  TStopwatch sw;
  sw.Start();

  binnedRDevents.resize(nBins);
  binnedMCevents.resize(nBins);
  binnedEtaAcc.resize(nBins);
  for (size_t iBin = 0; iBin < nBins; iBin++)
    {
      setBin(iBin);
      for (size_t iEvent = 0; iEvent < RDevents.size(); iEvent++)
	{
	  if (!RDevents[iEvent].accepted())
	    continue;
	  binnedRDevents[iBin].push_back(RDevents[iEvent]);
	}

      for (size_t iEvent = 0; iEvent < MCevents.size(); iEvent++)
	{
	  if (!flatMC && !MCevents[iEvent].accepted())
	    continue;
	  binnedMCevents[iBin].push_back(MCevents[iEvent]);
	}

      if (!flatMC)
	{
	  double countAllMC = 0;  // no of MC events generated in bin
	  for (size_t iEvent = 0; iEvent < MCallEvents.size(); iEvent++)
	    {
	      if (!MCallEvents[iEvent].accepted())
		continue;
	      countAllMC++;
	    }
	  binnedEtaAcc[iBin] = 1.*binnedMCevents[iBin].size() / countAllMC;
	}
      else
	{
	  binnedEtaAcc[iBin] = 1;  // No acceptance effects -> All accepted.
	}
    }

  sw.Stop();
  cout << "data binned after " << sw.CpuTime() << " s." << endl;
}


double
likelihood::probabilityDensity(const vector<double>& x, double theta, double phi) const
{
  return this->probabilityDensity(x, event(theta, phi));
}

double
likelihood::probabilityDensity(const vector<double>& x, const event& e) const
{
  double sum = 0;
  waveset::const_iterator it;
  size_t idx_base = 0;
  for (it = ws.begin(); it != ws.end(); it++)
    {
      // norm = abs^2
      sum += norm(it->sum(&x[idx_base], e));
      idx_base += 2*it->waves.size();
    }
  return sum;
}

double
likelihood::MCweight(int reflectivity, const wave& w1, const wave& w2) const
{
  int id = reflectivity + ((w1.l << 16) | (w1.m << 12)
			   | (w2.l << 8) | (w2.m << 4));
  if (weights.find(id) != weights.end())
    return weights[id];

  // Uses Kahan's summation
  double sum = 0;
  double c = 0;
  const vector<event>& pMCevents
    = flatMC ? MCevents : binnedMCevents[currentBin];
  for (size_t i = 0; i < pMCevents.size(); i++)
    {
      // Forming this sum as REAL sum and with no conjugate, because
      // the form of the decay amplitudes allows this.  This is not
      // the most general form!
      double y = pMCevents[i].MCweight(reflectivity,w1,w2) - c;
      double t = sum + y;
      c = (t - sum) - y;  // compensation term.
      sum = t;
    }

  /*
  cout << "calculated MCweight " << sum / pMCevents.size()
       << " from " << pMCevents.size() << " MC events for "
       << "(l1,m1,l2,m2) = "
       << "(" << w1.l << "," << w1.m << "," << w2.l << "," << w2.m << ")"
       << endl;
  */

  return weights[id] = sum / pMCevents.size();
}

double
likelihood::calc_mc_likelihood(const vector<double>& x) const
{
  double sumMC = 0;
  waveset::const_iterator it;
  size_t waveset_base = 0;
  for (it = ws.begin(); it != ws.end(); it++)
    {
      vector<wave>::const_iterator wave1, wave2;
      size_t idx1;
      for (wave1 = it->waves.begin(), idx1 = waveset_base;
	   wave1 != it->waves.end();
	   wave1++, idx1 += 2)
	{
	  complex<double> a1(x[idx1], x[idx1+1]);
	  size_t idx2;
	  for (wave2 = it->waves.begin(), idx2 = waveset_base;
	       wave2 != it->waves.end();
	       wave2++, idx2 += 2)
	    {
	      complex<double> conj_a2(x[idx2], -x[idx2+1]);
	      sumMC +=  real(a1*conj_a2
			     *MCweight(it->reflectivity, *wave1, *wave2));
	    }
	}
      waveset_base = idx1;
    }

  return sumMC * binnedEtaAcc[currentBin];
}

double
likelihood::calc_rd_likelihood(const vector<double>& x) const
{
  // Uses Kahan summation
  double sumRD = 0;
  double c = 0;
  const vector<event>& events = binnedRDevents[currentBin];
  for (size_t i = 0; i < events.size(); i++)
    {
      double y = log(probabilityDensity(x, events[i])) - c;
      double t = sumRD + y;
      c = (t - sumRD) - y;
      sumRD = t;
      //sumRD += y;
    }

  return sumRD;
}


double
likelihood::calc_likelihood(const vector<double>& x) const
{
  double lhMC = calc_mc_likelihood(x);
  double lhRD = calc_rd_likelihood(x);

  //cout << "likelihood = " << lhRD << " - " << lhMC << endl;
  return lhRD - lhMC;
}


double
likelihood::operator()(const vector<double>& x) const
{
  return -this->calc_likelihood(x);
}
