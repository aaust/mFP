#include <complex>
#include <vector>
#include <iostream>

#include "TH1.h"
#include "TStopwatch.h"

#include "control.h"
#include "likelihood.h"

using namespace std;


likelihood::likelihood(const waveset& ws_,
		       eventStream* RDevents_,
		       eventStream* MCevents_,
		       size_t nBins_, double threshold_, double binWidth_,
		       size_t idxBranching_)
  : ws(ws_),
    RDevents(RDevents_),
    MCevents(MCevents_),
    nBins(nBins_),
    threshold(threshold_),
    binWidth(binWidth_),
    idxBranching(idxBranching_),
    currentBin(0)
{
  // Set the indices for the respective waves.
  lastIdx = 0;
  for (waveset::iterator it = ws.begin(); it != ws.end(); it++)
    {
      for (vector<wave>::iterator itWave = it->getWaves().begin();
	   itWave != it->getWaves().end(); itWave++)
	{
	  itWave->setIndex(lastIdx);
	  lastIdx += 2;
	}
    }
}

void
likelihood::calcMCweights()
{
  cout << "calculating MCweight" << flush;
  TStopwatch sw;
  sw.Start();

  // Uses Kahan's summation
  MCweights.ResizeTo(lastIdx / 2, lastIdx / 2);
  TMatrixD sum(lastIdx / 2, lastIdx / 2);
  TMatrixD c(lastIdx / 2, lastIdx / 2);
  eventsAccepted = 0;
  eventsInBin = 0;
  for (Long_t i = 0; i < MCevents->nEvents(); i++)
    {
      // Forming this sum as REAL sum and with no conjugate, because
      // the form of the decay amplitudes allows this.  This is not
      // the most general form!
      const event& e = MCevents->getEvent(i);
      if (!e.inBin())
	continue;
      eventsInBin++;
      if (!e.accepted())
	continue;
      eventsAccepted++;

      for (waveset::iterator it = ws.begin(); it != ws.end(); it++)
	{
	  for (vector<wave>::iterator itWave1 = it->getWaves().begin();
	       itWave1 != it->getWaves().end(); itWave1++)
	    {
	      size_t idx1 = itWave1->getIndex() / 2;
	      for (vector<wave>::iterator itWave2 = it->getWaves().begin();
		   itWave2 != it->getWaves().end(); itWave2++)
		{
		  size_t idx2 = itWave2->getIndex() / 2;
		  double y = (e.MCweight(it->reflectivity,*itWave1,*itWave2)
			      - c(idx1, idx2));
		  double t = sum(idx1, idx2) + y;
		  //c(idx1, idx2) = (t - sum(idx1, idx2)) - y;  // compensation term.
		  sum(idx1,idx2) = t;
		}
	    }
	}
    }

  MCweights = sum;
  MCweights *= 1./eventsInBin;
  sw.Stop();
  cout << " done after " << sw.CpuTime() << "s CPU, " << sw.RealTime() << "s wall time" << endl;
  MCweights.Print();
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
  for (it = ws.begin(); it != ws.end(); it++)
    {
      // norm = abs^2
      sum += norm(it->sum(x, e));
    }

  return (idxBranching == 16 ? 1-x[idxBranching+1]:x[idxBranching])*sum;
}

double
likelihood::MCweight(int reflectivity, const wave& w1, const wave& w2) const
{
  return MCweights(w1.getIndex() / 2, w2.getIndex() / 2);

  int id = reflectivity + ((w1.l << 16) | (w1.m << 12)
			   | (w2.l << 8) | (w2.m << 4));
  if (weights.find(id) != weights.end())
    return weights[id];

  // Uses Kahan's summation
  double sum = 0;
  double c = 0;
  eventsAccepted = 0;
  eventsInBin = 0;
  for (Long_t i = 0; i < MCevents->nEvents(); i++)
    {
      // Forming this sum as REAL sum and with no conjugate, because
      // the form of the decay amplitudes allows this.  This is not
      // the most general form!
      const event& e = MCevents->getEvent(i);
      if (!e.inBin())
	continue;
      eventsInBin++;
      if (!e.accepted())
	continue;
      eventsAccepted++;
      double y = e.MCweight(reflectivity,w1,w2) - c;
      double t = sum + y;
      c = (t - sum) - y;  // compensation term.
      sum = t;
    }

  ///*
  cout << "calculated MCweight " << sum
       << " from " << MCevents->nEvents() << " MC events for "
       << "(l1,m1,l2,m2) = "
       << "(" << w1.l << "," << w1.m << "," << w2.l << "," << w2.m << ")"
       << endl;
//  */

  return weights[id] = sum;
}

#if 0
complex<double>
likelihood::MCmomentWeight(int L1, int M1, int L2, int M2)
{

  double spherical = ROOT::Math::sph_legendre(l, m, this->theta);

  // Uses Kahan summation
  complex<double> sum = 0;
  complex<double> comp = 0;
  const vector<event>& events = binnedMCevents[currentBin];
  for (size_t i = 0; i < events.size(); i++)
    {
      double YL1M1 = ROOT::Math::sph_legendre(L1, M1, events[i].theta);
      double YL2M2 = ROOT::Math::sph_legendre(L2, M2, events[i].theta);

      double phi = events[i].phi;
      double c = cos((M2-M1)*phi);
      double s = sin((M2-M1)*phi);
     
      complex<double> y = YL1M1*YL2M2*complex<double>(c,s) - comp;
      complex<double> t = sum + y;
      comp = (t - sum) - y;
      sum = t;
      //sumRD += y;
    }

  return (2*L2+1) / (4*M_PI) * binnedEtaAcc[currentBin] / events.size() * sum;
}
#endif


double
likelihood::calc_mc_likelihood(const vector<double>& x) const
{
  double sumMC = 0;
  waveset::const_iterator it;
  for (it = ws.begin(); it != ws.end(); it++)
    {
      vector<wave>::const_iterator wave1, wave2;
      for (wave1 = it->waves.begin();
	   wave1 != it->waves.end();
	   wave1++)
	{
	  complex<double> a1(x[wave1->getIndex()], x[wave1->getIndex()+1]);
	  for (wave2 = it->waves.begin();
	       wave2 != it->waves.end();
	       wave2++)
	    {
	      complex<double> conj_a2(x[wave2->getIndex()], -x[wave2->getIndex()+1]);
	      sumMC +=  real(a1*conj_a2
			     *MCweight(it->reflectivity, *wave1, *wave2));
	    }
	}
    }

  return  (idxBranching == 16 ? 1-x[idxBranching+1]:x[idxBranching]) * sumMC * eventsAccepted / eventsInBin;
}

double
likelihood::calc_rd_likelihood(const vector<double>& x) const
{
  // Uses Kahan summation
  double sumRD = 0;
  double c = 0;
  for (Long_t i = 0; i < RDevents->nEvents(); i++)
    {
      const event& e = RDevents->getEvent(i);
      if (!e.inBin())
	continue;
      double y = log(probabilityDensity(x, e)) - c;
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

  cout /*<< "nevents = " << binnedRDevents[currentBin].size()*/ << " likelihood = " << lhRD << " - " << lhMC << endl;
  return lhRD - lhMC;
}


double
likelihood::operator()(const vector<double>& x) const
{
  return -this->calc_likelihood(x);
}


complex<double>
likelihood::calcMoment(int L, int M) const
{
  // Uses Kahan's summation
  complex<double> sum = 0;
  complex<double> c = 0;
  for (Long_t i = 0; i < RDevents->nEvents(); i++)
    {
      const event& e = RDevents->getEvent(i);
      complex<double> y = e.momentWeight(L, M) - c;
      complex<double> t = sum + y;
      c = (t - sum) - y;  // compensation term.
      sum = t;
    }
  return sum;
}
