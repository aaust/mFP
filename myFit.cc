#include <complex>
#include <vector>
#include <map>

using namespace std;

#include <string>
#include <stdio.h>

#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TRandom1.h"
#include "Math/SpecFunc.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"
#include "TStopwatch.h"
#include "TFile.h"

#include "control.h"

bool flatMC = false;
#define NFLATMCEVENTS 100000
double massLow = 0;
double massHigh = 9999;
int iBin;


#define NWAVES 8

struct wave {
  string name;
  int l, m;

  wave() { }
  wave(int ll, int mm) { l = ll; m = mm; }
  wave(const char* name_, int ll, int mm) : name(name_), l(ll), m(mm) {}
  wave(const wave& o) { l = o.l; m = o.m; name = o.name; }
};

class event;

struct coherent_waves {
  int reflectivity;
  int spinflip;
  vector<wave> waves;

  coherent_waves() {}
  coherent_waves(const coherent_waves& o) { reflectivity = o.reflectivity; spinflip = o.spinflip; waves = o.waves; }

  complex<double> sum(const double* x, const event& e) const;
  size_t getNwaves() { return waves.size(); }
};

typedef vector<coherent_waves> waveset;

class event {
public:
  double mass;
  double tPrime;
  double theta;
  double phi;

  event() { mass = tPrime = theta = phi = 0; }
  /*
  event(const event& other) { theta = other.theta; phi = other.phi; }
  */
  event(double th, double ph) { mass = tPrime = 0; theta = th; phi = ph; }
  event(double mass_, double tPrime_, double th, double ph)
  { mass = mass_; tPrime = tPrime_; theta = th; phi = ph; }
  event(const event& o) { mass = o.mass; tPrime = o.tPrime; theta = o.theta; phi = o.phi; }

  double decayAmplitude(int reflectivity, const wave& w) const
  { return decayAmplitude(reflectivity, w.l, w.m); };
  double decayAmplitude(int reflectivity, int l, int m) const;
  double MCweight(int reflectivity, const wave& w1, const wave& w2) const
  { return MCweight(reflectivity, w1.l, w1.m, w2.l, w2.m); }
  double MCweight(int reflectivity, int l1, int m1, int l2, int m2) const;

  bool accepted() const {
    return (this->mass >= massLow && this->mass < massHigh
	    //&& this->tPrime > 0.1 && this->tPrime < 0.3);
	    //&& this->tPrime > 0.3);
            && 1);
  }
};


class likelihood : public ROOT::Minuit2::FCNBase {
  vector<coherent_waves> ws;
  vector<event> RDevents;
  vector<event> MCevents;
  vector<event> MCallEvents;

  vector<vector<event> > binnedRDevents;
  vector<vector<event> > binnedMCevents;
  vector<double> binnedEtaAcc;

  size_t nBins;
  double threshold;
  double binWidth;

  size_t currentBin;
public:
  likelihood(waveset ws_,
	     vector<event>& RDevents_,
	     vector<event>& MCevents_,
	     vector<event>& MCallEvents_,
	     size_t nBins_, double threshold_, double binWidth_);

  double Up() const { return 0.5; }

public:
  double
  decay(int reflectivity, int l, int m, double theta, double phi) const;

  double
  probabilityDensity(const vector<double>& x, double theta, double phi) const;

  double
  probabilityDensity(const vector<double>& x, const event& e) const;

  double
  MCweight(int reflectivity, const wave& w1, const wave& w2) const;

  double
  calc_mc_likelihood(const vector<double>& x) const;

  double
  calc_rd_likelihood(const vector<double>& x) const;

  double
  calc_likelihood(const vector<double>& x) const;

  double
  operator()(const vector<double>& x) const;

  void
  setBin(size_t iBin) { currentBin = iBin;
    massLow = threshold + iBin*binWidth;
    massHigh = threshold + (iBin+1)*binWidth;
  }
};


// Returns the real / imaginary part of the amplitude, depending on which one
// is non-zero.

// The amplitude is:  Ylm(theta, phi) - epsilon (-)^m Yl-m(theta, phi)
//                  = Ylm(theta, phi) - epsilon Ylm(theta, phi)*
//                  = N Plm(theta) (e^(i m phi) - epsilon e^(-i m phi))
//                  = N' Plm(theta) {cos,sin}(m phi)

// NOTE the phase is ignored, as different reflectivities don't interfere

double
event::decayAmplitude(int reflectivity, int l, int m) const
{
  double spherical = ROOT::Math::sph_legendre(l, m, this->theta);

  // This absorbs the factor 2 from e^i \pm e^-i = 2 {i sin, cos}
  double factor = 1; //sqrt(nrdevents);
  if (m != 0)
    factor *= sqrt(2.);    

  if (reflectivity == +1)
    return factor*spherical*sin(m*this->phi);
  else
    return factor*spherical*cos(m*this->phi);
}


double
event::MCweight(int reflectivity, int l1, int m1, int l2, int m2) const
{
  // This is real and no conjugate is employed because of the special
  // form of the two-pseudoscalar decay amplitudes.
  return (this->decayAmplitude(reflectivity,l1,m1)
	  * this->decayAmplitude(reflectivity,l2,m2));
}

complex<double>
coherent_waves::sum(const double* x, const event& e) const
{
  complex<double> result = 0;
  size_t idx = 0;
  vector<wave>::const_iterator it;
  for (it = this->waves.begin(); it != this->waves.end(); it++)
    {
      complex<double> a(x[idx], x[idx+1]);
      result += a * e.decayAmplitude(this->reflectivity,*it);
      idx += 2;
    }

  return result;
}



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
	  double countAllMC = 0;
	  for (size_t iEvent = 0; iEvent < MCallEvents.size(); iEvent++)
	    {
	      if (!MCallEvents[iEvent].accepted())
		continue;
	      countAllMC += 1;
	    }
	  binnedEtaAcc[iBin] = binnedMCevents[iBin].size() / countAllMC;
	}
      else
	{
	  binnedEtaAcc[iBin] = 1;  // No acceptance effects -> All accepted.
	}
    }
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

map<int, double> weights;

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
  cout << "calculated MCweight " << sum / (nmcevents - countRejected)
       << " from " << (nmcevents - countRejected) << " MC events "
       << "out of " << countGenerated << " for "
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



void
myFit()
{
  gRandom = new TRandom1;

  vector<wave> positive;
  positive.push_back(wave("D+", 2, 1));
  positive.push_back(wave("P+", 1, 1));
  positive.push_back(wave("F+", 4, 1));

  vector<wave> negative;
  negative.push_back(wave("S0", 0, 0));
  negative.push_back(wave("P0", 1, 0));
  negative.push_back(wave("P-", 1, 1));
  negative.push_back(wave("D0", 2, 0));
  negative.push_back(wave("D-", 2, 1));

  coherent_waves wsPos, wsNeg;
  wsPos.reflectivity = +1;
  wsPos.spinflip = +1;
  wsPos.waves = positive;

  wsNeg.reflectivity = -1;
  wsNeg.spinflip = +1;
  wsNeg.waves = negative;

  vector<coherent_waves> ws;
  ws.push_back(wsPos);
  ws.push_back(wsNeg);

  assert(NWAVES == negative.size() + positive.size());

  struct {
    const char* name;
    double value;
    bool fixed;
  }
  startingValues[2 * NWAVES] =
    { { "Rea(+,2,1)", gRandom->Uniform(5), false },
      { "Ima(+,2,1)", 0, true },
      { "Rea(+,1,1)", gRandom->Uniform(5), false },
      { "Ima(+,1,1)", gRandom->Uniform(5), false },
      { "Rea(+,4,1)", gRandom->Uniform(1), false },
      { "Ima(+,4,1)", gRandom->Uniform(1), false },
      { "Rea(-,0,0)", gRandom->Uniform(5), false },
      { "Ima(-,0,0)", 0, true },
      { "Rea(-,1,0)", gRandom->Uniform(5), false },
      { "Ima(-,1,0)", gRandom->Uniform(5), false },
      { "Rea(-,1,1)", gRandom->Uniform(5), false },
      { "Ima(-,1,1)", gRandom->Uniform(5), false },
      { "Rea(-,2,0)", gRandom->Uniform(5), false },
      { "Ima(-,2,0)", gRandom->Uniform(5), false },
      { "Rea(-,2,1)", gRandom->Uniform(5), false },
      { "Ima(-,2,1)", gRandom->Uniform(5), false },
    };

  TH2* hRD = new TH2D("hRD", "RD", 10, -1, 1, 10, -M_PI, M_PI);
  TH2* hMC = new TH2D("hMC", "MC", 10, -1, 1, 10, -M_PI, M_PI);
  TH1* hMass = new TH1D("hMass", "mass distribution",
			250, 0.5, 3);
  TH1* htprime = new TH1D("htprime", "t' distribution",
			  250, 0, 1);
  TH1* hMassMC = new TH1D("hMassMC", "MC mass distribution",
			250, 0.5, 3);
  TH1* htprimeMC = new TH1D("htprimeMC", "t' distribution",
			  250, 0, 1);

  //vector<event> RDevents(2500, event(0,0));

  FILE* fd = fopen(dataFile.c_str(), "r");
  if (!fd)
    {
      cerr << "Can't open input file '" << dataFile << "'." << endl;
      abort();
    }
  char line[99999];

  vector<event> RDevents;
  while (fgets(line, 99999, fd) /*&& nrdevents < 5000*/)
    {
      double m, tPr, theta, phi;
      sscanf(line, "%lf %lf %lf %lf", &m, &tPr, &theta, &phi);

      hMass->Fill(m);
      htprime->Fill(tPr);
      event e(m, tPr, theta, phi);
      RDevents.push_back(e);
      hRD->Fill(cos(theta), phi);
    }
  fclose(fd);
  cout << "read " << RDevents.size() << " RD events" << endl;

  vector<event> MCevents;
  vector<event> MCallEvents;
  if (!flatMC)
    {
      // Note that Max writes cos(theta) instead of theta
      fd = fopen(MCFile.c_str(), "r");
      if (!fd)
	{
	  cerr << "Can't open input file '" << dataFile << "'." << endl;
	  abort();
	}
      while (fgets(line, 99999, fd))
	{
	  int acc;
	  double m, tPr, theta, phi;
	  sscanf(line, "%d %lf %lf %lf %lf", &acc, &m, &tPr, &theta, &phi);

	  event e(m, tPr, acos(theta), phi);
	  if (acc)
	    {
	      hMassMC->Fill(m);
	      htprimeMC->Fill(tPr);
	      MCevents.push_back(e);
	    }
	  MCallEvents.push_back(e);
	}
      fclose(fd);
      cout << "read " << MCevents.size() << " accepted MC events out of "
	   << MCallEvents.size() << " total ("
	   << 100.*MCevents.size() / MCallEvents.size()
	   << "% overall acceptance)"
	   << endl;
    }
  else
    {
      for (int iMC = 0; iMC < NFLATMCEVENTS; iMC++)
	{
	  event e(acos(gRandom->Uniform(-1,1)), gRandom->Uniform(-M_PI, M_PI));
	  MCevents.push_back(e);
	}
    }

  likelihood myL(ws, RDevents, MCevents, MCallEvents,
		 nBins, threshold, binWidth);

  TStopwatch sw;

  TFitterMinuit* minuit = new TFitterMinuit();
  minuit->SetMinuitFCN(&myL);

  double lower = threshold;
  double upper = threshold + nBins*binWidth;
  TH1* hSwave = new TH1D("hSwave", "S wave intensity (#epsilon = -1)",
			 nBins, lower, upper);
  TH1* hD0wave = new TH1D("hD0wave", "D0 wave intensity (#epsilon = -1)",
			  nBins, lower, upper);
  TH1* hDwave = new TH1D("hDwave", "D wave intensity",
			 nBins, lower, upper);
  TH1* hPwave = new TH1D("hPwave", "P wave intensity",
			 nBins, lower, upper);
  TH1* hFwave = new TH1D("hFwave", "F wave intensity",
			 nBins, lower, upper);
  TH1* hPhaseDP = new TH1D("hPhaseDP", "D - P phase",
			 nBins, lower, upper);
  TH1* hPhaseDF = new TH1D("hPhaseDF", "D - F phase",
			 nBins, lower, upper);
  TH1* hPhaseD0S = new TH1D("hPhaseD0S", "D0 - S phase",
			    nBins, lower, upper);
  TH1* hIntensity = new TH1D("hIntensity", "total intensity as predicted",
			     nBins, lower, upper);
  //TH3* hPredict = new TH3D("hPredict", "prediction", nBins, 0, nBins, 100, -1, 1, 100, -M_PI, M_PI);
  TStopwatch fulltime;
  fulltime.Start();
  for (iBin = 0; iBin < nBins; iBin++)
    {
      sw.Start();

      myL.setBin(iBin);

      minuit->Clear();
      if (!flatMC)
	// MCweights are the same in different mass bins for the
	// two-body amplitudes.
	weights.clear();

      // The MC part of the likelihood function will evaluate to the
      // number of events in the fit.  In order to speed up the
      // calculation we scale the starting values such that this
      // condition obtains.
      vector<double> vStartingValues(2*NWAVES);
      for (size_t iSV = 0; iSV < 2*NWAVES; iSV++)
	{
	  vStartingValues[iSV] = startingValues[iSV].value;
	}
      size_t nGood = 0;
      for (size_t iEvent = 0; iEvent < RDevents.size(); iEvent++)
	{
	  if (!RDevents[iEvent].accepted())
	    continue;
	  nGood++;
	}
      if (nGood == 0)
	continue;
      double ratio = sqrt(nGood / myL.calc_mc_likelihood(vStartingValues));
      for (size_t iSV = 0; iSV < 2*NWAVES; iSV++)
	{
	  vStartingValues[iSV] *= ratio;
	}

      for (int j= 0; j < 2*NWAVES; j++)
	{
	  minuit->SetParameter(j, startingValues[j].name,
			       startingValues[j].value*ratio, 1, 0, 0);
	  if (startingValues[j].fixed)
	    minuit->FixParameter(j);
	}

      minuit->CreateMinimizer();
      int iret = minuit->Minimize();
      sw.Stop();
      cout << "iret = " << iret << " after " << sw.CpuTime() << " s." << endl;
      if (iret == 0)
	{
	  vector<double> vStartingValue(2*NWAVES);
	  for (int j = 0; j < 2*NWAVES; j++)
	    {
	      vStartingValues[j] = minuit->GetParameter(j);
	    }

	  for (int j = 0; j < minuit->GetNumberTotalParameters(); j++)
	    startingValues[j].value = minuit->GetParameter(j);

	  complex<double> aDwave(minuit->GetParameter(0),
				 minuit->GetParameter(1));
	  complex<double> aPwave(minuit->GetParameter(2),
				 minuit->GetParameter(3));
	  complex<double> aFwave(minuit->GetParameter(4),
				 minuit->GetParameter(5));
	  complex<double> aSwave(minuit->GetParameter(6),
				 minuit->GetParameter(7));
	  complex<double> aD0wave(minuit->GetParameter(12),
				  minuit->GetParameter(13));

	  hDwave->SetBinContent(iBin+1,norm(aDwave));
	  hDwave->SetBinError(iBin+1, 2*abs(aDwave)*minuit->GetParError(0));
	  hPwave->SetBinContent(iBin+1,norm(aPwave));
	  double error = 2*(sqrt(pow(real(aPwave) * minuit->GetParError(2), 2)
			       + pow(imag(aPwave) * minuit->GetParError(3), 2)
			       + (2*real(aPwave)*imag(aPwave)
				  * minuit->GetCovarianceMatrixElement(2, 3))));
	  hPwave->SetBinError(iBin+1, error);

	  double phase = arg(aDwave / aPwave);
	  if (iBin > 1)
	    {
	      double oldPhase = hPhaseDP->GetBinContent(iBin);
	      while (phase - oldPhase > M_PI)
		phase -= 2*M_PI;
	      while (oldPhase - phase > M_PI)
		phase += 2*M_PI;
	    }
	  hPhaseDP->SetBinContent(iBin+1,phase);
	  hPhaseDP->SetBinError(iBin+1, .2);

	  error = 2*(sqrt(pow(real(aFwave) * minuit->GetParError(4), 2)
			       + pow(imag(aFwave) * minuit->GetParError(5), 2)
			       + (2*real(aFwave)*imag(aFwave)
				  * minuit->GetCovarianceMatrixElement(4, 5))));
	  hFwave->SetBinContent(iBin+1,norm(aFwave));
	  hFwave->SetBinError(iBin+1, error);


	  error = 2*(sqrt(pow(real(aSwave) * minuit->GetParError(6), 2)
			       + pow(imag(aSwave) * minuit->GetParError(7), 2)
			       + (2*real(aSwave)*imag(aSwave)
				  * minuit->GetCovarianceMatrixElement(6, 7))));
	  hSwave->SetBinContent(iBin+1,norm(aSwave));
	  hSwave->SetBinError(iBin+1, error);

	  error = 2*(sqrt(pow(real(aD0wave) * minuit->GetParError(12), 2)
			       + pow(imag(aD0wave) * minuit->GetParError(13), 2)
			       + (2*real(aD0wave)*imag(aD0wave)
				  * minuit->GetCovarianceMatrixElement(12, 13))));
	  hD0wave->SetBinContent(iBin+1,norm(aD0wave));
	  hD0wave->SetBinError(iBin+1, error);

	  phase = arg(aD0wave / aSwave);
	  if (iBin > 1)
	    {
	      double oldPhase = hPhaseD0S->GetBinContent(iBin);
	      while (phase - oldPhase > M_PI)
		phase -= 2*M_PI;
	      while (oldPhase - phase > M_PI)
		phase += 2*M_PI;
	    }
	  hPhaseD0S->SetBinContent(iBin+1,phase);
	  hPhaseD0S->SetBinError(iBin+1, .2);

	  phase = arg(aDwave / aFwave);
	  if (iBin > 1)
	    {
	      double oldPhase = hPhaseDF->GetBinContent(iBin);
	      if (phase - oldPhase > M_PI)
		phase -= 2*M_PI;
	      else if (oldPhase - phase > M_PI)
		phase += 2*M_PI;
	    }
	  hPhaseDF->SetBinContent(iBin+1,phase);
	  hPhaseDF->SetBinError(iBin+1, .2);

	  vector<double> result;
	  for (int iPar = 0; iPar < minuit->GetNumberTotalParameters(); iPar++)
	    {
	      result.push_back(minuit->GetParameter(iPar));
	    }

	  double intensity = 0;
	  for (int ix = 0; ix < 100; ix++)
	    {
	      double x = -1 + 2./100*ix;
	      for (int iy = 0; iy < 100; iy++)
		{
		  double y = -M_PI + 2*M_PI/100*iy;
		  intensity += myL.probabilityDensity(result, acos(x), y);
		  //hPredict->SetBinContent(iBin+1, ix + 1, iy + 1, myL.probabilityDensity(result, acos(x), y));
		}
	    }
	  hIntensity->SetBinContent(iBin + 1, intensity / 10000);
	}
    }

  fulltime.Stop();
  cout << "took " << fulltime.CpuTime() << " s CPU time, " << fulltime.RealTime() << " s wall time" << endl;

#if 0
  TCanvas* c = new TCanvas();
  c->Divide(5, 5);
  for (int i = 1; i <= 25; i++)
    {
      hPredict->GetXaxis()->SetRange(i, i);
      char name[555];
      sprintf(name, "yz%d", i);
      TH2* hProjection = (TH2*)hPredict->Project3D(name);
      c->cd(i);
      hProjection->Draw("colz");
    }

  c = new TCanvas();
  c->Divide(2,2);
  c->cd(1); hPwave->Draw();
  c->cd(2); hDwave->Draw();
  c->cd(3); hFwave->Draw();
  c->cd(4); hPhaseDP->Draw();
#endif
}


int main()
{
  if (!readControlFile("control.txt"))
    {
      return 1;
    }

  for (int i = 0; i < nFits; i++)
    {
      char outFileName[999];
      snprintf(outFileName, 999, "out%2.2d.root", i);
      TFile* out = TFile::Open(outFileName, "RECREATE");
      out->cd();
      myFit();
      out->Write();
      delete out;
    }

  return 0;
}
