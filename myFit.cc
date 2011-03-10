#include <complex>
#include <vector>
#include <map>

using namespace std;

#include <stdio.h>

#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TRandom1.h"
#include "Math/SpecFunc.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"
#include "TStopwatch.h"

#define NRDEVENTS 100000
size_t nrdevents;
#define NMCEVENTS 10000 //(4*NRDEVENTS)

struct wave {
  int l, m;
  std::complex<double> a;

  wave() { }
  wave(int ll, int mm) { l = ll; m = mm; a = 0; }
  wave(const wave& o) { l = o.l; m = o.m; a = o.a; }
};

class event;

struct coherent_waves {
  int reflectivity;
  int spinflip;
  vector<wave> waves;

  coherent_waves() {}
  coherent_waves(const coherent_waves& o) { reflectivity = o.reflectivity; spinflip = o.spinflip; waves = o.waves; }

  complex<double> sum(const vector<double>& x, size_t idx_base, const event& e) const;
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
};

event RDevents[NRDEVENTS]; // not a reference because of Cint limitations
event MCevents[NMCEVENTS]; // idem

class likelihood : public ROOT::Minuit2::FCNBase {
  vector<coherent_waves> ws;
public:
  likelihood(waveset ws_) /*,
	     vector<event>& RDevents_,
	     vector<event>& MCevents_) */
    : ws(ws_) /*, RDevents(RDevents_), MCevents(MCevents_)*/ {}
  //double operator() (const std::vector<double>& x) const;
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
  calc_likelihood(const vector<double>& x) const;

public:
  double
  operator()(const vector<double>& x) const;

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
coherent_waves::sum(const vector<double>& x, size_t idx_base, const event& e) const
{
  complex<double> result = 0;
  size_t idx = idx_base;
  vector<wave>::const_iterator it;
  for (it = this->waves.begin(); it != this->waves.end(); it++)
    {
      complex<double> a(x[idx], x[idx+1]);
      result += a * e.decayAmplitude(this->reflectivity,*it);
      idx += 2;
    }

  return result;
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
      sum += norm(it->sum(x, idx_base, e));
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
  int countRejected = 0;
  for (size_t i = 0; i < NMCEVENTS /*MCevents.size()*/; i++)
    {
      if (0 && cos(MCevents[i].theta) < -.8)
	{
	  countRejected++;
	  continue;
	}
      // Forming this sum as REAL sum and with no conjugate, because
      // the form of the decay amplitudes allows this.  This is not
      // the most general form!
      double y = MCevents[i].MCweight(reflectivity,w1,w2) - c;
      double t = sum + y;
      c = (t - sum) - y;  // compensation term.
      sum = t;
    }

  return (weights[id] = sum / (NMCEVENTS - countRejected));
}

double massLow = 0;
double massHigh = 9999;

double
likelihood::calc_likelihood(const vector<double>& x) const
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

  // Uses Kahan summation
  double sumRD = 0;
  double c = 0;
  for (size_t i = 0; i < nrdevents; i++)
    {
      if (RDevents[i].mass < massLow || RDevents[i].mass > massHigh)
	continue;
      double y = log(probabilityDensity(x, RDevents[i])) - c;
      double t = sumRD + y;
      c = (t - sumRD) - y;
      sumRD = t;
      //sumRD += y;
    }

  //cout << "likelihood = " << sumRD << " - " << sumMC << endl;
  return sumRD - sumMC;
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
  positive.push_back(wave(2, 1));
  positive.push_back(wave(1, 1));
  positive.push_back(wave(4, 1));

  vector<wave> negative;
  negative.push_back(wave(0,0));
  negative.push_back(wave(1,0));
  negative.push_back(wave(1,1));
  negative.push_back(wave(2,0));
  negative.push_back(wave(2,1));

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

  TH2* hRD = new TH2D("hRD", "RD", 10, -1, 1, 10, -M_PI, M_PI);
  TH2* hMC = new TH2D("hMC", "MC", 10, -1, 1, 10, -M_PI, M_PI);
  TH1* hMass = new TH1D("hMass", "mass distribution",
			250, 0.5, 3);

  //vector<event> RDevents(2500, event(0,0));

  FILE* fd = fopen("dataEtaPpi.txt", "r");
  char line[99999];
  nrdevents = 0;
  while (fgets(line, 99999, fd))
    {
      double m, tPr, theta, phi;
      sscanf(line, "%lf %lf %lf %lf", &m, &tPr, &theta, &phi);

      hMass->Fill(m);
      event e(m, tPr, theta, phi);
      RDevents[nrdevents++] = e;
      hRD->Fill(cos(theta), phi);
    }
  fclose(fd);
  cout << "read " << nrdevents << " RD events" << endl;
  //nrdevents = nrdevents / 10;

  //vector<event> MCevents(10000, event(0,0));
  for (int i = 0; i < NMCEVENTS; i++)
    {
      event e(acos(gRandom->Uniform(-1,1)), gRandom->Uniform(-4*atan(1),4*atan(1)));
      MCevents[i] =  e;
      hMC->Fill(cos(e.theta), e.phi);
    }
  likelihood myL(ws); //, RDevents, MCevents);

  TStopwatch sw;

  TFitterMinuit* minuit = new TFitterMinuit();
  minuit->SetMinuitFCN(&myL);

  struct {
    const char* name;
    double value;
    bool fixed;
  }
  startingValues[16] = { { "Rea(+,2,1)", gRandom->Uniform(5), false },
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
			 { "Ima(-,2,1)", gRandom->Uniform(5), false } };

  int nBins = 40;
  double binWidth = 0.05;

  TH1* hDwave = new TH1D("hDwave", "D wave intensity",
			 nBins, 0.65, 0.65 + nBins*binWidth);
  TH1* hPwave = new TH1D("hPwave", "P wave intensity",
			 nBins, 0.65, 0.65 + nBins*binWidth);
  TH1* hFwave = new TH1D("hFwave", "F wave intensity",
			 nBins, 0.65, 0.65 + nBins*binWidth);
  TH1* hPhaseDP = new TH1D("hPhaseDP", "D - P phase",
			 nBins, 0.65, 0.65 + nBins*binWidth);
  TH1* hPhaseDF = new TH1D("hPhaseDF", "D - F phase",
			 nBins, 0.65, 0.65 + nBins*binWidth);
  TH1* hIntensity = new TH1D("hIntensity", "total intensity as predicted",
			     nBins, 0.65, 0.65 + nBins*binWidth);
  TH3* hPredict = new TH3D("hPredict", "prediction", nBins, 0, nBins, 100, -1, 1, 100, -M_PI, M_PI);
  TStopwatch fulltime;
  fulltime.Start();
  for (int i = 0; i < nBins; i++)
    {
      massLow = 0.65 + i*binWidth;
      massHigh = 0.65 + (i+1)*binWidth;

      cout << "mass bin [" << massLow << ", " << massHigh << "]" << endl;
      sw.Start();

      minuit->Clear();

      for (int j= 0; j < 16; j++)
	{
	  minuit->SetParameter(j, startingValues[j].name,
			       startingValues[j].value, 1, 0, 0);
	  if (startingValues[j].fixed)
	    minuit->FixParameter(j);
	}

      minuit->CreateMinimizer();
      int iret = minuit->Minimize();
      sw.Stop();
      cout << "iret = " << iret << " after " << sw.CpuTime() << " s." << endl;
      if (iret == 0)
	{
	  for (int j = 0; j < minuit->GetNumberTotalParameters(); j++)
	    startingValues[j].value = minuit->GetParameter(j);

	  complex<double> aDwave(startingValues[0].value,
				 startingValues[1].value);
	  complex<double> aPwave(startingValues[2].value,
				 startingValues[3].value);
	  complex<double> aFwave(startingValues[4].value,
				 startingValues[5].value);

	  hDwave->SetBinContent(i+1,norm(aDwave));
	  hDwave->SetBinError(i+1, 2*abs(aDwave)*minuit->GetParError(0));
	  hPwave->SetBinContent(i+1,norm(aPwave));
	  double error = 2*(sqrt(pow(real(aPwave) * minuit->GetParError(2), 2)
			       + pow(imag(aPwave) * minuit->GetParError(3), 2)
			       + (2*real(aPwave)*imag(aPwave)
				  * minuit->GetCovarianceMatrixElement(2, 3))));
	  hPwave->SetBinError(i+1, error);

	  double phase = arg(aDwave / aPwave);
	  if (i > 1)
	    {
	      double oldPhase = hPhaseDP->GetBinContent(i);
	      if (phase - oldPhase > M_PI)
		phase -= 2*M_PI;
	      else if (oldPhase - phase > M_PI)
		phase += 2*M_PI;
	    }
	  hPhaseDP->SetBinContent(i+1,phase);
	  hPhaseDP->SetBinError(i+1, .2);

	  error = 2*(sqrt(pow(real(aFwave) * minuit->GetParError(4), 2)
			       + pow(imag(aFwave) * minuit->GetParError(5), 2)
			       + (2*real(aFwave)*imag(aFwave)
				  * minuit->GetCovarianceMatrixElement(4, 5))));
	  hFwave->SetBinContent(i+1,norm(aFwave));
	  hFwave->SetBinError(i+1, error);

	  phase = arg(aDwave / aFwave);
	  if (i > 1)
	    {
	      double oldPhase = hPhaseDF->GetBinContent(i);
	      if (phase - oldPhase > M_PI)
		phase -= 2*M_PI;
	      else if (oldPhase - phase > M_PI)
		phase += 2*M_PI;
	    }
	  hPhaseDF->SetBinContent(i+1,phase);
	  hPhaseDF->SetBinError(i+1, .2);

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
		  hPredict->SetBinContent(i, ix + 1, iy + 1, myL.probabilityDensity(result, acos(x), y));
		}
	    }
	  hIntensity->SetBinContent(i + 1, intensity / 10000);
	}
    }

  fulltime.Stop();
  cout << "took " << fulltime.CpuTime() << " s CPU time, " << fulltime.RealTime() << " s wall time" << endl;


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
      
}
