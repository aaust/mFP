#include <complex>
#include <iostream>
using namespace std;

#include <math.h>

#include "Math/SpecFunc.h"

#include "TFile.h"
#include "TTree.h"
#include "TRandom1.h"

#include "gHist.h"
#include "bw.h"

double
getAcceptanceWeight(double mass, double t, double costh, double phi)
{
  return 1;
}


double
decayAmplitude(int reflectivity, int l, int m, double costh, double phi)
{
  double spherical;
  spherical = ROOT::Math::sph_legendre(l, m, acos(costh));
  double factor = .5; //sqrt(nrdevents);
  if (m != 0)
    factor *= sqrt(2.);

  double result;
  if (reflectivity == +1)
    result = factor*spherical*sin(m*phi);
  else
    result = factor*spherical*cos(m*phi);

  //lookupTable[id] = result;
  return result; // * this->tPrime*exp(-7.*this->tPrime);
}


complex<double>
getBWAmplitude(int J, double mass, double t)
{
  switch(J)
    {
    case 1: return BW(mass*mass, mEtaP, mPi, 1, 1.6, 0.4);
    case 2: return 0.5*BW_a2_pietap_coupled(mass*mass);
    case 4: return BW_a4_pietap(mass*mass);
    default: return 0;
    }
}

double
getEventWeight(double mass, double t, double costh, double phi)
{
  // 2*mass because the |BW|^2 are distributions in m^2.
  return norm(getBWAmplitude(1, mass, t) * decayAmplitude(1, 1, 1, costh, phi)
  	      + getBWAmplitude(2, mass, t) * decayAmplitude(1, 2, 1, costh, phi)
	      + getBWAmplitude(4, mass, t) * decayAmplitude(1, 4, 1, costh, phi));
}

int
main()
{
  TFile *f = new TFile("predict.root", "RECREATE");
  f->cd();

  gRandom = new TRandom1;

  // find maximum event weight, simply generate 1e6 events, and record
  // the maximum.
  double maxWeight = 0;
  for (int i = 0; i < 1000000; i++)
    {
      double mass = gRandom->Uniform(2) + mEtaP + mPi;
      double t = 1.;
      double costh = gRandom->Uniform(-1,1);
      double phi = gRandom->Uniform(2*M_PI);

      double w = getEventWeight(mass, t, costh, phi);
      gHist.Fill("hWeights", "weight of randomly generated events", 1000, 0, 1000, w);
      if (w > maxWeight)
	maxWeight = w;
    }
  cout << "maximum weight = " << maxWeight << endl;

  float mass;
  float t = 1;
  float costh;
  float phi;

  TTree *tree = new TTree("tPredict", "prediction");
  tree->Branch("m", &mass, "m/F");
  tree->Branch("t", &t, "t/F");
  tree->Branch("costh", &costh, "costh/F");
  tree->Branch("phi", &phi, "phi/F");

  int goodEvents = 0;
  do {
    do {
      mass = gRandom->Uniform(2) + mEtaP + mPi;
      t = 1.;
      costh = gRandom->Uniform(-1,1);
      phi = gRandom->Uniform(-M_PI, M_PI);
    } while(getEventWeight(mass, t, costh, phi) < gRandom->Uniform(maxWeight));

    gHist.Fill("hMvsCosth", "m vs cos th;m;cos th", 100, mEtaP + mPi, mEtaP + mPi + 2, 100, -1, 1,
	       mass, costh);
    tree->Fill();
  } while (++goodEvents < 50000);

  // Some illustrative plots ...
  const int nBins = 200;
  for (int i = 0; i < nBins; i++)
    {
      double mass = mEtaP + mPi + 2.*i/nBins;
      complex<double> Dwave = getBWAmplitude(2, mass, 1);
      complex<double> phaseD = Dwave / abs(Dwave);
      complex<double> Pwave = getBWAmplitude(1, mass, 1);
      complex<double> Gwave = getBWAmplitude(4, mass, 1);

      gHist.getHist("hReDwave", "Re(Dwave)", nBins, mEtaP + mPi, mEtaP + mPi + 2)->SetBinContent(i, real(Dwave / phaseD));
      gHist.getHist("hRePwave", "Re(Pwave)", nBins, mEtaP + mPi, mEtaP + mPi + 2)->SetBinContent(i, real(Pwave / phaseD));
      gHist.getHist("hImPwave", "im(Pwave)", nBins, mEtaP + mPi, mEtaP + mPi + 2)->SetBinContent(i, imag(Pwave / phaseD));
      gHist.getHist("hReGwave", "Re(Gwave)", nBins, mEtaP + mPi, mEtaP + mPi + 2)->SetBinContent(i, real(Gwave / phaseD));
      gHist.getHist("hImGwave", "Im(Gwave)", nBins, mEtaP + mPi, mEtaP + mPi + 2)->SetBinContent(i, imag(Gwave / phaseD));

      gHist.getHist("hIntD", "intensity D", nBins, mEtaP + mPi, mEtaP + mPi + 2)->SetBinContent(i, norm(Dwave));
      gHist.getHist("hIntP", "intensity P", nBins, mEtaP + mPi, mEtaP + mPi + 2)->SetBinContent(i, norm(Pwave));
      gHist.getHist("hIntG", "intensity G", nBins, mEtaP + mPi, mEtaP + mPi + 2)->SetBinContent(i, norm(Gwave));

      gHist.getHist("hPhaseDP", "arg(D / P)", nBins, mEtaP + mPi, mEtaP + mPi + 2)->SetBinContent(i, arg(Dwave / Pwave));
      gHist.getHist("hPhaseDG", "arg(D / G)", nBins, mEtaP + mPi, mEtaP + mPi + 2)->SetBinContent(i, arg(Dwave / Gwave));
      gHist.getHist("hPhasePG", "arg(P / G)", nBins, mEtaP + mPi, mEtaP + mPi + 2)->SetBinContent(i, arg(Pwave / Gwave));
    }


  f->Write();
  delete f;
  return 0;
}
