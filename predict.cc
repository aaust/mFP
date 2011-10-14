#include <complex>
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
  double factor = 1; //sqrt(nrdevents);
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

double
getEventWeight(double mass, double t, double costh, double phi)
{
  complex<double> BWa2 = BW_a2_pietap_coupled(mass*mass);
  complex<double> BWpi1 = BW(mass*mass, mEtaP, mPi, 1, 1.6, 0.4);
  complex<double> BWa4 = BW_a4_pietap(mass*mass);
  double relativeStrength = .5;

  // 1/mass accounts for the phase space, still needs verfication.
  return norm(relativeStrength * BWa2 * decayAmplitude(1, 2, 1, costh, phi)
	      + BWpi1 * decayAmplitude(1, 1, 1, costh, phi)
	      + BWa4 * decayAmplitude(1, 4, 1, costh, phi)) / mass;
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
      phi = gRandom->Uniform(2*M_PI);
    } while(getEventWeight(mass, t, costh, phi) < gRandom->Uniform(maxWeight));

    gHist.Fill("hMvsCosth", "m vs cos th;m;cos th", 100, mEtaP + mPi, mEtaP + mPi + 2, 100, -1, 1,
	       mass, costh);
    tree->Fill();
  } while (++goodEvents < 50000);

  f->Write();
  delete f;
  return 0;
}
