#include <complex>
#include <iostream>
using namespace std;

#include <math.h>

#include "Math/SpecFunc.h"

#include "TFile.h"
#include "TTree.h"
#include "TRandom1.h"

#include "gHist.h"
#include "fitInfo.h"


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


int
main()
{
  const char *fnResult = "out00.root";
  TFile *fResult = TFile::Open(fnResult, "READ");
  if (!fResult)
    {
      cerr << "Can't open input file '" << fnResult << "'" << endl;
      return 1;
    }
  TTree *treeResult;
  fResult->GetObject("tFitResults", treeResult);
  if (!treeResult)
    {
      cerr << "Can't find input tree 'tFitResults'" << endl;
      return 1;
    }

  fitInfo *info = 0;
  TIter next(treeResult->GetUserInfo());
  TObject *obj;
  while ((obj = next()))
    {
      if (obj->InheritsFrom("fitInfo"))
	{
	  info = (fitInfo*)obj;
	  break;
	}
    }
  if (!info)
    {
      cerr << "Tree without user info" << endl;
      return 1;
    }

  double *values = new double[info->getNvars()];
  double massLow, massHigh;
  treeResult->SetBranchAddress("massLow", &massLow);
  treeResult->SetBranchAddress("massHigh", &massHigh);
  treeResult->SetBranchAddress("values", values);

  vector<double> lowerBounds;
  vector<double> upperBounds;
  vector<vector<double> > allValues;
  for (Long_t i = 0; i < treeResult->GetEntries(); i++)
    {
      // This assumes the bins are sorted.
      treeResult->GetEntry(i);
      lowerBounds.push_back(massLow);
      upperBounds.push_back(massHigh);
      allValues.push_back(vector<double>(values,values+info->getNvars()));
    }

  const char *fnMC = "/data/zup/diefenbach/trees/EtaPr3Pi.root";

  TDirectory *oldDir = gDirectory;
  TFile *fMC = TFile::Open(fnMC, "READ");
  if (!fMC)
    {
      cerr << "Can't open input file '" << fnMC << "'." << endl;
      abort();
    }

  TTree *t;
  fMC->GetObject("trMC", t);
  if (!t)
    {
      cerr << "Can't find tree 'trMC' in file '" << fnMC << "'." << endl;
      abort();
    }
  oldDir->cd();

  int acc;
  float mX;
  float tPr;
  float costh;
  float phi;

  t->SetBranchAddress("accepted", &acc);
  if (t->GetBranch("mX"))
    t->SetBranchAddress("mX", &mX);
  else if (t->GetBranch("mKK"))
    t->SetBranchAddress("mKK", &mX);
  else if (t->GetBranch("mKpi"))
    t->SetBranchAddress("mKpi", &mX);
  else
    {
      cerr << "unknown format of MC input tree" << endl;
      abort();
    }
  t->SetBranchAddress("tPrime", &tPr);
  t->SetBranchAddress("costhGJ", &costh);
  t->SetBranchAddress("phiGJ", &phi);

  float pvx; t->SetBranchAddress("pvx", &pvx);
  float pvy; t->SetBranchAddress("pvy", &pvy);
  float pvz; t->SetBranchAddress("pvz", &pvz);
  float pg1X; t->SetBranchAddress("pg1X", &pg1X);
  float pg1Y; t->SetBranchAddress("pg1Y", &pg1Y);
  float pg1Z; t->SetBranchAddress("pg1Z", &pg1Z);
  float pg1E; t->SetBranchAddress("pg1E", &pg1E);
  float pg2X; t->SetBranchAddress("pg2X", &pg2X);
  float pg2Y; t->SetBranchAddress("pg2Y", &pg2Y);
  float pg2Z; t->SetBranchAddress("pg2Z", &pg2Z);
  float pg2E; t->SetBranchAddress("pg2E", &pg2E);
  float pPimX; t->SetBranchAddress("pPimX", &pPimX);
  float pPimY; t->SetBranchAddress("pPimY", &pPimY);
  float pPimZ; t->SetBranchAddress("pPimZ", &pPimZ);
  float pPimE; t->SetBranchAddress("pPimE", &pPimE);
  
  TFile *fOutput = TFile::Open("predict.root", "RECREATE");

  for (Long_t i = 0; i < t->GetEntries(); i++)
    {
      if ((i+1)*10LL/t->GetEntries() > i*10/t->GetEntries())
	cout << (i+1)*100/t->GetEntries() << "% " << flush;

      t->GetEntry(i);
      if (!acc)
	continue;
      event e(mX, tPr, acos(costh), phi);

      if (mX > upperBounds[upperBounds.size() - 1]
	  || mX < lowerBounds[0])
	continue;

      // Find bin containing this event.
      size_t idx = lowerBounds.size() / 2;
      size_t limitLow = 0;
      size_t limitHigh = lowerBounds.size();
      while (!(mX >= lowerBounds[idx] && mX < upperBounds[idx])
	     && limitLow < limitHigh)
	{
	  if (mX < lowerBounds[idx])
	    limitHigh = idx;
	  else
	    limitLow = idx+1;
	  idx = (limitHigh + limitLow) / 2;
	  //cout << mX << idx << " " << limitLow << " " << limitHigh
	  // << " " << lowerBounds[limitLow] << " " << upperBounds[limitHigh-1] << endl;
	}
      if (limitLow >= limitHigh)
	{
	  cout << "no bin found for " << mX << endl;
	  continue;
	}

      double weight = info->getWaveSet().getEventWeight(allValues[idx], e);
      gHist.Fill("hMvsCosth", "m vs cos th",
		 info->getNbins(), info->getLower(), info->getUpper(), 100, -1, 1,
		 mX, costh, weight);
      gHist.Fill("hpvz", "pvz", 500, -100, 0,
		 pvz, weight);
      gHist.Fill("htPr", "t'", 1000, 0, 1,
		 tPr, weight);
      gHist.Fill("hXYECAL", "XY ECAL", 200, -.04, .04, 200, -0.04, 0.04,
		 pg1X / pg1Z, pg1Y / pg1Z, weight);
      gHist.Fill("hXYECAL", "XY ECAL", 200, -.04, .04, 200, -0.04, 0.04,
		 pg2X / pg2Z, pg2Y / pg2Z, weight);
      gHist.Fill("hEvsAngle", "E vs angle", 200, 0, 200, 100, 0, 0.06,
		 pg1E, hypot(pg1X, pg1Y) / pg1Z, weight);
      gHist.Fill("hEvsAngle", "E vs angle", 200, 0, 200, 100, 0, 0.06,
		 pg2E, hypot(pg2X, pg2Y) / pg2Z, weight);
      if (mX < 2)
	{
	  gHist.Fill("hpPimE", "energy of non-eta' #pi-",
		     500, 0, 200,
		     pPimE, weight);
	}
    }
  cout << endl;
  fOutput->Write();
}
