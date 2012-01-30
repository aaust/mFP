#include <complex>
#include <iostream>
using namespace std;

#include <math.h>

#include "Math/SpecFunc.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "gHist.h"
#include "fitInfo.h"


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
  UInt_t NMCevents;
  treeResult->SetBranchAddress("massLow", &massLow);
  treeResult->SetBranchAddress("massHigh", &massHigh);
  treeResult->SetBranchAddress("values", values);
  treeResult->SetBranchAddress("NMCevents", &NMCevents);


  vector<double> lowerBounds;
  vector<double> upperBounds;
  vector<vector<double> > allValues;
  vector<UInt_t> NMCperBin;
  for (Long_t i = 0; i < treeResult->GetEntries(); i++)
    {
      // This assumes the bins are sorted.
      treeResult->GetEntry(i);
      lowerBounds.push_back(massLow);
      upperBounds.push_back(massHigh);
      allValues.push_back(vector<double>(values,values+info->getNvars()));
      NMCperBin.push_back(NMCevents);
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
  float pPim3X; t->SetBranchAddress("pPim3X", &pPim3X);
  float pPim3Y; t->SetBranchAddress("pPim3Y", &pPim3Y);
  float pPim3Z; t->SetBranchAddress("pPim3Z", &pPim3Z);
  float pPim3E; t->SetBranchAddress("pPim3E", &pPim3E);
  float pPip3X; t->SetBranchAddress("pPip3X", &pPip3X);
  float pPip3Y; t->SetBranchAddress("pPip3Y", &pPip3Y);
  float pPip3Z; t->SetBranchAddress("pPip3Z", &pPip3Z);
  float pPip3E; t->SetBranchAddress("pPip3E", &pPip3E);
  float pEta3X; t->SetBranchAddress("pEta3X", &pEta3X);
  float pEta3Y; t->SetBranchAddress("pEta3Y", &pEta3Y);
  float pEta3Z; t->SetBranchAddress("pEta3Z", &pEta3Z);
  float pEta3E; t->SetBranchAddress("pEta3E", &pEta3E);
  float pPimX; t->SetBranchAddress("pPimX", &pPimX);
  float pPimY; t->SetBranchAddress("pPimY", &pPimY);
  float pPimZ; t->SetBranchAddress("pPimZ", &pPimZ);
  float pPimE; t->SetBranchAddress("pPimE", &pPimE);
  
  TFile *fOutput = TFile::Open("predict.root", "RECREATE");

  Long_t nTotal = t->GetEntries();
  for (Long_t i = 0; i < nTotal; i++)
    {
      if ((i+1)*10LL/nTotal > i*10/nTotal)
	cout << (i+1)*100/nTotal << "% " << flush;

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

      TLorentzVector lvEta3(pEta3X, pEta3Y, pEta3Z, pEta3E);
      TLorentzVector lvPim3(pPim3X, pPim3Y, pPim3Z, pPim3E);
      TLorentzVector lvPip3(pPip3X, pPip3Y, pPip3Z, pPip3E);

      double weight = info->getWaveSet()[0].getEventWeight(allValues[idx], e) / NMCperBin[idx];
      gHist.Fill("hMvsCosth", "m vs cos th",
		 info->getNbins(), info->getLower(), info->getUpper(), 100, -1, 1,
		 mX, costh, weight);
      //if (1.9 < mX && mX < 2.1)
	{
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

	  gHist.Fill("hDalitz", "Dalitz Plot;pi+eta;pi-eta", 50, 0.45, 0.7, 50, 0.45, 0.7,
		     (lvEta3 + lvPip3).M2(), (lvEta3 + lvPim3).M2());
	  gHist.Fill("hDalitzFlipped", "Dalitz Plot;pi+eta;pi-pi+", 50, 0.45, 0.75, 50, 0.06, 0.18,
		     (lvEta3 + lvPip3).M2(), (lvPip3 + lvPim3).M2());
		     
	  //}
	  //      if (mX < 2)
	  //{
	  gHist.Fill("hpPimE", "energy of non-eta' #pi-",
		     500, 0, 200,
		     pPimE, weight);
	}
    }
  cout << endl;
  fOutput->Write();
}
