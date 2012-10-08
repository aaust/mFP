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
#include "control.h"

class plotter 
{
public:
  virtual Long_t getEntries() = 0;
  virtual void fillHistograms(const fitInfo* info, double weight) = 0;
  virtual void getEntry(Long_t n) = 0;
  virtual bool isAccepted() = 0;
  virtual event getEvent() = 0;
};

class tobiPlotter : public plotter
{
private:
  TTree *trMC;
  TFile *fMC;

  // General
  int acc;
  float mX;
  float tPr;
  float costh;
  float phi;

  // specific for event type
  float pvx;
  float pvy;
  float pvz;
  float pg1X;
  float pg1Y;
  float pg1Z;
  float pg1E;
  float pg2X;
  float pg2Y;
  float pg2Z;
  float pg2E;
  float pPim3X;
  float pPim3Y;
  float pPim3Z;
  float pPim3E;
  float pPip3X;
  float pPip3Y;
  float pPip3Z;
  float pPip3E;
  float pEta3X;
  float pEta3Y;
  float pEta3Z;
  float pEta3E;
  float pPimX;
  float pPimY;
  float pPimZ;
  float pPimE;

  TLorentzVector lvEta3;
  TLorentzVector lvPim3;
  TLorentzVector lvPip3;

public:
  tobiPlotter(const char* fnMC)
  {
    TDirectory *oldDir = gDirectory;
    fMC = TFile::Open(fnMC, "READ");
    if (!fMC)
      {
	cerr << "Can't open input file '" << fnMC << "'." << endl;
	abort();
      }

    fMC->GetObject("trMC", trMC);
    if (!trMC)
      {
	cerr << "Can't find tree 'trMC' in file '" << fnMC << "'." << endl;
	abort();
      }
    oldDir->cd();

    TTree *t = trMC;
    t->SetBranchAddress("accepted", &acc);
    t->SetBranchAddress("mX", &mX);
    t->SetBranchAddress("tPrime", &tPr);
    t->SetBranchAddress("costhGJ", &costh);
    t->SetBranchAddress("phiGJ", &phi);

    t->SetBranchAddress("pvx", &pvx);
    t->SetBranchAddress("pvy", &pvy);
    t->SetBranchAddress("pvz", &pvz);
    t->SetBranchAddress("pg1X", &pg1X);
    t->SetBranchAddress("pg1Y", &pg1Y);
    t->SetBranchAddress("pg1Z", &pg1Z);
    t->SetBranchAddress("pg1E", &pg1E);
    t->SetBranchAddress("pg2X", &pg2X);
    t->SetBranchAddress("pg2Y", &pg2Y);
    t->SetBranchAddress("pg2Z", &pg2Z);
    t->SetBranchAddress("pg2E", &pg2E);
    t->SetBranchAddress("pPim3X", &pPim3X);
    t->SetBranchAddress("pPim3Y", &pPim3Y);
    t->SetBranchAddress("pPim3Z", &pPim3Z);
    t->SetBranchAddress("pPim3E", &pPim3E);
    t->SetBranchAddress("pPip3X", &pPip3X);
    t->SetBranchAddress("pPip3Y", &pPip3Y);
    t->SetBranchAddress("pPip3Z", &pPip3Z);
    t->SetBranchAddress("pPip3E", &pPip3E);
    t->SetBranchAddress("pEta3X", &pEta3X);
    t->SetBranchAddress("pEta3Y", &pEta3Y);
    t->SetBranchAddress("pEta3Z", &pEta3Z);
    t->SetBranchAddress("pEta3E", &pEta3E);
    t->SetBranchAddress("pPimX", &pPimX);
    t->SetBranchAddress("pPimY", &pPimY);
    t->SetBranchAddress("pPimZ", &pPimZ);
    t->SetBranchAddress("pPimE", &pPimE);
  }
  
  ~tobiPlotter() {
    delete fMC;
  }
  
  
  Long_t getEntries()
  { return trMC->GetEntries(); }
  void getEntry(Long_t n)
  { trMC->GetEntry(n); }
  
  bool isAccepted() { return acc; }
  
  event getEvent() { return event(mX, tPr, acos(costh), phi); }
  
  void fillHistograms(const fitInfo* info, double weight)
  {
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
};

class alexPlotter : public plotter
{
private:
  TTree *trMC;
  TFile *fMC;

  // General
  int acc;
  float mX;
  float tPr;
  float costh;
  float phi;

  // specific for event type
  float pvx;
  float pvy;
  float pvz;
  float pPX;
  float pPY;
  float pPZ;
  float pPE;
  float pPipX;
  float pPipY;
  float pPipZ;
  float pPipE;
  float pPimX;
  float pPimY;
  float pPimZ;
  float pPimE;
  
public:
  alexPlotter(const char* fnMC)
  {
    TDirectory *oldDir = gDirectory;
    fMC = TFile::Open(fnMC, "READ");
    if (!fMC)
      {
	cerr << "Can't open input file '" << fnMC << "'." << endl;
	abort();
      }
    
    fMC->GetObject("trMC", trMC);
    if (!trMC)
      {
	cerr << "Can't find tree 'trMC' in file '" << fnMC << "'." << endl;
	abort();
      }
    oldDir->cd();
    
    TTree *t = trMC;
    t->SetBranchAddress("accepted", &acc);
    t->SetBranchAddress("mX", &mX);
    t->SetBranchAddress("tPrime", &tPr);
    t->SetBranchAddress("costhGJ", &costh);
    t->SetBranchAddress("phiGJ", &phi);
    
    t->SetBranchAddress("pvx", &pvx);
    t->SetBranchAddress("pvy", &pvy);
    t->SetBranchAddress("pvz", &pvz);
    t->SetBranchAddress("pPX", &pPX);
    t->SetBranchAddress("pPY", &pPY);
    t->SetBranchAddress("pPZ", &pPZ);
    t->SetBranchAddress("pPE", &pPE);
    t->SetBranchAddress("pPipX", &pPipX);
    t->SetBranchAddress("pPipY", &pPipY);
    t->SetBranchAddress("pPipZ", &pPipZ);
    t->SetBranchAddress("pPipE", &pPipE);
    t->SetBranchAddress("pPimX", &pPimX);
    t->SetBranchAddress("pPimY", &pPimY);
    t->SetBranchAddress("pPimZ", &pPimZ);
    t->SetBranchAddress("pPimE", &pPimE);
  }
  
  ~alexPlotter() {
    delete fMC;
  }
  
  Long_t getEntries()
  { return trMC->GetEntries(); }
  void getEntry(Long_t n)
  { trMC->GetEntry(n); }
  
  bool isAccepted() { return acc; }
  
  event getEvent() { return event(mX, tPr, acos(costh), phi); }
  
  void fillHistograms(const fitInfo* info, double weight)
  {
    gHist.Fill("hMvsCosth", "m(X) vs cos(#theta)",
	       info->getNbins(), info->getLower(), info->getUpper(), 100, -1, 1,
	       mX, costh, weight);
    //if ((fabs(costh)) < 0.85)
    //if ((fabs(costh)) < (-0.125*mX + 1.05))
    gHist.Fill("hMvsPhi", "m vs #phi",
	       info->getNbins(), info->getLower(), info->getUpper(), 100, -M_PI, M_PI,
	       mX, phi, weight);
    gHist.Fill("hM", "m",
	       info->getNbins(), info->getLower(), info->getUpper(),
	       mX, weight);
    gHist.Fill("hpvz", "Primary Vertex along Z", 500, -100, 0,
	       pvz, weight);
    gHist.Fill("htPr", "t'", 1000, 0, 1,
	       tPr, weight);
    gHist.Fill("hpPE", "energy of proton",
	       500, 0, 200,
	       pPE, weight);
    gHist.Fill("hpPipE", "energy of #pi+",
	       500, 0, 200,
	       pPipE, weight);
    gHist.Fill("hpPimE", "energy of #pi-",
	       500, 0, 200,
	       pPimE, weight);
  }
};
  
int
main()
{
  /*const char *controlfn = "control.txt";
    if (!readControlFile(controlfn))
    {
    return 1;
    }
  */
  
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

  //const char* fnMC = MCFiles[].c_str();
  //const char *fnMC = "/afs/cern.ch/work/a/aaust/public/200.2000.pf160.15.mc.root";
  const char *fnMC = "/afs/cern.ch/user/a/aaust/scratch0/KKtrees/950.2450.pf150.KpoKm.mc.root";
  plotter *currentPlotter = new alexPlotter(fnMC);
  
  TFile *fOutput = TFile::Open("predict.root", "RECREATE");

  Long_t nTotal = currentPlotter->getEntries();
  for (Long_t i = 0; i < nTotal; i++)
    {
      if ((i+1)*10LL/nTotal > i*10LL/nTotal)
	cout << (i+1)*100LL/nTotal << "% " << flush;

      currentPlotter->getEntry(i);
      if (!currentPlotter->isAccepted())
	continue;
      event e = currentPlotter->getEvent();

      if (e.mass > upperBounds[upperBounds.size() - 1]
	  || e.mass < lowerBounds[0])
	continue;

      // Find bin containing this event.
      size_t idx = lowerBounds.size() / 2;
      size_t limitLow = 0;
      size_t limitHigh = lowerBounds.size();
      while (!(e.mass >= lowerBounds[idx] && e.mass < upperBounds[idx])
	     && limitLow < limitHigh)
	{
	  if (e.mass < lowerBounds[idx])
	    limitHigh = idx;
	  else
	    limitLow = idx+1;
	  idx = (limitHigh + limitLow) / 2;
	  //cout << mX << idx << " " << limitLow << " " << limitHigh
	  // << " " << lowerBounds[limitLow] << " " << upperBounds[limitHigh-1] << endl;
	}
      if (limitLow >= limitHigh)
	{
	  //cout << "no bin found for " << e.mass << endl;
	  continue;
	}

      double weight = info->getWaveSet().getEventWeight(allValues[idx], e) / NMCperBin[idx];
      currentPlotter->fillHistograms(info, weight);
    }
  cout << endl;
  fOutput->Write();
}
