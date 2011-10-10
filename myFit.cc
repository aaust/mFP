#include <complex>
#include <vector>
#include <map>

using namespace std;

#include <string>
#include <stdio.h>
#include <omp.h>

#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TRandom1.h"
#include "TFitterMinuit.h"
#include "TStopwatch.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"

#include "control.h"
#include "wave.h"
#include "3j.h"
#include "event.h"
#include "likelihood.h"
#include "gHist.h"

#define NFLATMCEVENTS 100000
double massLow = 0;
double massHigh = 9999;
int iBin;

class combinedLikelihood : public ROOT::Minuit2::FCNBase {
public:
  vector<likelihood*> myLs;
  waveset ws;
  size_t nBins;
  double threshold;
  double binWidth;
public:
  combinedLikelihood(waveset ws_,
		     size_t nBins_, double threshold_, double binWidth_)
    : ws(ws_), nBins(nBins_), threshold(threshold_), binWidth(binWidth_)
  {
  }

  ~combinedLikelihood() { for (size_t i = 0; i < myLs.size(); i++) delete myLs[i]; }

  void
  addChannel(vector<event>& RDevents,
	     vector<event>& MCevents,
	     vector<event>& MCallEvents)
  {
    size_t idxBranching = /*2*NWAVES*/ 16 + this->getNChannels();
    myLs.push_back(new likelihood(ws, RDevents, MCevents, MCallEvents, nBins, threshold, binWidth, idxBranching));
  }


  double Up() const { return 0.5; }

  void setBin(size_t iBin) { for (size_t i = 0; i < myLs.size(); i++) myLs[i]->setBin(iBin); }
  size_t eventsInBin() const {
    size_t sum = 0;
    for (size_t i = 0; i < myLs.size(); i++) sum += myLs[i]->eventsInBin(); 
    return sum;
  }
  void clearWeights() { for (size_t i = 0; i < myLs.size(); i++) myLs[i]->clearWeights(); }

  double
  calc_mc_likelihood(const vector<double>& x) const
  {
    double result = 0;
    for (size_t i = 0; i < myLs.size(); i++)
      result += myLs[i]->calc_mc_likelihood(x);
    return result;
  }

  double
  operator()(const vector<double>& x) const
  {
    double result = 0;
    for (size_t i = 0; i < myLs.size(); i++)
      {
	double BR = 1;
	// We only take the Branching Ratios into account when there's
	// more than one channel.  It should work with a single
	// channel as long as the corresponding variable is present,
	// but who says that it always will?
	if (getNChannels() > 1)
	  {
	    if (i == 0)
	      {
		// First channel gets 1 - \sum BR
		for (size_t j = 0; j < getNChannels() - 1; j++)
		  BR -= x[x.size() - 1 - j];
	      }
	    else
	      {
		BR = x[x.size() - getNChannels() + i];
	      }
	  }
	// The BR enters the individual real data events by
	// multiplying the respective probability density whose
	// logarithm is then summed to yield the log-likelihood, and
	// it can thus be taken out of the logarithm contained in the
	// RD likelihood calculation as follows.
	result += (myLs[i]->eventsInBin()*log(BR) + myLs[i]->calc_rd_likelihood(x)
		   - BR*myLs[i]->calc_mc_likelihood(x));
      }
    return -result;
  }

  size_t getNChannels() const { return myLs.size(); }
};

struct tStartingValue {
  string name;
  Double_t value;
  bool fixed;

  tStartingValue(const string& n, Double_t v, bool f)
    : name(n), value(v), fixed(f)
  {}
};


void
fillRDhists(const event& e)
{
  double lower = threshold;
  double upper = threshold + nBins*binWidth;

  gHist.Fill("htprime", "t' distribution", 250, 0, 1,
	     e.tPrime);
  gHist.Fill("hRD", "RD", 10, -1, 1, 10, -M_PI, M_PI,
	     cos(e.theta), e.phi);
  if (e.mass < 1.5)
    gHist.Fill("hCosThVsPhiLow", "cos(#theta) vs. #phi for low mass;cos(#theta);#phi",
	       20, -1, 1, 20, -M_PI, M_PI,
	       cos(e.theta), e.phi);
  if (e.mass > 2.2)
    gHist.Fill("hCosThVsPhiHigh", "cos(#theta) vs. #phi for low mass;cos(#theta);#phi",
	       20, -1, 1, 20, -M_PI, M_PI,
	       cos(e.theta), e.phi);
  gHist.Fill("hMassFine", "mass distribution",
	     250, threshold, 3, e.mass);
  gHist.Fill("hMassRough", "expected MC likelihood of result",
	     nBins, lower, upper,
	     e.mass);

  gHist.Fill("hMvsCosth", "m vs costh Rd", nBins, lower, upper, 40, -1, 1,
	     e.mass, cos(e.theta));
}

void
fillMChists(const event& e, bool acc)
{
  double lower = threshold;
  double upper = threshold + nBins*binWidth;

  gHist.Fill("hThVsMgen", "generated cos(#theta_{#eta'}) vs M;cos(#theta);M/GeV",
	     100, -1, 1, nBins, lower, upper,
	     cos(e.theta), e.mass);
  gHist.Fill("hPhiVsMgen", "generated #phi vs M;#phi;M/GeV",
	     40, -M_PI, M_PI, nBins, lower, upper,
	     e.phi, e.mass);
  gHist.Fill("hMVsTgen", "generated M vs t';M/GeV;t'/GeV^{2}",
	     nBins, lower, upper, 40, 0.05, 0.45,
	     e.mass, e.tPrime);

  gHist.Fill("hMvsCosthGen", "m vs costh MC Gen", nBins, lower, upper, 40, -1, 1,
	     e.mass, cos(e.theta));

  if (acc)
    {
      gHist.Fill("hThVsMacc", "accepted cos(#theta_{#eta'}) vs M;cos(#theta);M/GeV",
		 100, -1, 1, nBins, lower, upper,
		 cos(e.theta), e.mass);
      gHist.Fill("hPhiVsMacc", "accepted #phi vs M;#phi;M/GeV",
		 40, -M_PI, M_PI, nBins, lower, upper,
		 e.phi, e.mass);
      gHist.Fill("hMVsTacc", "accepted M vs t';M/GeV;t'/GeV^{2}",
		 nBins, lower, upper, 40, 0.05, 0.45,
		 e.mass, e.tPrime);
      gHist.Fill("hMassMC", "MC mass distribution",
		 250, threshold, 3,
		 e.mass);
      gHist.Fill("htprimeMC", "t' distribution",
		 250, 0, 1,
		 e.tPrime);
      gHist.Fill("hMvsCosthAcc", "m vs costh MC Acc", nBins, lower, upper, 40, -1, 1,
		 e.mass, cos(e.theta));

    }
}


void __attribute((noinline))
myFit()
{
  gRandom = new TRandom1;

  double lower = threshold;
  double upper = threshold + nBins*binWidth;

  vector<wave> positive;
  positive.push_back(wave("D+", 2, 1, nBins, lower, upper, true));
  positive.push_back(wave("P+", 1, 1, nBins, lower, upper));
  //positive.push_back(wave("F+", 3, 1, nBins, lower, upper));
  positive.push_back(wave("G+", 4, 1, nBins, lower, upper));
  //positive.push_back(wave("D++", 2, 2, nBins, lower, upper));

  vector<wave> negative;
  negative.push_back(wave("S0", 0, 0, nBins, lower, upper, true));
  negative.push_back(wave("P0", 1, 0, nBins, lower, upper));
  negative.push_back(wave("P-", 1, 1, nBins, lower, upper));
  negative.push_back(wave("D0", 2, 0, nBins, lower, upper));
  negative.push_back(wave("D-", 2, 1, nBins, lower, upper));

  coherent_waves wsPos, wsNeg;
  wsPos.reflectivity = +1;
  wsPos.spinflip = +1;
  wsPos.waves = positive;

  wsNeg.reflectivity = -1;
  wsNeg.spinflip = +1;
  wsNeg.waves = negative;

  waveset ws;
  ws.push_back(wsPos);
  ws.push_back(wsNeg);

  size_t lastIdx = 0;
  for (size_t i = 0; i < ws.size(); i++)
    for (size_t j = 0; j < ws[i].waves.size(); j++, lastIdx += 2)
      ws[i].waves[j].setIndex(lastIdx);

  // Find the non-zero moments for the given waveset.
  std::vector<std::pair<size_t, size_t> > vecMom = listOfMoments(ws);

  // Prepare the moment histograms
  map<std::pair<size_t,size_t>, TH1*> mhMoments;
  for (std::vector<std::pair<size_t, size_t> >::const_iterator it = vecMom.begin();
       it != vecMom.end(); it++)
    {
      char name[999];
      char title[999];
      snprintf(name, 999, "hMoment%zd%zd", it->first, it->second);
      snprintf(title, 999, "Moment H(%zd,%zd)", it->first, it->second);
      mhMoments[*it] = gHist.getHist(name, title, nBins, lower, upper);
    }

  vector<tStartingValue> startingValues;
  for (size_t i = 0; i < positive.size(); i++)
    {
      const wave& w = positive[i];
      startingValues.push_back(tStartingValue((string("Re(") + w.name + string(")")),
					      gRandom->Uniform(5),
					      false));
      if (w.phaseLocked)
	startingValues.push_back(tStartingValue((string("Im(") + w.name + string(")")),
						0,
						true));
      else
	startingValues.push_back(tStartingValue((string("Im(") + w.name + string(")")),
						gRandom->Uniform(5),
						false));
    }
  for (size_t i = 0; i < negative.size(); i++)
    {
      const wave& w = negative[i];
      startingValues.push_back(tStartingValue(tStartingValue((string("Re(") + w.name + string(")")),
							     gRandom->Uniform(5),
							     false)));
      if (w.phaseLocked)
	startingValues.push_back(tStartingValue((string("Im(") + w.name + string(")")),
						0,
						true));
      else
	startingValues.push_back(tStartingValue((string("Im(") + w.name + string(")")),
						gRandom->Uniform(5),
						false));
    }

  startingValues.push_back(tStartingValue("BR1", 1, true));
  startingValues.push_back(tStartingValue("BR2", 0.57, false));

  combinedLikelihood myL(ws, nBins, threshold, binWidth);

  for (size_t iFile = 0; iFile < dataFiles.size(); iFile++)
    {
      vector<event> RDevents;

      const char* fn = dataFiles[iFile].c_str();
      size_t len = strlen(fn);

      if (len > 5
	  && fn[len - 5] == '.' && fn[len - 4] == 'r'
	  && fn[len - 3] == 'o' && fn[len - 2] == 'o'
	  && fn[len - 1] == 't')
	{
	  TDirectory *oldDir = gDirectory;
	  TFile *f = TFile::Open(fn, "READ");
	  if (!f)
	    {
	      cerr << "Can't open input file '" << fn << "'." << endl;
	      abort();
	    }

	  TTree *tree;
	  f->GetObject("events", tree);
	  oldDir->cd();
	  if (tree)
	    {
	      float mX;
	      float tPr;
	      float theta;
	      float phi;
	      float likeK, likePi;

	      tree->SetBranchAddress("mKpi", &mX);
	      tree->SetBranchAddress("tPrime", &tPr);
	      tree->SetBranchAddress("theta", &theta);
	      tree->SetBranchAddress("phi", &phi);
	      tree->SetBranchAddress("likeK", &likeK);
	      tree->SetBranchAddress("likePi", &likePi);

	      RDevents.reserve(tree->GetEntries());
	      for (Long_t i = 0; i < tree->GetEntries(); i++)
		{
		  tree->GetEntry(i);
		  if (likeK != -1 && likePi > 2*likeK)
		    continue;
		  event e(mX, tPr, theta, phi);
		  RDevents.push_back(e);
		  fillRDhists(e);
		}
	    }
	  else
	    {
	      f->GetObject("schluter/trRDEtap3pi", tree);
	      if (!tree)
		{
		  cerr << "can't find physics tree" << endl;
		  abort();
		}

	      float m;
	      float costh;
	      float phi;
	      float t;
	      float mCandEtaP1, mCandEtaP2;

	      tree->SetBranchAddress("m", &m);
	      tree->SetBranchAddress("t", &t);
	      tree->SetBranchAddress("costh", &costh);
	      tree->SetBranchAddress("phi", &phi);
	      tree->SetBranchAddress("mCandEtaP1", &mCandEtaP1);
	      tree->SetBranchAddress("mCandEtaP2", &mCandEtaP2);

	      RDevents.reserve(tree->GetEntries());
	      for (Long_t i = 0; i < tree->GetEntries(); i++)
		{
		  tree->GetEntry(i);
		  gHist.Fill("hmEtap", "m(#pi#pi#eta)", 1000, 0.7, 1.7, mCandEtaP1);
		  gHist.Fill("hmEtap", "m(#pi#pi#eta)", 1000, 0.7, 1.7, mCandEtaP2);
		  if (fabs(mCandEtaP1 - .958) > 0.02 && fabs(mCandEtaP2 - .958) > 0.02)
		    continue;
		  event e(m, -t, acos(costh), phi);
		  RDevents.push_back(e);
		  fillRDhists(e);
		}
	      
	    }
	  f->Close();
	}
      else
	{
	  FILE* fd = fopen(dataFiles[iFile].c_str(), "r");
	  if (!fd)
	    {
	      cerr << "Can't open input file '" << dataFiles[iFile] << "'."
		   << endl;
	      abort();
	    }
	  char line[99999];

	  while (fgets(line, 99999, fd) /*&& nrdevents < 5000*/)
	    {
	      double m, tPr, theta, phi;
	      sscanf(line, "%lf %lf %lf %lf", &m, &tPr, &theta, &phi);

	      event e(m, tPr, theta, phi);
	      RDevents.push_back(e);
	      fillRDhists(e);
	    }
	  fclose(fd);
	}
      cout << "read " << RDevents.size() << " RD events" << endl;

      vector<event> MCevents;
      vector<event> MCallEvents;
      if (!flatMC)
	{
	  const char* fn = MCFiles[iFile].c_str();
	  size_t len = strlen(fn);

	  if (len > 5
	      && fn[len - 5] == '.' && fn[len - 4] == 'r'
	      && fn[len - 3] == 'o' && fn[len - 2] == 'o'
	      && fn[len - 1] == 't')
	    {
	      TDirectory *oldDir = gDirectory;
	      TFile *f = TFile::Open(fn, "READ");
	      if (!f)
		{
		  cerr << "Can't open input file '" << fn << "'." << endl;
		  abort();
		}

	      TTree *t;
	      f->GetObject("trMC", t);
	      if (!t)
		{
		  cerr << "Can't find tree 'trMC' in file '" << fn << "'." << endl;
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

	      MCallEvents.reserve(t->GetEntries());
	      MCevents.reserve(t->GetEntries() / 10);  // Right order of magnitude.
	      for (Long_t i = 0; i < t->GetEntries(); i++)
		{
		  t->GetEntry(i);
		  event e(mX, tPr, acos(costh), phi);

		  fillMChists(e, acc);
		  if (acc)
		    MCevents.push_back(e);
		  MCallEvents.push_back(e);
		}
	      f->Close();
	    }
	  else
	    {
	      // Note that Max writes cos(theta) instead of theta
	      FILE* fd = fopen(MCFiles[iFile].c_str(), "r");
	      if (!fd)
		{
		  cerr << "Can't open input file '" << MCFiles[iFile] << "'." << endl;
		  abort();
		}
	      char line[99999];

	      while (fgets(line, 99999, fd))
		{
		  int acc;
		  double m, tPr, theta, phi;
		  sscanf(line, "%d %lf %lf %lf %lf", &acc, &m, &tPr, &theta, &phi);

		  event e(m, tPr, acos(theta), phi);

		  fillMChists(e, acc);
		  if (acc)
		    MCevents.push_back(e);
		  MCallEvents.push_back(e);
		}
	      fclose(fd);
	    }
	  cout << "read " << MCevents.size() << " accepted MC events out of "
	       << MCallEvents.size() << " total ("
	       << 100.*MCevents.size() / MCallEvents.size()
	       << "% overall acceptance)"
	       << endl;
	}
      else
	{
	  MCevents.reserve(NFLATMCEVENTS);
	  for (int iMC = 0; iMC < NFLATMCEVENTS; iMC++)
	    {
	      event e(acos(gRandom->Uniform(-1,1)), gRandom->Uniform(-M_PI, M_PI));
	      MCevents.push_back(e);
	    }
	}

      myL.addChannel(RDevents, MCevents, MCallEvents);
    }

  if (myL.getNChannels() == 0)
    {
      cerr << "no data." << endl;
      abort();
    }
  cout << "combined fit of " << myL.getNChannels() << " channels." << endl;

  TStopwatch sw;

  TFitterMinuit* minuit = new TFitterMinuit();
  minuit->SetMinuitFCN(&myL);

  //TH3* hPredict = new TH3D("hPredict", "prediction", nBins, 0, nBins, 100, -1, 1, 100, -M_PI, M_PI);

  TH1* hBR = new TH1D("hBR", "relative Branching Ratio",
		      nBins, lower, upper);
  hBR->SetMinimum(0);

  const size_t nParams = lastIdx + myL.getNChannels();
  TStopwatch fulltime;
  fulltime.Start();
  for (iBin = 0; iBin < nBins; iBin++)
    {
      sw.Start();
      myL.setBin(iBin);

      minuit->Clear();
      if (!flatMC)
	// flat MCweights are the same in different mass bins for the
	// two-body amplitudes.
	myL.clearWeights();

      vector<double> vStartingValues(nParams);

      // Use random starting values if the user didn't demand
      // continuity between bins.
      if (!continuous)
	for (size_t iSV = 0; iSV < nParams; iSV++)
	  {
	    if (!startingValues[iSV].fixed)
	      startingValues[iSV].value = gRandom->Uniform(1);
	  }
      for (size_t iSV = 0; iSV < nParams; iSV++)
	{
	  vStartingValues[iSV] = startingValues[iSV].value;
	}

      // The MC part of the likelihood function will evaluate to the
      // number of events in the fit.  In order to speed up the
      // calculation we scale the starting values such that this
      // condition obtains.
      size_t nGood = myL.eventsInBin();
      if (nGood == 0)
	continue;
      double ratio = sqrt(nGood / myL.calc_mc_likelihood(vStartingValues));
      for (size_t iSV = 0; iSV < nParams - myL.getNChannels(); iSV++)
	{
	  vStartingValues[iSV] *= ratio;
	}

      // Set starting values.
      for (size_t j= 0; j < nParams - myL.getNChannels(); j++)
	{
	  if (!startingValues[j].fixed)
	    {
	      minuit->SetParameter(j, startingValues[j].name.c_str(),
				   vStartingValues[j], vStartingValues[j].value*0.01, 0, 0);
	    }
	  else
	    {
	      minuit->SetParameter(j, startingValues[j].name.c_str(),
				   vStartingValues[j], 1, 0, 0);
	      minuit->FixParameter(j);
	    }
	}
      for (size_t j = nParams - myL.getNChannels(); j < nParams; j++)
	{
	  minuit->SetParameter(j, startingValues[j].name.c_str(),
			       vStartingValues[j].value, .1, 0, 1);
	  if (startingValues[j].fixed)
	    minuit->FixParameter(j);
	}

      // Run minimizer.
      minuit->CreateMinimizer();
      int iret = minuit->Minimize();
      sw.Stop();
      cout << "iret = " << iret << " after " << sw.CpuTime() << " s." << endl;

      if (iret == 0)
	{
	  vector<double> vStartingValue(nParams);
	  for (size_t j = 0; j < nParams; j++)
	    {
	      vStartingValues[j] = minuit->GetParameter(j);
	    }

	  for (int j = 0; j < minuit->GetNumberTotalParameters(); j++)
	    startingValues[j].value = minuit->GetParameter(j);

	  gHist.getHist("hMClikelihood", "MC likelihood of result", nBins, lower, upper)
	    ->SetBinContent(iBin+1, myL.calc_mc_likelihood(vStartingValues));

	  for (size_t iCoherent = 0; iCoherent < ws.size(); iCoherent++)
	    {
	      vector<wave>& waves = ws[iCoherent].getWaves();
	      for (size_t iWave1 = 0; iWave1 < waves.size(); iWave1++)
		{
		  waves[iWave1].fillHistIntensity(iBin, minuit);
		  if (iWave1 != waves.size()-1)
		    {
		      for (size_t iWave2 = iWave1 + 1; iWave2 < waves.size(); iWave2++)
			waves[iWave1].fillHistPhase(iBin, waves[iWave2], minuit);
		    }
		}
	    }

	  if (myL.getNChannels() == 2)
	    {
	      hBR->SetBinContent(iBin+1, minuit->GetParameter(16 + 1));
	      hBR->SetBinError(iBin+1, minuit->GetParError(16 + 1));
	    }

	  vector<double> result;
	  for (int iPar = 0; iPar < minuit->GetNumberTotalParameters(); iPar++)
	    {
	      result.push_back(minuit->GetParameter(iPar));
	    }


	  for (std::vector<std::pair<size_t, size_t> >::const_iterator it = vecMom.begin();
	       it != vecMom.end(); it++)
	    {
	      mhMoments[*it]->SetBinContent(iBin + 1, decomposeMoment(*it, ws, vStartingValues));
	      mhMoments[*it]->SetBinError(iBin + 1, decomposeMomentError(*it, ws, minuit));
	    }

#if 0
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
#endif
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


int main(int argc, char **argv)
{
  const char *controlfn;
  if (argc == 1)
    controlfn = "control.txt";
  else
    controlfn = argv[1];
  if (!readControlFile(controlfn))
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
      gHist.clear();
    }

  return 0;
}
