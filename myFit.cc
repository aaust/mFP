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
#include "TFitterMinuit.h"
#include "TStopwatch.h"
#include "TFile.h"

#include "control.h"
#include "wave.h"
#include "event.h"
#include "likelihood.h"

#define NFLATMCEVENTS 100000
double massLow = 0;
double massHigh = 9999;
int iBin;

#define NWAVES 8

class combinedLikelihood : public ROOT::Minuit2::FCNBase {
public:
  vector<likelihood> myLs;
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

  void
  addChannel(vector<event>& RDevents,
	     vector<event>& MCevents,
	     vector<event>& MCallEvents)
  {
    size_t idxBranching = 2*NWAVES + this->getNChannels();
    myLs.push_back(likelihood(ws, RDevents, MCevents, MCallEvents, nBins, threshold, binWidth, idxBranching));
  }


  double Up() const { return 0.5; }

  void setBin(size_t iBin) { for (size_t i = 0; i < myLs.size(); i++) myLs[i].setBin(iBin); }
  size_t eventsInBin() const {
    size_t sum = 0;
    for (size_t i = 0; i < myLs.size(); i++) sum += myLs[i].eventsInBin(); 
    return sum;
  }
  void clearWeights() { for (size_t i = 0; i < myLs.size(); i++) myLs[i].clearWeights(); }

  double
  calc_mc_likelihood(const vector<double>& x) const
  {
    double result = 0;
    for (size_t i = 0; i < myLs.size(); i++)
      result += myLs[i].calc_mc_likelihood(x);
    return result;
  }

  double
  operator()(const vector<double>& x) const
  {
    double result = 0;
    for (size_t i = 0; i < myLs.size(); i++)
      result -= myLs[i].calc_likelihood(x);
    return result;
  }

  size_t getNChannels() { return myLs.size(); }
};



void
myFit()
{
  gRandom = new TRandom1;

  vector<wave> positive;
  positive.push_back(wave("D+", 2, 1));
  positive.push_back(wave("P+", 1, 1));
  positive.push_back(wave("G+", 4, 1));

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
  startingValues[2 * NWAVES + 2] =
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
      { "BR1", 1, true },
      { "BR2", 0.6, false },
    };

  TH2* hRD = new TH2D("hRD", "RD", 10, -1, 1, 10, -M_PI, M_PI);
  TH1* hMassFine = new TH1D("hMassFine", "mass distribution",
			    250, 0.5, 3);
  TH1* htprime = new TH1D("htprime", "t' distribution",
			  250, 0, 1);
  TH1* hMassMC = new TH1D("hMassMC", "MC mass distribution",
			250, 0.5, 3);
  TH1* htprimeMC = new TH1D("htprimeMC", "t' distribution",
			  250, 0, 1);

  combinedLikelihood myL(ws, nBins, threshold, binWidth);

  for (size_t iFile = 0; iFile < dataFiles.size(); iFile++)
    {
      FILE* fd = fopen(dataFiles[iFile].c_str(), "r");
      if (!fd)
	{
	  cerr << "Can't open input file '" << dataFiles[iFile] << "'." << endl;
	  abort();
	}
      char line[99999];

      vector<event> RDevents;
      while (fgets(line, 99999, fd) /*&& nrdevents < 5000*/)
	{
	  double m, tPr, theta, phi;
	  sscanf(line, "%lf %lf %lf %lf", &m, &tPr, &theta, &phi);

	  hMassFine->Fill(m);
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
	  fd = fopen(MCFiles[iFile].c_str(), "r");
	  if (!fd)
	    {
	      cerr << "Can't open input file '" << MCFiles[iFile] << "'." << endl;
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
  TH1* hGwave = new TH1D("hGwave", "G wave intensity",
			 nBins, lower, upper);
  TH1* hPhaseDP = new TH1D("hPhaseDP", "D - P phase",
			 nBins, lower, upper);
  TH1* hPhaseDG = new TH1D("hPhaseDG", "D - G phase",
			 nBins, lower, upper);
  TH1* hPhaseD0S = new TH1D("hPhaseD0S", "D0 - S phase",
			    nBins, lower, upper);
  TH1* hIntensity = new TH1D("hIntensity", "total intensity as predicted",
			     nBins, lower, upper);
  //TH3* hPredict = new TH3D("hPredict", "prediction", nBins, 0, nBins, 100, -1, 1, 100, -M_PI, M_PI);

  TH1* hBR = new TH1D("hBR", "relative Branching Ratio",
		      nBins, lower, upper);

  int listOfMoments[][2] = { { 0, 0 },
			     { 1, 0 },
			     { 1, 1 },
			     { 2, 0 },
			     { 2, 1 },
			     { 2, 2 },
			     { 3, 0 },
			     { 3, 1 },
			     { 3, 2 },
			     { 4, 0 },
			     { 4, 1 },
			     { 4, 2 } };
  size_t nMoments = 12;
  vector<TH1D*> hMomentsRe(nMoments);
  vector<TH1D*> hMomentsIm(nMoments);
  vector<TH1D*> hMomentsPWA(nMoments);
  for (size_t i = 0; i < nMoments; i++)
    {
      int L = listOfMoments[i][0];
      int M = listOfMoments[i][1];

      char name[999], title[999];
      snprintf(name, 999, "hH%d%dRe", L, M);
      snprintf(title, 999, "Re(H(%d%d)) from experiment", L, M);
      hMomentsRe[i] = new TH1D(name, title, nBins, lower, upper);

      snprintf(name, 999, "hH%d%dIm", L, M);
      snprintf(title, 999, "Im(H(%d%d)) from experiment", L, M);
      hMomentsIm[i] = new TH1D(name, title, nBins, lower, upper);

      snprintf(name, 999, "hH%d%dPWA", L, M);
      snprintf(title, 999, "hH%d%d from PWA", L, M);
      hMomentsPWA[i] = new TH1D(name, title, nBins, lower, upper);
    }

  const size_t nParams = 2*NWAVES + myL.getNChannels();
  TStopwatch fulltime;
  fulltime.Start();
  for (iBin = 0; iBin < nBins; iBin++)
    {
      sw.Start();

      myL.setBin(iBin);

      // Calculate moments, fill in hist.
      for (size_t i = 0; i < nMoments; i++)
	{
	  int L = listOfMoments[i][0];
	  int M = listOfMoments[i][1];
	  complex<double> mom = myL.myLs[0].calcMoment(L, M);
	  hMomentsRe[i]->SetBinContent(iBin+1, real(mom));
	  hMomentsIm[i]->SetBinContent(iBin+1, imag(mom));
	}

      minuit->Clear();
      if (!flatMC)
	// flat MCweights are the same in different mass bins for the
	// two-body amplitudes.
	myL.clearWeights();

      // The MC part of the likelihood function will evaluate to the
      // number of events in the fit.  In order to speed up the
      // calculation we scale the starting values such that this
      // condition obtains.
      vector<double> vStartingValues(nParams);
      for (size_t iSV = 0; iSV < nParams; iSV++)
	{
	  vStartingValues[iSV] = startingValues[iSV].value;
	}
      size_t nGood = myL.eventsInBin();
      if (nGood == 0)
	continue;
      double ratio = sqrt(nGood / myL.calc_mc_likelihood(vStartingValues));
      for (size_t iSV = 0; iSV < nParams - myL.getNChannels(); iSV++)
	{
	  vStartingValues[iSV] *= ratio;
	}

      for (size_t j= 0; j < nParams - myL.getNChannels(); j++)
	{
	  if (!startingValues[j].fixed)
	    minuit->SetParameter(j, startingValues[j].name,
				 startingValues[j].value*ratio, 1, 0, 0);
	  else
	    {
	      minuit->SetParameter(j, startingValues[j].name,
				   startingValues[j].value, 1, 0, 0);
	      minuit->FixParameter(j);
	    }
	}
      for (size_t j = nParams - myL.getNChannels(); j < nParams; j++)
	{
	  minuit->SetParameter(j, startingValues[j].name,
			       startingValues[j].value, .1, 0, 0);
	  if (startingValues[j].fixed)
	    minuit->FixParameter(j);
	}

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

	  complex<double> aDwave(minuit->GetParameter(0),
				 minuit->GetParameter(1));
	  complex<double> aPwave(minuit->GetParameter(2),
				 minuit->GetParameter(3));
	  complex<double> aGwave(minuit->GetParameter(4),
				 minuit->GetParameter(5));
	  complex<double> aSwave(minuit->GetParameter(6),
				 minuit->GetParameter(7));
	  complex<double> aP0wave(minuit->GetParameter(8),
				  minuit->GetParameter(9));
	  complex<double> aPmwave(minuit->GetParameter(10),
				  minuit->GetParameter(11));
	  complex<double> aD0wave(minuit->GetParameter(12),
				  minuit->GetParameter(13));
	  complex<double> aDmwave(minuit->GetParameter(14),
				  minuit->GetParameter(15));

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

	  error = 2*(sqrt(pow(real(aGwave) * minuit->GetParError(4), 2)
			       + pow(imag(aGwave) * minuit->GetParError(5), 2)
			       + (2*real(aGwave)*imag(aGwave)
				  * minuit->GetCovarianceMatrixElement(4, 5))));
	  hGwave->SetBinContent(iBin+1,norm(aGwave));
	  hGwave->SetBinError(iBin+1, error);


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

	  phase = arg(aDwave / aGwave);
	  if (iBin > 1)
	    {
	      double oldPhase = hPhaseDG->GetBinContent(iBin);
	      if (phase - oldPhase > M_PI)
		phase -= 2*M_PI;
	      else if (oldPhase - phase > M_PI)
		phase += 2*M_PI;
	    }
	  hPhaseDG->SetBinContent(iBin+1,phase);
	  hPhaseDG->SetBinError(iBin+1, .2);

	  //hBR->SetBinContent(iBin+1, minuit->GetParameter(2*NWAVES + 1));
	  //hBR->SetBinError(iBin+1, minuit->GetParError(2*NWAVES + 1));

	  vector<double> result;
	  for (int iPar = 0; iPar < minuit->GetNumberTotalParameters(); iPar++)
	    {
	      result.push_back(minuit->GetParameter(iPar));
	    }

	  const double s2 = sqrt(2);
	  const double s3 = sqrt(3);
	  const double s5 = sqrt(5);
	  const double s6 = sqrt(6);
	  const double s10 = sqrt(10);
	  const double s15 = sqrt(15);
	  const double s21 = sqrt(21);
	  const double s30 = sqrt(30);
	  double H00 = (norm(aSwave) + norm(aP0wave) + norm(aPmwave)
			+ norm(aD0wave) + norm(aDmwave) + norm(aPwave)
			+ norm(aDwave));
	  double H10 = real(2/s3*aSwave*conj(aP0wave)
			    + 4/s15*aP0wave*conj(aD0wave)
			    + 2/s5*aPmwave*conj(aDmwave)
			    + 2/s5*aPwave*conj(aDwave));
	  double H11 = real(2/s6*aSwave*conj(aPmwave)
			    + 2/s10*aP0wave*conj(aDmwave)
			    - 2/s30*aPmwave*conj(aD0wave));
	  double H20 = (2/s5*real(aSwave*conj(aD0wave))
			+ .4*norm(aP0wave)
			- .2*norm(aPmwave)
			- .2*norm(aPwave)
			+ 2./7.*norm(aD0wave)
			+ 1./7.*norm(aDmwave)
			+ 1./7.*norm(aDwave));
	  double H21 = real(2/s10*aSwave*conj(aDmwave)
			    + .4*s3/s2*aP0wave*conj(aPmwave)
			    + 2./7./s2*aD0wave*conj(aDmwave));
	  double H22 = (.2*s3/s2*norm(aPmwave)
			+ .2*s3/s2*norm(aPwave)
			+ 1./7.*s3/s2*norm(aDmwave)
			- 1./7.*s3/s2*norm(aDwave));
	  double H30 = real(6./7.*s3/s5*aP0wave*conj(aD0wave)
			    - 6./7./s5*aPmwave*conj(aDmwave)
			    - 6./7./s5*aPwave*conj(aDwave));
	  double H31 = real(4./7.*s3/s5*aP0wave*conj(aDmwave)
			    + 6./7./s5*aPmwave*conj(aD0wave));
	  double H32 = real(2./7.*s3/s2*aPmwave*conj(aDmwave)
			    - 2./7.*s3/s2*aPwave*conj(aDwave));
	  double H40 = (2./7.*norm(aD0wave)
			- 4./21.*norm(aDmwave)
			- 4./21.*norm(aDwave));
	  double H41 = 2./7.*s5/s3*real(aD0wave*conj(aDmwave));
	  double H42 = s10/s21*norm(aDmwave) - s10/s21*norm(aDwave);

	  double vH[] = { H00, H10, H11, H20, H21, H22,
			  H30, H31, H32, H40, H41, H42, };
	  for (size_t iMom = 0; iMom < 12; iMom++)
	    {
	      hMomentsPWA[iMom]->SetBinContent(iBin+1, vH[iMom]);
	      hMomentsPWA[iMom]->SetBinError(iBin+1, .1);
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
