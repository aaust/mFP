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
  TH1* hMassFine = new TH1D("hMassFine", "mass distribution",
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
      fd = fopen(MCFile.c_str(), "r");
      if (!fd)
	{
	  cerr << "Can't open input file '" << MCFile << "'." << endl;
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
	// flat MCweights are the same in different mass bins for the
	// two-body amplitudes.
	myL.clearWeights();

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
