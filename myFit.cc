#include <complex>
#include <vector>
#include <map>
#include <algorithm>

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
#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "control.h"
#include "wave.h"
#include "3j.h"
#include "event.h"
#include "likelihood.h"
#include "startingValue.h"
#include "fitInfo.h"
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

  void
  fillPredict(const vector<double>& x, TH2* hth, TH2* hph ) const
  {
    cout << "Calculating weighted MC" << endl;

    for (size_t i = 0; i < myLs.size(); i++)
      myLs[i]->fillPredict(x, hth, hph);
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
	//cout << "nevents = " << myLs[i]->eventsInBin() << " likelihood = " << myLs[i]->calc_rd_likelihood(x) << " - " << myLs[i]->calc_mc_likelihood(x) << endl;
      }
    return -result;
  }

  size_t getNChannels() const { return myLs.size(); }
  size_t MCeventsInBin() const { size_t sum = 0; for (size_t i = 0; i < myLs.size(); i++) sum += myLs[i]->MCeventsInBin(); return sum; }
};


// from : Numerical Recipes in C: The Art of Scientific Computing (Second Edition), published by Cambridge University Press, SECTION 9.5

void laguer(const vector<complex<double> > &a, complex<double> &x, int &its)
{
  const int MR=8,MT=10,MAXIT=MT*MR;
	const double EPS=numeric_limits<double>::epsilon();
	const double frac[MR+1]={0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
	complex<double> dx,x1,b,d,f,g,h,sq,gp,gm,g2;
	int m=a.size()-1;
	for (int iter=1;iter<=MAXIT;iter++) {
		its=iter;
		b=a[m];
		double err=abs(b);
		d=f=0.0;
		double abx=abs(x);
		for (int j=m-1;j>=0;j--) {
			f=x*f+d;
			d=x*d+b;
			b=x*b+a[j];
			err=abs(b)+abx*err;
		}
		err *= EPS;
		if (abs(b) <= err) return;
		g=d/b;
		g2=g*g;
		h=g2-2.0*f/b;
		sq=sqrt(double(m-1)*(double(m)*h-g2));
		gp=g+sq;
		gm=g-sq;
		double abp=abs(gp);
		double abm=abs(gm);
		if (abp < abm) gp=gm;
		dx=max(abp,abm) > 0.0 ? double(m)/gp : polar(1+abx,double(iter));
		x1=x-dx;
		if (x == x1) return;
		if (iter % MT != 0) x=x1;
		else x -= frac[iter/MT]*dx;
	}
	throw("too many iterations in laguer");
};

void zroots(const vector<complex<double> > &a, vector<complex<double> > &roots, const bool &polish)
{
	const double EPS=1.0e-14;
	int i,its;
	complex<double> x,b,c;
	int m=a.size()-1;
	vector<complex<double> > ad(m+1);
	for (int j=0;j<=m;j++) ad[j]=a[j];
	for (int j=m-1;j>=0;j--) {
		x=0.0;
		vector<complex<double> > ad_v(j+2);
		for (int jj=0;jj<j+2;jj++) ad_v[jj]=ad[jj];
		laguer(ad_v,x,its);
		if (abs(imag(x)) <= 2.0*EPS*abs(real(x)))
		  x=complex<double>(real(x),0.0);
		roots[j]=x;
		b=ad[j+1];
		for (int jj=j;jj>=0;jj--) {
			c=ad[jj];
			ad[jj]=b;
			b=x*b+c;
		}
	}
	if (polish)
		for (int j=0;j<m;j++)
			laguer(a,roots[j],its);
	for (int j=1;j<m;j++) {
		x=roots[j];
		for (i=j-1;i>=0;i--) {
			if (real(roots[i]) <= real(x)) break;
			roots[i+1]=roots[i];
		}
		roots[i+1]=x;
	}
};

// Sort roots

void sortroots(vector<complex<double> > &roots)
{ 
  if ( real(roots[2]) > real(roots[0]) ) swap(roots[0], roots[2]);
  if ( real(roots[3]) < real(roots[1]) ) swap(roots[1], roots[3]);

  if ( imag(roots[0]) < 0 ) roots[0] = conj(roots[0]);
  if ( imag(roots[3]) < 0 ) roots[3] = conj(roots[3]);
  if ( imag(roots[2]) > 0 ) roots[2] = conj(roots[2]);
  if ( imag(roots[1]) > 0 ) roots[1] = conj(roots[1]);

  if ( imag(roots[0]) > imag(roots[3]) ) swap(roots[0], roots[3]);
  //if ( imag(roots[2]) < imag(roots[1]) ) swap(roots[2], roots[1]);

}
  
// from : Techniques of Amplitude Analysis for two-pseudoscalar systems, S.U. Chung, Phys.Rev.D 56, 1997, 7299

void waves2coeff(vector<double> &values, vector<complex<double> > &coeff)
{
   complex<double> S0(values[4], values[5]);
   complex<double> P0(values[6], values[7]);
   complex<double> PM(values[8], values[9]);
   complex<double> D0(values[10],values[11]);
   complex<double> DM(values[12],values[13]);

   complex<double> a0 = S0 + sqrt(3)*P0 + sqrt(5)*D0;
   complex<double> a1 = -2.*sqrt(3)*(PM + sqrt(5)*DM);
   complex<double> a2 = 2.*S0 - 4.*sqrt(5)*D0;
   complex<double> a3 = -2.*sqrt(3)*(PM - sqrt(5)*DM);
   complex<double> a4 = S0 - sqrt(3)*P0 + sqrt(5)*D0;
   
   coeff.push_back(a0);
   coeff.push_back(a1);
   coeff.push_back(a2);
   coeff.push_back(a3);
   coeff.push_back(a4);

};


void roots2waves(complex<double> &a4, vector<complex<double> > &roots, vector<double> &waves)
{
  complex<double> u1 = roots[0];
  complex<double> u2 = roots[1];
  complex<double> u3 = roots[2];
  complex<double> u4 = roots[3];
  
  complex<double> s0 = a4/6. * (2.*u1*u2*u3*u4 + u1*u2 + u1*u3 + u1*u4 + u2*u3 + u2*u4 + u3*u4 + 2.);
  complex<double> p0 = a4/(2.*sqrt(3)) * (u1*u2*u3*u4 - 1.);
  complex<double> pm = a4/(4.*sqrt(3)) * (u1*u2*u3 + u2*u3*u4 + u3*u4*u1 + u4*u1*u2 + u1 + u2 + u3 + u4);
  complex<double> d0 = a4/(6.*sqrt(5)) * (u1*u2*u3*u4 - u1*u2 - u1*u3 - u1*u4 - u2*u3 - u2*u4 - u3*u4 + 1.);
  complex<double> dm = a4/(4.*sqrt(15))* (u1*u2*u3 + u2*u3*u4 + u3*u4*u1 + u4*u1*u2 - u1 - u2 - u3 - u4);

  // rotate, such that im(s0) == 0
  // might change sign!
  complex<double> im (0.,-1.);
  //  cout << s0 << " " << abs(s0) << endl;
  //s0 *= exp(im*arg(s0));
  //cout << s0 << " " << abs(s0) << endl;
  /*  p0 *= exp(im*arg(s0));
  pm *= exp(im*arg(s0));
  d0 *= exp(im*arg(s0));
  dm *= exp(im*arg(s0));*/
  
  waves.push_back(real(s0));
  //waves.push_back(0);
  waves.push_back(imag(s0));
  waves.push_back(real(p0));
  waves.push_back(imag(p0));
  waves.push_back(real(pm));
  waves.push_back(imag(pm));
  waves.push_back(real(d0));
  waves.push_back(imag(d0));
  waves.push_back(real(dm));
  waves.push_back(imag(dm));
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
  gHist.Fill("hMassCoarse", "expected MC likelihood of result",
	     nBins, lower, upper,
	     e.mass);

  gHist.Fill("hMvsCosth", "m vs costh Rd", nBins, lower, upper, 100, -1, 1,
	     e.mass, cos(e.theta));
  if (e.accepted())
    gHist.Fill("hMvsPhi", "m vs phi Rd", nBins, lower, upper, 100, -M_PI, M_PI,
	       e.mass, e.phi);
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
  TRandom1* gRandom = new TRandom1();
  gRandom->SetSeed2(0); // truely random!!!

  double lower = threshold;
  double upper = threshold + nBins*binWidth;

  vector<wave> positive;
  positive.push_back(wave("P+", 1, 1, nBins, lower, upper, true));
  positive.push_back(wave("D+", 2, 1, nBins, lower, upper));
  //positive.push_back(wave("F+", 3, 1, nBins, lower, upper));
  //positive.push_back(wave("G+", 4, 1, 78, 1.7, upper));
  //positive.push_back(wave("D++", 2, 2, nBins, lower, upper));

  vector<wave> negative;
  negative.push_back(wave("S0", 0, 0, nBins, lower, upper, true));
  negative.push_back(wave("P0", 1, 0, nBins, lower, upper));
  negative.push_back(wave("P-", 1, 1, nBins, lower, upper));
  negative.push_back(wave("D0", 2, 0, nBins, lower, upper));
  negative.push_back(wave("D-", 2, 1, nBins, lower, upper));
  //negative.push_back(wave("G0", 4, 0, 78, 1.7, upper));
  //negative.push_back(wave("G-", 4, 1, 78, 1.7, upper));
 
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

  TTree *outTree = new TTree("tFitResults", "fit results tree");
  double *values = new double[lastIdx];
  TMatrixDSym *covMat = new TMatrixDSym(lastIdx);
  UInt_t NMCevents;
  char branchDesc[99];
  snprintf(branchDesc, 99, "values[%zd]/D", lastIdx);
  outTree->Branch("massLow", &massLow, "massLow/D");
  outTree->Branch("massHigh", &massHigh, "massHigh/D");
  outTree->Branch("values", values, branchDesc);
  outTree->Branch("covMat", "TMatrixDSym", &covMat);
  outTree->Branch("NMCevents", &NMCevents, "NMCevents/i");
  outTree->GetUserInfo()->Add(new fitInfo(modelName, ws, startingValues, nBins, threshold, binWidth, lastIdx));

  /*  TTree* predictTree = new TTree("predictTree","Weighted MC");
  double weight(0);
  float mX(0), mc_theta(0), mc_phi(0), mc_costh(0);
  if (mc){
    out->Branch("accepted", &mc_acc, "mc_acc/I");
    out->Branch("mX", &mc_two, "mc_two/F");
    out->Branch("tPrime", &tmax, "tmax/F");
    out->Branch("costhGJ", &mc_costh, "mc_costh/F");
    out->Branch("phiGJ", &mc_phi, "mc_phi/F");
  }
  */

  // for second fit with calculated solution, take positive result from last bin ('continuous')
  vector<double> positiveStart(4);
  positiveStart.push_back(gRandom->Uniform(5));
  positiveStart.push_back(gRandom->Uniform(5));
  positiveStart.push_back(gRandom->Uniform(5));
  positiveStart.push_back(gRandom->Uniform(5));

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
	      //tree->SetBranchAddress("likeK", &likeK);
	      //tree->SetBranchAddress("likePi", &likePi);

	      RDevents.reserve(tree->GetEntries());
	      for (Long_t i = 0; i < tree->GetEntries(); i++)
		{
		  tree->GetEntry(i);
		  //if (likeK != -1 && likePi > 2*likeK)
		  //continue;
		  event e(mX, tPr, theta, phi);
		  RDevents.push_back(e);
		  fillRDhists(e);
		}
	    }
	  else
	    {
	      f->GetObject("schluter/trRDEtap3pi", tree);
	      if (tree)
		{
		  float m;
		  float costh;
		  float phi;
		  float t;
		  float mCandEtaP1, mCandEtaP2;
		  float CL1, CL2;
		  float pG1X, pG1Y, pG1Z, pG1E;
		  float pG2X, pG2Y, pG2Z, pG2E;

		  tree->SetBranchAddress("m", &m);
		  tree->SetBranchAddress("t", &t);
		  tree->SetBranchAddress("costh", &costh);
		  tree->SetBranchAddress("phi", &phi);
		  tree->SetBranchAddress("mCandEtaP1", &mCandEtaP1);
		  tree->SetBranchAddress("mCandEtaP2", &mCandEtaP2);
		  tree->SetBranchAddress("CL1", (int*)&CL1); // mistakes happen -> int*
		  tree->SetBranchAddress("CL2", (int*)&CL2);
		  tree->SetBranchAddress("pG1X", &pG1X);
		  tree->SetBranchAddress("pG1Y", &pG1Y);
		  tree->SetBranchAddress("pG1Z", &pG1Z);
		  tree->SetBranchAddress("pG1E", &pG1E);
		  tree->SetBranchAddress("pG2X", &pG2X);
		  tree->SetBranchAddress("pG2Y", &pG2Y);
		  tree->SetBranchAddress("pG2Z", &pG2Z);
		  tree->SetBranchAddress("pG2E", &pG2E);

		  RDevents.reserve(tree->GetEntries());
		  for (Long_t i = 0; i < tree->GetEntries(); i++)
		    {
		      tree->GetEntry(i);

		      TLorentzVector lvEta(pG1X + pG2X, pG1Y + pG2Y, pG1Z + pG2Z, pG1E + pG2E);

		      gHist.Fill("hmEtap", "m(#pi#pi#eta)", 1000, 0.7, 1.7, mCandEtaP1);
		      gHist.Fill("hmEtap", "m(#pi#pi#eta)", 1000, 0.7, 1.7, mCandEtaP2);
		      gHist.Fill("hCL1pre", "CL1 pre mass cut", 200, 0, 1, CL1);
		      gHist.Fill("hCL2pre", "CL2 pre mass cut", 200, 0, 1, CL2);

		      gHist.Fill("hmEtaVsCL2pre", "m(#gamma#gamma) vs. CL2", 100, 0.5, 0.62, 400, 0, 1,
				 lvEta.M(), CL2);

		      gHist.Fill("hmEta", "m(#gamma#gamma)", 400, 0.5, 0.62, lvEta.M());
		      gHist.Fill("hmEtaVsE", "m(#gamma#gamma) vs E", 200, 0.5, 0.62, 100, 0, 200,
				 lvEta.M(), lvEta.E());
		      if (fabs(mCandEtaP1 - .958) < 0.02 || fabs(mCandEtaP2 - .958) < 0.02 /*CL2 > 0.01*/)
			{
			  gHist.Fill("hmEtaPostCL", "m(#gamma#gamma)", 400, 0.5, 0.62, lvEta.M());
			  gHist.Fill("hmEtaVsEPostCL", "m(#gamma#gamma) vs E", 200, 0.5, 0.62, 100, 0, 200,
				     lvEta.M(), lvEta.E());
			  //continue;
			}
		      gHist.Fill("hmEtappost", "m(#pi#pi#eta) post CL cut", 1000, 0.7, 1.7, mCandEtaP1);
		      gHist.Fill("hmEtappost", "m(#pi#pi#eta) post CL cut", 1000, 0.7, 1.7, mCandEtaP2);

		      if (fabs(mCandEtaP1 - .958) > 0.02 && fabs(mCandEtaP2 - .958) > 0.02)
			continue;
		      gHist.Fill("hCL1", "CL1 post mass cut", 200, 0, 1, CL1);
		      gHist.Fill("hCL2", "CL2 post mass cut", 200, 0, 1, CL2);
		      gHist.Fill("hmEtaVsCL2", "m(#gamma#gamma) vs. CL2 post mass cut", 100, 0.5, 0.62, 400, 0, 1,
				 lvEta.M(), CL2);
		      event e(m, -t, acos(costh), phi);
		      RDevents.push_back(e);
		      fillRDhists(e);
		    }
		}
	      else
		{
		  f->GetObject("tPredict", tree);
		  if (!tree)
		    {
		      cerr << "no known tree found" << endl;
		      abort();
		    }

		  float m;
		  float costh;
		  float phi;
		  float t;
		  tree->SetBranchAddress("m", &m);
		  tree->SetBranchAddress("t", &t);
		  tree->SetBranchAddress("costh", &costh);
		  tree->SetBranchAddress("phi", &phi);
 
		  RDevents.reserve(tree->GetEntries());
		  for (Long_t i = 0; i < tree->GetEntries(); i++)
		    {
		      tree->GetEntry(i);

		      event e(m, t, acos(costh), phi);
		      RDevents.push_back(e);
		      fillRDhists(e);
		    }
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
	      /*	      cout << "reserve " << t->GetEntries()/10 << " accepted MC events out of "
	       << t->GetEntries() << " total"
	       << endl;*/

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
	       << 100LL*MCevents.size() / MCallEvents.size()
	       << "% overall acceptance)"
	       << endl;
	}
      else
	{
	  MCevents.reserve(NFLATMCEVENTS);
	  for (int iMC = 0; iMC < NFLATMCEVENTS; iMC++)
	    {
	      double costh = gRandom->Uniform(-1,1);
	      double phi = gRandom->Uniform(-M_PI, M_PI);
	      //cout << costh << " " << acos(costh) << " " << phi << endl;
	      event e(acos(costh), phi);
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
  const size_t nParams = lastIdx + myL.getNChannels();

  TFitterMinuit* minuit = new TFitterMinuit();
  minuit->SetMinuitFCN(&myL);

  TFitterMinuit* fitamb = new TFitterMinuit();
  fitamb->SetMinuitFCN(&myL);
  
  // Do not show fit results:
  minuit->SetPrintLevel(0);
  fitamb->SetPrintLevel(0);

  //TH3* hPredict = new TH3D("hPredict", "prediction", nBins, 0, nBins, 100, -1, 1, 100, -M_PI, M_PI);

  TH2* hthpre = new TH2D("hthpre", "prediction for cos(#theta)", nBins, threshold, threshold + nBins*binWidth, 100, -1, 1);
  TH2* hphpre = new TH2D("hphpre", "prediction for #phi", nBins, threshold, threshold + nBins*binWidth, 100, -M_PI, M_PI);

  TH1* hBR = new TH1D("hBR", "relative Branching Ratio",
		      nBins, lower, upper);
  hBR->SetMinimum(0);

  TStopwatch fulltime;
  int failBins = 0;
  fulltime.Start();
  // possibiliy to start from high mass end
  //for (iBin = nBins-1; iBin >= 0; iBin--)
  for (iBin = 0; iBin < nBins; iBin++)
    {
      
      cout << "Bin " << iBin << " : " << lower+(iBin*binWidth) << "GeV" << endl;

      sw.Start();
      myL.setBin(iBin);

      minuit->Clear();
      fitamb->Clear();
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
	      startingValues[iSV].value = gRandom->Uniform(5);
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
      double ratio = (vStartingValues[0] < 0 ? -1 : 1)*sqrt(nGood / myL.calc_mc_likelihood(vStartingValues));
      for (size_t iSV = 0; iSV < nParams - myL.getNChannels(); iSV++)
	{
	  vStartingValues[iSV] *= ratio;
	}

      // Set starting values.
      for (size_t j= 0; j < nParams - myL.getNChannels(); j++)
	{
	  if (!startingValues[j].fixed)
	    {
	      //if (j != 0)
		minuit->SetParameter(j, startingValues[j].name.c_str(),
				   vStartingValues[j], 10/*vStartingValues[j]*0.01*/, 0, 0);
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
			       vStartingValues[j], .1, 0, 1);
	  if (startingValues[j].fixed)
	    minuit->FixParameter(j);
	}

      // Run minimizer.
      minuit->CreateMinimizer();
      int iret = minuit->Minimize();

      while ( iret != 0 ){
	// if fit did not converge, use random starting values until it does
	for (size_t iSV = 0; iSV < nParams; iSV++)
	  {
	    if (!startingValues[iSV].fixed)
	      startingValues[iSV].value = gRandom->Uniform(5);
	  }
	for (size_t iSV = 0; iSV < nParams; iSV++)
	  {
	    vStartingValues[iSV] = startingValues[iSV].value;
	  }

	// Set starting values.
	for (size_t j= 0; j < nParams - myL.getNChannels(); j++)
	  {
	    if (!startingValues[j].fixed)
	      {
		//if (j != 0)
		minuit->SetParameter(j, startingValues[j].name.c_str(),
				     vStartingValues[j], 10/*vStartingValues[j]*0.01*/, 0, 0);
	      }
	    else
	      {
		minuit->SetParameter(j, startingValues[j].name.c_str(),
				     vStartingValues[j], 1, 0, 0);
		minuit->FixParameter(j);
	      }
	  }
	
	iret = minuit->Minimize();
      }
      sw.Stop();
      //      cout << "iret = " << iret << " after " << sw.CpuTime() << " s." << endl;

      if (iret == 0)
	{
	  //	  vector<double> vStartingValue(nParams);
	  
	  for (size_t j = 0; j < nParams; j++)
	    {
	      vStartingValues[j] = minuit->GetParameter(j);
	    }

	  for (size_t i = 0; i < ws.size(); i++)
	    for (size_t j = 0; j < ws[i].waves.size(); j++)
	      {
		size_t idx = ws[i].waves[j].getIndex();
		values[idx] = minuit->GetParameter(idx);
		values[idx+1] = minuit->GetParameter(idx+1);
	      }

	  if (!ambiguous){
	    for (size_t i1 = 0; i1 < ws.size(); i1++)
	      for (size_t j1 = 0; j1 < ws[i1].waves.size(); j1++)
		{
		  const wave& w1 = ws[i1].waves[j1];
		  size_t idx1 = w1.getIndex();
		  size_t idx1Cov = w1.idxInCovariance(minuit);
		  
		  for (size_t i2 = 0; i2 < ws.size(); i2++)
		    for (size_t j2 = 0; j2 < ws[i2].waves.size(); j2++)
		      {
			const wave& w2 = ws[i2].waves[j2];
			size_t idx2 = w2.getIndex();
			size_t idx2Cov = w2.idxInCovariance(minuit);
			
			// Four possibilities:
			// 1) both free -> simply copy values to covMat
			if (!w1.phaseLocked && !w2.phaseLocked)
			  {
			    (*covMat)(idx1,idx2) = minuit->GetCovarianceMatrixElement(idx1Cov, idx2Cov);
			    (*covMat)(idx1,idx2+1) = minuit->GetCovarianceMatrixElement(idx1Cov, idx2Cov+1);
			    (*covMat)(idx1+1,idx2) = minuit->GetCovarianceMatrixElement(idx1Cov+1, idx2Cov);
			    (*covMat)(idx1+1,idx2+1) = minuit->GetCovarianceMatrixElement(idx1Cov+1, idx2Cov+1);
			  }
			// 2) first free, second fixed -> covariances with Im(w1) are zero
			else if (!w1.phaseLocked && w2.phaseLocked)
			  {
			    (*covMat)(idx1,idx2) = minuit->GetCovarianceMatrixElement(idx1Cov, idx2Cov);
			    (*covMat)(idx1,idx2+1) = 0;
			    (*covMat)(idx1+1,idx2) = minuit->GetCovarianceMatrixElement(idx1Cov+1, idx2Cov);
			    (*covMat)(idx1+1,idx2+1) = 0;
			  }
			// 3) vice versa
			else if (w1.phaseLocked && !w2.phaseLocked)
			  {
			    (*covMat)(idx1,idx2) = minuit->GetCovarianceMatrixElement(idx1Cov, idx2Cov);
			    (*covMat)(idx1,idx2+1) = minuit->GetCovarianceMatrixElement(idx1Cov, idx2Cov+1);
			    (*covMat)(idx1+1,idx2) = 0;
			    (*covMat)(idx1+1,idx2+1) = 0;
			  }
			// 4) both fixed (impossible)
			else
			  {
			    (*covMat)(idx1,idx2) = minuit->GetCovarianceMatrixElement(idx1Cov, idx2Cov);
			    (*covMat)(idx1,idx2+1) = 0;
			    (*covMat)(idx1+1,idx2) = 0;
			    (*covMat)(idx1+1,idx2+1) = 0;
			  }
		      }
		}
	    
	    NMCevents = myL.MCeventsInBin();
	    outTree->Fill();
	  }
	  
	  for (int j = 0; j < minuit->GetNumberTotalParameters(); j++)
	    startingValues[j].value = minuit->GetParameter(j);
	  
	  gHist.getHist("hMClikelihood", "MC likelihood of result", nBins, lower, upper)
	    ->SetBinContent(iBin+1, myL.calc_mc_likelihood(vStartingValues));
	  
	  if (!ambiguous){
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
	  }
	  
	  if (myL.getNChannels() == 2)
	    {
	      hBR->SetBinContent(iBin+1, minuit->GetParameter(nParams - 1));
	      hBR->SetBinError(iBin+1, minuit->GetParError(nParams - 1));
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
	  

	  // Fill predict histograms if wanted
	  // (relatively time consuming)
	  if (predict)
	    {
	      myL.fillPredict(result, hthpre, hphpre);	      
	    }
	  
	  // Ambiguities only for negative reflectivities
	  // for the moment hard-wired for this wave set: P+, D+, S0, P0, P-, D0, D-
	  if (ambiguous && lastIdx == 14)
	    {
	      // calculate coefficients
	      vector<double> amplitudes;
	      for (int k = 0; k < 14; k++)
		amplitudes.push_back(values[k]);
	    
	      vector<complex<double> > input;
	      waves2coeff(amplitudes, input);

	      // find 4 complex roots
	      vector<complex<double> > roots(4);
	      zroots(input, roots, true);
	      
	      sortroots(roots);

	      for (int l=0; l<4; l++){
		stringstream ss;         // conversion int to string
		ss << l;
		string num = ss.str();

		gHist.getHist(("hRe"+num).c_str(), "Real part of root",
			      nBins, threshold, threshold + nBins*binWidth)
		  ->SetBinContent(iBin+1, real(roots[l]));
		gHist.getHist(("hIm"+num).c_str(), "Imaginary part of root",
			      nBins, threshold, threshold + nBins*binWidth)
		  ->SetBinContent(iBin+1, /*TMath::Abs(*/imag(roots[l]));
	      }

	      // calculate 8 ambiguous solutions by complex conjugating
	      vector<double> ambiguous[8];
	      for (int j=0; j<8; j++){

		// first, fill positive reflectivity
		ambiguous[j].push_back(positiveStart[0]);
		ambiguous[j].push_back(positiveStart[1]);
		ambiguous[j].push_back(positiveStart[2]);
		ambiguous[j].push_back(positiveStart[3]);

		vector<complex<double> > solution(roots);      
		// 1 3 5 7
		if (j&0x1) solution[3] = conj(solution[3]);
		// 2 3 6 7
		if (j&0x2) solution[2] = conj(solution[2]);
		// 4 5 6 7
		if (j&0x4) solution[1] = conj(solution[1]);
	      
		//      vector<double> waves;
		roots2waves(input[4], solution, ambiguous[j]);
		// for (int m=0; m<14; m++)
		// cout << values[m] << " " << ambiguous[j][m] << endl;
	      }

	      // now run fit again with 8 solutions as starting value
	      for (int n = 0; n < 8; n++){
		// only run and plot result for solution 3
		if (n==3){
		  
		  for (size_t j= 0; j < nParams - myL.getNChannels(); j++)
		    {
		      if (!startingValues[j].fixed /*|| j < 4*/)
			{
			  fitamb->SetParameter(j, startingValues[j].name.c_str(),
					       ambiguous[n][j], 10/*ambiguous[n][j]*0.01*/, 0, 0);
			}
		      else
			{
			  fitamb->SetParameter(j, startingValues[j].name.c_str(),
					       ambiguous[n][j], 1, 0, 0);
			  fitamb->FixParameter(j);
			}
		    }
		  for (size_t j = nParams - myL.getNChannels(); j < nParams; j++)
		    {
		      fitamb->SetParameter(j, startingValues[j].name.c_str(),
					   vStartingValues[j], .1, 0, 1);
		      if (startingValues[j].fixed)
			fitamb->FixParameter(j);
		    }
		  
		  // Run minimizer.
		  fitamb->CreateMinimizer();
		  int amret = fitamb->Minimize();

		  // if fit did not converge, use new random starting values until it does
		  while ( amret != 0 ){
		    ambiguous[n][0] = gRandom->Uniform(10);
		    ambiguous[n][1] = gRandom->Uniform(10);
		    ambiguous[n][2] = gRandom->Uniform(10);
		    ambiguous[n][3] = gRandom->Uniform(10);

		    for (size_t j= 0; j < nParams - myL.getNChannels(); j++)
		      {
			if (!startingValues[j].fixed /*|| j < 4*/)
			  {
			    fitamb->SetParameter(j, startingValues[j].name.c_str(),
						 ambiguous[n][j], 10/*ambiguous[n][j]*0.01*/, 0, 0);
			  }
			else
			  {
			    fitamb->SetParameter(j, startingValues[j].name.c_str(),
						 ambiguous[n][j], 1, 0, 0);
			    fitamb->FixParameter(j);
			  }
		      }
		    for (size_t j = nParams - myL.getNChannels(); j < nParams; j++)
		      {
			fitamb->SetParameter(j, startingValues[j].name.c_str(),
					     vStartingValues[j], .1, 0, 1);
			if (startingValues[j].fixed)
			  fitamb->FixParameter(j);
		      }
		    amret = fitamb->Minimize();
		  }
		  //		  cout << "amret = " << amret << " after " << sw.CpuTime() << " s." << endl;
		  
		  if (amret == 0)
		    {
		      //	  vector<double> vStartingValue(nParams);
		      
		      for (int i = 0; i < 4; i++)
			{
			  positiveStart[i] = fitamb->GetParameter(i);
			}
		      
		      for (size_t i = 0; i < ws.size(); i++)
			for (size_t j = 0; j < ws[i].waves.size(); j++)
			  {
			    size_t idx = ws[i].waves[j].getIndex();
			    values[idx] = fitamb->GetParameter(idx);
			    values[idx+1] = fitamb->GetParameter(idx+1);
			  }
		      
		      for (size_t i1 = 0; i1 < ws.size(); i1++)
			for (size_t j1 = 0; j1 < ws[i1].waves.size(); j1++)
			  {
			    const wave& w1 = ws[i1].waves[j1];
			    size_t idx1 = w1.getIndex();
			    size_t idx1Cov = w1.idxInCovariance(fitamb);
			    
			    for (size_t i2 = 0; i2 < ws.size(); i2++)
			      for (size_t j2 = 0; j2 < ws[i2].waves.size(); j2++)
				{
				  const wave& w2 = ws[i2].waves[j2];
				  size_t idx2 = w2.getIndex();
				  size_t idx2Cov = w2.idxInCovariance(fitamb);
				  
				  // Four possibilities:
				  // 1) both free -> simply copy values to covMat
				  if (!w1.phaseLocked && !w2.phaseLocked)
				    {
				      (*covMat)(idx1,idx2) = fitamb->GetCovarianceMatrixElement(idx1Cov, idx2Cov);
				      (*covMat)(idx1,idx2+1) = fitamb->GetCovarianceMatrixElement(idx1Cov, idx2Cov+1);
				      (*covMat)(idx1+1,idx2) = fitamb->GetCovarianceMatrixElement(idx1Cov+1, idx2Cov);
				      (*covMat)(idx1+1,idx2+1) = fitamb->GetCovarianceMatrixElement(idx1Cov+1, idx2Cov+1);
				    }
				  // 2) first free, second fixed -> covariances with Im(w1) are zero
				  else if (!w1.phaseLocked && w2.phaseLocked)
				    {
				      (*covMat)(idx1,idx2) = fitamb->GetCovarianceMatrixElement(idx1Cov, idx2Cov);
				      (*covMat)(idx1,idx2+1) = 0;
				      (*covMat)(idx1+1,idx2) = fitamb->GetCovarianceMatrixElement(idx1Cov+1, idx2Cov);
				      (*covMat)(idx1+1,idx2+1) = 0;
				    }
				  // 3) vice versa
				  else if (w1.phaseLocked && !w2.phaseLocked)
				    {
				      (*covMat)(idx1,idx2) = fitamb->GetCovarianceMatrixElement(idx1Cov, idx2Cov);
				      (*covMat)(idx1,idx2+1) = fitamb->GetCovarianceMatrixElement(idx1Cov, idx2Cov+1);
				      (*covMat)(idx1+1,idx2) = 0;
				      (*covMat)(idx1+1,idx2+1) = 0;
				    }
				  // 4) both fixed (impossible)
				  else
				    {
				      (*covMat)(idx1,idx2) = fitamb->GetCovarianceMatrixElement(idx1Cov, idx2Cov);
				      (*covMat)(idx1,idx2+1) = 0;
				      (*covMat)(idx1+1,idx2) = 0;
				      (*covMat)(idx1+1,idx2+1) = 0;
				    }
				}
			  }
		    
		      NMCevents = myL.MCeventsInBin();
		      outTree->Fill();
		      
		      /*		      for (size_t iCoherent = 0; iCoherent < ws.size(); iCoherent++)
			{
			vector<wave>& waves = ws[iCoherent].getWaves();
			for (size_t iWave1 = 0; iWave1 < waves.size(); iWave1++)
			{
			complex<double> a(minuit->GetParameter(waves[iWave1].getIndex()),
			minuit->GetParameter(waves[iWave1].getIndex() + 1));
			
			gHist.getHist(("hAmbi"+waves[iWave1].name+"sol"+amb).c_str(),
			("fit results for this wave, solution "+amb).c_str(),
			nBins, threshold, threshold + nBins*binWidth)
			->SetBinContent(iBin+1, norm(a));
			}
			}*/
		      
		      for (size_t iCoherent = 0; iCoherent < ws.size(); iCoherent++)  
			{
			  vector<wave>& waves = ws[iCoherent].getWaves();
			  for (size_t iWave1 = 0; iWave1 < waves.size(); iWave1++)
			    {
			      waves[iWave1].fillHistIntensity(iBin, fitamb);
			      if (iWave1 != waves.size()-1)
				{
				  for (size_t iWave2 = iWave1 + 1; iWave2 < waves.size(); iWave2++)
				    waves[iWave1].fillHistPhase(iBin, waves[iWave2], fitamb);
				}
			    }  
			  
			}
		    }  
		}
	      }
	    }	      
	}
      else failBins++;
    }
  
  fulltime.Stop();
  cout << "Took " << fulltime.CpuTime() << " s CPU time, " << fulltime.RealTime() << " s wall time." << endl;
  cout << "In " << failBins << " bins TFitterMinuit::Minimization DID not converge !" << endl;
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
