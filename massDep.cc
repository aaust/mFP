#include <vector>
#include <complex>
using namespace std;

#include "TObject.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TFitterMinuit.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"

#include "bw.h"
#include "fitInfo.h"
#include "gHist.h"


// Abstract interface.
// "value" is supposed to contain phase-space factors.
class waveParametrization {
public:
  virtual complex<double> value(double m, const double* par) const = 0;
  virtual size_t getNpar() const = 0;
};


class simpleBreitWigner : public waveParametrization {
public:
  simpleBreitWigner(double m1_, double m2_, int J_)
    : m1(m1_), m2(m2_), J(J_) {};
  virtual complex<double> value(double m, const double* par) const;
  virtual size_t getNpar() const { return 2; }
private:
  double m1, m2;
  int J;
};

complex<double>
simpleBreitWigner::value(double m, const double* par) const
{
  return breakupMomentum(m*m, m1, m2)*BW(m*m, m1, m2, J, par[0], par[1]);
}

class BreitWignerSum : public waveParametrization {
public:
   BreitWignerSum(double m1_, double m2_, int J_)
    : m1(m1_), m2(m2_), J(J_) {};
  virtual complex<double> value(double m, const double* par) const;
  virtual size_t getNpar() const { return 6; }
private:
  double m1, m2;
  int J;
};

complex<double>
BreitWignerSum::value(double m, const double* par) const
{
  return breakupMomentum(m*m, m1, m2)*(BW(m*m, m1, m2, J, par[0], par[1])
				       + (complex<double>(par[4], par[5])
					  *BW(m*m, m1, m2, J, par[2], par[3])));
}


class chiSquare : public ROOT::Minuit2::FCNBase {
public:
  chiSquare(TTree* tree_);

  double Up() const { return 1.; }
  double operator()(const vector<double>& x) const;

private:
  TTree *tree;
  const fitInfo *info;
  double valueInBin(Long_t i, const vector<double>& x) const;

  mutable double massLow;
  mutable double massHigh;
  mutable double* values;
  mutable TMatrixDSym* covMat;
};

chiSquare::chiSquare(TTree* tree_) : tree(tree_)
{
  info = 0;
  TIter next(const_cast<TTree*>(tree)->GetUserInfo());
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
      abort();
    }

  values = new double[info->getNvars()];
  covMat = 0;
  tree->SetBranchAddress("massLow", &massLow);
  tree->SetBranchAddress("massHigh", &massHigh);
  tree->SetBranchAddress("values", values);
  tree->SetBranchAddress("covMat", &covMat);
}

double
chiSquare::valueInBin(Long_t i, const vector<double>& x) const
{
  tree->GetEntry(i);

  // Let's do this manually for the time being:
  // x[0] ... x[6]  : params of D-wave
  // x[7] ... x[10]  : params of P-wave
  // x[11] ... x[18] : params of G-wave

  double m = .5*(massLow + massHigh);
  double m1 = mEtaP;
  double m2 = mPi;
  double phaseSpace = pow(sqrt(breakupMomentum(m*m, m1, m2)),2);

  const double *par = &x[1];
  complex<double> Dwave = 0;
  if (par[0] > 0 && par[1] > 0 && par[2] > par[0] && par[3] > 0)
    Dwave = x[0]*phaseSpace*(BW_a2_pietap_coupled(m*m) //BW(m*m, m1, m2, 2, par[0], par[1])
			     + (complex<double>(par[4], par[5])
				*BW(m*m, m1, m2, 2, par[2], par[3])));
  complex<double> phaseD = 1;
  if (abs(Dwave) != 0)
    phaseD = Dwave / abs(Dwave);
  Dwave /= phaseD;

  par = &x[9];
  complex<double> Pwave;
  if (par[0] > 0 && par[1] > 0)
    Pwave = complex<double>(x[7],x[8])*phaseSpace*BW(m*m, m1, m2, 1, par[0], par[1]);
  Pwave /= phaseD;

  par = &x[13];
  complex<double> Gwave;
  if (par[0] > 0 && par[1] > 0 && par[2] > par[0] && par[3] > 0)
    Gwave = complex<double>(x[11],x[12])*phaseSpace*(BW(m*m, m1, m2, 4, par[0], par[1])
						     + (complex<double>(par[4], par[5])
							*BW(m*m, m1, m2, 4, par[2], par[3])));

  //cout << Dwave << " " << Pwave << " " << Gwave << endl;

  TVectorD eps(5);
  eps[0] = values[0] - real(Dwave);
  eps[1] = values[2] - real(Pwave);
  eps[2] = values[3] - imag(Pwave);
  if (m > 1.7)
    {
      eps[3] = values[4] - real(Gwave);
      eps[4] = values[5] - imag(Gwave);
    }
  else
    eps[3] = eps[4] = 0;

  TMatrixDSym cov(5);
  cov(0,0) = (*covMat)(0,0);
  for (int i = 1; i < 5; i++)
    cov(0,i) = cov(i,0) = (*covMat)(0,i+1);
  for (int i = 1; i < 5; i++)
    for (int j = 1; j < 5; j++)
      cov(i,j) = cov(j,i) = (*covMat)(i+1, j+1);

  return cov.Invert().Similarity(eps);
}


double
chiSquare::operator()(const vector<double>& x) const 
{
  double value = 0;
  for (Long_t i = 3; i < 30 /*(Long_t)info->getNbins()*/; i++)
    {
      value += valueInBin(i, x);
    }
  cout << value << endl;
#if 1
  cout << " x = {";
  for (size_t i = 0; i < x.size(); i++)
    cout << x[i] << " ";
  cout << " } " << endl;
#endif
  return value;
}


int main(int argc, char **argv)
{
  if (argc != 2)
    {
      cerr << "Please give input file on command line." << endl;
      return 1;
    }
  TFile *f = TFile::Open(argv[1], "READ");
  if (!f)
    {
      cerr << "Can't open input file '" << argv[1] << "'" << endl;
      return 1;
    }
  TTree *tree;
  f->GetObject("tFitResults", tree);
  if (!tree)
    {
      cerr << "Can't find input tree 'tFitResults'" << endl;
      return 1;
    }

  fitInfo *info = 0;
  TIter next(tree->GetUserInfo());
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

  TFile *fout = TFile::Open("outMassDep.root", "RECREATE");
  fout->cd();

  double *values = new double[info->getNvars()];
  TMatrixDSym *covMat = 0;
  double massLow, massHigh;
  tree->SetBranchAddress("massLow", &massLow);
  tree->SetBranchAddress("massHigh", &massHigh);
  tree->SetBranchAddress("values", values);
  tree->SetBranchAddress("covMat", &covMat);
  assert((Long_t)info->getNbins() == tree->GetEntries());
  vector<TH1*> vhPar;
  for (size_t j = 0; j < info->getNvars(); j++)
    {
      char name[99];
      char title[99];
      snprintf(name, 99, "h%zd", j);
      snprintf(title, 99, "param %zd (%s)", j, info->getParamNames()[j].c_str());
      vhPar.push_back(gHist.getHist(name, title, info->getNbins(), 0, 1));
    }
  for (Long_t i = 0; i < (Long_t)info->getNbins(); i++)
    {
      tree->GetEntry(i);
      cout << "bin " << i << " m = [" << massLow << "," << massHigh << "]" << endl;
      covMat->Print();
      for (size_t j = 0; j < info->getNvars(); j++)
	{
	  cout << values[j] << " ";
	  vhPar[j]->SetBinContent(i, values[j]);
	  vhPar[j]->SetBinError(i, sqrt((*covMat)[j][j]));
	}
      cout << endl;
    }

  chiSquare fitFunc(tree);

  TFitterMinuit* minuit = new TFitterMinuit();
  minuit->SetMinuitFCN(&fitFunc);

  minuit->SetParameter(0, "D strength", 57.4152, 1, 0, 10000);
  minuit->SetParameter(1, "a_2 mass", 1.3183, 0, 0, 0);
  minuit->FixParameter(1);
  minuit->SetParameter(2, "a_2 width", 0.107, 0, 0, 0);
  minuit->FixParameter(2);
  minuit->SetParameter(3, "a_2' mass", 1.7, 0.02, 1.4, 10);
  minuit->SetParameter(4, "a_2' width", .3, 0.01, 0, 1);
  minuit->SetParameter(5, "a_2' strength Re", 1., 0.1, 0, 0);
  minuit->SetParameter(6, "a_2' strength Im", 0, 0.1, 0, 0);

  minuit->SetParameter(7, "P strength Re", 100, 1, 0, 0);
  minuit->SetParameter(8, "P strength Im", 0, 1, 0, 0);
  minuit->SetParameter(9, "P-wave mass", 1.399, 0.02, 0, 10);
  minuit->SetParameter(10, "P-wave width", .3, 0.05, 0, 10);

  minuit->SetParameter(11, "G strength Re", 10, 0.1, 0, 0);
  minuit->SetParameter(12, "G strength Im", 0, 0.1, 0, 0);
  minuit->SetParameter(13, "a_4 mass", 2.001, 0, 0, 0);
  minuit->FixParameter(13);
  minuit->SetParameter(14, "a_4 width", 0.235, 0, 0, 0);
  minuit->FixParameter(14);
  minuit->SetParameter(15, "a_4' mass", 2.5, 0.01, 2.1, 10);
  minuit->SetParameter(16, "a_4' width", .6, 0.01, 0, 10);
  minuit->SetParameter(17, "a_4' strength Re", 1., 0.1, 0, 0);
  minuit->SetParameter(18, "a_4' strength Im", 0, 0.1, 0, 0);
  
  TStopwatch sw;
  sw.Start();
  minuit->CreateMinimizer();
  int iret = minuit->Minimize();
  sw.Stop();
  cout << "iret = " << iret << " after " << sw.CpuTime() << " s." << endl;

  //if (!iret)
    {
      vector<double> x;
      for (int i = 0; i < minuit->GetNumberTotalParameters(); i++)
	x.push_back(minuit->GetParameter(i));

      for (int i = 0; i < 50; i++)
	{
	  double m = 1.1 + 2./50*i;

	  double m1 = mEtaP;
	  double m2 = mPi;
	  const double *par = &x[1];
	  complex<double> Dwave = x[0]*breakupMomentum(m*m, m1, m2)*(BW_a2_pietap_coupled(m*m) //BW(m*m, m1, m2, 2, par[0], par[1])
								     + (complex<double>(par[4], par[5])
									*BW(m*m, m1, m2, 2, par[2], par[3])));
	  complex<double> phaseD = 1;
	  if (abs(Dwave) != 0)
	    phaseD = Dwave / abs(Dwave);
	  Dwave /= phaseD;

	  par = &x[9];
	  complex<double> Pwave;
	  if (par[0] > 0 && par[1] > 0)
	    Pwave = complex<double>(x[7],x[8])*breakupMomentum(m*m, m1, m2)*BW(m*m, m1, m2, 1, par[0], par[1]);
	  Pwave /= phaseD;

	  par = &x[13];
	  complex<double> Gwave;
	  if (par[0] > 0 && par[1] > 0 && par[2] > par[0] && par[3] > 0)
	    Gwave = complex<double>(x[11],x[12])*breakupMomentum(m*m, m1, m2)*(BW(m*m, m1, m2, 4, par[0], par[1])
								       + (complex<double>(par[4], par[5])
									  *BW(m*m, m1, m2, 4, par[2], par[3])));

	  gHist.getHist("hDwaveRe", "Dwave", 50, 1.1, 3.1)->SetBinContent(i,real(Dwave));
	  gHist.getHist("hPwaveRe", "Re Pwave", 50, 1.1, 3.1)->SetBinContent(i,real(Pwave));
	  gHist.getHist("hPwaveIm", "Im Pwave", 50, 1.1, 3.1)->SetBinContent(i,imag(Pwave));
	  gHist.getHist("hGwaveRe", "Re Gwave", 50, 1.1, 3.1)->SetBinContent(i,real(Gwave));
	  gHist.getHist("hGwaveIm", "Im Gwave", 50, 1.1, 3.1)->SetBinContent(i,imag(Gwave));

	  gHist.getHist("hDwaveInt", "int DWave", 50, 1.1, 3.1)->SetBinContent(i,norm(Dwave));
	  gHist.getHist("hPwaveInt", "int Pwave", 50, 1.1, 3.1)->SetBinContent(i,norm(Pwave));
	  gHist.getHist("hGwaveInt", "int Gwave", 50, 1.1, 3.1)->SetBinContent(i,norm(Gwave));

	  gHist.getHist("hPhaseDP", "phase D - P", 50, 1.1, 3.1)->SetBinContent(i, arg(Dwave / Pwave));
	  gHist.getHist("hPhaseDG", "phase D - G", 50, 1.1, 3.1)->SetBinContent(i, arg(Dwave / Gwave));
	  gHist.getHist("hPhasePG", "phase P - G", 50, 1.1, 3.1)->SetBinContent(i, arg(Pwave / Gwave));
	}
    }

  fout->Write();
  delete fout;
  delete f;
  return 0;
}
