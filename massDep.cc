#include <vector>
#include <complex>
using namespace std;

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TFitterMinuit.h"
#include "TMatrixDSym.h"

#include "bw.h"
#include "fitInfo.h"


// Abstract interface.
// This is supposed to contain phase-space factors.
class waveParametrization {
public:
  virtual complex<double> value(double m, const double* par) = 0;
  virtual size_t getNpar() const = 0;
};


class simpleBreitWigner : public waveParametrization {
public:
  simpleBreitWigner(double m1_, double m2_, int J_)
    : m1(m1_), m2(m2_), J(J_) {};
  virtual complex<double> value(double m, const double* par);
  virtual size_t getNpar() const { return 2; }
private:
  double m1, m2;
  int J;
};

complex<double>
simpleBreitWigner::value(double m, const double* par)
{
  return breakupMomentum(m*m, m1, m2)*BW(m*m, m1, m2, J, par[0], par[1]);
}

class chiSquare : public ROOT::Minuit2::FCNBase {
public:
  chiSquare(const fitInfo* info_) : info(info_) {}

  double Up() const { return 1.; }
  double operator()(const vector<double>& x) const;

  void addBin(double massLow, double massHigh,
	      const vector<double>& values,
	      const TMatrixDSym& covMat)
  {
    this->massLow.push_back(massLow);
    this->massHigh.push_back(massHigh);
    this->values.push_back(vector<double>(&values[0],&values[info->getNparams()]));
    this->covMats.push_back(TMatrixDSym(covMat));
  }
private:
  const fitInfo* info;

  vector<double> massLow;
  vector<double> massHigh;
  vector<vector<double> > values;
  vector<TMatrixDSym> covMats;
};

double
chiSquare::operator()(const vector<double>& x) const 
{
  return 0;
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

  double *values = new double[info->getNparams()];
  TMatrixDSym *covMat = 0;
  double massLow, massHigh;
  tree->SetBranchAddress("massLow", &massLow);
  tree->SetBranchAddress("massHigh", &massHigh);
  tree->SetBranchAddress("values", values);
  tree->SetBranchAddress("covMat", &covMat);
  assert((Long_t)info->getNbins() == tree->GetEntries());
  for (Long_t i = 0; i < (Long_t)info->getNbins(); i++)
    {
      tree->GetEntry(i);
      cout << "bin " << i << " m = [" << massLow << "," << massHigh << "]" << endl;
      for (size_t j = 0; j < info->getNparams(); j++)
	cout << values[j] << " ";
      cout << endl;
    }

  chiSquare fitFunc(info);

  delete f;
  return 0;
}
