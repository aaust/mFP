#include <vector>
using namespace std;

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TFitterMinuit.h"
#include "TMatrixDSym.h"

#include "bw.h"
#include "fitInfo.h"

class chiSquare : public ROOT::Minuit2::FCNBase {
public:
  double Up() const { return 1.; }
  double operator()(const vector<double>& x) const;
};

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

  delete f;
  return 0;
}
