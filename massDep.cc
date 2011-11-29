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
#include "TDecompChol.h"

#include "bw.h"
#include "fitInfo.h"
#include "gHist.h"
#include "fitModel.h"

const bool flipImag = false;

class chiSquare : public ROOT::Minuit2::FCNBase {
public:
  chiSquare(TTree* tree_);

  double Up() const { return 1.; }
  double operator()(const vector<double>& x) const;

  ~chiSquare() { delete model; }

private:
  TTree *tree;
  const fitInfo *info;
  fitModel *model;
  double valueInBin(Long_t i, const vector<double>& x) const;

  mutable double massLow;
  mutable double massHigh;
  mutable double* values;
  mutable TMatrixDSym* covMat;
};

chiSquare::chiSquare(TTree* tree_) : tree(tree_)
{
  info = 0;
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
      abort();
    }

  values = new double[info->getNvars()];
  covMat = 0;
  tree->SetBranchAddress("massLow", &massLow);
  tree->SetBranchAddress("massHigh", &massHigh);
  tree->SetBranchAddress("values", values);
  tree->SetBranchAddress("covMat", &covMat);

  model = fitModel::getFitModelForName(info->getModelName());
}

double
chiSquare::valueInBin(Long_t iBin, const vector<double>& x) const
{
  tree->GetEntry(iBin);

  // Let's do this manually for the time being:
  // x[0] ... x[6]  : params of D-wave
  // x[7] ... x[10]  : params of P-wave
  // x[11] ... x[18] : params of G-wave

  double m = .5*(massLow + massHigh);
  if (m < 0.77+mPi || m > 2.7)
    return 0;

  model->evaluateAt(m, x);
  complex<double> Dwave = model->valueForWave("D+");
  complex<double> Pwave = model->valueForWave("P+");
  complex<double> Gwave = model->valueForWave("G+");

  TVectorD eps(5);
  eps[0] = values[0] - real(Dwave);
  assert(values[1] == 0);

  if (m < 2.)
    {
      eps[1] = values[2] - real(Pwave);
      eps[2] = (flipImag ? -1 : 1) * values[3] - imag(Pwave);
    }
  else
    eps[1] = eps[2] = 0;
  if (m > 1.5 && m < 2.7)
    {
      eps[3] = values[4] - real(Gwave);
      eps[4] =  (flipImag ? -1 : 1) * values[5] - imag(Gwave);
    }
  else
    eps[3] = eps[4] = 0;
  //eps[1] = eps[2] = eps[3] = eps[4] = 0;

  gHist.Fill("hDwaveEvolution", "evolution of Dwave",
	     info->getNbins(), info->getLower(), info->getLower() + (info->getNbins()+1) * info->getBinWidth(), 1000, 500, -500,
	     m, real(Dwave));
  gHist.Fill("hPwaveReEvolution", "evolution of Re Pwave",
	     info->getNbins(), info->getLower(), info->getUpper(), 1000, 500, -500,
	     m, real(Pwave));
  gHist.Fill("hPwaveImEvolution", "evolution of Im Pwave",
	     info->getNbins(), info->getLower(), info->getUpper(), 1000, 500, -500,
	     m, (flipImag ? -1 : 1)*imag(Pwave));
  if (1 || m > 1.4)
    {
      gHist.Fill("hGwaveReEvolution", "evolution of Re Gwave",
		 info->getNbins(), info->getLower(), info->getUpper(), 1000, 500, -500,
		 m, real(Gwave));
      gHist.Fill("hGwaveImEvolution", "evolution of Im Gwave",
		 info->getNbins(), info->getLower(), info->getUpper(), 1000, 500, -500,
		 m, (flipImag ? -1 : 1) * imag(Gwave));
    }

  TMatrixDSym cov(14);
  int colSkipped = 0, rowSkipped = 0;
  for (int i = 0; i < 14; i++)
    {
      if (i == 1 || i == 6)
	rowSkipped++;
      colSkipped = 0;
      for (int j = 0; j < 14; j++)
	{
	  if (j == 1 || j == 6)
	    colSkipped++;
	  cov(i,j) = (*covMat)(i + rowSkipped, j + colSkipped);
	}
    }
  bool status;
  TMatrixDSym weight(TDecompChol(cov).Invert(status));
  if (!status)
    {
      cout << "non-positive-definite covariance matrix in bin " << iBin << endl;
      abort();
    }

  TMatrixDSym G(5);
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++)
      G(i,j) = weight(i,j);

  if (G.Similarity(eps) > 10000)
    {
      cout << "Gigantic value " << G.Similarity(eps) << " in bin " << iBin << endl;
      G.Print();
      eps.Print();
    }
  if (G.Similarity(eps) < 0)
    {
      cout << "negative value " << G.Similarity(eps) << " in bin " << iBin << endl;
      G.Print();
      eps.Print();
    }
  return G.Similarity(eps);
}


double
chiSquare::operator()(const vector<double>& x) const 
{
  double value = 0;
  for (Long_t i = 0; i < (Long_t)info->getNbins(); i++)
    {
      value += valueInBin(i, x);
    }
  if (value < 0)
    {
      for (Long_t i = 0; i < 25 /*(Long_t)info->getNbins()*/; i++)
	{
	  if (i == 4 || i == 23)
	    continue;
	  cout << valueInBin(i, x) << " ";
	}
      cout << endl;
    }
#if 0
  cout << " x = {";
  for (size_t i = 0; i < x.size(); i++)
    cout << x[i] << ", ";
  cout << " } " << endl;
#endif
  gHist.Fill("hchi2", "chi2s obtained during fitting steps", 1000, 0, 200000, value);
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
  cout << info->getNbins() << " " << tree->GetEntries() << endl;
  //assert((Long_t)info->getNbins() == tree->GetEntries());
  vector<TH1*> vhPar;
  for (size_t j = 0; j < info->getNvars(); j++)
    {
      char name[99];
      char title[99];
      snprintf(name, 99, "h%zd", j);
      snprintf(title, 99, "param %zd (%s)", j, info->getParamNames()[j].c_str());
      vhPar.push_back(gHist.getHist(name, title, info->getNbins(), info->getLower(), info->getUpper()));
    }
  for (Long_t i = 0; i < (Long_t)info->getNbins()-1; i++)
    {
      tree->GetEntry(i);
      cout << "bin " << i << " m = [" << massLow << "," << massHigh << "]" << endl;
      for (size_t j = 0; j < info->getNvars(); j++)
	{
	  vhPar[j]->SetBinContent(i, values[j]);
	  vhPar[j]->SetBinError(i, sqrt((*covMat)[j][j]));
	}
      cout << endl;
    }

  chiSquare fitFunc(tree);

  TFitterMinuit* minuit = new TFitterMinuit();
  minuit->SetMinuitFCN(&fitFunc);

  double vals[] = //{12.5161, 1.3183, 0.107, 1.90331, 0.99997, 10.2365, 29.0293, -192.712, 237.324, 1.37779, 0.400319, -1.42882, -0.294888, 2.001, 0.235, 2.71232, 0.153474, -1.84574, -12.2434,  }
  //{77.728, 1.3183, 0.107, 6.82604, 0.470435, -51.0442, 23.7011, -297.761, -23.8848, 1.39258, 0.400842, -0.965659, -0.348415, 2.001, 0.235, 4.64518, 0.0361006, -84.1068, -172.716,  } 
  //{77.688, 1.3183, 0.107, 7.71038, 0.415327, -61.8941, 28.4093, -297.45, -24.1899, 1.39193, 0.400172, -0.996119, -0.372143, 2.001, 0.235, 8.30321, 0.0033461, -1047.66, -2202.12,  } 
  //{13.1807, 1.3183, 0.107, 9.98529, 0.999804, -257.594, 1157.81, -1152.6, -1647.65, 3.07546, 7.26251, -0.965659, -0.348415, 2.001, 0.235, 4.64518, 0.0361006, -84.1068, -242.008,  } 
    //  {13.7333, 1.3183, 0.107, 2.99991, 0.782959, -152.037, 59.0287, -19237.9, 8934.96, 1.99998, 1.1686e-06, -0.965659, -0.348415, 2.001, 0.235, 4.64518, 0.0361006, -84.1068, -242.008,  }
  //{91.2756, 1.28142, 0.0888551, 2.99993, 0.00757554, -17.8072, 9.42176, -400282, 609404, 1.89095, 5.39743e-10, 46.211, -8.99411, 2.001, 0.235, 9.98809, 0.0106126, -575.415, -105.118,  }
  //{144.564, 1.31293, 0.129749, 1.76617, 5.55112e-17, -17.8247, -1.58795e+06, -213947, 313149, 1.92468, 1.52465e-09, -108.293, -20.1702, 2.001, 0.235, 9.99909, 0.010058, -568.593, -167.354,  } 
  //{112.217, 1.29986, 0.0999167, 1.76678, 6.99046e-12, -39651.9, 58861.1, 165077, -134251, 1.1064, 3.88578e-16, -0.00484132, -0.337098, 2.001, 0.235, 9.23719, 0.935933, 374.254, -85.3632,  } 
  //{134.733, 1.31361, 0.119147, 1.79072, 0.00399844, 7700.92, 5441.68, 17827.8, 31456.4, 1.10929, 0.0231332, 0.474292, 4.70543, 2.001, 0.235, 9.38998, 7.48521e-05, 817.705, -434.581,  }
  //{137.706, 1.31412, 0.132691, 1.78024, 0.00239979, -0.0639616, 0.793151, 25.5136, 5.09123, 1.6503, 0.307946,  }
  //{13.7333, 1.3183, 0.107, 2.99991, 0.782959, 0, 0, -19237.9, 8934.96, 1.99998, 1.1686e-06, -0.965659, -0.348415, 2.001, 0.235, 4.64518, 0.0361006, 0, 0,  }
    //{234, 1.3183, 0.107, 1.9, .5, 0, 0, 100, -18, 1.4, 0.4, 121.445, 43.7584, 2.0253, 0.28722, 5, 0.7, 0, 0 } // 2.001, 0.235, 5., 0.7, 0, 0 }
    //    { 493.863, 1.313, 0.1057, 1.8, 0.6, 0, 0, 142.454, -94.109, 1.35, 0.3, 130.12059638, 24.3762, 1.99, 0.236181, 5, 0.4, 1.3, -1.5, }
    { 315.601, 1.31176, 0.106659, 1.8, 0.6, 0, 0, 95.8984, -54.8463, 1.38494, 0.43218, 76.5587, 35.9412, 2.07055, 0.643189, 5, 0.4, 1.3, -1.5, }
  ;

  minuit->SetParameter(0, "D strength", vals[0], .1, 0, 0);
  //minuit->FixParameter(0);
  minuit->SetParameter(1, "a_2 mass", vals[1], 0.01, 0, 0);
  minuit->FixParameter(1);
  minuit->SetParameter(2, "a_2 width", vals[2], 0.01, 0.0, 0.);
  minuit->FixParameter(2);
  minuit->SetParameter(3, "a_2' mass", vals[3], 0.05, 1.5, 2.2);
  minuit->FixParameter(3);
  minuit->SetParameter(4, "a_2' width", vals[4], 0.2, 0.3, 1.);
  minuit->FixParameter(4);
  minuit->SetParameter(5, "a_2' strength Re", vals[5], 0.1, 0, 0);
  minuit->FixParameter(5);
  minuit->SetParameter(6, "a_2' strength Im", vals[6], 0., 0, 0);
  minuit->FixParameter(6);

  minuit->SetParameter(7, "P strength Re", vals[7], 1, 0, 0);
  //minuit->FixParameter(7);
  minuit->SetParameter(8, "P strength Im", vals[8], 1, 0, 0);
  //minuit->FixParameter(8);
  minuit->SetParameter(9, "P-wave mass", vals[9], 0.02, 0, 2);
  //minuit->FixParameter(9);
  minuit->SetParameter(10, "P-wave width", vals[10], 0.25, 0, 1.5);
  //minuit->FixParameter(10);

  minuit->SetParameter(11, "G strength Re", vals[11], 0.1, 0, 0);
  //minuit->FixParameter(11);
  minuit->SetParameter(12, "G strength Im", vals[12], 0.1, 0, 0);
  //minuit->FixParameter(12);
  minuit->SetParameter(13, "a_4 mass", vals[13], 0.01, 0, 0);
  //minuit->FixParameter(13);
  minuit->SetParameter(14, "a_4 width", vals[14], 0.02, 0, 0);
  //minuit->FixParameter(14);
  minuit->SetParameter(15, "a_4' mass", vals[15], 0.01, 2.1, 10);
  minuit->FixParameter(15);
  minuit->SetParameter(16, "a_4' width", vals[16], 0.01, 0, 10);
  minuit->FixParameter(16);
  minuit->SetParameter(17, "a_4' strength Re", vals[17], 0.1, 0, 0);
  minuit->FixParameter(17);
  minuit->SetParameter(18, "a_4' strength Im", vals[18], 0.1, 0, 0);
  minuit->FixParameter(18);

  minuit->SetParameter(19, "D BG exp width", 3.246672960364, 0.1, 0, 10);
  minuit->FixParameter(19);
  minuit->SetParameter(20, "D BG const", 0.3924742977682, 0.1, 0, 0);
  minuit->FixParameter(20);
  minuit->SetParameter(21, "D BG linear", 15.89705694124, 0.1, 0, 0);
  //minuit->FixParameter(21);
  minuit->SetParameter(22, "D BG quadratic", -13.36962331843, 0.1, 0, 0);
  minuit->FixParameter(22);

  /*
  minuit->SetParameter(19, "D BG exp width", 0, 0.1, 0, 10);
  minuit->FixParameter(19);
  minuit->SetParameter(20, "D BG const", 0, 0.1, 0, 0);
  minuit->FixParameter(20);
  minuit->SetParameter(21, "D BG linear", 0, 0.1, 0, 0);
  minuit->FixParameter(21);
  minuit->SetParameter(22, "D BG quadratic", 0, 0.1, 0, 0);
  minuit->FixParameter(22);
  */
  /*
  minuit->SetParameter(23, "G BG exp width", 9.02584, 0.1, 0, 10);
  minuit->FixParameter(23);
  minuit->SetParameter(24, "G BG const", -148.13, 0.1, 0, 0);
  minuit->FixParameter(24);
  minuit->SetParameter(25, "G BG linear", -56.0833, 0.1, 0, 0);
  minuit->FixParameter(25);
  minuit->SetParameter(26, "G BG quadratic", 8.77022, 0.1, 0, 0);
  minuit->FixParameter(26);
  */
  minuit->SetParameter(23, "G BG exp width", 0, 0.1, 0, 10);
  minuit->FixParameter(23);
  minuit->SetParameter(24, "G BG const", 0, 0.1, 0, 0);
  minuit->FixParameter(24);
  minuit->SetParameter(25, "G BG linear", 0, 0.1, 0, 0);
  minuit->FixParameter(25);
  minuit->SetParameter(26, "G BG quadratic", 0, 0.1, 0, 0);
  minuit->FixParameter(26);

  TStopwatch sw;
  sw.Start();
  minuit->CreateMinimizer();
  int iret = minuit->Minimize();
  sw.Stop();
  cout << "iret = " << iret << " after " << sw.CpuTime() << " s." << endl;
  if (iret)
    minuit->PrintResults(1, 0.);

  //if (!iret)
    {
      fitModel* model = fitModel::getFitModelForName(info->getModelName());

      vector<double> x;
      for (int i = 0; i < minuit->GetNumberTotalParameters(); i++)
	x.push_back(minuit->GetParameter(i));

      for (size_t i = 0; i < info->getNbins(); i++)
	{
	  double m = info->getThreshold() + info->getBinWidth()*(i + 0.5);

	  if (m < 0.77+mPi)
	    continue;

	  model->evaluateAt(m, x);
	  complex<double> phaseD = model->valueForWave("phaseD");
	  complex<double> Dwave = model->valueForWave("D+");
	  complex<double> DwaveBG = model->valueForWave("D+BG");
	  complex<double> Pwave = model->valueForWave("P+");
	  complex<double> Gwave = model->valueForWave("G+");
	  if (flipImag)
	    {
	      Pwave = conj(Pwave);
	      Gwave = conj(Gwave);
	    }

	  gHist.getHist("hPhaseD", "#phi(D)", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,arg(phaseD));
	  gHist.getHist("hDwaveRe", "Dwave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,real(Dwave));
	  gHist.getHist("hDwaveBG", "DwaveBG", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,abs(DwaveBG));
	  gHist.getHist("hPwaveRe", "Re Pwave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,real(Pwave));
	  gHist.getHist("hPwaveIm", "Im Pwave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,(flipImag ? -1 : 1) * imag(Pwave));
	  gHist.getHist("hGwaveRe", "Re Gwave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,real(Gwave));
	  gHist.getHist("hGwaveIm", "Im Gwave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,(flipImag ? -1 : 1) * imag(Gwave));

	  gHist.getHist("hDwaveInt", "int DWave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,norm(Dwave));
	  gHist.getHist("hPwaveInt", "int Pwave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,norm(Pwave));
	  gHist.getHist("hGwaveInt", "int Gwave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,norm(Gwave));

	  gHist.getHist("hPhaseDP", "phase P - D", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i, arg(Pwave));
	  gHist.getHist("hPhaseDG", "phase G - D", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i, arg(Gwave));
	  gHist.getHist("hPhasePG", "phase P - G", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i, arg(Pwave / Gwave));
	}
    }

  fout->Write();
  delete fout;
  delete f;
  return 0;
}
