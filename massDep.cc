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

const bool flipImag = true;

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
  //cout << " m = " << m << endl;
  double m1 = mEtaP;
  double m2 = mPi;
  double phaseSpace = 1; //sqrt(m); //1; //sqrt(breakupMomentum(m*m, m1, m2)); // / m; // sqrt(m);

  const double *par = &x[1];
  complex<double> Dwave = 0;
  if (par[0] > 0 && par[1] > 0 && par[2] > par[0] && par[3] > 0)
    Dwave = x[0]*phaseSpace*(BW(m*m, m1, m2, 2, par[0], par[1]) //BW(m*m, m1, m2, 2, par[0], par[1])
			     //+ (complex<double>(par[4], par[5])
			     //	*BW(m*m, m1, m2, 2, par[2], par[3]))
			     );
  //cout << BW_a2_pieta_coupled(m*m) << BW(m*m, m1, m2, 2, par[2], par[3]) << endl;
  complex<double> phaseD = 1;
  if (abs(Dwave) != 0)
    phaseD = Dwave / abs(Dwave);
  Dwave /= phaseD;

  par = &x[9];
  complex<double> Pwave = 0;
  if (par[0] > 0 && par[1] > 0 && par[0] > mEtaP + mPi)
    Pwave = complex<double>(x[7],x[8])*phaseSpace*BW(m*m, m1, m2, 1, par[0], par[1]);
  //cout << "BW(" << m << "," << m1 << "," << m2 << ",1,1," << par[0] << "," << par[1] << ") = " << Pwave << endl;
  Pwave /= phaseD;

  par = &x[13];
  complex<double> Gwave = 0;
  if (par[0] > 0 && par[1] > 0 && par[2] > par[0] && par[3] > 0)
    Gwave = complex<double>(x[11],x[12])*phaseSpace*(BW(m*m, m1, m2, 4, par[0], par[1])
						     //+ (complex<double>(par[4], par[5])
						     //	*BW(m*m, m1, m2, 4, par[2], par[3]))
						     );
  Gwave /= phaseD;

  //cout << "mDPG = " << m << " " << Dwave << " " << Pwave << " " << Gwave << endl;
  TVectorD eps(5);
  eps[0] = values[0] - real(Dwave);
  assert(values[1] == 0);
  eps[1] = values[2] - real(Pwave);
  eps[2] = (flipImag ? -1 : 1) * values[3] - imag(Pwave);
  if (m > 1.4)
    {
      eps[3] = values[4] - real(Gwave);
      eps[4] =  (flipImag ? -1 : 1) * values[5] - imag(Gwave);
    }
  else
    eps[3] = eps[4] = 0;

  gHist.Fill("hDwaveEvolution", "evolution of Dwave",
	     info->getNbins(), info->getLower(), info->getUpper(), 1000, -500, 500,
	     m, real(Dwave));
  gHist.Fill("hPwaveReEvolution", "evolution of Re Pwave",
	     info->getNbins(), info->getLower(), info->getUpper(), 1000, -500, 500,
	     m, real(Pwave));
  gHist.Fill("hPwaveImEvolution", "evolution of Im Pwave",
	     info->getNbins(), info->getLower(), info->getUpper(), 1000, -500, 500,
	     m, imag(Pwave));
  if (m > 1.4)
    {
      gHist.Fill("hGwaveReEvolution", "evolution of Re Gwave",
		 info->getNbins(), info->getLower(), info->getUpper(), 1000, -500, 500,
		 m, real(Gwave));
      gHist.Fill("hGwaveImEvolution", "evolution of Im Gwave",
		 info->getNbins(), info->getLower(), info->getUpper(), 1000, -500, 500,
		 m, imag(Gwave));
    }

  TMatrixDSym cov(5);
  cov(0,0) = (*covMat)(0,0);
  for (int i = 1; i < 5; i++)
    cov(0,i) = cov(i,0) = (flipImag && (i == 3 || i == 5) ? -1 : 1) * (*covMat)(0,i+1);
  for (int i = 1; i < 5; i++)
    for (int j = 1; j < 5; j++)
      cov(i,j) = cov(j,i) =  (flipImag && ((i == 3 && j != 3)
					   || (i == 5 && j != 3)) ? -1 : 1) * (*covMat)(i+1, j+1);

  bool status;
  TMatrixDSym G(TDecompChol(cov).Invert(status));
  if (!status)
    {
      cout << "non-positive-definite covariance matrix in bin " << iBin << endl;
      cov.Print();
      return 0;
    }

  if (0 && iBin == 4)
    {
      //cov.Invert().Print();
      cout << values[0] << " - " << real(Dwave) << endl;
      cout << G.Similarity(eps) << endl;
    }

  if (G.Similarity(eps) < 10000)
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
      //if (i == 4 || i == 23)
      //continue;
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
    {160, 1.3183, 0.107, 2, .5, 0, 0, 200, 0, 1.6, 0.4, 3.5, 0, 2.001, 0.235, 5., 0.7, 0, 0 }
;

  minuit->SetParameter(0, "D strength", vals[0], .1, 0, 0);
  //minuit->FixParameter(0);
  minuit->SetParameter(1, "a_2 mass", vals[1], 0.01, 0, 0);
  //minuit->FixParameter(1);
  minuit->SetParameter(2, "a_2 width", vals[2], 0.01, 0.0, 0.);
  //minuit->FixParameter(2);
  minuit->SetParameter(3, "a_2' mass", vals[3], 0.05, 1.4, 3);
  minuit->FixParameter(3);
  minuit->SetParameter(4, "a_2' width", vals[4], 0.2, 0, 1.7);
  minuit->FixParameter(4);
  minuit->SetParameter(5, "a_2' strength Re", vals[5], 0.1, 0, 0);
  minuit->FixParameter(5);
  minuit->SetParameter(6, "a_2' strength Im", vals[6], 0.1, 0, 0);
  minuit->FixParameter(6);

  minuit->SetParameter(7, "P strength Re", vals[7], 1, 0, 0);
  //minuit->FixParameter(7);
  minuit->SetParameter(8, "P strength Im", vals[8], 1, 0, 0);
  //minuit->FixParameter(8);
  minuit->SetParameter(9, "P-wave mass", vals[9], 0.02, 0, 2);
  //minuit->FixParameter(9);
  minuit->SetParameter(10, "P-wave width", vals[10], 0.25, 0, 1.5);
  //minuit->FixParameter(10);

  minuit->SetParameter(11, "G strength Re", vals[11], 0.1, -5, 5);
  //minuit->FixParameter(11);
  minuit->SetParameter(12, "G strength Im", vals[12], 0.1, 0, 0);
  //minuit->FixParameter(12);
  minuit->SetParameter(13, "a_4 mass", vals[13], 0.01, 0, 10);
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

  /*
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
  */
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
      vector<double> x;
      for (int i = 0; i < minuit->GetNumberTotalParameters(); i++)
	x.push_back(minuit->GetParameter(i));

      for (size_t i = 0; i < info->getNbins(); i++)
	{
	  double m = info->getThreshold() + info->getBinWidth()*i;

	  double m1 = mEtaP;
	  double m2 = mPi;
	  double phaseSpace = 1; //sqrt(m); //breakupMomentum(m*m,m1,m2));// sqrt(m); //breakupMomentum(m*m, m1, m2)) / m;

	  const double *par = &x[1];
	  complex<double> Dwave = 0;
	  if (par[0] > 0 && par[1] > 0 && par[2] > par[0] && par[3] > 0)
	    Dwave = x[0]*phaseSpace*(BW(m*m, m1, m2, 2, par[0], par[1]) //BWcoupled(m*m, m1, m2, mPi, 0.77, 1, par[0], par[1], 0.2) //BW(m*m, m1, m2, 2, par[0], par[1])
				     //+ (complex<double>(par[4], par[5])
				     //	*BW(m*m, m1, m2, 2, par[2], par[3]))
				     );
	  //cout << Dwave<< BW_a2_pieta_coupled(m*m) << BW(m*m, m1, m2, 2, par[2], par[3]) << endl;
	  complex<double> phaseD = 1;
	  if (abs(Dwave) != 0)
	    phaseD = Dwave / abs(Dwave);
	  Dwave /= phaseD;

	  par = &x[9];
	  complex<double> Pwave = 0;
	  if (par[0] > 0 && par[1] > 0 && par[0] > mEta + mPi)
	    Pwave = complex<double>(x[7],x[8])*phaseSpace*BW(m*m, m1, m2, 1, par[0], par[1]);
	  //cout << "BW(" << m << "," << m1 << "," << m2 << ",1,1," << par[0] << "," << par[1] << ") = " << Pwave << endl;
	  gHist.getHist("hPhaseP", "#phi(P)", 50, 0.725, 0.725+50*0.05)->SetBinContent(i,arg(Pwave));
	  Pwave /= phaseD;

	  par = &x[13];
	  complex<double> Gwave = 0;
	  if (par[0] > 0 && par[1] > 0 && par[2] > par[0] && par[3] > 0)
	    Gwave = complex<double>(x[11],x[12])*phaseSpace*(BW(m*m, m1, m2, 4, par[0], par[1])
							     //+ (complex<double>(par[4], par[5])
							     //	*BW(m*m, m1, m2, 4, par[2], par[3]))
							     );
	  Gwave /= phaseD;

	  gHist.getHist("hPhaseD", "#phi(D)", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,arg(phaseD));
	  gHist.getHist("hDwaveRe", "Dwave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,real(Dwave));
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
