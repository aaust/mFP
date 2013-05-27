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
#include "TRandom1.h"

#include "bw.h"
#include "fitInfo.h"
#include "gHist.h"
#include "fitModel.h"

const bool flipImag = true;

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
  if (m < 1.05 || m > 2.45)
    return 0;

  //if (iBin == 12)// || iBin<=4)
  //    return 0;

  model->evaluateAt(m, x);
  complex<double> Swave = model->valueForWave("S0");
  complex<double> Dwave = model->valueForWave("D0");

  TVectorD eps(3);
  eps[0] = values[2] - real(Swave);
  assert(values[3] == 0);

  eps[1] = values[4] - real(Dwave);
  eps[2] = fabs(values[5]) - imag(Dwave);
  
  gHist.Fill("hSwaveEvolution", "evolution of Swave",
	     info->getNbins(), info->getLower(), info->getLower() + (info->getNbins()+1) * info->getBinWidth(), 1000, 500, -500,
	     m, real(Swave));
  gHist.Fill("hDwaveReEvolution", "evolution of Re Dwave",
	     info->getNbins(), info->getLower(), info->getUpper(), 1000, 500, -500,
	     m, real(Dwave));
  gHist.Fill("hDwaveImEvolution", "evolution of Im Dwave",
	     info->getNbins(), info->getLower(), info->getUpper(), 1000, 500, -500,
	     m, (flipImag ? -1 : 1)*imag(Dwave));

  TMatrixDSym cov(6);
  int colSkipped = 0, rowSkipped = 0;
  for (int i = 0; i < 6; i++)
    {
      if ( i==1 || i==2 )
	rowSkipped++;
      colSkipped = 0;
      for (int j = 0; j < 6; j++)
	{
	  if ( j==1 || j==2 )
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
  
  TMatrixDSym G(3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      // i am looking only at S0 (real) and D0 (real, imaginary)
      G(i,j) = weight(i+1,j+1);

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

  TRandom1* gRandom = new TRandom1();
  gRandom->SetSeed2(0); // truely random!!!

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
	  if (j==5) vhPar[j]->SetBinContent(i, fabs(values[j]));
	  else vhPar[j]->SetBinContent(i, values[j]);
	  vhPar[j]->SetBinError(i, sqrt((*covMat)[j][j]));
	}
      cout << endl;
    }

  chiSquare fitFunc(tree);

  TFitterMinuit* minuit = new TFitterMinuit();
  minuit->SetMinuitFCN(&fitFunc);

  double vals[] = //{ 50, 50 , 1.275, 0.185, 50, 50, 1.515, 0.07, 50, 50, 2.13, .27, 200, 0, 1, .065, 50, 50, 1.5, 0.104, 0, 100, 1.730, 0.100, 50, 50, 1.37, .3, 50, 50, 0.5, 0.5, -0.5, 50, 50, 0.5, 0.5, -0.5} //WA102 
  //{ 50, 50, 1.275, 0.185, 50, 50, 1.32, .107, 50, 50, 1.515, 0.07, 200, 0, 1, .065, 100, 0, 1.5, 0.104, 100, 0, 1.730, 0.100, 50, 50, 0.5, -0.5, 0.5, 50, 50, 0.5, -0.5, 0.5} // f2/a2
  //{ 0, -80, 1.275, 0.185, 0, -80, 1.525, 0.073, 0, -50, 2.3, .150, /*0, 0, 0, 0, 0,*/ 200, 0, 1, .100, 100, 0, 1.505, 0.109, 100, 0, 1.720, 0.135/*, 0, 0, 0, 0, 0*/} //PDG
  //{-79, 30, 1.275, 0.185, -44, -50, 1.515, 0.07, 14, -31, 2.13, 0.27, -130, -539, 1.061, 0.200, -26, 90, 1.5, 0.104, -28, 49, 1.66, 0.196, -159, 395, 1.28, 0.358} //fit results without bg, fixed mass and width of f01500 and D
  //  {-79, 35, 1.275, 0.185, -44, -46, 1.515, 0.07, 12, -31, 2.13, 0.2, -133, -537, 1.062, 0.200, 22, 106, 1.486, 0.139, -12, 43, 1.667, 0.131, -182, 400, 1.290, 0.365, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}  //fit results without bg, fixed mass and width of D
  //{-138, -22, 1.275, 0.185, 32, -24, 1.515, 0.07, 40, 27, 2.5, 6.6, -50, -458, 1.09, 0.358, -62, -4, 1.486, 0.108, 8, -41, 1.712, 0.197, -182, 500, 1.290, 0.18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}  //fit results without bg, fixed mass and width of f and f'
    //{40, 27, 1.275, 0.185, 38, -95, 1.515, 0.07, 239, -249, 1.09, 0.358, -333, -14, 1.486, 0.108, 34, -24, 1.712, 0.197, 26, -44, 1.290, 0.18, 0.15, 6, -2, .5, 1, 0, 0, 0, 0, 0}  //fit results without bg, no f'', fixed all mass and width
    //    {-57, -22, 1.275, 0.185, 75, -124, 1.515, 0.07, 4, 2, -2.5, 0.6,
  // 239, -249, 1.09, 0.358, -33.3, -14, 1.486, 0.108, 34, -24, 1.712, 0.197, 26, -44, 1.290, 0.18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}// D with BG
  //   {21, -8, 1.32, 0.132, 50, 2, 1.515, 0.07, -0.5, 3, 0.4, -0.1, 0.2,//0.4, -2.5, 0.6,
  //-147, -283, 1.09, 0.358,
  //7370, 0.05, 2.5, -0.5,
  //39, 30, 1.486, 0.108, 21, -15, 1.712, 0.197, 32, -22, 1.290, 0.18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}// D with BG
    //-19, -16, 1.32, 0.132, -138, 91, 1.515, 0.07, 37, 75, 0.4, -0.1, 0.2,//0.4, -2.5, 0.6,
     //-147, -283, 1.09, 0.358,
    //664, 0.07, 4.2, -0.7,
  //, 302, 1.52, 0.108, -16, -50, 1.72, 0.100, 32, -22, 1.290, 0.18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}// D with BG
    //    {-60,  60, 1.29, 0.177, -4, -40, 1.515, 0.07, 2256, -50, 1.6, 1.9, //0.2,//0.4, -2.5, 0.6,
     //-147, -283, 1.09, 0.358,
    //     14925, 0.12, 3.2, //-0.7,
  //     -28, 91, 1.52, 0.108, -30, 30, 1.72, 0.100, 17, 384, 1.290, 0.18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}// D with BG
  //{-46,  56, 1.24, 0.177, -28, -14, 1.515, 0.12, 266, -187, -6.4, 4.2,//0.4, -2.5, 0.6,
  //-147, -283, 1.09, 0.358,
  //1027, 2.2, 3,
  //18, 57, 1.52, 0.108, -30, 12, 1.72, 0.100, 70, 377, 1.290, 0.18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}// shown on 20130109 with ex    
  //{-42,  52, 1.27, 0.216, -24, -18, 1.51, 0.1, 4, -3, -8.1, 5.1,//0.4, -2.5, 0.6,
  //1259, 2.2, 4,
  //37, 65, 1.48, 0.146, -34, 6, 1.72, 0.17, 69, 320, 1.29, 0.26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}// PS
    {-53,  52, 1.275, 0.185, -10, -30, 1.525, 0.073,
     11, -5, -6, 4,//0.4, -2.5, 0.6,
     //-147, -283, 1.09, 0.358,
     //450, 1.5,3.5,
     1259, 2.2, 2,
     13, 78, 1.48, 0.14, -27, 0, 1.72, 0.13, 43, 350, 1.29, 0.26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}// PS
  //{-60,  60, 1.29, 0.177, -4, -40, 1.515, 0.07, 708, -61, 2.2, 0.4, //0.2,//0.4, -2.5, 0.6,
  //-147, -283, 1.09, 0.358,
  //13807, 0.13, 3.2, //-0.7,
  //-28, 91, 1.52, 0.108, -30, 30, 1.72, 0.100, 17, 384, 1.290, 0.18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}//  with exp(-cm*m)
  //    {10,  10, 1.270, 0.185, 10, 10, 1.525, 0.073, 
  //708, -61, 1.2, 1.4,
  //13807, 0.13, 3.2,
  //10, 10, 1.505, 0.109, 10, 10, 1.72, 0.135, 10, 10, 1.300, 0.350, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}// start over with pdg
  ;
  
  minuit->SetParameter(0, "f_2 strength Re", vals[0], 1, 0, 0);
  //minuit->FixParameter(0);
  minuit->SetParameter(1, "f_2 strength Im", vals[1], 1, 0, 0);
  //minuit->FixParameter(1);
  minuit->SetParameter(2, "f_2 mass", vals[2], 1, 1.20, 1.35);
  //minuit->FixParameter(2);
  minuit->SetParameter(3, "f_2 width", vals[3], 1, 0.1, 0.5);
  //minuit->FixParameter(3);
  minuit->SetParameter(4, "f_2' strength Re", vals[4], 1, 0, 0);
  //minuit->FixParameter(4);
  minuit->SetParameter(5, "f_2' strength Im", vals[5], 1, 0, 0);
  //minuit->FixParameter(5);
  minuit->SetParameter(6, "f_2' mass", vals[6], 1, 1.45, 1.6);
  //minuit->FixParameter(6);
  minuit->SetParameter(7, "f_2' width", vals[7], 1, 0.05, 0.5);
  //minuit->FixParameter(7);
  
  minuit->SetParameter(8, "D BG strength Re", vals[8], 1, 0, 0);
  //minuit->FixParameter(8);
  minuit->SetParameter(9, "D BG strength Im", vals[9], 1, 0, 0);
  ///minuit->FixParameter(9);
  //minuit->SetParameter(10, "D BG const", vals[10], 1, 0, 0);
  //minuit->FixParameter(10);
  minuit->SetParameter(10, "D BG linear", vals[10], .1, 0, 0);
  //minuit->FixParameter(10);
  minuit->SetParameter(11, "D BG quadratic", vals[11], .1, 0, 0);
  //minuit->FixParameter(11);
  
  /*  
  minuit->SetParameter(12, "f_0(980) strength Re", vals[12], 1, 0, 0);
  //minuit->FixParameter(12);
  minuit->SetParameter(13, "f_0(980) strength Im", vals[13], 1, 0, 0);
  //minuit->FixParameter(13);
  minuit->SetParameter(14, "f_0(980) mass", vals[14], 0.01, 0, 0);
  minuit->FixParameter(14);
  minuit->SetParameter(15, "f_0(980) width", vals[15], 0.01, 0, 0);
  minuit->FixParameter(15);
  */

  //minuit->SetParameter(12, "S BG strength Re", vals[12], 1, 0, 0);
  //minuit->FixParameter(12);
  minuit->SetParameter(12, "S BG strength Re", vals[12], 1, 0, 0);
  //minuit->FixParameter(12);
  //minuit->SetParameter(12, "S BG const", vals[12], 1, 0, 0);
  //minuit->FixParameter(12);
  minuit->SetParameter(13, "S BG linear", vals[13], .1, 0, 0);
  //minuit->FixParameter(13);
  minuit->SetParameter(14, "S BG quadratic", vals[14], .1, 0, 0);
  //minuit->FixParameter(14);

  minuit->SetParameter(15, "f_0(1500) strength Re", vals[15], 1, 0, 0);
  //minuit->FixParameter(15);
  minuit->SetParameter(16, "f_0(1500) strength Im", vals[16], 1, 0, 0);
  //minuit->FixParameter(16);
  minuit->SetParameter(17, "f_0(1500) mass", vals[17], 1, 1.45, 1.55);
  //minuit->FixParameter(17);
  minuit->SetParameter(18, "f_0(1500) width", vals[18], 1, 0.05, 0.15);
  //minuit->FixParameter(18);
  minuit->SetParameter(19, "f_0(1710) strength Re", vals[19], 1, 0, 0);
  //minuit->FixParameter(19);
  minuit->SetParameter(20, "f_0(1710) strength Im", vals[20], 1, 0, 0);
  //minuit->FixParameter(20);
  minuit->SetParameter(21, "f_0(1710) mass", vals[21], 1, 1.65, 1.75);
  //minuit->FixParameter(21);
  minuit->SetParameter(22, "f_0(1710) width", vals[22], 1, 0.05, 0.2);
  //minuit->FixParameter(22);
  
  minuit->SetParameter(23, "f_0(1370) strength Re", vals[23], 1, 0, 0);
  //minuit->FixParameter(23);
  minuit->SetParameter(24, "f_0(1370) strength Im", vals[24], 1, 0, 0);
  //minuit->FixParameter(24);
  minuit->SetParameter(25, "f_0(1370) mass", vals[25], 1, 1.25, 1.40);
  //minuit->FixParameter(25);
  minuit->SetParameter(26, "f_0(1370) width", vals[26], 1, 0.05, 0.40);
  //minuit->FixParameter(26);
  
  /*
  minuit->SetParameter(24, "f_0(1370) strength Re", vals[24], 1, 0, 0);
  //minuit->FixParameter(24);
  minuit->SetParameter(25, "f_0(1370) strength Im", vals[25], 1, 0, 0);
  //minuit->FixParameter(25);
  minuit->SetParameter(26, "f_0(1370) mass", vals[26], 0.01, 0, 0);
  minuit->FixParameter(26);
  minuit->SetParameter(27, "f_0(1370) width", vals[27], 0.01, 0, 0);
  minuit->FixParameter(27);
  */
  
  
  
  /*
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
      fitModel* model = fitModel::getFitModelForName(info->getModelName());

      vector<double> x;
      for (int i = 0; i < minuit->GetNumberTotalParameters(); i++)
	x.push_back(minuit->GetParameter(i));

      for (size_t i = 0; i < info->getNbins(); i++)
	{
	  double m = info->getThreshold() + info->getBinWidth()*(i + 0.5);

	  //if (m < 2*mK || m > 2.45)
	  //continue;

	  model->evaluateAt(m, x);
	  //complex<double> phaseD = model->valueForWave("phaseD");
	  complex<double> Dwave = model->valueForWave("D0");
	  complex<double> DwaveBG = model->valueForWave("D0BG");
	  complex<double> Swave = model->valueForWave("S0");
	  complex<double> SwaveBG = model->valueForWave("S0BG");

	  complex<double> f2 = model->valueForWave("f2");
	  complex<double> f2prime = model->valueForWave("f2prime");

	  complex<double> f01370 = model->valueForWave("f01370");
	  complex<double> f01500 = model->valueForWave("f01500");
	  complex<double> f01710 = model->valueForWave("f01710");

	  //gHist.getHist("hPhaseD", "#phi(D)", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,arg(phaseD));
	  gHist.getHist("hSwaveRe", "Re(Swave)", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,real(Swave));
	  gHist.getHist("hSwaveIm", "Im(Swave)", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,imag(Swave));
	  gHist.getHist("hDwaveRe", "Re(Dwave)", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,real(Dwave));
	  gHist.getHist("hDwaveIm", "Im(Dwave)", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,imag(Dwave));
	  //ist.getHist("hDwaveBG", "DwaveBG", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,abs(DwaveBG));
	  /*	  gHist.getHist("hPwaveRe", "Re Pwave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,real(Pwave));
		  gHist.getHist("hPwaveIm", "Im Pwave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,(flipImag ? -1 : 1) * imag(Pwave));
	  gHist.getHist("hGwaveRe", "Re Gwave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,real(Gwave));
	  gHist.getHist("hGwaveIm", "Im Gwave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,(flipImag ? -1 : 1) * imag(Gwave));*/

	  gHist.getHist("hSwaveInt", "int SWave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,norm(Swave));
	  gHist.getHist("hSwaveBG", "int SWave BG", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,norm(SwaveBG));
	  
	  gHist.getHist("hDwaveInt", "int Dwave", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,norm(Dwave));
	  gHist.getHist("hDwaveBG", "int DwaveBG", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i,norm(DwaveBG));

	  gHist.getHist("hPhaseSD", "phase D - S", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i, arg(Swave/Dwave));
	  
	  gHist.getHist("hf2", "f2(1270)", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i, norm(f2));
	  gHist.getHist("hf2prime", "f2'(1525)", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i, norm(f2prime));

	  gHist.getHist("hf01370", "f0 1370", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i, norm(f01370));
	  gHist.getHist("hf01500", "f01500", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i, norm(f01500));
	  gHist.getHist("hf01710", "f01710", info->getNbins(), info->getLower(), info->getUpper())->SetBinContent(i, norm(f01710));
	}
    }

  fout->Write();
  delete fout;
  delete f;
  return 0;
}
