#include <complex>
#include <vector>
#include <map>

using namespace std;

#include "TRandom1.h"
#include "Math/SpecFunc.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"
#include "TStopwatch.h"

#define NRDEVENTS 1000
#define NMCEVENTS (4*NRDEVENTS)

struct wave {
  int l, m;
  std::complex<double> a;

  wave() { }
  wave(int ll, int mm) { l = ll; m = mm; a = 0; }
  wave(const wave& o) { l = o.l; m = o.m; a = o.a; }
};

struct coherent_waves {
  int reflectivity;
  int spinflip;
  vector<wave> waves;

  coherent_waves() {}
  coherent_waves(const coherent_waves& o) { reflectivity = o.reflectivity; spinflip = o.spinflip; waves = o.waves; }
};

typedef vector<coherent_waves> waveset;

struct event {
  double theta;
  double phi;

  event() { theta = phi = 0; }
  /*
  event(const event& other) { theta = other.theta; phi = other.phi; }
  */
  event(double th, double ph) { theta = th; phi = ph; }
};

event RDevents[NRDEVENTS]; // not a reference because of Cint limitations
event MCevents[NMCEVENTS]; // idem

class likelihood : public ROOT::Minuit2::FCNBase {
  vector<coherent_waves> ws;
public:
  likelihood(waveset ws_) /*,
	     vector<event>& RDevents_,
	     vector<event>& MCevents_) */
    : ws(ws_) /*, RDevents(RDevents_), MCevents(MCevents_)*/ {}
  //double operator() (const std::vector<double>& x) const;
  double Up() const { return 0.5; }

private:
  double
  decay(int reflectivity, int l, int m, double theta, double phi) const;

  complex<double>
  coherentSum(const coherent_waves& waves, const event& e) const;

  double
  probabilityDensity(const event& e) const;

  double
  mcWeight(int reflectivity, int l1, int m1, int l2, int m2) const;

  double
  calc_likelihood() const;

public:
  double
  operator()(const vector<double>& x) const;

};


// Returns the real / imaginary part of the amplitude, depending on which one
// is non-zero.

// The amplitude is:  Ylm(theta, phi) - epsilon (-)^m Yl-m(theta, phi)
//                  = Ylm(theta, phi) - epsilon Ylm(theta, phi)*
//                  = N Plm(theta) (e^(i m phi) - epsilon e^(-i m phi))
//                  = N' Plm(theta) {cos,sin}(m phi)

// NOTE the phase is ignored, as different reflectivities don't interfere

double
likelihood::decay(int reflectivity, int l, int m, double theta, double phi) const
{
  double spherical = ROOT::Math::sph_legendre(l, m, theta);

  double factor = .5;
  if (m != 0)
    factor = sqrt(.5);    

  if (reflectivity == +1)
    return 2.*factor*spherical*sin(m*phi);
  else
    return 2.*factor*spherical*cos(m*phi);
}

complex<double>
likelihood::coherentSum(const coherent_waves& waves, const event& e) const
{
  complex<double> sum = 0;
  vector<wave>::const_iterator it;
  for (it = waves.waves.begin(); it != waves.waves.end(); it++)
    {
      sum += it->a * decay(waves.reflectivity,it->l,it->m,e.theta,e.phi);
    }

  return sum;
}

double
likelihood::probabilityDensity(const event& e) const
{
  double sum = 0;
  waveset::const_iterator it;
  for (it = ws.begin(); it != ws.end(); it++)
    {
      // norm = abs^2
      sum += norm(coherentSum(*it, e));
    }
  return sum;
}

map<int, double> weights;

double
likelihood::mcWeight(int reflectivity, int l1, int m1, int l2, int m2) const
{
  int id = reflectivity + ((l1 << 16) | (m1 << 12) | (l2 << 8) | (m2 << 4));
  if (weights.find(id) != weights.end())
    return weights[id];

  // Uses Kahan's summation
  double sum = 0;
  double c = 0;
  for (size_t i = 0; i < NMCEVENTS /*MCevents.size()*/; i++)
    {
      // Forming this sum as REAL sum and with no conjugate, because
      // the form of the decay amplitudes allows this.  This is not
      // the most general form!
      double y = (decay(reflectivity,l1,m1,MCevents[i].theta,MCevents[i].phi)
		  * decay(reflectivity,l2,m2,MCevents[i].theta,MCevents[i].phi)) - c;
      double t = sum + y;
      c = (t - sum) - y;  // compensation term.
      sum = t;
    }

  return (weights[id] = sum / NMCEVENTS);
}

double
likelihood::calc_likelihood() const
{
  double sumMC = 0;
  waveset::const_iterator it;
  for (it = ws.begin(); it != ws.end(); it++)
    {
      vector<wave>::const_iterator wave1, wave2;
      for (wave1 = it->waves.begin(); wave1 != it->waves.end(); wave1++)
	for (wave2 = it->waves.begin(); wave2 != it->waves.end(); wave2++)
	  sumMC +=  real(wave1->a*conj(wave2->a)
			 *mcWeight(it->reflectivity,
				   wave1->l, wave1->m,
				   wave2->l, wave2->m));
    }

  // Uses Kahan summation
  double sumRD = 0;
  double c = 0;
  for (size_t i = 0; i < NRDEVENTS /*RDevents.size()*/; i++)
    {
      double y = log(probabilityDensity(RDevents[i])) - c;
      double t = sumRD + y;
      c = (t - sumRD) - y;
      sumRD = t;
      //sumRD += y;
    }

  cout << "likelihood = " << sumRD << " - " << sumMC << endl;
  return sumRD - sumMC;
}


double
likelihood::operator()(const vector<double>& x) const
{
  size_t idx = 0;
  vector<coherent_waves>::iterator it;
  for (it = const_cast<likelihood*>(this)->ws.begin(); it != const_cast<likelihood*>(this)->ws.end(); it++)
    {
      vector<wave>::iterator wave;
      for (wave = it->waves.begin(); wave != it->waves.end(); wave++)
	{
	  wave->a = complex<double>(x[idx], x[idx + 1]);
	  idx += 2;
	}
    }

  return -this->calc_likelihood();
}



void
myFit()
{
  gRandom = new TRandom1;

  vector<wave> positive;
  positive.push_back(wave(2, 1));
  positive.push_back(wave(1, 1));

  vector<wave> negative;
  negative.push_back(wave(0,0));
  /*
  negative.push_back(wave(1,0));
  negative.push_back(wave(1,1));
  negative.push_back(wave(2,0));
  negative.push_back(wave(2,1));
  */

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

  //vector<event> RDevents(2500, event(0,0));
  for (int i = 0; i < NRDEVENTS; i++)
    {
      event e(acos(gRandom->Uniform(-1,1)), gRandom->Uniform(-4*atan(1),4*atan(1)));
      RDevents[i] = e;
    }

  //vector<event> MCevents(10000, event(0,0));
  for (int i = 0; i < NMCEVENTS; i++)
    {
      event e(acos(gRandom->Uniform(-1,1)), gRandom->Uniform(-4*atan(1),4*atan(1)));
      MCevents[i] =  e;
    }
  likelihood myL(ws); //, RDevents, MCevents);

  TStopwatch sw;

  TFitterMinuit* minuit = new TFitterMinuit();
  minuit->SetMinuitFCN(&myL);

  bool allowPositive = false;

  for (int i = 0; i < 5; i++)
    {
      sw.Start();

      minuit->Clear();
      /*
      minuit->SetParameter(0,"Rea(+,2,1)", 0, 0, 0, 0);
      minuit->SetParameter(1,"Ima(+,2,1)", 0, 0, 0, 0);
      minuit->SetParameter(2,"Rea(+,1,1)", 0, 0, 0, 0);
      minuit->SetParameter(3,"Ima(+,1,1)", 0, 0, 0, 0);
      */
      minuit->SetParameter(0,"Rea(+,2,1)", gRandom->Uniform(-10,10), 0.1, 0, 0);
      minuit->SetParameter(1,"Ima(+,2,1)", 0, 0.1, 0, 0);
      minuit->FixParameter(1);
      minuit->SetParameter(2,"Rea(+,1,1)", gRandom->Uniform(-10,10), 0.1, 0, 0);
      minuit->SetParameter(3,"Ima(+,1,1)", gRandom->Uniform(-10,10), 0.1, 0, 0);

      minuit->SetParameter(4,"Rea(-,0,0)", gRandom->Uniform(-10,10), 0.5, 0, 0);
      minuit->SetParameter(5,"Ima(-,0,0)",   0, 0.1, 0, 0);
      minuit->FixParameter(5);
      /*
      minuit->SetParameter(6,"Rea(-,1,0)", gRandom->Uniform(-10,10), 0.1, 0, 0);
      minuit->SetParameter(7,"Ima(-,1,0)", gRandom->Uniform(-10,10), 0.1, 0, 0);
      minuit->SetParameter(6,"Rea(-,1,1)", gRandom->Uniform(-10,10), 0.1, 0, 0);
      minuit->SetParameter(7,"Ima(-,1,1)", gRandom->Uniform(-10,10), 0.1, 0, 0);
      minuit->SetParameter(6,"Rea(-,2,0)", gRandom->Uniform(-10,10), 0.1, 0, 0);
      minuit->SetParameter(7,"Ima(-,2,0)", gRandom->Uniform(-10,10), 0.1, 0, 0);
      minuit->SetParameter(6,"Rea(-,2,1)", gRandom->Uniform(-10,10), 0.1, 0, 0);
      minuit->SetParameter(7,"Ima(-,2,1)", gRandom->Uniform(-10,10), 0.1, 0, 0);
      */

      minuit->CreateMinimizer();
      int iret = minuit->Minimize();
      sw.Stop();
      cout << "iret = " << iret << " after " << sw.CpuTime() << " s." << endl;
    }
}

#if 0
max = 4;
(*            s eta    l  m *)
indices = { { 1,  1, {{1, 1},
                      {2, 1} (*,
                      {2, 2},
                      {4, 1} *) } },
            (* {-1,  1, {{2, 2}}}, *)
            { 1, -1, {{0, 0},
                      {1, 0},
                      {2, 0},
                      {1, 1},
                      {2, 1}}}};

thetam[0] = 1/2;
thetam[m_] := 1/Sqrt[2]

a[s_,eta_,l_,m_] := rea[s,eta,l,m] + I ima[s,eta,l,m]

decay[eta_, l_,m_,theta_,phi_] := thetam[m] *
       (WignerD[{l,0,m},theta,phi] - eta (-1)^m WignerD[{l, 0, -m}, theta,phi])

coherentSum[indices_, theta_, phi_] := \
     (* indices = {s, eta, {{l,m}, {l,m}, {l, m}, ...}} *)
     Sum[ a[indices[[1]], indices[[2]], wave[[1]], wave[[2]]] \
         * decay[indices[[2]], wave[[1]], wave[[2]], theta, phi],
         {wave, indices[[3]]}];

probabilityDensity[theta_,phi_] := \
     Sum[Abs[coherentSum[idx, theta, phi]]^2, {idx, indices}]

NMCEvent = 100000;
mctheta = ArcCos[RandomVariate[UniformDistribution[{-1, 1}], NMCEvent]];
mcphi = RandomVariate[UniformDistribution[{0, 2Pi}], NMCEvent];

(* <<mcweights_not_backwards *)
Get["mcweights"]

mcweight[eta_,l1_,m1_,l2_,m2_] := \
   mcweight[eta,l1,m1,l2,m2] = mcw[eta,l1,m1,l2,m2] (* = \
   1/NMCEvent*Sum[decay[eta,l1,m1,mctheta[[i]],mcphi[[i]]]*
                  Conjugate[decay[eta,l2,m2,mctheta[[i]],mcphi[[i]]]],
                  {i,NMCEvent}] *)

     //(* data = ReadList["data.txt", Number]; *) (* "Hauke_PmEta.txt", Number]; *) (*dataHauke.txt", Number];*)

data = ReadList["Hauke_PmEta.txt", Number];
{mass, tprime, thetaAll, phiAll} = Transpose[Partition[data, 4]];

NEvent = Length[theta];

(*
theta = ArcCos[RandomVariate[UniformDistribution[{-1, 1}], NEvent]];
phi = RandomVariate[UniformDistribution[{0, 2Pi}], NEvent];
*)

(* probabilityDensity[theta[[1]], phi[[1]]] *)



Print["binning " <> ToString[Length[mass]] <> " events"]

mBinned = {};
tpBinned = {};
thetaBinned = {};
phiBinned = {};

For[i=0, i<40, i++,
    massbinlow = 0.7 (* 1.1 *) + i*0.05;
    massbinhigh = 0.7 (* 1.1 *) + (i+1)*0.05;

    mBin = {};
    tpBin = {};
    thetaBin = {};
    phiBin = {};
    NEvent = 0;
    For[j = 0, j<Length[mass], j++,
        If[mass[[j]] <= massbinhigh && mass[[j]] > massbinlow
           (* && thetaAll[[j]] > -0.8 *)
           && tprime[[j]] > 0.6 (* && tprime[[j]] < 0.6 *),
           AppendTo[mBin, mass[[j]]];
           AppendTo[tpBin, tprime[[j]]];
           AppendTo[thetaBin, thetaAll[[j]]];
           AppendTo[phiBin, phiAll[[j]]]
          ]
       ];
    AppendTo[mBinned, mBin];
    AppendTo[tpBinned, tpBin];
    AppendTo[thetaBinned, thetaBin];
    AppendTo[phiBinned, phiBin];
]

Print["data binned in " <> ToString[Length[thetaBinned]] <> " bins"];

val = {};
sol = {};

For[i=1, i < Length[thetaBinned], i++,
    NEvent = Length[thetaBinned[[i]]];
    Print[Length[thetaBinned[[i]]]];
    mcpart = Re@Sum[a[block[[1]], block[[2]], wave1[[1]], wave1[[2]]]
                    * Conjugate[
                       a[block[[1]], block[[2]], wave2[[1]], wave2[[2]]] *
                       mcweight[block[[2]], wave1[[1]], wave1[[2]],
                                            wave2[[1]], wave2[[2]]]],
                 {block, indices}, {wave1, block[[3]]}, {wave2, block[[3]]}];

    likelihood = (Sum[Re @ Log[probabilityDensity[thetaBinned[[i]][[j]],
                                                  phiBinned[[i]][[j]]]],
                      {j, 1, Length[thetaBinned[[i]]]}]
                  - mcpart);

    result = NMaximize[{likelihood, ima[1,1,2,1] == 0, ima[1,-1,1,0] == 0},

                                 { rea[1,1,1,1], ima[1,1,1,1],
                                   rea[1,1,2,1], ima[1,1,2,1],
                                   (* rea[1,1,4,1], ima[1,1,4,1], *)
                                   rea[1,-1,0,0], ima[1,-1,0,0],
                                   rea[1,-1,1,0], ima[1,-1,1,0],
                                   rea[1,-1,2,0], ima[1,-1,2,0],
                                   rea[1,-1,1,1], ima[1,-1,1,1],
                                   rea[1,-1,2,1], ima[1,-1,2,1] }];
(*
Table[{rea[s,eta,l,m], ima[s,eta,l,m]},
                     {s, 1, 1, 2}, {eta, -1, 1, 2},
                     {l, 0, lmax}, {m, (eta + 1)/2, l}]//Flatten];
*)
    AppendTo[val, {NEvent, N[(# Log[#] - #)&[NEvent]], result[[1]]}];
    AppendTo[sol, result[[2]]];

]

Print[val]
Print[sol]

Gwave = ListPlot[Abs[rea[1,1,4,1] + I ima[1,1,4,1]] /. sol,
                DisplayFunction->Identity, ImageSize->{400,300}];

Dwave = ListPlot[Abs[rea[1,1,2,1] + I ima[1,1,2,1]] /. sol,
                DisplayFunction->Identity, ImageSize->{400,300}];

Dwave2 = ListPlot[Abs[rea[1,1,2,2] + I ima[1,1,2,2]] /. sol,
                DisplayFunction->Identity, ImageSize->{400,300}];

Pwave = ListPlot[Abs[rea[1,1,1,1] + I ima[1,1,1,1]] /. sol,
                DisplayFunction->Identity, ImageSize->{400,300}];

phase = ListPlot[Arg[(rea[1,1,2,1] + I ima[1,1,2,1]) / (rea[1,1,1,1] + I ima[1,1,1,1])] /. sol,
                DisplayFunction->Identity, ImageSize->{400,300}];

Export["plots.png", GraphicsGrid[{{Pwave, Dwave}, {Dwave2, phase}}]]

predict1 = Plot3D[probabilityDensity[ArcCos[cth], phi] /. {rea[_,+1,_,_] -> 0, ima[_,+1,_,_] -> 0} /. sol[[13]],
                 {cth, -1, 1}, {phi, -Pi, Pi},
                 DisplayFunction->Identity, ImageSize->{400,300}];

predict2 = Plot3D[probabilityDensity[ArcCos[cth], phi] /. {rea[_,+1,_,_] -> 0, ima[_,+1,_,_] -> 0} /. sol[[14]],
                 {cth, -1, 1}, {phi, -Pi, Pi},
                 DisplayFunction->Identity, ImageSize->{400,300}];

predict3 = Plot3D[probabilityDensity[ArcCos[cth], phi] /. {rea[_,-1,_,_] -> 0, ima[_,-1,_,_] -> 0} /. sol[[13]],
                  {cth, -1, 1}, {phi, -Pi, Pi},
                  DisplayFunction->Identity, ImageSize->{400,300}];

predict4 = Plot3D[probabilityDensity[ArcCos[cth], phi] /. {rea[_,-1,_,_] -> 0, ima[_,-1,_,_] -> 0} /. sol[[14]],
                  {cth, -1, 1}, {phi, -Pi, Pi},
                  DisplayFunction->Identity, ImageSize->{400,300}];

Export["predictions.png", GraphicsGrid[{{predict1, predict2}, {predict3, predict4}}]];

data1 = Histogram3D[{Cos[thetaBinned[[12]]], phiBinned[[12]]} //Transpose,
                  DisplayFunction->Identity, ImageSize->{400,300}];
data2 = Histogram3D[{Cos[thetaBinned[[13]]], phiBinned[[13]]} //Transpose,
                  DisplayFunction->Identity, ImageSize->{400,300}];
data3 = Histogram3D[{Cos[thetaBinned[[14]]], phiBinned[[14]]} //Transpose,
                  DisplayFunction->Identity, ImageSize->{400,300}];
data4 = Histogram3D[{Cos[thetaBinned[[15]]], phiBinned[[15]]} //Transpose,
                  DisplayFunction->Identity, ImageSize->{400,300}];

Export["data.png", GraphicsGrid[{{data1, data2}, {data3, data4}}]];


(* Save["mcweights_not_backwards", mcw] *)
(* Save["mcweights", mcw]*)

Save["solution.m", sol]

Quit[]
#endif
