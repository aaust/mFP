void format(TH1D* hist, int i=0){

  hist->SetTitleFont(132);
  hist->SetStats(0);
  hist->SetFillStyle(0);
  hist->GetXaxis()->SetTitleSize(0.04);
  //hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleFont(132);
  hist->GetXaxis()->SetLabelFont(132); 
  hist->GetYaxis()->SetTitleSize(0.04);
  //hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleFont(132);
  hist->GetYaxis()->SetLabelFont(132);
  TGaxis::SetMaxDigits(3);
  hist->GetXaxis()->SetTitle("Invariant Mass of K^{+}K^{-} (GeV/c^{2})");
  hist->GetXaxis()->SetRangeUser(1.05,2.44);
  if (i==0)
    {
      hist->GetYaxis()->SetTitle("Intensity / (10MeV/c^{2})     ");
    }
  if (i==1)
    {
      hist->GetYaxis()->SetTitle("Phase       ");
    }
  if (i==3)
    {
      hist->GetYaxis()->SetTitle("Amplitude / (10MeV/c^{2})     ");
    }

}

void addtext(int i=0) {
  TText *data = new TText(.635, .85, "COMPASS 2009");
  data->SetNDC(true);
  data->SetTextColor(1);
  data->SetTextSize(0.045);
  data->SetTextFont(132);
  data->Draw();
  TLatex *data0 = new TLatex(.658, .80, "p p #rightarrow p_{f}K^{+}K^{-} p_{s}");
  data0->SetNDC(true);
  data0->SetTextColor(1);
  data0->SetTextSize(0.035);
  data0->SetTextFont(132);
  data0->Draw();
  TText *text = new TText(.4, .3, "preliminary");
    text->SetNDC(true);
    text->SetTextColor(16);
    text->SetTextSize(0.12);
    text->SetTextFont(132);
    text->SetTextAngle(20.);
    text->Draw();
}

void massdep(){

  TCanvas* c1 = new TCanvas("c1","Mass Dependent Fit",1200,1000);
  c1->Divide(2,2);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  
  TFile *out = TFile::Open("out16.root"); 
  TH1D *S = (TH1D *) gROOT->FindObject("hIntensityS0");
  c1->cd(1);
  format(S);
  S->Draw();
  TH1D *D = (TH1D *) gROOT->FindObject("hIntensityD0");
  c1->cd(4);
  format(D);
  D->Draw(); 
  TH1D *SD = (TH1D *) gROOT->FindObject("hPhaseS0D0");
  c1->cd(2);
  format(SD, 1);
  SD->Scale(-180./TMath::Pi());
  SD->Draw(); 

  TFile *massdep = TFile::Open("outMassDep.root"); 
  TH1D *Sfit = (TH1D *) gROOT->FindObject("hSwaveInt");
  TH1D *SBG = (TH1D *) gROOT->FindObject("hSwaveBG");
  TH1D *hf01370 = (TH1D *) gROOT->FindObject("hf01370");
  TH1D *hf01500 = (TH1D *) gROOT->FindObject("hf01500");
  TH1D *hf01710 = (TH1D *) gROOT->FindObject("hf01710");
  c1->cd(1);
  Sfit->SetLineColor(2);
  Sfit->SetLineWidth(2);
  Sfit->Draw("same");
  SBG->SetLineColor(4);
  SBG->Draw("same");
  hf01370->SetLineColor(6);
  hf01370->Draw("same");
  hf01500->SetLineColor(7);
  hf01500->Draw("same");
  hf01710->SetLineColor(8);
  hf01710->Draw("same");
  TH1D *Dfit = (TH1D *) gROOT->FindObject("hDwaveInt");
  TH1D *DBG = (TH1D *) gROOT->FindObject("hDwaveBG");
  TH1D *hf2 = (TH1D *) gROOT->FindObject("hf2");
  //  TH1D *ha2 = (TH1D *) gROOT->FindObject("ha2");
  TH1D *hf2prime = (TH1D *) gROOT->FindObject("hf2prime");
  //TH1D *hf2primeprime = (TH1D *) gROOT->FindObject("hf2primeprime");
  addtext(0);
  c1->cd(4);
  Dfit->SetLineColor(2);
  Dfit->SetLineWidth(2);
  Dfit->Draw("same"); 
  DBG->SetLineColor(4);
  DBG->Draw("same"); 
  hf2->SetLineColor(6);
  hf2->Draw("same");
  //ha2->SetLineColor(8);
  //ha2->Draw("same");
  hf2prime->SetLineColor(7);
  hf2prime->Draw("same");
  //hf2primeprime->SetLineColor(8);
  //hf2primeprime->Draw("same");
  TH1D *SDfit = (TH1D *) gROOT->FindObject("hPhaseSD");
  SDfit->Scale(-180./TMath::Pi());
  addtext(0);
  c1->cd(2);
  SDfit->SetLineColor(2);
  SDfit->SetLineWidth(2);
  SDfit->Draw("same"); 
  addtext(0);

  TCanvas* c2 = new TCanvas("c2","Real and Imaginary",1200,1000);
  c2->Divide(2,2);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  
  c2->cd(1);
  TH1D *ReS = (TH1D *) gROOT->FindObject("h2");
  format(ReS,3);
  ReS->SetTitle("Re(S_{0})");
  ReS->Draw();
  TH1D *ReSfit = (TH1D *) gROOT->FindObject("hSwaveRe");
  ReSfit->SetLineColor(2);
  ReSfit->SetLineWidth(2);
  ReSfit->Draw("same"); 
  addtext(0);

  c2->cd(2);
  TH1D *ReD = (TH1D *) gROOT->FindObject("h4");
  format(ReD,3);
  ReD->SetTitle("Re(D_{0})");
  ReD->Draw();
  TH1D *ReDfit = (TH1D *) gROOT->FindObject("hDwaveRe");
  ReDfit->SetLineColor(2);
  ReDfit->SetLineWidth(2);
  ReDfit->Draw("same"); 
  addtext(0);

  c2->cd(3);
  TH1D *ImSfit = (TH1D *) gROOT->FindObject("hSwaveIm");
  format(ImSfit,3);
  ImSfit->SetTitle("Im(S_{0})");
  ImSfit->SetLineColor(2);
  ImSfit->SetLineWidth(2);
  ImSfit->GetYaxis()->SetRangeUser(-1,1);
  ImSfit->Draw(); 
  addtext(0);

  c2->cd(4);
  TH1D *ImD = (TH1D *) gROOT->FindObject("h5");
  format(ImD,0);
  ImD->SetTitle("Im(D_{0})");
  ImD->Draw();
  TH1D *ImDfit = (TH1D *) gROOT->FindObject("hDwaveIm");
  ImDfit->SetLineColor(2);
  ImDfit->SetLineWidth(2);
  ImDfit->Draw("same"); 
  addtext(0);

}
