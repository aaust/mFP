void format(TH1D* hist, int i=0){

  hist->SetTitle("");
  hist->SetTitleFont(132);
  hist->SetStats(0);
  hist->SetFillStyle(0);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleFont(132);
  hist->GetXaxis()->SetLabelFont(132); 
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleFont(132);
  hist->GetYaxis()->SetLabelFont(132);
  TGaxis::SetMaxDigits(3);
  if (i==0)
    {
      hist->GetXaxis()->SetTitle("Invariant Mass of K^{+}K^{-} (GeV/c^{2})");
      hist->GetYaxis()->SetTitle("Intensity / 20MeV/c^{2}       ");
    }
  else
    hist->SetLineColor(2);

}

void massdep(){

  TCanvas* c1 = new TCanvas("c1","Mass Dependent Fit",1200,1200);
  c1->Divide(2,2);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  
  TFile *out = TFile::Open("out02.root"); 
  TH1D *S = (TH1D *) gROOT->FindObject("hIntensityS0");
  c1->cd(1);
  S->Draw();
  TH1D *D = (TH1D *) gROOT->FindObject("hIntensityD0");
  c1->cd(4);
  D->Draw(); 
  TH1D *SD = (TH1D *) gROOT->FindObject("hPhaseS0D0");
  c1->cd(2);
  SD->Draw(); 

  TFile *massdep = TFile::Open("outMassDep.root"); 
  TH1D *Sfit = (TH1D *) gROOT->FindObject("hSwaveInt");
  TH1D *SBG = (TH1D *) gROOT->FindObject("hSwaveBG");
  TH1D *hf01370 = (TH1D *) gROOT->FindObject("hf01370");
  TH1D *hf01500 = (TH1D *) gROOT->FindObject("hf01500");
  TH1D *hf01710 = (TH1D *) gROOT->FindObject("hf01710");
  c1->cd(1);
  Sfit->SetLineColor(2);
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
  TH1D *hf2prime = (TH1D *) gROOT->FindObject("hf2prime");
  c1->cd(4);
  Dfit->SetLineColor(2);
  Dfit->Draw("same"); 
  DBG->SetLineColor(4);
  DBG->Draw("same"); 
  hf2->SetLineColor(6);
  hf2->Draw("same");
  hf2prime->SetLineColor(7);
  hf2prime->Draw("same");
  TH1D *SDfit = (TH1D *) gROOT->FindObject("hPhaseSD");
  c1->cd(2);
  SDfit->SetLineColor(2);
  SDfit->Draw("same"); 

  TCanvas* c2 = new TCanvas("c2","Real and Imaginary",1200,1200);
  c2->Divide(2,2);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  
  c2->cd(1);
  TH1D *ReS = (TH1D *) gROOT->FindObject("h2");
  ReS->Draw();
  TH1D *ReSfit = (TH1D *) gROOT->FindObject("hSwaveRe");
  ReSfit->SetLineColor(2);
  ReSfit->Draw("same"); 

  c2->cd(2);
  TH1D *ReD = (TH1D *) gROOT->FindObject("h4");
  ReD->Draw();
  TH1D *ReDfit = (TH1D *) gROOT->FindObject("hDwaveRe");
  ReDfit->SetLineColor(2);
  ReDfit->Draw("same"); 

  c2->cd(3);
  TH1D *ImSfit = (TH1D *) gROOT->FindObject("hSwaveIm");
  ImSfit->SetLineColor(2);
  ImSfit->GetYaxis()->SetRangeUser(-1,1);
  ImSfit->Draw(); 

  c2->cd(4);
  TH1D *ImD = (TH1D *) gROOT->FindObject("h5");
  ImD->Draw();
  TH1D *ImDfit = (TH1D *) gROOT->FindObject("hDwaveIm");
  ImDfit->SetLineColor(2);
  ImDfit->Draw("same"); 



}
