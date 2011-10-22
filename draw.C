{
  c= new TCanvas; c->Divide(3,2);
  c->cd(1); h0->Draw(); hDwaveRe->Draw("same");
  c->cd(2); h2->Draw(); hPwaveRe->Draw("same");
  c->cd(3); h3->Draw(); hPwaveIm->Draw("same");
  //c->cd(4); hPhaseD->Draw();
  c->cd(5); h4->Draw(); hGwaveRe->Draw("same");
  c->cd(6); h5->Draw(); hGwaveIm->Draw("same");

  c= new TCanvas; c->Divide(3,2);
  c->cd(1); h0->Draw(); hDwaveEvolution->Draw("samecol"); hDwaveRe->Draw("same");
  c->cd(2); h2->Draw(); hPwaveReEvolution->Draw("samecol"); hPwaveRe->Draw("same");
  c->cd(3); h3->Draw(); hPwaveImEvolution->Draw("samecolz"); hPwaveIm->Draw("same");
  c->cd(4); hPhaseDP->Draw();
  c->cd(5); h4->Draw(); hGwaveReEvolution->Draw("samecolz"); hGwaveRe->Draw("same");
  c->cd(6); h5->Draw(); hGwaveImEvolution->Draw("samecolz"); hGwaveIm->Draw("same");
}

