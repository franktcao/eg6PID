int derp(){

   TH2F * h2 = new TH2F("asdf","asdf",100,0,1,100,0,1);
   //h2->FillRandom("gaus",1000);
   h2->Fill(0.1,0.1);
   h2->Fill(0.1,0.1);
   h2->Fill(0.1,0.1);
   h2->Fill(0.1,0.1);
   h2->Fill(0.1,0.1);

TCanvas * c =new TCanvas();
   h2->Draw("colz");
 c->SaveAs("asdf.pdf");
   return 0;
}

