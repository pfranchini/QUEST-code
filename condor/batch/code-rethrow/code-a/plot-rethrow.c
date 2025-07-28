{

  char file_name[]="/home/franchinip/dataQUEST/QUEST/ND3/cosmics/output/cosmics-throw.root";

  ROOT::RDataFrame d1("Stepping", file_name);

  auto h_fxPos  = d1.Filter("fxPos").Histo1D({"h_fxPos", "x Pos", 100, -100., 100.},"fxPos");
  auto h_fyPos  = d1.Filter("fyPos").Histo1D({"h_fyPos", "y Pos", 100, -100., 100.},"fyPos");
  auto h_fzPos  = d1.Filter("fzPos").Histo1D({"h_fzPos", "z Pos", 100,    0., 150.},"fzPos");

  auto h_fxMom  = d1.Filter("fxMom").Histo1D({"h_fxMom", "x Mom", 100, -1.5, 1.5},"fxMom");
  auto h_fyMom  = d1.Filter("fyMom").Histo1D({"h_fyMom", "y Mom", 100, -1.5, 1.5},"fyMom");
  auto h_fzMom  = d1.Filter("fzMom").Histo1D({"h_fzMom", "z Mom", 100, -1.5, 1.5},"fzMom");

  TCanvas *c1 = new TCanvas("c1", "Cosmics", 1400 , 800);
  c1->Divide(3,2);

  c1->cd(1);
  gPad->SetLogy();
  h_fxPos->GetXaxis()->SetTitle("x [mm]");
  h_fxPos->Draw();
  c1->cd(2);
  gPad->SetLogy();
  h_fyPos->GetXaxis()->SetTitle("y [mm]");
  h_fyPos->Draw();
  c1->cd(3);
  gPad->SetLogy();
  h_fzPos->GetXaxis()->SetTitle("z [mm]");
  h_fzPos->Draw();

  c1->cd(4);
  gPad->SetLogy();
  h_fxMom->GetXaxis()->SetTitle("");
  h_fxMom->Draw();
  c1->cd(5);
  gPad->SetLogy();
  h_fyMom->GetXaxis()->SetTitle("");
  h_fyMom->Draw();
  c1->cd(6);
  gPad->SetLogy();
  h_fzMom->GetXaxis()->SetTitle("");
  h_fzMom->Draw();

  

  /*
  c1->cd(2);
  gPad->SetLogy();
  h_fEdep_low->GetXaxis()->SetTitle("Deposited energy [MeV]");
  h_fEdep_low->Draw();

  c1->cd(3);
  gPad->SetLogy();
  h_fEdep_ROI->GetXaxis()->SetTitle("Deposited energy [MeV]");
  h_fEdep_ROI->Draw();
  */
  //  c1->SaveAs("fEdep.png");  

  //h_fEdep_low->Draw();



}
