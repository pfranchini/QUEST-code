{

  char file_name[]="merge.root";

  ROOT::RDataFrame d1("Scoring", file_name);

  auto h_fEdep  = d1.Filter("fEdep>0&&fEdep<1.").Histo1D({"h_fEdep", "Full spectrum", 100, 0., 1.},"fEdep");
  auto h_fEdep_ROI  = d1.Filter("fEdep>25e-6&&fEdep<0.01").Histo1D({"h_fEdep_ROI", "ROI [25eV-10keV]", 50, 0., 0.01},"fEdep");
  auto h_fEdep_low  = d1.Filter("fEdep>0&&fEdep<=25e-6").Histo1D({"h_fEdep_low", "low energies [0-25eV]", 100, 0., 25e-6},"fEdep");

  
  std::cout << "Total events:         " << h_fEdep->GetEntries() << std::endl;
  std::cout << "ROI entries [-10keV]: " << h_fEdep_ROI->GetEntries() << std::endl;
  std::cout << "Low energy entries:   " << h_fEdep_low->GetEntries() << std::endl;

  TCanvas *c1 = new TCanvas("c1", "Cosmics", 1400 , 600);
  c1->Divide(3,1);

  c1->cd(1);
  gPad->SetLogy();
  h_fEdep->GetXaxis()->SetTitle("Deposited energy [MeV]");
  h_fEdep->Draw();

  c1->cd(2);
  gPad->SetLogy();
  h_fEdep_low->GetXaxis()->SetTitle("Deposited energy [MeV]");
  h_fEdep_low->Draw();

  c1->cd(3);
  gPad->SetLogy();
  h_fEdep_ROI->GetXaxis()->SetTitle("Deposited energy [MeV]");
  h_fEdep_ROI->Draw();

  c1->SaveAs("fEdep.png");  

  //h_fEdep_low->Draw();



}
