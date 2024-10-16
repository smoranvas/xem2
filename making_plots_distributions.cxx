{
  gStyle->SetOptStat(0);
  TString Case = "HMS";
  TString file;
  if (Case=="SHMS") file="/cache/hallc/xem2/analysis/ONLINE/REPLAYS/SHMS/PRODUCTION/shms_replay_production_17876_-1.root";
  //if (Case=="HMS")  file="/work/smoran/xem2/data/HMS/hms_replay_production_4757_-1.root"; //5.878
  //if (Case=="HMS")  file="/work/smoran/xem2/data/HMS/hms_replay_production_4690_-1.root"; //6.600
  //if (Case=="HMS")  file="/work/smoran/xem2/data/HMS/hms_replay_production_4805_-1.root"; //5.360
  //if (Case=="HMS")  file="/work/smoran/xem2/data/HMS/hms_replay_production_4847_-1.root"; //4.78 
  //if (Case=="HMS")  file="/work/smoran/xem2/data/HMS/hms_replay_production_4876_-1.root"; //4.27
  //if (Case=="HMS")  file="/work/smoran/xem2/data/HMS/hms_replay_production_4903_-1.root"; //3.81
  //if (Case=="HMS")  file="/work/smoran/xem2/data/HMS/hms_replay_production_4926_-1.root"; //3.40
  //if (Case=="HMS")  file="/work/smoran/xem2/data/HMS/hms_replay_production_4949_-1.root"; //3.040
  //if (Case=="HMS")  file="/work/smoran/xem2/data/HMS/hms_replay_production_4970_-1.root"; //2.71
  //  if (Case=="HMS")  file="/work/smoran/xem2/data/HMS/hms_replay_production_5005_-1.root"; //2.42

  TString run= "4903" ; 
  if (Case=="HMS")  file= (const char*)Form("/work/smoran/xem2/data/HMS/hms_replay_production_%s_-1.root", (const char*)run);

  const Double_t p_spec        = 3.81; // in GeV, it is the central momentum

  TFile *f = new TFile((const char*)file);
  TTree *TSP, *TSH;

  if (Case=="SHMS"){TSP= (TTree*)f->Get("TSP");} // P Scaler Data
  if (Case=="HMS"){ TSH= (TTree*)f->Get("TSH");} // H Scaler Data

  TTree *T   = (TTree*)f->Get("T");      // Hall A Analizer Output DST 
  TTree *E   = (TTree*)f->Get("E");      // Hall A Epics Data

  Double_t ptardp , ptarth , ptarph , Cal , ngcer , ptary, xptar, yptar, Q2, Xbj, ytar, W, nu , xtar , AvgCurr;
  Double_t xfp, yfp, xpfp, ypfp;  // FOCAL PLANE VARIABLES
  Double_t theta, Ep;
  Double_t theta_central ;
  if (Case=="SHMS") {theta_central= 8.02*TMath::Pi()/180.  ;}  //8 degrees, this value is in 'ecSHMS_Angle' variable
  if (Case=="HMS"){ theta_central = 20.0*TMath::Pi()/180.  ;   }

  const Double_t Mp            = 0.938; // Proton's mass
  const Double_t Eb            = 10.541650; //beam's energy

  //Set branches of the tree 
  if (Case=="SHMS"){
    T->SetBranchAddress("P.gtr.dp", &ptardp); //data delta (momentum acceptance)
    T->SetBranchAddress("P.gtr.th", &xptar); //data target theta (dx/dz)    THIS IS 'xptar'
    T->SetBranchAddress("P.gtr.ph", &yptar); //data target phi (dy/dz)      THIS IS 'yptar'
    T->SetBranchAddress("P.cal.etottracknorm", &Cal); //total deposited energy in the calorimeter normalized by track's momentum
    T->SetBranchAddress("P.ngcer.npeSum", &ngcer); //total # of photoelectrons SHMS Noble Gas cherenkov
    T->SetBranchAddress("P.kin.Q2"   , &Q2);
    T->SetBranchAddress("P.kin.x_bj" , &Xbj);
    T->SetBranchAddress("P.kin.W" , &W);
    T->SetBranchAddress("P.kin.nu" , &nu);
    T->SetBranchAddress("P.gtr.y" , &ytar); 
    T->SetBranchAddress("P.dc.x_fp" , &xfp);
    T->SetBranchAddress("P.dc.y_fp" , &yfp);
    T->SetBranchAddress("P.dc.xp_fp" , &xpfp);
    T->SetBranchAddress("P.dc.yp_fp" , &ypfp);
  }
  else if (Case=="HMS"){
    T->SetBranchAddress("H.gtr.dp", &ptardp); //data delta (momentum acceptance)
    T->SetBranchAddress("H.gtr.th", &xptar); //data target theta (dx/dz)    THIS IS 'xptar' or  'ptarth''
    T->SetBranchAddress("H.gtr.ph", &yptar); //data target phi (dy/dz)      THIS IS 'yptar' or  'ptarph'
    T->SetBranchAddress("H.cal.etottracknorm", &Cal); //total deposited energy in the calorimeter normalized by track's momentum
    T->SetBranchAddress("H.kin.Q2"   , &Q2);
    T->SetBranchAddress("H.kin.x_bj" , &Xbj);
    T->SetBranchAddress("H.kin.W" , &W);
    T->SetBranchAddress("H.kin.nu" , &nu);
    T->SetBranchAddress("H.gtr.y" , &ytar);
    T->SetBranchAddress("H.gtr.x" , &xtar);
    T->SetBranchAddress("H.dc.x_fp" , &xfp);
    T->SetBranchAddress("H.dc.y_fp" , &yfp);
    T->SetBranchAddress("H.dc.xp_fp" , &xpfp);
    T->SetBranchAddress("H.dc.yp_fp" , &ypfp);
    T->SetBranchAddress("H.cer.npeSum", &ngcer);
    T->SetBranchAddress("H.bcm.bcm4a.AvgCurrent", &AvgCurr);

//theta = TMath::ACos(( cos(theta_central) + yptar * sin(theta_central) ) / TMath::Sqrt( 1. + xptar * xptar + yptar * yptar)); //ThetaLab


}

 TCanvas *c= new TCanvas("c","assigment4");
 //c->Divide(2,1);
 Int_t nbins=170;

 TH1F *Q2_h = new TH1F("Q2", "Q2" , nbins, 0,2.5);
 TH1F *Xb_h = new TH1F("Xb", "Xb" , nbins, 0.1 , 2.5);
 TH1F *W_h  = new TH1F("W", "W"   , nbins, -0.5,2.5);
 TH1F *nu_h = new TH1F("nu","nu"  , nbins , -0.5, 2.5 );
 TH1F *tmp  = new TH1F("tmp","tmp"  , nbins , -10, 100.5 );
 TH1F *ecal_h = new TH1F("ecal","ecal"  , nbins , -0.5, 2.5 );
 TH1F *cer_h  = new TH1F("cer","cer"  , nbins , -0.5, 30.5 );
 TH1F *xfp_h  = new TH1F("xfp","xfp"  , nbins , -60.0, 60.0 );

 TH2F *Q2_comp = new TH2F("Q2_comp","theta v/s delta"       ,nbins,-9,9,nbins,17,23);
 TH2F *Xb_comp = new TH2F("Xb_comp","Xb cal v/s Xb tree"    ,nbins,0,10,nbins,0,10);
 TH2F *W_comp  = new TH2F("W_comp","W cal v/s W tree"       ,nbins,0,10,nbins,0,10);
 TH2F *nu_comp = new TH2F("nu_comp","nu cal v/s nu tree"    ,nbins,0,20,nbins,0,45);
 TH2F *xytar_comp = new TH2F("xytar_comp","xtar v/s ytar"    ,nbins,-15,15,nbins,-15,15);


 Q2_h->GetXaxis()->SetTitle("cal.etottracknorm");
 Xb_h->GetXaxis()->SetTitle("Xbj");
 W_h->GetXaxis()->SetTitle("W");
 nu_h->GetXaxis()->SetTitle("nu");
 tmp->GetXaxis()->SetTitle("AvgCurr");

 Q2_h->GetYaxis()->SetTitle("counts");
 Xb_h->GetYaxis()->SetTitle("counts");
 W_h->GetYaxis()->SetTitle("counts");
 nu_h->GetYaxis()->SetTitle("counts");
 tmp ->GetYaxis()->SetTitle("counts");

 Q2_comp->GetXaxis()->SetTitle("delta");
 Q2_comp->GetYaxis()->SetTitle("theta");
 Xb_comp->GetXaxis()->SetTitle("Xbj calc");
 Xb_comp->GetYaxis()->SetTitle("Xbj tree");
 W_comp->GetXaxis()->SetTitle("W calc");
 W_comp->GetYaxis()->SetTitle("W tree");
 nu_comp->GetXaxis()->SetTitle("nu calc");
 nu_comp->GetYaxis()->SetTitle("nu tree");
 xytar_comp->GetXaxis()->SetTitle("Ytar");
 xytar_comp->GetYaxis()->SetTitle("Xtar");

 Long64_t nentries = T->GetEntries();
 for (int i=0;i<nentries;i++){
   T->GetEntry(i);
   Double_t theta;
   if (Case =="SHMS") theta = TMath::ACos(( cos(theta_central) - yptar * sin(theta_central) ) / TMath::Sqrt( 1. + xptar * xptar + yptar * yptar)); //ThetaLab
   if (Case=="HMS")   theta = TMath::ACos(( cos(theta_central) + yptar * sin(theta_central) ) / TMath::Sqrt( 1. + xptar * xptar + yptar * yptar)); //ThetaLab

   Double_t Ep      = p_spec*(1.0+0.01*ptardp); // E prime
   Double_t Q2_calc = 4*Eb*Ep*sin(theta/2.)*sin(theta/2.);
   Double_t nu_calc = Eb-Ep;
   Double_t W_calc  = TMath::Sqrt(Mp*Mp + 2*Mp*nu_calc - Q2_calc)  ; //invariant mass electron-nucleon interaction
   Double_t Xb_calc = Q2_calc/2./Mp/nu_calc;
   Double_t delta   = ptardp ; 
   //   if (abs(ptardp) && Cal>0.7 &&  ngcer>2 && abs(xptar)<0.085 && abs(yptar)<0.032  ){
   if (true) {
     /*   ecal_h->Fill(Cal);
   cer_h->Fill(ngcer);
   xfp_h->Fill(xfp);
   Xb_h->Fill(Xb_calc);
   tmp->Fill(AvgCurr);*/
   xytar_comp->Fill( ytar,xtar );
 }
  }

 cout<< "RUN NUMBER: " << run <<endl;
 cout<< "#ECAL histo: " << endl;

 std::cout << "bin_centers_ecalRUN_"<<run<<" = np.array([";
 for (int i = 1; i <= nbins; ++i) {
   double binCenter = ecal_h->GetBinCenter(i);
   std::cout << binCenter;
   if (i < nbins) std::cout << ", ";
 }
 std::cout << "])\n\n";

 std::cout << "bin_contents_ecalRUN_"<<run<<" = np.array([";
 for (int i = 1; i <= nbins; ++i) {
   double binContent = ecal_h->GetBinContent(i);
   std::cout << binContent;
   if (i < nbins) std::cout << ", ";
 }
 std::cout << "])\n\n" << std::endl;


 cout<< "#CER histo: " << endl;

 std::cout << "bin_centers_cerRUN_"<<run<<" = np.array([";
 for (int i = 1; i <= nbins; ++i) {
   double binCenter = cer_h->GetBinCenter(i);
   std::cout << binCenter;
   if (i < nbins) std::cout << ", ";
 }
 std::cout << "])\n\n";

 std::cout << "bin_contents_cerRUN_"<<run<<" = np.array([";
 for (int i = 1; i <= nbins; ++i) {
   double binContent = cer_h->GetBinContent(i);
   std::cout << binContent;
   if (i < nbins) std::cout << ", ";
 }
 std::cout << "])\n\n" << std::endl;


 cout<< "#XFP histo: " << endl;

 std::cout << "bin_centers_xfpRUN_"<<run<<" = np.array([";
 for (int i = 1; i <= nbins; ++i) {
   double binCenter = xfp_h->GetBinCenter(i);
   std::cout << binCenter;
   if (i < nbins) std::cout << ", ";
 }
 std::cout << "])\n\n";

 std::cout << "bin_contents_xfpRUN_"<<run<<" = np.array([";
 for (int i = 1; i <= nbins; ++i) {
   double binContent = xfp_h->GetBinContent(i);
   std::cout << binContent;
   if (i < nbins) std::cout << ", ";
 }
 std::cout << "])\n\n" << std::endl;

 cout<< "#Xb histo: " << endl;

 std::cout << "bin_centers_xRUN_"<<run<<" = np.array([";
 for (int i = 1; i <= nbins; ++i) {
   double binCenter = Xb_h->GetBinCenter(i);
   std::cout << binCenter;
   if (i < nbins) std::cout << ", ";
 }
 std::cout << "])\n\n";

 std::cout << "bin_contents_xRUN_"<<run<<" = np.array([";
 for (int i = 1; i <= nbins; ++i) {
   double binContent = Xb_h->GetBinContent(i);
   std::cout << binContent;
   if (i < nbins) std::cout << ", ";
 }
 std::cout << "])\n\n" << std::endl;





 Double_t min=0.01;
 Q2_h->SetMinimum(min);
 Xb_h->SetMinimum(min);
 nu_h->SetMinimum(min);
 W_h->SetMinimum(min);
 tmp->SetMinimum(min);

 const Int_t nColors = 18;
 Int_t colors[nColors] = {
   TColor::GetColor("#ff0000"), // Red
   TColor::GetColor("#ff3300"), // Orange-Red
   TColor::GetColor("#ff6600"), // Orange
   TColor::GetColor("#ff9900"), // Yellow-Orange
   TColor::GetColor("#ffcc00"), // Yellow
   TColor::GetColor("#ffff00"), // Light Yellow
   TColor::GetColor("#ccff00"), // Yellow-Green
   TColor::GetColor("#99ff00"), // Light Green
   TColor::GetColor("#66ff00"), // Green
   TColor::GetColor("#33ff00"), // Spring Green
   TColor::GetColor("#00ff00"), // Bright Green
   TColor::GetColor("#00ff33"), // Sea Green
   TColor::GetColor("#00ff66"), // Aqua
   TColor::GetColor("#00ff99"), // Cyan
   TColor::GetColor("#00ffcc"), // Light Cyan
   TColor::GetColor("#00ffff"), // Sky Blue
   TColor::GetColor("#00ccff"), // Light Blue
   TColor::GetColor("#0099ff")  // Blue
 };
 const Int_t nXb = 18;
 const Double_t Xb_values[nXb] = {0.2, 0.28, 0.36, 0.44, 0.52, 0.60, 0.68, 0.76, 0.84, 0.92, 1.0, 1.08, 1.16, 1.24, 1.32, 1.40, 1.48, 1.56};
 const Int_t nPoints = 500;
 Double_t x_values[nPoints];
 Double_t y_values[nPoints];
 Double_t delta_values[nPoints];

 TLegend *legend = new TLegend(0.8, 0.1, 0.99, 0.9);
 legend->SetNColumns(1);

 for (Int_t i = 0; i < nPoints; i++) {
   x_values[i] = 0.1 + (Eb - 0.1 - 0.1) * i / (nPoints - 1);  // Avoid division by zero
   delta_values[i] = (x_values[i] / p_spec - 1) * 100;  
 }


 //  c->cd(1); c_1->SetLogy(); Q2_h->Draw(); 
 //  c->cd(2);  
 xytar_comp->Draw("COLZ");







 /*
Q2_comp->Draw("colz");

  for (Int_t j = 0; j < nXb; j++) {
    Double_t Xb = Xb_values[j];
    Double_t C = (4 * Eb) / (2 * Mp * Xb);

    for (Int_t i = 0; i < nPoints; i++) {
      y_values[i] = 2 * TMath::ASin(TMath::Sqrt((2 * Mp * Xb * (Eb- x_values[i])) / (4 * Eb * x_values[i])));
      y_values[i] *= 180 / TMath::Pi();  // Convert to degrees                                                                                                                                                                                          
    }

    TGraph *graph = new TGraph(nPoints,  delta_values , y_values);
    //    graph->SetLineColor(colors[j]);
    graph->SetLineColor(kBlack);
    graph->SetLineWidth(2);
    graph->SetTitle("Lines of Constant Xb;#delta;#theta (degrees)");

    if (j == 0) {
      graph->Draw("L same");
    } else {
      graph->Draw("L same");
    }

    legend->AddEntry(graph,(const char*) Form("Xb = %.2f", Xb), "l");
  }


  c->SetLogz();
 */

  /*
  c->cd(3); c_3->SetLogy(); Xb_h->Draw();
  c->cd(4);  Xb_comp->Draw("colz"); 
  c->cd(5); c_5->SetLogy(); W_h->Draw();
  c->cd(6);  W_comp->Draw("colz");
  c->cd(7); c_7->SetLogy(); nu_h->Draw();
  c->cd(8);  nu_comp->Draw("colz");

  c_2->SetLogz();
  c_4->SetLogz();
  c_6->SetLogz();
  c_8->SetLogz();

  */
  if (Case=="SHMS") {c->SaveAs("assigment4_SHMS.pdf");  c->SaveAs("assigment4_SHMS.root");}
  if (Case=="HMS")  {c->SaveAs("distributions.pdf");    c->SaveAs("distributions.root");}


}
