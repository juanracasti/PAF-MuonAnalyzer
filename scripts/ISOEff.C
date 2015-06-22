void ISOEff(TString Sample = "DR74X_50ns_MC_DY", float dR = 0.4) 
{ 

  TFile *f  = TFile::Open("../files/" + Sample + ".root");

  const int nVar = 3;
  const int nISO = 4;
  const int nMu  = 2;

  TH1F* h[nVar][nISO+1][nMu];

  f->cd();

  //if (dR == 0.4) {
    
    h[0][0][0] = (TH1F*)gDirectory->Get("h_Eff_pt_TightID_Mu1");
    h[0][1][0] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO04_Mu1");
    h[0][2][0] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO04dBeta_Mu1");
    h[0][3][0] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO04PFWeighted_Mu1");
    h[0][4][0] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO04PUPPI_Mu1");

    h[0][0][1] = (TH1F*)gDirectory->Get("h_Eff_pt_TightID_Mu2");
    h[0][1][1] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO04_Mu2");
    h[0][2][1] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO04dBeta_Mu2");
    h[0][3][1] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO04PFWeighted_Mu2");
    h[0][4][1] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO04PUPPI_Mu2");

    
    h[1][0][0] = (TH1F*)gDirectory->Get("h_Eff_eta_TightID_Mu1");
    h[1][1][0] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO04_Mu1");
    h[1][2][0] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO04dBeta_Mu1");
    h[1][3][0] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO04PFWeighted_Mu1");
    h[1][4][0] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO04PUPPI_Mu1");
    
    h[1][0][1] = (TH1F*)gDirectory->Get("h_Eff_eta_TightID_Mu2");
    h[1][1][1] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO04_Mu2");
    h[1][2][1] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO04dBeta_Mu2");
    h[1][3][1] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO04PFWeighted_Mu2");
    h[1][4][1] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO04PUPPI_Mu2");
    
    
    h[2][0][0] = (TH1F*)gDirectory->Get("h_Eff_npv_TightID_Mu1");
    h[2][1][0] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO04_Mu1");
    h[2][2][0] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO04dBeta_Mu1");
    h[2][3][0] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO04PFWeighted_Mu1");
    h[2][4][0] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO04PUPPI_Mu1");
    
    h[2][0][1] = (TH1F*)gDirectory->Get("h_Eff_npv_TightID_Mu2");
    h[2][1][1] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO04_Mu2");
    h[2][2][1] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO04dBeta_Mu2");
    h[2][3][1] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO04PFWeighted_Mu2");
    h[2][4][1] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO04PUPPI_Mu2"); 
    
    //}

  // else if (dR == 0.3) {
    
  //   h[0][0][0] = (TH1F*)gDirectory->Get("h_Eff_pt_TightID_Mu1");
  //   h[0][1][0] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO03_Mu1");
  //   h[0][2][0] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO03dBeta_Mu1");
  //   h[0][3][0] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO03PFWeighted_Mu1");
  //   h[0][4][0] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO03PUPPI_Mu1");

  //   h[0][0][1] = (TH1F*)gDirectory->Get("h_Eff_pt_TightID_Mu2");
  //   h[0][1][1] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO03_Mu2");
  //   h[0][2][1] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO03dBeta_Mu2");
  //   h[0][3][1] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO03PFWeighted_Mu2");
  //   h[0][4][1] = (TH1F*)gDirectory->Get("h_Eff_pt_TightISO03PUPPI_Mu2");

    
  //   h[1][0][0] = (TH1F*)gDirectory->Get("h_Eff_eta_TightID_Mu1");
  //   h[1][1][0] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO03_Mu1");
  //   h[1][2][0] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO03dBeta_Mu1");
  //   h[1][3][0] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO03PFWeighted_Mu1");
  //   h[1][4][0] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO03PUPPI_Mu1");
    
  //   h[1][0][1] = (TH1F*)gDirectory->Get("h_Eff_eta_TightID_Mu2");
  //   h[1][1][1] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO03_Mu2");
  //   h[1][2][1] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO03dBeta_Mu2");
  //   h[1][3][1] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO03PFWeighted_Mu2");
  //   h[1][4][1] = (TH1F*)gDirectory->Get("h_Eff_eta_TightISO03PUPPI_Mu2");
    
    
  //   h[2][0][0] = (TH1F*)gDirectory->Get("h_Eff_npv_TightID_Mu1");
  //   h[2][1][0] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO03_Mu1");
  //   h[2][2][0] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO03dBeta_Mu1");
  //   h[2][3][0] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO03PFWeighted_Mu1");
  //   h[2][4][0] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO03PUPPI_Mu1");
    
  //   h[2][0][1] = (TH1F*)gDirectory->Get("h_Eff_npv_TightID_Mu2");
  //   h[2][1][1] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO03_Mu2");
  //   h[2][2][1] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO03dBeta_Mu2");
  //   h[2][3][1] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO03PFWeighted_Mu2");
  //   h[2][4][1] = (TH1F*)gDirectory->Get("h_Eff_npv_TightISO03PUPPI_Mu2"); 
    
  // }
  
  TGraphAsymmErrors* r[nVar][nISO][nMu];
  TH1F*              h1[nVar][nISO+1][nMu];

  double ptbins1[24] = {10,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,100,
  		       110,120,130,140,150,170,200};
  double ptbins2[11] = {10,20,25,30,35,40,50,60,75,100,200};
  double ptbins3[6]  = {10,20,40,60,100,200};
		       
  double etabins[14] = {-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,-0.2,0.2,0.3,0.9,1.2,1.6,2.1,2.4};

  double npvbins[11] = {0,10,12,14,16,18,20,22,25,30,40};

  
  for (int iISO = 0; iISO < nISO+1; iISO++) {
    for (int iMu = 0; iMu < nMu; iMu++)     {
      
      h1[0][iISO][iMu] = (TH1F*) h[0][iISO][iMu]->Rebin(10, "h1", ptbins2);
      h1[1][iISO][iMu] = (TH1F*) h[1][iISO][iMu]->Rebin(13, "h1", etabins);
      h1[2][iISO][iMu] = (TH1F*) h[2][iISO][iMu]->Rebin(10, "h1", npvbins);
      
    }
  }
   
  for (int iVar = 0; iVar < nVar; iVar++)   {
    for (int iISO = 0; iISO < nISO; iISO++) {
      for (int iMu = 0; iMu < nMu; iMu++)   {

  	r[iVar][iISO][iMu] = new TGraphAsymmErrors(h1[iVar][iISO+1][iMu], h1[iVar][0][iMu]); 

      }
    }
  }
  
  
  ////////////////////////////////////
  /////////////// DRAW ///////////////
  ////////////////////////////////////

  //gStyle->SetOptStat("e");
  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.80);
  gStyle->SetStatY(0.80);
  gStyle->SetStatW(0.19);
  gStyle->SetStatH(0.10);

  TCanvas* c[nVar][nMu];
  TLegend* l2;
  TLegend* l3;
  TMultiGraph* mg[nVar][nMu];


  if (1) {

    c[0][0] = new TCanvas("Efficiencies vs pt Mu1",   "Efficiencies vs pt Mu1",   750, 750);
    c[0][1] = new TCanvas("Efficiencies vs pt Mu2",   "Efficiencies vs pt Mu2",   750, 750);
    c[1][0] = new TCanvas("Efficiencies vs eta Mu1",  "Efficiencies vs eta Mu1",  750, 750);
    c[1][1] = new TCanvas("Efficiencies vs eta Mu2",  "Efficiencies vs eta Mu2",  750, 750);
    c[2][0] = new TCanvas("Efficiencies vs npv Mu1",  "Efficiencies vs npv Mu1",  750, 750);
    c[2][1] = new TCanvas("Efficiencies vs npv Mu2",  "Efficiencies vs npv Mu2",  750, 750);


    for (int iVar=0; iVar<nVar; iVar++) {
      for (int iMu=0; iMu<nMu; iMu++)   {
	mg[iVar][iMu] = new TMultiGraph();
      }
    }
    

    for (int iVar=0; iVar<nVar; iVar++) {
      for (int iMu=0; iMu<nMu; iMu++)   {
	    
    	c[iVar][iMu]->Range(-0.8247861,-0.8247861,0.8247861,0.8247861);
    	c[iVar][iMu]->SetLeftMargin(0.15);
    	c[iVar][iMu]->SetRightMargin(0.15);
    	c[iVar][iMu]->SetTopMargin(0.15);
    	c[iVar][iMu]->SetBottomMargin(0.15);
    	//c[iVar][iMu]->SetLogy();
	
      }
    }
    

    l2 = new TLegend(0.60, 0.20, 0.83, 0.35);
    l2->SetBorderSize(1);
    l2->SetFillStyle(1001);
    l2->SetTextFont(42);
    l2->AddEntry(r[0][0][0], "Standard",              "lp");
    l2->AddEntry(r[0][1][0], "#Delta#beta corrected", "lp");
    l2->AddEntry(r[0][2][0], "PF-Weighted",           "lp");
    l2->AddEntry(r[0][3][0], "PUPPI",                 "lp");


    l3 = new TLegend(0.60, 0.68, 0.83, 0.83); 
    l3->SetBorderSize(1);
    l3->SetFillStyle(1001);
    l3->SetTextFont(42);
    l3->AddEntry(r[0][0][0], "Standard",              "lp");
    l3->AddEntry(r[0][1][0], "#Delta#beta corrected", "lp");
    l3->AddEntry(r[0][2][0], "PF-Weighted",           "lp");
    l3->AddEntry(r[0][3][0], "PUPPI",                 "lp");


    for (int iMu=0; iMu<nMu; iMu++)
      {
    	c[0][iMu]->cd();  
    	r[0][0][iMu]->GetXaxis()->SetTitle("p_{T}");
	
    	c[1][iMu]->cd();  
    	r[1][0][iMu]->GetXaxis()->SetTitle("#eta");
	
    	c[2][iMu]->cd();  
    	r[2][0][iMu]->GetXaxis()->SetTitle("NPV");
      }
    
    for (int iVar=0; iVar<nVar; iVar++)   {
      for (int iMu=0; iMu<nMu; iMu++)   {
	
    	c[iVar][iMu]->cd();
	
    	r[iVar][0][iMu]->SetLineColor(4);
    	r[iVar][0][iMu]->SetLineWidth(2);
    	r[iVar][0][iMu]->SetMarkerColor(4);
    	r[iVar][0][iMu]->SetMarkerSize(0.8);
    	r[iVar][0][iMu]->SetMarkerStyle(kFullTriangleUp);
    	r[iVar][0][iMu]->SetFillStyle(3001);
    	r[iVar][0][iMu]->SetFillColor(4);	   
    	r[iVar][1][iMu]->SetLineColor(6);
    	r[iVar][1][iMu]->SetLineWidth(2);
    	r[iVar][1][iMu]->SetMarkerColor(6);
    	r[iVar][1][iMu]->SetMarkerSize(0.8);
    	r[iVar][1][iMu]->SetMarkerStyle(kFullSquare);
    	r[iVar][1][iMu]->SetFillStyle(3001);
    	r[iVar][1][iMu]->SetFillColor(6);
    	r[iVar][2][iMu]->SetLineColor(8);
    	r[iVar][2][iMu]->SetLineWidth(2);
    	r[iVar][2][iMu]->SetMarkerColor(8);
    	r[iVar][2][iMu]->SetMarkerSize(1.0);
    	r[iVar][2][iMu]->SetMarkerStyle(kOpenSquare);
    	r[iVar][2][iMu]->SetFillStyle(3001);
    	r[iVar][2][iMu]->SetFillColor(8);
    	r[iVar][3][iMu]->SetLineColor(kOrange+7);
    	r[iVar][3][iMu]->SetLineWidth(2);
    	r[iVar][3][iMu]->SetMarkerColor(kOrange+7);
    	r[iVar][3][iMu]->SetMarkerSize(1.0);
    	r[iVar][3][iMu]->SetMarkerStyle(kOpenCircle);
    	r[iVar][3][iMu]->SetFillStyle(3001);
    	r[iVar][3][iMu]->SetFillColor(kOrange+7);
	
    	mg[iVar][iMu]->Add(r[iVar][0][iMu], "p");
    	mg[iVar][iMu]->Add(r[iVar][2][iMu], "p");
    	mg[iVar][iMu]->Add(r[iVar][1][iMu], "p");
    	mg[iVar][iMu]->Add(r[iVar][3][iMu], "p");	    

    	c[iVar][iMu] ->Clear();
    	mg[iVar][iMu]->Draw("a");
    	c[iVar][iMu] ->Modified();
    	c[iVar][iMu] ->Update();

    	if (iVar==0) {  
    	  mg[iVar][iMu]->GetXaxis()->SetTitle("p_{T}");
    	  mg[iVar][iMu]->GetXaxis()->SetRangeUser(0.0,200.0); 
    	}
    	if (iVar==1) {
    	  mg[iVar][iMu]->GetXaxis()->SetTitle("#eta");
    	  mg[iVar][iMu]->GetXaxis()->SetRangeUser(-2.5,2.5); 
    	}
    	if (iVar==2) {
    	  mg[iVar][iMu]->GetXaxis()->SetTitle("\#PV");
    	  mg[iVar][iMu]->GetXaxis()->SetRangeUser(0.0,45.0); 
    	}
	
    	mg[iVar][iMu]->GetXaxis()->SetLabelOffset(0.015);
    	mg[iVar][iMu]->GetXaxis()->SetLabelSize(0.025);
    	mg[iVar][iMu]->GetXaxis()->SetTitleSize(0.035);
    	mg[iVar][iMu]->GetXaxis()->SetTitleOffset(1.2);
    	mg[iVar][iMu]->GetYaxis()->SetTitle("Efficiency");
    	mg[iVar][iMu]->GetYaxis()->SetLabelOffset(0.015);
    	mg[iVar][iMu]->GetYaxis()->SetLabelSize(0.025);
    	mg[iVar][iMu]->GetYaxis()->SetTitleSize(0.035);
    	mg[iVar][iMu]->GetYaxis()->SetTitleOffset(1.6);
	
    	mg[iVar][iMu]->GetYaxis()->SetRangeUser(0.65,1.05);
    	l2->Draw(); 
	
    	//c[iVar]->Update();
	
      }
    }
    
  }
  
  
}
