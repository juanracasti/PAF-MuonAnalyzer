void IDEffAllMu(TString Sample = "DR74X_50ns_MC_DY") 
{ 

  TFile *f  = TFile::Open("../files/" + Sample + ".root");

  const int nVar = 3;
  const int nID = 6;

  TH1F* h[nVar][nID+1];

  TString Var[nVar] = {
    "pt",
    "eta",
    "npv"
  };

  TString ID[nID+1] = {
    "NoID",
    "TightID",
    "MediumID",
    "HWWID",
    "TightIDGoT",
    "TightIDipsHWW",
    "MediumIDipsHWW"
  };

  for (unsigned int iVar = 0; iVar < nVar; iVar++) {
    for (unsigned int iID = 0; iID < nID+1; iID++) {
      
	h[iVar][iID] = (TH1F*)f->Get("h_Eff_"+Var[iVar]+"_"+ID[iID]+"_AllMu");
      
    }
  }
  
  
  TGraphAsymmErrors* r[nVar][nID];
  TH1F*              h1[nVar][nID+1];

  double ptbins1[24] = {10,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,100,
  		       110,120,130,140,150,170,200};
  double ptbins2[12] = {10,20,25,30,35,40,50,60,75,100,150,200};
  double ptbins3[6]  = {10,20,40,60,100,200};
		       
  double etabins[14] = {-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,-0.2,0.2,0.3,0.9,1.2,1.6,2.1,2.4};

  double npvbins[11] = {0,10,12,14,16,18,20,22,25,30,40};

  
  for (unsigned int iID = 0; iID < nID+1; iID++) {
    
    h1[0][iID] = (TH1F*) h[0][iID]->Rebin(11, "h1", ptbins2);
    h1[1][iID] = (TH1F*) h[1][iID]->Rebin(13, "h1", etabins);
    h1[2][iID] = (TH1F*) h[2][iID]->Rebin(10, "h1", npvbins);
    
  }
  
  for (unsigned int iVar = 0; iVar < nVar; iVar++) {
    for (unsigned int iID = 0; iID < nID; iID++)   {
      
      r[iVar][iID] = new TGraphAsymmErrors(h1[iVar][iID+1], h1[iVar][0]); 
      
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

  TCanvas* c[nVar];
  TLegend* l2;
  TLegend* l3;
  TMultiGraph* mg[nVar];


  if (1) {
    
    c[0] = new TCanvas("Efficiencies vs pt",   "Efficiencies vs pt",   750, 750);
    c[1] = new TCanvas("Efficiencies vs eta",  "Efficiencies vs eta",  750, 750);
    c[2] = new TCanvas("Efficiencies vs npv",  "Efficiencies vs npv",  750, 750);
    

    for (unsigned int iVar=0; iVar<nVar; iVar++) {
      
      mg[iVar] = new TMultiGraph();
      
    }
    

    for (unsigned int iVar=0; iVar<nVar; iVar++) {
      
      c[iVar]->Range(-0.8247861,-0.8247861,0.8247861,0.8247861);
      c[iVar]->SetLeftMargin(0.15);
      c[iVar]->SetRightMargin(0.15);
      c[iVar]->SetTopMargin(0.15);
      c[iVar]->SetBottomMargin(0.15);
      //c[iVar]->SetLogy();
      
    }
    

    l2 = new TLegend(0.20, 0.20, 0.43, 0.35);
    l2->SetBorderSize(1);
    l2->SetFillStyle(1001);
    l2->SetTextFont(42);
    l2->AddEntry(r[0][0], "Tight ID",               "lp");
    l2->AddEntry(r[0][1], "Medium ID",              "lp");
    l2->AddEntry(r[0][2], "HWW ID",                 "lp");
    l2->AddEntry(r[0][3], "Tight ID, GLB or TRK",   "lp");
    l2->AddEntry(r[0][4], "Tight ID, HWW IP cuts",  "lp");
    l2->AddEntry(r[0][5], "Medium ID, HWW IP cuts", "lp");


    l3 = new TLegend(0.60, 0.68, 0.83, 0.83); //(0.60, 0.68, 0.83, 0.83)
    l3->SetBorderSize(1);
    l3->SetFillStyle(1001);
    l3->SetTextFont(42);
    l3->AddEntry(r[0][0], "Tight ID",               "lp");
    l3->AddEntry(r[0][1], "Medium ID",              "lp");
    l3->AddEntry(r[0][2], "HWW ID",                 "lp");
    l3->AddEntry(r[0][3], "Tight ID, GLB or TRK",   "lp");
    l3->AddEntry(r[0][4], "Tight ID, HWW IP cuts",  "lp");
    l3->AddEntry(r[0][5], "Medium ID, HWW IP cuts", "lp");


    c[0]->cd();  
    r[0][0]->GetXaxis()->SetTitle("p_{T}");
    
    c[1]->cd();  
    r[1][0]->GetXaxis()->SetTitle("#eta");
    
    c[2]->cd();  
    r[2][0]->GetXaxis()->SetTitle("NPV");
     
    for (unsigned int iVar=0; iVar<nVar; iVar++)   {
	
      c[iVar]->cd();
	
      r[iVar][0]->SetLineColor(4);
      r[iVar][0]->SetLineWidth(2);
      r[iVar][0]->SetMarkerColor(4);
      r[iVar][0]->SetMarkerSize(0.8);
      r[iVar][0]->SetMarkerStyle(kFullSquare);
      r[iVar][0]->SetFillStyle(3001);
      r[iVar][0]->SetFillColor(4);	   
      r[iVar][1]->SetLineColor(6);
      r[iVar][1]->SetLineWidth(2);
      r[iVar][1]->SetMarkerColor(6);
      r[iVar][1]->SetMarkerSize(0.8);
      r[iVar][1]->SetMarkerStyle(kOpenSquare);
      r[iVar][1]->SetFillStyle(3001);
      r[iVar][1]->SetFillColor(6);
      r[iVar][2]->SetLineColor(8);
      r[iVar][2]->SetLineWidth(2);
      r[iVar][2]->SetMarkerColor(8);
      r[iVar][2]->SetMarkerSize(1.0);
      r[iVar][2]->SetMarkerStyle(kFullTriangleUp);
      r[iVar][2]->SetFillStyle(3001);
      r[iVar][2]->SetFillColor(8);
      r[iVar][3]->SetLineColor(kOrange+7);
      r[iVar][3]->SetLineWidth(2);
      r[iVar][3]->SetMarkerColor(kOrange+7);
      r[iVar][3]->SetMarkerSize(1.0);
      r[iVar][3]->SetMarkerStyle(kOpenCircle);
      r[iVar][3]->SetFillStyle(3001);
      r[iVar][3]->SetFillColor(kOrange+7);
      r[iVar][4]->SetLineColor(kCyan+2);
      r[iVar][4]->SetLineWidth(2);
      r[iVar][4]->SetMarkerColor(kCyan+2);
      r[iVar][4]->SetMarkerSize(0.85);
      r[iVar][4]->SetMarkerStyle(kOpenCircle);
      r[iVar][4]->SetFillStyle(0);
      r[iVar][4]->SetFillColor(kCyan+2);
      r[iVar][5]->SetLineColor(kRed+1);
      r[iVar][5]->SetLineWidth(2);
      r[iVar][5]->SetMarkerColor(kRed+1);
      r[iVar][5]->SetMarkerSize(0.85);
      r[iVar][5]->SetMarkerStyle(kOpenTriangleUp);
      r[iVar][5]->SetFillStyle(0);
      r[iVar][5]->SetFillColor(kRed+1);
	
      mg[iVar]->Add(r[iVar][0], "p");
      mg[iVar]->Add(r[iVar][2], "p");
      mg[iVar]->Add(r[iVar][1], "p");
      mg[iVar]->Add(r[iVar][3], "p");
      mg[iVar]->Add(r[iVar][4], "p");
      mg[iVar]->Add(r[iVar][5], "p"); 

      c[iVar] ->Clear();
      mg[iVar]->Draw("a");
      c[iVar] ->Modified();
      c[iVar] ->Update();

      if (iVar==0) {  
	mg[iVar]->GetXaxis()->SetTitle("p_{T}");
	mg[iVar]->GetXaxis()->SetRangeUser(0.0,200.0); 
      }
      if (iVar==1) {
	mg[iVar]->GetXaxis()->SetTitle("#eta");
	mg[iVar]->GetXaxis()->SetRangeUser(-2.5,2.5); 
      }
      if (iVar==2) {
	mg[iVar]->GetXaxis()->SetTitle("\#PV");
	mg[iVar]->GetXaxis()->SetRangeUser(0.0,45.0); 
      }
	
      mg[iVar]->GetXaxis()->SetLabelOffset(0.015);
      mg[iVar]->GetXaxis()->SetLabelSize(0.025);
      mg[iVar]->GetXaxis()->SetTitleSize(0.035);
      mg[iVar]->GetXaxis()->SetTitleOffset(1.2);
      mg[iVar]->GetYaxis()->SetTitle("Efficiency");
      mg[iVar]->GetYaxis()->SetLabelOffset(0.015);
      mg[iVar]->GetYaxis()->SetLabelSize(0.025);
      mg[iVar]->GetYaxis()->SetTitleSize(0.035);
      mg[iVar]->GetYaxis()->SetTitleOffset(1.6);
	
      if (Sample.Contains("DY")) {
	mg[iVar]->GetYaxis()->SetRangeUser(0.88,1.03);
	l2->Draw(); 
      }
	
      else if (Sample.Contains("QCD")) {
	mg[iVar]->GetYaxis()->SetRangeUser(0.00,1.00);
	l2->Draw(); 
      }
	
      //c[iVar]->Update();
	

    }
    
  }
  
  
}
