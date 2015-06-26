void ISORocCurve(TString BaseID = "TightID", TString dR = "04")
{

  TFile *f  = TFile::Open("../files/DR74X_50ns_MC_DY.root");
  TFile *g  = TFile::Open("../files/DR74X_50ns_MC_QCD.root");

  const int nISO = 4;

  TH1F* h[2][nISO];

  TString ISO[nISO] = {
    "",
    "dBeta",
    "PFWeighted",
    "PUPPI"
  };

  for (unsigned int i = 0; i < nISO; i++) {
    
    h[0][i] = (TH1F*)f->Get("h_RC_"+BaseID+"_ISO"+dR+ISO[i]+"_AllMu");
    h[1][i] = (TH1F*)g->Get("h_RC_"+BaseID+"_ISO"+dR+ISO[i]+"_AllMu");
    
  }

  double Eff[2][nISO][25];
  double Err[2][nISO][25];

  for (unsigned int i = 0; i < 2; i++)        {
    for (unsigned int j = 0; j < nISO; j++)   {
      for (unsigned int k = 1; k <= 25; k++)  {
	
	Eff[i][j][k-1] = h[i][j]->GetBinContent(k)/h[i][j]->GetBinContent(0);
	Err[i][j][k-1] = 
	  (h[i][j]->GetBinError(k) + h[i][j]->GetBinError(0)*h[i][j]->GetBinContent(k)/h[i][j]->GetBinContent(0))/
	  h[i][j]->GetBinContent(0);
	
      }
    }
  }

  TGraphErrors* ROC[nISO];

  for (unsigned int i = 0; i < nISO; i++) {
    ROC[i] = new TGraphErrors(25, Eff[1][i], Eff[0][i], Err[1][i], Err[0][i]);
    //ROC[i] = new TGraphErrors(25, Eff[1][i], Eff[0][i]);
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

  TCanvas* c;
  TLegend* l;
  TMultiGraph* mg;

  if (1) {

    c = new TCanvas("Efficiency ROC Curve",   "Efficiency ROC Curve",   750, 750);

    c->Range(-0.8247861,-0.8247861,0.8247861,0.8247861);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.15);
    c->SetTopMargin(0.15);
    c->SetBottomMargin(0.15);

    mg = new TMultiGraph();

    for (unsigned int i = 0; i < nISO; i++) {
      mg->Add(ROC[i]);      
    }

    TLegend* l = new TLegend(0.60, 0.25, 0.83, 0.43);
    l->SetBorderSize(0);
    l->SetFillStyle(1001);
    l->SetTextFont(42);
    l->AddEntry(ROC[0], "Standard",              "lp");
    l->AddEntry(ROC[1], "#Delta#beta corrected", "lp");
    l->AddEntry(ROC[2], "PF-Weighted",           "lp");
    l->AddEntry(ROC[3], "PUPPI",                 "lp");

    ROC[0]->SetLineColor(4);
    ROC[0]->SetLineWidth(2);
    ROC[0]->SetMarkerColor(4);
    ROC[0]->SetMarkerSize(0.8);
    ROC[0]->SetMarkerStyle(kFullTriangleUp);
    ROC[0]->SetFillStyle(3001);
    ROC[0]->SetFillColor(4);	   
    ROC[1]->SetLineColor(6);
    ROC[1]->SetLineWidth(2);
    ROC[1]->SetMarkerColor(6);
    ROC[1]->SetMarkerSize(0.8);
    ROC[1]->SetMarkerStyle(kFullSquare);
    ROC[1]->SetFillStyle(3001);
    ROC[1]->SetFillColor(6);
    ROC[2]->SetLineColor(8);
    ROC[2]->SetLineWidth(2);
    ROC[2]->SetMarkerColor(8);
    ROC[2]->SetMarkerSize(1.0);
    ROC[2]->SetMarkerStyle(kOpenSquare);
    ROC[2]->SetFillStyle(3001);
    ROC[2]->SetFillColor(8);
    ROC[3]->SetLineColor(kOrange+7);
    ROC[3]->SetLineWidth(2);
    ROC[3]->SetMarkerColor(kOrange+7);
    ROC[3]->SetMarkerSize(1.0);
    ROC[3]->SetMarkerStyle(kOpenCircle);
    ROC[3]->SetFillStyle(3001);
    ROC[3]->SetFillColor(kOrange+7);

    mg->Draw("alp");
    l->Draw();

    mg->GetXaxis()->SetTitle("QCD Efficiency");
    mg->GetXaxis()->SetLabelOffset(0.015);
    mg->GetXaxis()->SetLabelSize(0.025);
    mg->GetXaxis()->SetTitleSize(0.035);
    mg->GetXaxis()->SetTitleOffset(1.2);
    mg->GetYaxis()->SetTitle("DY Efficiency");
    mg->GetYaxis()->SetLabelOffset(0.015);
    mg->GetYaxis()->SetLabelSize(0.025);
    mg->GetYaxis()->SetTitleSize(0.035);
    mg->GetYaxis()->SetTitleOffset(1.6);

    mg->GetXaxis()->SetRangeUser(0.00,0.30);
    mg->GetYaxis()->SetRangeUser(0.70,1.00);


  }
    
}

  
