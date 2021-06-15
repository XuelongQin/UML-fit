// #include "interface/utils.h"

string year[3] = {"2016","2017","2018"};

std::map<int,std::vector<float>> fM_sigmas =
  {
   {2016, {0.023, 0.015, 0.017, 0.013, 0.005, 0.010, 0.006, 0.013}},
   {2017, {0.018, 0.014, 0.015, 0.010, 0.004, 0.008, 0.005, 0.011}},
   {2018, {0.015, 0.010, 0.011, 0.008, 0.006, 0.006, 0.006, 0.008}},
  };


void plotSimFitComparison_manual(int q2Bin = 6)
{
  string shortString = Form("b%ip1",q2Bin);

  static const int nPars = 8;
  string parName[nPars] = {"Fl","P1","P2","P3","P4p","P5p","P6p","P8p"};

  static const int nMassPars = 3;
  string massParName[nMassPars] = {"mean_{RT}^{%s}",
				   "fsig_"+shortString+"_%s",
				   "slope^{%s}"};

  static const int nMassConstPars = 14;
  string massConstParName[nMassConstPars] = {"f_{M}^{%s}",
					     "#alpha_{RT1}^{%s}",
					     "#alpha_{RT2}^{%s}",
					     "#sigma_{RT1}^{%s}",
					     "#sigma_{RT2}^{%s}",
					     "f^{RT%s}",
					     "n_{RT1}^{%s}",
					     "n_{RT2}^{%s}",
					     "#alpha_{WT1}^{%s}",
					     "#alpha_{WT2}^{%s}",
					     "#sigma_{WT1}^{%s}",
					     "mean_{WT}^{%s}",
					     "n_{WT1}^{%s}",
					     "n_{WT2}^{%s}"};

  static const int nSPars = 6;
  string sParName[nSPars] = {"Fs","As","A4s","A5s","A7s","A8s"};

  double PDG_Fl  [2] = {0.571,0.463};
  double PDG_Fleh[2] = {0.007,0.028};
  double PDG_Flel[2] = {0.007,0.040};

  double bestParSim[nPars] = {5.5337e-01,-1.5926e-02,-1.2078e-03,2.3698e-01,-9.5575e-01,-6.3722e-03,1.8206e-03,-2.1799e-01};
  double lowParSim [nPars] = {-5.41e-04,-3.14e-03,-9.47e-04,-2.46e-03,-3.25e-03,-1.92e-03,-1.58e-03,-3.95e-03};
  double highParSim[nPars] = {5.41e-04,3.14e-03,9.47e-04,2.46e-03,3.25e-03,1.92e-03,1.58e-03,3.95e-03};
  double bestPar[3][nPars];
  double lowPar [3][nPars];
  double highPar[3][nPars];
  
  double bestSParSim[nSPars] = {1.0204e-01,-4.1464e-01,2.3050e-02,4.9952e-03,-6.4214e-03,-3.7205e-01};
  double lowSParSim [nSPars] = {-1.19e-03,-3.06e-03,-7.29e-03,-2.44e-03,-3.62e-03,-3.62e-03};
  double highSParSim[nSPars] = {1.19e-03,3.06e-03,7.29e-03,2.44e-03,3.62e-03,3.62e-03};
  double bestSPar[3][nSPars];
  double lowSPar [3][nSPars];
  double highSPar[3][nSPars];
  
  double bestMassParSim[3*nMassPars] = {5.2754e+00,5.2758e+00,5.2756e+00,8.0127e-01,7.9899e-01,8.0311e-01,-5.7326e+00,-5.5572e+00,-5.5962e+00};
  double errMassParSim [3*nMassPars] = {7.07e-05,5.92e-05,4.32e-05,1.00e-03,7.92e-04,6.26e-4,3.65e-02,2.75e-02,2.09e-02};
  double bestMassConstParSim[3*nMassConstPars] = {9.9999e-01,1.0000e+00,1.0000e+00,1.6961e+00,1.5924e+00,1.6432e+00,-2.1039e+00,-2.4412e+00,-2.4312e+00,2.6405e-02,2.6021e-02,2.6010e-02,4.9160e-02,4.8716e-02,4.9286e-02,5.3159e-01,5.4499e-01,5.3495e-01,2.3803e+00,3.4712e+00,3.0201e+00,4.0165e+00,4.2124e+00,4.3951e+00,8.9729e-01,8.6277e-01,8.3543e-01,1.1477e+00,1.1141e+00,1.0864e+00,5.2726e-02,5.0570e-02,5.2088e-02,5.2736e+00,5.2743e+00,5.2745e+00,8.5247e+00,1.2334e+01,1.3871e+01,1.3334e+01,1.8758e+01,2.0588e+01};
  double errMassConstParSim [3*nMassConstPars] = {1.49e-04,1.05e-04,1.12e-04,1.34e-02,1.16e-02,1.14e-02,1.67e-02,2.36e-02,2.58e-02,8.49e-05,7.32e-05,6.51e-05,1.67e-04,1.38e-04,1.26e-04,3.70e-03,3.27e-03,2.91e-03,4.99e-02,7.81e-02,8.15e-02,1.19e-01,1.57e-01,1.42e-01,1.26e-02,1.08e-02,9.94e-03,1.62e-02,1.41e-02,1.20e-02,3.49e-04,3.12e-04,3.11e-04,2.17e-04,2.01e-04,1.97e-04,6.20e-01,9.61e-01,1.35e+00,8.26e-01,1.18e+00,1.06e+00};
  double bestMassPar[3*nMassPars];
  double errMassPar [3*nMassPars];
  double bestMassConstPar[3*nMassConstPars];
  double errMassConstPar [3*nMassConstPars];

  double massConstVal[3*nMassConstPars];
  double massConstErr[3*nMassConstPars];

  string finName = "simFitResults4d/simFitResult_data_fullAngularMass_Swave_%s_MCStat_b%i_fullStat_noMeanCon.root";
  string finName2 = "/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%s_fM_";
  if (q2Bin==4) finName2 = finName2 + "Jpsi_newbdt.root";
  else if (q2Bin==6) finName2 = finName2 + "Psi_newbdt.root";
  else finName2 = finName2 + "newbdt.root";

  for (int iy=0; iy<3; ++iy) {
    auto fin = TFile::Open(Form(finName.c_str(),year[iy].c_str(),q2Bin));
    if ( !fin || fin->IsZombie() ) {
      cout<<year[iy]<<" file is problematic!"<<endl;
      return;
    }
    auto fin2 = TFile::Open(Form(finName2.c_str(),year[iy].c_str()));
    if ( !fin2 || fin2->IsZombie() ) {
      cout<<year[iy]<<" file2 is problematic!"<<endl;
      return;
    }

    auto tin= (TTree*)fin->Get("fitResultsTree");
    if ( !tin || tin->IsZombie() || tin->GetEntries() < 1 ) {
      cout<<year[iy]<<" tree is problematic!"<<endl;
      fin->Close();
      return;
    }
    for (int iPar=0; iPar<nPars; ++iPar) {
      tin->SetBranchAddress((parName[iPar]+"_best").c_str(),&bestPar[iy][iPar]);
      tin->SetBranchAddress((parName[iPar]+"_low" ).c_str(),&lowPar [iy][iPar]);
      tin->SetBranchAddress((parName[iPar]+"_high").c_str(),&highPar[iy][iPar]);
    }
    tin->GetEntry(0);

    auto ws = (RooWorkspace*)fin->Get("wsp_out");
    if ( !ws || ws->IsZombie() ) {
      cout<<year[iy]<<" workspace is problematic!"<<endl;
      return;
    }
    auto ws2 = (RooWorkspace*)fin2->Get("w");
    if ( !ws2 || ws2->IsZombie() ) {
      cout<<year[iy]<<" workspace2 is problematic!"<<endl;
      return;
    }

    for (int iPar=0; iPar<nSPars; ++iPar) {
      auto par = (RooRealVar*)ws->var(sParName[iPar].c_str());
      bestSPar[iy][iPar] = par->getValV();
      lowSPar [iy][iPar] = par->getErrorLo();
      highSPar[iy][iPar] = par->getErrorHi();
    }
    for (int iPar=0; iPar<nMassPars; ++iPar) {
      auto par = (RooRealVar*)ws->var(Form(massParName[iPar].c_str(),year[iy].c_str()));
      bestMassPar[iPar*3+iy] = par->getValV();
      errMassPar [iPar*3+iy] = par->getError();
    }
    ws2->loadSnapshot(Form("reference_fit_RT_%i",q2Bin));
    for (int iPar=0; iPar<nMassConstPars; ++iPar) {
      auto par = (RooRealVar*)ws->var(Form(massConstParName[iPar].c_str(),year[iy].c_str()));
      bestMassConstPar[iPar*3+iy] = par->getValV();
      errMassConstPar [iPar*3+iy] = par->getError();
      if (iPar==8) ws2->loadSnapshot(Form("reference_fit_WT_%i",q2Bin));
      if (iPar>0) {
	auto par2 = (RooRealVar*)ws2->var(Form(massConstParName[iPar].c_str(),Form("%i",q2Bin)));
	if ( !par2 || par2->IsZombie() ) {
	  cout<<Form(massConstParName[iPar].c_str(),year[iy].c_str())<<" constraint not found"<<endl;
	  return;
	}
	massConstVal[iPar*3+iy] = par2->getValV();
	massConstErr[iPar*3+iy] = par2->getError();
      } else {
	double nrt_mc = ws2->var(Form("nRT_%i",q2Bin))->getVal();
	double nwt_mc = ws2->var(Form("nWT_%i",q2Bin))->getVal();
	massConstVal[iy] = 1.;
	massConstErr[iy] = fM_sigmas[atoi(year[iy].c_str())][q2Bin] / ( nwt_mc / (nrt_mc + nwt_mc) );
	cout<<"==========="
	    <<iy<<"\t"<<nrt_mc<<"\t"<<nwt_mc<<"\t"
	    <<massConstErr[iy]<<"\t"
	    <<bestMassConstPar[iy]<<"\t"
	    <<errMassConstPar [iPar*3+iy]<<endl;
      }
    }

    fin->Close();
    fin2->Close();

  }

  // Plot P-wave parameters
  double x[3][nPars];
  double xe[nPars];
  double xSim[nPars];
  double xeSim[nPars];

  double ySim [nPars];
  double ylSim[nPars];
  double yhSim[nPars];
  double y [3][nPars];
  double yl[3][nPars];
  double yh[3][nPars];

  for (int iPar=0; iPar<nPars; ++iPar) {
    xSim[iPar]  = 0.5+iPar;
    xeSim[iPar] = 0.5;
    xe[iPar]    = 1.0/6;
    ySim [iPar] = 0;
    ylSim[iPar] = -1 * lowParSim [iPar];
    yhSim[iPar] = highParSim[iPar];
    for (int iy=0; iy<3; ++iy) {
      x [iy][iPar] = (2.0*iy+1)/6+iPar;
      y [iy][iPar] = bestPar[iy][iPar] - bestParSim[iPar];
      yl[iy][iPar] = -1 * lowPar [iy][iPar];
      yh[iy][iPar] = highPar[iy][iPar];
    }
  }

  TCanvas canv("canv","canv",1500,1000);
  canv.cd();

  gStyle->SetOptStat(0);
  // cCover.cd()->SetTicky(2);
  string plotTitle = Form("Comparison of simultaneous fit results - q2 bin %i;;(single-year result) - (simultaneous result)",q2Bin);
  auto hLab = new TH1S ("hLab",plotTitle.c_str(),nPars,0,nPars);
  for (int iBin=0; iBin<nPars; ++iBin)
    hLab->GetXaxis()->SetBinLabel(iBin+1,parName[iBin].c_str());
  double yRange = 0.1;
  if (q2Bin==4) yRange = 0.05;
  hLab->SetMinimum(-1*yRange);
  hLab->SetMaximum(yRange);
  // hLab->GetYaxis()->SetTickLength(0.006);
  hLab->Draw();

  auto grSim = new TGraphAsymmErrors(nPars,xSim,ySim,xeSim,xeSim,ylSim,yhSim);
  grSim->SetName("grSim");
  grSim->SetFillColor(30);
  grSim->Draw("2");

  TGraphAsymmErrors* grPDG;
  double yPDG[1];
  double ylPDG[1];
  double yhPDG[1];
  if (q2Bin==4 || q2Bin==6) {
    uint iCR = (q2Bin-4) / 2;
    yPDG [0] = PDG_Fl[iCR]-bestParSim[0];
    ylPDG[0] = PDG_Flel[iCR];
    yhPDG[0] = PDG_Fleh[iCR];
    grPDG = new TGraphAsymmErrors(1,xSim,yPDG,xeSim,xeSim,PDG_Flel,PDG_Fleh);
    grPDG->SetName("grPDG");
    grPDG->SetLineColor(15);
    grPDG->SetLineWidth(2);
    grPDG->Draw("P");
  }

  TLine line (0,0,nPars,0);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.SetLineColor(13);
  line.Draw();

  TGraphAsymmErrors* gr [3];
  int colours [3] = {2,kViolet+7,4};
  for (int iy=0; iy<3; ++iy) {
    gr[iy] = new TGraphAsymmErrors(nPars,x[iy],y[iy],xe,xe,yl[iy],yh[iy]);
    gr[iy]->SetName(Form("gr%i",iy));
    gr[iy]->SetLineColor(colours[iy]);
    gr[iy]->SetLineWidth(2);
    gr[iy]->Draw("P");
  }

  TLegend leg (0.15,0.66,0.4,0.9);
  leg.AddEntry(grSim,"Simultaneous result","f");
  if (q2Bin==4 || q2Bin==6) leg.AddEntry(grPDG,"PDG value","lep");
  for (int iy=0; iy<3; ++iy)
    leg.AddEntry(gr[iy],("Result "+year[iy]).c_str(),"lep");
  leg.Draw();

  canv.SaveAs(("plotSimFit4d_d/comparisonSimFit_result_"+shortString+".pdf").c_str());

  // Plot S-wave parameters
  double xS[3][nSPars];
  double xSe[nSPars];
  double xSSim[nSPars];
  double xSeSim[nSPars];

  double ySSim [nSPars];
  double ySlSim[nSPars];
  double yShSim[nSPars];
  double yS [3][nSPars];
  double ySl[3][nSPars];
  double ySh[3][nSPars];

  for (int iPar=0; iPar<nSPars; ++iPar) {
    xSSim[iPar]  = 0.5+iPar;
    xSeSim[iPar] = 0.5;
    xSe[iPar]    = 1.0/6;
    ySSim [iPar] = 0;
    ySlSim[iPar] = -1 * lowSParSim [iPar];
    yShSim[iPar] = highSParSim[iPar];
    for (int iy=0; iy<3; ++iy) {
      xS [iy][iPar] = (2.0*iy+1)/6+iPar;
      yS [iy][iPar] = bestSPar[iy][iPar] - bestSParSim[iPar];
      ySl[iy][iPar] = -1 * lowSPar [iy][iPar];
      ySh[iy][iPar] = highSPar[iy][iPar];
    }
  }

  TCanvas canvS("canvS","canvS",1500,1000);
  canvS.cd();

  string plotTitleS = Form("Comparison of simultaneous fit S-wave parameters - q2 bin %i;;(single-year result) - (simultaneous result)",q2Bin);
  auto hLabS = new TH1S ("hLabS",plotTitleS.c_str(),nSPars,0,nSPars);
  for (int iBin=0; iBin<nSPars; ++iBin)
    hLabS->GetXaxis()->SetBinLabel(iBin+1,sParName[iBin].c_str());
  double yRangeS = 0.2;
  if (q2Bin==4) yRangeS = 0.05;
  hLabS->SetMinimum(-1*yRangeS);
  hLabS->SetMaximum(yRangeS);
  hLabS->Draw();

  auto grSimS = new TGraphAsymmErrors(nSPars,xSSim,ySSim,xSeSim,xSeSim,ySlSim,yShSim);
  grSimS->SetName("grSimS");
  grSimS->SetFillColor(30);
  grSimS->Draw("2");

  TLine lineS (0,0,nSPars,0);
  lineS.SetLineWidth(2);
  lineS.SetLineStyle(2);
  lineS.SetLineColor(13);
  lineS.Draw();

  TGraphAsymmErrors* grS [3];
  for (int iy=0; iy<3; ++iy) {
    grS[iy] = new TGraphAsymmErrors(nSPars,xS[iy],yS[iy],xSe,xSe,ySl[iy],ySh[iy]);
    grS[iy]->SetName(Form("grS%i",iy));
    grS[iy]->SetLineColor(colours[iy]);
    grS[iy]->SetLineWidth(2);
    grS[iy]->Draw("P");
  }

  TLegend legS (0.65,0.73,0.9,0.9);
  legS.AddEntry(grSimS,"Simultaneous result","f");
  for (int iy=0; iy<3; ++iy)
    legS.AddEntry(grS[iy],("Result "+year[iy]).c_str(),"lep");
  legS.Draw();

  canvS.SaveAs(("plotSimFit4d_d/comparisonSimFit_Swave_"+shortString+".pdf").c_str());

  // Plot mass parameters
  double xMass[3*(nMassPars+nMassConstPars)];
  double xMasse[3*(nMassPars+nMassConstPars)];

  double yMassSim [3*(nMassPars+nMassConstPars)];
  double yMasseSim[3*(nMassPars+nMassConstPars)];
  double yMass [3*(nMassPars+nMassConstPars)];
  double yMasse[3*(nMassPars+nMassConstPars)];
  double yMassConst [3*(nMassPars+nMassConstPars)];
  double yMasseConst[3*(nMassPars+nMassConstPars)];

  for (int iPar=0; iPar<3*(nMassPars+nMassConstPars); ++iPar) {
    xMass[iPar]  = 0.5+iPar;
    xMasse[iPar] = 0.5;
    yMassSim [iPar] = 0.;
    yMasseSim[iPar] = 1.;
  }
  for (int iPar=0; iPar<3*nMassPars; ++iPar) {
    yMass [iPar] = (bestMassPar[iPar] - bestMassParSim[iPar]) / errMassParSim[iPar];
    yMasse[iPar] = errMassPar[iPar] / errMassParSim[iPar];
    yMassConst [iPar] = 1e99;
    yMasseConst[iPar] = 0;
  }
  for (int iPar=0; iPar<3*nMassConstPars; ++iPar) {
    yMass [iPar+3*nMassPars] = (bestMassConstPar[iPar] - bestMassConstParSim[iPar]) / errMassConstParSim[iPar];
    yMasse[iPar+3*nMassPars] = errMassConstPar[iPar] / errMassConstParSim[iPar];
    yMassConst [iPar+3*nMassPars] = (massConstVal[iPar] - bestMassConstParSim[iPar]) / errMassConstParSim[iPar];
    yMasseConst[iPar+3*nMassPars] = massConstErr[iPar] / errMassConstParSim[iPar];
  }

  TCanvas canvM("canvM","canvM",1500,1000);
  canvM.cd();

  string plotTitleM = Form("Comparison of simultaneous fit mass parameters - q2 bin %i;;pulls",q2Bin);
  auto hLabM = new TH1S ("hLabM",plotTitleM.c_str(),3*(nMassPars+nMassConstPars),0,3*(nMassPars+nMassConstPars));
  for (int iBin=0; iBin<3*nMassPars; ++iBin)
    hLabM->GetXaxis()->SetBinLabel(iBin+1,Form(massParName[iBin/3].c_str(),year[iBin%3].c_str()));
  for (int iBin=0; iBin<3*nMassConstPars; ++iBin)
    hLabM->GetXaxis()->SetBinLabel(iBin+1+3*nMassPars,Form(massConstParName[iBin/3].c_str(),year[iBin%3].c_str()));
  double yRangeM = 9;
  if (q2Bin==4) yRangeM = 30;
  hLabM->SetMinimum(-1*yRangeM);
  hLabM->SetMaximum(yRangeM);
  hLabM->Draw();

  auto grSimM = new TGraphErrors(3*(nMassPars+nMassConstPars),xMass,yMassSim,xMasse,yMasseSim);
  grSimM->SetName("grSimM");
  grSimM->SetFillColor(30);
  grSimM->Draw("2");

  TLine lineM (0,0,3*(nMassPars+nMassConstPars),0);
  lineM.SetLineWidth(2);
  lineM.SetLineStyle(2);
  lineM.SetLineColor(13);
  lineM.Draw();

  auto grM = new TGraphErrors(3*(nMassPars+nMassConstPars),xMass,yMass,xMasse,yMasse);
  grM->SetName("grM");
  grM->SetLineColor(2);
  grM->SetLineWidth(2);
  grM->Draw("P");

  auto grConst = new TGraphErrors(3*(nMassPars+nMassConstPars),xMass,yMassConst,xMasse,yMasseConst);
  grConst->SetName("grConst");
  grConst->SetLineColor(1);\
  grConst->SetLineWidth(2);
  grConst->SetFillStyle(0);
  grConst->SetMarkerStyle(20);
  grConst->Draw("5p");

  TLegend legM (0.1,0.73,0.35,0.9);
  legM.AddEntry(grSimM,"Simultaneous result","f");
  legM.AddEntry(grM,"Single-year result","lep");
  legM.AddEntry(grConst,"Constraint","fp");
  legM.Draw();

  canvM.SaveAs(("plotSimFit4d_d/comparisonSimFit_massParameters_"+shortString+".pdf").c_str());

}
