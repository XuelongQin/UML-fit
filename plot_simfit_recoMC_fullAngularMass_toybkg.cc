#include <TFile.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <list>
#include <map>

// #include <RooRealVar.h>
#include <RooAbsPdf.h>
// #include <RooWorkspace.h>
#include <RooCategory.h>
#include <RooSuperCategory.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooMinimizer.h>
#include <RooPlot.h>
#include <RooHistFunc.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooNumIntConfig.h>
#include <RooAddition.h>
#include <RooRandom.h>
// #include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooCBShape.h>
#include "RooDoubleCBFast.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"

#include "utils.h"
#include "PdfSigRTMass.h"
#include "PdfSigWTMass.h"
#include "PdfSigAngMass.h"
#include "ShapeSigAng.h"

#include "BoundCheck.h"
#include "BoundDist.h"
#include "Penalty.h"
#include "Fitter.h"
#include "RooBernsteinSideband.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;
std::map<int,float> scale_to_data;

TCanvas* c [4*nBins];

double power = 1.0;

static const int nRanges = 5;
double rangeVal[nRanges+1] = {5.00,5.10,5.20,5.35,5.45,5.60};
string rangeName[nRanges] = {"lsb","ltail","peak","rtail","rsb"};

void plot_simfit_recoMC_fullAngularMass_toybkgBin(int q2Bin, int parity, uint nSample, std::vector<int> years)
{

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  map< string, vector<int>* > rangeComp;
  for (int i=0; i<nRanges; ++i)
    rangeComp[rangeName[i]] = new vector<int> {i};
  rangeComp["sb"]        = new vector<int> {0,4};
  rangeComp["largepeak"] = new vector<int> {1,2,3};
  rangeComp["largesb"]   = new vector<int> {0,1,3,4};

  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  string all_years = "";
  string year = ""; 

  RooCategory sample ("sample", "sample");
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    all_years += year;
    sample.defineType(("data"+year+Form("_subs%i",nSample)).c_str());
  }

  string fileDir = "/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/";
  string fileName = fileDir + "simFitResults4d/simFitResult_recoMC_fullAngularMass_toybkg" + all_years + Form("_dataStat-%i_b%i.root", nSample, q2Bin);
  TFile* fin = new TFile(fileName.c_str(),"READ");
  auto wsp_out = (RooWorkspace*)fin->Get(("ws_"+shortString+Form("_s%i_pow1.0",nSample+1)).c_str());
  auto combData = (RooDataSet*)wsp_out->data("allcombData");
  auto simPdf = (RooSimultaneous*)wsp_out->pdf("simPdf");
  auto mass = (RooRealVar*)wsp_out->var("mass");
  auto ctK = (RooRealVar*)wsp_out->var("ctK");
  auto ctL = (RooRealVar*)wsp_out->var("ctL");
  auto phi = (RooRealVar*)wsp_out->var("phi");

  string plotString = shortString + "_" + all_years + Form("_s%i",nSample);
   
  int confIndex = 2*nBins*parity  + q2Bin;
  string longString  = "Fit to reconstructed events";
  longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);

  // plot fit projections 
  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to RECO-level MC - "+longString).c_str(),2000,1500);
  if (years.size()>1) c[confIndex]->Divide(4, years.size());
  else c[confIndex]->Divide(2,2);
  
  cout<<"plotting 4d canvas"<<endl;
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
  
    std::vector<RooPlot*> frames;
    frames.push_back( prepareFrame( mass ->frame(Title((longString+year).c_str()))));
    frames.push_back( prepareFrame( ctK ->frame(Title((longString+year).c_str())) ));
    frames.push_back( prepareFrame( ctL ->frame(Title((longString+year).c_str())) ));
    frames.push_back( prepareFrame( phi ->frame(Title((longString+year).c_str())) ));
    TLegend* leg = new TLegend (0.60,0.75,0.9,0.9);

    auto singleYearPdf = simPdf->getPdf(("data"+year+Form("_subs%d",nSample)).c_str());
    auto singleYearData = combData->reduce(("sample==sample::data"+year+Form("_subs%d",nSample)).c_str());

    cout<<"canvas ready"<<endl;
    for (unsigned int fr = 0; fr < frames.size(); fr++){
      cout<<"fr "<<fr<<" data"<<year<<Form("_subs%d",nSample)<<endl;
        singleYearData->plotOn(frames[fr], MarkerColor(kRed+1), LineColor(kRed+1), Binning(40), Name(("plData"+year).c_str()));
        // combData->plotOn(frames[fr], MarkerColor(kRed+1), LineColor(kRed+1), Binning(40), Cut(("sample==sample::data"+year+Form("_subs%d",nSample)).c_str()), Name(("plData"+year).c_str()));
        
        singleYearPdf->plotOn(frames[fr],
		       // Slice(sample, ("data"+year+Form("_subs%d",nSample)).c_str()), 
		       // ProjWData(RooArgSet(sample), *combData), 
		       LineWidth(1), 
		       Name(("plPDF"+year).c_str()), 
		       NumCPU(4));

        singleYearPdf->plotOn(frames[fr],
		       // Slice(sample, ("data"+year+Form("_subs%d",nSample)).c_str()), 
		       // ProjWData(RooArgSet(sample), *combData), 
		       LineWidth(1), 
		       Name(("plPDFbkg"+year).c_str()), 
		       NumCPU(4),
		       LineColor(8),
		       Components( ("bkg_pdf_"+year).c_str() ));

        singleYearPdf->plotOn(frames[fr],
		       // Slice(sample, ("data"+year+Form("_subs%d",nSample)).c_str()), 
		       // ProjWData(RooArgSet(sample), *combData), 
		       LineWidth(1), 
		       Name(("plPDFsig"+year).c_str()), 
		       NumCPU(4),
		       LineColor(880),
		       Components( ("PDF_sig_ang_mass_"+shortString+"_"+year).c_str() ));

        if (fr == 0) { 
          leg->AddEntry(frames[fr]->findObject(("plData"+year).c_str()),
			"Data",	"lep");
          leg->AddEntry(frames[fr]->findObject(("plPDF"+year ).c_str()),
			"Total PDF","l");
          leg->AddEntry(frames[fr]->findObject(("plPDFsig"+year ).c_str()),
			"Signal","l");
          leg->AddEntry(frames[fr]->findObject(("plPDFbkg"+year ).c_str()),
			"Background","l");
        }

        c[confIndex]->cd(iy*4+fr+1);
        gPad->SetLeftMargin(0.19); 
        frames[fr]->Draw();
        if (fr == 0) leg->Draw("same");
    }
  }

  c[confIndex]->SaveAs( ("plotSimFit4d_d/simFitResult_recoMC_fullAngularMass_toybkg_" + plotString +  ".pdf").c_str() );

  for (int i=0; i<nRanges; ++i)
    mass->setRange(rangeName[i].c_str(),rangeVal[i],rangeVal[i+1]);

  // plot fit projections in left sideband region 
  for (const auto& [name, comp] : rangeComp) {
    auto canv = new TCanvas (Form("canv%s",name.c_str()),
			     Form("canv%s",name.c_str()),
			     1500,500*years.size());
    canv->Divide(3, years.size());

    string massCut = Form("(mass>%.2f && mass<%.2f)",
			  rangeVal[comp->at(0)],
			  rangeVal[comp->at(0)+1]);
    string rangeList = rangeName[comp->at(0)];

    for (uint iCut=1; iCut<comp->size(); ++iCut) {
      massCut = massCut + Form("||(mass>%.2f && mass<%.2f)",
			  rangeVal[comp->at(iCut)],
			  rangeVal[comp->at(iCut)+1]);
      rangeList = rangeList + "," + rangeName[comp->at(iCut)];
    }

    for (unsigned int iy = 0; iy < years.size(); iy++) {
      year.clear(); year.assign(Form("%i",years[iy]));
  
      std::vector<RooPlot*> frames;
      frames.push_back( prepareFrame( ctK ->frame(Title((longString+year).c_str())) ));
      frames.push_back( prepareFrame( ctL ->frame(Title((longString+year).c_str())) ));
      frames.push_back( prepareFrame( phi ->frame(Title((longString+year).c_str())) ));
      TLegend* leg = new TLegend (0.45,0.7,0.9,0.9);

      auto singleYearPdf = simPdf->getPdf(("data"+year+Form("_subs%d",nSample)).c_str());
      auto singleYearData = combData->reduce(("sample==sample::data"+year+Form("_subs%d",nSample)).c_str());

      for (unsigned int fr = 0; fr < frames.size(); fr++){

	singleYearData->plotOn(frames[fr],
			 MarkerColor(kRed+1),
			 LineColor(kRed+1),
			 Binning(40),
			 Cut(massCut.c_str()),
			 // Cut(("("+massCut+")&&sample==sample::data"+year+Form("_subs%d",nSample)).c_str()),
			 Name(("plData"+name+year).c_str()));
        
	singleYearPdf->plotOn(frames[fr],
		       // Slice(sample, ("data"+year+Form("_subs%d",nSample)).c_str()), 
		       // ProjWData(RooArgSet(sample), *combData), 
		       ProjectionRange(rangeList.c_str()),
		       LineWidth(1), 
		       Name(("plPDF"+name+year).c_str()),
		       NumCPU(4));

	singleYearPdf->plotOn(frames[fr],
		       // Slice(sample, ("data"+year+Form("_subs%d",nSample)).c_str()), 
		       // ProjWData(RooArgSet(sample), *combData), 
		       ProjectionRange(rangeList.c_str()),
		       LineWidth(1), 
		       LineColor(8),
		       Name(("plPDFbkg"+name+year).c_str()),
		       Components( ("bkg_pdf_"+year).c_str() ),
		       NumCPU(4));

	singleYearPdf->plotOn(frames[fr],
		       // Slice(sample, ("data"+year+Form("_subs%d",nSample)).c_str()), 
		       // ProjWData(RooArgSet(sample), *combData), 
		       ProjectionRange(rangeList.c_str()),
		       LineWidth(1), 
		       LineColor(880),
		       Name(("plPDFsig"+name+year).c_str()),
		       Components( ("PDF_sig_ang_mass_"+shortString+"_"+year).c_str() ),
		       NumCPU(4));

	if (fr == 0) { 
          leg->AddEntry(frames[fr]->findObject(("plData"+name+year).c_str()),
			"Data",	"lep");
          leg->AddEntry(frames[fr]->findObject(("plPDF"+name+year ).c_str()),
			"Total PDF","l");
          leg->AddEntry(frames[fr]->findObject(("plPDFsig"+name+year ).c_str()),
			"Signal","l");
          leg->AddEntry(frames[fr]->findObject(("plPDFbkg"+name+year ).c_str()),
			"Background","l");
        }
        canv->cd(iy*3+fr+1);
        gPad->SetLeftMargin(0.16); 
        frames[fr]->Draw();
        if (fr == 0) leg->Draw("same");
      }
    }

    canv->SaveAs( ("plotSimFit4d_d/simFitResult_recoMC_fullAngularMass_toybkg_" + plotString + "_" + name + ".pdf").c_str() );
    
  }

  return;

}



void plot_simfit_recoMC_fullAngularMass_toybkgBin1(int q2Bin, int parity, uint nSample, std::vector<int> years)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      plot_simfit_recoMC_fullAngularMass_toybkgBin(q2Bin, parity, nSample, years);
  else
    plot_simfit_recoMC_fullAngularMass_toybkgBin(q2Bin, parity, nSample, years);
}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency
  //                [1] odd efficiency
  //                [-1] for each parity recursively

  int q2Bin   = -1;
  int parity  = -1; 

  if ( argc > 1 ) q2Bin   = atoi(argv[1]);
  if ( argc > 2 ) parity  = atoi(argv[2]);

  uint nSample = 0;
  if ( argc > 3 ) nSample = atoi(argv[3]);

  std::vector<int> years;
  if ( argc > 4 && atoi(argv[4]) != 0 ) years.push_back(atoi(argv[4]));
  else {
    cout << "No specific years selected, using default: 2016" << endl;
    years.push_back(2016);
  }
  if ( argc > 5 && atoi(argv[5]) != 0 ) years.push_back(atoi(argv[5]));
  if ( argc > 6 && atoi(argv[6]) != 0 ) years.push_back(atoi(argv[6]));

  cout <<  "q2Bin       " << q2Bin        << endl;
  cout <<  "parity      " << parity       << endl;
  cout <<  "nSample     " << nSample      << endl;
  cout <<  "years[0]    " << years[0]     << endl;
//   cout << years[1] << endl;
//   cout << years[2] << endl;


  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;

  if ( q2Bin==-1 )   cout << "Running all the q2 bins" << endl;
  if ( parity==-1 )  cout << "Running both the parity datasets" << endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      plot_simfit_recoMC_fullAngularMass_toybkgBin1(q2Bin, parity, nSample, years);
  else
    plot_simfit_recoMC_fullAngularMass_toybkgBin1(q2Bin, parity, nSample, years);

  return 0;

}
