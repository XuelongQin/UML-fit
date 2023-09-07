#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>
#include <list>
#include <map>

#include <RooAbsPdf.h>
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
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooCBShape.h>
#include "RooDoubleCBFast.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "utils.h"
#include "PdfSigRTMass.h"
#include "PdfSigWTMass.h"

#include "RooBernsteinSideband.h"

//tmp
#include "RooMCStudy.h"


using namespace RooFit;
using namespace std;

static const int nBins = 9;

double power = 1.0;


void cocktailsample(int q2Bin, int parity, bool multiSample, uint nSample, int XGBv, bool localFiles, bool plot, int save, std::vector<int> years)
{

    string XGBstr = "";
    if (XGBv>9) XGBstr = Form("-TMVAv%i",XGBv-10);
    else if (XGBv>0) XGBstr = Form("_XGBv%i",XGBv);
    string shortString = Form("b%ip%i",q2Bin,parity);
    cout<<"Conf: "<<shortString<<endl;
    string stat = nSample > 0 ? "_dataStat":"_MCStat";
    uint firstSample = ( multiSample || nSample==0 ) ? 0 : nSample-1;
    uint lastSample = nSample > 0 ? nSample-1 : 0;
    RooRealVar* ctK = new RooRealVar("ctK", "cos(#theta_{K})", -1  , 1  );
    RooRealVar* ctL = new RooRealVar("ctL", "cos(#theta_{l})", -1  , 1  );
    RooRealVar* phi = new RooRealVar("phi", "#phi", -3.14159, 3.14159  );
    RooArgList vars (* ctK,* ctL,* phi);
    RooRealVar* rand = new RooRealVar("rand", "rand", 0,1);
    RooRealVar* mass = new RooRealVar("mass","m(#mu#muK#pi)", 5.,5.6,"GeV");
    RooRealVar* wei  = new RooRealVar("weight","weight",1);
    RooArgSet reco_vars (*ctK, *ctL, *phi, *rand, *mass, *wei);
    RooArgSet observables (*ctK, *ctL, *phi, *mass);
    std::vector<RooWorkspace*> wsp, wsp_mcmass, wsp_sb;
    std::vector<std::vector<RooDataSet*>> data;
    string bkgpdfname = "bkg_pdf_%i";
    // Random generators
    TRandom gen_nevt;
    string foutname = Form("/afs/cern.ch/user/x/xuqin/eos/Kmumu/mom_cocktail/recoMCDataset_b%i.root",q2Bin);
    TFile *fout = new TFile(foutname.c_str(),"recreate");
    RooWorkspace *wsp_out = new RooWorkspace(Form("ws_b%ip%i",q2Bin,1-parity),Form("ws_b%ip%i",q2Bin,1-parity));
    for (unsigned int iy = 0; iy < years.size(); iy++) {
        string filename_data = Form("recoMCDataset_b%i_%i%s.root", q2Bin, years[iy], XGBstr.c_str());
        if (!localFiles) filename_data = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/MC-datasets/" + filename_data;
        retrieveWorkspace( filename_data, wsp, Form("ws_b%ip%i", q2Bin, 1-parity ));
    data.push_back( createDataset( nSample,  firstSample,  lastSample, wsp[iy],  
                                   q2Bin,  parity,  years[iy], 
                                   observables,  shortString  )); 
        cout<< " create dataset !" << endl;
        // now generate bkg events
        int nbkg_togen = nbkg_years[years[iy]][q2Bin];

        // Read angular pdf for sidebands from external file 
        string filename_sb = Form("/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/sidebands/savesb_%i_b%i.root",years[iy],q2Bin);
        retrieveWorkspace( filename_sb, wsp_sb, "wsb");

        RooBernsteinSideband* bkg_ang_pdf = (RooBernsteinSideband*) wsp_sb[iy]->pdf(Form("BernSideBand_bin%i_%i", q2Bin, years[iy]));
        //RooArgSet*  bkg_ang_params = (RooArgSet*) bkg_ang_pdf->getParameters(observables);

        // read mass pdf for background
        RooRealVar* slope       = new RooRealVar    (Form("slope^{%i}",years[iy]),  Form("slope^{%i}",years[iy]) , wsp_sb[iy]->var("slope")->getVal(), -10., 0.);
        RooExponential* bkg_exp = new RooExponential(Form("bkg_exp_%i",years[iy]),  Form("bkg_exp_%i",years[iy]) ,  *slope,   *mass  );
        cout << Form("exponential slope for %f", slope->getVal())  << endl;

    // retrieve sideband range from input file
        float max_lsb = wsp_sb[iy]->var(Form("max_sbl_bin%i_%i", q2Bin, years[iy]))->getVal();
        float min_rsb = wsp_sb[iy]->var(Form("min_sbr_bin%i_%i", q2Bin, years[iy]))->getVal();
        mass->setRange("sbleft",  5., max_lsb);
        mass->setRange("sbright", min_rsb, 5.6);
        cout << Form("sideband mass range: [5., %.2f] U [%.2f, 5.6]", max_lsb, min_rsb)  << endl;
        //wei->setRange("sbw",1.0,1.0);
        //RooPolynomial *bkg_weight = new RooPolynomial(Form("bkg_weight_%i",years[iy]),  Form("bkg_weight_%i",years[iy]) , *wei);
        // create 4D pdf  for background and import to workspace
        RooProdPdf* bkg_pdf = new RooProdPdf(Form(bkgpdfname.c_str(),years[iy]), Form(bkgpdfname.c_str(),years[iy]),
                        RooArgList(*bkg_ang_pdf,*bkg_exp)); 

        RooAbsPdf::GenSpec* genSpec = bkg_pdf->prepareMultiGen( observables, NumEvents(gen_nevt.Poisson(nbkg_togen)));

        // now generate toy bkg sample
        for (uint itoy = 0; itoy <= lastSample-firstSample; itoy++){
            // set the random generator seed for reproducibility (in case comparing multisample to single sample) 
            RooRandom::randomGenerator()->SetSeed(itoy+firstSample);
            RooDataSet *toy_bkg = bkg_pdf->generate(*genSpec) ;
            toy_bkg->Print();
            data[iy][itoy]->append(*toy_bkg);
            data[iy][itoy]->SetName(Form("data_y%i_b%i_subs%i",years[iy],q2Bin,itoy));
            wsp_out->import( *data[iy][itoy]);
            data[iy][itoy]->get(1);
            cout << "weight MC " << data[iy][itoy]->weight();
            toy_bkg->get(1);
            cout << "weight bkg " << toy_bkg->weight();
        }
        fout->cd();
        wsp_out->Write();
    }
    fout->Close();
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

  bool multiSample = false;
  uint nSample = 0;
  if ( argc > 3 && atoi(argv[3]) > 0 ) multiSample = true;
  if ( argc > 4 ) nSample = atoi(argv[4]);

  if (nSample==0) multiSample = false;

  int XGBv = 0; 
  if ( argc > 5 ) XGBv = atoi(argv[5]);

  bool localFiles = false;
  if ( argc > 6 && atoi(argv[6]) > 0 ) localFiles = true;

  bool plot = true;
  if ( argc > 7 && atoi(argv[7]) == 0 ) plot = false;

  int save = 0;
  if ( argc > 8 ) save = atoi(argv[8]);

  std::vector<int> years;
  if ( argc > 9 && atoi(argv[9]) != 0 ) years.push_back(atoi(argv[9]));
  else {
    cout << "No specific years selected, using default: 2016" << endl;
    years.push_back(2016);
  }
  if ( argc > 10 && atoi(argv[10]) != 0 ) years.push_back(atoi(argv[10]));
  if ( argc > 11 && atoi(argv[11]) != 0 ) years.push_back(atoi(argv[11]));

  cout <<  "q2Bin       " << q2Bin        << endl;
  cout <<  "parity      " << parity       << endl;
  cout <<  "multiSample " << multiSample  << endl;
  cout <<  "nSample     " << nSample      << endl;
  cout <<  "local files " << localFiles   << endl;
  cout <<  "plot        " << plot         << endl;
  cout <<  "save        " << save         << endl;
  cout <<  "years[0]    " << years[0]     << endl;
//   cout << years[1] << endl;
//   cout << years[2] << endl;


  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;

  if ( q2Bin==-1 )   cout << "Running all the q2 bins" << endl;
  if ( parity==-1 )  cout << "Running both the parity datasets" << endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      cocktailsample(q2Bin, parity, multiSample, nSample,XGBv, localFiles, plot, save, years);
  else
    cocktailsample(q2Bin, parity, multiSample, nSample, XGBv, localFiles, plot, save, years);

  return 0;

}



