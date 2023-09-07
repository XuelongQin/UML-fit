#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>
#include <list>
#include <map>
#include "RooStats/SPlot.h"
#include <string>
#include <TFile.h>
#include <RooDataSet.h>
#include <RooWorkspace.h>
#include <RooAbsPdf.h>
#include <RooCBShape.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooSimultaneous.h>
#include <RooDataHist.h>
#include <RooAbsReal.h>
#include <RooHistFunc.h>
#include <RooCategory.h>
#include "RooBernsteinSideband.h"
#include <RooProdPdf.h>
#include <RooRandom.h>
#include <RooAddPdf.h>
#include "RooDoubleCBFast.h"
#include <RooMinimizer.h>

#include "PdfSigMass.h"
#include "utils.h"
#include "PdfSigRTMass.h"
#include "PdfSigWTMass.h"



using namespace std;
using namespace ROOT;
using namespace RooFit;
using namespace RooStats;

double cos_theta_k;
double cos_theta_l;
double phi_kst_mumu;
double value;
double error;
double Fl = 0.;
double AFB = 0.;
double S3 = 0.;
double S4 = 0.;
double S5 = 0.;
double S6 = 0.;
double S7 = 0.;
double S8 = 0.;
double S9 = 0.;

double P1 = 0.;
double P2 = 0.;
double P3 = 0.;
double P4p = 0.;
double P5p = 0.;
double P6p = 0.;
double P8p = 0.;

double weight ;
double M_6s ;
double M_6c ;
double f_1s;
double f_3 ;
double f_4 ;
double f_5;
double f_6s ;
double f_6c;
double f_7 ;
double f_8 ;
double f_9 ;



string paralist[15] = {"FlS", "AFBS", "S3S", "S4S", "S5S", "S7S", "S8S", "S9S", "P1S", "P2S", "P3S", "P4pS", "P5pS", "P6pS", "P8pS"};
string momlist[9] = {"M1s","M6s","M6c","M3","M4","M5","M7","M8","M9"};
vector<double> vResult(9);
vector<double> vError(9);
vector<double> vPResult(15);
vector<double> vPError(15);

static const int nBins = 9;
TCanvas* c [4*nBins];


tuple<RooDataSet*,RooDataSet*,RooDataSet*> splot(RooWorkspace *ws_pdf ,int q2Bin, int year, int isample){
    RooAbsPdf *mass_pdf = ws_pdf->pdf(Form("PDF_mass_b%ip1_%i",q2Bin,year));
    RooRealVar *nsig = ws_pdf->var(Form("nsigb%ip1_%i",q2Bin,year));
    RooRealVar *nbkg = ws_pdf->var(Form("nbkgb%ip1_%i",q2Bin,year));  
    cout << "nsig is " << nsig->getVal() << endl;
    cout << "nbkg is " << nbkg->getVal() << endl; 
	//mass_pdf->Print();
    RooDataSet *data = (RooDataSet*)ws_pdf->data(Form("data_y%i_b%i_subs%i",year,q2Bin,isample));
	//data->Print();
    //RooWorkspace* ws_out = new RooWorkspace(Form("ws_Splotout_b%i_s100",q2Bin));
    RooStats::SPlot *sData = new RooStats::SPlot("sData", "An SPlot", *data, mass_pdf, RooArgSet(*nsig, *nbkg)); 
    //data->Print();

    std::cout << "\n\n------------------------------------------\n\nCheck SWeights:" << std::endl;

    std::cout << std::endl
                << "Yield of sig is\t" << nsig->getVal() << ".  From sWeights it is "
                << sData->GetYieldFromSWeight(Form("nsigb%ip1_%i",q2Bin,year)) << std::endl;

    std::cout << "Yield of bkg is\t" << nbkg->getVal() << ".  From sWeights it is "
                << sData->GetYieldFromSWeight(Form("nbkgb%ip1_%i",q2Bin,year)) << std::endl
                << std::endl;

    for (Int_t i = 0; i < 10; i++) {
        std::cout << "sig Weight for event " << i << std::right << std::setw(12) << sData->GetSWeight(i, Form("nsigb%ip1_%i",q2Bin,year)) << "  bkg Weight"
                    << std::setw(12) << sData->GetSWeight(i, Form("nbkgb%ip1_%i",q2Bin,year)) << "  Total Weight" << std::setw(12) << sData->GetSumOfEventSWeight(i)
                    << std::endl;
    }   

    ///ws_out->import(*totdata, Rename("momdataSplot"));
    RooDataSet *data_sig = new RooDataSet("SigSweight","Signal component",data,*data->get(),0,Form("nsigb%ip1_%i_sw",q2Bin,year));
    RooDataSet *data_bkg = new RooDataSet("BkgSweight","Background component",data,*data->get(),0,Form("nbkgb%ip1_%i_sw",q2Bin,year));
    return make_tuple(data,data_sig,data_bkg);
}

// the code used to constrct histo from efficiency 
tuple<TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*> Constructhisto(int q2Bin){
   	stringstream myString;
	myString.clear();
	myString.str("");
	myString << "/afs/cern.ch/user/x/xuqin/eos/Kmumu/momeff/1Deffb" << q2Bin << "y2016_od.root" ;
	//any year any parity is okay as we only need the binning results
	cout << "\n[Moment::ReCalValue]\tTry to open " << myString.str().c_str() << endl;
	TFile *efffile = new TFile(myString.str().c_str(), "READ");
	if (efffile!=0){
		cout<< "file is existed" << endl;
	}
	TH1D *f1seffW = (TH1D*)efffile->Get("f1sWTeff");
	TH1D *m6seffW = (TH1D*)efffile->Get("m6sWTeff");
	TH1D *m6ceffW = (TH1D*)efffile->Get("m6cWTeff");
	TH1D *f3effW = (TH1D*)efffile->Get("f3WTeff");
	TH1D *f4effW = (TH1D*)efffile->Get("f4WTeff");
	TH1D *f5effW = (TH1D*)efffile->Get("f5WTeff");
	TH1D *f7effW = (TH1D*)efffile->Get("f7WTeff");
	TH1D *f8effW = (TH1D*)efffile->Get("f8WTeff");
	TH1D *f9effW = (TH1D*)efffile->Get("f9WTeff");


	int bins = f1seffW->GetNbinsX();
	//TAxis *axis = m6seffW->GetXaxis();
	double xbinsf1s[bins+1], xbinsm6s[bins+1], xbinsm6c[bins+1], xbinsf3[bins+1], xbinsf4[bins+1], xbinsf5[bins+1], xbinsf7[bins+1], xbinsf8[bins+1], xbinsf9[bins+1];
	for (int i=0;i<=bins;i++){
		xbinsf1s[i] = f1seffW->GetBinLowEdge(i+1);
		xbinsm6s[i] = m6seffW->GetBinLowEdge(i+1);
		xbinsm6c[i] = m6ceffW->GetBinLowEdge(i+1);
		xbinsf3[i] = f3effW->GetBinLowEdge(i+1);
		xbinsf4[i] = f4effW->GetBinLowEdge(i+1);
		xbinsf5[i] = f5effW->GetBinLowEdge(i+1);
		xbinsf7[i] = f7effW->GetBinLowEdge(i+1);
		xbinsf8[i] = f8effW->GetBinLowEdge(i+1);
		xbinsf9[i] = f9effW->GetBinLowEdge(i+1);

		//    cout << xbins[i] << endl;
	}


	TH1D *f1s = new TH1D("f1s","",bins,xbinsf1s);
	TH1D *m6s = new TH1D("m6s","",bins,xbinsm6s);
	TH1D *m6c = new TH1D("m6c","",bins,xbinsm6c);
	TH1D *f3 = new TH1D("f3","",bins,xbinsf3);
	TH1D *f4 = new TH1D("f4","",bins,xbinsf4);
	TH1D *f5 = new TH1D("f5","",bins,xbinsf5);
	TH1D *f7 = new TH1D("f7","",bins,xbinsf7);
	TH1D *f8 = new TH1D("f8","",bins,xbinsf8);
	TH1D *f9 = new TH1D("f9","",bins,xbinsf9); 

    return make_tuple(f1s,m6s,m6c,f3,f4,f5,f7,f8,f9);
}

tuple<double,double> MoMValue(TH1D *his){
    double mvalue = his->GetMean();
    double merror = his->GetRMS()/TMath::Sqrt(his->GetEntries()-1);
    return make_tuple(mvalue,merror);
}


tuple<vector<double>,vector<double>,vector<double>,vector<double>> CalValue(TH1D *f1s, TH1D *m6s, TH1D *m6c, TH1D *f3, TH1D *f4, TH1D *f5, TH1D *f7, TH1D *f8, TH1D *f9){
    vector<double> momvalue(15);
    vector<double> momerror(15); 
    auto M1stuple = MoMValue(f1s);
    auto M6stuple = MoMValue(m6s);
    auto M6ctuple = MoMValue(m6c);
    auto M3tuple = MoMValue(f3);
    auto M4tuple = MoMValue(f4);
    auto M5tuple = MoMValue(f5);
    auto M7tuple = MoMValue(f7);
    auto M8tuple = MoMValue(f8);
    auto M9tuple = MoMValue(f9);

    momvalue[0] = get<0>(M1stuple);
    momerror[0] = get<1>(M1stuple);
    momvalue[1] = get<0>(M6stuple);
    momerror[1] = get<1>(M6stuple);
    momvalue[2] = get<0>(M6ctuple);
    momerror[2] = get<1>(M6ctuple);
    momvalue[3] = get<0>(M3tuple);
    momerror[3] = get<1>(M3tuple);
    momvalue[4] = get<0>(M4tuple);
    momerror[4] = get<1>(M4tuple);
    momvalue[5] = get<0>(M5tuple);
    momerror[5] = get<1>(M5tuple);
    momvalue[6] = get<0>(M7tuple);
    momerror[6] = get<1>(M7tuple);
    momvalue[7] = get<0>(M8tuple);
    momerror[7] = get<1>(M8tuple);
    momvalue[8] = get<0>(M9tuple);
    momerror[8] = get<1>(M9tuple);

    vector<double> parvalue(15);
    vector<double> parerror(15);

//Fl
    parvalue[0] = 2.0 - 2.5 * momvalue[0];
    parerror[0] = 2.5 * momerror[0];

//AFBS
    parvalue[1] = 3.0/4.0*5.0/2.0*momvalue[1];
    parerror[1] = 3.0/4.0*5.0/2.0*momerror[1];

//S3S
    parvalue[2] = 25.0/8.0*momvalue[3];
    parerror[2] = 25.0/8.0*momerror[3];

//S4S
    parvalue[3] = 25.0/8.0*momvalue[4];
    parerror[3] = 25.0/8.0*momerror[4];

//S5S
    parvalue[4] = 5.0/2.0*momvalue[5];
    parerror[4] = 5.0/2.0*momerror[5];

//S7S
    parvalue[5] = 5.0/2.0*momvalue[6];
    parerror[5] = 5.0/2.0*momerror[6];

//S8S
    parvalue[6] = 25.0/8.0*momvalue[7];
    parerror[6] = 25.0/8.0*momerror[7];

//S9S
    parvalue[7] = 25.0/8.0*momvalue[8];
    parerror[7] = 25.0/8.0*momerror[8];

//P1S
    parvalue[8] = 2*parvalue[2]/(1-parvalue[0]);
    parerror[8] = TMath::Sqrt(TMath::Power(2*parerror[2]/(1-parvalue[0]),2)+TMath::Power(2*parvalue[2]*parerror[0]/TMath::Power(1-parvalue[0],2),2));

//P2S
    parvalue[9] = 2.0/3.0*parvalue[1]/(1-parvalue[0]);
    //cout << "P2S is " << parvalue[9] << endl;
    parerror[9] = 2.0/3.0 * TMath::Sqrt(TMath::Power((1.0/(1.0-parvalue[0])*parerror[1]),2)+TMath::Power((parerror[0]*parvalue[1]/TMath::Power(1.0-parvalue[0],2)),2));
    //cout << "error of P2S is " << parerror[9] << endl;
//P3S
    parvalue[10] = -1.0*parvalue[7]/(1.0-parvalue[0]);
    parerror[10] = TMath::Sqrt(TMath::Power((1.0/(1.0-parvalue[0])*parerror[7]),2)+TMath::Power((parerror[0]*parvalue[7]/TMath::Power(1.0-parvalue[0],2)),2));

//P4pS
    parvalue[11] = 2.0*parvalue[3]/TMath::Sqrt(parvalue[0]*(1.0-parvalue[0]));
    parerror[11] = 2.0*TMath::Sqrt(TMath::Power(parerror[3]*1.0/TMath::Sqrt(parvalue[0]*(1.0-parvalue[0])),2)+TMath::Power((parerror[0]*(1.0-2.0*parvalue[0]*parvalue[3]/(2.0*TMath::Power((1.0-parvalue[0])*parvalue[0],3.0/2.0)))),2));

//P5pS
    parvalue[12] = parvalue[4]/TMath::Sqrt(parvalue[0]*(1.0-parvalue[0]));
    parerror[12] = TMath::Sqrt(TMath::Power(parerror[4]*1.0/TMath::Sqrt(parvalue[0]*(1.0-parvalue[0])),2)+TMath::Power((parerror[0]*(1.0-2.0*parvalue[0]*parvalue[4]/(2.0*TMath::Power((1.0-parvalue[0])*parvalue[0],3.0/2.0)))),2));

//P6pS
    parvalue[13] = -parvalue[5]/TMath::Sqrt(parvalue[0]*(1.0-parvalue[0]));
    parerror[13] = TMath::Sqrt(TMath::Power(parerror[5]*1.0/TMath::Sqrt(parvalue[0]*(1.0-parvalue[0])),2)+TMath::Power((parerror[0]*(1.0-2.0*parvalue[0]*parvalue[5]/(2.0*TMath::Power((1.0-parvalue[0])*parvalue[0],3.0/2.0)))),2));

//P8pS
    parvalue[14] = 2.0 * parvalue[6]/TMath::Sqrt(parvalue[0]*(1.0-parvalue[0]));
    parerror[14] = 2.0*TMath::Sqrt(TMath::Power(parerror[6]*1.0/TMath::Sqrt(parvalue[0]*(1.0-parvalue[0])),2)+TMath::Power((parerror[0]*(1.0-2.0*parvalue[0]*parvalue[6]/(2.0*TMath::Power((1.0-parvalue[0])*parvalue[0],3.0/2.0)))),2));  

    return make_tuple(momvalue,momerror,parvalue,parerror);


}





void Recocock_MOM(int q2Bin,std::vector<int> years, int iSample){
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
    int parity = 1;
    string all_years = "";
    string yearstring = ""; 
    string XGBstr = "";
    bool multiSample = true;
    bool localFiles = false;
    int XGBv = 8;
    if (XGBv>9) XGBstr = Form("-TMVAv%i",XGBv-10);
    else if (XGBv>0) XGBstr = Form("_XGBv%i",XGBv);
    RooArgSet c_vars_rt, c_pdfs_rt;
    RooArgSet c_vars_wt, c_pdfs_wt;
    RooArgSet c_vars; 
    RooWorkspace * ws_pars = new RooWorkspace("ws_pars");

    uint firstSample = iSample;
    uint lastSample = iSample;
    string shortString = Form("b%ip%i",q2Bin,parity);
    cout<<"Conf: "<<shortString<<endl;

    std::vector<RooWorkspace*> wsp, wsp_mcmass,wsp_sb;
    std::vector<std::vector<RooDataSet*>> data;
    double power = 1.0;
    gInterpreter->GenerateDictionary("std::pair<std::string, RooDataSet*>", "map;string;RooDataSet.h");
    gInterpreter->GenerateDictionary("std::map<std::string, RooDataSet*>",  "map;string;RooDataSet.h");
    gInterpreter->GenerateDictionary("std::pair<std::map<string,RooDataSet*>::iterator, bool>", "map;string;RooDataSet.h");
    std::map<std::string, RooDataSet*> map;
    // create workspace to save results to out file (or for internal usage)
    /*RooWorkspace* wksp = new RooWorkspace(((multiSample?"wsMulti_":"ws_")+shortString+Form("_s%i_pow%.1f",nSample,power)).c_str(),
                        (multiSample?"Workspace with set of RECO subsample mom results":
                        (nSample>0?"Workspace with RECO subsample mom result":
                        "Workspace with full RECO mom result")));*/

    RooRealVar* ctK = new RooRealVar("ctK", "cos(#theta_{K})", -1  , 1  );
    RooRealVar* ctL = new RooRealVar("ctL", "cos(#theta_{l})", -1  , 1  );
    RooRealVar* phi = new RooRealVar("phi", "#phi", -3.14159, 3.14159  );
    RooRealVar* mass = new RooRealVar("mass","m(#mu#muK#pi)", 5.,5.6,"GeV");
    RooArgList vars (* ctK,* ctL,* phi);
    RooRealVar* rand = new RooRealVar("rand", "rand", 0,1);
    RooArgSet reco_vars ( *mass);
    RooArgSet observables (*mass,*ctK,*ctL,*phi);
    string isample = "";
    RooCategory sample ("sample", "sample");
    string bkgpdfname = "bkg_pdf_%i";
    string year = ""; 
    TRandom gen_nevt;
    std::vector<RooAbsPdf*> PDF_sig_mass(0), pdf_sig_mass, pdf_bkg_mass, final_PDF;
    std::vector<RooGaussian*> c_deltaPeaks, c_fm;
    for (unsigned int iy = 0; iy < years.size(); iy++) {
        year.clear(); year.assign(Form("%i",years[iy]));
        all_years += year;
        for (uint is = firstSample; is <= lastSample; is++) {
            isample.clear(); isample.assign( Form("%i",is) );
            sample.defineType(("data"+year+"_subs"+isample).c_str());
        }
    }

    RooSimultaneous* simPdf = new RooSimultaneous("simPdf", "simultaneous pdf", sample);
    cout << "construct simpdf" << endl;

    string outname = "./Momresults/cocktail/cocktail_" ;
    outname = outname + all_years + Form("_b%i_subs%i.root",q2Bin,iSample);
    TFile *fout = new TFile(outname.c_str(),"RECREATE");
    cout << "output tree will be saved in " << outname << endl;
    TTree* momResultsTree = new TTree("momResultsTree","momResultsTree");

    for (unsigned int iMom=0; iMom<vResult.size(); ++iMom){
        momResultsTree->Branch(Form("%s_value",momlist[iMom].c_str()),&vResult[iMom]);
        momResultsTree->Branch(Form("%s_error",momlist[iMom].c_str()),&vError[iMom]);
    }

    TTree* parResultsTree = new TTree("parResultsTree","parResultsTree");
    for (unsigned int iPar=0; iPar<vPResult.size(); ++iPar){
        parResultsTree->Branch(Form("%s_value",paralist[iPar].c_str()),&vPResult[iPar]);
        parResultsTree->Branch(Form("%s_error",paralist[iPar].c_str()),&vPError[iPar]);
    }

    string filename_data = Form("recoMCDataset_b%i.root",q2Bin);
    filename_data = "/afs/cern.ch/user/x/xuqin/eos/Kmumu/mom_cocktail/" + filename_data;
    TFile *inputdata = new TFile(filename_data.c_str(),"read");
    RooWorkspace *wsp_cock = (RooWorkspace*)inputdata->Get(Form("ws_b%ip0",q2Bin));
    for (unsigned int iy = 0; iy < years.size(); iy++){
        year.clear(); year.assign(Form("%i",years[iy]));

        std::vector<RooDataSet*> singyeardata;
        for (uint is = firstSample; is <= lastSample; is++) {
            RooDataSet *datasub = (RooDataSet*)wsp_cock->data(Form("data_y%i_b%i_subs%i",years[iy],q2Bin,is))->reduce("mass>5.0 && mass <5.6");
            singyeardata.push_back(datasub);
        }
        data.push_back(singyeardata);

        // Signal Mass Component
        // import mass PDF from fits to the MC
        //mass->setRange("sbleft",  5., 5.6);
        // Mass Component
        // import mass PDF from fits to the MC
        string filename_mc_mass = Form("/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/massFits-%s-XGBv8/%i.root",q2Bin==4?"Jpsi":(q2Bin==6?"Psi":"LMNR"),years[iy]);
        if (localFiles)
            filename_mc_mass = Form("results_fits_%i_fM_%s.root",years[iy], q2Bin==4?"Jpsi":(q2Bin==6?"Psi":"lmnr"));
        if (!retrieveWorkspace( filename_mc_mass, wsp_mcmass, "w"))  return;

        //wsp_mcmass[iy]->Print();
        wsp_mcmass[iy]->loadSnapshot(Form("reference_fit_RT_%i",q2Bin));
        RooRealVar* mean_rt       = new RooRealVar (Form("mean_{RT}^{%i}",years[iy])    , "massrt"      , wsp_mcmass[iy]->var(Form("mean_{RT}^{%i}",q2Bin))->getVal()     ,      5,    6, "GeV");
        RooRealVar* sigma_rt      = new RooRealVar (Form("#sigma_{RT1}^{%i}",years[iy] ), "sigmart1"    , wsp_mcmass[iy]->var(Form("#sigma_{RT1}^{%i}",q2Bin))->getVal()  ,      0,    1, "GeV");
        RooRealVar* alpha_rt1     = new RooRealVar (Form("#alpha_{RT1}^{%i}",years[iy] ), "alphart1"    , wsp_mcmass[iy]->var(Form("#alpha_{RT1}^{%i}", q2Bin))->getVal() ,      0,   10 );
        RooRealVar* n_rt1         = new RooRealVar (Form("n_{RT1}^{%i}",years[iy])      , "nrt1"        , wsp_mcmass[iy]->var(Form("n_{RT1}^{%i}", q2Bin))->getVal()      ,      0.01,  100.);

        RooAbsPdf* dcb_rt;
        RooRealVar* alpha_rt2 = new RooRealVar (Form("#alpha_{RT2}^{%i}",years[iy] ), "alphart2"    , 0,    -10,   10 );
        RooRealVar* n_rt2     = new RooRealVar (Form("n_{RT2}^{%i}",years[iy])      , "nrt2"        , 0.01,   0.01,  100.);
        RooRealVar* sigma_rt2 = new RooRealVar (Form("#sigma_{RT2}^{%i}",years[iy] ), "sigmaRT2"  ,   0 , 0,   0.12, "GeV");
        RooRealVar* f1rt      = new RooRealVar (Form("f^{RT%i}",years[iy])          , "f1rt"      ,   0 , 0.,  1.);
        if (q2Bin != 7){
        alpha_rt2 -> setVal(wsp_mcmass[iy]->var(Form("#alpha_{RT2}^{%i}", q2Bin))->getVal() );
        n_rt2     -> setVal(wsp_mcmass[iy]->var(Form("n_{RT2}^{%i}", q2Bin)     )->getVal() );
        }
        if (q2Bin >= 4){
        sigma_rt2-> setVal(wsp_mcmass[iy]->var(Form("#sigma_{RT2}^{%i}",q2Bin))->getVal() );
        f1rt     -> setVal(wsp_mcmass[iy]->var(Form("f^{RT%i}", q2Bin))->getVal() );
        if (q2Bin < 7) 
            dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, alpha_rt2, n_rt1, n_rt2 ,f1rt, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt );
        else
            dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, n_rt1 ,f1rt, q2Bin, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt );
        
        } 
        else{
            alpha_rt2->setRange(0,10);
            dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, alpha_rt1, alpha_rt2, n_rt1, n_rt2 , wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt  );
        }
        
        /// create constrained PDF for RT mass
        RooArgList constr_rt_list = RooArgList(c_pdfs_rt);
        constr_rt_list.add(*dcb_rt);
        RooProdPdf * c_dcb_rt = new RooProdPdf(("c_dcb_rt_"+year).c_str(), ("c_dcb_rt_"+year).c_str(), constr_rt_list );
        c_vars.add(c_vars_rt);
    
        /// create WT component
        wsp_mcmass[iy]->loadSnapshot(Form("reference_fit_WT_%i",q2Bin));

    //     RooRealVar* mean_wt     = new RooRealVar (Form("mean_{WT}^{%i}",years[iy])      , "masswt"     ,  wsp_mcmass[iy]->var(Form("mean_{WT}^{%i}", q2Bin))->getVal()    ,      5,    6, "GeV");
        RooRealVar* sigma_wt    = new RooRealVar (Form("#sigma_{WT1}^{%i}",years[iy])   , "sigmawt"    ,  wsp_mcmass[iy]->var(Form("#sigma_{WT1}^{%i}", q2Bin))->getVal() ,      0,    1, "GeV");
        RooRealVar* alpha_wt1   = new RooRealVar (Form("#alpha_{WT1}^{%i}",years[iy] )  , "alphawt1"   ,  wsp_mcmass[iy]->var(Form("#alpha_{WT1}^{%i}", q2Bin))->getVal() ,      0,   10 );
        RooRealVar* alpha_wt2   = new RooRealVar (Form("#alpha_{WT2}^{%i}",years[iy] )  , "alphawt2"   ,  wsp_mcmass[iy]->var(Form("#alpha_{WT2}^{%i}", q2Bin))->getVal() ,      0,   10 );
        RooRealVar* n_wt1       = new RooRealVar (Form("n_{WT1}^{%i}",years[iy])        , "nwt1"       ,  wsp_mcmass[iy]->var(Form("n_{WT1}^{%i}", q2Bin))->getVal()      ,      0.01, 100.);
        RooRealVar* n_wt2       = new RooRealVar (Form("n_{WT2}^{%i}",years[iy])        , "nwt2"       ,  wsp_mcmass[iy]->var(Form("n_{WT2}^{%i}", q2Bin))->getVal()      ,      0.01, 100.);
        double mean_wt_val = wsp_mcmass[iy]->var(Form("mean_{WT}^{%i}", q2Bin))->getVal(); 
        double mean_wt_err = wsp_mcmass[iy]->var(Form("mean_{WT}^{%i}", q2Bin))->getError(); 

        RooAbsPdf* dcb_wt;
        double deltaPeakValue = mean_rt->getVal()-mean_wt_val;
        double mean_rt_err = wsp_mcmass[iy]->var(Form("mean_{RT}^{%i}", q2Bin))->getError(); 
        double deltaPeakError = sqrt(mean_rt_err*mean_rt_err+mean_wt_err*mean_wt_err); 

        RooRealVar* deltaPeakVar = new RooRealVar ( Form("deltaPeakVar^{%i}", years[iy]),Form("deltaPeakVar^{%i}", years[iy]), deltaPeakValue, 0., 0.2) ;
        RooGaussian* c_deltaPeaks = new RooGaussian(Form("deltaPeaks^{%i}"  , years[iy]) , "c_deltaPeaks", *deltaPeakVar, RooConst( deltaPeakValue ), RooConst(deltaPeakError )); // value to be checked
        RooFormulaVar*mWT_data = new  RooFormulaVar(Form("mWT_data^{%i}",years[iy]), "@0 + @1", RooArgList(*mean_rt, *deltaPeakVar));
        RooRealVar* mean_wt     = (RooRealVar*)mWT_data;
        
        c_pdfs_wt.add(*c_deltaPeaks);
        c_vars_wt.add(*deltaPeakVar);

        dcb_wt = createWTMassShape(q2Bin, mass, mean_wt, sigma_wt, alpha_wt1, alpha_wt2, n_wt1, n_wt2 , wsp_mcmass[iy], years[iy], true, c_vars_wt, c_pdfs_wt );

        RooArgList constr_wt_list = RooArgList(c_pdfs_wt);
        constr_wt_list.add(*dcb_wt);
        RooProdPdf * c_dcb_wt = new RooProdPdf(("c_dcb_wt_"+year).c_str(), ("c_dcb_wt_"+year).c_str(), constr_wt_list );
        c_vars.add(c_vars_wt);

        cout << "deltaPeak constraint built, of value " << deltaPeakValue << " +/- " << deltaPeakError << endl;


        // original
        //     RooAbsPdf* dcb_wt = createWTMassShape(q2Bin, mass, mean_wt, sigma_wt, alpha_wt1, alpha_wt2, n_wt1, n_wt2 , wsp_mcmass[iy], years[iy], true, c_vars_wt, c_pdfs_wt );
        // 
        //     /// create constrained PDF for WT mass
        //     RooArgList constr_wt_list = RooArgList(c_pdfs_wt);
        //     constr_wt_list.add(*dcb_wt);
        //     RooProdPdf * c_dcb_wt = new RooProdPdf(("c_dcb_wt_"+year).c_str(), ("c_dcb_wt_"+year).c_str(), constr_wt_list );
        //     c_vars.add(c_vars_wt);



        RooRealVar* mFrac = new RooRealVar(Form("mFrac^{%i}",years[iy]),"mistag fraction",0.13, 0, 1);

        if (q2Bin < 5){
            PDF_sig_mass.push_back( new PdfSigMass(("PDF_sig_mass_"+shortString+"_"+year).c_str(),
                                                ("PDF_sig_mass_"+year).c_str(),
                                                *mass,
                                                *mean_rt, *sigma_rt, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2,
                                                *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,
                                            *mFrac,
                                            *c_dcb_rt,
                                            *c_dcb_wt  
                                            ));}
        else if (q2Bin==7){
            PDF_sig_mass.push_back( new PdfSigMass(("PDF_sig_mass_"+shortString+"_"+year).c_str(),
                                                ("PDF_sig_mass_"+year).c_str(),
                                                *mass,
                                                *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *n_rt1, *f1rt,
                                                *mean_wt, *sigma_wt,             *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,
                                            *mFrac,
                                            *c_dcb_rt ,
                                            *c_dcb_wt
                                            ));}
        else{
        PDF_sig_mass.push_back( new PdfSigMass(("PDF_sig_mass_"+shortString+"_"+year).c_str(),
                                            ("PDF_sig_mass_"+year).c_str(),
                                            *mass,
                                            *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2, *f1rt,
                                            *mean_wt, *sigma_wt,             *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,
                                        *mFrac,
                                        *c_dcb_rt ,
                                        *c_dcb_wt
                                        ));
        }


        /// create constraint on mFrac (here there is no efficiency, therefore value set to measured value on MC)
        double nrt_mc   =  wsp_mcmass[iy]->var(Form("nRT_%i",q2Bin))->getVal(); 
        double nwt_mc   =  wsp_mcmass[iy]->var(Form("nWT_%i",q2Bin))->getVal(); 
        //double fraction = nwt_mc / (nrt_mc + nwt_mc);
        double fraction = nwt_mc/nrt_mc;
        c_fm.push_back(new RooGaussian(Form("c_fm^{%i}",years[iy]) , "c_fm" , *mFrac,  
                                        RooConst(fraction) , 
                                        RooConst(fM_sigmas[years[iy]][q2Bin]*(nrt_mc + nwt_mc)/nrt_mc)
                                        ) );
        cout << fraction << "   " << fM_sigmas[years[iy]][q2Bin]*(nrt_mc + nwt_mc)/nrt_mc << endl;                                    
        c_vars.add(*mFrac);  

        pdf_sig_mass.push_back(new RooProdPdf(("pdf_sig_mass_"+year).c_str(), 
                                            ("pdf_sig_mass_"+year).c_str(),
                                            RooArgList(*PDF_sig_mass[iy],*c_fm[iy])));



        string filename_sb = Form("/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/sidebands/savesb_%i_b%i.root", years[iy],q2Bin);
        cout << "Read sideband file " << filename_sb << endl;
        retrieveWorkspace( filename_sb, wsp_sb, "wsb");
        // read mass pdf for background
        RooRealVar* slope       = new RooRealVar    (Form("slope^{%i}",years[iy]),  Form("slope^{%i}",years[iy]) , wsp_sb[iy]->var("slope")->getVal(), -10., 0.);
        // RooRealVar* slope       = new RooRealVar    (Form("slope^{%i}",years[iy]),  Form("slope^{%i}",years[iy]) , wsp_sb[iy]->var("slope")->getVal(), -10., 0.);
        pdf_bkg_mass.push_back(new RooExponential(Form("pdf_bkg_mass_%i",years[iy]),  Form("pdf_bkg_mass_%i",years[iy]) ,  *mass,   *slope  ));
        int nbkg_ex = nbkg_years[years[iy]][q2Bin];
        int nsig_ex = nsig_years[years[iy]][q2Bin];
        RooRealVar *nsig = new RooRealVar(("nsig"+shortString+"_"+year).c_str(),"number of signal",nsig_ex,0,nsig_ex+nbkg_ex);
        RooRealVar *nbkg = new RooRealVar(("nbkg"+shortString+"_"+year).c_str(),"number of bkg",nbkg_ex,0,nsig_ex+nbkg_ex);
        final_PDF.push_back( new RooAddPdf( ("PDF_mass_"+shortString+"_"+year).c_str(),
                                            ("PDF_mass_"+shortString+"_"+year).c_str(),
                                            RooArgList(*pdf_sig_mass[iy], *pdf_bkg_mass[iy]),
                                            RooArgList(*nsig,*nbkg)
                                        ));


        for (uint is = firstSample; is <= lastSample; is++) {
            if ( !data[iy][0] || data[iy][0]->IsZombie() ) {
                cout<<"Dataset " << is  << " not found in file: "<<filename_data<<endl;
                return;
            }
            cout << " xixi " << endl;
            ws_pars->import(*data[iy][0]);
            cout << "entries of is " << is << " is " << data[iy][0]->sumEntries() <<endl;   
            map.insert( map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year+Form("_subs%d",is)).c_str(), data[iy][0]) );
            simPdf->addPdf(*final_PDF[iy], ("data"+year+Form("_subs%d",is)).c_str());
        }
    }

    //TFile* fout = new TFile(("simFitResults4d/simFitResult_recoMC_fullMass" + all_years + stat + Form("_b%i.root", q2Bin)).c_str(),"UPDATE");

    // save initial par values into a workspace 
    ws_pars->import(*simPdf,RecycleConflictNodes());
    RooArgSet *params = (RooArgSet *)simPdf->getParameters(*mass);
    // The kTRUE flag imports the values of the objects in (*params) into the workspace
    // If not set, the present values of the workspace parameters objects are stored
    ws_pars->saveSnapshot("initial_pars", *params, kTRUE);


    // Construct combined dataset in (x,sample)
    RooDataSet allcombData ("allcombData", "combined data", 
                            reco_vars,
                            Index(sample), 
                            Import(map)); 

    RooDataSet* combData = 0;
    RooAbsReal* nll = 0;
 
    for (uint is = firstSample; is <= lastSample; is++) {

        string the_cut = Form("sample==sample::data%d_subs%d", years[0], is);
        if (years.size() > 1){
            for (unsigned int iy=1; iy < years.size(); iy++){
                the_cut = the_cut + Form("|| sample==sample::data%d_subs%d", years[iy], is);
        }
        }
        combData = (RooDataSet*)allcombData.reduce(Cut(the_cut.c_str()));
        cout<<"Fitting subsample "<<is<<" with "<<combData->numEntries()<<" entries"<<endl;
        ws_pars->loadSnapshot("initial_pars");

        TStopwatch subTime;
        nll = ws_pars->pdf("simPdf")->createNLL(*combData,
                                                RooFit::Extended(kFALSE),
                                                RooFit::Constrain(c_vars),
                                                RooFit::NumCPU(1)
                                                );
            
        RooMinimizer m(*nll) ;   
        m.optimizeConst (kTRUE); // do not recalculate constant terms
        m.setOffsetting(kTRUE);  //  Enable internal likelihood offsetting for enhanced numeric precision.
        m.setPrintLevel(-1);
        m.setPrintEvalErrors(-1);
        m.setMinimizerType("Minuit2");
        //     m.setVerbose(kTRUE); 

        subTime.Start(true);
        m.setStrategy(0);
        m.migrad() ;
        m.hesse() ;
        m.setStrategy(2);
        m.migrad() ;
        m.hesse() ;
        subTime.Stop();
        cout << "fitting done  " << subTime.CpuTime() << endl;
        
        RooFitResult* fitResult = m.save(("result_" + shortString + Form("subs%d",is)).c_str()) ; 
        fitResult->Print("v");
        fout->cd();
        fitResult->Write();


        auto histo = Constructhisto(q2Bin);
        TH1D *f1s = get<0>(histo);
        TH1D *m6s = get<1>(histo);
        TH1D *m6c = get<2>(histo);
        TH1D *f3 = get<3>(histo);
        TH1D *f4 = get<4>(histo);
        TH1D *f5 = get<5>(histo);
        TH1D *f7 = get<6>(histo);
        TH1D *f8 = get<7>(histo);
        TH1D *f9 = get<8>(histo);
        cout << "begin to fill histograms " << endl;
        for (unsigned int iy = 0; iy < years.size(); iy++){
            
            cout << "\n[Moment::ReCalValue]\tTry to fill year" << years[iy] << "reco mc sample" << endl; 
            stringstream myString;
            myString.clear();
            myString.str("");
            //myString << "/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/moment/1Deff/momentfinal/eff/1Deffb" << q2Bin << "y"<< year <<".root" ;
            myString << "/afs/cern.ch/user/x/xuqin/eos/Kmumu/momeff/1Deffb" << q2Bin << "y" << years[iy] << "_" << "od" << ".root" ;
            //odd efficiency for even event; even efficiency for odd event
            cout << "\n[Moment::ReCalValue]\tTry to open" << myString.str().c_str() << endl;
            TFile *efffile = new TFile(myString.str().c_str(), "READ");
            if (efffile!=0){
                cout<< "file is existed" << endl;
            }

            TH1D *f1sweight;
            TH1D *m6sweight;
            TH1D *m6cweight;
            TH1D *f3weight;
            TH1D *f4weight;
            TH1D *f5weight;
            TH1D *f7weight;
            TH1D *f8weight;
            TH1D *f9weight;

            TH1D *m6sdiffsignweight;
            TH1D *m6cdiffsignweight;
            TH1D *f5diffsignweight;
            TH1D *f8diffsignweight;
            TH1D *f9diffsignweight;

            f1sweight = (TH1D*)efffile->Get("f1stotweight");
            m6sweight = (TH1D*)efffile->Get("m6ssamesignweight");
            m6cweight = (TH1D*)efffile->Get("m6csamesignweight");
            f3weight = (TH1D*)efffile->Get("f3totweight");
            f4weight = (TH1D*)efffile->Get("f4totweight");
            f5weight = (TH1D*)efffile->Get("f5samesignweight");
            f7weight = (TH1D*)efffile->Get("f7totweight");
            f8weight = (TH1D*)efffile->Get("f8samesignweight");
            f9weight = (TH1D*)efffile->Get("f9samesignweight");

            m6sdiffsignweight = (TH1D*)efffile->Get("m6sdiffsignweight");
            m6cdiffsignweight = (TH1D*)efffile->Get("m6cdiffsignweight");
            f5diffsignweight = (TH1D*)efffile->Get("f5diffsignweight");
            f8diffsignweight = (TH1D*)efffile->Get("f8diffsignweight");
            f9diffsignweight = (TH1D*)efffile->Get("f9diffsignweight");

            cout << "Read efficiency successfully" << endl;
            //ws_pars->Print();
            auto datasplot = splot(ws_pars,q2Bin,years[iy],is);
            RooDataSet *tot_data = get<0>(datasplot);
            RooDataSet *sig_data = get<1>(datasplot);
            RooDataSet *bkg_data = get<2>(datasplot);


            int totentries = sig_data->numEntries();
            //sig_data->Print();
            cout << "total entries in total data is " << totentries << endl;

            for (int i=0;i<totentries;i++){

                cos_theta_k = sig_data->get(i)->getRealValue("ctK");
                cos_theta_l = sig_data->get(i)->getRealValue("ctL");
                phi_kst_mumu = sig_data->get(i)->getRealValue("phi");
                weight = sig_data->weight();            
                f_1s = 1 - cos_theta_k * cos_theta_k;
                M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
                M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
                f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
                f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);	                

                double weightf1s,weightm6s,weightm6c,weightf3,weightf4,weightf5,weightf7,weightf8,weightf9;
                double weightdiffm6s, weightdiffm6c,weightdifff5,weightdifff8,weightdifff9;

                weightf1s = f1sweight->GetBinContent(f1sweight->GetXaxis()->FindBin(f_1s));
                weightm6s = m6sweight->GetBinContent(m6sweight->GetXaxis()->FindBin(M_6s));
                weightm6c = m6cweight->GetBinContent(m6cweight->GetXaxis()->FindBin(M_6c));
                weightf3 = f3weight->GetBinContent(f3weight->GetXaxis()->FindBin(f_3));
                weightf4 = f4weight->GetBinContent(f4weight->GetXaxis()->FindBin(f_4));
                weightf5 = f5weight->GetBinContent(f5weight->GetXaxis()->FindBin(f_5));
                weightf7 = f7weight->GetBinContent(f7weight->GetXaxis()->FindBin(f_7));
                weightf8 = f8weight->GetBinContent(f8weight->GetXaxis()->FindBin(f_8));
                weightf9 = f9weight->GetBinContent(f9weight->GetXaxis()->FindBin(f_9));
                f1s->Fill(f_1s,weight*weightf1s);
                f3->Fill(f_3,weight*weightf3);
                f4->Fill(f_4,weight*weightf4);
                f7->Fill(f_7,weight*weightf7);


                m6s->Fill(M_6s, weight*weightm6s);
                m6c->Fill(M_6c, weight*weightm6c);
                f5->Fill(f_5, weight*weightf5);
                f8->Fill(f_8, weight*weightf8);
                f9->Fill(f_9,weight*weightf9);

                weightdiffm6s = m6sdiffsignweight->GetBinContent(m6sdiffsignweight->GetXaxis()->FindBin(M_6s));
                weightdiffm6c = m6cdiffsignweight->GetBinContent(m6cdiffsignweight->GetXaxis()->FindBin(M_6c));
                weightdifff5 = f5diffsignweight->GetBinContent(f5diffsignweight->GetXaxis()->FindBin(f_5));
                weightdifff8 = f8diffsignweight->GetBinContent(f8diffsignweight->GetXaxis()->FindBin(f_8));
                weightdifff9 = f9diffsignweight->GetBinContent(f9diffsignweight->GetXaxis()->FindBin(f_9));


                m6s->Fill(-M_6s, weight*weightdiffm6s);
                m6c->Fill(-M_6c, weight*weightdiffm6c);
                f5->Fill(-f_5, weight*weightdifff5);
                f8->Fill(-f_8, weight*weightdifff8);
                f9->Fill(-f_9,weight*weightdifff9);

            }
        }
        auto resulttuple = CalValue(f1s,m6s,m6c,f3,f4,f5,f7,f8,f9);
        vResult.assign(get<0>(resulttuple).begin(),get<0>(resulttuple).end());
        vError.assign(get<1>(resulttuple).begin(),get<1>(resulttuple).end());
        momResultsTree->Fill();
        vPResult.assign(get<2>(resulttuple).begin(),get<2>(resulttuple).end());
        vPError.assign(get<3>(resulttuple).begin(),get<3>(resulttuple).end());
        parResultsTree->Fill();

        string plotString = shortString + "_" + Form("subs%i_",is) + all_years;
        //if (nSample>0) plotString = plotString + Form("_s%i",nSample);

        int confIndex = 0;
        string longString  = "Fit to reconstructed toy events ";
        longString = longString + Form(parity==1?" (sub%i_q2-bin %i even)":" (sub%i_q2-bin %i odd)",is,q2Bin);

        // plot fit projections 
        c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to RECO-level toy MC - "+longString).c_str(),3600,1200);
        c[confIndex]->Divide(years.size());
        cout<<"plotting 4d canvas: "<< years.size() << endl;

        for (unsigned int iy = 0; iy < years.size(); iy++) {
            c[confIndex]->cd(iy+1);
            year.clear(); year.assign(Form("%i",years[iy]));
        
            std::vector<RooPlot*> frames;
            frames.push_back( prepareFrame( mass ->frame(Title((longString+year).c_str()))));
            TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);

            auto singleYearPdf = simPdf->getPdf(("data"+year+Form("_subs%d",is)).c_str());
            auto singleYearData = combData->reduce(("sample==sample::data"+year+Form("_subs%d",is)).c_str());
        
            for (unsigned int fr = 0; fr < frames.size(); fr++){
                cout << "year is " << year << endl;
                cout<<"fr "<<fr<<" data"<<year<<Form("_subs%d",is)<<endl;
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
		       Components( ("pdf_bkg_mass_"+year).c_str() ));

                singleYearPdf->plotOn(frames[fr],
		       // Slice(sample, ("data"+year+Form("_subs%d",nSample)).c_str()), 
		       // ProjWData(RooArgSet(sample), *combData), 
		       LineWidth(1), 
		       Name(("plPDFsig"+year).c_str()), 
		       NumCPU(4),
		       LineColor(880),
		       Components( ("pdf_sig_mass_"+year).c_str() ));

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
        c[confIndex]->SaveAs( ("plotMom_Fit/cocktail/simFitResult_cocktail_fullMass_" + plotString +  ".pdf").c_str() );
    }
    fout->cd();
    ws_pars->Write();
    momResultsTree->Write();
    parResultsTree->Write();
    fout->Close();
}


int main(int argc, char **argv) {

    int q2Bin = atoi(argv[1]);
    int iSample = atoi(argv[2]);
    std::vector<int> years;
    if ( argc > 3 && atoi(argv[3]) != 0 ) years.push_back(atoi(argv[3]));
    else {
        cout << "No specific years selected, using default: 2016" << endl;
        years.push_back(2016);
    }
    if ( argc > 4  && atoi(argv[4])  != 0 ) years.push_back(atoi(argv[4]));
    if ( argc > 5 && atoi(argv[5]) != 0 ) years.push_back(atoi(argv[5]));
    Recocock_MOM(q2Bin,years,iSample);
    return 0;
}




 