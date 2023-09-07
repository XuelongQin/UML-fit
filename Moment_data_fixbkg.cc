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
#include "RooBernstein.h"

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
double f_FS;
double f_AS;
double f_AS4;
double f_AS5;
double f_AS7;
double f_AS8;



string paralist[21] = {"FlS", "AFBS", "S3S", "S4S", "S5S", "S7S", "S8S", "S9S", "P1S", "P2S", "P3S", "P4pS", "P5pS", "P6pS", "P8pS", "FS", "AS", "AS4", "AS5", "AS7", "AS8"};
string momlist[15] = {"M1s","M6s","M6c","M3","M4","M5","M7","M8","M9","MFS","MAS","MAS4","MAS5","MAS7","MAS8"};
vector<double> vResult(15);
vector<double> vError(15);
vector<double> vPResult(21);
vector<double> vPError(21);

static const int nBins = 9;
TCanvas* c [4*nBins];


tuple<RooDataSet*,RooDataSet*,RooDataSet*> splot(RooWorkspace *ws_pdf ,int q2Bin, int year, int isample){
    RooAbsPdf *mass_pdf = ws_pdf->pdf(Form("PDF_mass_b%ip1_%i",q2Bin,year));
    RooRealVar *nsig = ws_pdf->var(Form("nsigb%ip1_%i",q2Bin,year));
    RooRealVar *nbkg = ws_pdf->var(Form("nbkgb%ip1_%i",q2Bin,year));  
    cout << "nsig is " << nsig->getVal() << endl;
    cout << "nbkg is " << nbkg->getVal() << endl; 
	//mass_pdf->Print();
    RooDataSet *data = (RooDataSet*)ws_pdf->data(Form("data%i_b%ip1_subs%i",year,q2Bin,isample));
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
tuple<TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*> Constructhisto(int q2Bin){
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
    TH1D *fFSeffW = (TH1D*)efffile->Get("fFSWTeff");
    TH1D *fASeffW = (TH1D*)efffile->Get("fASWTeff");
    TH1D *fAS4effW = (TH1D*)efffile->Get("fAS4WTeff");
    TH1D *fAS5effW = (TH1D*)efffile->Get("fAS5WTeff");
    TH1D *fAS7effW = (TH1D*)efffile->Get("fAS7WTeff");
    TH1D *fAS8effW = (TH1D*)efffile->Get("fAS8WTeff");


	int bins = f1seffW->GetNbinsX();
	//TAxis *axis = m6seffW->GetXaxis();
	double xbinsf1s[bins+1], xbinsm6s[bins+1], xbinsm6c[bins+1], xbinsf3[bins+1], xbinsf4[bins+1], xbinsf5[bins+1], xbinsf7[bins+1], xbinsf8[bins+1], xbinsf9[bins+1], xbinsfFS[bins+1],xbinsfAS[bins+1],xbinsfAS4[bins+1],xbinsfAS5[bins+1],xbinsfAS7[bins+1], xbinsfAS8[bins+1] ;
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
        xbinsfFS[i] = fFSeffW->GetBinLowEdge(i+1);
        xbinsfAS[i] = fASeffW->GetBinLowEdge(i+1);
        xbinsfAS4[i] = fAS4effW->GetBinLowEdge(i+1);
        xbinsfAS5[i] = fAS5effW->GetBinLowEdge(i+1);
        xbinsfAS7[i] = fAS7effW->GetBinLowEdge(i+1);
        xbinsfAS8[i] = fAS8effW->GetBinLowEdge(i+1);

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
    TH1D *fFS = new TH1D("fFS","",bins,xbinsfFS);
    TH1D *fAS = new TH1D("fAS","",bins,xbinsfAS);
    TH1D *fAS4 = new TH1D("fAS4","",bins,xbinsfAS4);
    TH1D *fAS5 = new TH1D("fAS5","",bins,xbinsfAS5);
    TH1D *fAS7 = new TH1D("fAS7","",bins,xbinsfAS7);
    TH1D *fAS8 = new TH1D("fAS8","",bins,xbinsfAS8);

    return make_tuple(f1s,m6s,m6c,f3,f4,f5,f7,f8,f9,fFS,fAS,fAS4,fAS5,fAS7,fAS8);
}

tuple<double,double> MoMValue(TH1D *his){
    double mvalue = his->GetMean();
    double merror = his->GetRMS()/TMath::Sqrt(his->GetEntries()-1);
    return make_tuple(mvalue,merror);
}


tuple<vector<double>,vector<double>,vector<double>,vector<double>> CalValue(TH1D *f1s, TH1D *m6s, TH1D *m6c, TH1D *f3, TH1D *f4, TH1D *f5, TH1D *f7, TH1D *f8, TH1D *f9, TH1D *fFS, TH1D *fAS, TH1D *fAS4, TH1D *fAS5, TH1D* fAS7, TH1D* fAS8){
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
    auto MFStuple = MoMValue(fFS);
    auto MAStuple = MoMValue(fAS);
    auto MAS4tuple = MoMValue(fAS4);
    auto MAS5tuple = MoMValue(fAS5);
    auto MAS7tuple = MoMValue(fAS7);
    auto MAS8tuple = MoMValue(fAS8);

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

    momvalue[9] = get<0>(MFStuple);
    momerror[9] = get<1>(MFStuple);

    momvalue[10] = get<0>(MAStuple);
    momerror[10] = get<1>(MAStuple);

    momvalue[11] = get<0>(MAS4tuple);
    momerror[11] = get<1>(MAS4tuple);

    momvalue[12] = get<0>(MAS5tuple);
    momerror[12] = get<1>(MAS5tuple);

    momvalue[13] = get<0>(MAS7tuple);
    momerror[13] = get<1>(MAS7tuple);

    momvalue[14] = get<0>(MAS8tuple);
    momerror[14] = get<1>(MAS8tuple);

    vector<double> parvalue(21);
    vector<double> parerror(21);

//Fl
    parvalue[0] = (10.0*momvalue[9]+15.0*momvalue[0]-18.0)/(30.0*momvalue[9]+15.0*momvalue[0]-34.0);
    parerror[0] = TMath::Sqrt( TMath::Power((60.0*(-4.0+5.0*momvalue[9])/TMath::Power((-34.0+15.0*momvalue[0]+30.0*momvalue[9]),2)*momerror[9]),2) + TMath::Power((100.0*(-2.0+3.0*momvalue[0])/TMath::Power((-34.0+15.0*momvalue[0]+30.0*momvalue[9]),2)*momerror[0]),2));

//FS
    parvalue[15] = 15.0/4.0*(momvalue[0]+2.0*momvalue[9]-2.0);
    parerror[15] = 15.0/4.0*TMath::Sqrt(TMath::Power(momerror[0],2)+4*TMath::Power(momerror[9],2));

//AFBS
    parvalue[1] = 3.0/4.0*5.0/2.0*momvalue[1]/(1.0-parvalue[15]);;
    parerror[1] = 3.0/4.0*5.0/2.0*TMath::Sqrt(TMath::Power((1.0/(1.0-parvalue[15])*momerror[1]),2)+TMath::Power((parerror[15]*momvalue[1]/TMath::Power(1-parvalue[15],2)),2));

//S3S
    parvalue[2] = 25.0/8.0*momvalue[3]/(1-parvalue[15]);
    parerror[2] = 25.0/8.0*TMath::Sqrt(TMath::Power((1.0/(1.0-parvalue[15])*momerror[3]),2)+TMath::Power((parerror[15]*momvalue[3]/TMath::Power(1.0-parvalue[15],2)),2));

//S4S
    parvalue[3] = 25.0/8.0*momvalue[4]/(1-parvalue[15]);
    parerror[3] = 25.0/8.0*TMath::Sqrt(TMath::Power((1.0/(1.0-parvalue[15])*momerror[4]),2)+TMath::Power((parerror[15]*momvalue[4]/TMath::Power(1.0-parvalue[15],2)),2));

//S5S
    parvalue[4] = 5.0/2.0*momvalue[5]/(1-parvalue[15]);
    parerror[4] = 5.0/2.0*TMath::Sqrt(TMath::Power((1.0/(1.0-parvalue[15])*momerror[5]),2)+TMath::Power((parerror[15]*momvalue[5]/TMath::Power(1.0-parvalue[15],2)),2));   

//S7S
    parvalue[5] = 5.0/2.0*momvalue[6]/(1-parvalue[15]);
    parerror[5] = 5.0/2.0*TMath::Sqrt(TMath::Power((1.0/(1.0-parvalue[15])*momerror[6]),2)+TMath::Power((parerror[15]*momvalue[6]/TMath::Power(1.0-parvalue[15],2)),2));

//S8S
    parvalue[6] = 25.0/8.0*momvalue[7]/(1-parvalue[15]);
    parerror[6] = 25.0/8.0*TMath::Sqrt(TMath::Power((1.0/(1.0-parvalue[15])*momerror[7]),2)+TMath::Power((parerror[15]*momvalue[7]/TMath::Power(1.0-parvalue[15],2)),2));

//S9S
    parvalue[7] = 25.0/8.0*momvalue[8]/(1-parvalue[15]);
    parerror[7] = 25.0/8.0*TMath::Sqrt(TMath::Power((1.0/(1.0-parvalue[15])*momerror[8]),2)+TMath::Power((parerror[15]*momvalue[8]/TMath::Power(1.0-parvalue[15],2)),2));

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

//AS
    parvalue[16] = 15.0/4.0 * momvalue[10];
    parerror[16] = 15.0/4.0 * momerror[10];

//AS4
    parvalue[17] = 15.0/4.0 * momvalue[11];
    parerror[17] = 15.0/4.0 * momerror[11];

//AS5
    parvalue[18] = 3.0 * momvalue[12];
    parerror[18] = 3.0 * momerror[12];

//AS7
    parvalue[19] = 3.0 * momvalue[13];
    parerror[19] = 3.0 * momerror[13];

//AS8
    parvalue[20] = 15.0/4.0 * momvalue[14];
    parerror[20] = 15.0/4.0 * momerror[14];    

    return make_tuple(momvalue,momerror,parvalue,parerror);
}

void data_MOM(int q2Bin,std::vector<int> years){
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
    int parity = 1;
    string all_years = "";
    string yearstring = ""; 
    string XGBstr = "";
    int nSample=0;
    int q2stat=0;

    bool localFiles=false;

    int XGBv = 8;
    if (XGBv>9) XGBstr = Form("-TMVAv%i",XGBv-10);
    else if (XGBv>0) XGBstr = Form("-XGBv%i",XGBv);
    RooArgSet c_vars_rt, c_pdfs_rt;
    RooArgSet c_vars_wt, c_pdfs_wt;
    RooArgSet c_vars; 
    RooWorkspace * ws_pars = new RooWorkspace("ws_pars");

    string shortString = Form("b%ip%i",q2Bin,parity);
    cout<<"Conf: "<<shortString<<endl;

    std::vector<RooWorkspace*> wsp, wsp_mcmass,wsp_sb;
    std::vector<std::vector<RooDataSet*>> data;
    double power = 1.0;
    int fitOption = 1;
    gInterpreter->GenerateDictionary("std::pair<std::string, RooDataSet*>", "map;string;RooDataSet.h");
    gInterpreter->GenerateDictionary("std::map<std::string, RooDataSet*>",  "map;string;RooDataSet.h");
    gInterpreter->GenerateDictionary("std::pair<std::map<string,RooDataSet*>::iterator, bool>", "map;string;RooDataSet.h");
    std::map<std::string, RooDataSet*> map;

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

    uint iSample = 0;
    uint firstSample = iSample;
    uint lastSample = iSample;
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

    string outname = "./Momresults/data/data_fixbkg_" ;
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


    for (unsigned int iy = 0; iy < years.size(); iy++){
        year.clear(); year.assign(Form("%i",years[iy]));

        string filename_data = Form("recoDATADataset_b%i_%i.root",q2Bin,years[iy]);
        filename_data = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/data-datasets/" + filename_data;
        TFile *inputdata = new TFile(filename_data.c_str(),"read");
        RooWorkspace *wsp_data = (RooWorkspace*)inputdata->Get(Form("ws_b%ip0",q2Bin));

        std::vector<RooDataSet*> singyeardata;
        for (uint is = firstSample; is <= lastSample; is++) {
            RooDataSet* datasub;
            RooDataSet* singlesample = new RooDataSet((Form("data%i_",years[iy])+shortString + "_subs0").c_str(), 
				       ((Form("data%i_",years[iy])+shortString + "_subs0").c_str()), 
				       RooArgSet(observables));

            datasub = (RooDataSet*)wsp_data->data(Form("data_b%i",q2Bin))->reduce("mass>5.0 && mass <5.6");
            singlesample->append(*datasub);
            singyeardata.push_back (singlesample);
        }
        data.push_back(singyeardata);

        //// Signal part of the pdf ////
        // Mass Component - import mass PDF from fits to the MC
        string filename_mc_mass = Form("/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/massFits-%s%s/%i.root",q2Bin==4?"Jpsi":(q2Bin==6?"Psi":"LMNR"),XGBstr.c_str(), years[iy]);
        if (localFiles)
        filename_mc_mass = Form("results_fits_%i_fM_%s%s.root",years[iy], q2Bin==4?"Jpsi":(q2Bin==6?"Psi":"lmnr"), XGBv>0 ? Form("_MCw_xgbv%i",XGBv):"");
        if (!retrieveWorkspace( filename_mc_mass, wsp_mcmass, "w")) {
        cout<<"Workspace not found in mass-fit file: "<<filename_mc_mass<<endl;
        return;
        }

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
        if (q2Bin < 7) {
            if (nSample==0)  dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, alpha_rt2, n_rt1, n_rt2 ,f1rt, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt );
            else             dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, alpha_rt2, n_rt1, n_rt2 ,f1rt, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt, q2stat );
        } else {
            if (nSample==0)  dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, n_rt1 ,f1rt, q2Bin, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt );
            else             dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, n_rt1 ,f1rt, q2Bin, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt, q2stat );
        }
        } 
        else{
            alpha_rt2->setRange(0,10);
            if (nSample==0) dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, alpha_rt1, alpha_rt2, n_rt1, n_rt2 , wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt  );
            else dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, alpha_rt1, alpha_rt2, n_rt1, n_rt2 , wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt, q2stat );
        }
        
        // create constrained PDF for RT mass
        RooArgList constr_rt_list = RooArgList(c_pdfs_rt);
        constr_rt_list.add(*dcb_rt);
        RooProdPdf * c_dcb_rt = new RooProdPdf(("c_dcb_rt_"+year).c_str(), ("c_dcb_rt_"+year).c_str(), constr_rt_list );
        c_vars.add(c_vars_rt);
    
        /// create WT component
        wsp_mcmass[iy]->loadSnapshot(Form("reference_fit_WT_%i",q2Bin));

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
        RooRealVar* mean_wt = (RooRealVar*) mWT_data;
        cout << "deltaPeak constraint built, of value " << deltaPeakValue << " +/- " << deltaPeakError << endl;

        c_pdfs_wt.add(*c_deltaPeaks);
        c_vars_wt.add(*deltaPeakVar);

        if (nSample==0) dcb_wt = createWTMassShape(q2Bin, mass, mean_wt, sigma_wt, alpha_wt1, alpha_wt2, n_wt1, n_wt2 , wsp_mcmass[iy], years[iy], true, c_vars_wt, c_pdfs_wt );
        else dcb_wt = createWTMassShape(q2Bin, mass, mean_wt, sigma_wt, alpha_wt1, alpha_wt2, n_wt1, n_wt2 , wsp_mcmass[iy], years[iy], true, c_vars_wt, c_pdfs_wt, q2stat );

        // create constrained PDF for WT mass
        RooArgList constr_wt_list = RooArgList(c_pdfs_wt);
        constr_wt_list.add(*dcb_wt);
        RooProdPdf * c_dcb_wt = new RooProdPdf(("c_dcb_wt_"+year).c_str(), ("c_dcb_wt_"+year).c_str(), constr_wt_list );
        c_vars.add(c_vars_wt);


        RooRealVar* mFrac = new RooRealVar(Form("mFrac^{%i}",years[iy]),"mistag fraction",0.13, 0, 1);

        /*if (q2Bin < 5){
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
        */
        /// create constraint on mFrac (here there is no efficiency, therefore value set to measured value on MC)
        double nrt_mc   =  wsp_mcmass[iy]->var(Form("nRT_%i",q2Bin))->getVal(); 
        double nwt_mc   =  wsp_mcmass[iy]->var(Form("nWT_%i",q2Bin))->getVal(); 
        double fraction = nwt_mc / (nrt_mc + nwt_mc);
        //double fraction = nwt_mc/nrt_mc;
        //double frac_sigma = fM_sigmas[years[iy]][q2Bin]*(nrt_mc + nwt_mc)/nrt_mc;
        double frac_sigma = fM_sigmas[years[iy]][q2Bin];
        c_fm.push_back(new RooGaussian(Form("c_fm^{%i}",years[iy]) , "c_fm" , *mFrac,  
                                        RooConst(fraction) , 
                                        RooConst(frac_sigma))
                                        );
        cout << fraction << "   " << frac_sigma << endl;                                    
        c_vars.add(*mFrac);  

        PDF_sig_mass.push_back(new RooAddPdf(("signal_PDF_"+year).c_str(),
				 ("signal_PDF_"+year).c_str(),
				  RooArgList(*c_dcb_wt,*c_dcb_rt), *mFrac));
        pdf_sig_mass.push_back(new RooProdPdf(("pdf_sig_mass_"+year).c_str(), 
                                            ("pdf_sig_mass_"+year).c_str(),
                                            RooArgList(*PDF_sig_mass[iy],*c_fm[iy])));
        //c_vars.add(*c_fm[iy]);

        //// Background components ////
        // Read angular pdf for sidebands from external file 
        string filename_sb = Form("savesb_%i_b%i.root", years[iy], q2Bin);
        if (!localFiles) filename_sb = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/sidebands/" + filename_sb;
        if (!(retrieveWorkspace(filename_sb, wsp_sb, "wsb"))) return;

    // create exponential mass pdf for background
        // in Jpsi bin, if fitOption > 0 -> bernstein polynbomial
        //RooRealVar* slope = new RooRealVar(Form("slope^{%i}",years[iy]),  Form("slope^{%i}",years[iy]) , wsp_sb[iy]->var("slope")->getVal(), -10., 0.);
        RooRealVar* slope = new RooRealVar(Form("slope^{%i}",years[iy]),  Form("slope^{%i}",years[iy]) ,-5., -10., 0.);
        RooAbsPdf* bkg_mass_pdf = 0;
        double pol_bmax =1.;
        RooRealVar* b0_bkg_mass = new RooRealVar(Form("b0_bkg_mass-%i",years[iy]) , Form("b0_bkg_mass-%i",years[iy])  ,  pol_bmax  );
        RooRealVar* b1_bkg_mass = new RooRealVar(Form("b1_bkg_mass-%i",years[iy]) , Form("b1_bkg_mass-%i",years[iy])  ,  0.1,  0., pol_bmax);
        RooRealVar* b2_bkg_mass = new RooRealVar(Form("b2_bkg_mass-%i",years[iy]) , Form("b2_bkg_mass-%i",years[iy])  ,  0.1,  0., pol_bmax);
        RooRealVar* b3_bkg_mass = new RooRealVar(Form("b3_bkg_mass-%i",years[iy]) , Form("b3_bkg_mass-%i",years[iy])  ,  0.0 );
        RooRealVar* b4_bkg_mass = new RooRealVar(Form("b4_bkg_mass-%i",years[iy]) , Form("b4_bkg_mass-%i",years[iy])  ,  0.1 , 0., pol_bmax);
        if (q2Bin==4 && fitOption>0){
            b0_bkg_mass->setConstant(true);
            b3_bkg_mass->setConstant(true);
            if (iy==0){
                b1_bkg_mass->setVal(0.21402);
                b2_bkg_mass->setVal(0.39604);
                b4_bkg_mass->setVal(0.030159);
            }
            else if (iy==1){
                b1_bkg_mass->setVal(0.25924);
                b2_bkg_mass->setVal(0.37700);
                b4_bkg_mass->setVal(0.022907);
            }
            else {
                b1_bkg_mass->setVal(0.27240);
                b2_bkg_mass->setVal(0.40341);
                b4_bkg_mass->setVal(0.025352);
            }
            b1_bkg_mass->setConstant(true);
            b2_bkg_mass->setConstant(true);
            b4_bkg_mass->setConstant(true);
            bkg_mass_pdf = new RooBernstein(Form("bkg_mass1_%i",years[iy]),Form("bkg_mass1_%i",years[iy]),  *mass, RooArgList(*b0_bkg_mass, *b1_bkg_mass, *b2_bkg_mass, *b3_bkg_mass,* b4_bkg_mass));
        } else { 
            if (q2Bin==6){
                if (iy==0){
                    slope->setVal(-5.9601);
                }
                else if (iy==1){
                    slope->setVal(-5.8669);
                }
                else{
                    slope->setVal(-5.9322);
                }
            }
            slope->setConstant(true);
            bkg_mass_pdf = new RooExponential(Form("bkg_mass1_%i",years[iy]),  Form("bkg_mass1_%i",years[iy]) ,  *mass,   *slope  );
        } 
        pdf_bkg_mass.push_back(bkg_mass_pdf);


        // RooRealVar* slope       = new RooRealVar    (Form("slope^{%i}",years[iy]),  Form("slope^{%i}",years[iy]) , wsp_sb[iy]->var("slope")->getVal(), -10., 0.);

        double fullentries = data[iy][0]->sumEntries();
        cout << "number of fullentries is " << fullentries << endl;
        double nbkg_ex = 0.3 * fullentries;
        double nsig_ex = 0.7 * fullentries;
        RooRealVar *nsig = new RooRealVar(("nsig"+shortString+"_"+year).c_str(),"number of signal",nsig_ex,0,0.8*fullentries);
        RooRealVar *nbkg = new RooRealVar(("nbkg"+shortString+"_"+year).c_str(),"number of bkg",nbkg_ex,0,0.4*fullentries);
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
                                                RooFit::Extended(kTRUE),
                                                RooFit::Constrain(c_vars),
                                                RooFit::NumCPU(8)
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
        //m.minos();
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
        TH1D *fFS = get<9>(histo);
        TH1D *fAS = get<10>(histo);
        TH1D *fAS4 = get<11>(histo);
        TH1D *fAS5 = get<12>(histo);
        TH1D *fAS7 = get<13>(histo);
        TH1D *fAS8 = get<14>(histo);

        for (unsigned int iy = 0; iy < years.size(); iy++){
            int year = years[iy];
            cout << "\n[Moment::ReCalValue]\tTry to fill year" << year << "reco mc sample" << endl; 
            stringstream myString;
            myString.clear();
            myString.str("");
            //myString << "/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/moment/1Deff/momentfinal/eff/1Deffb" << q2Bin << "y"<< year <<".root" ;
            myString << "/afs/cern.ch/user/x/xuqin/eos/Kmumu/momeff/1Deffb" << q2Bin << "y" << year << "_" << "od" << ".root" ;
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
            TH1D *fFSweight;
            TH1D *fASweight;
            TH1D *fAS4weight;
            TH1D *fAS5weight;
            TH1D *fAS7weight;
            TH1D *fAS8weight;

            TH1D *m6sdiffsignweight;
            TH1D *m6cdiffsignweight;
            TH1D *f5diffsignweight;
            TH1D *f8diffsignweight;
            TH1D *f9diffsignweight;
            TH1D *fASdiffsignweight;
            TH1D *fAS4diffsignweight;
            TH1D *fAS7diffsignweight;

            f1sweight = (TH1D*)efffile->Get("f1stotweight");
            m6sweight = (TH1D*)efffile->Get("m6ssamesignweight");
            m6cweight = (TH1D*)efffile->Get("m6csamesignweight");
            f3weight = (TH1D*)efffile->Get("f3totweight");
            f4weight = (TH1D*)efffile->Get("f4totweight");
            f5weight = (TH1D*)efffile->Get("f5samesignweight");
            f7weight = (TH1D*)efffile->Get("f7totweight");
            f8weight = (TH1D*)efffile->Get("f8samesignweight");
            f9weight = (TH1D*)efffile->Get("f9samesignweight");
            fFSweight = (TH1D*)efffile->Get("fFStotweight");
            fASweight = (TH1D*)efffile->Get("fASsamesignweight");
            fAS4weight = (TH1D*)efffile->Get("fAS4samesignweight");
            fAS5weight = (TH1D*)efffile->Get("fAS5totweight");
            fAS7weight = (TH1D*)efffile->Get("fAS7samesignweight");
            fAS8weight = (TH1D*)efffile->Get("fAS8totweight");

            m6sdiffsignweight = (TH1D*)efffile->Get("m6sdiffsignweight");
            m6cdiffsignweight = (TH1D*)efffile->Get("m6cdiffsignweight");
            f5diffsignweight = (TH1D*)efffile->Get("f5diffsignweight");
            f8diffsignweight = (TH1D*)efffile->Get("f8diffsignweight");
            f9diffsignweight = (TH1D*)efffile->Get("f9diffsignweight");
            fASdiffsignweight = (TH1D*)efffile->Get("fASdiffsignweight");
            fAS4diffsignweight = (TH1D*)efffile->Get("fAS4diffsignweight");
            fAS7diffsignweight = (TH1D*)efffile->Get("fAS7diffsignweight");

            cout << "Read efficiency successfully" << endl;
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
                f_FS = 1 - cos_theta_l * cos_theta_l;
                f_AS = (1 - cos_theta_l * cos_theta_l) * cos_theta_k;
                f_AS4 = 2 * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l)) * cos_theta_l * TMath::Cos(phi_kst_mumu);
                f_AS5 = TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l)) * TMath::Cos(phi_kst_mumu);
                f_AS7 = TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l)) * TMath::Sin(phi_kst_mumu);
                f_AS8 = 2 * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l)) * cos_theta_l * TMath::Sin(phi_kst_mumu);

                double weightf1s,weightm6s,weightm6c,weightf3,weightf4,weightf5,weightf7,weightf8,weightf9,weightfFS,weightfAS,weightfAS4,weightfAS5,weightfAS7,weightfAS8;
                double weightdiffm6s, weightdiffm6c,weightdifff5,weightdifff8,weightdifff9,weightdifffAS,weightdifffAS4,weightdifffAS7;

                weightf1s = f1sweight->GetBinContent(f1sweight->GetXaxis()->FindBin(f_1s));
                weightm6s = m6sweight->GetBinContent(m6sweight->GetXaxis()->FindBin(M_6s));
                weightm6c = m6cweight->GetBinContent(m6cweight->GetXaxis()->FindBin(M_6c));
                weightf3 = f3weight->GetBinContent(f3weight->GetXaxis()->FindBin(f_3));
                weightf4 = f4weight->GetBinContent(f4weight->GetXaxis()->FindBin(f_4));
                weightf5 = f5weight->GetBinContent(f5weight->GetXaxis()->FindBin(f_5));
                weightf7 = f7weight->GetBinContent(f7weight->GetXaxis()->FindBin(f_7));
                weightf8 = f8weight->GetBinContent(f8weight->GetXaxis()->FindBin(f_8));
                weightf9 = f9weight->GetBinContent(f9weight->GetXaxis()->FindBin(f_9));
                weightfFS = fFSweight->GetBinContent(fFSweight->GetXaxis()->FindBin(f_FS));
                weightfAS = fASweight->GetBinContent(fASweight->GetXaxis()->FindBin(f_AS));
                weightfAS4 = fAS4weight->GetBinContent(fAS4weight->GetXaxis()->FindBin(f_AS4));
                weightfAS5 = fAS5weight->GetBinContent(fAS5weight->GetXaxis()->FindBin(f_AS5));
                weightfAS7 = fAS7weight->GetBinContent(fAS7weight->GetXaxis()->FindBin(f_AS7));
                weightfAS8 = fAS8weight->GetBinContent(fAS8weight->GetXaxis()->FindBin(f_AS8));
                f1s->Fill(f_1s,weight*weightf1s);
                f3->Fill(f_3,weight*weightf3);
                f4->Fill(f_4,weight*weightf4);
                f7->Fill(f_7,weight*weightf7);
                fFS->Fill(f_FS, weight*weightfFS);
                fAS5->Fill(f_AS5, weight*weightfAS5);
                fAS8->Fill(f_AS8, weight*weightfAS8);

                m6s->Fill(M_6s, weight*weightm6s);
                m6c->Fill(M_6c, weight*weightm6c);
                f5->Fill(f_5, weight*weightf5);
                f8->Fill(f_8, weight*weightf8);
                f9->Fill(f_9,weight*weightf9);
                fAS->Fill(f_AS, weight*weightfAS);
                fAS4->Fill(f_AS4,weight*weightfAS4);
                fAS7->Fill(f_AS7, weight*weightfAS7);

                weightdiffm6s = m6sdiffsignweight->GetBinContent(m6sdiffsignweight->GetXaxis()->FindBin(M_6s));
                weightdiffm6c = m6cdiffsignweight->GetBinContent(m6cdiffsignweight->GetXaxis()->FindBin(M_6c));
                weightdifff5 = f5diffsignweight->GetBinContent(f5diffsignweight->GetXaxis()->FindBin(f_5));
                weightdifff8 = f8diffsignweight->GetBinContent(f8diffsignweight->GetXaxis()->FindBin(f_8));
                weightdifff9 = f9diffsignweight->GetBinContent(f9diffsignweight->GetXaxis()->FindBin(f_9));
                weightdifffAS = fASdiffsignweight->GetBinContent(fASdiffsignweight->GetXaxis()->FindBin(f_AS));
                weightdifffAS4 = fAS4diffsignweight->GetBinContent(fAS4diffsignweight->GetXaxis()->FindBin(f_AS4));
                weightdifffAS7 = fAS7diffsignweight->GetBinContent(fAS7diffsignweight->GetXaxis()->FindBin(f_AS7));

                m6s->Fill(-M_6s, weight*weightdiffm6s);
                m6c->Fill(-M_6c, weight*weightdiffm6c);
                f5->Fill(-f_5, weight*weightdifff5);
                f8->Fill(-f_8, weight*weightdifff8);
                f9->Fill(-f_9,weight*weightdifff9);
                fAS->Fill(-f_AS, weight*weightdifffAS);
                fAS4->Fill(-f_AS4,weight*weightdifffAS4);
                fAS7->Fill(-f_AS7, weight*weightdifffAS7);
            }
        }
        auto resulttuple = CalValue(f1s,m6s,m6c,f3,f4,f5,f7,f8,f9,fFS,fAS,fAS4,fAS5,fAS7,fAS8);
        //cout << "haha1" << endl;
        vResult.assign(get<0>(resulttuple).begin(),get<0>(resulttuple).end());
        vError.assign(get<1>(resulttuple).begin(),get<1>(resulttuple).end());
        momResultsTree->Fill();
        //fout->cd();
        //momResultsTree->Write();
        vPResult.assign(get<2>(resulttuple).begin(),get<2>(resulttuple).end());
        vPError.assign(get<3>(resulttuple).begin(),get<3>(resulttuple).end());
        parResultsTree->Fill();

        string plotString = shortString + "_" + Form("subs%i_",is) + all_years;
        //if (nSample>0) plotString = plotString + Form("_s%i",nSample);

        int confIndex = 0;
        string longString  = "Fit to data events ";
        longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);

        // plot fit projections 
        c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to data - "+longString).c_str(),3600,1200);
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
		       NumCPU(4),
               Normalization(1.0, RooAbsReal::RelativeExpected));

                singleYearPdf->plotOn(frames[fr],
		       // Slice(sample, ("data"+year+Form("_subs%d",nSample)).c_str()), 
		       // ProjWData(RooArgSet(sample), *combData), 
		       LineWidth(1), 
		       Name(("plPDFbkg"+year).c_str()), 
		       NumCPU(4),
		       LineColor(8),
		       Components( ("bkg_mass1_"+year).c_str() ),
               Normalization(1.0, RooAbsReal::RelativeExpected));

                singleYearPdf->plotOn(frames[fr],
		       // Slice(sample, ("data"+year+Form("_subs%d",nSample)).c_str()), 
		       // ProjWData(RooArgSet(sample), *combData), 
		       LineWidth(1), 
		       Name(("plPDFsig"+year).c_str()), 
		       NumCPU(4),
		       LineColor(880),
		       Components( ("pdf_sig_mass_"+year).c_str() ),
               Normalization(1.0, RooAbsReal::RelativeExpected));

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
        c[confIndex]->SaveAs( ("plotMom_Fit/data/simFitResult_data_fullMass_fixbkg" + plotString +  ".pdf").c_str() );
    }
    fout->cd();
    ws_pars->Write();
    momResultsTree->Write();
    parResultsTree->Write();
    fout->Close();
}


int main(int argc, char **argv) {

    int q2Bin = atoi(argv[1]);
    std::vector<int> years;
    if ( argc > 2 && atoi(argv[2]) != 0 ) years.push_back(atoi(argv[2]));
    else {
        cout << "No specific years selected, using default: 2016" << endl;
        years.push_back(2016);
    }
    if ( argc > 3  && atoi(argv[3])  != 0 ) years.push_back(atoi(argv[3]));
    if ( argc > 4 && atoi(argv[4]) != 0 ) years.push_back(atoi(argv[4]));
    data_MOM(q2Bin,years);
    return 0;
}




