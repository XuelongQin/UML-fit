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

#include "RooDoubleCBFast.h"

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


void Gen_MOM(int q2Bin,std::vector<int> years){
    string all_years = "";
    string yearstring = ""; 

    for (unsigned int iy = 0; iy < years.size(); iy++) {
        yearstring.clear(); yearstring.assign(Form("%i",years[iy]));
        all_years += yearstring;
    }

    string outname = "./Momresults/GenMC_" ;
    outname = outname + all_years + Form("_b%i.root",q2Bin);
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

    TFile *Gen_file = new TFile(Form("/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/effDataset_b%i_2016.root",q2Bin));
    RooWorkspace *wsp_gen = (RooWorkspace*)Gen_file->Get(Form("ws_b%ip0",q2Bin));
    RooDataSet* data_gen = (RooDataSet*)wsp_gen->data(Form("data_genDen_ev_b%i",q2Bin));

    int totentries = data_gen->numEntries();
    cout << "total entries in gen data is " << totentries << endl;
    //break;

    for (int i=0;i<totentries;i++){
        cos_theta_k = data_gen->get(i)->getRealValue("ctK");
        cos_theta_l = data_gen->get(i)->getRealValue("ctL");
        phi_kst_mumu = data_gen->get(i)->getRealValue("phi");
        weight = data_gen->weight();            
        f_1s = 1 - cos_theta_k * cos_theta_k;
        M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
        M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
        f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
        f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
        f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
        f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
        f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
        f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);	                

        f1s->Fill(f_1s,weight);
        f3->Fill(f_3,weight);
        f4->Fill(f_4,weight);
        f7->Fill(f_7,weight);

        m6s->Fill(M_6s, weight);
        m6c->Fill(M_6c, weight);
        f5->Fill(f_5, weight);
        f8->Fill(f_8, weight);
        f9->Fill(f_9,weight);
    }

    fout->cd();
    f1s->Write();
    m6s->Write();
    m6c->Write();
    f3->Write();
    f4->Write();
    f5->Write();
    f7->Write();
    f8->Write();
    f9->Write();

    auto resulttuple = CalValue(f1s,m6s,m6c,f3,f4,f5,f7,f8,f9);
    //cout << "haha1" << endl;
    vResult.assign(get<0>(resulttuple).begin(),get<0>(resulttuple).end());
    vError.assign(get<1>(resulttuple).begin(),get<1>(resulttuple).end());
    momResultsTree->Fill();
    //fout->cd();
    //momResultsTree->Write();
    vPResult.assign(get<2>(resulttuple).begin(),get<2>(resulttuple).end());
    vPError.assign(get<3>(resulttuple).begin(),get<3>(resulttuple).end());
    parResultsTree->Fill();
    fout->cd();
    momResultsTree->Write();
    parResultsTree->Write();
    fout->Close();
}


int main(int argc, char **argv) {

    int q2Bin = atoi(argv[1]);
    std::vector<int> years;
    if ( argc > 2 && atoi(argv[2]) != 2 ) years.push_back(atoi(argv[2]));
    else {
        cout << "No specific years selected, using default: 2016" << endl;
        years.push_back(2016);
    }
    if ( argc > 3  && atoi(argv[3])  != 0 ) years.push_back(atoi(argv[3]));
    if ( argc > 4 && atoi(argv[4]) != 0 ) years.push_back(atoi(argv[4]));
    Gen_MOM(q2Bin,years);
    return 0;
}


