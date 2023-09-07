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

void Reco_MOM(int q2Bin,std::vector<int> years,int tag){
    //tag:
    //0: single mistag events
    //1: single correct tag events
    //2: mistag+correct tag events

    string all_years = "";
    string yearstring = ""; 

    for (unsigned int iy = 0; iy < years.size(); iy++) {
        yearstring.clear(); yearstring.assign(Form("%i",years[iy]));
        all_years += yearstring;
    }

    string outname = "./Momresults/RecoMC_" ;
    outname = outname + all_years + Form("_b%i_tag%i.root",q2Bin,tag);
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

        TH1D *m6sdiffsignweight;
        TH1D *m6cdiffsignweight;
        TH1D *f5diffsignweight;
        TH1D *f8diffsignweight;
        TH1D *f9diffsignweight;
        
        if (tag==0){
            f1sweight = (TH1D*)efffile->Get("f1sWTweight");
            m6sweight = (TH1D*)efffile->Get("m6sWTweight");
            m6cweight = (TH1D*)efffile->Get("m6cWTweight");
            f3weight = (TH1D*)efffile->Get("f3WTweight");
            f4weight = (TH1D*)efffile->Get("f4WTweight");
            f5weight = (TH1D*)efffile->Get("f5WTweight");
            f7weight = (TH1D*)efffile->Get("f7WTweight");
            f8weight = (TH1D*)efffile->Get("f8WTweight");
            f9weight = (TH1D*)efffile->Get("f9WTweight");
        }

        else if (tag==1){
            f1sweight = (TH1D*)efffile->Get("f1sCTweight");
            m6sweight = (TH1D*)efffile->Get("m6sCTweight");
            m6cweight = (TH1D*)efffile->Get("m6cCTweight");
            f3weight = (TH1D*)efffile->Get("f3CTweight");
            f4weight = (TH1D*)efffile->Get("f4CTweight");
            f5weight = (TH1D*)efffile->Get("f5CTweight");
            f7weight = (TH1D*)efffile->Get("f7CTweight");
            f8weight = (TH1D*)efffile->Get("f8CTweight");
            f9weight = (TH1D*)efffile->Get("f9CTweight");
        }

        else if (tag==2){
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

        }
        else{
            cout << "wrong tag configuration!" << endl;
            exit(0);
        }

        cout << "Read efficiency successfully" << endl;

        TFile *Reco_file = new TFile(Form("/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/MC-datasets/recoMCDataset_b%i_%i_XGBv8.root",q2Bin,year));
        RooWorkspace *wsp_reco = (RooWorkspace*)Reco_file->Get(Form("ws_b%ip0",q2Bin));
        RooDataSet* data_reco;
        if (tag==0){
            data_reco = (RooDataSet*)wsp_reco->data(Form("data_wtRECO_ev_b%i",q2Bin));
        }
        else if (tag==1){
            data_reco = (RooDataSet*)wsp_reco->data(Form("data_ctRECO_ev_b%i",q2Bin));
        }
        else if (tag==2){
            data_reco = (RooDataSet*)wsp_reco->data(Form("data_ctRECO_ev_b%i",q2Bin));
            RooDataSet *data_reco1 = (RooDataSet*)wsp_reco->data(Form("data_wtRECO_ev_b%i",q2Bin));
            data_reco->append(*data_reco1);
        }
        else{
            cout << "wrong tag configuration!" << endl;
            exit(0);
        }
        int totentries = data_reco->numEntries();
        cout << "total entries in total data is " << totentries << endl;

        for (int i=0;i<totentries;i++){

            cos_theta_k = data_reco->get(i)->getRealValue("ctK");
            cos_theta_l = data_reco->get(i)->getRealValue("ctL");
            phi_kst_mumu = data_reco->get(i)->getRealValue("phi");
            weight = data_reco->weight();            
            f_1s = 1 - cos_theta_k * cos_theta_k;
            M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
            M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
            f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
            f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);	                

            //cout<< "xixi" << endl;
            if (tag==0){
                M_6s = -M_6s;
                M_6c = -M_6c;
                f_5 = -f_5;
                f_8 = -f_8;
                f_9 = -f_9;
            }

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


            if (tag==2){
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
    int tag = atoi(argv[2]);
    std::vector<int> years;
    if ( argc > 3 && atoi(argv[3]) != 0 ) years.push_back(atoi(argv[3]));
    else {
        cout << "No specific years selected, using default: 2016" << endl;
        years.push_back(2016);
    }
    if ( argc > 4  && atoi(argv[4])  != 0 ) years.push_back(atoi(argv[4]));
    if ( argc > 5 && atoi(argv[5]) != 0 ) years.push_back(atoi(argv[5]));
    Reco_MOM(q2Bin,years,tag);
    return 0;
}
