#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>
#include <list>
#include <map>
#include <TMatrixT.h>
#include <TMatrixTSym.h>

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
#include <RooProduct.h>
#include <RooCBShape.h>
#include "RooDoubleCBFast.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
#include "RooBernstein.h"

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

TCanvas* cnll;
TCanvas* cZoom;
TCanvas* cPen;
TCanvas* c [4*nBins];

double power = 1.0;

void Generate_toy(int q2Bin, int parity, bool multiSample, uint nSample, uint q2stat, int XGBv, bool localFiles, float Fsvalue , std::vector<int> years)
{

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
    int fitOption=1;
    string shortString = Form("b%ip%i",q2Bin,parity);
    cout<<"Conf: "<<shortString<<endl;

    // Load variables and dataset
    // importing the complementary dataset, to fit with statistically uncorrelated efficiency

    string effCString = Form("effCHist_b%ip%i",q2Bin,parity);
    string effWString = Form("effWHist_b%ip%i",q2Bin,parity);
    string intCHistString = "MCint_"+shortString + "t1";
    string intWHistString = "MCint_"+shortString + "t0";
    string all_years = "";
    string year = ""; 
    string isample = ""; 
    uint firstSample = ( multiSample || nSample==0 ) ? 0 : nSample-1;
    uint lastSample = nSample > 0 ? nSample-1 : 0;
    string stat = "";
    if (nSample>0)   stat = Form("_b%istat-%i",q2stat,firstSample);
    if (multiSample) stat = stat + Form("-%i",lastSample);

    string sigpdfname = "PDF_sig_ang_mass_"+shortString+"_%i";
    string bkgpdfname = "bkg_pdf_%i";

    string XGBstr = "";
    if (XGBv>9) XGBstr = Form("-TMVAv%i",XGBv-10);
    else if (XGBv>0) XGBstr = Form("-XGBv%i",XGBv);

    std::vector<TFile*> fin_eff;
    std::vector<RooWorkspace*> wsp, wsp_mcmass, wsp_sb, wsp_Z4430;
    std::vector<std::vector<RooDataSet*>> data;
    std::vector<RooAbsReal*> effC, effW;
    std::vector<TH3D*> effCHist, effWHist;
    std::vector<TH1D*> intCHist, intWHist;
    std::vector< std::vector<double> > intCVec(years.size(), std::vector<double>(0));
    std::vector< std::vector<double> > intWVec(years.size(), std::vector<double>(0));
    std::vector<RooAbsPdf*> PDF_sig_ang_mass_bkg(0);
    std::vector<RooAbsPdf*> PDF_sig_ang_mass_bkg_penalty(0);
    std::vector<RooGaussian*> c_deltaPeaks, c_fm;
    RooArgSet c_vars_rt, c_pdfs_rt;
    RooArgSet c_vars_wt, c_pdfs_wt;
    RooArgSet c_vars; 

    //// from https://root-forum.cern.ch/t/combining-roodatasets-using-std-map-in-pyroot/16471/20
    gInterpreter->GenerateDictionary("std::pair<std::string, RooDataSet*>", "map;string;RooDataSet.h");
    gInterpreter->GenerateDictionary("std::map<std::string, RooDataSet*>",  "map;string;RooDataSet.h");
    gInterpreter->GenerateDictionary("std::pair<std::map<string,RooDataSet*>::iterator, bool>", "map;string;RooDataSet.h");
    std::map<std::string, RooDataSet*> map;

    RooRealVar* ctK = new RooRealVar("ctK", "cos(#theta_{K})", -1  , 1  );
    RooRealVar* ctL = new RooRealVar("ctL", "cos(#theta_{l})", -1  , 1  );
    RooRealVar* phi = new RooRealVar("phi", "#phi", -3.14159, 3.14159  );
    RooArgList vars (* ctK,* ctL,* phi);
    RooRealVar* mass = new RooRealVar("mass","mass", 5.,5.6,"GeV");
    RooRealVar* rand = new RooRealVar("rand", "rand", 0,1);
    RooArgSet reco_vars (*ctK, *ctL, *phi, *mass, *rand);
    RooArgSet observables (*ctK, *ctL, *phi, *mass);

    double fsig_year[3] = {0.75,0.78,0.78};
    // define angular parameters with ranges from positiveness requirements on the decay rate
    RooRealVar* Fl    = new RooRealVar("Fl","F_{L}",0.5,0,1);
    RooRealVar* P1    = new RooRealVar("P1","P_{1}",0,-1,1);   
    RooRealVar* P2    = new RooRealVar("P2","P_{2}",0,-0.5,0.5);
    RooRealVar* P3    = new RooRealVar("P3","P_{3}",0,-0.5,0.5);
    RooRealVar* P4p   = new RooRealVar("P4p","P'_{4}",0,-1*sqrt(2),sqrt(2));
    RooRealVar* P5p   = new RooRealVar("P5p","P'_{5}",0,-1*sqrt(2),sqrt(2));
    RooRealVar* P6p   = new RooRealVar("P6p","P'_{6}",0,-1*sqrt(2),sqrt(2));
    RooRealVar* P8p   = new RooRealVar("P8p","P'_{8}",0,-1*sqrt(2),sqrt(2));

    RooRealVar* Fs  = new RooRealVar("Fs","F_{S}",0.05,0,1);
    RooRealVar* As  = new RooRealVar("As","A^{S}",0,-1,1);
    RooRealVar* A4s = new RooRealVar("A4s","A_{4}^{S}",0,-1,1);
    RooRealVar* A5s = new RooRealVar("A5s","A_{5}^{S}",0,-1,1);
    RooRealVar* A7s = new RooRealVar("A7s","A_{7}^{S}",0,-1,1);
    RooRealVar* A8s = new RooRealVar("A8s","A_{8}^{S}",0,-1,1);

    Fl->setVal(0.58);
    P1->setVal(-0.04);
    P2->setVal(0.00);
    P3->setVal(0.26);
    P4p->setVal(-0.96);
    P5p->setVal(0.00);
    P6p->setVal(0.00);
    P8p->setVal(-0.24);

    //Fs value from fit: 0.085
    //Fs value from mom: 0.042
    Fs->setVal(Fsvalue);
    As->setVal(-0.26);
    A4s->setVal(0.06);
    A5s->setVal(0.01);
    A7s->setVal(0.00);
    A8s->setVal(-0.40);

    if (Fsvalue==0){
        As->setVal(0.00);
        A4s->setVal(0.00);
        A5s->setVal(0.00);
        A7s->setVal(0.00);
        A8s->setVal(0.00);
    }


    TRandom gen_nevt;

    string foutname = Form("/afs/cern.ch/user/x/xuqin/eos/Kmumu/mom_toy_Jpsi_Fs/controltoyDataset_Fs%i_q2st%i.root",int(100*Fsvalue),q2stat);
    TFile *fout = new TFile(foutname.c_str(),"recreate");
    RooWorkspace *wsp_out = new RooWorkspace("ws_toy","ws_toy");
    for (unsigned int iy = 0; iy < years.size(); iy++) {
        year.clear(); year.assign(Form("%i",years[iy]));

        // import KDE efficiency histograms and partial integral histograms
        string filename = Form((parity==0 ? "KDEeff_b%i_ev_%i.root" : "KDEeff_b%i_od_%i.root"),q2Bin,years[iy]);
        if (!localFiles) {
        if (XGBv<1) filename = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/eff/" + filename;
        else filename = Form("/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/eff-XGBv%i/",XGBv) + filename;
        }
        fin_eff.push_back( new TFile( filename.c_str(), "READ" ));
        if ( !fin_eff[iy] || !fin_eff[iy]->IsOpen() ) {
        cout<<"File not found: "<<filename<<endl;
        return;   
        }

        effCHist.push_back( (TH3D*)fin_eff[iy]->Get(effCString.c_str()));
        effWHist.push_back( (TH3D*)fin_eff[iy]->Get(effWString.c_str()));
        if ( !effCHist[iy] || effCHist[iy]->IsZombie() || !effWHist[iy] || effWHist[iy]->IsZombie() ) {
        cout<<"Efficiency histogram "<< effCString <<" or " << effWString << " not found in file: "<< filename <<endl;
        return;
        }

        // create efficiency functions
        RooDataHist* effCData = new RooDataHist(("effCData_"+shortString+"_"+year).c_str(),"effCData",vars,effCHist[iy]);
        RooDataHist* effWData = new RooDataHist(("effWData_"+shortString+"_"+year).c_str(),"effWData",vars,effWHist[iy]);
        effC.push_back( new RooHistFunc(("effC_"+shortString+"_"+year).c_str(),
                                        ("effC"+year).c_str() ,
                                        vars,
                                        *effCData,
                                        1));
        effW.push_back( new RooHistFunc(("effW_"+shortString+"_"+year).c_str(),
                                        ("effW"+year).c_str() ,
                                        vars,
                                        *effWData,
                                        1));

        // import precomputed integrals and fill a std::vector
        intCHist.push_back( (TH1D*)fin_eff[iy]->Get(intCHistString.c_str()));
        intWHist.push_back( (TH1D*)fin_eff[iy]->Get(intWHistString.c_str()));
        intCVec.push_back (vector<double> (0));
        intWVec.push_back (vector<double> (0));
        if ( !intCHist[iy] || intCHist[iy]->IsZombie() || !intWHist[iy] || intWHist[iy]->IsZombie() ) {
        cout << "Integral histogram " << intCHistString <<" or " << intWHistString << " not found in file: "<< filename << endl << "Abort" << endl;
        return;
        } else if ( strcmp( intCHist[iy]->GetTitle(), effCHist[iy]->GetTitle() ) || strcmp( intWHist[iy]->GetTitle(), effWHist[iy]->GetTitle() )) {
        // if the eff_config tag is different between efficiency and precomputed-integral means that they are inconsistent
        cout << "Integral histograms are incoherent with efficiency in file: " << filename << endl;
        cout << "Efficiency (CT) conf: " << effCHist[iy]->GetTitle() <<endl;
        cout << "Integral (CT) conf: "   << intCHist[iy]->GetTitle() <<endl;
        cout << "Efficiency (WT) conf: " << effWHist[iy]->GetTitle() <<endl;
        cout << "Integral (WT) conf: "   << intWHist[iy]->GetTitle() <<endl;
        cout << "Abort"<<endl;
        return;
        } 
        else {
        for (int i=1; i<=intCHist[iy]->GetNbinsX(); ++i) {
            intCVec[iy].push_back(intCHist[iy]->GetBinContent(i));
        }
        for (int i=1; i<=intWHist[iy]->GetNbinsX(); ++i) {
            intWVec[iy].push_back(intWHist[iy]->GetBinContent(i));
        }
        }

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

        // As the signal PDF is written as [ CT + mFrac * MT ] (see the PdfSigAngMass class),
        // the mFrac parameter represents the fraction between mistagged and correctly-tagged events (n_MT/n_CT)
        // Also, since the integral of the efficiencies contains the information on the mistag fraction in MC
        // this parameter represents the ratio between the fitted mFrac and the MC-based one ( n_MT_data / n_CT_data * n_CT_MC / n_MT_MC )
        RooRealVar* mFrac = new RooRealVar(Form("f_{M}^{%i}",years[iy]),"mistag fraction",1, 0, 15);
        // The values of fM_sigmas are computed on data-like MC subsamples as the fluctuation of the fraction of mistagged events ( n_MT/(n_MT+n_CT) )
        // this fluctuation needs to be propagated on the quantity of mFrac, to define a Gaussian contraint
        double nrt_mc   =  wsp_mcmass[iy]->var(Form("nRT_%i",q2Bin))->getVal(); 
        double nwt_mc   =  wsp_mcmass[iy]->var(Form("nWT_%i",q2Bin))->getVal(); 
        double fraction = nwt_mc / (nrt_mc + nwt_mc);
        // Propagation: sigma(mFrac) = sigma(n_MT/n_CT) * n_CT/n_MT
        //                           = sigma(fraction)/(1-fraction)^2 * (1-fraction)/fraction
        //                           = sigma(fraction) / fraction / (1-fraction)
        double frac_sigma = fM_sigmas[years[iy]][q2Bin]/fraction/(1-fraction);
        if (nSample!=0) frac_sigma = fM_sigmas[years[iy]][q2stat]/fraction/(1-fraction);

        RooGaussian* c_fm = new RooGaussian(Form("c_fm^{%i}",years[iy]) , "c_fm" , *mFrac,  
                                            RooConst(1.) , 
                                            RooConst(frac_sigma)
                                            );

        cout << "mFrac = " << fraction << " +/- " << fM_sigmas[years[iy]][q2Bin] << " ---> R = 1 +/- " << frac_sigma << endl;
        c_vars.add(*mFrac); 

        // Angular Component
        RooAbsReal* ang_rt = new ShapeSigAng(("PDF_sig_ang_rt_"+shortString+"_"+year).c_str(),
                                            ("PDF_sig_ang_rt_"+year).c_str(),
                                    *ctK,*ctL,*phi,
                                    *Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,
                        *Fs,*As,*A4s,*A5s,*A7s,*A8s,
                                    *effC[iy], intCVec[iy],
                                    true
                                    );
        
        RooAbsReal* ang_wt = new ShapeSigAng(("PDF_sig_ang_wt_"+shortString+"_"+year).c_str(),
                                            ("PDF_sig_ang_wt_"+year).c_str(),
                                    *ctK,*ctL,*phi,
                                    *Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,
                        *Fs,*As,*A4s,*A5s,*A7s,*A8s,
                                    *effW[iy], intWVec[iy],
                                    false 
                                    );
        
        // Angular * mass component
        PdfSigAngMass* pdf_sig_ang_mass = nullptr;
        //PdfSigAngMass* pdf_sig_ang_mass_penalty = nullptr;
        if (q2Bin < 5)  {
        pdf_sig_ang_mass = new PdfSigAngMass( Form(sigpdfname.c_str(),years[iy]),
                            Form(sigpdfname.c_str(),years[iy]),
                            *ctK,*ctL,*phi,*mass,
                            *mean_rt, *sigma_rt, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2,
                            *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
                            *mFrac, *c_fm,
                            *ang_rt, *ang_wt,
                            *c_dcb_rt, *c_dcb_wt
                            );
        
        /*pdf_sig_ang_mass_penalty = new PdfSigAngMass( Form((sigpdfname+"_penalty").c_str(),years[iy]),
                                Form((sigpdfname+"_penalty").c_str(),years[iy]),
                                *ctK,*ctL,*phi,*mass,
                                *mean_rt, *sigma_rt, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2,
                                *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
                                *mFrac, *c_fm,
                                *penTerm,
                                *ang_rt, *ang_wt,
                                *c_dcb_rt, *c_dcb_wt
                                );*/
        }      		                                           
        else if (q2Bin==7){
            pdf_sig_ang_mass = new PdfSigAngMass( Form(sigpdfname.c_str(),years[iy]),
                                                Form(sigpdfname.c_str(),years[iy]),
                                        *ctK,*ctL,*phi,*mass,
                                                *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *n_rt1, *f1rt,
                                                *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
                                        *mFrac, *c_fm,
                                        *ang_rt, *ang_wt,
                                        *c_dcb_rt, *c_dcb_wt
                                        );
        
            /*pdf_sig_ang_mass_penalty =  new PdfSigAngMass( Form((sigpdfname+"_penalty").c_str(),years[iy]),
                                                        Form((sigpdfname+"_penalty").c_str(),years[iy]),
                                                    *ctK,*ctL,*phi,*mass,
                                                        *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *n_rt1, *f1rt,
                                                        *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
                                                    *mFrac, *c_fm,
                                                *penTerm,
                                                    *ang_rt, *ang_wt,
                                                    *c_dcb_rt, *c_dcb_wt
                                                    );*/
        }      		                                           
        else {
        pdf_sig_ang_mass = new PdfSigAngMass( Form(sigpdfname.c_str(),years[iy]),
                            Form(sigpdfname.c_str(),years[iy]),
                            *ctK,*ctL,*phi,*mass,
                            *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2, *f1rt,
                            *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
                            *mFrac, *c_fm,
                            *ang_rt, *ang_wt,
                            *c_dcb_rt, *c_dcb_wt
                            );
        
        /*pdf_sig_ang_mass_penalty = new PdfSigAngMass( Form((sigpdfname+"_penalty").c_str(),years[iy]),
                                Form((sigpdfname+"_penalty").c_str(),years[iy]),
                                *ctK,*ctL,*phi,*mass,
                                *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2, *f1rt,
                                *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,                        
                                *mFrac, *c_fm,
                                *penTerm,
                                *ang_rt, *ang_wt,
                                *c_dcb_rt, *c_dcb_wt
                                );*/
        } 

        // Add constraint on mFrac to the pdf
        auto pdf_sig_ang_mass_mfc = new RooProdPdf(("PDF_sig_ang_mass_mfc_"+shortString+"_"+year).c_str(),
                            ("PDF_sig_ang_mass_mfc_"+year).c_str(),
                            *pdf_sig_ang_mass,
                            *c_fm);
        /*auto pdf_sig_ang_mass_penalty_mfc = new RooProdPdf(("PDF_sig_ang_mass_penalty_mfc_"+shortString+"_"+year).c_str(),
                            ("PDF_sig_ang_mass_penalty_mfc_"+year).c_str(),
                            *pdf_sig_ang_mass_penalty,
                            *c_fm);*/
        

        //// Background components ////
        // Read angular pdf for sidebands from external file 
        string filename_sb = Form("savesb_%i_b%i.root", years[iy], q2Bin);
        if (!localFiles) filename_sb = "/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/sidebands/" + filename_sb;
        if (!(retrieveWorkspace(filename_sb, wsp_sb, "wsb"))) return;

        RooBernsteinSideband* bkg_ang_pdf = (RooBernsteinSideband*) wsp_sb[iy]->pdf(Form("BernSideBand_bin%i_%i", q2Bin, years[iy]));
        RooArgSet* bkg_ang_params = (RooArgSet*)bkg_ang_pdf->getParameters(observables);
        auto iter = bkg_ang_params->createIterator();
        RooRealVar* ivar =  (RooRealVar*)iter->Next();
        while (ivar) {
            ivar->setConstant(true);
            ivar = (RooRealVar*) iter->Next();
        }

        // create exponential mass pdf for background
        // in Jpsi bin, if fitOption > 0 -> bernstein polynbomial
        RooRealVar* slope = new RooRealVar(Form("slope^{%i}",years[iy]),  Form("slope^{%i}",years[iy]) , -5., -10., 0.);
        RooAbsPdf* bkg_mass_pdf = 0;
        double pol_bmax =1.;
        RooRealVar* b0_bkg_mass = new RooRealVar(Form("b0_bkg_mass-%i",years[iy]) , Form("b0_bkg_mass-%i",years[iy])  ,  pol_bmax  );
        RooRealVar* b1_bkg_mass = new RooRealVar(Form("b1_bkg_mass-%i",years[iy]) , Form("b1_bkg_mass-%i",years[iy])  ,  0.1,  0., pol_bmax);
        RooRealVar* b2_bkg_mass = new RooRealVar(Form("b2_bkg_mass-%i",years[iy]) , Form("b2_bkg_mass-%i",years[iy])  ,  0.1,  0., pol_bmax);
        RooRealVar* b3_bkg_mass = new RooRealVar(Form("b3_bkg_mass-%i",years[iy]) , Form("b3_bkg_mass-%i",years[iy])  ,  0.0 );
        RooRealVar* b4_bkg_mass = new RooRealVar(Form("b4_bkg_mass-%i",years[iy]) , Form("b4_bkg_mass-%i",years[iy])  ,  0.1 , 0., pol_bmax);

        b0_bkg_mass->setVal(wsp_sb[iy]->var("pol_b0")->getVal());
        b1_bkg_mass->setVal(wsp_sb[iy]->var("pol_b0")->getVal());
        b2_bkg_mass->setVal(wsp_sb[iy]->var("pol_b0")->getVal());
        b3_bkg_mass->setVal(wsp_sb[iy]->var("pol_b0")->getVal());
        b4_bkg_mass->setVal(wsp_sb[iy]->var("pol_b0")->getVal());
        if (q2Bin==4 && fitOption>0){
            b0_bkg_mass->setConstant(true);
            b3_bkg_mass->setConstant(true);
            //b1_bkg_mass-> setVal(wsp_sb[iy]->var(Form("#sigma_{RT2}^{%i}",q2Bin))->getVal() );
            bkg_mass_pdf = new RooBernstein(Form("bkg_mass1_%i",years[iy]),Form("bkg_mass1_%i",years[iy]),  *mass, RooArgList(*b0_bkg_mass, *b1_bkg_mass, *b2_bkg_mass, *b3_bkg_mass,* b4_bkg_mass));

        } else { 
            bkg_mass_pdf = new RooExponential(Form("bkg_mass1_%i",years[iy]),  Form("bkg_mass1_%i",years[iy]) ,  *mass,   *slope  );
        } 

        RooProdPdf* bkg_pdf = new RooProdPdf(Form(bkgpdfname.c_str(),years[iy]),
                        Form(bkgpdfname.c_str(),years[iy]),
                        RooArgList(*bkg_ang_pdf,*bkg_mass_pdf));


        //// Finally build full pdf, including or not Z component ////
        RooAddPdf* full_pdf = 0;
        //RooAddPdf* full_pdf_penalty = 0;

        // fraction of signal and bkg pdfs
        RooRealVar *fsig = new RooRealVar( ("fsig_"+shortString+"_"+year).c_str(), ("fsig_"+shortString+"_"+year).c_str(),0,1 );
        fsig->setVal(fsig_year[iy]);
        full_pdf = new RooAddPdf(("PDF_sig_ang_fullAngularMass_bkg_"+shortString+"_"+year).c_str(),
                            ("PDF_sig_ang_fullAngularMass_bkg_"+shortString+"_"+year).c_str(),
                                RooArgList(*pdf_sig_ang_mass_mfc, *bkg_pdf),
                                RooArgList(*fsig)
                            );
        int n_togen = nbkg_years[years[iy]][q2stat] + nsig_years[years[iy]][q2stat];
        cout << "togen " << n_togen << endl;
        RooAbsPdf::GenSpec* genSpec = full_pdf->prepareMultiGen( observables, NumEvents(gen_nevt.Poisson(n_togen)));
        for (uint itoy = 0; itoy <= lastSample-firstSample; itoy++){
            RooRandom::randomGenerator()->SetSeed(itoy+firstSample);
            RooDataSet *toy_mom = full_pdf->generate(*genSpec) ;
            toy_mom->SetName(Form("data_y%i_b%i_subs%i",years[iy],q2Bin,itoy));
            wsp_out->import(*toy_mom);
        }
    }
    fout->cd();
    wsp_out->Write();
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

  int q2stat = 0;
    if ( argc > 5 ) q2stat = atoi(argv[5]);

  int XGBv = 0; 
  if ( argc > 6 ) XGBv = atoi(argv[6]);

  bool localFiles = false;
  if ( argc > 7 && atoi(argv[7]) > 0 ) localFiles = true;

  float Fsvalue = 0.04;
  if ( argc > 8 ) Fsvalue = atof(argv[8]);


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
  cout << "Fs value     " << Fsvalue      << endl;
  //cout <<  "plot        " << plot         << endl;
  //cout <<  "save        " << save         << endl;
  cout <<  "years[0]    " << years[0]     << endl;
//   cout << years[1] << endl;
//   cout << years[2] << endl;


  /*if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;

  if ( q2Bin==-1 )   cout << "Running all the q2 bins" << endl;
  if ( parity==-1 )  cout << "Running both the parity datasets" << endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      cocktailsample(q2Bin, parity, multiSample, nSample, q2stat,XGBv, localFiles,years);
  else*/
    Generate_toy(q2Bin, parity, multiSample, nSample,q2stat, XGBv, localFiles, Fsvalue ,years);

  return 0;

}