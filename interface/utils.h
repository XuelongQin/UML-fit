#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <RooConstVar.h>
#include <map>

using namespace std;
using namespace RooFit;

extern std::map<int, float> scale_to_data;

// redo jpsi and psi2s
std::map<int,std::vector<float>> frt_sigmas = {
  {2016, {0.021, 0.015, 0.016, 0.011, 0.009, 0.009, 0.008, 0.011}},
  {2017, {0.018, 0.013, 0.014, 0.010, 0.008, 0.008, 0.007, 0.011}},
  {2018, {0.013, 0.010, 0.011, 0.007, 0.006, 0.006, 0.006, 0.007}},
};

std::map<int,std::vector<float>> fM_sigmas = {
  {2016, {0.023, 0.015, 0.017, 0.013, 0.005, 0.010, 0.006, 0.013}},
  {2017, {0.018, 0.014, 0.015, 0.010, 0.004, 0.008, 0.005, 0.011}},
  {2018, {0.015, 0.010, 0.011, 0.008, 0.006, 0.006, 0.006, 0.008}},
};


std::map<int,std::vector<float>> nbkg_years = {
  {2016, {162, 535, 462,  810, 0.005, 1342, 0.006, 467}},
  {2017, {185, 496, 441,  711, 0.004, 1363, 0.005, 379}},
  {2018, {288, 842, 734, 1270, 0.002, 2954, 0.003, 779}},
};

std::map<int,std::vector<float>> nsig_years = {
  {2016, {205, 454, 391,  689, 0.005, 1174, 0.006,  704}},
  {2017, {307, 581, 495, 1013, 0.004, 1524, 0.005,  835}},
  {2018, {500, 981, 821, 1608, 0.002, 3018, 0.003, 1836}},
};



RooPlot* prepareFrame(RooPlot* frame){
    frame->GetYaxis()->SetTitleOffset(1.8);
    frame->SetMaximum(frame->GetMaximum()*1.15);
    frame->SetMinimum(0);
    return frame;
}


RooGaussian* constrainVar(RooRealVar* var, 
                          string inVarName,  
                          RooWorkspace *w, 
                          int year,
                          bool addToList,
                          RooArgSet &c_vars,
                          RooArgSet &c_pdfs
                          ){
    RooGaussian* gauss_constr = new RooGaussian(  Form("c_%s_%i", inVarName.c_str(), year) , 
                                                  Form("c_%s_%i", inVarName.c_str(), year) , 
                                                  *var,  
                                                  RooConst( w->var(inVarName.c_str())->getVal()  ), 
                                                  RooConst( w->var(inVarName.c_str())->getError())
                                                 ); 
//     gauss_constr ->Print();
    if (addToList){
      c_vars.add(*var);    
      c_pdfs.add(*gauss_constr);
    }
    return gauss_constr;                 
}

void constrainVar2(RooRealVar* var, 
                          string inVarName,  
                          RooWorkspace *w, 
                          int year,
                          bool addToList,
                          RooArgSet &c_vars,
                          RooArgSet &c_pdfs
                          ){
    RooGaussian* gauss_constr = new RooGaussian(  Form("c_%s_%i", inVarName.c_str(), year) , 
                                                  Form("c_%s_%i", inVarName.c_str(), year) , 
                                                  *var,  
                                                  RooConst( w->var(inVarName.c_str())->getVal()  ), 
                                                  RooConst( w->var(inVarName.c_str())->getError())
                                                 ); 
    if (addToList){
      c_vars.add(*var);    
      c_pdfs.add(*gauss_constr);
    }
}


bool retrieveWorkspace(string filename, std::vector<RooWorkspace*> &ws, std::string ws_name){

    TFile* f =  TFile::Open( filename.c_str() ) ;
    if ( !f || !f->IsOpen() ) {
      cout << "File not found: " << filename << endl;
      return false;
    }
    RooWorkspace* open_w = (RooWorkspace*)f->Get(ws_name.c_str());
    if ( !open_w || open_w->IsZombie() ) {
      cout<<"Workspace "<< ws_name <<  "not found in file: " << filename << endl;
      return false;
    }
    ws.push_back( open_w );
    f->Close();
    return true;
}


std::vector<RooDataSet*> createDataset(int nSample, uint firstSample, uint lastSample, RooWorkspace *ws, 
                                       int q2Bin, int parity, int year, //std::map<int,float> scale_to_data,
                                       RooArgSet reco_vars, RooArgSet vars, std::string shortString  ){

    RooDataSet* dataCT, *dataWT;
    std::vector<RooDataSet*> datasample;

    if (nSample>0){  
      for (uint is = firstSample; is <= lastSample; is++) {

	RooDataSet* isample = new RooDataSet(("data_"+shortString + Form("_subs%i_%i", is,year)).c_str(), 
					     ("data_"+shortString + Form("_subs%i_%i", is,year)).c_str(), 
					     RooArgSet(reco_vars));

        dataCT = (RooDataSet*)ws->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin))
          ->reduce( RooArgSet(vars), Form("rand > %f && rand < %f", is*scale_to_data[year], (is+1)*scale_to_data[year] )) ;
        dataWT = (RooDataSet*)ws->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin))
          ->reduce( RooArgSet(vars), Form("rand > %f && rand < %f", is*scale_to_data[year], (is+1)*scale_to_data[year] )) ;

        isample->append(*dataCT);
        isample->append(*dataWT);
        datasample.push_back (isample);
      }
    }
    else{
      RooDataSet* isample = new RooDataSet(("data_"+shortString + "_subs0").c_str(), 
					   ("data_"+shortString + "_subs0").c_str(), 
					   RooArgSet(vars));
      dataCT = (RooDataSet*)ws->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin)) ;
      dataWT = (RooDataSet*)ws->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin)) ;
      isample->append(*dataCT);
      isample->append(*dataWT);
      /* isample->append(*((RooDataSet*)dataCT->reduce(EventRange(0,0.1*dataCT->numEntries())))); */
      /* isample->append(*((RooDataSet*)dataWT->reduce(EventRange(0,0.1*dataWT->numEntries())))); */
      datasample.push_back (isample);
    }
    return datasample;
}

std::vector<RooDataSet*> createnullDataset(int nSample, uint firstSample, uint lastSample, RooWorkspace *ws, 
                                       int q2Bin, int parity, int year, //std::map<int,float> scale_to_data,
                                       RooArgSet reco_vars, RooArgSet vars, std::string shortString  ){

    RooDataSet* dataCT, *dataWT;
    std::vector<RooDataSet*> datasample;

    if (nSample>0){  
      for (uint is = firstSample; is <= lastSample; is++) {

	RooDataSet* isample = new RooDataSet(("data_"+shortString + Form("_subs%i_%i", is,year)).c_str(), 
					     ("data_"+shortString + Form("_subs%i_%i", is,year)).c_str(), 
					     RooArgSet(reco_vars));
        datasample.push_back (isample);
      }
    }
    else{
      RooDataSet* isample = new RooDataSet(("data_"+shortString + "_subs0").c_str(), 
					   ("data_"+shortString + "_subs0").c_str(), 
					   RooArgSet(vars));
      dataCT = (RooDataSet*)ws->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin)) ;
      dataWT = (RooDataSet*)ws->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin)) ;
      /* isample->append(*((RooDataSet*)dataCT->reduce(EventRange(0,0.1*dataCT->numEntries())))); */
      /* isample->append(*((RooDataSet*)dataWT->reduce(EventRange(0,0.1*dataWT->numEntries())))); */
      datasample.push_back (isample);
    }
    return datasample;
    }


std::vector<RooDataSet*> createDataset(int nSample, uint firstSample, uint lastSample, RooWorkspace *ws, 
                                       int q2Bin, int parity, int year, int scale, //std::map<int,float> scale_to_data,
                                       RooArgSet reco_vars, RooArgSet vars, std::string shortString  ){

    RooDataSet* dataCT, *dataWT;
    std::vector<RooDataSet*> datasample;

    if (nSample>0){  
      for (uint is = firstSample; is <= lastSample; is++) {

	RooDataSet* isample = new RooDataSet(("data_"+shortString + Form("_subs%i_%i", is,year)).c_str(), 
					     ("data_"+shortString + Form("_subs%i_%i", is,year)).c_str(), 
					     RooArgSet(reco_vars));

        dataCT = (RooDataSet*)ws->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin))
          ->reduce( RooArgSet(vars), Form("rand > %f && rand < %f", scale*is*scale_to_data[year], scale*(is+1)*scale_to_data[year] )) ;
        dataWT = (RooDataSet*)ws->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin))
          ->reduce( RooArgSet(vars), Form("rand > %f && rand < %f", scale*is*scale_to_data[year], scale*(is+1)*scale_to_data[year] )) ;

        isample->append(*dataCT);
        isample->append(*dataWT);
        datasample.push_back (isample);
      }
    }
    else{
      RooDataSet* isample = new RooDataSet(("data_"+shortString + "_subs0").c_str(), 
					   ("data_"+shortString + "_subs0").c_str(), 
					   RooArgSet(vars));
      dataCT = (RooDataSet*)ws->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin)) ;
      dataWT = (RooDataSet*)ws->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin)) ;
      isample->append(*dataCT);
      isample->append(*dataWT);
      /* isample->append(*((RooDataSet*)dataCT->reduce(EventRange(0,0.1*dataCT->numEntries())))); */
      /* isample->append(*((RooDataSet*)dataWT->reduce(EventRange(0,0.1*dataWT->numEntries())))); */
      datasample.push_back (isample);
    }
    return datasample;
}


std::vector<RooDataSet*> createDatasetInctMC(int nSample, uint firstSample, uint lastSample, RooWorkspace *ws, 
                                       int q2Bin, int parity, int year,int q2stat, //std::map<int,float> scale_to_data,
                                       RooArgSet reco_vars, RooArgSet vars, std::string shortString ){

    RooDataSet* dataCT, *dataWT;
    std::vector<RooDataSet*> datasample;

    if (nSample>0){  
      for (uint is = firstSample; is <= lastSample; is++) {

	RooDataSet* isample = new RooDataSet(("data_"+shortString + Form("_subs%i_%i", is,year)).c_str(), 
					     ("data_"+shortString + Form("_subs%i_%i", is,year)).c_str(), 
					     RooArgSet(reco_vars));
        int totentries = ws->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin)) -> numEntries() + ws->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin)) -> numEntries();
        double scale = (double)nsig_years[year][q2stat]/totentries;
        dataCT = (RooDataSet*)ws->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin))
          ->reduce( RooArgSet(vars), Form("rand > %f && rand < %f", scale*is, scale*(is+1) )) ;
        dataWT = (RooDataSet*)ws->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin))
          ->reduce( RooArgSet(vars), Form("rand > %f && rand < %f", scale*is, scale*(is+1))) ;
        isample->append(*dataCT);
        isample->append(*dataWT);
        datasample.push_back (isample);
      }
    }
    else{
      RooDataSet* isample = new RooDataSet(("data_"+shortString + "_subs0").c_str(), 
					   ("data_"+shortString + "_subs0").c_str(), 
					   RooArgSet(vars));
      dataCT = (RooDataSet*)ws->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin)) ;
      dataWT = (RooDataSet*)ws->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin)) ;
      isample->append(*dataCT);
      isample->append(*dataWT);
      /* isample->append(*((RooDataSet*)dataCT->reduce(EventRange(0,0.1*dataCT->numEntries())))); */
      /* isample->append(*((RooDataSet*)dataWT->reduce(EventRange(0,0.1*dataWT->numEntries())))); */
      datasample.push_back (isample);
    }
    return datasample;

}


std::vector<RooDataSet*> createDatasetInData(RooWorkspace *ws, int q2Bin, RooArgSet vars, std::string shortString ){

  RooDataSet* data;
  std::vector<RooDataSet*> datasample;

  RooDataSet* isample = new RooDataSet(("data_"+shortString + "_subs0").c_str(), 
				       ("data_"+shortString + "_subs0").c_str(), 
				       RooArgSet(vars));
  data = (RooDataSet*)ws->data(Form("data_b%i",q2Bin)) ;
  isample->append(*data);
  /* isample->append(*((RooDataSet*)data->reduce(EventRange(0,0.1*data->numEntries())))); */
  datasample.push_back (isample);

  return datasample;
}

std::vector<RooDataSet*> createDatasetInData(int nSample, uint firstSample, uint lastSample, RooWorkspace *ws,
					     int q2Bin, int year, RooArgSet vars, std::string shortString, int q2st )
{

  std::vector<RooDataSet*> datasample;

  RooDataSet* data = (RooDataSet*)ws->data(Form("data_b%i",q2Bin)) ;

  if (nSample>0){

    int nSigCR[3] = {627365, 787149, 1694952};
    int subSampStat = (int)(data->sumEntries()*nsig_years[year][q2st]/nSigCR[year-2016]);

    for (uint is = firstSample; is <= lastSample; is++) {

      RooDataSet* isample = new RooDataSet(("data_"+shortString + Form("_subs%i", is)).c_str(),
					   ("data_"+shortString + Form("_subs%i", is)).c_str(),
					   RooArgSet(vars));

      isample->append(*((RooDataSet*)data->reduce(EventRange(is*subSampStat+1,
							     (is+1)*subSampStat))));
      datasample.push_back (isample);
    }
  }
  else{
    RooDataSet* isample = new RooDataSet(("data_"+shortString + "_subs0").c_str(),
					 ("data_"+shortString + "_subs0").c_str(),
					 RooArgSet(vars));
    isample->append(*data);
    datasample.push_back (isample);

  }


  return datasample;
}

std::vector<RooDataSet*> createDatasetInToy(RooWorkspace *ws, int q2Bin, int year, RooArgSet vars, std::string shortString ){

  RooDataSet* data;
  std::vector<RooDataSet*> datasample;

  RooDataSet* isample = new RooDataSet(("data_"+shortString + Form("_subs0_%i",year)).c_str(), 
				       ("data_"+shortString + Form("_subs0_%i",year)).c_str(), 
				       RooArgSet(vars));
  data = (RooDataSet*)ws->data(Form("data_b%ip1_subs0_%i",q2Bin,year)) ;
  isample->append(*data);
  /* isample->append(*((RooDataSet*)data->reduce(EventRange(0,0.1*data->numEntries())))); */
  datasample.push_back (isample);

  return datasample;
}

std::vector<double> generatemF(int q2Bin, int year){
  string filename_mc_mass = "";
  if (q2Bin==4) filename_mc_mass = Form("/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%i_fM_Jpsi_newbdt.root",year);
  else if (q2Bin==6) filename_mc_mass = Form("/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%i_fM_Psi_newbdt.root",year);
  else filename_mc_mass = Form("/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%i_fM_newbdt.root",year);
  TFile *file = new TFile(filename_mc_mass.c_str());
  RooWorkspace* w = (RooWorkspace*)file->Get("w");
  double nrt_mc   =  w->var(Form("nRT_%i",q2Bin))->getVal(); 
  double nwt_mc   =  w->var(Form("nWT_%i",q2Bin))->getVal(); 

  std::vector<double> mFset;
  double mF;
  int iy = year-2016;
  double mFv = nwt_mc/(nwt_mc+nrt_mc);
  double mFsigmav = 0.008;
  double frac_sigma = mFsigmav/mFv;
  RooRealVar* mFrac = new RooRealVar(Form("f_{M}^{%i}",year),"mistag fraction",1, 0.5, 1.5);
  RooGaussian* c_fm = new RooGaussian(Form("c_fm^{%i}",year) , "c_fm" , *mFrac,
                                      RooConst(1.) ,
                                      RooConst(frac_sigma));
  RooDataSet *mistaggauss = c_fm->generate(*mFrac, 100); 
  for (int i =0;i<100;i++){
    mF= mistaggauss->get(i)->getRealValue(Form("f_{M}^{%i}",year));
    mFset.push_back(mF);
  }
  return mFset;                                      
}

