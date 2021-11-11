#include "RooDoubleCBFast.h"
#include <RooRealVar.h>
#include <RooWorkspace.h>


RooDoubleCBFast* createWTMassShape(int q2Bin,
                                   RooRealVar* x,
                                   RooRealVar* mean_wt,
                                   RooRealVar* sigma_wt,
                                   RooRealVar* alpha_wt1,
                                   RooRealVar* alpha_wt2,
                                   RooRealVar* n_wt1,
                                   RooRealVar* n_wt2,
                                   RooWorkspace *w,
                                   int year,
                                   bool constrainVars, 
                                   RooArgSet &c_vars,
                                   RooArgSet &c_pdfs,
				   int contraintStat = -1
                                   ){

    RooDoubleCBFast* dcb_wt = new RooDoubleCBFast ( Form("dcb_wt_%i", year)  , 
                                                   "dcb_wt"      , 
                                                   *x, 
                                                   *mean_wt, *sigma_wt, *alpha_wt1, *n_wt1, *alpha_wt2, *n_wt2
                                                   );

    if (constrainVars){
      int scaleWidth = -1;
      if (q2Bin==4 && contraintStat>=0) scaleWidth = 2*contraintStat+1;

      constrainVar2(mean_wt  , Form("mean_{WT}^{%i}",q2Bin)   , w, year, true, c_vars, c_pdfs, scaleWidth);
      constrainVar2(sigma_wt , Form("#sigma_{WT1}^{%i}",q2Bin), w, year, true, c_vars, c_pdfs, scaleWidth);
      constrainVar2(alpha_wt1, Form("#alpha_{WT1}^{%i}",q2Bin), w, year, true, c_vars, c_pdfs, scaleWidth);
      constrainVar2(alpha_wt2, Form("#alpha_{WT2}^{%i}",q2Bin), w, year, true, c_vars, c_pdfs, scaleWidth);
      constrainVar2(n_wt1    , Form("n_{WT1}^{%i}",q2Bin)     , w, year, true, c_vars, c_pdfs, scaleWidth);
      constrainVar2(n_wt2    , Form("n_{WT2}^{%i}",q2Bin)     , w, year, true, c_vars, c_pdfs, scaleWidth);
    }

    return dcb_wt;                                                   
}




