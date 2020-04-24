/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

#include "Riostream.h" 

#include "Penalty.h" 

ClassImp(Penalty) 

Penalty::Penalty(const char *name, const char *title, 
		 RooAbsReal& _P1,
		 RooAbsReal& _P2,
		 RooAbsReal& _P3,
		 RooAbsReal& _P4p,
		 RooAbsReal& _P5p,
		 RooAbsReal& _P6p,
		 RooAbsReal& _P8p,
		 Double_t _power,
		 Double_t _coeff1,
		 Double_t _coeff4,
		 Double_t _coeff5,
		 bool _verbose) :
  RooAbsReal(name,title), 
  P1("P1","P1",this,_P1),
  P2("P2","P2",this,_P2),
  P3("P3","P3",this,_P3),
  P4p("P4p","P4p",this,_P4p),
  P5p("P5p","P5p",this,_P5p),
  P6p("P6p","P6p",this,_P6p),
  P8p("P8p","P8p",this,_P8p)
{

  power = _power;
  coeff1 = _coeff1;
  coeff4 = _coeff4;
  coeff5 = _coeff5;

  verbose = _verbose;

}


Penalty::Penalty(const Penalty& other, const char* name) :  
  RooAbsReal(other,name), 
  P1("P1",this,other.P1),
  P2("P2",this,other.P2),
  P3("P3",this,other.P3),
  P4p("P4p",this,other.P4p),
  P5p("P5p",this,other.P5p),
  P6p("P6p",this,other.P6p),
  P8p("P8p",this,other.P8p)
{

  power = other.power;
  coeff1 = other.coeff1;
  coeff4 = other.coeff4;
  coeff5 = other.coeff5;

  verbose = other.verbose;

}



Double_t Penalty::evaluate() const 
{

  double ctL4phi1 = P4p*P4p + P5p*P5p + P6p*P6p + P8p*P8p - 2 + 2*fabs( 2*P2 - P4p*P5p +P6p*P8p );

  double ret4 = pow(1.0/(1+exp(coeff4*ctL4phi1)),power);

  if (verbose && ctL4phi1>0) std::cout<<"[OUT] 4 "<<ctL4phi1<<std::endl;

  double a0 = 1 - P1*P1 - P6p*P6p*(1+P1) - P8p*P8p*(1-P1) - 4*P2*P2 - 4*P2*P6p*P8p; 
  double a4 = 1 - P1*P1 - P4p*P4p*(1+P1) - P5p*P5p*(1-P1) - 4*P2*P2 + 4*P2*P4p*P5p; 

  double a1 = 4*P3*P8p*P8p - 4*P3*P6p*P6p - 8*P1*P3 + 2*P5p*P6p*(1+P1) - 2*P4p*P8p*(1-P1) - 4*P2*P4p*P6p + 4*P2*P5p*P8p;
  double a3 = 4*P3*P4p*P4p - 4*P3*P5p*P5p + 8*P1*P3 + 2*P5p*P6p*(1-P1) - 2*P4p*P8p*(1+P1) - 4*P2*P4p*P6p + 4*P2*P5p*P8p;

  double a2 = 2 + 2*P1*P1 - 8*P2*P2 - 16*P3*P3 - (P4p*P4p+P6p*P6p)*(1-P1) - (P5p*P5p+P8p*P8p)*(1+P1) + 4*P2*P4p*P5p - 4*P2*P6p*P8p + 8*P3*P4p*P8p + 8*P3*P5p*P6p;

  double b0 = P8p*P8p - 1 - P1 + 2*P2 + P6p*P8p; 
  double b2 = P4p*P4p - 1 + P1 + 2*P2 - P4p*P5p; 
  double b1 = P4p*P8p - 2*P3 + 0.5 * ( P4p*P6p - P5p*P8p );
  
  double c0 = P8p*P8p - 1 - P1 - 2*P2 - P6p*P8p; 
  double c2 = P4p*P4p - 1 + P1 - 2*P2 + P4p*P5p; 
  double c1 = P4p*P8p - 2*P3 - 0.5 * ( P4p*P6p - P5p*P8p );
  
  int nSteps = 500;
  double phi, sin2, sincos, cos2;
  double ctL1, ctL5p, ctL5m, ctL5;
  double loc_penalty;
  double ret = 1.;
  double variable;

  for (int step = 0; step<nSteps; ++step) {
    phi = 3.14159 * step / nSteps;
    sin2 = sin(phi)*sin(phi);
    sincos = sin(phi)*cos(phi);
    cos2 = cos(phi)*cos(phi);

    ctL5p = b0*sin2 + b1*sincos + b2*cos2;

    ctL5m = c0*sin2 + c1*sincos + c2*cos2;

    ctL5 = TMath::Max(ctL5p,ctL5m);

    ctL1 = a0*sin2*sin2 + a1*sin2*sincos + a2*sin2*cos2 + a3*sincos*cos2 + a4*cos2*cos2;

    if (verbose && ctL5<0 && ctL1<0) std::cout<<"[OUT] 1 "<<ctL1<<" 5p "<<ctL5p<<" 5m "<<ctL5m<<" phi "<<phi<<std::endl;

    if ( ctL1 >= 0 ) {
      if ( ctL5 <= 0 ) variable = coeff1*ctL1;
      else variable = sqrt( coeff1*coeff1*ctL1*ctL1 + coeff5*coeff5*ctL5*ctL5 );
    } else {
      if ( ctL5 >= 0 ) variable = coeff5*ctL5;
      else variable = -1.0 / sqrt( 1.0/coeff1/coeff1/ctL1/ctL1 + 1.0/coeff5/coeff5/ctL5/ctL5 );
    }

    loc_penalty = pow(1.0/(1+exp(-1.0*variable)),power);
    if ( ret > loc_penalty ) ret = loc_penalty;

  }
  
  return ret * ret4;

}
