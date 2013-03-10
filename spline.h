#ifndef SPLINE_H
#define SPLINE_H

double //Generates spline using equations (25) and (26) in Lecture 1
	SplineGenerator(int d //dimension of spline
	, int k, double t);
double //Deriv of spline corresponding to equation (30) in Lecture 1
	SplineDeriv(int d, int k, double t, int n);
double //Integral of spline corresponding to equation (31) in Lecture 1
	SplineIntegrate(int d, int k, double upper);
double //Bounded integral of spline corresponding to equation (32) in Lecture 1 
	SplineIntegrateBounded(int d, int k, double a, double b);
double  //Bounded integral of second derivative of cubic splines corresponding to equation (33) in Lecture 1  
	SplineSecondOrderIntegral( int k, int l, double a, double b);
double // Tikhonov regularization calculation corresponding to equation (40)
	Tikhonov(double *r, int m, double lower, double higher, double lambda ); 
double //Fwd rate calculation corresponding to equation (39)
	ForwardRate(double *l, int m, double S, double T, double offset );
double //Term rate calculation  
	TermRate(double *r, int m, double S, double T, double offset );
double //Instaneous rate calculation corresponding to equations (35) and (36) 
	InstaneousRate(double *r, int m, double t );
double //Discount factor calculation corresponding to equation (37)
	DiscountFactor(double *f, int m, double S, double T, double offset );
double //Annuity calculation corresponding to equation (16)
	AnnuityValuation(double *l, int m, double S, double T, double offset, double dayCount );
double //Valuation formula for swap floating leg corresponding to equation (17) 
	SwapFloatingValuation(double *l, int m, double S, double T, double offset, double dayCount );
double //Break even swap rate calculation corresponding to equation (19) 
	ImpliedSwapRate(double *f, int m, double *l, int n,  double T, double dayCountAnnuity, double dayCountSwap );
double //Break even basis spread calculation corresponding to equation (22)
	ImpliedBasisSpread(double *f, int m,  double *l, int n, double T, double dayCount );
double //Swap PV by calculating difference between fixed and floating
	SwapPV(double *f, int m, double S, double T, double offset, double dayCountFixed, double dayCountFloating, double couponRate );
#endif
