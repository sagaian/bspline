#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <math.h>
#include "spline.h"
#include "swaption.h"

double bCall(double T, double K, double Fo, double s)
{
	double total = 0;
	double dPlus  = (log(Fo/K)+.5*pow(s,2)*T)/(s*sqrt(T));
	double dMinus = (log(Fo/K)-.5*pow(s,2)*T)/(s*sqrt(T));
	total += Fo*gsl_cdf_ugaussian_P(dPlus);
	total -= K*gsl_cdf_ugaussian_P(dMinus);
	return total;
}

double bnCall(double T, double K, double Fo, double s)
{
	double total = 0;
	double d = (Fo-K)/(s*sqrt(T));
	total += d*gsl_cdf_ugaussian_P(d)+gsl_ran_ugaussian_pdf(d);
	total *= s*sqrt(T);
	return total;
}

double bPut(double T, double K, double Fo, double s)
{
	double total = 0;
	double dPlus  = (log(Fo/K)+.5*pow(s,2)*T)/(s*sqrt(T));
	double dMinus = (log(Fo/K)-.5*pow(s,2)*T)/(s*sqrt(T));
	total -= Fo*gsl_cdf_ugaussian_P(-dPlus);
	total += K*gsl_cdf_ugaussian_P(-dMinus);
	return total;
}

double bnPut(double T, double K, double Fo, double s)
{
	double total = 0;
	double d = -(Fo-K)/(s*sqrt(T));
	total += d*gsl_cdf_ugaussian_P(d)+gsl_ran_ugaussian_pdf(d);
	total *= s*sqrt(T);
	return total;
}

double pPay(double *l, double m, double To, double T, double K, double Fo, double s, model_type type)
{
	double p = AnnuityValuation(l, m, To, T, To,180);
	if(type == LOGNORMAL)
	{
		p*=bCall(T,K,Fo,s);
	}else{
		p*=bnCall(T,K,Fo,s);
	}	
	return p;
}

double pRec(double *l, double m, double To, double T, double K, double Fo, double s, model_type type)
{
	double p = AnnuityValuation(l, m, To, T, To,180);
	if(type == LOGNORMAL)
	{
		p*=bPut(T,K,Fo,s);
	}else{
		p*=bnPut(T,K,Fo,s);
	}
	return p;
}

int main()
{
	double l[] = {0.05784418, 0.06702633, -0.1040073, 0.02846515, -0.1164175, 0.000658484, 0.01335344, 
		     0.008509746, 0.003786817, 0.004015788, 0.000295547, 0.00124092, -0.000106187, -0.00000230};

	double f[] = {-0.004874, 0.008821, 0.011692, -0.007248, -0.068989, -0.007224, 0.017296, 0.007138,0.003684, 
		       0.003943, 0.000446, 0.000976, 0.000920, -0.000020};
	
	double a = pRec(&f[0], 14, 90, 450, .1, .09, .1, LOGNORMAL );
	double b = pRec(&f[0], 14, 90, 450, .1, .09, .1, NORMAL );

	double c = pPay(&f[0], 14, 90, 450, .1, .09, .1, LOGNORMAL );
	double d = pPay(&f[0], 14, 90, 450, .1, .09, .1, NORMAL );
	
	printf("%f\n", a);
}
