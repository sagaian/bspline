#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include "spline.h"
#include "swaption.h"
#include<assert.h>

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
	double p = annuityFromBasisFunction(l, m, To, To, T, .5, 3);
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
	double p = annuityFromBasisFunction(l, m, To, To, T, .5, 3);
	if(type == LOGNORMAL)
	{
		p*=bPut(T,K,Fo,s);
	}else{
		p*=bnPut(T,K,Fo,s);
	}
	return p;
}


double bnCallDeriv(double T, double K, double Fo, double s)
{
	double term1 = -(K-Fo)*gsl_sf_erfc(-(-K+Fo)/(sqrt(2*T)*s))/(2*s);
	double term2 = sqrt(T)*(exp(-pow(-K+Fo,2)/(2*T*pow(s,2)))/sqrt(2*M_PI))-term1;	
	return term1+term2;
}


double f(double T, double K, double Fo, double s, double sn)
{
	return bnCall(T,K,Fo,sn)-bCall(T,K,Fo,s);
}

double f1(double T,double K, double Fo, double s)
{
	return bnCallDeriv(T,K,Fo,s);
}

double lognormalToNormalVol(double T, double K, double Fo, double s, double approx, double epsilon)
{
	double xPrev = approx;
	double x = 0;
	while(1)
	{
		x = xPrev-f(T,K,Fo,s,xPrev)/f1(T,K,Fo,xPrev);
		if(abs(x-xPrev)<epsilon || x < epsilon) break;
	};
	return x;
}

double sabrNormalVol(double maturity, double strike, double forwardInitialValue, double volatilityInitialValue, double alpha, double beta, double rho)
{
        double zeta = alpha*(pow(forwardInitialValue, 1 - beta) - pow(strike, 1 - beta))/volatilityInitialValue/(1 - beta);
        double delta = log((sqrt(1 - rho*zeta + zeta*zeta) + zeta - rho)/(1 - rho));
        double fMid = (forwardInitialValue + strike)/2;
        double gamma1 = beta/fMid;
        double gamma2 = beta*(beta - 1)/fMid/fMid;
        double epsilon = maturity*alpha*alpha;
        double c = pow(fMid, beta);
        double normalVol = alpha*(forwardInitialValue - strike)/delta*(1 + ((2*gamma2 - gamma1*gamma1)*volatilityInitialValue*volatilityInitialValue*c*c/alpha/alpha/24 +
                                rho*gamma1*volatilityInitialValue*c/alpha/4 + (2 - 3*rho*rho)/24)*epsilon);
        return normalVol;
}

int main()
{
	double l[] = 
	{
		0.03190845, 0.002814546, 0.006497451, 0.006858175, 0.007298168, 0.006660141, 0.008085883,
		0.02793383, 0.03311066, 0.03791496, -0.00926921, 0.4582085, -394.6686, -2.3e-06 
	}; 

	double f[] =
	{	
		-0.01310881, 0.0012291, 0.001026269, 0.001448263, 0.001419661, 0.001422399, 0.003151944, 
		0.02461342, 0.0323508, 0.03535055, -0.006182185, 0.4203752, -361.9657, -2e-05
	};

	assert( sizeof(l)/sizeof(double) == 14);
	assert( sizeof(f)/sizeof(double) == 14);

	double nVol = lognormalToNormalVol(450,.1,.09,.1, .2,.000001);
	
	double a = pRec(&f[0], 14, 90, 450, .1, .09, .1, LOGNORMAL );
	double b = pRec(&f[0], 14, 90, 450, .1, .09, nVol, NORMAL );

	double c = pPay(&f[0], 14, 90, 450, .1, .09, .1, LOGNORMAL );
	double d = pPay(&f[0], 14, 90, 450, .1, .09, nVol, NORMAL );
	
	double diff = bCall(450, .1, .09, .1) - bnCall(450, .1, .09,nVol);
	double sabrVol = sabrNormalVol(450,.1,.09,nVol,.1,.1,.1); 
	printf("%f\n", diff);
	printf("%f,%f\n", a,b);
	printf("%f,%f\n", c,d);
	printf("%f\n", sabrVol);
}
