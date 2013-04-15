#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include "spline.h"
#include "swaption.h"
#include<assert.h>
#include"levmar.h"

#define ERROR .005

double bsDelta(double T, double K, double Fo, double s)
{
	double d1 = (log(Fo/K) + (pow(s,2) / 2) * T) / (s * sqrt(T));
	return gsl_cdf_ugaussian_P(d1);
}
 
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

double pPay(double *f, double m, double To, double T, double K, double Fo, double s, model_type type)
{
	double p = annuityFromBasisFunction(f, m, To, To, T, .5, 3);
	if(type == LOGNORMAL)
	{
		p*=bCall(T,K,Fo,s);
	}else{
		p*=bnCall(T,K,Fo,s);
	}	
	return p;
}

double pRec(double *f, double m, double To, double T, double K, double Fo, double s, model_type type)
{
	double p = annuityFromBasisFunction(f, m, To, To, T, .5, 3);
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

static inline double f(double T, double K, double Fo, double s, double sn)
{
	return bnCall(T,K,Fo,sn)-bCall(T,K,Fo,s);
}

static inline double f1(double T,double K, double Fo, double s)
{
	return bnCallDeriv(T,K,Fo,s);
}

double lognormalToNormalVol(double T, double K, double Fo, double s, double approx)
{
	double xPrev = approx;
	double x = 0;
	while(1)
	{
		x = xPrev - f(T,K,Fo,s,xPrev)/f1(T,K,Fo,xPrev);
		if(fabs(x - xPrev) < ERROR/10 || fabs(x) < ERROR/10) break;
		xPrev = x;
	};
	return x;
}

double sabrATMNormalVol(double maturity, double strike, double volatilityInitialValue, double alpha,  double rho, double beta)
{
	double alphaSq = alpha * alpha;
        double epsilon = maturity * alphaSq;
	double term1 = alpha / volatilityInitialValue * pow(strike, -beta);
	double term2 = (pow(beta, 2) - 2 * beta)/24 * pow(volatilityInitialValue, 2) / alphaSq * pow(strike, 2 * beta - 2);
	double term3 = beta * rho / (4 * alpha ) * volatilityInitialValue * pow(strike, beta - 1);
	double term4 = (2 - 3 * pow(rho,2) )/ 24;
	double atmNormalVol = term1 + term1 * epsilon * (term2 + term3 + term4);
	return atmNormalVol;
}


double sabrNormalVol(double maturity, double strike, double forwardInitialValue, double volatilityInitialValue, double alpha, double rho, double beta)
{
	double normalVol;
	if( strike == forwardInitialValue )
	{
		normalVol = sabrATMNormalVol(maturity, strike, volatilityInitialValue, alpha, rho, beta);	
	} else {
		double zeta = alpha*(pow(forwardInitialValue, 1 - beta) - pow(strike, 1 - beta))/volatilityInitialValue/(1 - beta);
		double delta = log((sqrt(1 - rho*zeta + zeta*zeta) + zeta - rho)/(1 - rho));
		double fMid = (forwardInitialValue + strike)/2;
		double gamma1 = beta/fMid;
		double gamma2 = beta*(beta - 1)/fMid/fMid;
		double epsilon = maturity*alpha*alpha;
		double c = pow(fMid, beta);
		normalVol = alpha*(forwardInitialValue - strike)/delta*(1 + ((2*gamma2 - gamma1*gamma1)*volatilityInitialValue*volatilityInitialValue*c*c/alpha/alpha/24 + rho*gamma1*volatilityInitialValue*c/alpha/4 + (2 - 3*rho*rho)/24)*epsilon);
		}
	return normalVol;
}

void getNormalVolData(double *lognormalVols, double *Fo, double * T, double *shifts, int numMaturities, int nShifts, double approx, double *normalVols)
{
	int i, j;
	printf("==== Normal Vols ====\n");
	for( i = 0; i < numMaturities; i++)
	{
		for( j = 0; j < nShifts; j++ )
		{
			int volIndex = i*nShifts + j;
			double lognormalVol = lognormalVols[volIndex];
			double normalVol = lognormalToNormalVol(T[i], Fo[i] + shifts[j], Fo[i],lognormalVol, approx);
			normalVols[volIndex] = normalVol;
			printf("%f, ", normalVol);
		}
		printf("\n");
	}
}

void sabrOptimization(double *p, double *x, int m, int n, void *data)
{
	double shifts[] = {-.01, -.005, 0, .005, .01 };
	double Fo[] = { 0.02177, 0.02245, 0.02382, 0.02656, 0.02907, 0.03205 };
	double T [] = {.25, .5, 1, 2, 3, 5};
	double numMaturities = sizeof(T)/sizeof(double);
 	int nShifts = sizeof(shifts)/sizeof(double);
	int i, j;	
	for( i = 0; i < numMaturities; i++)
	{
		for( j = 0; j < nShifts; j++ )
		{
			x[i * nShifts + j] = sabrNormalVol(T[i], Fo[i] + shifts[j], Fo[i], p[0], p[1], p[2], p[3]);
		}
	}
}

void getSABRCoeffs(double *coeffs, double *normalVols, double numCoeffs, double numVols, void *data)
{
	printf("==== SABR Coeffs (sigma, alpha, rho, beta) ====\n");
	int numIterations, i;
	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	opts[0] = LM_INIT_MU;
	opts[1] = 1E-15;
	opts[2] = 1E-15;
	opts[3] = 1E-20;
	opts[4] = LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing
	numIterations = dlevmar_dif(sabrOptimization, coeffs, normalVols, numCoeffs, numVols, 1000, opts, info, NULL, NULL, NULL); // no Jacobian
	assert( numIterations > 0);
	//printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", numIterations, info[5], info[6]);
	for(i = 0; i < numCoeffs + 1; ++i)
		printf("%.7g ", coeffs[i]); //include beta in output
	printf("\n");
}
 
void testOptimizationSabr(double *p, double *x, int m, int n )
{
	printf("Testing optimization\n");
	double shifts[] = {-.01, -.005, 0, .005, .01 };
	double Fo[] = { 0.02177, 0.02245, 0.02382, 0.02656, 0.02907, 0.03205 };
	double T [] = {.25, .5, 1, 2, 3, 5};
	double numMaturities = sizeof(T)/sizeof(double);
 	int nShifts = sizeof(shifts)/sizeof(double);
	int i, j;	
	for( i = 0; i < numMaturities; i++)
	{
		for( j = 0; j < nShifts; j++ )
		{
			double mkt =  *(x + i*nShifts + j);
			double model = sabrNormalVol(T[i], Fo[i] + shifts[j], Fo[i], p[0], p[1], p[2], p[3]);
			assert(fabs(mkt - model) < ERROR);
		}
	}
}

void testLognormalToNormalConversion(double *f, double *lognormalVols, double *normalVols, double numVols)
{
	printf("Testing lognormal to normal vol conversion\n");
	int i;
	for(i = 0; i < numVols; i++)
	{
		double a = pRec(&f[0], 14, 90, 450, .1, .09, lognormalVols[i], LOGNORMAL );
		double b = pRec(&f[0], 14, 90, 450, .1, .09, normalVols[i], NORMAL );
		assert( fabs(a - b) < ERROR );
		
		double c = pPay(&f[0], 14, 90, 450, .1, .09, lognormalVols[i], LOGNORMAL );
		double d = pPay(&f[0], 14, 90, 450, .1, .09, normalVols[i], NORMAL );
		assert( fabs(c - d) < ERROR );
	}	
}


int main()
{
	/*double l[] = 
	{
		0.03190845, 0.002814546, 0.006497451, 0.006858175, 0.007298168, 0.006660141, 0.008085883,
		0.02793383, 0.03311066, 0.03791496, -0.00926921, 0.4582085, -394.6686, -2.3e-06 
	};*/ 

	double f[] =
	{	
		-0.01310881, 0.0012291, 0.001026269, 0.001448263, 0.001419661, 0.001422399, 0.003151944, 
		0.02461342, 0.0323508, 0.03535055, -0.006182185, 0.4203752, -361.9657, -2e-05
	};

	double numBasisCoeffs = sizeof(f)/sizeof(double);
	
	double lognormalVols[] = {
		0.5731, 0.5172, 0.4767, 0.4572, 0.4492, 
		0.5509,	0.4993,	0.4629,	0.443,	0.4328,
		0.5152,	0.4687,	0.4373,	0.418,	0.4061,
		0.4608,	0.4208,	0.3936,	0.3758,	0.3642,
		0.4146,	0.3823,	0.3595,	0.3441,	0.3336,
		0.3733,	0.3455,	0.3255,	0.3115,	0.3021
	};
	
	double shifts[] = {-.01, -.005, 0, .005, .01 };
	double Fo[] = { 0.02177, 0.02245, 0.02382, 0.02656, 0.02907, 0.03205 };
	double T [] = {.25, .5, 1, 2, 3, 5};
	double approx = .1;

	int numMaturities = sizeof(T)/sizeof(double);
 	int numShifts = sizeof(shifts)/sizeof(double);
	int numVols = sizeof(lognormalVols)/sizeof(double);
	assert(numMaturities*sizeof(shifts)==sizeof(lognormalVols));
	assert(numShifts*numMaturities == numVols);	
	
	//normal vol conversion
	double *normalVols = malloc(sizeof(lognormalVols));
	getNormalVolData(lognormalVols, Fo, T, shifts, numMaturities, numShifts, approx, normalVols);
	//testLognormalToNormalConversion(f, lognormalVols, normalVols, numVols);
	
	//SABR optimization
	double sabrCoeffs[] = {.1, .3, 0, .5}; //sigma, alpha, rho, beta
	int numSABRCoeffs = sizeof(sabrCoeffs) / sizeof(double) - 1; //exclude beta from calibration
	getSABRCoeffs(sabrCoeffs, normalVols, numSABRCoeffs, numVols, NULL);
	//testOptimizationSabr(sabrCoeffs, normalVols, numSABRCoeffs, numVols );
	
	//shift in curve	
	double normalVol = sabrNormalVol(10, Fo[2], Fo[0], sabrCoeffs[0], sabrCoeffs[1], sabrCoeffs[2], sabrCoeffs[3]);		
	double recNoShift = pRec(f, numBasisCoeffs, 1, 10, Fo[2], Fo[2], normalVol, NORMAL);

	int i;
	for(i = 0; i < numBasisCoeffs; i++)
		f[i] -= .0001;


	double recShift = pRec(f, numBasisCoeffs, 1, 10, Fo[2], Fo[2], normalVol, NORMAL);
	double deltaSABR = (recNoShift - recShift)/.0001;
	double deltaBS = bsDelta(10, Fo[2], Fo[2], lognormalVols[12]);
	printf("==== Curve Shift ====\n");
	printf("Pre shift: %f, Post shift: %f, Delta from SABR: %f, Delta from BS: %f\n", recNoShift, recShift, deltaSABR, deltaBS );
	if(fabs(deltaSABR - deltaBS) > ERROR)
		printf("BS and SABR deltas are different. BS delta does not take into account volatility smile and will incorreclty hedge our delta risk.\n");
	free(normalVols);
}
