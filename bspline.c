#include"spline.h"
#include"levmar.h"
#include<assert.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<float.h>
#include<stdbool.h>

#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

#define SPLINE_DEGREE 3//2
#define T_MIN 0
#define T_MAX 30 //.1
#define EPSILON .0001
#define YEARS_TO_PRINT 18

//must be sorted
const double KNOTS[] =
{
	-3.0, -2.0, -1.0,
	0, .3, .6, 1, 1.5,
	2, 4, 8, 16, 28,
	40, 52, 64, 80, 100
};

/*const double KNOTS[] = {-2,-1,
	0,.25,.5,.75,1
};*/

const int N_KNOTS = sizeof(KNOTS) / sizeof(double);

double basisFunctionDegree0(double t, int knotIndex)
{
	assert(knotIndex < N_KNOTS - 1);
	double val = 0;
	if(t >= KNOTS[knotIndex] && t < KNOTS[knotIndex + 1])
		val = 1;
	return val;
}

double basisFunctionDegreeN(double t, int knotIndex, int degree)
{
	if(degree == 0)
		return basisFunctionDegree0(t, knotIndex);

	double frac1 = (t - KNOTS[knotIndex]) / (KNOTS[knotIndex + degree] - KNOTS[knotIndex]);
	double frac2 = (KNOTS[knotIndex + degree + 1] - t) / (KNOTS[knotIndex + degree + 1] - KNOTS[knotIndex + 1]);
	double term1 = frac1 * basisFunctionDegreeN(t, knotIndex, degree - 1);
	double term2 = frac2 * basisFunctionDegreeN(t, knotIndex + 1, degree - 1);
	return term1 + term2;

}

double basisFunctionDeriv(double t, int knotIndex, int degree)
{
	double frac1 = degree / (KNOTS[knotIndex + degree] - KNOTS[knotIndex]);
	double frac2 = degree / (KNOTS[knotIndex + degree + 1] - KNOTS[knotIndex + 1]);
	double term1 = frac1 * basisFunctionDegreeN(t, knotIndex, degree - 1);
	double term2 = frac2 * basisFunctionDegreeN(t, knotIndex + 1, degree - 1);
	return term1 - term2;
}

double basisFunctionIntegral(double t, int startKnotIndex, int degree)
{
	double sum = 0;
	double frac = (KNOTS[startKnotIndex + degree + 1] - KNOTS[startKnotIndex]) / (degree + 1);

	int knotIndex;
	for(knotIndex = startKnotIndex; knotIndex < N_KNOTS - degree - 2; knotIndex++)
		sum += frac * basisFunctionDegreeN(t, knotIndex, degree + 1);
	return sum;
}

double basisFunctionIntegralBounded(double ta, double tb, int knotIndex, int degree)
{
	double term1 = basisFunctionIntegral(tb, knotIndex, degree);
	double term2 = basisFunctionIntegral(ta, knotIndex, degree);
	return term1 - term2;
}

double sumBasisFunctionIntegrals(double *basisCoeffs, int numBasisCoeffs, double ta, double tb, int degree)
{
	double sum = 0;
	int knotIndex;
	assert(numBasisCoeffs == N_KNOTS - degree - 1);
	for(knotIndex = 0; knotIndex < numBasisCoeffs; knotIndex++)
	{
		sum += basisCoeffs[knotIndex] * basisFunctionIntegralBounded(ta, tb, knotIndex, degree);
	}
	return sum;
}

double instantaneousRateFromBasisFunction(double *basisCoeffs, int numBasisCoeffs, double t, int degree)
{
	double sum = 0;
	int knotIndex;
	assert(numBasisCoeffs == N_KNOTS - degree - 1);
	for(knotIndex = 0; knotIndex < numBasisCoeffs; knotIndex++)
	{
		sum += basisCoeffs[knotIndex] * basisFunctionDegreeN(t, knotIndex, degree);
	}
	return sum;
}

double discountFactorFromBasisFunction(double *basisCoeffs, int numBasisCoeffs, double ta, double tb, int degree)
{
	double rt = sumBasisFunctionIntegrals(basisCoeffs, numBasisCoeffs, ta, tb, degree);
	double discountFactor = exp(-rt);
	return discountFactor;
}

double forwardRateFromBasisFunction(double *basisCoeffs, int numBasisCoeffs, double ta, double tb, double dayCountFraction, int degree)
{
	double rt = sumBasisFunctionIntegrals(basisCoeffs, numBasisCoeffs, ta, tb, degree);
	double forwardRate = 1 / dayCountFraction * (exp(rt) - 1);
	return forwardRate;
}

double annuityFromBasisFunction(double *basisCoeffs, int numBasisCoeffs, double tValuation, double tStart, double tEnd, double dayCountFraction, int degree)
{
	double annuity = 0;
	double couponDate = tStart + dayCountFraction; //no coupon paid at tStart
	for( ; couponDate <= tEnd; couponDate += dayCountFraction)
		annuity += discountFactorFromBasisFunction(basisCoeffs, numBasisCoeffs, tValuation, couponDate, degree);
	annuity *= dayCountFraction;
	return annuity;
}

double floatingRateFromBasisFunction(double *forwardBasisCoeffs, double *discountBasisCoeffs, int numBasisCoeffs, double tValuation, double tStart, double tEnd, double dayCountFraction, int degree)
{
	double floatingRate = 0;
	double couponDate = tStart + dayCountFraction; //no coupon paid at tStart
	for( ; couponDate <= tEnd; couponDate += dayCountFraction)
	{
		double forwardRate = forwardRateFromBasisFunction(forwardBasisCoeffs, numBasisCoeffs, couponDate - dayCountFraction, couponDate, dayCountFraction, degree);
		double discountFactor =  discountFactorFromBasisFunction(discountBasisCoeffs, numBasisCoeffs, tValuation, couponDate, degree);
		floatingRate += forwardRate * discountFactor;
	}
	floatingRate *= dayCountFraction;
	return floatingRate;
}

double breakEvenSwapRate(double *forwardBasisCoeffs, double *discountBasisCoeffs, int numBasisCoeffs, double tValuation, double tStart, double tEnd, double fixedDayCountFraction, double floatingDayCountFraction, int degree)
{
	double fixed = annuityFromBasisFunction(discountBasisCoeffs, numBasisCoeffs, tValuation, tStart, tEnd, fixedDayCountFraction, degree);
	double floating = floatingRateFromBasisFunction(forwardBasisCoeffs, discountBasisCoeffs, numBasisCoeffs, tValuation, tStart, tEnd, floatingDayCountFraction, degree);
	return floating / fixed;
}

double breakEvenBasisSpread(double *forwardBasisCoeffs, double *discountBasisCoeffs, int numBasisCoeffs, double tValuation, double tStart, double tEnd, double dayCountFraction, int degree)
{
	double floating = floatingRateFromBasisFunction(forwardBasisCoeffs, discountBasisCoeffs, numBasisCoeffs, tValuation, tStart, tEnd, dayCountFraction, degree);
	double floatingAdjustment = floatingRateFromBasisFunction(discountBasisCoeffs, discountBasisCoeffs, numBasisCoeffs, tValuation, tStart, tEnd, dayCountFraction, degree); //use discount coeffs for forward rate extrapolation
	double nonStandardFixed = annuityFromBasisFunction(discountBasisCoeffs, numBasisCoeffs, tValuation, tStart, tEnd, dayCountFraction, degree); //not using fixed day count as in breakEvenSwapRate
	return (floating - floatingAdjustment) / nonStandardFixed;
}

void bSplineOptimization(double *p, double *x, int m, int n, void *data)
{
	double T[] =
	{
		0, 0, 6, 97, 188, 279, 370, 461, 552, 643,
		731, 1096, 1461, 1827, 2557, 3653, 4383, 5479, 7305, 9132, 10958,
		91, 183, 275, 366, 548, 731, 1096, 1461, 1827, 2557, 3653, 4383, 5479, 7305, 9132, 10958
	};

	int i = 0;
	m /= 2;

	x[0] = instantaneousRateFromBasisFunction(&p[0], m, i, 3);
	x[1] = instantaneousRateFromBasisFunction(&p[m], m, i, 3);

	for(i = 2; i < 10; i++ )
		x[i] = forwardRateFromBasisFunction(&p[0], m, T[i] / 360, (T[i] + 90) / 360, .25, 3);

	for(i = 10; i < 21; i++ )
		x[i] = breakEvenSwapRate(&p[0], &p[m], m, 0, 0, T[i] / 360, .5, .25, 3);

	for(i = 21; i < n; i++ )
		x[i] = breakEvenBasisSpread(&p[0], &p[m], m, 0, 0, T[i] / 360, .25, 3);
}

void testOptimization(double *p, double *x, int m, int n )
{
	printf("Testing optimization\n");
	
	double T[] =
	{
		0, 0, 6, 97, 188, 279, 370, 461, 552, 643,
		731, 1096, 1461, 1827, 2557, 3653, 4383, 5479, 7305, 9132, 10958,
		91, 183, 275, 366, 548, 731, 1096, 1461, 1827, 2557, 3653, 4383, 5479, 7305, 9132, 10958
	};

	
	int i = 0;
	m /= 2;
		
	assert( abs( x[0] - instantaneousRateFromBasisFunction(&p[0], m, i, 3) ) < EPSILON );
	assert( abs( x[1] - instantaneousRateFromBasisFunction(&p[m], m, i, 3) ) < EPSILON );

	for(i = 2; i < 10; i++ )
		assert( abs( x[i] - forwardRateFromBasisFunction(&p[0], m, T[i] / 360, (T[i] + 90) / 360, .25, 3) ) < EPSILON );

	for(i = 10; i < 21; i++ )
		assert( abs(x[i] - breakEvenSwapRate(&p[0], &p[m], m, 0, 0, T[i] / 360, .5, .25, 3) ) < EPSILON );

	for(i = 21; i < n; i++ )
		 assert( abs(x[i] - breakEvenBasisSpread(&p[0], &p[m], m, 0, 0, T[i] / 360, .25, 3) ) < EPSILON );

}

void testKnots()
{
	int tMinIndex = 0;
	int tMaxIndex = 0;

	int knotIndex;
	for(knotIndex = 0; knotIndex < N_KNOTS; knotIndex++)
	{
		if(KNOTS[knotIndex] == T_MIN)
			tMinIndex = knotIndex;

		if(KNOTS[knotIndex] > T_MAX)
		{
			tMaxIndex = knotIndex;
			break;
		}
	}

	assert(tMinIndex == SPLINE_DEGREE);
	assert(tMaxIndex == N_KNOTS - SPLINE_DEGREE - 2);
}

void testBasisFunction()
{
	int degree = SPLINE_DEGREE;
	double t;
	for(t = T_MIN; t <= T_MAX; t += .01)
	{
		double total = 0;
		for(int knotIndex = 0; knotIndex < N_KNOTS - degree - 2; knotIndex++)
		{
			double basis = basisFunctionDegreeN(t, knotIndex, degree);
			assert( basis >= 0 );
			total += basis;
		}

		assert(abs(total - 1.0) < EPSILON);
	}
}

void testBSpline()
{
	printf("Testing knots\n");
	testKnots();
	printf("Testing basis\n");
	testBasisFunction();
}

/*int main()
{
	double pOld[] = 
	{
		1.074308, -0.02751216, 0.008772677, 0.006251237, 0.007496152, 0.006584657, 0.008090704, 0.02769865,
 		0.03344964, 0.03718972, -0.006809754, 0.4342494, -373.7077, -2.3e-06, 1.139072, -0.03288946, 0.003648743,
		0.0007312583, 0.001670991, 0.001326829, 0.003167714, 0.02436871, 0.03270124, 0.03460087, -0.003641275, 
		0.3956395, -340.384, -2e-05 
	};


	double p[] =
	{
		0.03190845, 0.002814546, 0.006497451, 0.006858175, 0.007298168, 0.006660141, 0.008085883, 0.02793383, 0.03311066,
		0.03791496, -0.00926921, 0.4582085, -394.6686, -2.3e-06, -0.01310881, 0.0012291, 0.001026269, 0.001448263, 
		0.001419661, 0.001422399, 0.003151944, 0.02461342, 0.0323508, 0.03535055, -0.006182185, 0.4203752, -361.9657, -2e-05
	};


	double x[] =
	{
		0.0054588, 0.0007,
		0.005574, 0.006439, 0.006878, 0.007015,	0.007048, 0.006928, 0.007006, 0.007279,
		0.006900, 0.007980, 0.010080, 0.012480, 0.016900, 0.021060, 0.022980, 0.024780,
		0.025990, 0.026595, 0.026940,
		0.004588, 0.004925, 0.005113, 0.005250, 0.005350, 0.005300, 0.005163, 0.004925,
		0.004625, 0.004138, 0.003450,	0.003175, 0.002888, 0.002550, 0.002338, 0.002188
	};

	int m = sizeof(p) / sizeof(double);
	int n = sizeof(x) / sizeof(double);

	//optimization
	int ret;
	int i;
	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	opts[0] = LM_INIT_MU;
	opts[1] = 1E-15;
	opts[2] = 1E-15;
	opts[3] = 1E-20;
	opts[4] = LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing

	ret = dlevmar_dif(bSplineOptimization, p, x, m, n, 1000, opts, info, NULL, NULL, NULL); // no Jacobian

	printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", ret, info[5], info[6]);
	for(i = 0; i < m; ++i)
		printf("%.7g ", p[i]);
	printf("\n\nMinimization info:\n");
	for(i = 0; i < LM_INFO_SZ; ++i)
		printf("%g ", info[i]);
	printf("\n");
	
	
	//testing
	//testBSpline();	
	//testOptimization(p,x, m, n );

	int i;
	printf("=== Instantaneous 3M Libor Forward Curve ===\n" );
	for(i = 0; i < YEARS_TO_PRINT; i++)
	{
		double val= instantaneousRateFromBasisFunction(&p[0], m/2, i, 3);
		printf("%f\n", val);
	}
	
	printf("=== LIBOR/OIS Spread ===\n" );
	for(i = 0; i < YEARS_TO_PRINT; i++)
	{
		double val= instantaneousRateFromBasisFunction(&p[0], m/2, i, 3)-instantaneousRateFromBasisFunction(&p[14], m/2, i, 3);
		printf("%f\n", val);
	}
}*/


