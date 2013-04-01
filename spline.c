#include"spline.h"
#include"levmar.h"
#include<assert.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<float.h>
#include<stdbool.h> 

/*#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a > _b ? _a : _b; })*/
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
const double T[] = {
	-50.0, -30.0,-10.0, 0.0,
	20.0,
	50.0,
	120.0,
	250.0,
	500.0,
	2000.0,
	3500.0,
	5000.0,
	7500.0,
	11160.0,
	13680.0, 16200.0, 19800.0, 25200.0
};
const int K[] = { -3, -2, -1,
	0, 1, 2, 3,
	4, 5, 6, 7,
	8, 9, 10, 11,
	12, 13, 14
};

const int N = 10; 


int getEndIndex(double t){
	int index = 0;
	int maxIndex = 18 - 1; //TODO
	while( T[index] < t && index < maxIndex ){
		index++;
	}
	return index;
}

bool isWithinRange(int k, double t, int d){
	//assert( k >= T.size() );
	if( T[k] <= t && T[k+d+1] > t) return true;
	return false;
}

double SplineGenerator(int d, int k, double t){
	//assert( !(d < 0 || k < K[0] || k > K[N+6-d] ) ) ;
	if( d== 0 ){
		if( isWithinRange(k, t, d) ){
			return 1.0;
		}
		return 0.0;
	}

	double bPrev1 = SplineGenerator(d-1, k, t);
	double bPrev2 = SplineGenerator( d-1, k+1, t);
	double frac1 = (t-T[k])/(T[k+d]-T[k]);
	double frac2 = (T[k+d+1] - t)/(T[k+d+1]-T[k+1]);
	return frac1*bPrev1 + frac2*bPrev2;
}
// d - order of spline
// n - order of derivative
double SplineDeriv(int d, int k, double t, int n){
	//assert( !(d < 0 || k < K[0] || k > K[N+6-d] || d > 3) && n > 0 ) ;

	double frac1 = d/(T[k+d]-T[k]);
	double frac2 = d/(T[k+d+1]-T[k+1]);

	if( n == 1 )
	{
		if( !isWithinRange(k, t, d) || d== 0 ){
			return 0.0;
		}

		double bPrev1 = SplineGenerator(d-1, k, t);
		double bPrev2 = SplineGenerator( d-1, k+1, t);

		return frac1*bPrev1 - frac2*bPrev2;
	}


	double bPrev1 = SplineDeriv(d-1, k, t, n-1);
	double bPrev2 = SplineDeriv( d-1, k+1,t, n-1);
	return frac1*bPrev1 - frac2*bPrev2;
}




double SplineIntegrate(int d, int k, double upper){
	//assert( !(d < 0 || k < K[0] || k > K[N+6-d] || d > 2) ) ;
	double total = 0.0;
	if( upper > T[k] ){
		int endIndex = getEndIndex(upper);
		for( int i = k; i <= endIndex; i ++) {
			double frac = (T[k+d+1] - T[k])/(d+1); 
			double b = SplineGenerator(d+1, i, upper);
			total += frac*b;
		}
	}
	return total;
}

double SplineIntegrateBounded(int d, int k, double a, double b){
	return SplineIntegrate(d, k, b) - SplineIntegrate( d, k,a);
}

double SplineSecondOrderIntegral( int k, int l, double a, double b){
	double termOne = 0;

	if( k != l ){

		double bkOrderOne = SplineDeriv(3, k, b, 1); 
		double blOrderTwo = SplineDeriv(3, l, b, 2);
		double akOrderOne = SplineDeriv(3, k, a, 1); 
		double alOrderTwo = SplineDeriv(3, l, a, 2);

		termOne = bkOrderOne*blOrderTwo - akOrderOne*alOrderTwo;
	}

	int endIndex = getEndIndex(b) - 1;
	int startIndex = getEndIndex(a); 
	double termTwo = SplineDeriv(3, l, a, 3)*(SplineGenerator(3, k, T[startIndex]) - SplineGenerator(3, k, a) );
	termTwo += SplineDeriv(3, l, T[endIndex], 3)*(SplineGenerator(3, k, b) - SplineGenerator(3, k, T[endIndex]) );

	for(int i = startIndex; i <= endIndex; i++){
		double blOrderThree = SplineDeriv(3, l, T[i-1], 3);
		double bkDiff = SplineGenerator(3, k, T[i]) - SplineGenerator(3, k, T[i-1] );
		termTwo += blOrderThree*bkDiff;
	}

	return termOne - termTwo;
}

double Tikhonov(double *r, int m, double lower, double higher, double lambda ) 
{
	double total = 0;
	for(int i = 0; i < N+4; i++ ){
		int start = MAX(0,i-4); 
		for(int j = start; j < i; j++){
			total += 2*r[i]*r[j]*SplineSecondOrderIntegral(i, j, lower, higher); 
		}
		total += r[i]*r[i]*SplineSecondOrderIntegral(i, i, lower, higher);
	}
	total *= .5*lambda;
	return total;
}

double ForwardRate(double *l, int m, double S, double T, double offset )
{
	S -= offset;
	T -= offset;
	double total = 0;
	double delta = (T - S)/360; //T, S IN DAYS!!!
	for(int i = 0; i < m; i++){
		total += l[i]*SplineIntegrateBounded(3, i, S, T); 		
	}
	total = 1/delta*(exp(total)-1);
	return total; 
}


double TermRate(double *r, int m, double S, double T, double offset )
{
	S -= offset;
	T -= offset;
	double total = 0;
	for(int i = 0; i < m; i++){
		total += r[i]*SplineIntegrateBounded(3, i, S, T);
	}
	return total;
}

//both libor and ois
double InstaneousRate(double *r, int m, double t )
{
	double total = 0;
	for(unsigned long i = 0; i < m; i++){
		total += r[i]*SplineGenerator(3, i, t);
	}
	return total;
}

double DiscountFactor(double *f, int m, double S, double T, double offset )
{
	S -= offset;
	//T -= offset;
	double total = 0;
	for(int i = 0; i < m; i++){
		total -= f[i]*SplineIntegrateBounded(3, i, S, T); 		
	}
	return exp(total); 
}

double AnnuityValuation(double *f, int m, double S, double T, double offset, double dayCount )
{
	double alpha = dayCount/360; //day count fraction 180/360
	double total = 0;
	int nCoupons = (int)( (T-S) / dayCount );
	for( int i = 1; i <= nCoupons; i++ )
	{
		total += DiscountFactor( f, m, S, S + i*dayCount, offset );
	}
	total *= alpha;
	return total;
}

double SwapFloatingValuation(double *l, int m, double S, double T, double offset, double dayCount )
{
	double delta = dayCount/360;
	double total = 0;
	int nCoupons = ( int )( (T-S) / dayCount );

	for( int i = 1; i <= nCoupons; i++ )
	{
		double fwdRate = ForwardRate( l, m, S + (i-1)*dayCount, S + i*dayCount, offset );
		total += fwdRate* DiscountFactor( l, m, S + i*dayCount, T, offset );
	}
	total *= delta;
	return total;
}

double ImpliedSwapRate(double *f, int m, double *l, int n, double T, double dayCountAnnuity, double dayCountSwap )
{
	double swapFloatingValuation = SwapFloatingValuation( l, n, 0, T, 0, dayCountSwap );
	double annuityValuation = AnnuityValuation( f, m, 0, T, 0, dayCountAnnuity );
	return swapFloatingValuation / annuityValuation;
}

double ImpliedBasisSpread(double *f, int m, double *l, int n, double T, double dayCount )
{
	double num = SwapFloatingValuation(l, n, 0, T, 0, dayCount ) - SwapFloatingValuation(f, m, 0, T, 0, dayCount ) ;
	double den = AnnuityValuation( f, m, 0, T, 0, dayCount );
	assert( den != 0 );
	return num / den;
}


double SwapPV(double *f, int m, double S, double T, double offset, double dayCountFixed, double dayCountFloating, double couponRate )
{
	double fixed = AnnuityValuation( f, m, S, T, offset, dayCountFixed ) * couponRate;
	double floating = SwapFloatingValuation(f, m, S, T, offset, dayCountFloating );
	return fixed - floating;
}


/*int main()
  {
  printf("%f", SplineGenerator( 1,1, .1) );

  }*/

void area(double *p, double *x, int m, int n, void *data)
{
	register int i;

	for(i=0; i<n; i++)
		x[i]=-exp(p[0]-2)*(p[1]-3.5);
}

void spline(double *p, double *x, int m, int n, void *data)
{

	double d[35] = 
	{
		40898,
		40989,
		41080,
		41171,
		41262,
		41353,
		41444,
		41535,

		731,
		1096,
		1461,
		1827,
		2557,
		3653,
		4383,
		5479,
		7305,
		9132,
		10958,

		91,
		183,
		275,
		366,
		548,
		731,
		1096,
		1461,
		1827,
		2557,
		3653,
		4383,
		5479,
		7305,
		9132,
		10958

	};

	int i = 0;

	for(i = 0; i < 8; i++ )
	{
		x[i] = exp(TermRate(&p[0], 14, d[i], d[i]+90, 40892)-1)/4; 
	}
	for(i = 8; i < 19; i++ )
	{ 
		x[i] = ImpliedSwapRate( &p[14], 14, &p[0], 14, d[i], 180, 90 );
	}

	for(i = 19; i < n; i++ )
	{
		x[i] = ImpliedBasisSpread(&p[14], 14, &p[0], 14, d[i], 30 );
	}
}

/*int main()
{
	double p[28] =
	{ 0.00003672,
		0.00001374,
		0.00001978,
		0.00002284,
		0.00000829,
		0.00001560,
		0.00001770,
		0.00000978,
		0.00001010,
		0.00000639,
		0.00000433,
		-0.00000071,
		0.00000066,
		-0.00000230, 
		0.0000108642,
		-0.0000018685,
		-0.0000099828,
		-0.0000152516,
		-0.0000205267,
		-0.0000064600,
		0.0000070258,
		0.0000049045,
		0.0000054201,
		0.0000033175,
		0.0000006853,
		-0.0000135626,
		0.0000056523,
		-0.0000204759,
	};

	double x[35] = 
	{
		0.005574,
		0.006439,
		0.006878,
		0.007015,
		0.007048,
		0.006928,
		0.007006,
		0.007279,

		0.006900,
		0.007980,
		0.010080,
		0.012480,
		0.016900,
		0.021060,
		0.022980,
		0.024780,
		0.025990,
		0.026595,
		0.026940,

		0.004588,
		0.004925,
		0.005113,
		0.005250,
		0.005350,
		0.005300,
		0.005163,
		0.004925,
		0.004625,
		0.004138,
		0.003450,
		0.003175,
		0.002888,
		0.002550,
		0.002338,
		0.002188
	};


	double ans[28] =
	{
		0.01971084, 
		0.06077028,
		-0.02294467, 
		0.008100311, 
		-0.01779537, 
		-0.01561064, 
		-0.01227107, 
		0.001172863, 
		0.00458794, 
		0.002216043, 
		0.006467909, 
		0.003309302, 
		0.001374064, 
		-2.3e-06, 
		0.0007718124, 
		0.007924474, 
		0.02064914, 
		-0.0005977567, 
		-0.0156603, 
		-0.01834477, 
		0.02695873,
		0.01369224, 
		0.009672566, 
		0.003861522, 
		0.001193224, 
		0.0002345253,
		0.0001780809, 
		-2.04759e-05 
	}; 

	double v[10];
	int i;
	int startCount = 0;
	for(i = 0; i < 10; i++)
	{
		double val =(exp(TermRate(ans, 14,startCount, startCount + 90, 0)) - 1) * 4; 
		startCount += 100;
		printf("%f\n", val);
	}

*/
	/*int ret;
	  int i;
	  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	  opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	  opts[4]= LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 

	  int m=28; int n=35;
	  ret=dlevmar_dif(spline, p, x, m, n, 1000, opts, info, NULL, NULL, NULL);  // no Jacobian

	  printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", ret, info[5], info[6]);
	  for(i=0; i<m; ++i)
	  printf("%.7g ", p[i]);
	  printf("\n\nMinimization info:\n");
	  for(i=0; i<LM_INFO_SZ; ++i)
	  printf("%g ", info[i]);
	  printf("\n");
	 */






	/*
	   int n = 35; // num data points
	   int i = 0;

	   for(i = 0; i < 8; i++ )
	   {
	   double val = exp(TermRate(&p[0], 14, d[i], d[i]+90, 40892)-1)/4; 
	   printf("%f\n", val);
	   }
	   for(i = 8; i < 19; i++ )
	   { 
	   double val = ImpliedSwapRate( &p[14], 14, &p[0], 14, d[i], 180, 90 );
	   printf("%f\n", val);
	   }

	   for(i = 19; i < n; i++ )
	   {
	   double val = ImpliedBasisSpread(&p[14], 14, &p[0], 14, d[i], 30 );
	   printf("%f\n", val);
	   }  
	 */
	//return 0;
//}


