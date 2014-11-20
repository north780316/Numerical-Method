/*  */
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

/* Newton interpolation
 * PN(x) denotes polynomial of degree N
 * PN(x) = f(x0) + f[x1, x0]*(x - x0) + f[x2, x1, x0]*(x - x1)*(x - x0) + ...
 * PN(x) = f(x0) + summation from i = 1 to i = dataNum: f[xi, ..., x0]*(x - xi-1)*...*(x - x0)
 * */
#define MAX_DERIVATIVE_ORDER 3
#define DIFFERENTIAL 1e-1
#define DIFFERENTIAL_INVERSE 1e12

void computeDividedDifference ( int n, long double *x, long double *y, long double *d );
void computeDividedDifference ( int n, long double *x, long double *y, long double *d );
void printPolynomialFunction ( int n, long double *x, long double *y, long double *d );
void commandLineInteraction ( int n, long double *x, long double *y, long double *d );
void getInterpolationResult ( int n, long double *x, long double *y, long double *d, long double inputX, long double *intrResult );
void computeHighOrderDerivative ( int n, long double *x, long double *y, long double *d, long double inputX, long double *derivative, int order );
void computeHighOrderDerivativeUseDividedDifference ( int n, long double *x, long double *y, long double *d, long double inputX, long double *derivative, int order );

int
main ( int argc, char *argv[] ) 
{
	/* dataNum denotes the data quantity
	 * dataNum + 1 denotes the polynomial degree N 
	 * */
	int n = atoi(argv[1]);

	long double *x = (long double *) malloc (n * sizeof(long double));
	for ( int i = 0; i < n; ++i )
		x[i] = atof(argv[2 + i]);

	long double *y = (long double *) malloc (n * sizeof(long double));
	for ( int i = 0; i < n; ++i )
		y[i] = atof(argv[2 + i + n]);

	long double *d = (long double *) malloc (n * sizeof(long double));

	/* Compute the Newton divided difference */
	computeDividedDifference( n, x, y, d );

	/* Show the PN(x), the equation can be plotted by MATLAB */
	printPolynomialFunction( n, x, y, d );
	
	/* Interaction with user and get the interpolation value */
	commandLineInteraction( n, x, y, d );
	
	return 0;
}

void
computeDividedDifference ( int n, long double *x, long double *y, long double *d )
{
	/* Compute each of Newton divided difference. 
	 * d[0] is f(x0)
	 * d[1] is first divided difference: f[x1, x0] = (f(x1) - f(x0)) / (x1 - x0)
	 * d[n] is 'n' divided difference: f[xn, ..., x0] = (f[xn, ..., x1] - f[xn-1, ..., x0]) / (xn - x0)
	 * */
	/* Initialization */
	for ( int i = 0; i < n; ++i )
		d[i] = y[i];

	/* i+1 denotes the divided-difference length: [xi, xi-1, ..., x0] */
	for ( int i = 1; i < n; ++i )
		for ( int j = n - 1; j >= i; --j )
			d[j] = (d[j] - d[j - 1]) / (x[j] - x[j - i]);
}

void
printPolynomialFunction ( int n, long double *x, long double *y, long double *d )
{
	printf( "P%d(x) = %.5Le", n - 1, d[0] );
	for ( int i = 1; i < n; ++i )
	{
		if ( d[i] >= 0 )
			printf( " + %.5Le .* (", d[i] );
		else
			printf( " - %.5Le .* (", -d[i] );
		for ( int j = 0; j < i; ++j )
		{
			if ( x[j] >= 0 )
				printf( "x - %.5Le", x[j] );
			else
				printf( "x + %.5Le", x[j] );
			if ( j + 1 != i )
				printf( ") .* (" );
		}
		printf( ")" );
	}
	printf( "\n\n" );
	for ( int i = 1; i < n; ++i )
		printf( "f[x%d, ..., x0] = %.5Le\n", i, d[i] );
	printf( "\n" );
}

void
commandLineInteraction( int n, long double *x, long double *y, long double *d )
{
	printf( "Enter the Newton interpolation process\n" );

	long double max = x[0],
		   min = x[0];
	for ( int i = 1; i < n; ++i )
	{
		if ( x[i] > max ) max = x[i];
		if ( x[i] < min ) min = x[i];
	}

	/* Command line interaction with user */
	printf( "input x between x0 = %.5Le and x%d = %.5Le (0 to exit): ", min, n - 1, max );
	char buf[BUFSIZ];
	long double inputX, 
		   intrResult, derivative;
	while ( gets(buf) )
	{
		/* Exit */
		if ( strcmp(buf, "0") == 0 )
			break;
		inputX = atof(buf);

		/* Compute PN(x) by Horner's method */
		getInterpolationResult( n, x, y, d, inputX, &intrResult );
		printf( "interpolation result 'PN(%.5Le) = %.5Le'\n", inputX, intrResult );
		printf( "Use interpolation to evaluate derivatives:\n" );
		for ( int i = 1; i <= MAX_DERIVATIVE_ORDER && i < n; ++i )
		{
			computeHighOrderDerivative( n, x, y, d, inputX, &derivative, i );
			printf( "f^(%d)(%.5Le) = %.5Le'\n", i, inputX, derivative );
		}
		printf( "\nUse newton divided difference to evaluate derivatives:\n" );
		for ( int i = 1; i <= MAX_DERIVATIVE_ORDER && i < n; ++i )
		{
			computeHighOrderDerivativeUseDividedDifference( n, x, y, d, inputX, &derivative, i );
			printf( "f^(%d)(%.5Le) = %.5Le'\n", i, inputX, derivative );
		}

		printf( "\ninput x between x0 = %.5Le and x%d = %.5Le (0 to exit): ", min, n - 1, max );
	}
	
	printf( "\nExit the Newton interpolation process\n" );
}

void
getInterpolationResult ( int n, long double *x, long double *y, long double *d, long double inputX, long double *intrResult )
{
	/* Compute PN(x) by Horner's method */
	*intrResult = d[n - 1];
	for ( int i = n - 2; i >= 0; --i )
		*intrResult = *intrResult * (inputX - x[i]) + d[i];
}

void
computeHighOrderDerivative ( int n, long double *x, long double *y, long double *d, long double inputX, long double *derivative, int order )
{
	long double deltaX = DIFFERENTIAL,
		   coefficentTable[MAX_DERIVATIVE_ORDER][MAX_DERIVATIVE_ORDER + 1] = {
			{1, -1, 0, 0},
			{1, -2, 1, 0},
			{1, -3, 3, -1}
		   };

	long double dNf = 0, 
		   intrResult = 0;

	for ( int i = 0; i < order + 1; ++i )
	{
		getInterpolationResult( n, x, y, d, inputX - deltaX * i, &intrResult );
		dNf += coefficentTable[order - 1][i] * intrResult;
	}
	for ( int i = 0; i < order; ++i )
		dNf /= deltaX;

	*derivative = dNf;
}

void
computeHighOrderDerivativeUseDividedDifference ( int n, long double *x, long double *y, long double *d, long double inputX, long double *derivative, int order )
{
	long double factorial = 1;
	for ( int i = 1; i <= order; ++i )
		factorial *= i;
	*derivative = factorial * d[order];
}
