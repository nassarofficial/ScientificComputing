#include "Interpolation.h"
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;


double Interpolation::Spline(double x[], double y[], double a, int n){
	int interval, i;
	double h[n], alpha[n], l[n], u[n], z[n], c[n], b[n], d[n], interpolant,aj,bj,cj,dj;
	//setting length of the intervals 
	for (i=0; i<n-1; i++){	
		h[i] = x[i+1]-x[i];  	
	}
	// 
	for (i=0; i<n-1; i++){
		alpha[i] = ((3/h[i]) * (y[i+1] - y[i]) )- ( (3/h[i-1]) * (y[i] - y[i-1]) );
	}
	//Tridiagonal matrix formation
	//z[0], z[n] = 0
	l[0] = 1;
 		u[0] = 0;
	z[0] = 0;
	for (i = 1; i<n-1; i++){
		l[i] = (2 * (x[i+1] - x[i-1])) - (h[i-1] * u[i-1]);
		u[i] = h[i] / l[i]; 
		z[i] = (alpha[i] - (h[i-1] * z[i-1]) ) / l[i]; //second derivative of f(x[i])
	}
	l[n] = 1;
	z[n] = 0;
	c[n] = 0;
	for (i=n-1; i>=0; i--){
		c[i] = z[i] - (u[i] * c[i+1]) ;
		b[i] = ( (y[i+1] - y[i]) / (h[i]) ) -  ((h[i]/3) * (c[i+1] + (2*c[i])));
		d[i] = (c[i+1] - c[i]) / (3*h[i]);
	}
	sort(x,x+n);
	if (a < x[0] || a > x[n-1]){
		cout << "Not in range";
		return 0;
	}
	else{
		interval = findingInterval(x, a, n);
		//calculated coefficients of Sj(x) = aj + bj(x-xj) + cj(x-xj)^2 + dj(x-xj) ^3
		aj = y[interval];
		bj = b[interval];
		cj = c[interval];
		dj = d[interval];
		interpolant = aj + bj + pow(  cj * (a - x[interval]) , 2) + pow( dj * (a - x[interval]) , 3);
		//cout << "\naj = " << aj << "\nbj = " << bj << "\ncj = " << cj << "\ndj = " << dj << "\n" << "x[interval] = " << x[interval] << "\n	" << "a: " << a << "\n";
		return interpolant;
	}
}
double Interpolation::findingInterval(double x[], double a, int n) {
		for (int i = n-1; i > 0; i--){
			if (a < x[i] && a > x[i-1])	{
				return i-1;	
			}
		}
    return 0;
}
double Interpolation::NewtInt(double x[], double y[], double a, int n)
{
    double fdd[n];   //initializng an array of Fdds 
	double sum,mult;
	int i,j;
    //initializing FDDs with Ys 
    for (i=0; i<n; i++){
        fdd[i] = y[i];
    }
    //calculating FDDs for each iteration 
    for (j=0; j<n-1; j++){
        for (i=n-1; i>j;i--){	
            	fdd[i] = ( fdd[i]-fdd[i-1] )/ ( x[i]-x[i-j-1]) ;
        }
    }
    //evaluating fn(x) while calculating coefficients 
    for(i=n-1;i>=0;i--)
    {
        mult = 1;
        for(j=0;j<i;j++){
            mult *= (a-x[j]);
        }
        mult *= fdd[i];
        sum  += mult;
        cout << "Iteration " << (n-1)-i << ": " << sum << "\n";
    }
    return sum;
}


