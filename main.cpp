//
//  main.cpp
//  SC Proj 1
//
//  Created by Ahmed Nassar on 5/1/16.
//  Copyright © 2016 Ahmed Nassar. All rights reserved.
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <vector>
#include <time.h>
#include <complex>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <algorithm>

using namespace std;
using namespace std::chrono;


vector<double> get_row(vector<vector< double > > C, int row){
    int i, n;
    vector<double> temp;
    n = int(C.size());
    
    for (i=0;i<n; i++) {
        temp.push_back(C[row][i]);
        // cout << C[row][i] << endl;
    }
    return temp;
}

vector<double> get_cols(vector<vector< double > > C, int col){
    int i, n;
    vector<double> temp;
    n = int(C.size());
    
    for (i=0;i<n; i++) {
        temp.push_back(C[i][col]);
    }
    return temp;
}


vector<double> linear_regression(vector<vector<double>> points){
    double sumx=0,sumxy=0, sumy=0,sumx2=0;
    int n;
    vector<double> result;
    n = int(points.size());
    
    for (int i =0 ; i<n ; i++ ){
            sumx=sumx+points[i][0];
            sumy=sumy+points[i][1];
        
            sumxy=sumxy+points[i][0]*points[i][1];
            sumx2=sumx2+points[i][0]*points[i][0];
    }
    
    double xm=sumx/n;
    double ym=sumy/n;
    
    double a1 = (n*sumxy-sumx*sumy)/(n*sumx2-sumx*sumx); //slope
    double a0 = ym-a1*xm;	//bias

    result.push_back(a0);
    result.push_back(a1);
    
    return result;
}

vector<vector<double> > polynomial_regression (vector<vector<double> > points,int n){
    int i,j,N;
    N = int(points.size());

    vector<double> X,Y;
    X.reserve(2*n+1);
    Y.reserve(n+1);
    
    vector<vector<double>> B(n+1,vector<double>(n+2));

    
    for (i=0;i<2*n+1;i++)
    {
        X[i]=0;
        for (j=0;j<N;j++)
            X[i]=X[i]+pow(points[j][0],i);
    }

    for (i=0;i<=n;i++)
        for (j=0;j<=n;j++)
            B[i][j]=X[i+j];
    //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)

    for (i=0;i<n+1;i++)
    {
        Y[i]=0;
        for (j=0;j<N;j++)
            Y[i]=Y[i]+pow(points[j][0],i)*points[j][1];
    }
    for (i=0;i<=n;i++)
        B[i][n+1]=Y[i];
    n=n+1;
    
    cout<<"------- The Normal Matrix -------"<<endl;
    for (i=0;i<n;i++)
    {
        for (j=0;j<=n;j++)
            cout<<B[i][j]<<"  ";
        cout<<"\n";
        
    }

    return B;
}

vector<double> gauss_seidel_sor(vector<vector< double > > ls,vector< double > rs, int maxit, double es, double lambda){
    int n,i, iter,j,k,l;
    double sentinel =1,mul,maxea;
    n = int(ls.size());
    vector< double > x,xold,temp,ea;
    vector<vector< double > > C;
    x.reserve(n);
    xold.reserve(n);

    temp.reserve(n);
    ea.reserve(n);

    C = ls;
    for (k=0; k<n; k++) {
        C[k][k] = 0.0;
    }
    
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            // divide each val by diagonal
            C[i][j] = C[i][j]/ls[i][i];
            
        }
        // divide right side by diagonal
        rs[i]=rs[i]/ls[i][i];
    }
    
    iter = 0;

    while (sentinel == 1 and iter <= maxit) {
        for (int t=0; t<n; t++) {
            xold[t] = x[t];
        }
        
        for (l=0; l<n; l++) {
            temp = get_row(C,l);
            mul = inner_product(temp.begin(), temp.end(), x.begin(), 0.0);
            x[l] = rs[l] - mul;
            
            // multiply the right hand side by the relaxation

            x[l] = lambda*x[l]+(1.0-lambda)*xold[l];
            // Checking convergence by calculating approximate error
            if (x[l] != 0) {
                ea[l] = fabs((x[l]-xold[l])/x[l]) * 100;
            }
        }
        maxea = *max_element(ea.begin(), ea.end());
        if (maxea <= es){     // exit loop if
            sentinel = 0;
        }
        iter=iter+1;
    }
//    cout << iter << endl;

    return x;
    
}

vector<double> gauss_seidel(vector<vector< double > > ls,vector< double > rs,int maxit){
    
    int n = 0, i = 0, j = 0, sentinel = 1, iter =0;
    double es = 0.00001, maxea;
    n = int(ls.size());
    
    vector< double > temp,ea,xold,y;
    temp.reserve(n);
    ea.reserve(n);
    xold.reserve(n);
    y.reserve(n);
    
    iter = 0;
    while (sentinel==1 and iter <= maxit)
    {
        for (int t=0; t<n; t++) {
            xold[t] = y[t];
        }

        for (i = 0; i < n; i++)
        {
            y[i] = (rs[i] / ls[i][i]); // right hand side divided by diagonal
            for (j = 0; j < n; j++)
            {
                if (j == i)
                    continue;
                y[i] = y[i] - ((ls[i][j] / ls[i][i]) * temp[j]); // substituting the calculated values into the equations
                
                // Store for output
                temp[i] = y[i];

            }
            if (y[i] != 0) {
                ea[i] = fabs((y[i]-xold[i])/y[i]) * 100;
            }
        }
        maxea = *max_element(ea.begin(), ea.end());
        if (maxea <= es){     // exit loop if
            sentinel = 0;
        }

        iter=iter+1;
    }
//    cout << iter << endl;

    return temp;
}

vector<double> gauss_elimination(vector<vector< double > > a) {

    double cur,total,mat, maxval;
    int i,j,k,n,p,maxr;
    n = int(a.size());
    vector<double> temp;
    temp.reserve(n);
    
    
    for( i=0; i < (n); i++){
        cur = a[i][i];
        
        p = i;
        // find biggest val in column
        maxval = abs(a[i][i]);
        maxr = i;
        for (int k=i+1; k<n; k++) {
            if (abs(a[k][i]) > maxval) {
                maxval = abs(a[k][i]);
                maxr = k;
            }
        }
        
        // Swap maximum row tih column
        for (k=i; k<n+1;k++) {
            double temp = a[maxr][k];
            a[maxr][k] = a[i][k];
            a[i][k] = temp;
        }
        
        // Get the Triangular Matrix // Forward Substitution
        for(j = i+1; j < n; j++){
            mat = a[j][i] / a[i][i];
            for(k=0; k <= n; k++)
                a[j][k] -= mat * a[i][k];
        }
    
    }
    
    // Backgrond Substitution

    for(i = n-1; i >= 0; i--)
    {
        total = 0;
        for(j=i+1; j<n; j++)
            total += a[i][j] * temp[j];
        temp[i] = (a[i][n] - total) / a[i][i];
    }


    return temp;
}


int main() {
  
    int i,n,k;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Section 1
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    vector<vector<double>> polynomial_points {{0,1},{1,1.8},{2,1.3},{3,2.5},{4,6.3}};
    vector<vector<double>> poly_out_normal;
    
    vector<vector<double>> linear_points {{95,85},{85,95},{80,70},{70,65},{60,70}};
    vector<double> linear_out;
    k = int(polynomial_points.size());
    
    linear_out = linear_regression(linear_points);

    cout << "Bias: " << linear_out[0] <<endl;
    cout << "Slope: " << linear_out[1] <<endl;

    poly_out_normal = polynomial_regression(polynomial_points,2);
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Section 2
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Section 3
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    // Equations in a vector
    
    // 4 Equations
    // vector<vector< double > > a { {10, -1, 2, 0, 6},{-1, 11, -1, 3, 25},{2, -1, 10, -1, -11},{0, 3, -1, 8, 15} };
    
    // 3 Equations
    vector<vector< double > > a { { {3, -0.1, -0.2, 7.85},{0.1, 7, -0.3, -19.3},{0.3, -0.2, 10, 71.4} } };

    // 2 Equations
    //vector<vector< double > > a { {3, 2, 18},{-1, 2, 2}};
    
    //n = int(a.size());
    
    n = int(poly_out_normal.size());
    
    // Initialize vector for result
    vector <double> res_ele;
    vector <double> res_seidal, res_seidal2;
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    cout<<"\n";
    cout << "------- Gauss Elimination -------" << endl;

    // Get result as vector
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    //res_ele = gauss_elimination(a);
    
    res_ele = gauss_elimination(poly_out_normal);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << "Time:" << duration << " microseconds"<< endl;

    // Print Out Result
    for(i=0; i<n; i++)
        cout << "x" << i+1 << ": " << res_ele[i] << endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    cout<<"\n";
    cout << "------- Gauss Seidel : Iteration -------" << endl;
    
    
    res_ele = gauss_elimination(poly_out_normal);
    
    // 3 Equations
//    vector <vector<double>> ls = { {3, -0.1, -0.2},{0.1, 7, -0.3},{0.3, -0.2, 10} };
//    vector <double> rs = {7.85,-19.3,71.4};

    
    vector <vector<double>> ls {{5,10,30},{10,30,100,37.1},{30,100,354,130.3}};
    vector <double> rs = {12.9,37.1,130.3};
    
    // 4 Equations
//    vector<vector<double>> ls = { {10, -1, 2, 0},{-1, 11, -1, 3},{2, -1, 10, -1},{0, 3, -1, 8} };
//    vector<double> rs = {6, 25, -11, 15};

    
    //Polynomials
    
    high_resolution_clock::time_point t3 = high_resolution_clock::now();
    res_seidal = gauss_seidel(ls,rs,600);
    high_resolution_clock::time_point t4 = high_resolution_clock::now();
    auto duration2 = duration_cast<microseconds>( t4 - t3 ).count();
    cout << "Time:" << duration2 << " microseconds"<< endl;

    // Print Out Result
    for(i=0; i<n; i++)
        cout << "x" << i+1 << ": " << res_seidal[i] << endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    cout<<"\n";
    cout << "------- Gauss Seidel : SOR -------" << endl;

    high_resolution_clock::time_point t5 = high_resolution_clock::now();
    res_seidal2 = gauss_seidel_sor(ls,rs,600, 0.00001, 1.0);
    high_resolution_clock::time_point t6 = high_resolution_clock::now();
    auto duration3 = duration_cast<microseconds>( t6 - t5 ).count();
    cout << "Time:" << duration3 << " microseconds"<< endl;

    // Print Out Result
    for(i=0; i<n; i++)
        cout << "x" << i+1 << ": " << res_seidal2[i] << endl;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    return 0;
}