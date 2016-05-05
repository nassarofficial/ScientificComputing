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
#include <chrono>
using namespace std;
using namespace std::chrono;

vector<float> gauss_seidel(vector<vector< float > > ls,vector< float > rs,int iter){

    int n = 0, i = 0, j = 0;
    
    n = int(ls.size());
    
    vector< float > temp;
    temp.reserve(n);
    
    vector< float > y ;
    y.reserve(n);
    
    while (iter > 0)
    {
        for (i = 0; i < n; i++)
        {
            y[i] = (rs[i] / ls[i][i]);
            for (j = 0; j < n; j++)
            {
                if (j == i)
                    continue;
                y[i] = y[i] - ((ls[i][j] / ls[i][i]) * temp[j]);
                
                // Store for output
                temp[i] = y[i];
            }
        }
        iter--;
    }

    return temp;
}

vector<float> gauss_elimination(vector<vector< float > > a) {

    float cur,total,mat;
    int i,j,k,n,p,z,y;
    n = int(a.size());
    vector<float> temp;
    temp.reserve(n);
    
    
    for( i=0; i < (n); i++){
        cur = a[i][i];
        
        p = i;
        
        // Find largest in the columns
        for(k = i+1; k < n; k++)
            if(fabs(cur) < fabs(a[k][i])){
                cur = a[k][i] ;
                p = k;
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
  
    int i, n;
    
    // Equations in a vector
    //vector<vector< float > > a { {3, -0.1, -0.2, 7.85},{0.1, 7, -0.3, -19.3},{0.3, -0.2, 10, 71.4} };
    
    // 6 Equations
    //vector<vector< float > > a { {0.6, 1, -2, 1, 3, -1, 4},{2, -1, 1, 2, 1, -3, 20},{1, 3, -3, 1, 2, 1, -15},{5, 2, -1, -1, 2, 1, -3},{-3, -1, 2, 3, 1, 3, 16},{4, 3, 1, -6, -3, -2, -27} };

    // 4 Equations
    //vector<vector< float > > a { {10, -1, 2, 0, 6},{-1, 11, -1, 3, 25},{2, -1, 10, -1, -11},{0, 3, -1, 8, 15} };
    
    // 3 Equations
    vector<vector< float > > a { { {3, -0.1, -0.2, 7.85},{0.1, 7, -0.3, -19.3},{0.3, -0.2, 10, 71.4} } };

    // 2 Equations
    //vector<vector< float > > a { {3, 2, 18},{-1, 2, 2}};
    
    n = int(a.size());
    
    // Initialize vector for result
    vector <float> res_ele;
    vector <float> res_seidal;
    
    
    cout << "------- Gauss Elimination -------" << endl;

    // Get result as vector
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    res_ele = gauss_elimination(a);
    
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << "Time:" << duration << " microseconds"<< endl;

    // Print Out Result
    for(i=0; i<n; i++)
        cout << "x" << i+1 << ": " << res_ele[i] << endl;

    
    cout << "------- Gauss Seidel -------" << endl;
    
    // 3 Equations
    vector <vector<float>> ls = { {3, -0.1, -0.2},{0.1, 7, -0.3},{0.3, -0.2, 10} };
    
    vector <float> rs = {7.85,-19.3,71.4};
    // 4 Equations
    //vector<vector<float>> ls = { {10, -1, 2, 0},{-1, 11, -1, 3},{2, -1, 10, -1},{0, 3, -1, 8} };
    //vector<float> rs = {6, 25, -11, 15};

    high_resolution_clock::time_point t3 = high_resolution_clock::now();
    res_seidal = gauss_seidel(ls,rs,7);
    high_resolution_clock::time_point t4 = high_resolution_clock::now();
    auto duration2 = duration_cast<microseconds>( t4 - t3 ).count();
    cout << "Time:" << duration2 << " microseconds"<< endl;

    // Print Out Result
    for(i=0; i<n; i++)
        cout << "x" << i+1 << ": " << res_seidal[i] << endl;

    

    return 0;
}