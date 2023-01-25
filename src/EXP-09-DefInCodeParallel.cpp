//
//  01-DefInCode.cpp
//  matmul
//
//  Created by Nisanth M P on 18/01/23.
//

#include <iostream>
#include <sys/time.h>

using namespace std;
#define A(i,j) a[ (i)*n + (j) ]
#define B(i,j) b[ (i)*n + (j) ]
#define C(i,j) c[ (i)*n + (j) ]

#define N 2048 
#define NUM_THREADS 4
double AS[N][N], BS[N][N], CS[N][N];

int main() {
    
    int m, n, p;
    p = n = m = N;

    // Initialize Matrices
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
        {
            AS[i][j] = (double)rand()/ (double)RAND_MAX;
	    BS[i][j] = (double)rand()/ (double)RAND_MAX;
            CS[i][j] = 0;
        }

    struct timeval start, end;
 
    // start timer.
    gettimeofday(&start, NULL);
  
    // Matrix multiplication
    int i,j,k;
    double *a = &AS[0][0];
    double *b = &BS[0][0];
    double *c = &CS[0][0];
    #pragma omp parallel for
    for(i = 0; i < m; i += 1)
        for(j = 0; j < n; j += 1)            
            for(k = 0; k < p; ++k)
            {
		C(i,j) += A(i,k) * B(k,j);
	    }

    gettimeofday(&end, NULL);
  
    // Calculating total time taken by the program.
    double time_taken;
  
    time_taken = (end.tv_sec - start.tv_sec) * 1e6;
    time_taken = (time_taken + (end.tv_usec -
                              start.tv_usec)) * 1e-6;
  
    cout << "Time taken by program is : " << fixed
         << time_taken;
    cout << " sec" << endl;


    return 0;
}
