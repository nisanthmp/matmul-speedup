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
    for(i = 0; i < m; i += 4)
        for(k = 0; k < p; ++k)
            for(j = 0; j < n; j += 16)            
            {
                C(i,j) += A(i,k) * B(k,j);
                C(i,j+1) += A(i,k) * B(k,j+1);
                C(i,j+2) += A(i,k) * B(k,j+2);
                C(i,j+3) += A(i,k) * B(k,j+3);
                C(i,j+4) += A(i,k) * B(k,j+4);
                C(i,j+5) += A(i,k) * B(k,j+5);
                C(i,j+6) += A(i,k) * B(k,j+6);
                C(i,j+7) += A(i,k) * B(k,j+7);
                C(i,j+8) += A(i,k) * B(k,j+8);
                C(i,j+9) += A(i,k) * B(k,j+9);
                C(i,j+10) += A(i,k) * B(k,j+10);
                C(i,j+11) += A(i,k) * B(k,j+11);
                C(i,j+12) += A(i,k) * B(k,j+12);
                C(i,j+13) += A(i,k) * B(k,j+13);
                C(i,j+14) += A(i,k) * B(k,j+14);
                C(i,j+15) += A(i,k) * B(k,j+15);

                C(i+1,j) += A(i+1,k) * B(k,j);
                C(i+1,j+1) += A(i+1,k) * B(k,j+1);
                C(i+1,j+2) += A(i+1,k) * B(k,j+2);
                C(i+1,j+3) += A(i+1,k) * B(k,j+3);
                C(i+1,j+4) += A(i+1,k) * B(k,j+4);
                C(i+1,j+5) += A(i+1,k) * B(k,j+5);
                C(i+1,j+6) += A(i+1,k) * B(k,j+6);
                C(i+1,j+7) += A(i+1,k) * B(k,j+7);
                C(i+1,j+8) += A(i+1,k) * B(k,j+8);
                C(i+1,j+9) += A(i+1,k) * B(k,j+9);
                C(i+1,j+10) += A(i+1,k) * B(k,j+10);
                C(i+1,j+11) += A(i+1,k) * B(k,j+11);
                C(i+1,j+12) += A(i+1,k) * B(k,j+12);
                C(i+1,j+13) += A(i+1,k) * B(k,j+13);
                C(i+1,j+14) += A(i+1,k) * B(k,j+14);
                C(i+1,j+15) += A(i+1,k) * B(k,j+15);

                C(i+2,j) += A(i+2,k) * B(k,j);
                C(i+2,j+1) += A(i+2,k) * B(k,j+1);
                C(i+2,j+2) += A(i+2,k) * B(k,j+2);
                C(i+2,j+3) += A(i+2,k) * B(k,j+3);
                C(i+2,j+4) += A(i+2,k) * B(k,j+4);
                C(i+2,j+5) += A(i+2,k) * B(k,j+5);
                C(i+2,j+6) += A(i+2,k) * B(k,j+6);
                C(i+2,j+7) += A(i+2,k) * B(k,j+7);
                C(i+2,j+8) += A(i+2,k) * B(k,j+8);
                C(i+2,j+9) += A(i+2,k) * B(k,j+9);
                C(i+2,j+10) += A(i+2,k) * B(k,j+10);
                C(i+2,j+11) += A(i+2,k) * B(k,j+11);
                C(i+2,j+12) += A(i+2,k) * B(k,j+12);
                C(i+2,j+13) += A(i+2,k) * B(k,j+13);
                C(i+2,j+14) += A(i+2,k) * B(k,j+14);
                C(i+2,j+15) += A(i+2,k) * B(k,j+15);

                C(i+3,j) += A(i+3,k) * B(k,j);
                C(i+3,j+1) += A(i+3,k) * B(k,j+1);
                C(i+3,j+2) += A(i+3,k) * B(k,j+2);
                C(i+3,j+3) += A(i+3,k) * B(k,j+3);
                C(i+3,j+4) += A(i+3,k) * B(k,j+4);
                C(i+3,j+5) += A(i+3,k) * B(k,j+5);
                C(i+3,j+6) += A(i+3,k) * B(k,j+6);
                C(i+3,j+7) += A(i+3,k) * B(k,j+7);
                C(i+3,j+8) += A(i+3,k) * B(k,j+8);
                C(i+3,j+9) += A(i+3,k) * B(k,j+9);
                C(i+3,j+10) += A(i+3,k) * B(k,j+10);
                C(i+3,j+11) += A(i+3,k) * B(k,j+11);
                C(i+3,j+12) += A(i+3,k) * B(k,j+12);
                C(i+3,j+13) += A(i+3,k) * B(k,j+13);
                C(i+3,j+14) += A(i+3,k) * B(k,j+14);
                C(i+3,j+15) += A(i+3,k) * B(k,j+15);
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
