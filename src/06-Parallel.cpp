//
//  main.cpp
//  matmul
//
//  Created by Nisanth M P on 18/01/23.
//

#include <iostream>
//#include <bits/stdc++.h>
#include <sys/time.h>
#include <omp.h>

using namespace std;
#define A(i,j) a[ (i)*n + (j) ]
#define B(i,j) b[ (i)*n + (j) ]
#define C(i,j) c[ (i)*n + (j) ]

#define n 2048 
//4096
#define NUM_THREADS 4
double AS[n][n], BS[n][n], CS[n][n];
int pc = 128;
int rc = 128;

int main() {
    
    int p, q, r;
    p = q = r = n;

    // Initialize Matrices
    
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
        {
            AS[i][j] = (double)rand()/ (double)RAND_MAX;
        BS[i][j] = (double)rand()/ (double)RAND_MAX;
        CS[i][j] = 0;
        }

    struct timeval start, end;
 
    //omp_set_num_threads(NUM_THREADS); 
    // start timer.
    gettimeofday(&start, NULL);
  
    // unsync the I/O of C and C++.
    //ios_base::sync_with_stdio(false);

    // Matrix multiplication
// #pragma omp parallel for
    int i,j,k, pi, rk;
    double *a;
    double *b;
    double *c;
    /*register double Cij00, Cij01, Cij02, Cij03,
	            Cij04, Cij05, Cij06, Cij07,
	            Cij08, Cij09, Cij10, Cij11,
	            Cij12, Cij13, Cij14, Cij15;
    register double Aik;*/
    //omp_set_nested(1);
    #pragma omp parallel for
    for(rk = 0; rk < r; rk += rc)
    //#pragma omp parallel for
    for(pi = 0; pi < p; pi += pc)
    {
        a = &(AS[pi][rk]);
	b = &(BS[rk][0]);
	c = &(CS[pi][0]);
    //#pragma omp parallel for private(i,j,k) shared(A,B,C)
    //#pragma omp parallel for
    for(j = 0; j < q; j += 16)            
        for(i = 0; i < pc; i += 4)
	{
    	    //Cij00 = Cij01 = Cij02 = Cij03 = Cij04 = Cij05 = Cij06 = Cij07 = Cij08 = Cij09 = Cij10 = Cij11 = Cij12 = Cij13 = Cij14 = Cij15 = 0.0;
            for(k = 0; k < rc; ++k)
            {
		/*Aik = A[i][k];
                Cij00 += Aik * B[k][j];
                Cij01 += Aik * B[k][j+1];
                Cij02 += Aik * B[k][j+2];
                Cij03 += Aik * B[k][j+3];
                Cij04 += Aik * B[k][j+4];
                Cij05 += Aik * B[k][j+5];
                Cij06 += Aik * B[k][j+6];
                Cij07 += Aik * B[k][j+7];
                Cij08 += Aik * B[k][j+8];
                Cij09 += Aik * B[k][j+9];
                Cij10 += Aik * B[k][j+10];
                Cij11 += Aik * B[k][j+11];
                Cij12 += Aik * B[k][j+12];
                Cij13 += Aik * B[k][j+13];
                Cij14 += Aik * B[k][j+14];
                Cij15 += Aik * B[k][j+15];*/
                
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
            /*C[i][j] = Cij00;
            C[i][j+1] = Cij01;
            C[i][j+2] = Cij02;
            C[i][j+3] = Cij03;
            C[i][j+4] = Cij04;
            C[i][j+5] = Cij05;
            C[i][j+6] = Cij06;
            C[i][j+7] = Cij07;
            C[i][j+8] = Cij08;
            C[i][j+9] = Cij09;
            C[i][j+10] = Cij10;
            C[i][j+11] = Cij11;
            C[i][j+12] = Cij12;
            C[i][j+13] = Cij13;
            C[i][j+14] = Cij14;
            C[i][j+15] = Cij15;*/
	}
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
