// yy_next_combination.c
// inputs:
// output:

#include <math.h>
#include "mex.h"

/* mex api */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )   
{
    double *prev_comb,*n;  // inputs' pointers
    double *next_comb;     // output
    int i, K, done=1;        // temporary variables
    
    /* input relating */
    prev_comb = mxGetPr(prhs[0]);
    K = mxGetN(prhs[0]);
    n = mxGetPr(prhs[1]);
    
    /* applying output */
    plhs[0] = mxCreateDoubleMatrix(1,K,mxREAL);
    next_comb = mxGetPr(plhs[0]);
    
    /* computing outputs */
    // first combnation
    if (prev_comb[0]<0.5){
        for (i=0;i<K;){next_comb[i]=++i;}
        return;
    }
    // the last element of prev_comb is smaller than n
    if (prev_comb[K-1]<n[0]){
        for (i=0;i<K;i++){
            next_comb[i] = prev_comb[i];
        }
        next_comb[K-1]++;
        return;
    }
    // other scenarios
    i = 0;
    while (done){
        if ((int)(n[0]-prev_comb[i]) > K-i-1){
            next_comb[i]=prev_comb[i];
            i++;
        } else {done = 0;} 
    }
    next_comb[i-1]++; 
    for (;i<K;i++){next_comb[i] = next_comb[i-1]+1;}
}


