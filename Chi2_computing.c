// Chi2_computing.c
// inputs: data statistics, node sizes of the tested dual
// output: Chi2 value

#include <math.h>
#include "mex.h"

/* mex api */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )   
{
    double *p_stcs, *p_ns;    // inputs' pointers
    mwSize n_stcs;            // inputs' dimensions
    double *p_out;            // output
    mwSize i, j, k;           // temporary variables
    double Ns, *Nxs, Nys;
    mwSize Txy, nsx;
    
    /* input relating */
    p_stcs = mxGetPr(prhs[0]);
    p_ns = mxGetPr(prhs[1]);
    
    nsx = (mwSize)p_ns[0];
    Txy = nsx*(mwSize)p_ns[1];
    Nxs = mxGetPr(mxCreateDoubleMatrix(nsx,1,mxREAL));
//     mexPrintf("%d %d",nsx,Txy);
    /* read useful dimensianlity */
    n_stcs = mxGetM(prhs[0]);
    
    /* applying output */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    p_out = mxGetPr(plhs[0]);
    
    /* computing outputs */
    p_out[0] = 0.0;    // initial values
    Ns = 0.0;
    Nys = 0.0;
    for (j=0;j<nsx;j++) {Nxs[j] = 0.0;}
    if (p_stcs[0]>0.5){
        p_out[0] += p_stcs[0]*log(p_stcs[0]);
        Nxs[0] += p_stcs[0];
        Nys += p_stcs[0];
    }
    j = 2;
    k = 2;
    for (i=1;i<n_stcs;i++){
        if (p_stcs[i]>0.5){
            p_out[0] += p_stcs[i]*log(p_stcs[i]);
            Nxs[j-1] += p_stcs[i];
            Nys += p_stcs[i];
        }
        if (j==nsx){
            if (Nys>0.5) {
                p_out[0] -= Nys*log(Nys);
                Ns += Nys;
                Nys = 0.0;
            }
            j = 1;
        } else {j++;}
        if (k==Txy) {
            if (Ns>0.5) {
                for (k=0;k<nsx;k++) {
                    if (Nxs[k]>0.5) {p_out[0] += Nxs[k]*log(Ns/Nxs[k]);}
                    Nxs[k] = 0.0;
                }
                Ns = 0.0;
            }
            k = 1;
        } else {k++;}
    }
    p_out[0] += p_out[0];
}


