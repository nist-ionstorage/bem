#ifndef _fastlap_h_included_
#define _fastlap_h_included_

#include "mulStruct.h"
#include "mulGlobal.h"

void matSetup(ssystem *sys, int size, snglrty *snglist);
snglrty *loadSnglrty(int numSing, double *px, int *pshape, int *plhsType, int *prhsType);
fieldpt *loadFieldpt(int size, double *px, int *dtype, double *xnrm);
void matSetup(ssystem *sys, int size, snglrty *snglist);
int ONsolve(ssystem *sys, snglrty *snglist, fieldpt *fptlist, int size, int maxiter, double *tol);
int gmres(ssystem *sys, snglrty *snglist, fieldpt *fptlist, double *p, double *r,
          double *ap, double *z, double **bv, double **bh, int size, int maxiter, double *tol);
void apply(ssystem *sys, double *q, double *p, int size);

#endif