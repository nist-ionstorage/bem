#ifndef _mulmulti_h_included_
#define _mulmulti_h_included_

#include "mulStruct.h"
#include "mulGlobal.h"

int multerms(int order);
int costerms(int order);
int sinterms(int order);
void xyz2sphere(double x, double y, double z, double x0, double y0, double z0,
                double *rho, double *cosA, double *beta);
int iindex(int n, int m); //this function should be renamed
int sindex(int n, int m, int cterms);
double iPwr(int e);
double fact(int x);
void evalFactFac(double **array, int order);
void mulMultiAlloc(int maxsngs, int order, int depth);
void evalLegendre(double cosA, double *vector, int order);
double **mulQ2Multi(snglrty **sngs, int numsngs, double x, double y, double z,
                    int didthis, int order, double **mat);
void restartMulti(void);

#endif
