#ifndef _calp_h_included_
#define _calp_h_included_

#include "mulStruct.h"
#include "mulGlobal.h"

void initcalcp(snglrty *sing_list, int order, double* pAreas);
void centroid(snglrty *pp, double x2);
void Cross_Product(double vector1[], double vector2[], double result_vector[]);
void ComputeMoments(snglrty *panel, int order);
void calcSphericalPoint(snglrty *pp, int order);
void calcSphericalPanel(snglrty *pp, double **I, int order);
void zrotate(double cosa, double sina, double **M, double **Mt, int order);
void dp(snglrty *panel);
int planarize(snglrty *pq);
double jacobid(int m, int mp, int n, double x, double **binom, double *fact);
double **createBinom(int order);
double *createFactorial(int size);
double **createPmn(int order);
double calcp(snglrty *sing, double x, double y, double z, int deriv, double *normal, double *pfd);
void dumpnums(int flag, int size);
void testCalcp(snglrty *pp);
void restartCalcp(void);

#endif
