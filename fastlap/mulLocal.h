#ifndef _mullocal_h_included_
#define _mullocal_h_included_

#include "mulStruct.h"
#include "mulGlobal.h"

void evalFacFra(double **array, int order);
void evalSqrtFac(double **arrayout, double **arrayin, int order);
void evalSinCos(double beta, int order);
double sinB(int sum);
double cosB(int sum);
double **mulMulti2Local(double x, double y, double z,
                        double xp, double yp, double zp, int order, double **mat);
double **mulLocal2Local(double x, double y, double z,
                        double xc, double yc, double zc, int order);
void restartLocal(void);
double **mulQ2Local(snglrty **sngs, int numsngs,
                    double x, double y, double z, int order, int didthis, double **mat);
double **mulLocal2P(double x, double y, double z, fieldpt **fpts, int numfpts, int order);

#endif