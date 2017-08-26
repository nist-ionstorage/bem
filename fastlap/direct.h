#ifndef _direct_h_included_
#define _direct_h_included_

#include "mulStruct.h"
#include "mulGlobal.h"

void coeffSwap(snglrty **sngs, int numsngs, int rows, double **mat);
void invert(double **mat, int size, int *reorder);

#endif