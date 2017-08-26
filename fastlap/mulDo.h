#ifndef _muldo_h_included_
#define _muldo_h_included_

#include "mulStruct.h"
#include "mulGlobal.h"

void mulDirect(ssystem *sys);
void mulPrecond(ssystem *sys, int size);
void spmulPrecond(ssystem *sys, double *work, int size);
void mulUp(ssystem *sys);
void mulEval(ssystem *sys);
void mulDown(ssystem *sys);

#endif