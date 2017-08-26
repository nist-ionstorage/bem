#ifndef _mulmats_h_included_
#define _mulmats_h_included_

#include "mulStruct.h"
#include "mulGlobal.h"

void mulMatDirect(ssystem *sys);
void olmulMatPrecond(ssystem *sys);
void mulMatUp(ssystem *sys);
void mulMatEval(ssystem *sys);
void mulMatDown(ssystem *sys);

#endif
