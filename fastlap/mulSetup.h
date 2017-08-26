#ifndef _mulsetup_h_included_
#define _mulsetup_h_included_

#include "mulStruct.h"
#include "mulGlobal.h"

static int placeq(int flag, ssystem *sys, snglrty *snglrtys, int lhsSize, int rhsSize, fieldpt *fieldpts);
static void indexkid(ssystem *sys, cube *dad, int *qindex, int *pindex, int *pcindex);
static void getnbrs(ssystem *sys);
static void linkcubes(ssystem *sys);
static void setMaxq(ssystem *sys);
static void getAllInter(ssystem *sys);

void getrelations(ssystem *sys);
void setPosition(ssystem *sys);
void setExactUp(ssystem *sys, int numterms);
void setExactDown(ssystem *sys, int numterms);
int cntDwnwdChg(cube *cp, int depth);
void setTranslation(ssystem *sys, int val);

#endif
