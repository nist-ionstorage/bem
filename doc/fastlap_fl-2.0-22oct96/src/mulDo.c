/*
*
* This software is being provided to you, the LICENSEE, by the Massachusetts
* Institute of Technology (M.I.T.) under the following license.  By
* obtaining, using and/or copying this software, you agree that you have
* read, understood, and will comply with these terms and conditions:  
*
* Permission to use, copy, modify and distribute this software and its
* documentation for any purpose and without fee or royalty is hereby granted,
* provided that you agree to comply with the following copyright notice and
* statements, including the disclaimer, and that the same appear on ALL
* copies of the software and documentation, including modifications that you
* make for internal use or for distribution:
*
* Copyright 1992 by the Massachusetts Institute of Technology.  All rights
* reserved.  
*
* THIS SOFTWARE IS PROVIDED "AS IS", AND M.I.T. MAKES NO REPRESENTATIONS OR
* WARRANTIES, EXPRESS OR IMPLIED.  By way of example, but not limitation,
* M.I.T. MAKES NO REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS
* FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE OR
* DOCUMENTATION WILL NOT INFRINGE ANY THIRD PARTY PATENTS, COPYRIGHTS,
* TRADEMARKS OR OTHER RIGHTS.   
*
* The name of the Massachusetts Institute of Technology or M.I.T. may NOT
* be used in advertising or publicity pertaining to distribution of the
* software.  Title to copyright in this software and any associated
* documentation shall at all times remain with M.I.T., and USER agrees to
* preserve same.  
*
* Written by: K. Nabors, T. Korsmeyer, and J. White
*
*/

/* # ***** sort to /src/main
   # ***** */
#include <stdio.h>
#include <math.h>
#include "mulStruct.h"
#include "mulGlobal.h"

#if OPCNT == ON
static directops=0, upops=0, downops=0, evalops=0;
#endif

/* 
Compute the direct piece. 
*/
mulDirect(sys)
ssystem *sys;
{
int i, j, k, psize, qsize;
double pc, *p, *qn, **mat;
cube *nextc;
/* Assumes the potential vector has been zero'd!!!! */
  for(nextc=sys->directlist; nextc != NULL; nextc = nextc->dnext) {
    psize = nextc->numfieldpts;       /* Equals number of field pts. */
    /* c2.0 You have to be in the eval list to be in the direct list. */ 
    p = nextc->eval;
    /* c2.0 The loop on nearest neighbors now may include the self cube. */
    for(i=nextc->directnumvects - 1; i >= 0; i--) {
      mat = nextc->directmats[i];
      qn = nextc->directq[i];
      for(j = psize - 1; j >= 0; j--) {
	for(k = nextc->directnumeles[i] - 1; k >= 0; k--) {
	  p[j] += mat[j][k] * qn[k];
#if OPCNT == ON
	  directops++;
#endif
	}
      }
    }
  }
}

/*
mulPrecond is the overlapped preconditioner.  The raw guess vector x, in
the linear system ACx=p is turned into Cx.  This is a bit confusing because
C is an operator which turns fields into charges (the inverse of what A does,
of course) so x is put into p and Cx appears in q, namely Cx is the new set
of charges to work with.  
*/
mulPrecond(sys, size)
ssystem *sys;
int size;
{
  int i, j, k, qsize;
  double *p, *q, *pn, **mat;
  cube *nextc, *pnbr;

  /* Assumes the charge vector q has been zeroed!!!! */
  for(nextc=sys->precondlist; nextc != NULL; nextc = nextc->pnext) {
    qsize = nextc->upnumeles[0];       /* Equals number of singularities. */
    q = nextc->upvects[0];
    /* Through all nearest neighbors with fieldpoints (pnbrs). */
    for(i=nextc->numpnbrs - 1; i >= 0; i--) {
      mat = nextc->precondmats[i];
      pnbr = nextc->pnbrs[i];
      pn = pnbr->eval;
      for(j = qsize - 1; j >= 0; j--) {
	for(k = pnbr->numfieldpts - 1; k >= 0; k--) {
	  q[j] += mat[j][k] * pn[k];
	}
      }
    }
  }
}

/*
spmulPrecond is the *special* overlapped preconditioner.  The raw guess vector x, 
in the linear system ACx=p is turned into Cx.  This is a bit confusing because
C is an operator which turns fields into charges (the inverse of what A does,
of course) so x is put into p and Cx appears in q, namely Cx is the new set
of charges to work with.  
*/
spmulPrecond(sys, work, size)
ssystem *sys;
double *work;
int size;
{
  int i, j, k, qsize;
  double *p, *q, **mat;
  cube *nextc, *pnbr;

  /* Assumes the charge vector q has been zeroed!!!! */
  for(nextc=sys->precondlist; nextc != NULL; nextc = nextc->pnext) {
    qsize = nextc->upnumeles[0];       /* Equals number of singularities. */
    q = nextc->upvects[0];
    /* Through all nearest neighbors with singularities as that is how
       we selected the field points.  So we call them pnbrs. */
    for(i=nextc->numqnbrs - 1; i >= 0; i--) {
      mat = nextc->precondmats[i];
      pnbr = nextc->qnbrs[i];
      for(j = qsize - 1; j >= 0; j--) {
	for(k = pnbr->upnumeles[0] - 1; k >= 0; k--) {
	  q[j] += mat[j][k] * work[((pnbr->sngs[k])->myfpt)->index];
	}
      }
    }
  }
}

/* 
Loop through upward pass. 
*/
mulUp(sys)
ssystem *sys;
{
int i, j, k, l;
int msize;
double *multi, *rhs, **mat;
cube *nc;

  if(sys->depth < 2) return;	/* ret if upward pass not possible/worth it */

/* Through all the depths, starting from the bottom and not doing top. */
  for(i = sys->depth; i > 0; i--) {  
  /* Through all the cubes at depth. */
    for(nc=sys->multilist[i]; nc != NULL; nc = nc->mnext) {
      msize = nc->multisize;
      multi = nc->multi;
      for(j=0; j < msize; j++) multi[j] = 0;
    /* Through all the nonempty children of cube. */
      for(j=nc->upnumvects - 1; j >= 0; j--) {
	mat = nc->upmats[j];
	rhs = nc->upvects[j];
	for(k = nc->upnumeles[j] - 1; k >= 0; k--) {
	  for(l = msize - 1; l >= 0; l--) {
	    multi[l] += mat[l][k] * rhs[k];
#if OPCNT == ON
	    upops++;
#endif
	  }
	}
      }
    }
  }
}


/*
 * mulEval does the evaluation pass. It can be used after mulDown or alone. 
 */
void mulEval(sys)
ssystem *sys;
{
  int i, j, k, size;
  cube *nc;
  double *eval, **mat, *vec;

  if(sys->depth < 2) return;	/* ret if upward pass not possible/worth it */

  for(nc = sys->evallist; nc != NULL; nc = nc->enext) {
    size = nc->numfieldpts;     /* number of field points in cube */
    eval = nc->eval;		/* vector of evaluation pnt potentials */

    /* do the evaluations */
    for(i = nc->evalnumvects - 1; i >= 0; i--) {
      mat = nc->evalmats[i];
      vec = nc->evalvects[i];
      for(j = size - 1; j >= 0; j--) {
	for(k = nc->evalnumeles[i] - 1; k >= 0; k--) {
	  eval[j] += mat[j][k] * vec[k];
#if OPCNT == ON
	  evalops++;
#endif
	}
      }
    }
  }
}

/* 
 * mulDown loops through the downward pass. 
 */
mulDown(sys)
ssystem *sys;
{
  cube *nc;
  int depth, i, j, k, lsize;
  double **mat, *rhs, *local;

  if(sys->depth < 2) return;  /* ret if downward pass not possible/worth it */

  for(depth=2; depth <= sys->depth; depth++) {
    for(nc=sys->locallist[depth]; nc != NULL; nc = nc->lnext) {
      lsize = nc->localsize;
      local = nc->local;
      for(j=0; j < lsize; j++) local[j] = 0;
  /* Through all the locals for the cube. */
      for(i=nc->downnumvects - 1; i >= 0; i--) {
	mat = nc->downmats[i];
	rhs = nc->downvects[i];
	for(j = lsize - 1; j >= 0; j--) {
	  for(k = nc->downnumeles[i] - 1; k >= 0; k--) {
	    local[j] += mat[j][k] * rhs[k];
#if OPCNT == ON
	    downops++;
#endif
	  }
	}
      }
    }
  }
}


#if OPCNT == ON
printops()
{
  printf("Number of Direct Multi-Adds = %d\n", directops);
  printf("Number of Upward Pass Multi-Adds = %d\n", upops);
  printf("Number of Downward Pass Multi-Adds = %d\n", downops);
  printf("Number of Evaluation Pass Multi-Adds = %d\n", evalops);
  printf("Total Number of Multi-Adds = %d\n", directops+upops+downops+evalops);
}
#endif


