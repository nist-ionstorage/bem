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

/* # ***** sort to /src/direct
   # ***** */
#include <stdio.h>
#include "mulStruct.h"
#include "mulGlobal.h"

/* Allocates the direct matrix (for q to p) and fills it with calcp
 * coefficients for the left- and right-hand sides of the problem
 */
double **Q2P(sngs,numsngs,fpts,numfpts,swapOnly,mat)
snglrty **sngs;
fieldpt **fpts;
int numsngs, numfpts, swapOnly;
double **mat;
{
  int i, j;
  snglrty *pq;
  fieldpt *pf;
  double temp, calcp();
  if(swapOnly == TRUE)
      /* We've already allocated and loaded the matrix during the setup 
	 of the RHS, so we just have to recover the hidden LHS coefficients.*/
      coeffSwap(sngs, numsngs, numfpts, mat);
  else {
    CALLOC(mat, numfpts, double*, ON, AQ2P);
    for(i=0; i < numfpts; i++) {
      CALLOC(mat[i], 2 * numsngs, double, ON, AQ2P);
	/* We are doing a RHS or field computation so call calcp.  The first
	   half of the double-sized array will have the RHS coeff, the second
	   half the (hidden) LHS. */
      pf = fpts[i];
      for(j=0; j < numsngs; j++) {
	pq = sngs[j];
                          /* jt: new parameters for calcp() */
	mat[i][j] = calcp(pq, pf->x, pf->y, pf->z, pf->deriv, 
			  pf->nrm, &(mat[i][j+numsngs]));
      }
    }
  }
  return(mat);
}

/* Only allocates a matrix of identical dimension to half of
 * the direct matrix, as is needed for preconditioning.
 */
double **Q2PAlloc(nrows,ncols)
int nrows, ncols;
{
  double **mat;
  int i;
  CALLOC(mat, nrows, double*, ON, AQ2P);
  for(i=0; i < nrows; i++) {
    CALLOC(mat[i], ncols, double, ON, AQ2P);
  }
  return(mat);
}

/*
  - returned matrix has L below the diagonal, U above (GVL1 pg 58)
  - if allocate == TRUE ends up storing P and LU (could be a lot)
*/
double **ludecomp(matin, size, allocate)
double **matin;
int size;
int allocate;
{
  double factor, **mat;
  int i, j, k;

  if(allocate == TRUE) {
    /* allocate for LU matrix and copy A */
    MALLOC(mat, size, double*, ON, AMSC);
    for(i = 0; i < size; i++) {
      MALLOC(mat[i], size, double, ON, AMSC);
      for(j = 0; j < size; j++) mat[i][j] = matin[i][j];
    }
  }
  else mat = matin;

  for(k = 0; k < size-1; k++) {	/* loop on rows */
    if(mat[k][k] == 0.0) {
      fprintf(stderr, "FLE-ludecomp: zero piovt\n");
      exit(0);
    }
    for(i = k+1; i < size; i++) { /* loop on remaining rows */
      factor = (mat[i][k] /= mat[k][k]);
      for(j = k+1; j < size; j++) { /* loop on remaining columns */
	mat[i][j] -= (factor*mat[k][j]);
      }
    }
  }
  return(mat);
}

/*
  For direct solution of Pq = psi, used if preconditioning.
*/
void solve(mat, x, b, size)
double **mat, *x, *b;
int size;
{
  int i, j;

  /* copy rhs */
  if(x != b) for(i = 0; i < size; i++) x[i] = b[i];

  /* forward elimination */
  for(i = 0; i < size; i++) {	/* loop on pivot row */
    for(j = i+1; j < size; j++) { /* loop on elimnation row */
      x[j] -= mat[j][i]*x[i];
    }
  }

  /* back substitution */
  for(i--; i > -1; i--) {		/* loop on rows */
    for(j = i+1; j < size; j++) { /* loop on columns */
      x[i] -= mat[i][j]*x[j];
    }
    x[i] /= mat[i][i];
  }
}

/* 
  In-place inverts a matrix using guass-jordan.
*/
void invert(mat, size, reorder)
double **mat;
int size, *reorder;
{
  int i, j, k, best;
  double normal, multiplier, bestval, nextbest;
/*
  matlabDump(mat,size,"p");
*/
  for(i=0; i < size; i++) {

/* You should not ever need to pivot, but just in case. */
#if PIVOT == ON 
    best = i;
    bestval = ABS(mat[i][i]);
    for(j = i+1; j < size; j++) {
      nextbest = ABS(mat[i][j]);
      if(nextbest > bestval) {
	best = j;
	bestval = nextbest;
      }
    }

    * If reordering, find the best pivot. *
    if(reorder != NULL) {
      reorder[i] = best;
      if(best != i) {
	for(k=0; k < size; k++) {
	  bestval = mat[k][best];
	  mat[k][best] = mat[k][i];
	  mat[k][i] = bestval;
	}
      }
    }
#endif 

    /* First i^{th} column of A. */
    normal = 1.0 / mat[i][i];
    for(j=0; j < size; j++) {
      mat[j][i] *= normal;
    }
    mat[i][i] = normal;

    /* Fix the backward columns. */
    for(j=0; j < size; j++) {
      if(j != i) {
        multiplier = -mat[i][j];
        for(k=0; k < i; k++) {
          mat[k][j] += mat[k][i] * multiplier;
        }
        mat[i][j] = mat[i][i] * multiplier;
        for(k=i+1; k < size; k++) {
          mat[k][j] += mat[k][i] * multiplier;
        }
      }
    }
  }

  /* Unravel the reordering, starting with the last column. */
  if(reorder != NULL) {
    for(i=size-2; i >= 0; i--) {
      if(reorder[i] != i) {
	for(k=0; k < size; k++) {
	  bestval = mat[k][i];
	  mat[k][reorder[i]] = mat[k][i];
	  mat[k][i] = bestval;
	}
      }
    }
  }
/*
  matlabDump(mat,size,"c");
*/

}


/* 
Checks to see if the matrix has the M-matrix sign pattern and if
it is diagonally dominant. 
*/
matcheck(mat, rows, size)
double **mat;
int rows, size;
{
  double rowsum;
  int i, j;

  for(i = rows - 1; i >= 0; i--) {
    for(rowsum = 0.0, j = size - 1; j >= 0; j--) {
      if((i != j)  && (mat[i][j] > 0.0)) {
	printf("violation mat[%d][%d] =%g\n", i, j, mat[i][j]);
      }
      if(i != j) rowsum += ABS(mat[i][j]);
    }
    printf("row %d diag=%g rowsum=%g\n", i, mat[i][i], rowsum);
    if(rowsum > mat[i][i]) {
      for(j = size - 1; j >= 0; j--) {
	printf("col%d = %g ", j, mat[i][j]);
      }
      printf("\n");
    }
  }
}


matlabDump(mat, size, name)
double **mat;
int size;
char *name;
{
FILE *foo;
int i,j;
char fname[100];

  sprintf(fname, "%s.m", name);
  foo = fopen(fname, "w");
  fprintf(foo, "%s = [\n", name);
  for(i=0; i < size; i++) {
    for(j=0; j < size; j++) {
      fprintf(foo, "%.10e  ", mat[i][j]);
    }
    fprintf(foo, "\n");
  }
  fprintf(foo, "]\n");
}
/*
* This routine puts the hidden LHS coefficients in the part of the matrix 
* we actually use in the computations.
*/
coeffSwap(sngs, numsngs, rows, mat)
snglrty **sngs;
int numsngs, rows;
double **mat;
{
  int i, j;
  snglrty *pq;
  double temp;

  for(j=0; j < numsngs; j++) {
    pq = sngs[j];
    for(i=rows-1; i >= 0; i--) {
      temp = mat[i][j];
      mat[i][j] = mat[i][j+numsngs];
      mat[i][j+numsngs] = temp;
    }
  }
}








