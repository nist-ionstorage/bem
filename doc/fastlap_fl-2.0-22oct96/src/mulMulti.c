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

/*cdel extern int nr_mulLocal_calls, nr_mulMulti_calls;*/

/*
  Globals used for temporary storage.
*/
double *Irn, *Mphi;		/* (1/r)^n+1, m*phi vect's */
double *Ir, *phi;		/* 1/r and phi arrays, used to update above */
double *Rho, *Rhon;		/* rho and rho^n array */
double *Beta, *Betam;		/* beta and beta*m array */
double *tleg;		        /* Temporary Legendre storage. */
double **factFac;		/* factorial factor array: (n-m+1)...(n+m) */

/* 
   Used various places.  Returns number of coefficients in the multipole 
   expansion. 
*/
int multerms(order)
int order;
{
  return(costerms(order) + sinterms(order));
}

/*
  returns number of cos(m*phi)-weighted terms in the real (not cmpx) multi exp
*/
int costerms(order)
int order;
{
  return(((order + 1) * (order + 2)) / 2);
}

/*
  returns number of sin(m*phi)-weighted terms in the real (not cmpx) multi exp
*/
int sinterms(order)
int order;
{
  return((((order + 1) * (order + 2)) / 2) - (order+1));
}


/*
  takes two sets of cartesian absolute coordinates; finds rel. spherical coor.
*/
void xyz2sphere(x, y, z, x0, y0, z0, rho, cosA, beta)
double x, y, z, x0, y0, z0, *rho, *cosA, *beta;
{
  /* get relative coordinates */
  x -= x0;			/* "0" coordinates play the role of origin */
  y -= y0;
  z -= z0;
  /* get spherical coordinates */
  *rho = sqrt(x*x + y*y + z*z);

  if(*rho == 0.0) *cosA = 1.0;
  else *cosA = z/(*rho);

  if(x == 0.0 && y == 0.0) *beta = 0.0;
  else *beta = atan2(y, x);

}

/*
  gives the linear index into vector from n and m used by all routines dealing
  with moments (cosine parts) and Leg. function evals, e.g. Mn^m and Pn^m(cosA)
  used for all cases except for the sine part of arrays (use sindex() instead)
  assumed entry order: (n,m) = (0,0) (1,0) (1,1) (2,0) (2,1) (2,2) (3,0)...
*/
int index(n, m)
int n, m;
{
#if NOABRT == OFF
  if(m > n) {
    fprintf(stderr, "FLE-index: m = %d > n = %d\n", m, n);
    exit(0);
  }
  if(n < 0 || m < 0) {
    fprintf(stderr, "FLE-index: n = %d or m = %d negative\n", n, m);
    exit(0);
  }
#endif
  return(m + (n*(n+1))/2);
}

/*
  gives the linear index into vector from n and m used by all routines dealing
  with moments (sine parts), e.g. Mn^m 
  assumes an array with all m = 0 (Mn^0) entries omitted to save space
  assumed entry order: (n,m) = (1,1) (2,1) (2,2) (3,1) (3,2) (3,3) (4,1)...
*/
int sindex(n, m, cterms)
int n, m, cterms;		/* cterms is costerms(order) */
{
#if NOABRT == OFF
  if(m > n) {
    fprintf(stderr, "FLE-sindex: m = %d > n = %d\n", m, n);
    exit(0);
  }
  if(n < 0 || m < 0) {
    fprintf(stderr, "FLE-sindex: n = %d or m = %d negative\n", n, m);
    exit(0);
  }
  if(m == 0) {
    fprintf(stderr, "FLE-sindex: attempt to index M%d^0\n", n);
    exit(0);
  }
#endif
  return(cterms + m + (n*(n+1))/2 - (n+1));
}

/*
  returns i = sqrt(-1) to the power of the argument
*/
double iPwr(e)
int e;				/* exponent, computes i^e */
{
  if(e == 0) return(1.0);
  if(e % 2 != 0) {
    fprintf(stderr, "FLE-iPwr: odd exponent %d\n", e);
    exit(0);
  }
  else {
    e = e/2;			/* get power of negative 1 */
    if(e % 2 == 0) return(1.0);
    else return(-1.0);
  }
}

/*
  returns factorial of the argument (x!)
*/
double fact(x)
int x;
{
  double ret = 1.0;
  if(x == 0 || x == 1) return(1.0);
  else if(x < 0) {
    fprintf(stderr, "FLE-fact: attempt to take factorial of neg no. %d\n", x);
    exit(0);
  }
  else {
    while(x > 1) {
      ret *= x;
      x--;
    }
    return(ret);
  }
}

/*
  produces factorial factor array for mulMulti2P
*/
void evalFactFac(array, order)
int order;
double **array;
{
  int n, m;			/* array[n][m] = (m+n)!/(n-m)! */

  /* do first column of lower triangular part - always 1's */
  for(n = 0; n < order+1; n++) array[n][0] = 1.0;

  /* do remaining columns of lower triangular part */
  /*  use (n-m)!/(n+m)! = 1/(n-m+1)...(n+m) since m \leq n */
  /*  (array entry is divided into number to effect mul by (n-m)!/(n+m)!) */
  for(n = 1; n <= order; n++) {
    for(m = 1; m <= n; m++) {
      array[n][m] = (n-(m-1)) * array[n][m-1] * (n+m);
    }
  }

#if DISFAF == ON
  fprintf(stdout, "FACTORIAL FACTOR ARRAY:\n");
  dumpMat(array, order+1, order+1);
#endif

}

/*
  Allocates space for temporary vectors.
*/
void mulMultiAlloc(maxsngs, order, depth)
/*c2.0 maxsngs now means max # of singularities or fieldpoints in a cube.*/
int maxsngs, order, depth;
{
  int x;

#if DISSYN == ON
  extern int *multicnt, *localcnt, *evalcnt;
#endif
#if DMTCNT == ON
  extern int **Q2Mcnt, **Q2Lcnt, **Q2Pcnt, **L2Lcnt;
  extern int **M2Mcnt, **M2Lcnt, **M2Pcnt, **L2Pcnt;
#endif
  extern double *sinmkB, *cosmkB, **facFrA;


  if(order <= MAXORDER) order = MAXORDER;
  else {
#if NOWARN == OFF
    printf("FLW-mulMultiAlloc: larger than maximum order requested.\n");
#endif
  }


  CALLOC(Rho, maxsngs, double, ON, AMSC); /* rho array */
  CALLOC(Rhon, maxsngs, double, ON, AMSC); /* rho^n array */
  CALLOC(Beta, maxsngs, double, ON, AMSC); /* beta array */
  CALLOC(Betam, maxsngs, double, ON, AMSC); /* beta*m array */
  CALLOC(Irn, maxsngs, double, ON, AMSC);	/* (1/r)^n+1 vector */
  CALLOC(Ir, maxsngs, double, ON, AMSC); /* 1/r vector */
  CALLOC(Mphi, maxsngs, double, ON, AMSC); /* m*phi vector */
  CALLOC(phi, maxsngs, double, ON, AMSC); /* phi vector */
  CALLOC(tleg, costerms(2*order), double, ON, AMSC);
	       	/* temp legendre storage (2*order needed for local exp) */
  CALLOC(factFac, order+1, double*, ON, AMSC);
  for(x = 0; x < order+1; x++) {
    CALLOC(factFac[x], order+1, double, ON, AMSC);
  }
  evalFactFac(factFac, order);	/* get factorial factors for mulMulti2P */

#if DISSYN == ON
  /* for counts of local/multipole expansions and eval mat builds by level */
  CALLOC(localcnt, depth+1, int, ON, AMSC);
  CALLOC(multicnt, depth+1, int, ON, AMSC);
  CALLOC(evalcnt, depth+1, int, ON, AMSC);
#endif

#if DMTCNT == ON
  /* for counts of transformation matrices by level  */
  CALLOC(Q2Mcnt, depth+1, int*, ON, AMSC);
  CALLOC(Q2Lcnt, depth+1, int*, ON, AMSC);
  CALLOC(Q2Pcnt, depth+1, int*, ON, AMSC);
  CALLOC(L2Lcnt, depth+1, int*, ON, AMSC);
  CALLOC(M2Mcnt, depth+1, int*, ON, AMSC); 
  CALLOC(M2Lcnt, depth+1, int*, ON, AMSC); 
  CALLOC(M2Pcnt, depth+1, int*, ON, AMSC); 
  CALLOC(L2Pcnt, depth+1, int*, ON, AMSC);
  for(x = 0; x < depth+1; x++) {
    CALLOC(Q2Mcnt[x], depth+1, int, ON, AMSC);
    CALLOC(Q2Lcnt[x], depth+1, int, ON, AMSC);
    CALLOC(Q2Pcnt[x], depth+1, int, ON, AMSC);
    CALLOC(L2Lcnt[x], depth+1, int, ON, AMSC);
    CALLOC(M2Mcnt[x], depth+1, int, ON, AMSC);
    CALLOC(M2Lcnt[x], depth+1, int, ON, AMSC);
    CALLOC(M2Pcnt[x], depth+1, int, ON, AMSC);
    CALLOC(L2Pcnt[x], depth+1, int, ON, AMSC);
  }
#endif

  /* from here down could be switched out when the fake dwnwd pass is used */
  CALLOC(facFrA, 2*order+1, double*, ON, AMSC);
  for(x = 0; x < 2*order+1; x++) {
    CALLOC(facFrA[x], 2*order+1, double, ON, AMSC);
  }
  /* generate table of factorial fraction evaluations (for M2L and L2L) */
  evalFacFra(facFrA, order);
  CALLOC(sinmkB, 2*order+1, double, ON, AMSC); /* sin[(m+-k)beta] */
  CALLOC(cosmkB, 2*order+1, double, ON, AMSC); /* cos[(m+-k)beta] */
  cosmkB[0] = 1.0;		/* look up arrays used for local exp */
  /* generate array of sqrt((n+m)!(n-m)!)'s for L2L
  evalSqrtFac(sqrtFac, factFac, order); */
}

/*
 * evalLegendre returns a vector of Legendre function evaluations of the 
 * form Pn^m(cosA) n and m have maximum value = order.
 * Vector entries correspond to (n,m) = (0,0) (1,0) (1,1) (2,0) (2,1)...
 */
void evalLegendre(cosA, vector, order)
double cosA, *vector;
int order;
{
  int x;
  int n, m;			/* as in Pn^m, both <= order */
  double sinMA;			/* becomes sin^m(alpha) in higher order P's */
  double fact;			/* factorial factor */

  /* do evaluations of first four functions separately w/o recursions */
  vector[index(0, 0)] = 1.0;	/* P0^0 */
  if(order > 0) {
    vector[index(1, 0)] = cosA;	/* P1^0 */
    vector[index(1, 1)] = sinMA = -sqrt(1-cosA*cosA); /* P1^1 = -sin(alpha) */
  }
  if(order > 1) vector[index(2, 1)] = 3*sinMA*cosA; /* P2^1 = -3sin()cos() */

  /* generate remaining evaluations by recursion on lower triangular array */
  fact = 1.0;
  for(m = 0; m < order+1; m++) {
    if(m != 0 && m != 1) {	/* set up first two evaluations in row */
      fact *= (2*m - 1); /* (2(m-1)-1)!! -> (2m-1)!! */
      /* use recursion on m */
      if(vector[index(1, 1)] == 0.0) {
	vector[index(m, m)] = 0.0;
	if(m != order) vector[index(m+1, m)] = 0.0;	/* if not last row */
      }
      else {
	cosA = vector[index(1,0)]/vector[index(1,1)]; /* cosA= -cot(theta) */
	sinMA *= vector[index(1,1)]; /*(-sin(alpha))^(m-1)->(-sin(alpha))^m*/
	vector[index(m, m)] = fact * sinMA;
	if(m != order) {		/* do if not on last row */
	  vector[index(m+1, m)] = vector[index(1, 0)]*(2*m+1)
	      *vector[index(m, m)];
	}
      }
    }
    for(x = 2; x < order-m+1; x++) { /* generate row of evals recursively */
      vector[index(x+m, m)] = 
	  ((2*(x+m)-1)*vector[index(1, 0)]*vector[index(x+m-1, m)]
	   - (x + 2*m - 1)*vector[index(x+m-2, m)])/x;
    }
  }
}

static double **smat=NULL;
restartMulti()
{
  smat = NULL;
}

/* 
  Used for the upward pass. 
*/
double **mulQ2Multi(sngs, numsngs, x, y, z, didthis, order, mat)
double x, y, z, **mat; 
snglrty **sngs;
int numsngs, order, didthis;
{
  double **mulMulti2Multi();
  int i, j, k, l, terms = multerms(order);
  snglrty *pq;

  /* Just exchange RHS and LHS translation matrices if subseq pass. */
  if(didthis == TRUE) coeffSwap(sngs, numsngs, terms, mat);
  else { /* Calculate translation matrix for each panel. */

    if(mat != NULL) {
      ASSERT(mat == NULL);
    }
    /* Allocate the matrix for this set of panels. */
    CALLOC(mat, terms, double*, ON, AQ2M);
    for(i=0; i < terms; i++) CALLOC(mat[i], (2 * numsngs), double, ON, AQ2M);

    for(j = 0; j < numsngs; j++) { /* for each panel */
      pq = sngs[j];

      /* Compute Shift matrix, multipole coeffs for panels already done. */
      /* Note, mulMulti2Multi allocates space. */
      smat = mulMulti2Multi(pq->x, pq->y, pq->z, x, y, z, order, smat);

      /* Multiply the shift matrix by multipole coeffs for mono and dipole. */
      for(i = terms-1; i >= 0; i--) {
	mat[i][j] = 0.0;
	mat[i][j+numsngs] = 0.0;
	for(k = terms-1; k >= 0; k--) {
	  mat[i][j] += smat[i][k] * pq->multipoleR[k];
	  mat[i][j+numsngs] += smat[i][k] * pq->multipoleL[k];
	}
      }
    }
  }

#if DALQ2M
  dispQ2M(mat, sngs, numsngs, x, y, z, order);
#endif

  return(mat);
}
  

double **mulMulti2Multi(x, y, z, xp, yp, zp, order, mat)
double x, y, z, xp, yp, zp;	/* cube center, parent cube center */
int order;
double **mat;
{
  double rho, rhoPwr, cosA, beta, mBeta, temp1, temp2; 
  double iPwr(), fact();
  int r, j, k, m, n, c;
  int cterms = costerms(order), sterms = sinterms(order);
  int terms = cterms + sterms;

  /* Allocate the matrix (terms x terms ) if needed. */
  if(mat == NULL) {
    CALLOC(mat, terms, double*, ON, AM2M);
    for(r=0; r < terms; r++) CALLOC(mat[r], terms, double, ON, AM2M);
  }
  else {
    for(r = 0; r < terms; r++)
	for(c = 0; c < terms; c++) mat[r][c] = 0.0;
  }

  /* get relative distance in spherical coordinates */
  xyz2sphere(x, y, z, xp, yp, zp, &rho, &cosA, &beta);

  /* get the requisite Legendre function evaluations */
  evalLegendre(cosA, tleg, order);

  /* for each new moment (Nj^k) stuff the appropriate matrix entries */
  /* done completely brute force, one term at a time; uses exp in nb 12, p29 */
  for(j = 0; j <= order; j++) {
    for(k = 0; k <= j; k++) {
      for(n = 0, rhoPwr = 1.0; n <= j; n++, rhoPwr *= rho) {
	for(m = 0, mBeta = 0.0; m <= n; m++, mBeta += beta) {

	  if(k == 0) {		/* figure terms for Nj^0, ie k = 0 */
	    if(m <= j-n) {	/* if O moments are nonzero */
	      temp1 = fact(j)*rhoPwr*iPwr(2*m)*tleg[index(n, m)];
	      temp1 /= (fact(j-n+m)*fact(n+m));
	      mat[index(j, k)][index(j-n, m)] += temp1*cos(mBeta);
	      if(m != 0)		/* if sin term is non-zero */
		  mat[index(j, k)][sindex(j-n, m, cterms)] += temp1*sin(mBeta);
	    }
	  }
	  else {		/* figure terms for Nj^k, k != 0 */
	    temp1 = fact(j+k)*rhoPwr*tleg[index(n, m)]/fact(n+m);
	    temp2 = temp1*iPwr(2*m)/fact(j-n+k+m);
	    temp1 = temp1*iPwr(k-m-abs(k-m))/fact(j-n+abs(k-m));

	    /* write the cos(kPhi) coeff, bar(N)j^k */
	    if(m != 0) {
	      if(k-m < 0 && abs(k-m) <= j-n) {	/* use conjugates here */
		mat[index(j, k)][index(j-n, m-k)] += temp1*cos(mBeta);
		mat[index(j, k)][sindex(j-n, m-k, cterms)] += temp1*sin(mBeta);
	      }
	      else if(k-m == 0) {	/* double to compensate for 2Re sub. */
		mat[index(j, k)][index(j-n, k-m)] += 2*temp1*cos(mBeta);
		/* sin term is always zero */
	      }
	      else if(k-m > 0 && k-m <= j-n) {
		mat[index(j, k)][index(j-n, k-m)] += temp1*cos(mBeta);
		mat[index(j, k)][sindex(j-n, k-m, cterms)] -= temp1*sin(mBeta);
	      }
	      if(k+m <= j-n) {
		mat[index(j, k)][index(j-n, k+m)] += temp2*cos(mBeta);
		mat[index(j, k)][sindex(j-n, k+m, cterms)] += temp2*sin(mBeta);
	      }
	    }			/* do if m = 0 and O moments not zero */
	    else if(k <= j-n) mat[index(j, k)][index(j-n, k)] += temp2;

	    /* write the sin(kPhi) coeff, dblbar(N)j^k, if it is non-zero */
	    if(m != 0) {
	      if(k-m < 0 && abs(k-m) <= j-n) {	/* use conjugates here */
		mat[sindex(j, k, cterms)][index(j-n, m-k)] += temp1*sin(mBeta);
		mat[sindex(j, k, cterms)][sindex(j-n, m-k, cterms)] 
		    -= temp1*cos(mBeta);
	      }
	      else if(k-m == 0) {/* double to compensate for 2Re sub */
		mat[sindex(j, k, cterms)][index(j-n, k-m)] 
		    += 2*temp1*sin(mBeta);
		/* sine term is always zero */
	      }
	      else if(k-m > 0 && k-m <= j-n) {
		mat[sindex(j, k, cterms)][index(j-n, k-m)] += temp1*sin(mBeta);
		mat[sindex(j, k, cterms)][sindex(j-n, k-m, cterms)] 
		    += temp1*cos(mBeta);
	      }
	      if(k+m <= j-n) {
		mat[sindex(j, k, cterms)][index(j-n, k+m)] -= temp2*sin(mBeta);
		mat[sindex(j, k, cterms)][sindex(j-n, k+m, cterms)] 
		    += temp2*cos(mBeta);
	      }
	    }			/* do if m = 0 and moments not zero */
	    else if(k <= j-n) mat[sindex(j, k,cterms)][sindex(j-n, k, cterms)] 
		+= temp2;
	  }
	}
      }
    }
  }
#if DISM2M == ON
  dispM2M(mat, x, y, z, xp, yp, zp, order);
#endif
  return(mat);
}

/* 
  builds multipole evaluation matrix; used only for fake downward pass
  x,y,z is the multipole coord, fpts[]->x,y,z are the evaluation points.
*/
double **mulMulti2P(x, y, z, fpts, numfpts, order)
double x, y, z;			/* multipole expansion origin */
fieldpt **fpts;
int numfpts, order;
{
  double **mat;
  double cosTh;			/* cosine of elevation coordinate */
  double factorial;		/* 1/factorial = (n-m)!/(n+m)! */
  double tmp, matM1, matM2, matC1, matC2, matP1, matP2;
  int i, j, k, m, n, kold, start;
  int cterms = costerms(order), sterms = sinterms(order);
  int terms = cterms + sterms;

  CALLOC(mat, numfpts, double*, ON, AM2P);
  for(i=0; i < numfpts; i++) 
      CALLOC(mat[i], terms, double, ON, AM2P);

  /* get Legendre function evaluations, one set for each snglrty */
  /*   also get fieldpoint coordinates to set up rest of matrix */
  for(i = 0; i < numfpts; i++) { /* for each fieldpt, do a legendre eval set */
    xyz2sphere(fpts[i]->x, fpts[i]->y, fpts[i]->z,
	       x, y, z, &(Ir[i]), &cosTh, &(phi[i]));

    Irn[i] = Ir[i]; /* initialize (1/r)^n+1 vec. */
    Mphi[i] = phi[i];		/* initialize m*phi vector */

    evalLegendre(cosTh, mat[i], order);	/* wr moms to 1st (cos) half of row */

  }

#if DALM2P == ON
  fprintf(stdout,
	  "\nM2P MATRIX BUILD:\n    AFTER LEGENDRE FUNCTION EVALUATON\n");
  dumpMat(mat, numfpts, terms);
#endif

  /* add the (1/r)^n+1 factors to the left (cos(m*phi)) half of the matrix */
  for(j = 0, k = kold = 1; j < cterms; j++) { /* loop on columns of matrix */
    for(i = 0; i < numfpts; i++) mat[i][j] /= Irn[i]; /* divide by r^n+1 */
    k -= 1;
    if(k == 0) {		/* so that n changes as appropriate */
      kold = k = kold + 1;
      for(i = 0; i < numfpts; i++) Irn[i] *= Ir[i]; /* r^n -> r^n+1 */
    }
  }

#if DALM2P == ON
  fprintf(stdout,
	  "    AFTER ADDITION OF (1/R)^N+1 FACTORS\n");
  dumpMat(mat, numfpts, terms);
#endif

  /* add the factorial fraction factors to the left (cos(m*phi)) part of mat */
  /*  note that (n-m)!/(n+m)! = 1/(n-m+1)...(n+m) since m \leq n */
  for(n = 1; n <= order; n++) {
    for(m = 1; m <= n; m++) {
      for(i = 0; i < numfpts; i++) mat[i][index(n, m)] /= factFac[n][m];
    }
  }

#if DALM2P == ON
  fprintf(stdout,
	  "    AFTER ADDITION OF FACTORIAL FRACTION FACTORS\n");
  dumpMat(mat, numfpts, terms);
#endif

  /* copy left half of matrix to right half for sin(m*phi) terms */
  for(i = 0; i < numfpts; i++) { /* loop on rows of matrix */
    for(n = 1; n <= order; n++) { 
      for(m = 1; m <= n; m++) {	/* copy a row */
	mat[i][sindex(n, m, cterms)] = mat[i][index(n, m)];
      }
    }
  }

#if DALM2P == ON
  fprintf(stdout,
	  "    AFTER COPYING SINE (RIGHT) HALF\n");
  dumpMat(mat, numfpts, terms);
#endif

  /* add factors of cos(m*phi) and sin(m*phi) to left and right halves resp. */
  for(m = 1; m <= order; m++) {	/* lp on m in Mn^m (no m=0 since cos(0)=1) */
    for(n = m; n <= order; n++) { /* loop over cols with same m */
      for(i = 0; i < numfpts; i++) { /* add factors to a column */
	mat[i][index(n, m)] *= cos(Mphi[i]);
	mat[i][sindex(n, m, cterms)] *= sin(Mphi[i]);
      }
    }
    for(i = 0; i < numfpts; i++) Mphi[i] += phi[i]; /* (m-1)*phi->m*phi */
  }


  /*  
   * jt: modifications for normal derivative evaluations.
   */

/*cdel nr_mulMulti_calls++;*/

  for(i = 0; i < numfpts; i++)
    if(fpts[i]->deriv == 1) {
      for(n = 0; n < order; n++) {
        /* m=0 */
        factorial = (n+2)*(n+1);
        matP1 = factorial*mat[i][index(n+1, 1)];
        matP2 = factorial*mat[i][sindex(n+1, 1, cterms)];
        matC1 = (n+1)*mat[i][index(n+1, 0)];
        tmp  = fpts[i]->nrm[0]*matP1;
        tmp += fpts[i]->nrm[1]*matP2;
        tmp -= fpts[i]->nrm[2]*matC1;
        mat[i][index(n, 0)] = tmp;
        for(m = 1; m <= n; m++) {
          matM1 = mat[i][index(n+1, m-1)];
          if (m == 1)
            matM2 = 0.0;
          else
            matM2 = mat[i][sindex(n+1, m-1, cterms)];

          factorial = (n+m+2)*(n+m+1);
          matP1 = factorial*mat[i][index(n+1, m+1)];
          matP2 = factorial*mat[i][sindex(n+1, m+1, cterms)];
          matC1 = (n+1+m)*mat[i][index(n+1, m)];
          matC2 = (n+1+m)*mat[i][sindex(n+1, m, cterms)];            

          tmp  = fpts[i]->nrm[0]*0.5*(matP1-matM1);
          tmp += fpts[i]->nrm[1]*0.5*(matP2+matM2);
          tmp -= fpts[i]->nrm[2]*matC1;
          mat[i][index(n, m)] = tmp;
          tmp  = fpts[i]->nrm[0]*0.5*(matP2-matM2);
          tmp -= fpts[i]->nrm[1]*0.5*(matP1+matM1);
          tmp -= fpts[i]->nrm[2]*matC2;
          mat[i][sindex(n, m, cterms)] = tmp;
        }
      }
      mat[i][index(order, 0)] = 0.0;
      for(m = 1; m <= order ; m++) {
        mat[i][index(order, m)] = 0.0;
        mat[i][sindex(order, m, cterms)] = 0.0;
      }
    }


#if DISM2P == ON
  dispM2P(mat, x, y, z, fpts, numfpts, order);
#endif

  return(mat);
}
