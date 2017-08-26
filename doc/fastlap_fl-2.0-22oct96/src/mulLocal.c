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

/*cdel extern int nr_mulLocal_calls, nr_mulMulti_calls; */

/*
  globals used for temporary storage
*/
double **facFrA;		/* array of factorial fractions */
double *cosmkB;			/* array used to look up cos[(m+-k)beta] */
double *sinmkB;			/* array used to look up sin[(m+-k)beta] */

/*
  initializes the factorial fraction array used in M2L, L2L matrix calculation
*/
void evalFacFra(array, order)
int order;			/* array is 2*order+1 x 2*order+1 */
double **array;			/* array[num][den] = num!/den! */
{
  int d, i;
  for(i = 0; i <= 2*order; i++) {
    array[i][i] = 1.0; /* do main diagonal */
    if(i > 0 && i < 2*order) array[i+1][i] = i+1; /* do first sub diagonal */
  }
  for(d = 3; d <= 2*order; d++) { /* loop on lower triangular rows */
    for(i = 1; i < d-1; i++) {	/* loop on columns */
      array[d][i] = array[d-1][i]*array[d][d-1];
    }
  }
  /* invert lower part entries and copy to top */
  for(d = 2; d <= 2*order; d++) {
    for(i = 1; i <= d-1; i++) {
      array[i][d] = 1/array[d][i];
    }
  }
  /* copy 1st row and column from computed values */
  for(d = 1; d <= 2*order; d++) {
    array[0][d] = array[1][d];
    array[d][0] = array[d][1];
  }

#if DISFAC == ON
  fprintf(stdout, "FACTORIAL FRACTION ARRAY:\n");
  dumpMat(array, 2*order+1, 2*order+1);
#endif

}

/*
  initializes sqrt((m+n)!/(n-m)!) lookup table (for L2L)
*/
void evalSqrtFac(arrayout, arrayin, order)
int order;
double **arrayout, **arrayin;
{
  int n, m;			/* arrayout[n][m] = sqrt((m+n)!/(n-m)!) */

  /* set up first column, always ones */
  for(n = 0; n < order+1; n++) arrayout[n][0] = 1.0;

  /* set up lower triangular (n+m)!/(n-m)! */
  for(n = 1; n <= order; n++) {
    for(m = 1; m <= n; m++) {
      arrayout[n][m] = sqrt(arrayin[n][m]);
    }
  }

#if DISSFA == ON
  fprintf(stdout, "SQUARE ROOT FACTORIAL ARRAY:\n");
  dumpMat(arrayout, order+1, order+1);
#endif

}


/*
  initializes cos[(m+-k)beta] and sin[(m+-k)beta] lookup tables (M2L and L2L)
*/
void evalSinCos(beta, order)
int order;
double beta;
{
  int i;
  double temp = beta;

  for(i = 1; i <= 2*order; beta += temp, i++) {
    sinmkB[i] = sin(beta);
    cosmkB[i] = cos(beta);
  }
}

/*
  looks up sin[(m+-k)beta]
*/
double sinB(sum)
int sum;
{
  if(sum < 0) return(-sinmkB[abs(sum)]);
  else return(sinmkB[sum]);
}

/*
  looks up cos[(m+-k)beta]
*/
double cosB(sum)
int sum;
{
  return(cosmkB[abs(sum)]);
}

/* 
  Used for all but no local downward pass. 
*/
double **mulMulti2Local(x, y, z, xp, yp, zp, order, mat)
int order;
double x, y, z, xp, yp, zp;	/* multipole and local cube centers */
double **mat;			/* Place to store matrix, typically null. */
{
  int i, j, k, n, m;
  int terms = multerms(order);	/* the number of non-zero moments */
  int ct = costerms(order);	/* the number of non-zero cos (bar) moments */
  double rho, cosA, beta;	/* spher. position of multi rel to local */
  double rhoJ, rhoN;		/* rho^j and (-1)^n*rho^(n+1) in main loop */
  double rhoFac;		/* = rhoJ*rhoN intermediate storage */
  double temp1, temp2, temp3;
  double iPwr(), sinB(), cosB();
  extern double *tleg, *Ir, *Irn, *phi, *Mphi; /* external temporary storage */

  /* allocate a zeroed multi to local transformation matrix */
  if(mat == NULL) {
    CALLOC(mat, terms, double*, ON, AM2L);
    for(i = 0; i < terms; i++)
	CALLOC(mat[i], terms, double, ON, AM2L);
  }
  else {
    for(i = 0; i < terms; i++)
	for(j = 0; j < terms; j++) mat[i][j] = 0.0;
  }

  /* find relative spherical coordinates */
  xyz2sphere(x, y, z, xp, yp, zp, &rho, &cosA, &beta);

  /* generate legendre function evaluations */
  evalLegendre(cosA, tleg, 2*order); /* multi->loc needs 2x legendres */

  /* generate sin[(m+-k)beta] and cos[(m+-k)beta] look up arrays */
  /*  other lookup arrays generated in mulMultiAlloc() */
  evalSinCos(beta, order);

  /* generate multi to local transformation matrix; uses NB12 pg30 */
  /*  rhoFac factor divides could be reduced to once per loop */
  for(j = 0, rhoJ = 1.0; j <= order; rhoJ *= rho, j++) {
    for(k = 0; k <= j; k++) {	/* loop on Nj^k's, local exp moments */
      for(n = 0, rhoN = rho; n <= order; rhoN *= (-rho), n++) {
	rhoFac = 1.0 / (rhoJ*rhoN);  /* gives (-1)^n/rho^(j+n+1) factor */
	for(m = 0; m <= n; m++) {  /* loop on On^m's, multipole moments */

	  /* generate a bar(N)j^k and dblbar(N)j^k entry */
	  if(k == 0) {	       /* use abbreviated formulae in this case */

	    /* generate only bar(N)j^0 entry (dblbar(N)j^0 = 0 always) */
	    if(m != 0) {
	      temp1 = tleg[index(j+n, m)]*facFrA[j+n-m][n+m];
	      mat[index(j, 0)][index(n, m)] += temp1*cosB(m)*rhoFac;
	      mat[index(j, 0)][sindex(n, m, ct)] += temp1*sinB(m)*rhoFac;
	    }
	    else mat[index(j, 0)][index(n, 0)] 
		+= tleg[index(j+n, 0)]*facFrA[j+n][n]*rhoFac;
	  }
	  else {
	    temp1 = tleg[index(j+n, abs(m-k))]
		*facFrA[j+n-abs(m-k)][n+m]*iPwr(abs(k-m)-k-m);
	    temp2 = tleg[index(j+n, m+k)]*facFrA[j+n-m-k][n+m];
	    temp3 = tleg[index(j+n, k)]*facFrA[j+n-k][n]*2;

	    /* generate bar(N)j^k entry */
	    if(m != 0) {
	      mat[index(j, k)][index(n, m)] 
		  += (temp1*cosB(m-k)+temp2*cosB(m+k))*rhoFac;
	      mat[index(j, k)][sindex(n, m, ct)] 
		  += (temp1*sinB(m-k)+temp2*sinB(m+k))*rhoFac;
	    }
	    else mat[index(j, k)][index(n, 0)] += temp3*cosB(k)*rhoFac;

	    /* generate dblbar(N)j^k entry */
	    if(m != 0) {
	      mat[sindex(j, k, ct)][index(n, m)] 
		  += (-temp1*sinB(m-k)+temp2*sinB(m+k))*rhoFac;
	      mat[sindex(j, k, ct)][sindex(n, m, ct)] 
		  += (temp1*cosB(m-k)-temp2*cosB(m+k))*rhoFac;
	    }
	    else mat[sindex(j, k, ct)][index(n, 0)] += temp3*sinB(k)*rhoFac;
	  }
	}
      }
    }
  }

#if DISM2L == ON
  dispM2L(mat, x, y, z, xp, yp, zp, order);
#endif

  return(mat);
}

/* 
  Used only for true (Greengard) downward pass - similar to Multi2Local
*/
double **mulLocal2Local(x, y, z, xc, yc, zc, order)
int order;
double x, y, z, xc, yc, zc;	/* parent and child cube centers */
{
  int i, j, k, n, m;
  int terms = multerms(order);	/* the number of non-zero moments */
  int ct = costerms(order);	/* the number of non-zero cos (bar) moments */
  double **mat;			/* the transformation matrix */
  double rho, cosA, beta;	/* spher. position of multi rel to local */
  double rhoJ, rhoN;		/* rho^j and (-1)^n*rho^(n+1) in main loop */
  double rhoFac;		/* = rhoJ*rhoN intermediate storage */
  double temp1, temp2, temp3;
  double iPwr(), sinB(), cosB();
  extern double *tleg, *Ir, *Irn, *phi, *Mphi; /* external temporary storage */

  /* allocate the local to local transformation matrix */
  CALLOC(mat, terms, double*, ON, AL2L);
  for(i = 0; i < terms; i++)
      CALLOC(mat[i], terms, double, ON, AL2L);

  /* find relative spherical coordinates */
  xyz2sphere(x, y, z, xc, yc, zc, &rho, &cosA, &beta);

  /* generate legendre function evaluations */
  evalLegendre(cosA, tleg, 2*order); /* local->local needs 2x legendres */

  /* generate sin[(m+-k)beta] and cos[(m+-k)beta] look up arrays */
  /*  other lookup arrays generated in mulMultiAlloc() */
  evalSinCos(beta, order);

  /* generate local to local transformation matrix; uses NB12 pg36Y */
  /*  rhoFac factor divides could be reduced to once per loop */
  for(j = 0, rhoJ = 1.0; j <= order; rhoJ *= (-rho), j++) {
    for(k = 0; k <= j; k++) {	/* loop on Nj^k's, local exp moments */
      for(n = j, rhoN = rhoJ; n <= order; rhoN *= (-rho), n++) {
	for(m = 0; m <= n; m++) { /* loop on On^m's, old local moments */

	  /* generate a bar(N)j^k and dblbar(N)j^k entry */
	  rhoFac = rhoN/rhoJ;	/* divide to give (-rho)^(n-j) factor */
	  if(k == 0 && n-j >= m) {  /* use abbreviated formulae in this case */

	    /* generate only bar(N)j^0 entry (dblbar(N)j^0 = 0 always) */
	    if(m != 0) {
	      temp1 = tleg[index(n-j, m)]*facFrA[0][n-j+m]*rhoFac;
	      mat[index(j, 0)][index(n, m)] += temp1*cosB(m);
	      mat[index(j, 0)][sindex(n, m, ct)] += temp1*sinB(m);
	    }
	    else mat[index(j, 0)][index(n, 0)] += tleg[index(n-j, 0)]
		    *facFrA[0][n-j]*rhoFac;
	  }
	  else {
	    if(n-j >= abs(m-k)) temp1 = tleg[index(n-j, abs(m-k))]
		*facFrA[0][n-j+abs(m-k)]*iPwr(m-k-abs(m-k))*rhoFac;
	    if(n-j >= m+k) temp2 = tleg[index(n-j, m+k)]
		*facFrA[0][n-j+m+k]*iPwr(2*k)*rhoFac;
	    if(n-j >= k) temp3 = 2*tleg[index(n-j, k)]
		*facFrA[0][n-j+k]*iPwr(2*k)*rhoFac;

	    /* generate bar(N)j^k entry */
	    if(m != 0) {
	      if(n-j >= abs(m-k)) {
		mat[index(j, k)][index(n, m)] += temp1*cosB(m-k);
		mat[index(j, k)][sindex(n, m, ct)] += temp1*sinB(m-k);
	      }
	      if(n-j >= m+k) {
		mat[index(j, k)][index(n, m)] += temp2*cosB(m+k);
		mat[index(j, k)][sindex(n, m, ct)] += temp2*sinB(m+k);
	      }
	    }
	    else if(n-j >= k) mat[index(j, k)][index(n, 0)] += temp3*cosB(k);

	    /* generate dblbar(N)j^k entry */
	    if(m != 0) {
	      if(n-j >= abs(m-k)) {
		mat[sindex(j, k, ct)][index(n, m)] += (-temp1*sinB(m-k));
		mat[sindex(j, k, ct)][sindex(n, m, ct)] += temp1*cosB(m-k);
	      }
	      if(n-j >= m+k) {
		mat[sindex(j, k, ct)][index(n, m)] += (-temp2*sinB(m+k));
		mat[sindex(j, k, ct)][sindex(n, m, ct)] += temp2*cosB(m+k);
	      }
	    }
	    else if(n-j >= k) 
		mat[sindex(j, k, ct)][index(n, 0)] += temp3*sinB(k);
	  }
	}
      }
    }
  }

#if DISL2L == ON
  dispL2L(mat, x, y, z, xc, yc, zc, order);
#endif

  return(mat);
}

/* We loop over the panels in the subject cube and find the shift 
   matrix (via mulMulti2Local) to shift the snglrty multipole expansions 
   to a local expansion centered at x,y,z. The vectors for the panel
   with both RHS and LHS coefficients are returned in mat. */

static double **shft_mat=NULL;
restartLocal()
{
  shft_mat = NULL;
}


double **mulQ2Local(sngs, numsngs, x, y, z, order, didthis, mat)
double x, y, z;     /* This is the center for the local expansion. */
snglrty **sngs;     /* Array of pointers to the snglrty structures. */
int numsngs, order; /* Number of snglrtys in this cube, order of expansions. */
int didthis;        /* mat contains translation. */
double **mat;       /* contains translation, usually NULL. */
{
  int i, j, k;
  int terms = multerms(order);
  snglrty *pq;

  if(didthis == TRUE) {
    coeffSwap(sngs, numsngs, terms, mat);
  }
  else {
    /* Allocate the matrix we return. */
    CALLOC(mat, terms, double*, ON, AQ2L);
    for (i = 0; i < terms; i++) {
      CALLOC(mat[i], (2*numsngs), double, ON, AQ2L); /* It's 2* as we provide
							for both RHS and LHS
							coefficients. */
    }

    for (j = 0; j < numsngs; j++) { /* loop over each panel.*/
      pq = sngs[j];

      /* Get the matrix. Note, multi2local allocates shft_mat first time. */
      shft_mat = mulMulti2Local(pq->x, pq->y, pq->z, x, y, z, order, shft_mat);

      /* Multiply the shift matrix by multipole coeffs for mono and dipole. */
      for(i = terms-1; i >= 0; i--) {
	mat[i][j] = 0.0;
	mat[i][j+numsngs] = 0.0;
	for(k = terms-1; k >= 0; k--) {
	  mat[i][j] += shft_mat[i][k] * pq->multipoleR[k];
	  mat[i][j+numsngs] += shft_mat[i][k] * pq->multipoleL[k];
	}
      }
    }
  }
#if DALQ2L == ON
  fprintf(stdout,"    Q2L MATRIX\n");
  dumpMat(mat, terms, numsngs);
#endif
  
#if DISQ2L == ON
  dispQ2L(mat, sngs, numsngs, x, y, z, order);
#endif

  return(mat);
}

/*
  builds local expansion evaluation matrix; not used for fake dwnwd pass
  follows NB10 equation marked circle(2A) except roles of j,k and n,m 
  switched very similar to mulMulti2P()
*/
double **mulLocal2P(x, y, z, fpts, numfpts, order)
double x, y, z;
fieldpt **fpts;
int numfpts, order;
{
  double **mat;
  double cosTh;			/* cosine of elevation coordinate */
  double fact();
  extern double *Irn, *Mphi, *phi, *Ir;
  int i, j, k, m, n, kold, start;
  int cterms = costerms(order), terms = multerms(order);
  double tmp, matM1, matP1, matC1, matM2, matP2, matC2;

  CALLOC(mat, numfpts, double*, ON, AL2P);
  for(i = 0; i < numfpts; i++) 
      CALLOC(mat[i], terms, double, ON, AL2P);

  /* get Legendre function evaluations, one set for each snglrty */
  /*   also get snglrty coordinates to set up rest of matrix */
  for(i = 0; i < numfpts; i++) { /* for each fieldpt, do a legendre eval set */
    xyz2sphere(fpts[i]->x, fpts[i]->y, fpts[i]->z,
	       x, y, z, &(Ir[i]), &cosTh, &(phi[i]));
    Irn[i] = 1.0; /* initialize r^n vec. */
    Mphi[i] = phi[i];		/* initialize m*phi vector */
    evalLegendre(cosTh, mat[i], order);	/* wr moms to 1st (cos) half of row */
  }

#if DALL2P == ON
  fprintf(stdout,
	  "\nL2P MATRIX BUILD:\n    AFTER LEGENDRE FUNCTION EVALUATON\n");
  dumpMat(mat, numfpts, terms);
#endif

  /* add the r^n factors to the left (cos(m*phi)) half of the matrix */
  for(j = 0, k = kold = 1; j < cterms; j++) { /* loop on columns of matrix */
    for(i = 0; i < numfpts; i++) mat[i][j] *= Irn[i]; /* multiply by r^n */
    k -= 1;
    if(k == 0) {		/* so that n changes as appropriate */
      kold = k = kold + 1;
      for(i = 0; i < numfpts; i++) Irn[i] *= Ir[i]; /* r^n -> r^n+1 */
    }
  }

#if DALL2P == ON
  fprintf(stdout,"    AFTER ADDITION OF R^N FACTORS\n");
  dumpMat(mat, numfpts, terms);
#endif

  /* add the factorial fraction factors to the left (cos(m*phi)) part of mat */
  for(n = 0; n <= order; n++) {
    for(m = 0; m <= n; m++) {
      for(i = 0; i < numfpts; i++) mat[i][index(n, m)] /= fact(n+m);
    }
  }

#if DALL2P == ON
  fprintf(stdout,"    AFTER ADDITION OF FACTORIAL FRACTION FACTORS\n");
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

#if DALL2P == ON
  fprintf(stdout,"    AFTER COPYING SINE (RIGHT) HALF\n");
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
   *  jt: modifications for normal derivative evaluations. 
   */

/*cdel nr_mulLocal_calls++;*/

  for(i = 0; i < numfpts; i++) {
    if(fpts[i]->deriv == 1) {
      for(n = order; n >= 2; n--) {
        /* m=0 */
        matP1 = mat[i][index(n-1, 1)];
        matP2 = mat[i][sindex(n-1, 1, cterms)];
        matC1 = mat[i][index(n-1, 0)];
        tmp  = fpts[i]->nrm[0]*matP1;
        tmp += fpts[i]->nrm[1]*matP2;
        tmp += fpts[i]->nrm[2]*matC1;
        mat[i][index(n, 0)] = tmp;
        for(m = 1; m <= n ; m++) {
          matM1 = mat[i][index(n-1, m-1)];
          if (m == 1)
            matM2 = 0.0;            
          else
            matM2 = mat[i][sindex(n-1, m-1, cterms)];

          if (m == n) {
            matC1 = 0.0;
            matC2 = 0.0;
            matP1 = 0.0;
            matP2 = 0.0;
          }
          else {
            matC1 = mat[i][index(n-1, m)];
            matC2 = mat[i][sindex(n-1, m, cterms)];
            if (m == n-1) {
              matP1 = 0.0;
              matP2 = 0.0;
            }
            else {
              matP1 = mat[i][index(n-1, m+1)];
              matP2 = mat[i][sindex(n-1, m+1, cterms)];
            }
          }
          tmp  = fpts[i]->nrm[0]*0.5*(matP1-matM1);
          tmp += fpts[i]->nrm[1]*0.5*(matP2+matM2);
          tmp += fpts[i]->nrm[2]*matC1;
          mat[i][index(n, m)] = tmp;
          tmp  = fpts[i]->nrm[0]*0.5*(matP2-matM2);
          tmp -= fpts[i]->nrm[1]*0.5*(matP1+matM1);
          tmp += fpts[i]->nrm[2]*matC2;
          mat[i][sindex(n, m, cterms)] = tmp;
          }
      }
      /* n = 1 */
      matC1 = mat[i][index(0, 0)];
      mat[i][index(1, 0)] = fpts[i]->nrm[2]*matC1;
      mat[i][index(1, 1)] = -0.5*fpts[i]->nrm[0]*matC1;
      mat[i][sindex(1, 1, cterms)] = -0.5*fpts[i]->nrm[1]*matC1;
      /* n = 0 */
      mat[i][index(0,0)] = 0.0;
    }
  }    

        
#if DISL2P == ON
  dispL2P(mat, x, y, z, fpts, numfpts, order);
#endif

  return(mat);
}

