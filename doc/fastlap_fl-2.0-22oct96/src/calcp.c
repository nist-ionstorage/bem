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
* Modifications for normal derivative evaluations in calcp() by Johannes Tausch
*
*/
#define JTDEBUG 0
/*cdel extern int nDiagOnDiel; */       /* my own variable for m-dielectrics */  
/*cdel extern double lambda; */         /* my own variable for m-dielectrics */  
#include <stdio.h>
#include <math.h>

#include "mulStruct.h"
#include "mulGlobal.h"

#define FALSE 0
#define TRUE 1
#define XI 0
#define YI 1
#define ZI 2
#define EQUIV_TOL 1.0e-9
#define PLANAR_TOL 1.0e-3
#define MAX_REAL 1.0e+20;

/* Obvious Constants. */
#define PI 3.1415927
#define TWOPI 6.2831853
#define HALFPI 1.5707963E+00

/* Constants to save typing. */
#define FIVE3 1.666666666667
#define SEVEN3 2.3333333333333
#define ONE6 0.16666666666667
#define ONE3 0.3333333333333
#define FT3 4.66666667


/* Defines breakpoints panel multipoles. */
#define LIMITFOURTH 9.0
#define LIMITSECOND 36.0

/* Constants used for Hart arctan approximation. */
#define B1 0.24091197
#define B2 3.7851122
#define B3 5.6770721
#define B4 5.6772854
#define B5 5.6770747

#define Dot_Product(V1,V2) V1[XI]*V2[XI]+V1[YI]*V2[YI]+V1[ZI]*V2[ZI]
#define DotP_Product(V1,R,S,T) (V1[XI])*(R)+(V1[YI])*(S)+(V1[ZI])*(T)

#if JTDEBUG
extern int printit;
int jtflag, calcflag;
#endif


static int num2nd=0, num4th=0, numexact=0;
static int num2ndsav=0, num4thsav=0, numexactsav=0;

/*
initcalcp does a lot.  In particular, it generates the panel's coordinate
system, and re-labels the panel corners in that coordinate system.  Care
must be taken in exactly how this coordinate system is generated to make
sure it is a right-handed coordinate system.  In particular:

The panel's normal is assumed to be pointing at your nose if you are
looking down on a panel and its corners 0, 1, 2, 3 are labeled clockwise.

*/
initcalcp(sing_list, order)
  snglrty *sing_list;
  int order;
{
  snglrty *pq;
  double vtemp[3];
  double length, maxlength, minlength, length20, length31, sum, sum2, delta;
  double normalize();
  int i, j, next;

  pq = sing_list; 
  while (pq != NULL) {
    /* Differentiate between point sings and panels. */
    if (pq->shape == POINT ) {
      /* Load the centroid which for a point sing is the first vertex. */
      pq->x = pq->corner[0][XI];
      pq->y = pq->corner[0][YI];
      pq->z = pq->corner[0][ZI];
      /* Allocate and compute the expansion coefficients. */
      calcSphericalPoint(pq, order);
    }
    else {
      /* Calculate edge lengths. */
      maxlength = 0.0;
      minlength = MAX_REAL;
      for(i=0; i < pq->shape; i++) {    
	if(i == (pq->shape -1)) next = 0;
	else next = i + 1;
	for(sum= 0, j = 0; j < 3; j++) {
	  delta = pq->corner[next][j] - pq->corner[i][j];
	  sum += delta * delta;
	}
	pq->length[i] = length = sqrt(sum);
	maxlength = MAX(maxlength, length);
	minlength = MIN(minlength, length);
      }
      
      /* Get diags and lengths. */
      for(sum= 0.0, sum2 = 0.0, i = 0; i < 3; i++) {     
	pq->X[i] = delta = pq->corner[2][i] - pq->corner[0][i];
	sum += delta * delta;
	if(pq->shape == 3) pq->Y[i] = pq->corner[1][i] - pq->corner[0][i];      
	else {
	  pq->Y[i] = delta = pq->corner[1][i] - pq->corner[3][i];      
	  sum2 += delta * delta;
	}
      }
      length20 = sqrt(sum);
      length31 = sqrt(sum2);
      
      /* Check on lengths for quad. */
      if(pq->shape == 3) {
	pq->max_diag = maxlength;
	pq->min_diag = minlength;
      }
      else {
	length = MAX(length20, length31);
	/*
	  if( maxlength > length) {
	  fprintf(stderr, "FLE-initcalcp: Convex panel in input.\n");
	  exit(0);
	  }
	  */
	pq->max_diag = length;
	pq->min_diag = MIN(length20, length31);
      }
      
      /* Z-axis is normal to two diags. */
      Cross_Product(pq->X, pq->Y, pq->Z);
      pq->area = 0.5 * normalize(pq->Z);
      normalize(pq->X);
      
      /* Real Y-axis is normal to X and Z. */
      Cross_Product(pq->Z, pq->X, pq->Y);
      
      /* Project the corner points into the plane defined by edge midpoints. */
      if(planarize(pq) == FALSE) {     
	/* Planarization exceeds  PLANAR_TOL, write diagnostic. */
#if NOWARN == OFF
	fprintf(stdout, "FLW-initcalcp: Panel skewed beyond the PLANAR_TOL tolerance.\n");
#endif
      }

      /* Calculate the centroid. */
      centroid(pq, length20);      


      /* jt: define the directions pq->X and pq->Y a la Newman;
         this is necessary for calculating the gradient of the potential 
         in calcp()
         */
      if ( pq->shape == 3 ) {
        for ( i=0; i<3; i++ )
          pq->X[i] = pq->corner[0][i] + pq->corner[0][i]
            - pq->corner[1][i] - pq->corner[2][i] ;
      }
      else {
        for ( i=0; i<3; i++ )
          pq->X[i] = pq->corner[1][i] + pq->corner[0][i]
            - pq->corner[3][i] - pq->corner[2][i];
      }

      normalize(pq->X);
      /* Y-axis is normal to two X and Z direction */
      Cross_Product(pq->Z, pq->X, pq->Y);      

      /* Put corners in the newly defined coord system. */
      for(i=0; i < pq->shape; i++) {
	pq->corner[i][XI] -= pq->x;
	pq->corner[i][YI] -= pq->y;
	pq->corner[i][ZI] -= pq->z;
      }
      for(i=0; i < pq->shape; i++) {
	vtemp[XI] = Dot_Product(pq->corner[i], pq->X);
	vtemp[YI] = Dot_Product(pq->corner[i], pq->Y);
	vtemp[ZI] = Dot_Product(pq->corner[i], pq->Z);
	VCOPY(pq->corner[i], vtemp);
	if(fabs(pq->corner[i][ZI]) > (EQUIV_TOL * pq->min_diag)) {
	  fprintf(stderr, "FLE-initcalcp: renormalized z=%g\n", 
		  pq->corner[i][ZI]);
	  exit(0);
	}
	pq->corner[i][ZI] = 0.0;
      }
      /* Finally, calculate the moments and multipole coefficients. 
	 ComputeMoments calls calcSphericalPanel. */
      ComputeMoments(pq, order);
    }
    /* Iterate for the next singularity. */
    pq = pq->next;
  }
}




/*
Changes the corner points so that they lie in the plane defined by the
panel diagonals and any midpoint of an edge.
*/
planarize(pq)
snglrty *pq;
{
  double origin[3], corner[3], delta[4][3], px, py, dx, dy, dz;
  int i, j, numcorners = pq->shape;
  double tolsq = PLANAR_TOL * pq->min_diag; 

  tolsq *= tolsq;

  /* Triangular panels are just fine already. */
  if(numcorners != 4) return(TRUE);

  /* Pick edge midpoint as origin. */
  for(i=0; i < 3; i++) origin[i] = 0.5 * (pq->corner[1][i] + pq->corner[0][i]);

  for(i=0; i < numcorners; i++) {
    for(j=0; j < 3; j++) corner[j] = pq->corner[i][j] - origin[j];
    px = Dot_Product(corner, pq->X);
    py = Dot_Product(corner, pq->Y);

    dx = px * pq->X[XI] + py * pq->Y[XI] + origin[XI] - pq->corner[i][XI];
    dy = px * pq->X[YI] + py * pq->Y[YI] + origin[YI] - pq->corner[i][YI];
    dz = px * pq->X[ZI] + py * pq->Y[ZI] + origin[ZI] - pq->corner[i][ZI];

    delta[i][XI] = dx;
    delta[i][YI] = dy;
    delta[i][ZI] = dz;
  }
  
  for(i=0; i < numcorners; i++) {
    for(j=0; j < 3; j++) {
      pq->corner[i][j] += delta[i][j];
    }
  }
  /* If moved beyond the tolerance, set flag. */
  if((dx * dx + dy * dy + dz * dz) > tolsq) {
    return(FALSE);
  }
  else {
    return(TRUE);
  }
}


/* 
Determines centroid of a panel (defined as the point which makes the
first moments vanish.  Calculation begins by projection into the
coordinate system defined by the panel normal as the z-axis and
edge02 as the x-axis.
*/
centroid(pp, x2)
snglrty *pp;
double x2;
{
  double vertex1[3], vertex3[3];
  double sum, dl, x1, y1, x3, y3, xc, yc;
  int i;

  /* Use vertex 0 as the origin. */
  for(i=0; i< 3; i++) {
    vertex1[i] = pp->corner[1][i] - pp->corner[0][i];
    if(pp->shape == 4) vertex3[i] = pp->corner[3][i] - pp->corner[0][i];
    else vertex3[i] = pp->corner[2][i] - pp->corner[0][i];
  }

  /* Project into the panel axes. */
  y1 = Dot_Product(vertex1, pp->Y);
  y3 = Dot_Product(vertex3, pp->Y);
  x1 = Dot_Product(vertex1, pp->X);
  x3 = Dot_Product(vertex3, pp->X);

  yc = ONE3 * (y1 + y3);
  xc = ONE3 * (x2 + ((x1 * y1 - x3 * y3)/(y1 - y3)));

  pp->x = pp->corner[0][XI] + xc * pp->X[XI] + yc * pp->Y[XI];
  pp->y = pp->corner[0][YI] + xc * pp->X[YI] + yc * pp->Y[YI];
  pp->z = pp->corner[0][ZI] + xc * pp->X[ZI] + yc * pp->Y[ZI];

}

double normalize(vector)
double vector[3];
{
  double length;
  int i;

  length = sqrt( vector[0]*vector[0] 
		+ vector[1]*vector[1] 
		+ vector[2]*vector[2]);
    
  for (i=0; i<3; i++) vector[i] = vector[i] / length;

  return length;
}

/* Assumes the vectors are normalized. */
int If_Equal(vector1, vector2)
  double vector1[3], vector2[3];
{
  int i;

  for (i=0; i<3; i++)
    if (fabs(vector1[i] - vector2[i]) > EQUIV_TOL) return FALSE;   
  return TRUE;
}

/* Calculates the cross product between two vectors. */
Cross_Product(vector1, vector2, result_vector)
  double vector1[], vector2[], result_vector[];
{
  result_vector[XI] = vector1[YI]*vector2[ZI] - vector1[ZI]*vector2[YI];
  result_vector[YI] = vector1[ZI]*vector2[XI] - vector1[XI]*vector2[ZI];
  result_vector[ZI] = vector1[XI]*vector2[YI] - vector1[YI]*vector2[XI];
}


double tilelength(nq)
snglrty *nq;
{
  return nq->max_diag;
}

/* 
Used to make sure that arrays for computeMoments and Calcspherical are
recomputed. 
*/
static int maxorderMom = -1;
static int maxorderSphere = -1;
restartCalcp()
{
  maxorderMom = -1;
  maxorderSphere = -1;
}

/*
ComputeMoments evaluates moments of quadrilateral surface relative to
local system, array S(15).  Note that S(2)=S(6)=0 because panel is planar.
*/

ComputeMoments(panel,order)
snglrty *panel;
int order;
{
  int i, j, nside,  N, M, N1, M1, M2, MN1, MN2, momOrder;
  double dx, dy, dxdy, dydx, x, y, z, SI, *xp, *yp, *xpn, *ypn;
  double ypp[3];
  static double *XP[4], *YP[4], **I;
  static double CS[16] = { 0.0, 1.0, 1.0, 1.5, 1.5, 3.75, 1.0, 3.0, 
			   1.5, 7.5, 1.5, 1.5, 3.75, 1.5, 7.5, 3.75 };

  /* An order of at least 4 is required. */
  momOrder = MAX(order,4);
  
  /* Allocate temporary storage and initialize arrays. */
  if(momOrder > maxorderMom) {
    for(i = 0; i < 4; i++) {
      CALLOC(XP[i], (momOrder+3), double, ON, AQ2P);
      CALLOC(YP[i], (momOrder+3), double, ON, AQ2P);
    }
    /* Allocate the euclidean moments matrix, Imn. */
    CALLOC(I, (momOrder+1), double*, ON, AQ2P);
    for(i=0; i <= momOrder; i++) CALLOC(I[i], (momOrder+1), double, ON, AQ2P);
    maxorderMom = momOrder;
  }

  /* First zero out the Moments matrix. */
  for(i = 0; i <= momOrder; i++) { 
    for(j = 0; j <= momOrder; j++) I[i][j] = 0.0; 
  }
    
  /* Compute powers of x and y at corner pts. */
  for(i = 0; i < panel->shape; i++) {
    xp = XP[i];
    yp = YP[i];
    xp[1] = panel->corner[i][XI];
    yp[1] = panel->corner[i][YI];
    for(j = 2; j <= momOrder+2; j++) {
      xp[j] = xp[j-1] * xp[1];
      yp[j] = yp[j-1] * yp[1];
    }
  }

  /* First moment, easy, just the panel area. */
  I[0][0] = panel->area;

  /* By using centroid, (1,0) and (0,1) are zero, so begin with (2,0). */
  for(nside = 0; nside < panel->shape; nside++) {
    xp = XP[nside];
    yp = YP[nside];
    if(nside == (panel->shape - 1)) {
      xpn = XP[0];
      ypn = YP[0];  
    }
    else {
      xpn = XP[nside + 1];
      ypn = YP[nside + 1];
    }

    dx = xpn[1] - xp[1];
    dy = ypn[1] - yp[1];

    if(fabs(dx) >= fabs(dy)) {
      dydx = dy/dx;
      for(M = 2; M <= momOrder; M++) {
	M1 = M + 1;
	M2 = M + 2;

	SI = ((xpn[M1] * ypn[1]) - (xp[M1] * yp[1])) / M1
	     + dydx * (xp[M2] - xpn[M2]) / (M1 * M2);
	I[M][0] += SI;

	for(N = 1; N <= M; N++) {
	  N1 = N + 1;
	  MN1 = M - N + 1;
	  SI = (xpn[MN1] * ypn[N1] - xp[MN1] * yp[N1]) / (MN1 * N1)
		     - (dydx * N * SI) / MN1;
	  I[M-N][N] += SI;
	}
      }
    }
    else {
      dxdy = dx/dy;
      for(M = 2; M <= momOrder; M++) {
	M1 = M + 1;
	M2 = M + 2;
	SI = (dxdy / (M1 * M2)) * (ypn[M2] - yp[M2]);
	I[0][M] += SI;
	for(N = 1; N <= M; N++) {
	  MN1 = M - N + 1;
	  MN2 = MN1 + 1;
          SI = dxdy * ((xpn[N] * ypn[MN2] - xp[N] * yp[MN2]) / (MN1 * MN2) 
			- (N * SI / MN1));
	  I[N][M-N] += SI;
	}
      }
    }
  }

  /* Now Create the S vector for calcp. */
  for(i = 0, M = 0; M <= 4; M++) {
    for(N = 0; N <= (4 - M); N++) {
      i++;
      panel->moments[i] = I[M][N] * CS[i];
    }
  }
  calcSphericalPanel(panel, I, order);
}


/*
  calcSphericalPoint computes the multipole coefficients for a point 
  source.  Note that the point source is unique among point singularities
  in that it has no orientation.  All the rest, dipole, vortex blob, etc.
  may need quite a bit of the code found in calcSphericalPanel becuase they 
  will require the rotation of their coefficients if they are defined in 
  local coordinates.
*/
calcSphericalPoint(pp, order)
snglrty *pp;
int order;
{
  int numterms;
  double *multiM;

  /* Allocate space for the multipole vector. */
  /* cftk Since pp->order exists, does that mean that individual sings can have
     their own order?  If so, here it is such that the number of terms is 
     unity. */
  pp->order = order;
  /* I allocate numterms elements, all but the first will be zero. Otherwise
     we will seg fault when we do matrix vector multiplies later.*/

  numterms = multerms(order);

  /*cftk Here you need to set up for whatever singularities you need.*/
  /* Allocate space for the source. */
  CALLOC(multiM, numterms, double, ON, AQ2P);
  if (pp->rhsType == POINT_SOURCE) {
    pp->multipoleR = multiM;
  }
  else {
    fprintf("FLE-calcpsph: unknown point singularity type: %d\n",pp->lhsType);
  }
  if (pp->lhsType == POINT_SOURCE) {
    pp->multipoleL = multiM;
  }
  else {
    fprintf("FLE-calcpsph: unknown point singularity type: %d\n",pp->lhsType);
  }
/* cftk fill in the point source definition here.*/
  multiM[index(0,0)] = 1.0;
}

/*
  calcSphericalPanel computes the multipole coefficients for a constant 
  source or dipole distribution on a planar panel.  First they are computed 
  in the panel coordinate system, then rotated into the global coord system. 
*/
calcSphericalPanel(pp, I, order)
double **I;
snglrty *pp;
int order;
{
  /* These are the temporary spherical harmonic coefficient matrices
     containing both the source terms. */
  static double **Mmn, **Mrmn, **Mtmn, **Mrtmn;
  /* These are the temporary spherical harmonic coefficient matrices
     containing both the dipole terms. */
  static double **Nmn, **Nrmn, **Ntmn, **Nrtmn;
  static double **Pmn, **Binom;
  static double *Fact;
  double *multiM, *multiD;
  double msumc, msums, sumc, sums, sign, sqfac;
  double cosa, sina, cosb, signb, cosg, sing;
  double dplus, dminus, coeff, rcoeff, icoeff;
  double **createBinom(), **createPmn();
  double *createFactorial();
  int i, m, mp, n, r, halfn, halfm, flrm, ceilm, numterms, rterms;
  double jacobid();

  if(order > maxorderSphere) {
    /* Allocate a temporary Multipole moments matrix, Mmn. */
    CALLOC(Mmn, (order+1), double*, ON, AQ2P);
    for(i = 0; i <= order; i++) CALLOC(Mmn[i], (order+1), double, ON, AQ2P);
    CALLOC(Mtmn, (order+1), double*, ON, AQ2P);
    for(i = 0; i <= order; i++) CALLOC(Mtmn[i], (order+1), double, ON, AQ2P);

    /* Allocate a rotated temporary Multipole matrix. */
    CALLOC(Mrmn, (order+1), double*, ON, AQ2P);
    for(i = 0; i <= order; i++) CALLOC(Mrmn[i], (order+1), double, ON, AQ2P);
    CALLOC(Mrtmn, (order+1), double*, ON, AQ2P);
    for(i = 0; i <= order; i++) CALLOC(Mrtmn[i], (order+1), double, ON, AQ2P);

    /* Allocate temporary spherical harmonic coefficient matrices for dipole 
       terms, N*mn. */
    CALLOC(Nmn, (order+1), double*, ON, AQ2P);
    for(i = 0; i <= order; i++) CALLOC(Nmn[i], (order+1), double, ON, AQ2P);
    CALLOC(Ntmn, (order+1), double*, ON, AQ2P);
    for(i = 0; i <= order; i++) CALLOC(Ntmn[i], (order+1), double, ON, AQ2P);

    CALLOC(Nrmn, (order+1), double*, ON, AQ2P);
    for(i = 0; i <= order; i++) CALLOC(Nrmn[i], (order+1), double, ON, AQ2P);
    CALLOC(Nrtmn, (order+1), double*, ON, AQ2P);
    for(i = 0; i <= order; i++) CALLOC(Nrtmn[i], (order+1), double, ON, AQ2P);

    Binom = createBinom((order+1) * 2);
    Pmn = createPmn((order+1) * 2);
    Fact = createFactorial((order+1) * 2 + 1);

    maxorderSphere = order;
  }

  /* Allocate space for the multipole vector. */
  pp->order = order;
  numterms = multerms(order);
  rterms = costerms(order);

  /*cftk Here you need to set up for whatever singularities you need.*/
  /* Allocate both source and dipole space */
  CALLOC(multiM, numterms, double, ON, AQ2P);
  CALLOC(multiD, numterms, double, ON, AQ2P);
  if (pp->rhsType == CONSTANT_SOURCE) {
    pp->multipoleR = multiM;
  }
  else if(pp->rhsType == CONSTANT_DIPOLE) {
    pp->multipoleR = multiD;
  }
  else {
    fprintf("FLE-calcpsph: unknown singularity type: %d\n",pp->lhsType);
  }
  if (pp->lhsType == CONSTANT_SOURCE) {
    pp->multipoleL = multiM;
  }
  else if(pp->lhsType == CONSTANT_DIPOLE) {
    pp->multipoleL = multiD;
  }
  else {
    fprintf("FLE-calcpsph: unknown singularity type: %d\n",pp->lhsType);
  }

  /* Create the M^0_n multipole coefficients first. */
  Nmn[0][0] = 0;
  for(n = 0; n <= order; n++) {
    halfn = n / 2;
    if(2 * halfn == n) { /* n even. */
      for(sumc = 0.0, i=0; i <= halfn; i++) {
	sumc += Binom[halfn][halfn - i] * I[n - 2 * i][2 * i];
      }
      Mmn[0][n] = sumc * Pmn[0][n];
      /* The dipole assignment is added. */
      if (n < order) Nmn[0][n+1] = (n+1) * sumc * Pmn[0][n];
    }
    else {
      Mmn[0][n] = 0.0;
      if (n < order) Nmn[0][n+1] = 0.0;
    }
  }

  /* Create the rest. */
  for(n = 1; n <= order; n++) {
    for(m = 1; m <= n; m++) {
      halfn = (n - m) / 2;
      if((2 * halfn) == (n - m)) { /* n-m even. */
	halfm = m / 2;
	if((2 * halfm) != m) { flrm = halfm; ceilm = halfm + 1;	}
	else { flrm = halfm; ceilm = halfm; }

	for(sumc = 0.0, sums = 0.0, i=0; i <= halfn; i++) {
	  sign = -1.0;
	  for(msumc = 0.0, r=0; r <= flrm; r++) {
	    sign *= -1.0;
	    msumc += sign * Binom[m][m - 2*r] * I[n - 2*(r+i)][2*(r+i)];
	  }

	  sign = 1.0;
	  for(msums = 0.0, r=0; r <= (ceilm-1); r++) {
	    sign *= -1.0;
	    msums += 
		sign * Binom[m][m - (2*r+1)] * I[n - (2*(r+i)+1)][2*(r+i) + 1];
	  }
	  sumc += Binom[halfn][halfn - i] * msumc;
	  sums += Binom[halfn][halfn - i] * msums;
	}
	/* The source assignment. */
        Mmn[m][n]  = sumc *  2.0 * Pmn[m][n];
	Mtmn[m][n] = sums * -2.0 * Pmn[m][n];
	/* The dipole assignment. */
	if (n < order) {
	  sqfac = sqrt( (double) (((n+1) * (n+1)) + (m*m)));
	  Nmn[m][n+1]  = sqfac * sumc *  2.0 * Pmn[m][n];
	  Ntmn[m][n+1] = sqfac * sums * -2.0 * Pmn[m][n];
	}
      }
    }
  }

  /* Compute the Euler angles between the panel and global coords. */
  eulerAngles(pp->X, pp->Y, pp->Z, &cosa, &sina, &cosb, &signb, &cosg, &sing);


  /* First rotation about cosa, sina (we are doing both source and dipole). */
  zrotate(cosa, sina, Mmn, Mtmn, order);
  zrotate(cosa, sina, Nmn, Ntmn, order);
  /* Rotate theta by 2pi if z prime is in the -xpp half space. */
  /* This operates on all coefficients, so N's are just included. */

  if(signb < 0) {
    cosb *= -1.0;
    for(n = 0; n <= order; n++) {
      if((n % 2) == 0) { rcoeff = 1.0; icoeff = -1.0; }
      else { rcoeff = -1.0; icoeff = 1.0; }
      
      for(m = 0; m <= n; m++) {
	Mmn[m][n] *= rcoeff;
	Mtmn[m][n] *= icoeff;
	Nmn[m][n] *= rcoeff;
	Ntmn[m][n] *= icoeff;
      }
    }
  }

  /* Now finish the rotation, assuming z prime is in + xpp half-space. */
  /* Prevously this rotation skiped the (n-mp) non-even terms as these
    were assumed to be zero. [Note that a rotation about z (above) does not
    introduce any terms not already non-zero.]  We now include all the terms,
    as the dipole has non zero terms where the source does not and 
    vice-versa.*/

  for(n = 0; n <= order; n++) {
    for(m = 0; m <= n; m++) {
      Mrmn[m][n] = 0.0;
      Mrtmn[m][n] = 0.0;
      Nrmn[m][n] = 0.0;
      Nrtmn[m][n] = 0.0;
      for(mp = 0; mp <= n; mp++) {
	coeff = Fact[n - mp] / Fact[n - m];
	dplus = jacobid(mp, m, n, cosb, Binom, Fact);
	dminus = jacobid(-mp, m, n, cosb, Binom, Fact);
	if(m == 0) { /* This test has to do with conjugation properties. */
	  icoeff = 0.0;
	  if((mp % 2) == 0) rcoeff = 0.5 * coeff * (dplus + dminus);
	  else rcoeff = 0.5 * coeff * (dplus - dminus);
	  Mrmn[m][n] += rcoeff * Mmn[mp][n];
	  Nrmn[m][n] += rcoeff * Nmn[mp][n];
	}
	else {
	  if((mp % 2) == 0) {
	    rcoeff = coeff * (dplus + dminus);
	    icoeff = coeff * (dplus - dminus);
	  }
	  else {
	    rcoeff = coeff * (dplus - dminus);
	    icoeff = coeff * (dplus + dminus);
	  }
	  Mrmn[m][n] += rcoeff * Mmn[mp][n];
	  Mrtmn[m][n] += icoeff * Mtmn[mp][n];
	  Nrmn[m][n] += rcoeff * Nmn[mp][n];
	  Nrtmn[m][n] += icoeff * Ntmn[mp][n];
	}
      }
    }
  }

  /* Final rotation about z by cosg, sing. (we are doing both source 
     and dipole). */

  zrotate(cosg, -sing, Mrmn, Mrtmn, order);
  zrotate(cosg, -sing, Nrmn, Nrtmn, order);

  /* Copy the Multipole matrix into the multi vector for this panel. */
  for(n = 0; n <= order; n++) {
    multiM[index(n,0)] = Mrmn[0][n];
    /*flip dipole sign here to agree with the panel normal convention.*/
    multiD[index(n,0)] = -Nrmn[0][n];
    for(m = 1; m <= n; m++) {
      multiM[index(n,m)] = Mrmn[m][n];
      multiM[sindex(n,m,rterms)] = Mrtmn[m][n];
      /*flip dipole sign here to agree with the panel normal convention.*/
      multiD[index(n,m)] = -Nrmn[m][n];
      multiD[sindex(n,m,rterms)] = -Nrtmn[m][n];
    }
  } 
}


/* Calculate the Euler Angles between X, Y, Z and global coordinate system. */
static eulerAngles(X, Y, Z, pcosa, psina, pcosb, psignb, pcosg, psing)
double X[3], Y[3], Z[3], *pcosa, *psina, *pcosb, *psignb, *pcosg, *psing;
{
  double ypp[3], xpp[3], norm;

/* First calculate ypp, the intersection of X-Y plane with global XY plane. */
  ypp[ZI] = 0.0;
  if(Y[ZI] != 0.0) {
    ypp[XI] = Y[ZI] * X[XI] - X[ZI] * Y[XI];
    ypp[YI] = Y[ZI] * X[YI] - X[ZI] * Y[YI];
  }
  else {
    ypp[XI] = Y[XI];
    ypp[YI] = Y[YI];
  }

/* Force ypp to be in the positive global Y plane. */
  if(ypp[YI] < 0) {
    ypp[YI] *= -1.0;
    ypp[XI] *= -1.0;
  }

/* Normalize the ypp vector. Use fact that ypp[ZI] = 0. */
  norm = sqrt(ypp[XI] * ypp[XI] + ypp[YI] * ypp[YI]);
  ypp[XI] /= norm;
  ypp[YI] /= norm;

  Cross_Product(ypp, Z, xpp);
  

/* Now calculate angles. */
  *pcosa = ypp[XI] * Y[XI] + ypp[YI] * Y[YI];
  *psina = ypp[XI] * X[XI] + ypp[YI] * X[YI];
  *pcosb = Z[ZI];
  if(xpp[ZI] <= 0) *psignb = 1.0;
  else *psignb = -1.0;
/*  *psignb = 1.0; for debugging. */
  *pcosg = ypp[YI];
  *psing = ypp[XI];
}

/* Rotate multipole coeffs about angle defined by cosa and sina. */
zrotate(cosa, sina, M, Mt, order)
double cosa, sina, **M, **Mt;
int order;
{
  double temp, cosm, sinm;
  int m, n;

  cosm = cosa;
  sinm = sina;
  for(m = 1; m <= order; m++) {
    for(n = m; n <= order; n++) {
      temp = M[m][n];
      M[m][n] = cosm * temp - sinm *  Mt[m][n];
      Mt[m][n] = cosm * Mt[m][n] + sinm * temp;
    }
    /* Compute  exp((m+1)a) = exp(a) * exp(ma). */
    temp = cosm;
    cosm = temp * cosa - sinm * sina;
    sinm = temp * sina + sinm * cosa;
  }
}

/*
More than just the jacobi polynomial, is actually a generalized 
spherical harmonic, which is kind of a funny name if you ask me.
For the nonstandard use of multipoles,

jacobid = 1/(2^m) * (n + |m|)!/(n + |mp|)!
              * (1-x)^(0.5 * (m-mp)) * (1+x)^(0.5 * (m+mp))
	         * jacobi(m-mp, m+mp, n-m, x);
*/

double jacobid(m, mp, n, x, binom, fact)
int n, m, mp;
double x, *fact, **binom;
{
  double ax, coeff = 1.0, val, part1, part2;
  int a, b, o, i;
  int temp, half, twosupm;

  /* First change m and mp to make sure m-mp > 0 and m+mp > 0. */
  if((m + mp) < 0) {
    m *= -1;
    mp *= -1;
    if(((mp - m) % 2) != 0) coeff *= -1.0;  
  }
  if((m - mp) < 0) {  /* Must swap m and mp. */
    temp = m;
    m = mp;
    mp = temp;
    coeff = (fact[n-m] * fact[n+m])
	    / (fact[n-mp] * fact[n+mp]);
    if(((m - mp) % 2) != 0) coeff *= -1.0;
  }

  /* Compute 2^m. */
  for(twosupm = 1, i=1; i <= m; i++) twosupm *= 2;

  if((m - mp) == 0) part1 = 1.0;
  else part1 = pow(1.0 - x, 0.5 * (m - mp));
  if((m + mp) == 0) part2 = 1.0;
  else part2 = pow(1.0 + x, 0.5 * (m + mp));

  val = coeff * part1 * part2 / ((double) twosupm);

  /* The jacobi polynomial. */
  a = m - mp;
  b = m + mp;
  o = n - m;
  for(ax = 1.0, i=o; i > 0; i--) {
    ax = 1 - ax * (1.0 - x) * ((double) (o - i + 1) * (a + b + i + o)) 
	     / ((double) (2 * i * (a + i)));
  }
  val *= binom[o + a][o] * ax;

/*
  printf("x=%g m=%d mp=%d, n=%d a=%d b=%d val=%g\n", x, m, mp, n, a, b, val);
*/
  return(val);
}

/* Create Binomial Array. */
double **createBinom(order)
int order;
{
  int m, n;
  double **Binom;

  CALLOC(Binom, (order+1), double *, ON, AQ2P);  
  for(m=0; m <= order; m++) {
    CALLOC(Binom[m], (order+1), double, ON, AQ2P);
  }

  for(m=0; m <= order; m++) {
    Binom[m][0] = 1.0;
    for(n = 1; n <= order; n++) {
      Binom[m][n] = ((m - (n-1))  * Binom[m][n-1]) / n;
    }
  }
  return(Binom);
}

/* Create Factorial Array. */
double *createFactorial(size)
int size;
{
  double *Fact;
  int i;

  CALLOC(Fact, size + 1, double, ON, AQ2P);  

  Fact[0] = 1.0;
  for(i=1; i <= size; i++) {
    Fact[i] = i * Fact[i-1];
  }
  return(Fact);
}


/* Create Legendre(0) array. */
double **createPmn(order)
int order;
{
  int n, m;
  double pmm;
  double **Pmn;

  CALLOC(Pmn, (order+1), double *, ON, AQ2P);  
  for(n=0; n <= order; n++) {
    CALLOC(Pmn[n], (order+2), double, ON, AQ2P);
  }

  pmm = 1.0;

  for(m=0; m <= order; m++) {
    Pmn[m][m] = pmm;
    Pmn[m][m+1] = 0.0;
    for(n = (m+2); n <= order; n++) {
      Pmn[m][n] = -((double) (m + n - 1)) / ((double) (n - m)) * Pmn[m][n-2];
    }
    pmm *= -(2.0 * m + 1.0);
  }
  return(Pmn);
}

/*
calcp computes the potential at x, y, z due to various singularities. 
At present the POINT_SOURCE, CONSTANT_SOURCE, and CONSTANT_DIPOLE are 
implemented.  Note that this program structure, where calcp is called
for any singularity type, is probably not the most efficient for any 
particular application.  For production codes, one is likely to want to 
tailor calcp to a particular set of singularity types and get rid of 
any others.  One might also think to move the if up a level and call
different routines for each singularity.  Furthermore, initcalcp is
calculating a whole bunch of panel attributes which are needed for 
the particular distribution algorithms, so different algorithms may not
need this call to procede the call to calcp.  A trivial example is 
the desingularized source method which needs no panel info at all
(although that is not quite the way I do it here).  In any case, for 
the most straightforward path to computing coefficients for new singularity 
types in FastLap, add the computation here with an if-block on the type.


A pointer to the coefficients for the type of singularity on the RHS 
is the return, a pointer to the coefficient for the type of singularity 
on the LHS is the last item in the arg list.

Note, the code for the constant strength distributions is subtle because 
there are 5 cases depending on the placement of the collocation point:

    CASE1: evaluation point projection strictly outside of panel (happens most
      of the time).
    CASE2: eval pnt proj. strictly inside panel (happens >= #panels times,
      at least for the self terms)
    CASE3: eval pnt proj. on a panel side, not on a corner (happens
      when paneled faces meet at right angles, also possible other ways)
    CASE4: eval pnt proj. on a panel corner (happens rarely to never)
    CASE5: eval pnt proj. on side extension (happens when paneled 
      faces meet at right angles, also possible other ways).

jt: modified this routine to include directional derivatives of the
    source/dipole potential. For the full gradient call this routine
    thrice (i.e. once for each direction normal = e_1, ... ). This is
    certainly not efficient, but more compatible with the rest of fastlap.

*/
double calcp(sing, x, y, z, deriv, normal, pfd)  /* jt: new variables */
snglrty *sing;
int deriv;                           /* jt: deriv=0/1 for potential/derivative */
double x, y, z, *normal, *pfd;       /* jt: normal: direction of derivative */
{
  double r[4], fe[4], xmxv[4], ymyv[4];
  double xc, yc, zc, zsq, xn, yn, zn, znabs, xsq, ysq, rsq, diagsq, dtol;
  double v, arg, st, ct, length, s1, c1, s2, c2, s12, c12, val;
  double *s;
  double rInv, r2Inv, r3Inv, r5Inv, r7Inv, r9Inv, zr2Inv;
  double ss1, ss3, ss5, ss7, ss9;
  double s914, s813, s411, s512, s1215;
  double fs, fd, fdsum;
  int okay, i, next;
  struct edge *edge;
  double *corner;
/* jt: variables for gradient evaluations */
  double nrmx, nrmy, nrmz;
  double rss3, ssx3, ssy3, ssx5, ssy5, rss7, ssx7, ssy7;
  double txy, rss9, ssx9, ssy9, rss11, ssx11, ssy11;  
  double fsx, fsy, fsz, fdx, fdy, fdz;
  double xri[4], yri[4], fln, fac, fh1, fh2, u1, u2, rr, zz;
  double nDrvSrc, nDrvDip;


  /* c2.0 Branch here for desingularized source formulation and return before
     doing any un-needed panel computations. */
  if(sing->lhsType == POINT_SOURCE || sing->rhsType == POINT_SOURCE) {
    xc = x - sing->x;
    yc = y - sing->y;
    zc = z - sing->z;
    /* cwarn: You might want to trap for r = 0. here. */
    fs = 1.0 / sqrt(xc*xc + yc*yc + zc*zc);
    if ( deriv == 1 ) {  /* jt: return gradient */
      rss3 = fs*fs*fs;
      fsx = (sing->x - x)*rss3; 
      fsy = (sing->y - y)*rss3; 
      fsz = (sing->z - z)*rss3;
      fs = fsx*normal[0] + fsy*normal[1] + fsz*normal[2];
    }
    if(pfd != NULL) *pfd = fs;
    return (fs);
  }

  /* Put the evaluation point into this panel's coordinates. */

  xc = x - sing->x;
  yc = y - sing->y;
  zc = z - sing->z;

  xn = DotP_Product(sing->X, xc, yc, zc);
  yn = DotP_Product(sing->Y, xc, yc, zc);
  zn = DotP_Product(sing->Z, xc, yc, zc);

  if ( deriv == 1 ) {   /* jt: put normal in panel's coordinates */
    nrmx = DotP_Product(sing->X, normal[0], normal[1], normal[2]);
    nrmy = DotP_Product(sing->Y, normal[0], normal[1], normal[2]);
    nrmz = DotP_Product(sing->Z, normal[0], normal[1], normal[2]);
  }

  zsq = zn * zn;
  xsq = xn * xn;
  ysq = yn * yn;
  rsq = zsq + xsq + ysq;

#if JTDEBUG
calcflag = 0;
#endif

  diagsq = sing->max_diag * sing->max_diag;
  dtol = EQUIV_TOL * sing->min_diag; /* jt: moved this statement out of the if block */
  znabs = fabs(zn);                  /* jt: moved this statement out of the if block */

  if(rsq > (LIMITFOURTH * diagsq)) { 

#if JTDEBUG
calcflag += 1;
#endif
    fs = 0.0; fd = 0.0;
    s = sing->moments;
    /* First, second moments. */
    r2Inv = 1.0 / rsq;
    rInv = sqrt(r2Inv);
    r3Inv = r2Inv * rInv;
    r5Inv = r3Inv * r2Inv;
    zr2Inv = zn * r2Inv;
    ss1 = s[1] * rInv;
    ss3 = -(s[3] + s[10]) * r3Inv;
    ss5 = (xsq * s[10] + (xn * yn * s[7]) + ysq * s[3]) * r5Inv;
    fs = ss1 + ONE3 * ss3 + ss5;
    fdsum = ss1 + ss3 + 5.0 * ss5;
    fd = zr2Inv * fdsum;
    if ( deriv == 1 ) {
      rss3 = r2Inv*ss1;
      ssx3 = -xn*rss3;
      ssy3 = -yn*rss3;
      ssx5 = (xn*(s[3]+3.0*s[10])+yn*s[7])*r5Inv;
      ssy5 = (yn*(s[10]+3.0*s[3])+xn*s[7])*r5Inv;
      rss7 = -5.0*r2Inv*ss5;
      ssx7 = xn*rss7;
      ssy7 = yn*rss7;
      fsx = ssx3 + ssx5 + ssx7;
      fsy = ssy3 + ssy5 + ssy7;
      fdx = zr2Inv*(3.0*ssx3+5.0*ssx5+7.0*ssx7);
      fdy = zr2Inv*(3.0*ssy3+5.0*ssy5+7.0*ssy7);
      fdz = r2Inv*fdsum - zr2Inv*zr2Inv*(3.0*ss1 + 5.0*ss3 + 35.0*ss5);
    }


    if(rsq < (LIMITSECOND * diagsq)) {

#if JTDEBUG
calcflag += 1;
#endif
    /* Third and fourth moments added for diagsq/r2 between 40 and 150. */
      s914 = s[9] + s[14];
      s813 = s[8] + s[13];
      s411 = s[4] + s[11];
      s512 = s[5] + s[12];
      s1215 = s[12] + s[15];
      r7Inv = r5Inv * r2Inv;
      r9Inv = r7Inv * r2Inv;
      ss5 = (-xn * s813 - yn * s411 + 0.1 * (s512 + s1215)) * r5Inv;

      ss7 = (FIVE3 *((xn * xsq * s[13] + yn * ysq * s[4]) 
		     + 3.0 * xn * yn * (xn * s[11]  +  yn * s[8]))
		     - xsq * s1215 - ysq * s512 - xn * yn * s914) * r7Inv;

      ss9 = (7.0 * (ONE6 * (xsq * xsq * s[15] + ysq * ysq * s[5])
		    + xsq * ysq * s[12])
	     + SEVEN3 * xn * yn * (xsq * s[14] + ysq * s[9])) * r9Inv;

      fs += ss5 + ss7 + ss9;
      fdsum = 5.0 * ss5 + 7.0 * ss7 + 9.0 * ss9;
      fd += zr2Inv * fdsum;
      if ( deriv == 1 ) {
      txy = 2*xn*yn;
      ssx5 = -s813*r5Inv;
      ssy5 = -s411*r5Inv;
      rss7 =  5.0*r2Inv*ss5;

      ssx7 = (5.0*(xn*xn*s[13] + txy*s[11] + yn*yn*s[8]) - s1215*(xn+xn) -
                       yn*s914)*r7Inv - xn*rss7;
      ssy7 = (5.0*(yn*yn*s[4] + xn*xn*s[11] + txy*s[8]) - s512*(yn+yn) -
                       xn*s914)*r7Inv - yn*rss7;
      rss9 = 7.0*ss7*r2Inv;
      ssx9 = (FT3*xn*xsq*s[15] + 14.0*xn*ysq*s[12] + 49.0*yn*(xsq*s[14] +
                       ONE3*ysq*s[9]))*r9Inv - xn*rss9;
      ssy9 = (FT3*yn*ysq*s[5] + 14.0*yn*xsq*s[12] + 49.0*xn*(ysq*s[9] +
                       ONE3*xsq*s[14]))*r9Inv - yn*rss9;
      rss11 = 9.0*ss9*r2Inv;
      ssx11 = -xn*rss11;
      ssy11 = -yn*rss11;

      fsx += ssx5+ssx7+ssx9+ssx11;
      fsy += ssy5+ssy7+ssy9+ssy11;
      fdx += zr2Inv*(5.0*ssx5 + 7.0*ssx7 + 9.0*ssx9 + 11.0*ssx11);
      fdy += zr2Inv*(5.0*ssy5 + 7.0*ssy7 + 9.0*ssy9 + 11.0*ssy11);
      fdz += r2Inv*fdsum - zr2Inv*zr2Inv*(35.0*ss5 + 63.0*ss7 + 99.0*ss9);
    }

      num4th++;
    }
    else num2nd++;
  }

  else {

    if ( deriv == 1 ) {
      fsx = fsy = 0.0;
      fdx = fdy = fdz = 0.0;
    }

    /* Always move the evaluation point a little bit off the panel. */
    if(znabs < dtol) { 
      zn = 0.5 * dtol;  /* Half of dtol insures detection for zero dipole. */
      znabs = 0.5 * dtol;
    }

    /* Once per corner computations. */
    for(okay = TRUE, i=0; i < sing->shape; i++) {
      corner = sing->corner[i]; 
      xmxv[i] = xc = xn - corner[XI];
      ymyv[i] = yc = yn - corner[YI];
      zc = zn - corner[ZI];
      fe[i] = xc * xc + zc * zc;
      r[i] = sqrt(yc * yc + fe[i]);
      if(r[i] < (1.005 * znabs)) {  /* If r almost z, on vertex normal. */
	okay = FALSE;
      }
      if ( deriv == 1 ) {
        xri[i] = xmxv[i]/r[i];
        yri[i] = ymyv[i]/r[i];
      }
    }

    /* Once per edge computations. */
    fs = 0.0; fd = 0.0;
    fsx = fsy = fdx = fdy = fdz = 0.0;   /* jt: initialize gradients */

    for(i=0; i < sing->shape; i++) {
      if(i == (sing->shape - 1)) next = 0;
      else next = i + 1;

      /* Now calculate the edge contributions to a panel. */
      length = sing->length[i];
      ct = (sing->corner[next][XI] - sing->corner[i][XI]) / length;
      st = (sing->corner[next][YI] - sing->corner[i][YI]) / length;

      /* v is projection of eval-i edge onto perpend to next-i edge. */
      /* Exploits the fact that corner points in panel coordinates. */
      v = xmxv[i] * st - ymyv[i] * ct;

      /* arg == zero if eval on next-i edge, but then v = 0. */
      arg = (r[i] + r[next] - length)/(r[i] + r[next] + length);
      fln = -log(arg);              /* jt: need this later */
      if(arg > 0.0) fs += v*fln;    /* jt: was: fs -= v*log(arg) */
      if ( deriv == 1 )
        if ( arg > 0.0 ) {
          fac = (r[i] + r[next] - length)*(r[i] + r[next] + length);
          fac = v*(length+length)/fac;
          fsx += fln*st - fac*(xri[i] + xri[next]);
          fsy -= fln*ct + fac*(yri[i] + yri[next]);
          fdz -= fac*( 1.0/r[i] + 1.0/r[next] );
        }

      /* Okay means eval not near a vertex normal, Use Hess-Smith. */
      if(okay) {
	s1 = v * r[i];
	c1 = znabs * (xmxv[i] * ct + ymyv[i] * st);
	s2 = v * r[next];
	c2 = znabs * (xmxv[next] * ct + ymyv[next] * st);
      }
      /* Near a vertex normal, use Newman. */
      else {
	s1 = (fe[i] * st) - (xmxv[i] * ymyv[i] * ct);
	c1 = znabs * r[i] * ct;
	s2 = (fe[next] * st) - (xmxv[next] * ymyv[next] * ct);
	c2 = znabs * r[next] * ct;
      }    

      s12 = (s1 * c2) - (s2 * c1);
      c12 = (c1 * c2) + (s1 * s2);
      val = atan2(s12, c12);
      fd += val;
      if ( deriv == 1 ) {
        u1   = xmxv[i]*ct + ymyv[i]*st;
        u2   = xmxv[next]*ct+ymyv[next]*st;
        if (!okay) {
          rr   = r[i]*r[i];
          fh1  = xmxv[i]*ymyv[i];
          fh2  = xmxv[next]*ymyv[next];
          fac  = c1/((c1*c1+s1*s1)*rr );
          fdx += (rr*v+fh1*u1)*fac;
          fdy -= fe[i]*u1*fac;
          rr   = r[next]*r[next];
          fac  = c2/((c2*c2+s2*s2)*rr);
          fdx -= (rr*v+fh2*u2)*fac;
          fdy += fe[next]*u2*fac;
        }
        else {
          fac  = zn/(c1*c1+s1*s1);
          fdx += (u1*v*xri[i]+r[i]*ymyv[i])*fac;
          fdy += (u1*v*yri[i]-r[i]*xmxv[i])*fac;
          fac  = zn/(c2*c2+s2*s2);
          fdx -= (u2*v*xri[next]+r[next]*ymyv[next])*fac;
          fdy -= (u2*v*yri[next]-r[next]*xmxv[next])*fac;
        }
      }
    }

    /* Adjust the computed values. */

    if(fd < 0.0) fd += TWOPI;
    if(zn < 0.0) fd *= -1.0;

    /* If in the same plane as panel, fd = 0. */
    if(znabs < dtol) {
      if(rsq < (dtol * dtol)) fd = -TWOPI; /* This is temporary, see sign flip below. */
      else fd = 0.0;
    }
    fs -= zn * fd;
    if ( deriv == 1 ) {
      fsx -= zn*fdx;
      fsy -= zn*fdy;
    }
    numexact++;
  }
/*  jt:
 *  fsx, fsy, fsz and fdx, fdy, fdz  are the partials of the source and the
 *  dipole potential in the panel coordinate system. If the field point is in the
 *  panel plane set fsz to zero; if it is in the panel itself, set fsz to +/- 2 PI.
 *  The sign depends on whether we approach from the interior (+) or the exterior
 *  (-) of the solid. 
 */
  if ( deriv == 1 ) {
    if(rsq < (dtol * dtol)) nDrvSrc = TWOPI;
    else nDrvSrc = nrmx*fsx + nrmy*fsy - nrmz*fd;
    nDrvDip = nrmx*fdx + nrmy*fdy + nrmz*fdz;
  }

  /* Return values of the source and dipole dist's.*/
  fd *= -1.0;  /* Because panel normal points out of computational domain. */ 

  if(fs < 0.0) {
    fprintf(stderr, "FLE-clacp: fs (exact source) is less than zero = %g\n", fs);
    fprintf(stderr, "FLE-calcp: Okay = %d Evaluation Point = %g %g %g\n", okay, x, y, z);
    fprintf(stderr, "FLE-calcp: Evaluation Point in local coords = %g %g %g\n", xn, yn, zn);
    fprintf(stderr, "FLE-calcp: Panel Description Follows\n");
    dp(sing);
    exit(0);
  }
  /*Put the correct coefficient in place for the LHS.*/
  if( sing->lhsType == CONSTANT_SOURCE ) {
  if( pfd != NULL )
    if ( deriv == 1 ) /* jt: normal derivative or potential */
      *pfd = nDrvSrc;
    else
      *pfd = fs;
  }
  else if( sing->lhsType == CONSTANT_DIPOLE ) {
  if( pfd != NULL ) 
    if( deriv == 1 ) /* jt: normal derivative or potential */
      *pfd = nDrvDip;
    else
      *pfd = fd;
  }
  else {
    fprintf("FLE-calcp: unknown singularity type: %d\n",sing->lhsType);
  }
  /*Put the correct coefficient in place for the RHS.*/
  if( sing->rhsType == CONSTANT_SOURCE ) {
    if ( deriv == 1 ) /* jt: normal derivative or potential */
      return (nDrvSrc);
    else
      return (fs);
  }
  else if( sing->rhsType == CONSTANT_DIPOLE ) {
    if( deriv == 1 ) /* jt: normal derivative or potential */
      return(nDrvDip);
    else
      return (fd);
  }
  else {
    fprintf("FLE-calcp: unknown singularity type: %d\n",sing->rhsType);
  }
}


/* Debugging Print Routines follow. */

dumpnums(flag, size)
int flag, size;
{
  double total;

  if(flag == ON) {		/* if first call */
    num2ndsav = num2nd;
    num4thsav = num4th;
    numexactsav = numexact;
  }
  else {
    total = num2ndsav + num4thsav + numexactsav;
#if MULDAT == ON
    fprintf(stdout, "Potential coefficient counts\n multipole only:\n");
    fprintf(stdout, 
	    "  2nd order: %d %.3g%%; 4th: %d %.3g%%; Integral: %d %.3g%%\n",
	    num2nd, 100*(num2ndsav/total), num4th, 100*(num4thsav/total), 
	    numexact, 100*(numexactsav/total));
#endif
    total = num2nd + num4th + numexact;
#if MULDAT == ON
    fprintf(stdout, " multipole plus adaptive:\n");
    fprintf(stdout, 
	    "  2nd order: %d %.3g%%; 4th: %d %.3g%%; Integral: %d %.3g%%\n",
	    num2nd, 100*(num2nd/total), num4th, 100*(num4th/total), 
	    numexact, 100*(numexact/total));
    fprintf(stdout, "Percentage of multiplies done by multipole: %.3g%%\n",
	    100*(size*size - total)/(size*size));
#endif
  }
}

dp(panel)
snglrty *panel;
{
int i, numterms;

  printf("shape=%d maxdiag=%g mindiag=%g area=%g\n", panel->shape, 
	 panel->max_diag, panel->min_diag, panel->area);

  printf("x=%g y=%g z=%g\n", panel->x, panel->y, panel->z);
  printf("X= %g %g %g\n", panel->X[0], panel->X[1], panel->X[2]);
  printf("Y= %g %g %g\n", panel->Y[0], panel->Y[1], panel->Y[2]);
  printf("Z= %g %g %g\n", panel->Z[0], panel->Z[1], panel->Z[2]);

  for(i=0; i < panel->shape; i++)
      printf("corner%d = %g %g %g\n", 
	     i, panel->corner[i][0], panel->corner[i][1], panel->corner[i][2]);

  for(i=0; i < panel->shape; i++)
      printf("length%d = %g\n", i, panel->length[i]);

  printf("moment coeffs:  ");
  for(i=0; i < 16; i++) {
    printf("%g  ", panel->moments[i]);
    if( (i % 6) == 0) printf("\n");
  }

  printf("\nmono multipole coeffs:  ");
  numterms = multerms(panel->order);
  for(i=0; i < numterms; i++) {
    printf("%g  ", panel->multipoleL[i]);
    if( (i % 6) == 0) printf("\n");
  }

  printf("\ndipole multipole coeffs:  ");
  numterms = multerms(panel->order);
  for(i=0; i < numterms; i++) {
    printf("%g  ", panel->multipoleR[i]);
    if( (i % 6) == 0) printf("\n");
  }

  printf("\n");
}



#define DIS 2
#define SCALE 2.0

testCalcp(pp)
snglrty *pp;
{

  double offx, offy, offz, x, y, z, mult;
  int i, j, k;

  offx = pp->x;
  offy = pp->y;
  offz = pp->z;

  mult = 0.5 * pp->max_diag;

  dp(pp);

  printf("\n\nCenter Point %g %g %g\n", offx, offy, offz);
  for(i=0; i < DIS; i++) {
    for(j=0; j < DIS; j++) {
      for(k=0; k < DIS; k++) {
	x = offx + i * mult * SCALE;
	y = offy + j * mult * SCALE;
	z = offz + k * mult * SCALE;
	printf("\nEval pt = %g %g %g\n", x, y, z);
        calcp(pp, x, y, z, NULL);
      }
    }
  }
}

