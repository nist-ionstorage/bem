/*
 *
 * extNeum.c:  Exterior Neumann Problem with fastlap
 *  sample program
 *          
 *
 *  usage:  extNeum  [-j=%d -m=%d -t=%d -r=%f -s]  panelfile
 *      where:
 *        -j        job:  0: density,  1: gradient,  2: both
 *        -m        multipole expansion order
 *        -t        tree depth
 *        -s        switch orientation of the panel normal
 *        -r        radius for points where the gradient is calculated
 *        panelfile file with the panel description
 *
 *
 *  solves the Exterior Neumann Problem for the geometry given by the
 *  panelfile.
 *
 *  Written by Johannes Tausch
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mulStruct.h"
#include "mulGlobal.h"


/* fastlap routines used by the driver */
int fastlap();
double normalize();
void Cross_Product();

#define MAX(A,B)  ( (A) > (B) ? (A) : (B) )
#define size 4000
#define ZERO 0
#define POINT_SOURCE 1
#define CONSTANT_SOURCE 11
#define CONSTANT_DIPOLE 12
#define LINEAR_SOURCE 21
#define LINEAR_DIPOLE 22
#define DIMEN 3

int read_panels( panelfile, shapep, xp, typep, lhsvectp, rhsvectp, sOrient ) 
char *panelfile;
double **xp, **lhsvectp, **rhsvectp;
int **shapep,  **typep, sOrient;
/*
 * Read in panel data from a file and returns the number of panels read.
 * If sOrient is nonzero then the panel vertices are read in reversed order.
 */
{
double *xq, *x, *lhsvect, *rhsvect;
int *shape, *type;
char cshape;
char buf[512];
int nPanels, i;
FILE *fp, *fopen();

if ( ( fp = fopen( panelfile, "r") ) == NULL ) {
  printf("\ncould not open %s\n",panelfile);
  exit(1);
}
fgets(buf, 511, fp);
printf("The header in this file is\n%s\n", buf);
i = 0;
while ( fgets(buf, 511, fp) != NULL )  i++;
nPanels = i;

if ( (x = (double*)malloc( nPanels*12*sizeof(double) ) ) == NULL ) exit(0);
if ( (lhsvect = (double*)malloc( nPanels*sizeof(double) ) ) == NULL ) exit(0);
if ( (rhsvect = (double*)malloc( nPanels*sizeof(double) ) ) == NULL ) exit(0);
if ( (shape = (int*)malloc( nPanels*sizeof(int) ) ) == NULL ) exit(0);
if ( (type  = (int*)malloc( nPanels*sizeof(int) ) ) == NULL ) exit(0);

rewind(fp);
fgets(buf, 511, fp);
for ( i=0; i<nPanels; i++ ) {
  fgets(buf, 511, fp);
  xq = x+12*i;
  if ( (buf[0] == 'q') || (buf[0] == 'Q') || (buf[0] == '4') ) {     /* quadrilat */
    if ( sOrient )
    sscanf(buf,"%c %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
           &cshape, xq+9, xq+10, xq+11, xq+6, xq+7, xq+8, xq+3, xq+4,
           xq+5, xq, xq+1, xq+2, rhsvect+i, lhsvect+i, type+i);

    else
    sscanf(buf,"%c %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
           &cshape, xq, xq+1, xq+2, xq+3, xq+4, xq+5, xq+6, xq+7, xq+8,
           xq+9, xq+10, xq+11, rhsvect+i, lhsvect+i, type+i);
    shape[i] = 4;
  }
  else {
    if ( (buf[0] == 't') || (buf[0] == 'T') || (buf[0] == '3') ) {   /* triangle */
      if ( sOrient )
      sscanf(buf,"%c %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
             &cshape, xq+6, xq+7, xq+8, xq+3, xq+4, xq+5, xq, xq+1, xq+2,
             rhsvect+i, lhsvect+i, type+i);
      else
      sscanf(buf,"%c %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
             &cshape, xq, xq+1, xq+2, xq+3, xq+4, xq+5, xq+6, xq+7, xq+8, 
             rhsvect+i, lhsvect+i, type+i);
      shape[i] = 3;
    }
    else {
      printf("\nERROR reading panels: Panel shape %c is not supported\n", buf[0]);
      exit(1);
    }
  }
}
printf("%d panels read\n",nPanels);
fclose(fp);
*xp = x;
*lhsvectp = lhsvect;
*rhsvectp = rhsvect;
*shapep = shape;
*typep = type;
return nPanels;
}




test_orientation( nPanels, x, xcoll )
int nPanels;
double *x, *xcoll;
{
double a1[3], a2[3], nrm[3];
double *xp, inner;
int i,k;

for ( i=0; i<nPanels; i++ ) {
  xp = x+12*i;
  for (k=0; k<3; k++ ) {
    a1[k] = xp[3+k] - xp[k];
    a2[k] = xp[6+k] - xp[k];
  }
  Cross_Product(a1, a2, nrm);
  normalize(nrm);
  inner = nrm[0]*xcoll[3*i] + nrm[1]*xcoll[3*i+1] + nrm[2]*xcoll[3*i+2];
/*  if ( inner < 0 ) {
    printf("\nError: panel %d has incorrect orientation inner = %lf\n",i,inner);
    exit(1);
  }*/
}

}


void get_collocPoints( nsing, x, dtype, shape, xcoll, xnrm, area )
int nsing, *dtype, *shape; 
double *x, *xcoll, *xnrm, *area;
/*
 * The field points are the collocation points (= centroids of panels).
 * The directions for the derivatives are the panel normals.
 * Mimics the calculation of the panel centroids in calcp().
 * This is important, because otherwise calcp() may not recognize a panel
 * with the field point in it.
 */
{
snglrty pq;
double *xp, delta, sum, sum2, length, maxlength, minlength, length20, length31;
int next,i,j,k;

for ( i=0; i<nsing; i++ ) {
  xp = x+12*i;
  pq.shape = shape[i];
  dtype[i] = 1;
  for(j=0;j<shape[i];j++) {
    for(k=0;k<3;k++) {
      (pq.corner[j])[k] = xp[(j*DIMEN)+k];
    }
  }
  /* Calculate edge lengths. */
  maxlength = 0.0;
  minlength = 1.0e20;
  for(k=0; k<shape[i]; k++) {    
    if(k == (shape[i] -1)) next = 0;
    else next = k + 1;
    for(sum= 0, j = 0; j < 3; j++) {
      delta = pq.corner[next][j] - pq.corner[k][j];
      sum += delta * delta;
    }
    pq.length[k] = length = sqrt(sum);
    maxlength = MAX(maxlength, length);
    minlength = MIN(minlength, length);
  }

  /* Get diags and lengths. */
  for(sum= 0.0, sum2 = 0.0, j = 0; j < 3; j++) {     
    pq.X[j] = delta = pq.corner[2][j] - pq.corner[0][j];
    sum += delta * delta;
    if(shape[i] == 3) pq.Y[j] = pq.corner[1][j] - pq.corner[0][j];      
    else {
      pq.Y[j] = delta = pq.corner[1][j] - pq.corner[3][j];      
      sum2 += delta * delta;
    }
  }
  length20 = sqrt(sum);
  length31 = sqrt(sum2);
      
  /* Check on lengths for quad. */
  if(shape[i] == 3) {
    pq.max_diag = maxlength;
    pq.min_diag = minlength;
  }
  else {
    length = MAX(length20, length31);
    pq.max_diag = length;
    pq.min_diag = MIN(length20, length31);
  }
      
  /* Z-axis is normal to two diags. */
  Cross_Product(pq.X, pq.Y, pq.Z);
  area[i] = pq.area = 0.5 * normalize(pq.Z);
  normalize(pq.X);
      
  /* Real Y-axis is normal to X and Z. */
  Cross_Product(pq.Z, pq.X, pq.Y);
      
  /* Project the corner points into the plane defined by edge midpoints. */
  planarize(&pq);    

  /* Calculate the centroid. */
  centroid(&pq, length20);

  xcoll[3*i]   = pq.x;
  xcoll[3*i+1] = pq.y;
  xcoll[3*i+2] = pq.z;

  /* the local Z-axis points into the domain, the normal points in 
   * the opposite direction 
   */
  xnrm[3*i]   = -pq.Z[0];
  xnrm[3*i+1] = -pq.Z[1];
  xnrm[3*i+2] = -pq.Z[2];
}

}/* get_collocPoints */






void main(int nargs, char *argv[]) 
{
double *x, *xcoll, *xnrm, *lhsvect, *rhsvect, *area;
double l8_error[4], l2_error[4], fourPi, radius=1.0, one=1.0, rad3, dummy;
int *type, *shape, *dtype, *rhstype, *lhstype, *rhsindex, *lhsindex;
char panelfile[80];
int nPanels, nummom, numlev, pos, i, j, k, nlhs, nrhs, nsing;
int fljob, job=2, switchOrientation=0;

double tolpar = 0.000001, tol;
int numitr, maxitr=50;

FILE *fpout0, *fpout1, *fpout2, *fpout3, *fopen();


tol = tolpar;

/* parse the command line */
panelfile[0] = 0;
nummom = numlev = -1;
for ( i=1; i<nargs; i++ )
  if ( argv[i][0] == '-' )
    switch ( argv[i][1] ) {
    case 'j': job = atoi( argv[i]+3 );
      break;
    case 't': numlev = atoi( argv[i]+3 );
      break;
    case 'm': nummom = atoi( argv[i]+3 );
      break;
    case 's': switchOrientation = 1;
      break;
    case 'r': radius = atof( argv[i]+3 );
    }
  else
    strcpy(panelfile,argv[i]);

if ( panelfile[0] == 0 ) {
  printf("\nName of the panel file > ");
  scanf("%s",panelfile);
}
if ( nummom < 0 ) {
  printf("Select expansion order > ");
  scanf("%d",&nummom);
}
if ( numlev < 0 ) {
  printf("Select tree depth > ");
  scanf("%d",&numlev);
}
if ( panelfile[0] == 0 ) {
  printf("Filename for input data > ");
  scanf("%s", panelfile);
}
rad3 = radius*radius*radius;
fourPi = 4*M_PI;
/*
 * get panels and collocation points
 */
nPanels = read_panels( panelfile, &shape, &x, &type, &lhsvect, &rhsvect, 
                       switchOrientation );

if ( (dtype  = (int*)malloc( nPanels*sizeof(int) ) ) == NULL ) exit(0);
if ( (xcoll = (double*)malloc( 3*nPanels*sizeof(double) ) ) == NULL ) exit(0);
if ( (xnrm  = (double*)malloc( 3*nPanels*sizeof(double) ) ) == NULL ) exit(0);
if ( (area  = (double*)malloc( nPanels*sizeof(double) ) ) == NULL ) exit(0);
if ( (lhstype  = (int*)malloc( nPanels*sizeof(int) ) ) == NULL ) exit(0);
if ( (rhstype  = (int*)malloc( nPanels*sizeof(int) ) ) == NULL ) exit(0);
if ( (lhsindex = (int*)malloc( nPanels*sizeof(int) ) ) == NULL ) exit(0);
if ( (rhsindex = (int*)malloc( nPanels*sizeof(int) ) ) == NULL ) exit(0);


nsing = nPanels;
get_collocPoints( nsing, x, dtype, shape, xcoll, xnrm, area );
/* test_orientation( nPanels, x, xcoll );*/

/*
 *  Set up the parameters for fastlap to solve the indirect integral formulation. 
 *  Here we ignore some of the information in the panelfile.
 */
if ( job != 1 ) {
  for ( i=0; i<nPanels; i++ ) {
    rhstype[i] = CONSTANT_SOURCE;
    lhstype[i] = CONSTANT_SOURCE;
    rhsvect[i] = -fourPi;
    rhsindex[i] = i;
    lhsindex[i] = i;
    dtype[i] = 1;
  }

  fljob = INDIRECT;
  nrhs = nsing;
  nlhs = nsing;
  tol = tolpar;
  numitr = fastlap(&nlhs,&nrhs,&nsing,x,shape,dtype,lhstype,rhstype,lhsindex,
                   rhsindex,lhsvect,rhsvect,xcoll,xnrm,&numlev,&nummom,&maxitr,
                   &tol,&fljob);

  printf("%d iterations knocked down the residual to: %lf",numitr, tol);

/*
 * Calculate the error for the sphere with boundary condition f=1
 * Does NOT work for other boundary conditions/geometries
 */
  fpout0 = fopen("density", "w"); 
  printf("Have opened density for density vector.\n");

  l2_error[0] = l8_error[0] = 0.0;
  for ( i=0; i<nlhs; i++ ) {
    dummy = fabs(lhsvect[i] - 1.0);
    fprintf(fpout0, "%d %lf %lf %lf\n",i,one,lhsvect[i],dummy);
    l2_error[0] += dummy*dummy*area[i];
    if ( dummy > l8_error[0] ) l8_error[0] = dummy;
  }
  (void)fflush(fpout0);
  (void)fclose(fpout0);
  l2_error[0] = sqrt( l2_error[0] );

}

/*
 * Calculate the gradient of the solution; call fastlap once for each component.
 * Assesment of the solution works only for the sphere.
 */
if ( job > 0 ) {
  fljob = FIELD;

  for ( i=0; i<nPanels; i++ ) {
    if ( job == 1 )                /* the known stuff is on the r.h.s */
      rhsvect[i] = 1.0;
    else
      rhsvect[i]  = lhsvect[i];
    xcoll[3*i]   *= radius;
    xcoll[3*i+1] *= radius;
    xcoll[3*i+2] *= radius;
  }

  /* Partial w.r.t. to x */
  for ( i=0, j=0; i<nPanels; i++, j+=3 ) {
    rhstype[i] = CONSTANT_SOURCE;
    lhstype[i] = CONSTANT_SOURCE;
    rhsindex[i] = i;
    lhsindex[i] = i;
    xnrm[j]   = 1.0;
    xnrm[j+1] = 0.0;
    xnrm[j+2] = 0.0;
    dtype[i] = 1;
  }

  nrhs = nsing;
  nlhs = nsing;
  tol = tolpar;
  numitr = fastlap(&nlhs,&nrhs,&nsing,x,shape,dtype,lhstype,rhstype,lhsindex,
                   rhsindex,lhsvect,rhsvect,xcoll,xnrm,&numlev,&nummom,&maxitr,
                   &tol,&fljob);

/*  fpout1 = fopen("d-dx", "w"); 
  printf("Have opened d-dx for d-dx vector.\n");*/

  l2_error[1] = l8_error[1] = 0.0;
  for ( i=0, j=0; i<nlhs; i++, j+=3 ) {
    dummy = fabs(lhsvect[i] + fourPi*xcoll[j]/rad3);
/*    fprintf(fpout1, "%d %lf %lf %lf\n",i,fourPi*xcoll[j]/rad3,lhsvect[i],dummy);*/
    l2_error[1] += dummy*dummy*area[i];
    if ( dummy > l8_error[1] ) l8_error[1] = dummy;
  }
/*  (void)fclose(fpout1);*/
  l2_error[1] = sqrt( l2_error[1] );

  /* Partial w.r.t. to y */
  for ( i=0, j=0; i<nPanels; i++, j+=3 ) {
    rhstype[i] = CONSTANT_SOURCE;
    lhstype[i] = CONSTANT_SOURCE;
    rhsindex[i] = i;
    lhsindex[i] = i;
    xnrm[j]   = 0.0;
    xnrm[j+1] = 1.0;
    xnrm[j+2] = 0.0;
    dtype[i] = 1;
  }

  nrhs = nsing;
  nlhs = nsing;
  tol = tolpar;
  numitr = fastlap(&nlhs,&nrhs,&nsing,x,shape,dtype,lhstype,rhstype,lhsindex,
                   rhsindex,lhsvect,rhsvect,xcoll,xnrm,&numlev,&nummom,&maxitr,
                   &tol,&fljob);

/*  fpout2 = fopen("d-dy", "w"); 
  printf("Have opened d-dy for d/dy vector.\n");*/

  l2_error[2] = l8_error[2] = 0.0;
  for ( i=0, j=0; i<nlhs; i++, j+=3 ) {
    dummy = fabs(lhsvect[i] + fourPi*xcoll[j+1]/rad3);
/*    fprintf(fpout2, "%d %lf %lf %lf\n",i,fourPi*xcoll[j+1]/rad3,lhsvect[i],dummy);*/
    l2_error[2] += dummy*dummy*area[i];
    if ( dummy > l8_error[2] ) l8_error[2] = dummy;
  }
/*  (void)fclose(fpout2);*/
  l2_error[2] = sqrt( l2_error[2] );

  /* Partial w.r.t. to z */
  for ( i=0, j=0; i<nPanels; i++, j+=3 ) {
    rhstype[i] = CONSTANT_SOURCE;
    lhstype[i] = CONSTANT_SOURCE;
    rhsindex[i] = i;
    lhsindex[i] = i;
    xnrm[j]   = 0.0;
    xnrm[j+1] = 0.0;
    xnrm[j+2] = 1.0;
    dtype[i] = 1;
  }

  nrhs = nsing;
  nlhs = nsing;
  tol = tolpar;
  numitr = fastlap(&nlhs,&nrhs,&nsing,x,shape,dtype,lhstype,rhstype,lhsindex,
                   rhsindex,lhsvect,rhsvect,xcoll,xnrm,&numlev,&nummom,&maxitr,
                   &tol,&fljob);

/*  fpout3 = fopen("d-dz", "w"); 
  printf("Have opened d-dz for d/dz vector.\n");*/

  l2_error[3] = l8_error[3] = 0.0;
  for ( i=0, j=0; i<nlhs; i++, j+=3 ) {
    dummy = fabs(lhsvect[i] + fourPi*xcoll[j+2]/rad3);
/*    fprintf(fpout3, "%d %lf %lf %lf\n",fourPi*xcoll[j+2]/rad3,lhsvect[i],dummy);*/
    l2_error[3] += dummy*dummy*area[i];
    if ( dummy > l8_error[3] ) l8_error[3] = dummy;
  }
/*  (void)fclose(fpout3);*/
  l2_error[3] = sqrt( l2_error[3] );

}

/*
 * Print errors
 */
if ( job != 1 )
  printf("\nDensity  : l2-error = %lf l8-error = %lf"  , l2_error[0], l8_error[0] );
if ( job > 0 ) {
  printf("\nd/dx Psi : l2-error = %lf l8-error = %lf", l2_error[1], l8_error[1] );
  printf("\nd/dy Psi : l2-error = %lf l8-error = %lf", l2_error[2], l8_error[2] );
  printf("\nd/dz Psi : l2-error = %lf l8-error = %lf", l2_error[3], l8_error[3] );
}
printf("\n");

}




