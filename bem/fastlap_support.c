/*
# -*- coding: utf-8 -*-
#
#   bem: triangulation and fmm/bem electrostatics tools 
#   fastlap support code originally from fastlap, (C) fastlap authors
*/

#include <math.h>
#include "mulStruct.h"
#include "mulGlobal.h"
#include "fastlap_support.h"

#define XI 0
#define YI 1
#define ZI 2
#define ONE3 0.3333333333333
#define Dot_Product(V1,V2) V1[XI]*V2[XI]+V1[YI]*V2[YI]+V1[ZI]*V2[ZI]
#define DotP_Product(V1,R,S,T) (V1[XI])*(R)+(V1[YI])*(S)+(V1[ZI])*(T)

void Dcentroid(int shape, double *pc, double *xcout)
{
  double corner[4][3], X[3], Y[3], Z[3], vertex1[3], vertex3[3];
  double sum, delta, x1, y1, x2, x3, y3, xc, yc;
  int i, j;
  double normalize();
  /* Load the corners. */
  for(i=0; i<4; i++) { 
      for(j=0; j<3; j++) { 
      corner[i][j] = *(pc++);
      }
  }

  /* Use vertex 0 as the origin and get diags and lengths. */
  for(sum=0, i=0; i<3; i++) {
    X[i] = delta = corner[2][i] - corner[0][i];
    sum += delta * delta;
    vertex1[i] = corner[1][i] - corner[0][i];
    if(shape == QUADRILAT) {
      vertex3[i] = corner[3][i] - corner[0][i];
      Y[i] = corner[1][i] - corner[3][i];
    }
    else if(shape == TRIANGLE) {
      vertex3[i] = corner[2][i] - corner[0][i];
      Y[i] = corner[1][i] - corner[0][i];
    }
    else {
      printf("Dcentroid FE: Shape indicator is neither triangle nor quadrilateral");
      exit(0);
    }
  }
  x2 = sqrt(sum);

  /* Z-axis is normal to two diags. */
  Cross_Product(X, Y, Z);
  normalize(X);
  normalize(Z);

  /* Real Y-axis is normal to X and Z. */
  Cross_Product(Z, X, Y);

  /* Project into the panel axes. */
  y1 = Dot_Product(vertex1, Y);
  y3 = Dot_Product(vertex3, Y);
  x1 = Dot_Product(vertex1, X);
  x3 = Dot_Product(vertex3, X);

  yc = ONE3 * (y1 + y3);
  xc = ONE3 * (x2 + ((x1 * y1 - x3 * y3)/(y1 - y3)));

  *(xcout+0) = corner[0][XI] + xc * X[XI] + yc * Y[XI];
  *(xcout+1) = corner[0][YI] + xc * X[YI] + yc * Y[YI];
  *(xcout+2) = corner[0][ZI] + xc * X[ZI] + yc * Y[ZI];
}
