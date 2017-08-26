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

/* # ***** sort to /src/io
   # ***** */
#include <stdio.h>
#include <math.h>
#include "mulStruct.h"
#include "mulGlobal.h"

disExtrasimpcube(pc)
cube *pc;
{
  printf("cubes[%d][%d][%d][%d]\n", pc->level, pc->j, pc->k, pc->l);
}

disExParsimpcube(pc)
cube *pc;
{
  cube *pa = pc->parent;
  printf("cubes[%d][%d][%d][%d], ", pc->level, pc->j, pc->k, pc->l);
  printf("parent = cubes[%d][%d][%d][%d]\n", pa->level, pa->j, pa->k, pa->l);
}

dissimpcube(pc)
cube *pc;
{
int i;
  printf("cube center: x=%g y=%g z=%g\n", pc->x, pc->y, pc->z);
  printf("index=%d level=%d exactUp=%d exactDown=%d numkids=%d\n",
	 pc->index,pc->level,pc->exactUp,pc->exactDown,pc->numkids);
  printf("numnbrs=%d upnumvects=%d directnumvects=%d downnumvects=%d\n",
	 pc->numnbrs, pc->upnumvects, pc->directnumvects, pc->downnumvects);
}

discube(pc)
cube *pc;
{
int i;
  printf("cube center: x=%g y=%g z=%g\n", pc->x, pc->y, pc->z);
  printf("index=%d level=%d exactUp=%d exactDown=%d numkids=%d\n",
	 pc->index,pc->level,pc->exactUp,pc->exactDown,pc->numkids);
  printf("numnbrs=%d upnumvects=%d directnumvects=%d downnumvects=%d\n",
	 pc->numnbrs, pc->upnumvects, pc->directnumvects, pc->downnumvects);
  if(pc->directnumvects > 0) {
    printf("num of elements in ");
    for(i=0; i < pc->directnumvects; i++) {
      printf("v%d = %d ", i, pc->directnumeles[i]);
    }
    printf("\nsngs\n");
    for(i=0; i < pc->directnumeles[0]; i++) {
      dischg(pc->sngs[i]);
    }
  }
  if(pc->downnumvects > 0) {
    printf("num of down elements in ");
    for(i=0; i < pc->downnumvects; i++) {
      printf("v%d = %d ", i, pc->downnumeles[i]);
    }
  }
}

disupcube(pc)
cube *pc;
{


}

disdirectcube(pc)
cube *pc;
{
int i;
  for(i=0; i < pc->directnumvects; i++) {
    printf("matrix %d\n", i);
    dismat(pc->directmats[i], pc->directnumeles[0], pc->directnumeles[i]);
  }
}


dissys(sys)
ssystem *sys;
{
int i, j, k, l, side;
  printf("side=%d depth=%d order=%d\n",
	 sys->side, sys->depth, sys->order);
  printf("Cube corner is x=%g y=%g z=%g\n", sys->minx, sys->miny, sys->minz);
  printf("Cube side length= %g\n", sys->length);
  printf("Printing all the cubes\n");
  for(i = 0, side = 1; i <= sys->depth; i++, side *= 2) {
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
	  fprintf(stdout, "\ncubes[%d][%d][%d][%d]\n", i, j, k, l);
	  dissimpcube(&(sys->cubes[i][j][k][l]));
/*	  disdirectcube(&(sys->cubes[i][j][k][l])); */
	}
      }
    }
  }
}



dismat(mat, rows, cols)
double **mat;
int rows, cols;
{
int i,j;
  if(cols != 0) {
    for(i=0; i < rows; i++) {
      printf("\n i=%d\n", i);
      for(j=0; j < cols; j++) {
        printf("%d %g  ", j, mat[i][j]);
        if(((j+1) % 5) == 0) printf("\n");
      }
    }
    printf("\n");
  }
}


disvect(v, size)
double *v;
int size;
{
int i;
  for(i=0; i < size; i++) {
    printf("i=%d %g ", i, v[i]);
    if(((i+1) % 5) == 0) printf("\n");
  }
  printf("\n");
}

dischg(pq)
snglrty *pq;
{
  printf("cond=%d index=%d\n", pq->cond, pq->index);
}

disallchg(pq) 
snglrty *pq;
{
snglrty *nq;
  for(nq = pq; nq != NULL; nq = nq->next) disfchg(pq);
}

disfchg(pq) 
snglrty *pq;
{
/*
  printf("Cond=%d Corners\n", pq->cond);
  printf("x0=%g y0=%g z0=%g\n", pq->x0, pq->y0, pq->z0);
  printf("x1=%g y1=%g z1=%g\n", pq->x1, pq->y1, pq->z1);
  printf("x2=%g y2=%g z2=%g\n", pq->x2, pq->y2, pq->z2);
  printf("x3=%g y3=%g z3=%g\n", pq->x3, pq->y3, pq->z3);
  printf("Center\n");
  printf("x=%g y=%g z=%g\n", pq->x, pq->y, pq->z);
*/
}

/*
dumpMat dumps a rows x cols matrix of doubles; assumes indices from zero. 
*/
void dumpMat(mat, rows, cols)
int rows, cols;
double **mat;
{
  int i, j;
  for(i = 0; i < rows; i++) {
    fprintf(stdout, "    row%d ", i);
    for(j = 0; j < cols; j++) {
      if(mat[i][j] < 0.0) fprintf(stdout, "%.5e ", mat[i][j]);
      else fprintf(stdout, " %.5e ", mat[i][j]);
    }
    fprintf(stdout, "\n");
  }
}

/*
dumpChgs dumps the relative coordinates of an array of snglrtys or evaluation 
points.
*/
void dumpChgs(sngs, numsngs, x, y, z)
int numsngs;
double x, y, z;
snglrty **sngs;
{
  int i;
  double rho, cosA, beta;
  for(i = 0; i < numsngs; i++) {
    xyz2sphere(sngs[i]->x, sngs[i]->y, sngs[i]->z,
	       x, y, z, &rho, &cosA, &beta);
    fprintf(stdout, "    %d %d ", i, sngs[i]->cond);
    if(rho < 0) fprintf(stdout, "(%.5e ", rho);
    else fprintf(stdout, "( %.5e ", rho);
    if(cosA < 0) fprintf(stdout, "%.5e ", cosA);
    else fprintf(stdout, " %.5e ", cosA);
    if(beta < 0) fprintf(stdout, "%.5e) ", beta);
    else fprintf(stdout, " %.5e) ", beta);
    if(x < 0) fprintf(stdout, "(%.5e ", sngs[i]->x);
    else fprintf(stdout, "( %.5e ", sngs[i]->x);
    if(y < 0) fprintf(stdout, "%.5e ", sngs[i]->y);
    else fprintf(stdout, " %.5e ", sngs[i]->y);
    if(z < 0) fprintf(stdout, "%.5e)\n", sngs[i]->z);
    else fprintf(stdout, " %.5e)\n", sngs[i]->z);
  }
}

/*
  display the matrix built for a given snglrty to multipole transformation
*/
void dispQ2M(mat, sngs, numsngs, x, y, z, order)
int numsngs, order;
double **mat, x, y, z;
snglrty **sngs;
{
  fprintf(stdout, "\nQ2M MATRIX: cube at (%.5e %.5e %.5e)\n", x, y, z);
  dumpMat(mat, multerms(order), numsngs);
  fprintf(stdout, 
	  "    CHARGES IN CUBE # cond (rho_i cos(alpha_i) beta_i) (x y z):\n");
  dumpChgs(sngs, numsngs, x, y, z);
}

/*
  display the matrix built for a given multipole to local transformation
*/
void dispM2L(mat, x, y, z, xp, yp, zp, order)
int order;
double **mat, x, y, z, xp, yp, zp;
{
  fprintf(stdout, 
   "\nM2L MATRIX: multi at (%.5e %.5e %.5e) -> local at (%.5e %.5e %.5e)\n",
	  x, y, z, xp, yp, zp);
  dumpMat(mat, multerms(order), multerms(order));
}

/*
  display the matrix built for a given snglrty to local transformation
*/
void dispQ2L(mat, sngs, numsngs, x, y, z, order)
int numsngs, order;
double **mat, x, y, z;
snglrty **sngs;
{
  fprintf(stdout, "\nQ2L MATRIX: cube at (%.5e %.5e %.5e)\n", x, y, z);
  dumpMat(mat, multerms(order), numsngs);
  fprintf(stdout, 
	  "    CHARGES IN CUBE # cond (rho_i cos(alpha_i) beta_i) (x y z):\n");
  dumpChgs(sngs, numsngs, x, y, z);
}

/*
  display the matrix built for a given multipole to multipole transformation
*/
void dispM2M(mat, x, y, z, xp, yp, zp, order)
int order;
double **mat, x, y, z, xp, yp, zp;
{
  fprintf(stdout, 
      "\nM2M MATRIX: cube at (%.5e %.5e %.5e) shifted to (%.5e %.5e %.5e)\n", 
	  x, y, z, xp, yp, zp);
  dumpMat(mat, multerms(order), multerms(order));
}

/*
  display the matrix built for a given local to local transformation
*/
void dispL2L(mat, x, y, z, xp, yp, zp, order)
int order;
double **mat, x, y, z, xp, yp, zp;
{
  fprintf(stdout, 
      "\nL2L MATRIX: cube at (%.5e %.5e %.5e) shifted to (%.5e %.5e %.5e)\n", 
	  x, y, z, xp, yp, zp);
  dumpMat(mat, multerms(order), multerms(order));
}

/*
  display the matrix built for a given multipole to potential transformation
*/
void dispM2P(mat, x, y, z, sngs, numsngs, order)
int numsngs, order;
double **mat, x, y, z;
snglrty **sngs;
{
  fprintf(stdout, "\nM2P MATRIX: cube at (%.5e %.5e %.5e)\n", x, y, z);
  dumpMat(mat, numsngs, multerms(order));
  fprintf(stdout, 
	  "    EVAL PNTS IN CUBE # cond (rho_i, cos(alpha_i), beta_i):\n");
  dumpChgs(sngs, numsngs, x, y, z);
}

/*
  display the matrix built for a given local to potential transformation
*/
void dispL2P(mat, x, y, z, sngs, numsngs, order)
int numsngs, order;
double **mat, x, y, z;
snglrty **sngs;
{
  fprintf(stdout, "\nL2P MATRIX: cube at (%.5e %.5e %.5e)\n", x, y, z);
  dumpMat(mat, numsngs, multerms(order));
  fprintf(stdout, 
	  "    EVAL PNTS IN CUBE # cond (rho_i, cos(alpha_i), beta_i):\n");
  dumpChgs(sngs, numsngs, x, y, z);
}

/*
  displays upward pass and moment vectors associated with a cube - debug only
*/
void dumpUpVecs(pc)
cube *pc;
{
  int i, j;
  fprintf(stdout, 
    "\nUPWARD PASS/MOMENT VECTORS, LEVEL %d CUBE AT (%.5e %.5e %.5e):\n",
	  pc->level, pc->x, pc->y, pc->z);
  for(i = 0; i < pc->upnumvects; i++) {
    fprintf(stdout, "%d", i);
    for(j = 0; j < pc->upnumeles[i]; j++) {
      if(pc->upvects[i][j] < 0.0) 
	  fprintf(stdout, " %.5e", pc->upvects[i][j]);
      else fprintf(stdout, "  %.5e", pc->upvects[i][j]);
    }
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "M");
  for(j = 0; j < pc->multisize; j++) {
    if(pc->multi[j] < 0.0) fprintf(stdout, " %.5e", pc->multi[j]);
    else fprintf(stdout, "  %.5e", pc->multi[j]);
  }
  fprintf(stdout, "\n");
}

/*
  displays the upward pass vectors for the eight level 1 cubes - debug only
*/
void dumpLevOneUpVecs(sys)
ssystem *sys;
{
  int i, j, k;
  cube *****cubes = sys->cubes;
  for(i = 0; i < 2; i++) {
    for(j = 0; j < 2; j++) {
      for(k = 0; k < 2; k++) {
	if(cubes[1][i][j][k] != NULL) dumpUpVecs(cubes[1][i][j][k]);
      }
    }
  }
}

/*
  checks a cube (direct, eval) list for bad cube structs - debug only
*/
void chkList(sys, listtype)
ssystem *sys;
int listtype;			/* DIRECT or EVAL */
{
  int cnt[BUFSIZ];		/* # of cubes processed by level */
  int depth = sys->depth;
  int lev, nn;
  int i, j, k;
  cube *nc;
  for(i = 0; i <= depth; i++) cnt[i] = 0;
  if(listtype == DIRECT) {
    fprintf(stderr, "\nChecking Direct list.");
    nc = sys->directlist;
  }
  else if(listtype == EVAL) {
    fprintf(stderr, "\nChecking Evaluation list.\n");
    nc = sys->evallist;
  }
  else {
    fprintf("FW-chkList: Bogus list type indicator.\n");
  }
  while(nc != NULL) {
    /* check number and level of neighbors */
    lev = nc->level;
    nn = nc->numnbrs;
    for(i = 0; i < nn; i++) {
      if(lev != ((nc->nbrs)[i])->level) {
	fprintf(stderr, "chkList: level %d cube has a level %d nbr\n", lev,
		((nc->nbrs)[i])->level);
	fprintf(stderr, " ok cubes ");
	for(j = 0; j <= depth; j++) fprintf(stderr, "lev%d: %d ", j, cnt[j]);
	fprintf(stderr, "\n");
	exit(0);
      }
    }
    /* check number of kids */
    if(lev == depth && nc->numkids != 0) {
      fprintf(stderr, "chkList: level %d cube has children\n", lev);
      fprintf(stderr, " ok cubes ");
      for(j = 0; j <= depth; j++) fprintf(stderr, "lev%d: %d ", j, cnt[j]);
      fprintf(stderr, "\n");
      exit(0);
    }
    /* if lowest level, check status of eval or direct vects */
    if(lev == depth) {
      if(listtype == DIRECT) {
	if(nc->directnumeles == NULL) {
	  fprintf(stderr, "chkList: level %d cube has bad direct info\n", lev);
	  fprintf(stderr, " ok cubes ");
	  for(j = 0; j <= depth; j++) fprintf(stderr, "lev%d: %d ", j, cnt[j]);
	  fprintf(stderr, "\n");
	  exit(0);
	}
      }
      if(listtype == EVAL) {
	if(nc->evalnumvects == 0) {
	  fprintf(stderr, "chkList: level %d cube has no eval info\n", lev);
	  fprintf(stderr, " ok cubes ");
	  for(j = 0; j <= depth; j++) fprintf(stderr, "lev%d: %d ", j, cnt[j]);
	  fprintf(stderr, "\n");
	  exit(0);
	}
      }
    }
    cnt[lev]++;
    if(listtype == DIRECT) nc = nc->dnext;
    else if(listtype == EVAL) nc = nc->enext;
    else {
      fprintf(stderr, "chkList: bad flag\n");
      exit(0);
    }
  }
  if(listtype == DIRECT) fprintf(stdout, "\nDirect ");
  else if(listtype == EVAL) fprintf(stdout, "\nEvaluation ");
  fprintf(stdout, "list ok: ");
  for(j = 0; j <= depth; j++) fprintf(stdout, "lev%d: %d ", j, cnt[j]);
  fprintf(stdout, "\n\n");
}

/*
  dumpList dumps info about the cubes in a list.
*/
void dumpList(sys, listtype)
ssystem *sys;
int listtype;			/* DIRECT, EVAL, MULTIL, or LOCAL */
{
  int cnt[BUFSIZ];		/* # of cubes processed by level */
  int depth;
  int lev, nn;
  int i, j, k;
  cube *nc;

  if(listtype == LOCAL) {
    fprintf(stderr, "\n*******Checking Local list.*******\n");
    for(depth=2; depth <= sys->depth; depth++) {
      for(nc=sys->locallist[depth]; nc != NULL; nc = nc->lnext) {
	fprintf(stderr, "ll: ind %d %d %d %d, dwn %d, exDwn %d, fpts %d, lsz %d, isz %d\n", nc->level, nc->j, nc->k, nc->l, nc->downcube, nc->exactDown, nc->numfieldpts, nc->localsize, nc->interSize);
      }
    }
  }
  else if(listtype == MULTIL) {
    fprintf(stderr, "\n*******Checking Multi list.*******\n");
    for(depth=2; depth <= sys->depth; depth++) {
      for(nc=sys->multilist[depth]; nc != NULL; nc = nc->mnext) {
	fprintf(stderr, "ml: ind %d %d %d %d, up %d, exUp %d, sngs %d, msz %d\n", nc->level, nc->j, nc->k, nc->l, nc->upcube, nc->exactUp, nc->upnumeles[0], nc->multisize);
      }
    }
  }
  else  if(listtype == DIRECT) {
    fprintf(stderr, "\n*******Checking Direct list.*******\n");
    for(nc=sys->directlist; nc != NULL; nc = nc->dnext) {
      fprintf(stderr,"dl: ind %d %d %d %d, up %d, dwn %d, exU %d, exD %d, fpt %d, lsz %d, nbrs %d\n", nc->level, nc->j, nc->k, nc->l, nc->upcube, nc->downcube, nc->exactUp, nc->exactDown, nc->numfieldpts, nc->localsize, nc->numnbrs);
    }
  }
  else if(listtype == EVAL) {
    fprintf(stderr, "\n*******Checking Evaluation list.*******\n");
    for(nc=sys->evallist; nc != NULL; nc = nc->enext) {
      fprintf(stderr,"el: ind %d %d %d %d, up %d, dwn %d, exU %d, exD %d, fpt %d, sng %d, evalnv %d, evalnels[0] %d,lsz %d, isz %d\nx %g, y %g, z%g\n", nc->level, nc->j, nc->k, nc->l, nc->upcube, nc->downcube, nc->exactUp, nc->exactDown, nc->numfieldpts, nc->upnumeles[0], nc->evalnumvects, nc->evalnumeles[0], nc->localsize, nc->interSize, nc->x, nc->y, nc->z);
      for (i=0; i < nc->numfieldpts; i++) {
	fprintf(stderr, "i %d, x %g, y %g, z %g,\n", i, nc->fpts[i]->x, nc->fpts[i]->y, nc->fpts[i]->z);
      }
    }
  }
  else {
    fprintf("FW-chkList: Bogus list type indicator.\n");
  }
}

/*
  chks a cube for bad cube struct (direct, local or eval) entries - debug only
*/
void chkCube(sys, nc, listtype)
ssystem *sys;
cube *nc;
int listtype;			/* DIRECT, LOCAL, or EVAL */
{
  int depth = sys->depth;
  int lev, nn;
  int i, j, k;
  if(nc != NULL) {
    /* check number and level of neighbors */
    lev = nc->level;
    nn = nc->numnbrs;
    for(i = 0; i < nn; i++) {
      if(lev != ((nc->nbrs)[i])->level) {
	fprintf(stdout, "chkCube: level %d cube has a level %d nbr\n", lev,
		((nc->nbrs)[i])->level);
/*	exit(0);*/
      }
    }
    /* check number of kids */
    if(lev == depth && nc->numkids != 0) {
      fprintf(stdout, "chkCube: level %d cube has children\n", lev);
/*      exit(0);*/
    }
    /* if lowest level, check status of eval and direct vects */
    if((lev == depth) && (listtype == DIRECT)) {
      if(nc->directnumeles == NULL) {
	fprintf(stdout, 
		"chkCube: level %d cube has NULL directnumeles\n", lev);
/*	exit(0);*/
      }
    }
    if((lev == depth) && (listtype == EVAL)) {
      if(nc->evalnumvects == 0) {
	fprintf(stdout, "chkCube: level %d cube has no eval info\n", lev);
/*	exit(0);*/
      }
      if(nc->eval == NULL) {
	fprintf(stdout, "chkCube: level %d cube has no eval pntr\n", lev);
      }
    }
  }
}

/*
  checks the lowest level cubes for trouble using chkCube - debug only
*/
void chkLowLev(sys, listtype)
ssystem *sys;
int listtype;			/* DIRECT, LOCAL or EVAL */
{
  int i, j, k, l, side, depth = sys->depth, cnt = 0;
  cube *nc, *****cubes = sys->cubes;
  for(i = 1, side = 1; i <= depth; i++, side *= 2);
  for(j=0; j < side; j++) {	/* loop through all cubes at level depth */
    for(k=0; k < side; k++) {
      for(l=0; l < side; l++) {
	nc = cubes[depth][j][k][l];
	if(nc != NULL) {
	  chkCube(sys, nc, listtype);
	  cnt++;
	}
      }
    }
  }
  fprintf(stdout,"Total lowest level (level %d) cubes checked = %d\n", 
	  depth, cnt);
}
  

/*
  core display routine used below
*/
void dumpSynCore1(str, depth, fcnt, exupcnt, exdowncnt, emcnt, tcnt)
int depth, *fcnt, *exupcnt, *exdowncnt, *emcnt, *tcnt;
char *str;
{
  int i;
  fprintf(stdout, "%-13s", str);
  for(i = 0; i <= depth; i++) {
    sprintf(str, "%d/%d/%d/%d/%d ", fcnt[i], exupcnt[i], exdowncnt[i], emcnt[i], tcnt[i]);
    if(i < 2) fprintf(stdout, "%8s", str);
    else if(i == 2) fprintf(stdout, "%12s", str);
    else if(i == 3) fprintf(stdout, "%16s", str);
    else if(i == 4) fprintf(stdout, "%20s", str);
    else if(i == 5) fprintf(stdout, "%24s", str);
    else fprintf(stdout, "%28s", str);
  }
  fprintf(stdout, "\n");
}
/*
  core display rtn used below
*/
dumpSynCore2(str, depth, cnt)
int depth, *cnt;
char *str;
{
  int i;

  fprintf(stdout, "%-13s", str);
  for(i = 0; i <= depth; i++) {
    sprintf(str, "%d    ", cnt[i]);
    if(i < 2) fprintf(stdout, "%8s", str);
    else if(i == 2) fprintf(stdout, "%12s", str);
    else if(i == 3) fprintf(stdout, "%16s", str);
    else if(i == 4) fprintf(stdout, "%20s", str);
    else if(i == 5) fprintf(stdout, "%24s", str);
    else fprintf(stdout, "%28s", str);
  }
  fprintf(stdout, "\n");
}

/*
  displays number of exactUp, exactDown, full, empty and total cubes per 
  level in all cubes, and eval, direct, multi and local lists
*/
void dumpSynop(sys)
ssystem *sys;
{
  int i, j, k, l, side, depth = sys->depth, lev;
  int exupcnt[BUFSIZ], exdowncnt[BUFSIZ], fcnt[BUFSIZ], emcnt[BUFSIZ], tcnt[BUFSIZ];
  extern int *multicnt, *localcnt;
  char str[BUFSIZ];
  cube *****cubes = sys->cubes, *nc;

  for(i = 0; i <= depth; i++) exupcnt[i] = exdowncnt[i] = fcnt[i] = emcnt[i] = tcnt[i] = 0;

  fprintf(stdout, 
	  "\nCUBE AND EXPANSION SYNOPSIS (full/exactUP/exactDown/empty/ttl):\n");
  fprintf(stdout, "             ");
  for(i = 0; i <= depth; i++) {
    sprintf(str, "level %d ", i);
    if(i < 2) fprintf(stdout, "%8s", str);
    else if(i == 2) fprintf(stdout, "%12s", str);
    else if(i == 3) fprintf(stdout, "%16s", str);
    else if(i == 4) fprintf(stdout, "%20s", str);
    else if(i == 5) fprintf(stdout, "%24s", str);
    else fprintf(stdout, "%28s", str);
  }
  fprintf(stdout, "\n");
  /* dump cube usage by level */
  for(i = 0, side = 1; i <= depth; i++, side *= 2) {
    for(j=0; j < side; j++) {	/* loop through all cubes at levels >= 0 */
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
	  nc = cubes[i][j][k][l];
	  tcnt[i]++;
	  if(nc != NULL) {
	    lev = nc->level;
	    fcnt[i]++;
	    if(nc->exactUp == TRUE) exupcnt[i]++;
	    if(nc->exactDown == TRUE) exdowncnt[i]++;
	  }
	  else emcnt[i]++;
	}
      }
    }
  }
  sprintf(str, "All cubes");
  dumpSynCore1(str, depth, fcnt, exupcnt, exdowncnt, emcnt, tcnt);
  
  for(i = 0; i <= depth; i++) exupcnt[i] = exdowncnt[i] = fcnt[i] = emcnt[i] = tcnt[i] = 0;
  /* dump cube direct list by level */
  for(nc = sys->directlist; nc != NULL; nc = nc->dnext) {
    lev = nc->level;
    tcnt[lev]++;
    if(nc->upnumvects > 0) fcnt[lev]++;
    else emcnt[lev]++;
    if(nc->exactUp == TRUE) exupcnt[lev]++;
    if(nc->exactDown == TRUE) exdowncnt[lev]++;
  }
  sprintf(str, "Direct list");
  dumpSynCore1(str, depth, fcnt, exupcnt, exdowncnt, emcnt, tcnt);

  for(i = 0; i <= depth; i++) exupcnt[i] = exdowncnt[i] = fcnt[i] = emcnt[i] = tcnt[i] = 0;
  /* dump cube local list by level */
  for(i = 0; i <= depth; i++) {
    for(nc = sys->locallist[i]; nc != NULL; nc = nc->lnext) {
      lev = nc->level;
      tcnt[lev]++;
      if(nc->upnumvects > 0) fcnt[lev]++;
      else emcnt[lev]++;
      if(nc->exactUp == TRUE) exupcnt[lev]++;
      if(nc->exactDown == TRUE) exdowncnt[lev]++;
    }
  }
  sprintf(str, "Local list");
  dumpSynCore1(str, depth, fcnt, exupcnt, exdowncnt, emcnt, tcnt);
    
  for(i = 0; i <= depth; i++) exupcnt[i] = exdowncnt[i] = fcnt[i] = emcnt[i] = tcnt[i] = 0;
  /* dump cube multipole list by level */
  for(i = 0; i <= depth; i++) {
    for(nc = sys->multilist[i]; nc != NULL; nc = nc->mnext) {
      lev = nc->level;
      tcnt[lev]++;
      if(nc->upnumvects > 0) fcnt[lev]++;
      else emcnt[lev]++;
      if(nc->exactUp == TRUE) exupcnt[lev]++;
      if(nc->exactDown == TRUE) exdowncnt[lev]++;
    }
  }
  sprintf(str, "Multi list");
  dumpSynCore1(str, depth, fcnt, exupcnt, exdowncnt, emcnt, tcnt);


  sprintf(str, "Multis built");
  dumpSynCore2(str, depth, multicnt);

  sprintf(str, "Locals built");
  dumpSynCore2(str, depth, localcnt);

}

/*
  like dumpMat but different formating and row labels (for dumpMatBldCnts)
*/
void dumpMatCnts(mat, depth, type)
int **mat, depth;
char *type;
{
  int i, j;
  char str[BUFSIZ];

  fprintf(stdout,
	  "\n%s MATRIX BUILD TOTALS (row = from cube, col = to cube):\n", 
	  type);

  for(i = 0; i <= depth; i++) {
    sprintf(str, " to %d ", i);
    if(i == 0) fprintf(stdout, "%13s", str);
    else if(i < 10) fprintf(stdout, "%6s", str);
    else fprintf(stdout, "%7s", str);
  }
  fprintf(stdout, "\n");

  for(i = 0; i <= depth; i++) {
    sprintf(str, "from %d ", i);
    fprintf(stdout, "%-7s", str); /* print row label */
    for(j = 0; j <= depth; j++) {
      sprintf(str, "%d ", mat[i][j]);
      if(j < 10) fprintf(stdout, "%6s", str);
      else fprintf(stdout, "%7s", str);
    }
    fprintf(stdout, "\n");
  }

}

/*
  display matrix build count totals
*/
void dumpMatBldCnts(sys)
ssystem *sys;
{
  int i;
  char type[BUFSIZ];
  extern int **Q2Mcnt, **Q2Lcnt, **Q2Pcnt, **L2Lcnt;
  extern int **M2Mcnt, **M2Lcnt, **M2Pcnt, **L2Pcnt;

  sprintf(type, "Q2M");
  dumpMatCnts(Q2Mcnt, sys->depth, type);

  sprintf(type, "Q2L");
  dumpMatCnts(Q2Lcnt, sys->depth, type);

  sprintf(type, "Q2P");
  dumpMatCnts(Q2Pcnt, sys->depth, type);

  sprintf(type, "M2M");
  dumpMatCnts(M2Mcnt, sys->depth, type);

  sprintf(type, "M2L");
  dumpMatCnts(M2Lcnt, sys->depth, type);

  sprintf(type, "M2P");
  dumpMatCnts(M2Pcnt, sys->depth, type);

  sprintf(type, "L2L");
  dumpMatCnts(L2Lcnt, sys->depth, type);

  sprintf(type, "L2P");
  dumpMatCnts(L2Pcnt, sys->depth, type);

}

/* 
  dumps state of important compile flags and arg list parameters
*/
void dumpConfig(fp)
FILE *fp;
{

  fprintf(fp, "\n\nFastLap CONFIGURATION FLAGS:\n");

  fprintf(fp, "  General Configuration\n");
  fprintf(fp, "    NOWARN");
  if(NOWARN == ON) 
      fprintf(fp, " == ON (no warnings to stdout)\n");
  else if(NOWARN == OFF) 
      fprintf(fp, " == OFF (warnings written on stdout)\n");
  fprintf(fp, "    NOABRT");
  if(NOABRT == ON) 
      fprintf(fp, " == ON (slow fatal error traps disabled)\n");
  else if(NOABRT == OFF) 
      fprintf(fp, " == OFF (slow fatal error traps enabled)\n");

  fprintf(fp, "  Multipole Configuration\n");

  fprintf(fp, "    DNTYPE");
  if(DNTYPE == NOLOCL) 
      fprintf(fp, " == NOLOCL (no locals in dwnwd pass)\n");
  else if(DNTYPE == NOSHFT) 
      fprintf(fp, " == NOSHFT (no local2local shift dwnwd pass)\n");
  else if(DNTYPE == GRENGD) 
      fprintf(fp, " == GRENGD (full Greengard dwnwd pass)\n");
  fprintf(fp, "    MULTI");
  if(MULTI == ON) fprintf(fp, " == ON (include multipole part of P*q)\n");
  else fprintf(fp, " == OFF (don't use multipole part of P*q)\n");
  fprintf(fp, "    RADINTER");
  if(RADINTER == ON) 
      fprintf(fp," == ON (allow parent level interaction list entries)\n");
  else 
   fprintf(fp," == OFF (use only cube level interaction list entries)\n");
  fprintf(fp, "    NNBRS == %d (max distance to a nrst neighbor)\n", NNBRS);
  fprintf(fp, "    ADAPT");
  if(ADAPT == ON) 
      fprintf(fp, " == ON (adaptive - no expansions in exact cubes)\n");
  else fprintf(fp, " == OFF (not adaptive - expansions in all cubes)\n");
  fprintf(fp, "    OPCNT");
  if(OPCNT == ON) 
      fprintf(fp, " == ON (count P*q ops - exit after mat build)\n");
  else fprintf(fp, " == OFF (no P*q op count - iterate to convergence)\n");
  fprintf(fp, "    MAXDEP");
  fprintf(fp, 
	  " == %d (assume no more than %d partitioning levels are needed)\n",
	  MAXDEP, MAXDEP);

  fprintf(fp, "  Linear System Solution Configuration\n");

  fprintf(fp, "    ITRTYP");
  if(ITRTYP == GMRES)
      fprintf(fp, " == GMRES (generalized minimum residuals)\n");
  else fprintf(fp, " == %d (not implemented - use GMRES)\n", ITRTYP);

  fprintf(fp, "    PRECOND");
  if(PRECOND == OL) {
    fprintf(fp, " == OL (use overlapped preconditioner)\n");
  }
  else if(PRECOND == SP) {
    fprintf(fp, " == SP (use special preconditioner)\n");
  }
  else fprintf(fp, " == NONE (no preconditioner)\n");

}

/*
  pads a string on the right up to a given length, truncates if too long
*/
char *padName(tostr, frstr, len)
char *tostr, *frstr;
int len;
{
  int i;

  for(i = 0; frstr[i] != '\0'; i++) tostr[i] = frstr[i];
  if(i > len) tostr[len] = '\0';		/* truncate */
  else {			/* pad */
    for(; i < len; i++) tostr[i] = ' ';
    tostr[len] = '\0';
  }
  return(tostr);
}

/*
  returns a string of spaces (doesn't stdio have this somewhere?)
*/
char *spaces(str, num)
char *str;
int num;
{
  int i;

  for(i = 0; i < num; i++) str[i] = ' ';
  str[num] = '\0';
  return(str);
}
    
/*
  dumps brief information about multipole set up
*/
void dumpMulSet(sy, size, numLev, order)
ssystem *sy;
int size, numLev, order;
{
  int numcubes, numsides, i, multerms();

  for(numcubes = 1, i = 0; i < numLev; numcubes *= 8, i++);
  for(numsides = 1, i = 0; i < numLev; numsides *= 2, i++);

  fprintf(stdout, "\nMULTIPOLE SETUP SUMMARY:\n");
  fprintf(stdout, "  Level 0 cube extremal coordinates\n");
  fprintf(stdout, "    x: %g to %g\n", 
	  sy->minx, sy->minx + numsides * (sy->length));
  fprintf(stdout, "    y: %g to %g\n", 
	  sy->miny, sy->miny + numsides * (sy->length));
  fprintf(stdout, "    z: %g to %g\n", 
	  sy->minz, sy->minz + numsides * (sy->length));
  fprintf(stdout, "  Level %d (lowest level) cubes\n    %d total\n", 
	  numLev, numcubes);
  fprintf(stdout, 
	  "    side length = %g\n  Maximum number of singularities in any = %d\n",
	  sy->length, sy->maxq);
  fprintf(stdout, 
	  "  Maximum number of singularities treated exactly = %d (limit = %d)\n",
	  sy->maxq, multerms(order));
}
