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

double **Q2P(), **Q2PAlloc();
double **mulMulti2P(), **mulQ2Multi(), **mulMulti2Multi();
double **mulLocal2Local(), **mulLocal2P(), **mulQ2Local(), **mulMulti2Local();

int *localcnt, *multicnt, *evalcnt;	/* counts of builds done by level */

int **Q2Mcnt, **Q2Lcnt, **Q2Pcnt, **L2Lcnt; /* counts of xformation mats */
int **M2Mcnt, **M2Lcnt, **M2Pcnt, **L2Pcnt;

/*
 * mulMatDirect creates the matrices for the piece of the problem that is done
 * directly.
 */
mulMatDirect(sys)
ssystem *sys;
{
  cube *nc, *nnbr;
  int i, nummats, swapOnly;
  extern double dirtime;

  /* c2.0 The direct list contains those cubes which have field points and at 
     least one nearest neighbor.  A cube MAY be its own nearest neighbor, but
     since a cube cannot be a nearest neighbor unless it has singularities, 
     there can be cubes with field points which have NO nearest neighbors and
     therefore are not in the direct list.  This, of course, requires the 
     existence of a unique evaluation list. */
  starttimer;
  for(nc=sys->directlist; nc != NULL; nc = nc->dnext) {
    if(nc->TranslationComputed == FALSE) {
      nummats = nc->numnbrs;
      /* Allocate space for the vects and mats. */
      nc->directnumvects = nummats;
      CALLOC(nc->directq, nummats, double*, ON, AMSC);
      CALLOC(nc->directnumeles, nummats, int, ON, AMSC);
      CALLOC(nc->directmats, nummats, double**, ON, AMSC);
    }
  }

  /* Now place in the matrices. */

  /*cftk Why don't we embed this in the above loop? */
  for(nc=sys->directlist; nc != NULL; nc = nc->dnext) {
    for(nummats=0, i=0; i < nc->numnbrs; i++) {
      nnbr = nc->nbrs[i];
      ASSERT(nnbr->upnumvects > 0);
      nc->directq[nummats] = nnbr->upvects[0];
      nc->directnumeles[nummats] = nnbr->upnumeles[0];

      if((nc->TranslationComputed==TRUE) && (nnbr->TranslationComputed==TRUE))
	  swapOnly = TRUE;
      else swapOnly = FALSE;

      nc->directmats[nummats] = Q2P(nnbr->sngs, nnbr->upnumeles[0], 
				    nc->fpts, nc->numfieldpts,
				    swapOnly, nc->directmats[nummats]);
      nummats++;

#if DMTCNT == ON
      Q2Pcnt[nc->level][nnbr->level]++;
#endif
    }
  }
  stoptimer;
  dirtime += dtime;
}

/* This near picks up only the hamming distance one cubes. */    
#define HNEAR(nbr, nj, nk, nl) \
((ABS((nbr)->j - (nj)) + ABS((nbr)->k - (nk)) + ABS((nbr)->l - (nl))) <= 1)

/* This near picks up all 27 neighboring cubes (including the given cube). */
#define NEAR(nbr, nj, nk, nl) \
((ABS((nbr)->j - (nj)) <= 1) && \
 (ABS((nbr)->k - (nk)) <= 1) && \
 (ABS((nbr)->l - (nl)) <= 1))

/* This near picks only the diagonal, for testing. */
#define DNEAR(nbr, nj, nk, nl) \
(((nbr)->j == (nj)) && \
 ((nbr)->k == (nk)) && \
 ((nbr)->l == (nl)) )
	

/*
 * olmulMatPrecond is an overlapping block preconditioner for the general
 * case of non-square direct interaction matrices.  It expects the subject 
 * cube to be the zeroth member of the nearest neighbors list.  However, 
 * the logic used here is that even though the self cube is in the direct 
 * list because it has field points, it may not be in its own nearest neighbor 
 * list as it may have no singularities. 
 */
olmulMatPrecond(sys)
ssystem *sys;
{
  cube *nc, *pnbr, *qnbr;
  double **mat, **matT, **matW, **nmat;
  int i, j, k, l, m, jj, nummats;
  int maxsize, nrows, ncols, rowOffset, colOffset;

  /* Allocate space for the precond mats and figure out the max number 
     of singularities or field points in any local problem so we can allocate
     work space. */
  for(maxsize=0, nc=sys->precondlist; nc != NULL; nc = nc->pnext) {
    nummats = nc->numpnbrs;
    CALLOC(nc->precondmats, nummats, double**, ON, AMSC);
    nrows = 0;
    /*cftk what do we need nummats here for? just use k*/
    for(nummats=0, k=0; k < nc->numpnbrs; k++) {
      pnbr = nc->pnbrs[k];
      nc->precondmats[nummats++] = Q2PAlloc(nc->upnumeles[0],pnbr->numfieldpts);
      nrows += pnbr->numfieldpts;
    }
    ncols = 0;
    for(k = 0; k < nc->numqnbrs; k++) {
      qnbr = nc->qnbrs[k];
      ncols += qnbr->upnumeles[0];
    }
    maxsize = MAX(maxsize,MAX(ncols, nrows));
  }
#if MULDAT == ON
  printf("  Preconditioning matrix for inversion size = %d\n", maxsize);
#endif

/* Allocate 3 matrices big enough for any set of 27. */
  MALLOC(mat, maxsize, double*, ON, AMSC);
  MALLOC(matT, maxsize, double*, ON, AMSC);
  MALLOC(matW, maxsize, double*, ON, AMSC);
  for(i=0; i < maxsize; i++) {
    MALLOC(mat[i], maxsize, double, ON, AMSC);
    MALLOC(matT[i], maxsize, double, ON, AMSC);
    MALLOC(matW[i], maxsize, double, ON, AMSC);
  }

/* Now go fill-in a local-problem matrix. */
  for(nc=sys->precondlist; nc != NULL; nc = nc->pnext) {
    rowOffset = 0;
      /* You must zero the column offset here as you may never enter the if
	 block below and thus exit the loop with the previous value. */
    colOffset = 0;

    /* This loop is for each block of rows of the local-problem matrix. 
     To generate rows, a neighbor cube must have field points.*/
    for(k=0; k < nc->numpnbrs; k++) {
      pnbr = nc->pnbrs[k];
      nrows = pnbr->numfieldpts;
      /* This loop is for each block of columns within the row block. */
      for(colOffset = 0, l = 0; l < nc->numqnbrs; l++) {
	qnbr = nc->qnbrs[l];
	for(m=0; m < pnbr->numnbrs; m++) {
	  if(pnbr->nbrs[m] == qnbr) break;
	}
	nmat = pnbr->directmats[m];
	ncols = pnbr->directnumeles[m];
	ASSERT (qnbr->upnumeles[0] == pnbr->directnumeles[m])
	for(i = nrows - 1; i >= 0; i--) {
	  for(j = ncols - 1; j >= 0; j--) {
	    mat[rowOffset+i][colOffset+j] = nmat[i][j];
	    matT[colOffset+j][rowOffset+i] = nmat[i][j];
	  }
	}
	colOffset += ncols;
      }
      rowOffset += nrows;
    }
	
    if( colOffset == rowOffset ) {
    /* For square systems, just invert. */
      invert(mat, colOffset, NULL);
    }
    else {
    /* For non-square systems form the pseudo inverse. */
      for(i = 0; i < colOffset; i++) {
	for(j = 0; j < colOffset; j++) {
	  matW[i][j] = 0.0;
	  for(jj = 0; jj < rowOffset; jj++) {
	    matW[i][j] += matT[i][jj] * mat[jj][j];
	  }
	}
      }

      invert(matW, colOffset, NULL);

      for(i = 0; i < colOffset; i++) {
	for(j = 0; j < rowOffset; j++) {
	  mat[i][j] = 0.0;
	  for(jj = 0; jj < colOffset; jj++) {
	    mat[i][j] += matW[i][jj] * matT[jj][j];
	  }
	}
      }
    }
    /* Copy out the preconditioner to the precond matrices. */
    nrows = nc->upnumeles[0];
    colOffset = 0;
    for(k=0; k < nc->numpnbrs; k++) {
      pnbr = nc->pnbrs[k];
      ncols = pnbr->numfieldpts;
      nmat = nc->precondmats[k];
      for(i = nrows - 1; i >= 0; i--) {
        for(j = ncols - 1; j >= 0; j--) {
          nmat[i][j] = mat[i][colOffset + j];
	}
      }
      colOffset += ncols;
    }
  }
}

/*
 * spmulMatPrecond is a special overlapping block preconditioner for 
 * desingularized problems.  It avoids non-square direct interaction matrices by 
 * finding the "correct" field point for each source point in the local problem.
 * Consequently the field points are operated on in source point order so the 
 * precond mats built by this function do not make the ordering transposition 
 * that one would expect (field point order --> source point order).  
 */
spmulMatPrecond(sys)
ssystem *sys;
{
  cube *nc, *pnbr, *qnbr;
  snglrty *pq, *kq;
  fieldpt *pf, myfpt;
  double realdummy, calcp();
  double **nmat, *colmat;
  double *rwork;
  int lwork, info, *ipiv;
  int i, j, k, ks, kk, kks;
  int maxsize, rank, nrows, ncols, colOffset;
  /* Allocate space for the precond work mats by counting the max number 
     of singularities any local problem. */
  for(maxsize=0, nc=sys->precondlist; nc != NULL; nc = nc->pnext) {
    CALLOC(nc->precondmats, nc->numqnbrs, double**, ON, AMSC);
    rank = 0;
    for( k=0; k < nc->numqnbrs; k++ ) {
      qnbr = nc->qnbrs[k];
      nc->precondmats[k] = Q2PAlloc(nc->upnumeles[0],qnbr->upnumeles[0]);
      rank += qnbr->upnumeles[0];
    }
    maxsize = MAX(maxsize,rank);
  }
  lwork = maxsize*8;
#if MULDAT == ON
  printf(" Preconditioning work space matrix size = %d\n", maxsize);
#endif

/* Allocate a column "matrix" big enough for any set of 27. */
  MALLOC(colmat, maxsize*maxsize, double, ON, AMSC);
  MALLOC(rwork, lwork, double, ON, AMSC);
  MALLOC(ipiv, maxsize, int, ON, AMSC);

/* Now go fill-in a local-problem matrix. */
  for(nc=sys->precondlist; nc != NULL; nc = nc->pnext) {
    /* Get a source point. */
    rank = 0;
    j = 0;
    for( kk=0; kk < nc->numqnbrs; kk++ ) {
      qnbr = nc->qnbrs[kk];
      for( kks=0; kks<qnbr->upnumeles[0]; kks++ ) {
        pq = qnbr->sngs[kks];
	/* Get a field point. */
	for( k=0; k < nc->numqnbrs; k++ ) {
          pnbr = nc->qnbrs[k];
	  for( ks=0; ks<pnbr->upnumeles[0]; ks++ ) {
            pf = (pnbr->sngs[ks])->myfpt;
	    realdummy = calcp(pq, pf->x, pf->y, pf->z, pf->deriv, 
				pf->nrm, &colmat[j++]);
	  }
	}
	rank++;
      }
    }

    /* If Lapack not available, load mat[i][j] for: invert(mat, j, NULL);*/
    dgetrf_(&rank,&rank,colmat,&rank,ipiv,&info);
    dgetri_(&rank,colmat,&rank,ipiv,rwork,&lwork,&info);

    /* Copy out the preconditioner to the precond matrices. */
    nrows = nc->upnumeles[0];
    colOffset = 0;
    for(k=0; k<nc->numqnbrs; k++) {
      pnbr = nc->qnbrs[k];
      ncols = pnbr->upnumeles[0];
      nmat = nc->precondmats[k];
      for(j=0; j<ncols; j++) {
        for(i=0; i<nrows; i++) {
          nmat[i][j] = colmat[(j+colOffset)*rank+i];
	}
      }
      colOffset += ncols;
    }
  }
}

/* 
 * MulMatUp computes the multipole to multipole or snglrty to
 * multipole matrices that map to a parent's multipole coeffs from its
 * children's multipoles or snglrtys. Note that only one set of
 * multipole to multipole matrices is computed per level by exploiting the
 * uniform break-up of three-space (ie many shifts have similar geometries).  
 */
mulMatUp(sys) 
ssystem *sys; 
{
  cube *nextc, *kid;
  int i, j, numterms, depth, order = sys->order;
  double **multimats[8];

  for(i=0; i < 8; i++) multimats[i] = NULL;

  numterms = multerms(order);

  if(sys->depth < 2) {
#if NOWARN == OFF
    fprintf(stdout, "FLW-mulMatUp: no multipole acceleration at all\n");
#endif
    return;	/* return if upward pass not possible */
  }

  /* Handle the lowest level cubes first (set up Q2M's). */
  for(nextc=sys->multilist[sys->depth]; nextc != NULL; nextc = nextc->mnext) {
    if(nextc->TranslationComputed == FALSE) {
      CALLOC(nextc->upmats, 1, double**, ON, AMSC);
    }
    nextc->upmats[0] = mulQ2Multi(nextc->sngs, nextc->upnumeles[0],
				  nextc->x, nextc->y, nextc->z, 
				  nextc->TranslationComputed,
				  order, nextc->upmats[0]);

#if DISSYN == ON
    multicnt[nextc->level]++;
#endif

#if DMTCNT == ON
    Q2Mcnt[nextc->level][nextc->level]++;
#endif

  }

#if MULDAT == ON
  if(sys->multilist[sys->depth] == NULL) {
    fprintf(stdout, 
	    "  No multipole expansions at level %d (lowest)\n", sys->depth);
  }
#endif
#if NOWARN == OFF
  if((sys->multilist[sys->depth] == NULL) && (sys->depth < 3)) {
    fprintf(stdout, "FLW-mulMatUp: no multipole acceleration at all.\n");
  }
#endif

  /* Allocate the vectors and matrices for the cubes. */
  /* No multipoles over root cube or its kids (would not be used if made). */
  for(depth = (sys->depth - 1); depth > 1; depth--) {
    /* Set up M2M's and Q2M's to compute multipoles needed for this level. */
#if MULDAT == ON
    if(sys->multilist[depth] == NULL) {
      fprintf(stdout, 
	      "  No multipole expansions at level %d\n", depth); 
    }
#endif
#if NOWARN == OFF
    if((sys->multilist[depth] == NULL) && (depth < 3)) {
      fprintf(stdout, "FLW-mulMatUp: no multipole acceleration at all.\n");
    }
#endif

    /* NULL out pointers to same-geometry M2M mats for this level */
    for(i=0; i < 8; i++) multimats[i] = NULL;

    /* Hit nonempty cubes at this level assigning ptrs to precomputed   */
    /* M2M mats (for this lev), or if kid is exact, computing Q2M matrices. */
    for(nextc=sys->multilist[depth]; nextc != NULL; nextc = nextc->mnext) {
      
#if DISSYN == ON
      multicnt[nextc->level]++;
#endif

      /* Save space for upvector sizes, upvect ptrs, and upmats. */
      if(nextc->TranslationComputed == FALSE) {
	CALLOC(nextc->upnumeles, nextc->upnumvects, int, ON, AMSC);
	CALLOC(nextc->upvects, nextc->upnumvects, double*, ON, AMSC);
	CALLOC(nextc->upmats, nextc->upnumvects, double**, ON, AMSC);
      }
      /* Go through nonempty kids and fill in upvectors and upmats. */
      for(i=0, j=0; j < nextc->numkids; j++) {
	/*2.0 Not only must kid exist, but kid must be in upward pass tree.*/
	if(((kid = nextc->kids[j])!=NULL)&&(nextc->kids[j]->upcube==TRUE)) {
	  if((kid->exactUp==FALSE)&&(nextc->TranslationComputed==FALSE)) {
	    /* if kid has a multi, compute translation on first pass only. */
	    nextc->upvects[i] = kid->multi;
	    nextc->upnumeles[i] = kid->multisize;
	    if(multimats[j] == NULL) { /* Build needed matrix only once. */
	      multimats[j] = mulMulti2Multi(kid->x,kid->y, kid->z,nextc->x, 
					    nextc->y,nextc->z,order,NULL);
	    }
	    nextc->upmats[i] = multimats[j];

#if DMTCNT == ON
	    M2Mcnt[kid->level][nextc->level]++; /* cnts use, ~computation */
#endif				/* only at most 8 mats really built/level */
	  }
	  else if(kid->exactUp == TRUE) { /* If kid is exact, it has no multi */
	    nextc->upvects[i] = kid->upvects[0];
	    nextc->upnumeles[i] = kid->upnumeles[0];
	    nextc->upmats[i] = mulQ2Multi(kid->sngs, kid->upnumeles[0],
					  nextc->x, nextc->y, nextc->z, 
					  nextc->TranslationComputed,
					  order, nextc->upmats[i]);
#if DMTCNT == ON
	    Q2Mcnt[kid->level][nextc->level]++;
#endif

	  }
	  i++;			/* only increments if kid is not empty */
	}
      }
    }
  }
}

/*
 * mulMatEval builds the transformation matrices for the final eval pass 
 * (M2P,L2P) for all cubes in the evaluation list generated by linkcubes().
 * These are the lowest level cubes which contain evaluation points.
 *
 * For each cube A in the evaluation list:
 * 1) if A is not exactDown (always the case if ADAPT = OFF),
 *    a) and if DNTYPE = GRENGD, build an L2P matrix from A to A,
 *    b) and if DNTYPE = NOSHFT, build an L2P matrix from each of A's 
 *       ancestors with level > 1 (including A) to A),
 *    c) and if DNTYPE = NOLOCL build an M2P matrix from each of A's fake 
 *       ilist entries to A (same action as 2b);
 * 2) if A is exactDown, find the 1st ancestor of A, cube B, 
 *    which either is not exactDown and is at level 2,3,4... or is at level 1
 *    a) if B is at level 2,3,4...
 *       i) if DNTYPE = GRENGD, construct an L2P from B to A and M2P's
 *          from the cubes in the true interaction lists of A and all its
 *	   ancestors up to and including B (a partial fake interaction list)
 *	j) if DNTYPE = NOSHFT, find cube C, the ancestor of B at level 1;
 *	   construct L2P's from the ancestors of B (including B but not C)
 *	   to A and Q- or M2P's from the cubes in the true interaction lists 
 *	   of A and all its ancestors up to and including B (a partial fake 
 *	   interaction list)
 *	k) if DNTYPE = NOLOCL, do 2b
 *    b) if B is at level 1 construct M2P's for all the cubes in A's
 *       fake interaction list
 *
 * True interaction list - RADINTER = OFF, those sibling (same level) cubes 
 * of a given cube who are children of the neighbors of the given cube's 
 * parent and are not neighbors of the given cube 
 * - ie those cubes required to cover snglrtys well separated from the given
 * cube but not accounted for in the parent's local expansion 
 * - the flag NNBRS is the number of sibling cube "shells" taken as neighbors 
 *  
 * fake interaction list - RADINTER = OFF, the combined true interaction lists
 * of a given cube and all its ancestors at levels 2,3,4...
 *
 * if RADINTER = ON, any 8 siblings of the given cube which form a well 
 * separated cube one level up are included in the lists as a single higher
 * level cube
 *  
 * if ADAPT = OFF, no cube is exact so step 2 is never done
 *
 * this routine is used alone if compiled with DNTYPE = NOLOCL or after
 * mulMatDown, which produces M2L and L2L matrices (DNTYPE = GRENGD) or
 * just M2L matrices (DNTYPE = NOSHFT) --  DNTYPE = GRENGD does the full
 * Greengard hiearchical downward pass
 *
 */
void mulMatEval(sys)
ssystem *sys;
{
  int i, j, k, ttlvects, vects;
  cube *na, *nc, *nexti;

  if(sys->depth < 2) return;	/* ret if upward pass not possible/worth it */

  /* c2.0 We work with the evaluation list here =/= direct list in general. */
  for(nc = sys->evallist; nc != NULL; nc = nc->enext) {

    ASSERT(nc->level == sys->depth);
/*c2.0 Can't assert this in general:    ASSERT(nc->upnumvects > 0);*/

    /* First count the number of transformations to do.  At the very least, 
       if the cube is not exactDown, and we are doing locals and shifting 
       locals, there will be the cube's own L2P.  If the cube is exactDown, 
       and we are doing locals, then we go up the tree thru the cube's 
       parents counting up their interaction list contributions which will 
       be evaluated in the cube. Then allocate based on this count. */ 
    for(na = nc, ttlvects = 0; na->level > 1; na = na->parent) { 
      if((na->exactDown == FALSE) && (DNTYPE != NOLOCL)) {
	ttlvects++;  /* allow for na to na local expansion (L2P) */
	if(DNTYPE == GRENGD) break; /* Only one local exp if shifting. */
      }
      else {
	ttlvects += na->interSize; /* room for Q2P and M2P xformations */
      }
    }

    if(nc->TranslationComputed == FALSE) {
      nc->evalnumvects = ttlvects; /* save ttl # of transformations to do */
      CALLOC(nc->evalvects, ttlvects, double*, ON, AMSC);
      CALLOC(nc->evalnumeles, ttlvects, int, ON, AMSC);
      CALLOC(nc->evalmats, ttlvects, double**, ON, AMSC);
    }

#if DILIST == ON
    fprintf(stdout, "\nInteraction list (%d entries) for ", ttlvects);
    disExParsimpcube(nc);
#endif

    /* Having counted and allocated everything, we now load the matrices
       based on the same logic as above.  Note that we check the ilist 
       entries to see if they are exact up, in which case we use Q2P, 
       rather than M2P. */
    for(j=0, na = nc, ttlvects = 0; na->level > 1; na = na->parent) { 
      if((na->exactDown == FALSE) && (DNTYPE != NOLOCL)) {  
	if(nc->TranslationComputed == FALSE) {
	  /* Build matrices for local expansion evaluation. */
	  nc->evalmats[j] = mulLocal2P(na->x, na->y, na->z, nc->fpts,
				       nc->numfieldpts, sys->order);
	  nc->evalnumeles[j] = na->localsize;
	  nc->evalvects[j] = na->local;
	  
#if DMTCNT == ON
	  L2Pcnt[na->level][nc->level]++;
#endif
	
#if DILIST == ON
	  fprintf(stdout, "L2P: ");
	  disExtrasimpcube(na);
#endif
	}
	j++; 
	/* Only one local exp if shifting. */
	if(DNTYPE == GRENGD) break; 
      }
      /* Build matrices for ancestor's (or cube's if 1st time) ilist */
      else {
	for(i=0; i < na->interSize; i++) {
	  nexti = na->interList[i];
	  /*c2.0 Here we check if the ilist member is exactUp. */
	  if(nexti->exactUp == TRUE) {
	    nc->evalvects[j] = nexti->upvects[0];
	    nc->evalmats[j] = Q2P(nexti->sngs, nexti->upnumeles[0], 
				  nc->fpts, nc->numfieldpts,
				  nc->TranslationComputed, nc->evalmats[j]);
	    nc->evalnumeles[j] = nexti->upnumeles[0];
	    
#if DMTCNT == ON
	    Q2Pcnt[nexti->level][nc->level]++;
#endif

#if DILIST == ON
	    fprintf(stdout, "Q2P: ");
	    disExtrasimpcube(nexti);
#endif
	    j++;
	  }
	  else {
	    if(nc->TranslationComputed == FALSE) {
	      nc->evalvects[j] = nexti->multi;
	      nc->evalmats[j] = mulMulti2P(nexti->x, nexti->y, nexti->z, 
					   nc->fpts, nc->numfieldpts, 
					   sys->order);
	      nc->evalnumeles[j] = nexti->multisize;
	    
#if DMTCNT == ON
	      M2Pcnt[nexti->level][nc->level]++;
#endif

#if DILIST == ON
	      fprintf(stdout, "M2P: ");
	      disExtrasimpcube(nexti);
#endif
	    }
	    j++;
	  }
	}
      }
    }
  }
}


/* 
 * mulMatDown sets up the matrices for the downward pass.
 * For each cube in the local list (parents always in list before kids):
 * 1) parent's local to child's local unless DNTYPE=NOSHFT or no parent local,
 * 2) multipoles for (Parent+parent's nbrs - child nbrs) to child's local.
 * -eval is sum of ancestral local evals for each lowest lev cube if NOSHFT
 * otherwise only lowest level local is evaluated (see mulMatEval).
 * -With ADAPT = OFF no cube is exactDown so local list is all cubes at lev>1
 * in which there are fieldpoints.
 * -mats that give potentials (M2P, L2P, Q2P) are calculated in mulMatEval()
 * -this routine makes only L2L, M2L and Q2L matrices.
 */
mulMatDown(sys)
ssystem *sys;
{
  int i, j, vects;
  cube *nc, *parent, *ni;
  int depth;

  ASSERT(DNTYPE != NOLOCL);	/* use mulMatEval() alone if NOLOCL */

  for(depth = 2; depth <= sys->depth; depth++) { /* no locals before level 2 */
    for(nc=sys->locallist[depth]; nc != NULL; nc = nc->lnext) {

      /* Allocate for interaction list, include one for parent if needed. */
      if((depth <= 2) || (DNTYPE == NOSHFT)) vects = nc->interSize;
      else vects = nc->interSize + 1;
      if(nc->TranslationComputed == FALSE) {
	nc->downnumvects = vects;
	CALLOC(nc->downvects, vects, double*, ON, AMSC);
	CALLOC(nc->downnumeles, vects, int, ON, AMSC);
	CALLOC(nc->downmats, vects, double**, ON, AMSC);
      }
      parent = nc->parent;
      /* A cube processed here cannot have an exact parent because if it did 
	 then it would have to be exact and would not be in the list. (Its 
	 assets would have been promoted.)*/
      ASSERT(parent->exactDown == FALSE); /* cube has >= #fpts of any of its kids*/

#if DISSYN == ON
      localcnt[nc->level]++;
#endif

      if((depth <= 2) || (DNTYPE == NOSHFT)) i = 0; /* No parent local. */
      else { /* Create the mapping matrix for the parent to kid. */
	i = 1;
	if(nc->TranslationComputed == FALSE) {
	  nc->downmats[0] = mulLocal2Local(parent->x, parent->y, parent->z,
					   nc->x, nc->y, nc->z, sys->order);
	  nc->downnumeles[0] = parent->localsize;
	  nc->downvects[0] = parent->local;

#if DMTCNT == ON
	  L2Lcnt[parent->level][nc->level]++;
#endif
	}
      }

      /* Go through the interaction list and create mapping matrices. */
      for(j = 0; j < nc->interSize; j++, i++) {
	ni = nc->interList[j];
	if(ni->exactUp == TRUE) { 
	  nc->downvects[i] = ni->upvects[0];
	  nc->downmats[i] = mulQ2Local(ni->sngs, ni->upnumeles[0], 
				       nc->x, nc->y, nc->z, sys->order, 
				       nc->TranslationComputed, 
				       nc->downmats[i]);
	  nc->downnumeles[i] = ni->upnumeles[0];
#if DMTCNT == ON
	  Q2Lcnt[ni->level][nc->level]++;
#endif
	}
	else {
	  if(nc->TranslationComputed == FALSE) {
	    nc->downvects[i] = ni->multi;
	    nc->downmats[i] = mulMulti2Local(ni->x, ni->y, ni->z, nc->x, 
					     nc->y, nc->z, sys->order, NULL);
	    nc->downnumeles[i] = ni->multisize;
#if DMTCNT == ON
	    M2Lcnt[ni->level][nc->level]++;
#endif
	  }
	}
      }
    }
  }
}
