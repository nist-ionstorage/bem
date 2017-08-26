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

/*
 * ssystem sets up the spatial hierarchy for snglrtys and expansions.
 */
ssystem *mulInit(autom, depth, order, snglrtys, lhsSize, rhsSize, fieldpts)
int autom;
int depth, order, lhsSize, rhsSize;
snglrty *snglrtys;
fieldpt *fieldpts;
{
  ssystem *sys;
  int qindex=1, pindex=1, cindex=1;

  CALLOC(sys, 1, ssystem, ON, AMSC);
  sys->depth = depth;	       /* Overwritten below if autom = ON*/
  sys->order = order;

  /* Create cubes, put singularities and field points in leaves. */
  sys->depth = placeq(autom, sys, snglrtys, lhsSize, rhsSize, fieldpts); 

  /* Get all the parents and kids for each cube. */
  getrelations(sys);		

  /* Figures out position of cube center. */
  setPosition(sys);		

  /* Index singularities, field points, and cubes.*/
  indexkid(sys, sys->cubes[0][0][0][0], &qindex, &pindex, &cindex); 

  /* Find diagonal relationship, field points to sings. */
#if PRECOND == SP
  getPairs(snglrtys, fieldpts);
#endif

  /* Note cubes done exactly and find # of nonempty kids for each cube. */  
  /*cftk note that the method of turning off adaptive, i.e. calling this with
    the number of terms = 0 requires that the real number of terms be computed
    inside setExact..., this is pretty clumsy, perhaps a flag? */
#if ADAPT == ON
  setExactUp(sys, multerms(sys->order)); 
  setExactDown(sys, multerms(sys->order));
#else
  setExactUp(sys, 0);
  setExactDown(sys, 0); 
#endif

  /* Get nbrs-type nearest neighbors with sing's of the downward-pass cubes. 
     This includes cubes due to exact upward-pass ancestor. */
  getnbrs(sys);
				   
  /* Make linked-lists of direct, multis, and locals to do at each level. */
  linkcubes(sys);

  /* If preconditioning, get lowest level 1st neighbors with field points 
     (pnbrs) for OL or with singularities (qnbrs) for SP. */
#if PRECOND == OL
  getqnbrs(sys);
  getpnbrs(sys);
#endif
#if PRECOND == SP
  getqnbrs(sys);
#endif




  /* Get max # singularities or fieldpoints (2.0) in cubes treated exactly 
     so that we can allocate for expansion quantities in L2P and Q2P, etc. */
  setMaxq(sys);                 

  /* Get the interaction lists at all levels. */
  getAllInter(sys);		

  return(sys);
}

/*
  placeq defines the level 0 cube and places the snglrtys and field points in
  it.  If auto depth is on (flag==ON) it will try to use the best number of 
  levels but the algorithm is dumb.  The return is number of levels used.
*/
static int placeq(flag, sys, snglrtys, lhsSize, rhsSize, fieldpts)
int flag, rhsSize, lhsSize;
ssystem *sys;
snglrty *snglrtys;
fieldpt *fieldpts;
{
  int i, j, k, l, side, isexact, multerms(), depth, maxSize;
  int xindex, yindex, zindex, limit = multerms(sys->order);
  double length0, length, maxTileLength;
  double minx, maxx, miny, maxy, minz, maxz, tilelength();
  snglrty *nextq, *compq;
  fieldpt *nextf, *nf;
  cube *****cubes, *nc;
  int oversize;

  /* Figure out the corner coordinate x<0, y<0, z<0 so that all panel 
     centroids are at x_c>x, y_c>y, z_c>z. This establishes the panel 
     constraint on the "lower left" corner.*/
  nextq = snglrtys;
  minx = maxx = nextq->x;
  miny = maxy = nextq->y;
  minz = maxz = nextq->z;

  for(nextq = nextq->next; nextq != NULL; nextq = nextq->next) {
    maxx = MAX(nextq->x, maxx);
    minx = MIN(nextq->x, minx);
    maxy = MAX(nextq->y, maxy);
    miny = MIN(nextq->y, miny);
    maxz = MAX(nextq->z, maxz);
    minz = MIN(nextq->z, minz);
  }

  /* Figure out the corner coordinate x<0, y<0, z<0 so that all field 
     points are at x_field>x, y_field>y, z_field>z. This establishes the 
     field point constraint on the "lower left" corner.*/
  nextf = fieldpts;

  for(nextf = nextf->next; nextf != NULL; nextf = nextf->next) {
    maxx = MAX(nextf->x, maxx);
    minx = MIN(nextf->x, minx);
    maxy = MAX(nextf->y, maxy);
    miny = MIN(nextf->y, miny);
    maxz = MAX(nextf->z, maxz);
    minz = MIN(nextf->z, minz);
  }

  sys->minx = minx;
  sys->miny = miny;
  sys->minz = minz;

  /* Now dimension level 0 cube so that it includes all panels. */

  length0 = MAX((maxx - minx), (maxy - miny));
  length0 = MAX((maxz - minz), length0);

  /* Create the vectors for storing the trial strengths and field values. */
  /*cftk Why do we need to allocate size + 1?*/
  /* We have to stay flexible and take the max on these as they do different
     things depending on the job. */
  maxSize = MAX(rhsSize, lhsSize);
  CALLOC(sys->q, maxSize+1, double, ON, AMSC);
  CALLOC(sys->p, maxSize+1, double, ON, AMSC);

  /* Auto depth is not supported and so the setting for it is trapped here.
   * Page down a bit to find the pre-set depth section.
   */
  if(flag == ON) {		/* set depth of partitions automatically */
    /*c2.0 Auto depth not supported at present.*/
    fprintf(stderr, 
	    "FLE-placeq: Sorry, but auto-depth is not supported.");
    exit(0);
    /* alloc spine of cube pntr array - leave enough room for depth = MAXDEP */
    CALLOC(cubes, MAXDEP+1, cube****, ON, AMSC); 

    /* allocate for levels 0, 1, and 2 (always used) */
    for(side = 1, i=0; i <= 2; side *= 2, i++) {
      CALLOC(cubes[i], side, cube***, ON, AMSC);
      for(j=0; j < side; j++) {
	CALLOC(cubes[i][j], side, cube**, ON, AMSC);
	for(k=0; k < side; k++) {
	  CALLOC(cubes[i][j][k], side, cube*, ON, AMSC);
	}
      }
    }
    /* side /= 2; */

    /* for each level > 2: allocate for full cubes, count sing.'s in each,
       quit loop if all lowest level cubes are exact */
    for(isexact = FALSE; isexact == FALSE; side *= 2, i++) {

      if(i > MAXDEP) {
	fprintf(stderr, 
		"FLE-placeq: levels required beyond MAXDEP == %d\n", 
		MAXDEP);
	exit(0);
      }

      length = (1.01 * length0)/side;

      CALLOC(cubes[i], side, cube***, OFF, AMSC);
      if(cubes[i] == NULL) {
	fprintf(stderr, "FLE-placeq: %d levels set up\n", i-1);
	exit(0);
      }
      for(j=0; j < side; j++) {
	CALLOC(cubes[i][j], side, cube**, OFF, AMSC);
	if(cubes[i][j] == NULL) {
	  fprintf(stderr, "FLE-placeq: %d levels set up\n", i-1);
	  exit(0);
	}
	for(k=0; k < side; k++) {
	  CALLOC(cubes[i][j][k], side, cube*, OFF, AMSC);
	  if(cubes[i][j][k] == NULL) {
	    fprintf(stderr, "FLE-placeq: %d levels set up\n", i-1);
	    exit(0);
	  }
	}
      }

      /* Count the number of sngularities per cube and allocate if needed */
      for(nextq = snglrtys; nextq != NULL; nextq = nextq->next) {
	xindex = (nextq->x - minx) / length;
	yindex = (nextq->y - miny) / length;
	zindex = (nextq->z - minz) / length;
	nc = cubes[i][xindex][yindex][zindex];
	if(nc == NULL) {
	  CALLOC(nc, 1, cube, OFF, AMSC);
	  if(nc == NULL) {
	    fprintf(stderr, "FLE-placeq: %d levels set up\n", i-1);
	    exit(0);
	  }
	  cubes[i][xindex][yindex][zindex] = nc;
	  nc->upnumvects = 1;
	  CALLOC(nc->upnumeles, 1, int, OFF, AMSC);
	  if(nc->upnumeles == NULL) {
	    fprintf(stderr, "FLE-placeq: %d levels set up\n", i-1);
	    exit(0);
	  }
	  nc->upnumeles[0] = 1;
	}
	else {
	  nc->upnumeles[0]++;
	}
      }

      /* if the current lowest level is not exact, loop back until it is */
      isexact = TRUE;
      for(j = 0; j < side; j++) {
	for(k = 0; k < side; k++) {
	  for(l = 0; l < side; l++) {
	    if(cubes[i][j][k][l] != NULL) {
	      if(cubes[i][j][k][l]->upnumeles[0] > limit) isexact = FALSE;
	    }
	  }
	}
      }
      /* clean up cube structs if need to go down another level */
      if(isexact == FALSE) {
	for(j = 0; j < side; j++) {
	  for(k = 0; k < side; k++) {
	    for(l = 0; l < side; l++) {
	      if(cubes[i][j][k][l] != NULL) {
		cubes[i][j][k][l]->upnumeles[0] = 0;
		cubes[i][j][k][l]->upnumvects = 0;
	      }
	    }
	  }
	}
      }
    }
    depth = i - 1;		/* the automatically set depth */
    side /= 2;
  }

  /* Here is where the preset depth part starts, using the preset depth, we 
   * try to insert singularities and field points.
   */
  else {
    /* Allocate the cubes, note calloc used because zeros everything. */
    depth = sys->depth;
    CALLOC(cubes, sys->depth+1, cube****, ON, AMSC);
    for(side = 1, i=0; i <= depth; side *= 2, i++) {
      CALLOC(cubes[i], side, cube***, ON, AMSC);
      for(j=0; j < side; j++) {
	CALLOC(cubes[i][j], side, cube**, ON, AMSC);
	for(k=0; k < side; k++) {
	  CALLOC(cubes[i][j][k], side, cube*, ON, AMSC);
	}
      }
    }
    side /= 2;
    length = (1.01 * length0)/side;

    /* Count the number of snglrtys per cube. */
    for(nextq = snglrtys; nextq != NULL; nextq = nextq->next) {
      xindex = (nextq->x - minx) / length;
      yindex = (nextq->y - miny) / length;
      zindex = (nextq->z - minz) / length;
      nc = cubes[depth][xindex][yindex][zindex];
      if(nc == NULL) {
	/* We've never seen this cube, allo. for cube & singularity counter. */
	CALLOC(nc, 1, cube, ON, AMSC);
	cubes[depth][xindex][yindex][zindex] = nc;
	CALLOC(nc->upnumeles, 1, int, ON, AMSC);
	/* c2.0 This leaf cube in upward, but maybe not downward pass tree. */ 
	nc->upcube = TRUE;
	nc->downcube = FALSE;
	nc->upnumeles[0] = 1;
/*cftk seems to me upnumvects is somewhat redundant and could be purged.*/
	nc->upnumvects = 1;
	/* c2.0 Initialize field point counter to zero. */
	nc->numfieldpts = 0;
      }
      else {
	/* We've seen this cube for sings. already. */
	nc->upnumeles[0]++;
      }
    }

    /* Count the number of field points per cube. */
    for(nextf = fieldpts; nextf != NULL; nextf = nextf->next) {
      xindex = (nextf->x - minx) / length;
      yindex = (nextf->y - miny) / length;
      zindex = (nextf->z - minz) / length;
      nc = cubes[depth][xindex][yindex][zindex];
      if(nc == NULL) {
	/* We've never seen this cube, allo. for cube. */
	CALLOC(nc, 1, cube, ON, AMSC);
	cubes[depth][xindex][yindex][zindex] = nc;
	/*2.0 This leaf cube in downward, and NOT in upward pass tree. */ 
	nc->downcube = TRUE;
	nc->upcube = FALSE;
	nc->numfieldpts = 1;
      }
      else {
	/* We've seen this cube either for sing's. or f'p'ts. already. */
	nc->downcube = TRUE;
	nc->numfieldpts++;
      }
    }
  }
  sys->length = length;
  sys->side = side;
  sys->cubes = cubes;

  /* Allocate space for the singularities and the field points. */
  for(j=0; j < side; j++) {
    for(k=0; k < side; k++) {
      for(l=0; l < side; l++) {
        nc = sys->cubes[depth][j][k][l];
        if(nc != NULL) {  /* Only fill out nonempty cubes. */
	  /* Allocate for the snglrty ptrs. */
	  if(nc->upnumeles != NULL) { /* Only allocate if cube has sings. */
	    CALLOC(nc->sngs, nc->upnumeles[0], snglrty*, ON, AMSC);
	    /* Zero the numsngs to use as index. */
	    nc->upnumeles[0] = 0;
	  }
	  /* Allocate for the field point ptrs, */
	  if(nc->numfieldpts != 0) {  /* Only allocate if cube has fpts. */
	    CALLOC(nc->fpts, nc->numfieldpts, fieldpt*, ON, AMSC);
	    /* Zero the numfieldpts to use as index. */
	    nc->numfieldpts = 0;
	  }
	}
      }
    }
  }

  /* Put the singularities in cubes; check their size against cube size. */
  for(oversize = 0, nextq = snglrtys; nextq != NULL; nextq = nextq->next) {
#if NOWARN == OFF
    if(tilelength(nextq) > length) {
      if(oversize++ < 10) {
	printf("FLW-placeq: oversized panel, cube length=%g panel length=%g\n",length, tilelength(nextq));
      }
    } 
#endif
    xindex = (nextq->x - minx) / length;
    yindex = (nextq->y - miny) / length;
    zindex = (nextq->z - minz) / length;
    nc = cubes[depth][xindex][yindex][zindex];
    nc->sngs[nc->upnumeles[0]++] = nextq;
  }

  /* Put the field points into the cubes.*/
  for(nextf = fieldpts; nextf != NULL; nextf = nextf->next) {
    xindex = (nextf->x - minx) / length;
    yindex = (nextf->y - miny) / length;
    zindex = (nextf->z - minz) / length;
    nc = cubes[depth][xindex][yindex][zindex];
    nc->fpts[nc->numfieldpts++] = nextf;
  }
  return(depth);
}
      

/*
  getRelations allocates parents links the children.  Note that the kids are
  put in their unique position as if all kids might be allocated.  This is 
  done in order to match a kid with the correct translation operator.  There
  are eight such operators and they are parent independent, rather than 
  redundantly storing one for each kid.  So its a trade-off, and one might
  like to explore the posibility of condensing the kid pointers and storing
  pointers to the correct translation operator. 
*/

getrelations(sys)
ssystem *sys;
{
cube *nextc, *parent, *****cubes = sys->cubes;
int i, j, k, l, side;
  for(i = sys->depth, side = sys->side; i >= 0; i--, side /= 2) {
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
          nextc = cubes[i][j][k][l];
	  if(nextc != NULL) {
	/* Get the parents and children pointers of nonempty cubes. */
	    if(i < sys->depth) {
	      nextc->numkids = 8; /* all cubes, even empties, are counted */
	      CALLOC(nextc->kids, nextc->numkids, cube*, ON, AMSC);
	      nextc->kids[0] = cubes[i+1][2*j][2*k][2*l]; /* empties get */
	      nextc->kids[1] = cubes[i+1][2*j][2*k][2*l+1]; /* null pointers */
	      nextc->kids[2] = cubes[i+1][2*j][2*k+1][2*l];
	      nextc->kids[3] = cubes[i+1][2*j][2*k+1][2*l+1];
	      nextc->kids[4] = cubes[i+1][2*j+1][2*k][2*l];
	      nextc->kids[5] = cubes[i+1][2*j+1][2*k][2*l+1];
	      nextc->kids[6] = cubes[i+1][2*j+1][2*k+1][2*l];
	      nextc->kids[7] = cubes[i+1][2*j+1][2*k+1][2*l+1];
	    }
	    if(i > 0) {
	      parent = cubes[i-1][j/2][k/2][l/2];
	      if(parent == NULL) {
		CALLOC(parent, 1, cube, ON, AMSC);
		cubes[i-1][j/2][k/2][l/2] = parent;
		/*2.0 By default, a parent is in neither tree. */
		parent->upcube = FALSE;
		parent->downcube = FALSE;
	      }
	      /*2.0 Then it joins the trees its children are in. */
	      nextc->parent = parent;
	      if (nextc->upcube == TRUE) {
		parent->upcube = TRUE;
	      }
	      if (nextc->downcube == TRUE) {
		parent->downcube = TRUE;
	      }
	    }
	  }
	}
      }
    }
  }
}

/*
 setPosition sets the position coordinates of the cubes.
*/
setPosition(sys)
ssystem *sys;
{
int i, j, k, l;
int side = sys->side;
double length = sys->length;
cube *nextc;

/* Mark the position of the lowest level cubes. */
  for(i=sys->depth; i >= 0; i--, side /= 2, length *= 2.0) {
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
	  nextc = sys->cubes[i][j][k][l];
	  if(nextc != NULL) {
	    nextc->x = length * ((double) j + 0.5) + sys->minx;
	    nextc->y = length * ((double) k + 0.5) + sys->miny;
	    nextc->z = length * ((double) l + 0.5) + sys->minz;
	    nextc->level = i;
	    nextc->j = j;
	    nextc->k = k;
	    nextc->l = l;
	  }
	}
      }
    }
  }
}

/*
  Recursive routine to give indexes to the snglrtys (& field points) so that 
  those in each cube are contiguous. In addition, insure that the snglrtys 
  (& field points) in each parent at each level are numbered contiguously.  
  This is used to support a psuedo-adaptive scheme.  Also get the pointer to 
  the appropriate section of the singularity strength and field vectors.  
  Uses the eval vector for the field coeffs at the lowest level.  Also index 
  the lowest level cubes.
*/
static indexkid(sys, dad, qindex, pindex, pcindex)
ssystem *sys;
cube *dad;
int *qindex, *pindex, *pcindex;
{
  int i;
  
  if(dad != NULL) {
    if(dad->numkids == 0) {
      if(dad->upnumeles != NULL) {
	CALLOC(dad->upvects, 1, double*, ON, AMSC);
	dad->upvects[0] = &(sys->q[*qindex]);
	for(i=0; i < dad->upnumeles[0]; i++) {
	  (dad->sngs[i])->index[0] = (*qindex)++;
	}
      }
      if(dad->numfieldpts != 0) {
	dad->eval = &(sys->p[*pindex]);
	for(i=0; i < dad->numfieldpts; i++) {
	  (dad->fpts[i])->index = (*pindex)++;
	}
      }
      dad->index = (*pcindex)++;
    }
    else {
      for(i=0; i < dad->numkids; i++) {
	indexkid(sys, dad->kids[i], qindex, pindex, pcindex);
      }
    }
  }
}
/* 
  getPairs records the field point / singularity pairs that form
  the diagonal in the linear system.  It depends on the input field
  points and singularities being in this order.  In principal this 
  information is not needed, but in practice building effective 
  preconditioning matrices for desingularized problems depends on it.
  */
getPairs(snglist, fptlist)
snglrty *snglist;
fieldpt *fptlist;
{
  snglrty *nq;
  fieldpt *fp;
  for( nq=snglist, fp=fptlist; nq!=NULL, fp!=NULL; nq=nq->next, fp=fp->next) { 
    nq->myfpt = fp;
    fp->mysng = nq;
  }
}


/*
 * setExactUp marks as exact those cubes containing fewer than numterms
 * number of sngularities.  If the number of sngularities in the kids is 
 * less than numterms, the box is marked as exact and the snglrtys are copied 
 * up.  Otherwise, the number of nonzero kids is counted and put in upnumvects 
 * as usual.  
 */
setExactUp(sys, numterms)
ssystem *sys;
int numterms;
{
  int i, j, k, l, m, n;
  int side = sys->side;
  int depth = sys->depth;
  int realNumterms, numsngs, first;
  cube *nc, *nkid, *****cubes = sys->cubes;
  int allexact;

  /* We have passed numterms, but we may be forcing non-adaptive by passing 0,
     so we need to be sure we really know the number of terms for setting
     the multisize below. */
  realNumterms = multerms(sys->order); 
  for(i=depth; i > 0; i--, side /= 2) {
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
          nc = cubes[i][j][k][l];
 	  if((nc != NULL) && (nc->upcube == TRUE)) {
	    if(i == depth) {
/*cdel	      printf("exact up cube: %d %d %d %d \n",i,j,k,l);*/
	      ASSERT(nc->upnumvects != 0);
	      if(nc->upnumeles[0] <= numterms) {
		nc->exactUp = TRUE;
		nc->multisize = nc->upnumeles[0];
	      }
	      else {
		nc->exactUp = FALSE;
		nc->multisize = realNumterms;
	      }
	    }
	    else {
	      /* Count the number of non-empty kids and their sngularities. */
	      for(allexact=TRUE, m=0, numsngs=0, nc->upnumvects=0; 
		  m < nc->numkids; m++) {
		nkid = nc->kids[m];
		if((nkid != NULL) && (nkid->upcube == TRUE)) {
		  nc->upnumvects += 1;
		  if(nkid->exactUp == FALSE) allexact = FALSE;
		  else numsngs += nkid->upnumeles[0];
		}
	      }
	      /* If all nonempty kids exact, and # sngs <= # terms, mark 
		 this cube exact too, copy sngs, and promote pointers to 
		 strength vector.  Note this EXPLOITS special ordering of 
		 the strength vector. */
	      if((allexact == FALSE) || (numsngs > numterms)) { 
		nc->exactUp = FALSE;
		nc->multisize = realNumterms;
	      }
/*cftk this test is redundant.*/
	      else if((allexact == TRUE) && (numsngs <= numterms)) { 
		nc->exactUp = TRUE;
		nc->upnumvects = 1;
		CALLOC(nc->upvects, 1, double*, ON, AMSC);
		CALLOC(nc->upnumeles, 1, int, ON, AMSC);
		nc->upnumeles[0] = numsngs;
		nc->multisize = numsngs;
		CALLOC(nc->sngs, numsngs, snglrty*, ON, AMSC);
		for(m=0, first=TRUE, numsngs=0; m < nc->numkids; m++) {
		  nkid = nc->kids[m]; 
		  if((nkid != NULL) && (nkid->upcube == TRUE)) {
		    if(first == TRUE) {
		      nc->upvects[0] = nkid->upvects[0];
		      first = FALSE;
		    }
		    for(n=0; n < nkid->upnumeles[0]; n++) {
		      nc->sngs[numsngs++] = nkid->sngs[n];
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

/* 
  setExactDown marks as exact those cubes containing fewer than numterms
  number of field points.  If the number of field points in the kids is less 
  than numterms, the box is marked as exact and the vector of local expansion 
  coefficients is set to the field vector.  
*/
setExactDown(sys, numterms)
ssystem *sys;
int numterms;
{
  int i, j, k, l, m, n;
  int side = sys->side;
  int depth = sys->depth;
  int realNumterms, numfpts, first;
  cube *nc, *nkid, *****cubes = sys->cubes;
  int allexact;

  /* We have passed numterms, but we may be forcing non-adaptive by passing 0,
     so we need to be sure we really know the number of terms for setting
     the multisize below. */
  realNumterms = multerms(sys->order); 
  for(i=depth; i > 0; i--, side /= 2) {
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
          nc = cubes[i][j][k][l];
 	  if((nc != NULL) && (nc->downcube == TRUE)) {
	    if(i == depth) {
	      if(nc->numfieldpts <= numterms) {
		nc->exactDown = TRUE;
		nc->localsize = nc->numfieldpts;
	      }
	      else {
		nc->exactDown = FALSE;
		nc->localsize = realNumterms;
	      }
	    }
	    else {
	      /* Count the number of field points in nonempty kids. */
	      for(allexact=TRUE, m=0, numfpts=0; m < nc->numkids; m++) {
		nkid = nc->kids[m];
		if((nkid != NULL) && (nkid->downcube == TRUE)) {
		  if(nkid->exactDown == FALSE) allexact = FALSE;
		  else numfpts += nkid->numfieldpts;
		}
	      }
	      /* If all nonempty kids exact, # sngs <= # terms, mark exact, 
		 copy sngs, and promote pointers to snglrty and potential.  
		 Note EXPLOITS special ordering of the pot and snglrty 
		 vectors. */
	      if((allexact == FALSE) || (numfpts > numterms)) { 
		nc->exactDown = FALSE;
		nc->localsize = realNumterms;
	      }
	      else if((allexact == TRUE) && (numfpts <= numterms)) { 
		nc->exactDown = TRUE;
		CALLOC(nc->fpts, numfpts, fieldpt*, ON, AMSC);
		nc->localsize = numfpts;
		nc->numfieldpts = numfpts;
		for(m=0, first=TRUE, numfpts=0; m < nc->numkids; m++) {
		  nkid = nc->kids[m]; 
		  if((nkid != NULL) && (nkid->downcube == TRUE)) {
		    if(first == TRUE) {
		      nc->local = nkid->local;
		      first = FALSE;
		    }
		    for(n=0; n < nkid->numfieldpts; n++) {
		      nc->fpts[numfpts++] = nkid->fpts[n];
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

/*
* getnbrs finds all the nearest neighbors to downward-pass cubes which 
* contain singularities.  These may include the subject cube itself, as 
* well as collegues which are the descendants of an exactUp ancestor.  
* We require these nbrs-type nearest neighbors at all levels becuase we 
* use them to figure out the interaction lists and these exist at all 
* levels. 
*/
static getnbrs(sys)
ssystem *sys;
{
  cube *nc, *np, *****cubes = sys->cubes;
  int depth = sys->depth;
  int i, j, k, l, m, n, p, side, es;
  int numnbrs;
  int flag;
  /* Return if depth = 0, no neighbors. */
  if(depth == 0) return;
  
  /* Every level, get the nearest nbrs & nbrs due to parents being exact. */
  for(i = 1, side = 2; i <= depth; i++, side *= 2) {
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
	  nc = cubes[i][j][k][l];
/*cftk we ask is nc not NULL, we loop over cubes this way a lot, it must 
  stink on pipeline arch.  Why don't we make a linked list of cubes at 
  each level?*/
	  if((nc != NULL) && (nc->downcube==TRUE)) {
	    /* Find sidelength of exact cube. We want to know if the parent 
	       is an exact up-cube.*/
	    flag = FALSE;
/*	    for(es=1, np=nc->parent; np->exactUp==TRUE; np = np->parent,es *= 2);*/
	    for(es=1, np=nc->parent; np != NULL; np = np->parent,es *= 2) {
	      if(np->exactUp == TRUE) flag = TRUE;
	      if((np->exactUp == FALSE) && (flag == TRUE)){
		es = es/2;
		break;
	      }
	    }
	    if (flag == FALSE) es=1;
	    /* Stack up the nearest nbrs plus nbrs in exact cube. */
	    numnbrs = 0;
	    for(m = MIN((j-NNBRS), es * (j/es));
		m < MAX((j+NNBRS+1), es * (1 + (j / es))); m++) {
	      for(n = MIN((k-NNBRS), es * (k/es));
		  n < MAX((k+NNBRS+1), es * (1 + (k/es))); n++) {
		for(p = MIN((l-NNBRS), es * (l/es));
		    p < MAX((l+NNBRS+1), es * (1+(l/es))); p++) {
		  if( (m >= 0) && (n >= 0) && (p >= 0)
		     && (m < side) && (n < side) && (p < side)
		     /* c2.0 The cube must exist in the upward pass tree and
			you may be your own nearest neighbor. */
		     && (cubes[i][m][n][p] != NULL) 
		     && (cubes[i][m][n][p]->upcube == TRUE) ) {
		       cstack[numnbrs++] = cubes[i][m][n][p];
		  }
		}
	      }
	    }
	    nc->numnbrs = numnbrs;
	    CALLOC(nc->nbrs, numnbrs, cube*, ON, AMSC);
	    /* c2.0 subject cube is being placed at start of nbrs list if
	       it is in the list at all (it is an upcube). */
/*cdel we get what we want without this test	    if( nc->upcube == TRUE) {*/
	      n = numnbrs-1;
	      for(m=numnbrs-1; m >= 0; m--) {
		if(cstack[m] == nc) {
		  nc->nbrs[0] = cstack[m];
		}
		else {
		  nc->nbrs[n] = cstack[m];
		  n--;
		}
	      }
	    /* }*/
            /* Otherwise we just copy out the list.
	    else {
	      for(m=numnbrs; m >= 0; m--) {
		nc->nbrs[m] = cstack[m];
	      }
	    }*/
	  }
	}
      }
    }
  }
}
/*
* getpnbrs finds all the nearest neighbors to lowest level upward-pass 
* cubes which contain field points so that we can construct complete 
* local problems to invert for the precond mats.
*/
static getpnbrs(sys)
ssystem *sys;
{
  cube *nc, *****cubes = sys->cubes;
  int depth = sys->depth;
  int i, j, k, l, m, n, p, side;
  int numnbrs, warn_count=0;
  
  /* At the lowest level only, get the nearest neighbors. */
  side = 1;
/*cftk diag, need int pow*/
  for ( i=1; i<=depth; i++ ) {
    side *= 2; 
  }

  for(nc=sys->precondlist; nc != NULL; nc = nc->pnext) {
    j = nc->j;
    k = nc->k;
    l = nc->l;
    numnbrs = 0;
    for(m = MAX((j-1),0); m < MIN((j+2),side); m++) {
      for(n = MAX((k-1),0); n < MIN((k+2), side); n++) {
        for(p = MAX((l-1),0); p < MIN((l+2), side); p++) {
	  if( (cubes[depth][m][n][p] != NULL) && (cubes[depth][m][n][p]->downcube == TRUE) ) {
	    cstack[numnbrs++] = cubes[depth][m][n][p];
	  }
	}
      }
    }
    nc->numpnbrs = numnbrs;
/*cftk cdel probably don't need this test.*/
    if( numnbrs != 0 ) {
      CALLOC(nc->pnbrs, numnbrs, cube*, ON, AMSC);
      /* c2.0 subject cube is being placed at start of pnbrs list 
	 if it is in the list at all (it is a down cube). */
/*      if( nc->downcube == TRUE) {*/
	n = numnbrs-1;
	for(m=numnbrs-1; m >= 0; m--) {
	  if(cstack[m] == nc) {
	    nc->pnbrs[0] = cstack[m];
	  }
	  else {
	    nc->pnbrs[n] = cstack[m];
	    n--;
	  }
	}
	/*     }*/
      /* Otherwise we just copy out the list. */
/*cdel      else {
	for(m=numnbrs; m >= 0; m--) {
	  nc->pnbrs[m] = cstack[m];
	}
      }*/
    }
#if NOWARN == OFF
    else {
      if( warn_count++ < 10) printf("FLW-getpnbrs:  Found precondlist cube with no pnbrs. Index: %d %d %d %d\n",i,j,k,l);
    }
#endif
  }
}

/*
* getqnbrs finds all the nearest neighbors to lowest level upward-pass 
* cubes which contain singularities so that we can construct complete 
* local problems to invert for the precond mats.
*/
static getqnbrs(sys)
ssystem *sys;
{
  cube *nc, *****cubes = sys->cubes;
  int depth = sys->depth;
  int i, j, k, l, m, n, p, side;
  int numnbrs, warn_count=0;
  
  /* At the lowest level only, get the nearest neighbors. */
  side = 1;
/*cftk diag, need int pow*/
  for ( i=1; i<=depth; i++ ) {
    side *= 2; 
  }

  for(nc=sys->precondlist; nc != NULL; nc = nc->pnext) {
    j = nc->j;
    k = nc->k;
    l = nc->l;
    numnbrs = 0;
    for(m = MAX((j-1),0); m < MIN((j+2),side); m++) {
      for(n = MAX((k-1),0); n < MIN((k+2), side); n++) {
        for(p = MAX((l-1),0); p < MIN((l+2), side); p++) {
	  if( (cubes[depth][m][n][p] != NULL) && (cubes[depth][m][n][p]->upcube == TRUE) ) {
	    cstack[numnbrs++] = cubes[depth][m][n][p];
	  }
	}
      }
    }
    nc->numqnbrs = numnbrs;
/*cftk cdel probably don't need this test.*/
    if( numnbrs != 0 ) {
      CALLOC(nc->qnbrs, numnbrs, cube*, ON, AMSC);
      /* c2.0 subject cube is being placed at start of pnbrs list 
	 if it is in the list at all (it is a down cube). */
      /* cdel    if( nc->upcube == TRUE) {*/
	n = numnbrs-1;
	for(m=numnbrs-1; m >= 0; m--) {
	  if(cstack[m] == nc) {
	    nc->qnbrs[0] = cstack[m];
	  }
	  else {
	    nc->qnbrs[n] = cstack[m];
	    n--;
	  }
	}
      /*  } */
      /* Otherwise we just copy out the list. */
/*cdel      else {
	for(m=numnbrs; m >= 0; m--) {
	  nc->qnbrs[m] = cstack[m];
	}
      }*/
    }
#if NOWARN == OFF
    else {
      if( warn_count++ < 10) printf("FLW-getqnbrs:  Found precondlist cube with no qnbrs. Index: %d %d %d %d\n",i,j,k,l);
    }
#endif
  }
}

/*
 * cntDwnwdChg returns number of snglrtys in lowest level cubes contained 
 *  in "cube."
 */
int cntDwnwdChg(cp, depth)
int depth;			/* number of lowest level */
cube *cp;
{
  int i;
  int cnt;
  cube *kidc;

  if(cp->level == depth) return(cp->upnumeles[0]);
  else for(i = 0; i < cp->numkids; i++) 
      cnt += cntDwnwdChg(cp->kids[i], depth);
  return(cnt);
}

/*
 * linkcubes sets up the links between cubes requiring multipole expansions, 
 * local expansions, direct evaluations and other evaluations at field 
 * points, and inclusion in the preconditioner.
 * Note, upnumvects and exactUp/exactDown must be set!!!
 */
static linkcubes(sys)
ssystem *sys;
{
  cube *nc, **plnc, **pdnc, **pmnc, **ppnc, **penc, *****cubes = sys->cubes;
  int i, j, k, l, cnt = 0;
  int side, depth=sys->depth, numterms=multerms(sys->order);

  /* Allocate the vector of heads of cubelists. */
  CALLOC(sys->multilist, sys->depth+1, cube*, ON, AMSC);
  CALLOC(sys->locallist, sys->depth+1, cube*, ON, AMSC);

  pdnc = &(sys->directlist);
#if PRECOND != NONE
  ppnc = &(sys->precondlist);
#endif
  penc = &(sys->evallist);
  for(i=0, side = 1; i <= sys->depth; i++, side *= 2) {
    pmnc = &(sys->multilist[i]);
    plnc = &(sys->locallist[i]);
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
	  nc = cubes[i][j][k][l];
	  if(nc != NULL) {
	    /* Do multi exp. if cube not exact and in upward pass tree. */
	    if((i > 1) && (nc->upcube == TRUE) && (nc->exactUp == FALSE)) {
	      CALLOC(nc->multi, numterms, double, ON, AMSC);
	      *pmnc = nc;
	      pmnc = &(nc->mnext);
	    }
	    /* Do local exp. if cube not exact and in downward pass tree. */
	    if((i > 1) && (nc->downcube == TRUE) && (nc->exactDown == FALSE)) {
	      CALLOC(nc->local, numterms, double, ON, AMSC);
	      *plnc = nc;
	      plnc = &(nc->lnext);
	    }
	    /* c2.0 Add to evaluation list if a downward-pass finest-level
	       cube, add to direct list if a downward-pass finest-level
	       cube WITH neighbors (nbrs), and add to the precond list if an
               upward-pass finest-level cube with neighbors with 
               field points (pnbrs). */ 
	    if( i == depth ) {
	      if( nc->downcube == TRUE) {
	        *penc = nc;    /* eval list. */
		penc = &(nc->enext);
		if (nc->numnbrs > 0) { 
		  *pdnc = nc;  /* direct list. */
		  pdnc = &(nc->dnext);
		}
	      }
#if PRECOND != NONE
	      if( nc->upcube == TRUE ) {
		  *ppnc = nc;    /* precond list. */
		  ppnc = &(nc->pnext);
	      }
#endif
	    }
	  }
	}
      }
    }
  }
}

/*
 * setMaxq Determines maximum number of singularities OR fieldpoints (2.0) c
 * ontained in a single cube.  This is needed for allocation of expansion 
 * quantities.
 */
static setMaxq(sys)
ssystem *sys;
{
  int i, j, k, l, side;
  int maxq=0, maxlq=0;
  cube *nc, *****cubes = sys->cubes;

  for(i = 1, side = 2; i <= sys->depth; i++, side *= 2) {
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
	  nc = cubes[i][j][k][l];
	  if(nc != NULL) {
	    if((nc->upcube == TRUE) && (nc->upnumeles != NULL)) {
	      maxq = MAX(maxq, nc->upnumeles[0]);
	    }
	    if(nc->downcube == TRUE) maxq = MAX(maxq, nc->numfieldpts);
	  }
	}
      }
    }
  }
  sys->maxq = maxq;
/*cftk I think maxlq is vestigial. Note that we exit with it set to 
zero here.  In the call to mulMultiAlloc, MAX(maxlq, maxq), will always
return maxq.*/
  sys->maxlq = maxlq;
}

/* 
 * markup sets the flag to "flag" in the child and its nearest nbrs
 */
static markUp(child, flag)
cube *child;
int flag;
{
  int i,j;
  cube *nc, *np;

/* c2.0 child need not be in the neighbors list, but if it is an upcube then
   it is
  child->flag = flag;*/
  for(i = 0; i < child->numnbrs; i++) {
    child->nbrs[i]->flag = flag;
  }
}

/* 
 * getInter forms the true interaction list (see also comment at mulMatEval())
 * for cube "child," excluding interaction cubes which are not in the upward 
 * pass tree.  The call to this routine has the if which selects only downward
 * pass cubes for "child."
 * The interaction list pointer is saved in the interList cube struct field.
 */
static getInter(child)
cube *child;
{
  int i, j, vects, usekids, lc, jc, kc, ln, jn, kn;
  int numnbr = (child->parent)->numnbrs; /* number of neighbors */
  cube **nbrc = (child->parent)->nbrs; /* list of neighbor pointers */
  cube *sib;			/* pointer to sibling (same level as child) */
  cube **pstack = &(cstack[0]); /* temporary storage pointer */

  /* Mark the neighbors of the child cube (possibly including the child),
     i.e. cube->flag = TRUE. */
  markUp(child, TRUE);

  /* Unmarked children of child's parent's neighbors become the ilist */
  for(i = 0; i < numnbr; i++) { /* loop on neighbors */
    /* Check nbr's kids for a marked kid. */
    for(usekids = FALSE, j = 0; j < nbrc[i]->numkids; j++) { 
      sib = (nbrc[i]->kids)[j];
      if((sib != NULL) && (sib->flag == TRUE)) { 
	usekids = TRUE; 
	break; 
      }
    }
    /* Use nbr if no kids marked. */
    /* ...and it's really not a 1st nrst nbr of the parent 
       - this stops parent-sized cubes from getting into the ilist
         when they have empty child-sized cubes that are 2nd or 1st
	 nrst nbrs of the child cube 
       - should work with NNBRS = 1 (never allows parent-sized in list)
         and NNBRS > 2 (but cannot allow greater than parent-sized)
       (29May90) */
#if ON == ON
    lc = (child->parent)->l;
    jc = (child->parent)->j;
    kc = (child->parent)->k;
    ln = nbrc[i]->l;
    jn = nbrc[i]->j;
    kn = nbrc[i]->k;
    /*c2.0 The parent's neighbor cube must be in the upward-pass tree now.*/
    /* cftk the test on upcube is not necessary, all neighbors are upcubes.
       So I changed the test to an assert and it never fails.
       Also, this next test can be simplified.*/ 
    ASSERT(nbrc[i]->upcube = TRUE);
    if((RADINTER == ON) && (usekids == FALSE) &&
       ((lc-1 != ln && lc+1 != ln && lc != ln)
       || (jc-1 != jn && jc+1 != jn && jc != jn)
       || (kc-1 != kn && kc+1 != kn && kc != kn))) {  
      *pstack = nbrc[i];
      pstack++;
    }
#else				/* USE THIS PART FOR TESTING ONLY */
    if(RADINTER && (usekids == FALSE)) { /* PRODUCES INCORRECT ILISTS!!! */
      *pstack = nbrc[i];
      pstack++;
    }
#endif
    else for(j = 0; j < nbrc[i]->numkids; j++) { /* use nbr's kids. */
      sib = (nbrc[i]->kids)[j];	/* get sib of child cube of interest */
      /*c2.0 The sib must be in the upward-pass tree now.*/
      if((sib != NULL) && (sib->flag == FALSE) && (sib->upcube == TRUE)) { 
	*pstack = sib;
	pstack++;
      }
    }
  }

  /* clear all the flags */
  markUp(child, FALSE);

  /* Allocate and save the interaction list. */
  child->interSize = vects = pstack - &(cstack[0]);
  CALLOC(child->interList, vects, cube*, ON, AMSC);
  for(j = 0; j < vects; j++) child->interList[j] = cstack[j];

  return(vects);		/* return number of interaction elements */
}

/*
 * getAllInter generates explicit, true interaction lists for all 
 * non-empty cubes w/lev > 1.
 */
static getAllInter(sys)
ssystem *sys;
{
  int i, j, k, l, side, depth = sys->depth;
  cube *nc, *****cubes = sys->cubes;
  for(i = 2, side = 4; i <= depth; i++, side *= 2) {
    for(j=0; j < side; j++) {	/* loop through all cubes at levels > 1 */
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
	  nc = cubes[i][j][k][l];
	  /*c2.0 The interaction list relates the downward-pass tree to the 
	    upward-pass tree. */
	  if((nc != NULL) && (nc->downcube == TRUE)) getInter(nc);
	}
      }
    }
  }
}


/*
 * setTranslation sets the translationComputed flag in every cube. 
 */
setTranslation(sys, val)
ssystem *sys;
int val;
{
  cube *nc, *****cubes = sys->cubes;
  int i, j, k, l, side;
  for(i = sys->depth, side = sys->side; i >= 0; i--, side /= 2) {
    for(j=0; j < side; j++) {
      for(k=0; k < side; k++) {
	for(l=0; l < side; l++) {
          nc = cubes[i][j][k][l];
	  if(nc != NULL) nc->TranslationComputed = val;
	}
      }
    }
  }
}











