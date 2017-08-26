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


struct snglrty {		/* Distribution or point singularity */
  struct snglrty *next;		/* Next singularity in linked list. */
  struct fieldpt *myfpt;        /* Pointer to corresponding field point.*/
  double corner[4][3];		/* Corner point coordinates. */
  int shape;                    /* 4=quad panel, 3=tri panel, 1=point. */
  int index[4];			/* Singularity value index in q array */
  double X[3], Y[3], Z[3];	/* Def coord system, Z is normal direction. */
  double max_diag;		/* Longest diagonal of panel. */
  double min_diag;		/* Shortest diagonal. */
  double length[4];		/* Edge lengths. */
  double area;			/* Area of the singularity. */
  double x, y, z;		/* Centroid of the singularity.  */
  double moments[16];		/* Moments of the panel. */
  double *multipoleR;		/* Stores Q2M for RHS. */
  double *multipoleL;		/* Stores Q2M for LHS. */
  int lhsType;     		/* Indicates which sng type is required LHS.*/
  int rhsType;     		/* Indicates which sng type is required RHS.*/
  int activeType;	        /* Indicates which sng type is required now.*/
  int cond;                     /* Sub-surface number. */
  int order;			/* Multipole order. */
};

typedef struct snglrty snglrty;

struct fieldpt {
  struct fieldpt *next;		/* Next fieldpt in linked list. */
  struct snglrty *mysng;        /* Pointer to corresponding singularity. */
  int index;			/* Field value index in p array */
  double x, y, z;               /* Position of field point. */
  int deriv;                    /* jt: Indicates derivative or potential */
  double nrm[3];                /* jt: direction of derivative */
};

typedef struct fieldpt fieldpt;

struct cube {		
/* Definition variables. */
  int index;			/* unique index. */
  int upcube;                   /* TRUE == cube in upward pass tree. */
  int downcube;                 /* TRUE == cube in downward pass tree. */
  int level;			/* 0 => root. */
  double x, y, z;		/* Position of cube center. */
  int j, k, l;			/* cube is cubes[level][j][k][l]. */
  int exactUp;                  /* Used to indicate that no multipole used. */
  int exactDown;                /* Used to indicate that no local used. */
  int flag;			/* used for marking for tree walks. */
  int TranslationComputed;	/* Indicates translation mats are computed. */
				/* Set to false if false in any of children. */

/* Upward Pass variables. */
  struct cube *mnext;		/* Ptr to next cube on which to do multi. */
  int upnumvects;		/* 0 if empty,  1 on bot level if not empty, 
				   else # nonempty kids. */ 
  int *upnumeles;		/* numeles[0] = # sngs on bot level, else
				   number of terms in kid's expansion. */
  double **upvects;             /* vects[0] = sngs on bot level, else vectors
                                   of kids' expansion terms. */
  int multisize;                /* Number of terms in the expansion. */
  double *multi;                /* Vector of multi coefficients. */
  double ***upmats;             /* Matrices for singularity to multi or multi 
                                   to multi. upmats[i] is multisize x 
                                   upnumeles[i]. */

/* Downward Pass variables. */
  struct cube *lnext;           /* Ptr to next cube on which to do local. */
  int downnumvects;		/* Number of cubes in iteraction list. */
  int *downnumeles;		/* # of eles in interact cube's expansion. */
  double **downvects;		/* Vects of interact cube's expansion. */

  int localsize;                /* Size of the local expansion */
  double *local;                /* Vector of local field coefs */
  double ***downmats;           /* Matrices for multi to field point, or multi
                                   to local or local to local.  downmats[i] is
                                   downnumele x localsize. */

  struct cube **interList;	/* explicit interaction list 
				   - for fake dwnwd passes and eval pass */
  int interSize;		/* number of elements in interList
				   - often != downnumvects nor evalnumvects */

  /* evaluation pass variables */
  struct cube *enext;		/* Pntr to next cube to evaluate */
  int evalnumvects;		/* For exact = #in inter list, o.w. = 1 */
  int *evalnumeles;		/* Num of elements in inter list entry exp */
  double **evalvects;		/* Multi, local, or sngs of ilist entry */

  int numfieldpts;              /* Number of fieldpoints, lowest level only.*/ 
  double *eval;			/* vector of field point fields in cube */
  double ***evalmats;		/* matrices for multi to potential, local to
				   potential or snglrty to potential */

/* Direct portion variables. */
  struct cube *dnext;		/* Ptr to next cube on which to do direct. */
  struct cube *pnext;		/* Ptr to next cube on which to do precond. */
/*  struct cube *rpnext;*/		/* Reverse ptr to next cube to do precond. */
  int directnumvects;		/* Number of vects, self plus nbrs. */
  int *directnumeles;		/* # of elements in the nbrs chg vect. 
				   directnumeles[0] = numsngs in cube. */
  double **directq;		/* Vecs of chg vecs, directq[0] this cube's. */
  double ***directmats;		/* Potential Coeffs in cube and neighbors. */
  double ***precondmats;	/* Precond Coeffs in cube and neighbors. */

/* Cube structure variables. */
  snglrty **sngs;               /* Array of singularity ptrs. Only used lowest 
                                   level. */
  fieldpt **fpts;               /* Array of field point ptrs. Only used lowest 
                                   level. */
  struct cube **nbrs;           /* Array of ptrs to neighbors with
                                   singularities for direct comps. */
  int numnbrs;                  /* Number of nbrs-type neighbors. */
  struct cube **pnbrs;          /* Array of ptrs to neighbors with 
                                   field pts for precond comps. */
  struct cube **qnbrs;          /* Array of ptrs to neighbors with 
                                   sing's for precond comps. */
  int numpnbrs;                 /* Number of pnbrs-type neighbors. */
  int numqnbrs;                 /* Number of qnbrs-type neighbors. */
  struct cube **kids;           /* Array of children ptrs. */
  int numkids;                  /* Number of kids. */
  struct cube *parent;          /* Ptr to parent cube. */

};
typedef struct cube cube;

struct ssystem {
  int side;                     /* # cubes per side on lowest level. */
  int depth;			/* # of levels of cubes. */
  int order;			/* # of levels of cubes. */
  double length;		/* Length per cube on lowest level. */
  double minx, miny, minz;	/* Coordinates of one corner of the domain. */
  int maxq;			/* Maximum number of sngs in a cube. */
  int maxlq;			/* Maximum # of sngs in lowest level cube. */
  double *q;                    /* The vector of lowest level snglrtys. */
  double *p;                    /* The vector of lowest level potentials. */
  cube *****cubes;		/* The array of cube pointers. */
  cube **multilist;             /* Array of ptrs to first cube in linked list
				   of cubes to do multi at each level. */
  cube **locallist;             /* Array of ptrs to first cube in linked list
				   of cubes to do local at each level. */
  cube *directlist;		/* Head of linked lst of leaf cubes w/ field
				   points and neighbors. */
  cube *evallist;		/* Head of linked lst of leaf cubes w/ f'p'ts. */
  cube *precondlist;		/* Head of linked lst of precond blks. */
};
typedef struct ssystem ssystem;
