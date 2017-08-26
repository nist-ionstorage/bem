/*
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

/* # ***** sort to /src/header
   # ***** */
#include <malloc.h>
#include "resusage.h"

#define VERSION 1.9
#define RELEASE "18 July 1995"

extern cube *cstack[];		/* defined in mulGlobal.c */

/*********************************************************************** 
  macros for allocation with checks for NULL pntrs and 0 byte requests
  - also keep an allocated memory count
  - CALLOC() is used when the memory must be zeroed
    its core should be either calloc() or ualloc() (not as fast but
    more space efficient, no free list - uses sbrk() and never frees)
  - MALLOC() used when memory can be anything
    core should be malloc() or ualloc()
***********************************************************************/
/* #define CALCORE(NUM, TYPE) calloc((unsigned)(NUM),sizeof(TYPE)) */
#define CALCORE(NUM, TYPE) ualloc((unsigned)(NUM)*sizeof(TYPE))
/* #define MALCORE malloc */
#define MALCORE ualloc

extern long memcount;
extern long memQ2M;
extern long memQ2L;
extern long memQ2P;
extern long memL2L;
extern long memM2M;
extern long memM2L;
extern long memM2P;
extern long memL2P;
extern long memMSC;

#define AQ2M 0
#define AQ2L 1
#define AQ2P 2
#define AL2L 3
#define AM2M 4
#define AM2L 5
#define AM2P 6
#define AL2P 7
#define AMSC 8

#define DUMPALLOCSIZ                                                   \
{                                                                      \
  (void)fprintf(stderr,                                                \
		"Total Memory Allocated: %d kilobytes (brk = 0x%x)\n", \
		memcount/1000, sbrk(0));                               \
/* # ***** awked out for release */                                    \
  (void)fprintf(stderr, " Q2M  matrix memory allocated: %7.d kilobytes\n",\
		memQ2M/1000);                                          \
  memcount = memQ2M;                                                   \
  (void)fprintf(stderr, " Q2L  matrix memory allocated: %7.d kilobytes\n",\
		memQ2L/1000);                                          \
  memcount += memQ2L;                                                  \
  (void)fprintf(stderr, " Q2P  matrix memory allocated: %7.d kilobytes\n",\
		memQ2P/1000);                                          \
  memcount += memQ2P;                                                  \
  (void)fprintf(stderr, " L2L  matrix memory allocated: %7.d kilobytes\n",\
		memL2L/1000);                                          \
  memcount += memL2L;                                                  \
  (void)fprintf(stderr, " M2M  matrix memory allocated: %7.d kilobytes\n",\
		memM2M/1000);                                          \
  memcount += memM2M;                                                  \
  (void)fprintf(stderr, " M2L  matrix memory allocated: %7.d kilobytes\n",\
		memM2L/1000);                                          \
  memcount += memM2L;                                                  \
  (void)fprintf(stderr, " M2P  matrix memory allocated: %7.d kilobytes\n",\
		memM2P/1000);                                          \
  memcount += memM2P;                                                  \
  (void)fprintf(stderr, " L2P  matrix memory allocated: %7.d kilobytes\n",\
		memL2P/1000);                                          \
  memcount += memL2P;                                                  \
  (void)fprintf(stderr, " Miscellaneous mem. allocated: %7.d kilobytes\n",\
		memMSC/1000);                                          \
  memcount += memMSC;                                                  \
  (void)fprintf(stderr, " Total memory (check w/above): %7.d kilobytes\n",\
		memcount/1000);                                        \
/* # ***** awked out for release */                                    \
}

#define CALLOC(PNTR, NUM, TYPE, FLAG, MTYP)                                 \
{                                                                           \
     if((NUM)*sizeof(TYPE)==0);			                            \
     else if(((PNTR)=(TYPE*)CALCORE(NUM, TYPE))==NULL) {                    \
       (void)fprintf(stderr,                                                \
	 "FLE-%s: out of memory at line %d\n",                              \
	       __FILE__, __LINE__);                                         \
       (void)fprintf(stderr, " (NULL pointer on %d byte request)\n",        \
		     (NUM)*sizeof(TYPE));                                   \
       DUMPALLOCSIZ;                                                        \
       DUMPRSS;                                                             \
       (void)fflush(stderr);                                                \
       (void)fflush(stdout);                                                \
       if(FLAG == ON) exit(0);                                              \
     }                                                                      \
     else {                                                                 \
       memcount += ((NUM)*sizeof(TYPE));                                    \
       if(MTYP == AQ2M) memQ2M += ((NUM)*sizeof(TYPE));                     \
       else if(MTYP == AQ2L) memQ2L += ((NUM)*sizeof(TYPE));                \
       else if(MTYP == AQ2P) memQ2P += ((NUM)*sizeof(TYPE));                \
       else if(MTYP == AL2L) memL2L += ((NUM)*sizeof(TYPE));                \
       else if(MTYP == AM2M) memM2M += ((NUM)*sizeof(TYPE));                \
       else if(MTYP == AM2L) memM2L += ((NUM)*sizeof(TYPE));                \
       else if(MTYP == AM2P) memM2P += ((NUM)*sizeof(TYPE));                \
       else if(MTYP == AL2P) memL2P += ((NUM)*sizeof(TYPE));                \
       else if(MTYP == AMSC) memMSC += ((NUM)*sizeof(TYPE));                \
       else {                                                               \
         (void)fprintf(stderr,"FLE-CALLOC: unknown memory type %d\n", MTYP);\
         exit(0);                                                           \
       }                                                                    \
     }                                                                      \
}

#define MALLOC(PNTR, NUM, TYPE, FLAG, MTYP)                                  \
{                                                                            \
     if((NUM)*sizeof(TYPE)==0);			                             \
     else if(((PNTR)=(TYPE*)MALCORE((unsigned)((NUM)*sizeof(TYPE))))==NULL) {\
       (void)fprintf(stderr,                                                 \
	 "FLE-%s: out of memory at line %d\n",                               \
	       __FILE__, __LINE__);                                          \
       (void)fprintf(stderr, " (NULL pointer on %d byte request)\n",         \
		     (NUM)*sizeof(TYPE));                                    \
       DUMPALLOCSIZ;                                                         \
       DUMPRSS;                                                              \
       (void)fflush(stderr);                                                 \
       (void)fflush(stdout);                                                 \
       if(FLAG == ON) exit(0);                                               \
     }                                                                       \
     else {                                                                  \
       memcount += ((NUM)*sizeof(TYPE));                                    \
       if(MTYP == AQ2M) memQ2M += ((NUM)*sizeof(TYPE));                     \
       else if(MTYP == AQ2L) memQ2L += ((NUM)*sizeof(TYPE));                \
       else if(MTYP == AQ2P) memQ2P += ((NUM)*sizeof(TYPE));                \
       else if(MTYP == AL2L) memL2L += ((NUM)*sizeof(TYPE));                \
       else if(MTYP == AM2M) memM2M += ((NUM)*sizeof(TYPE));                \
       else if(MTYP == AM2L) memM2L += ((NUM)*sizeof(TYPE));                \
       else if(MTYP == AM2P) memM2P += ((NUM)*sizeof(TYPE));                \
       else if(MTYP == AL2P) memL2P += ((NUM)*sizeof(TYPE));                \
       else if(MTYP == AMSC) memMSC += ((NUM)*sizeof(TYPE));                \
       else {                                                               \
         (void)fprintf(stderr,"FLE-MALLOC: unknown memory type %d\n", MTYP);\
         exit(0);                                                           \
       }                                                                    \
     }                                                                      \
}

/*****************************************************************************

misc. global macros

*****************************************************************************/
#define NOT !
#define  ABORT()						      \
{   (void)fflush(stdout);					      \
    (void)fprintf(stderr, "FLE-%s: panic at line %d.\n",              \
	    __FILE__, __LINE__);				      \
    (void)fflush(stderr);					      \
    abort();							      \
}

#define ASSERT(condition) if(NOT(condition)) ABORT()

#define INNER(pap,p,ap,size) for(pap=0.0,i=1; i<=size; i++) pap += p[i]*ap[i];

#ifndef MAX
#define MAX(A,B)  ( (A) > (B) ? (A) : (B) )
#endif

#ifndef MIN
#define MIN(A,B)  ( (A) > (B) ? (B) : (A) )
#endif

#define ABS(A) ( ( (A) > 0 ) ? (A) : (-(A)) )

#define VCOPY(A, B) A[0] = B[0]; A[1] = B[1]; A[2] = B[2];

#define TRUE 1
#define FALSE 0

#define ON 1
#define OFF 0

#ifndef M_PI
/* pi constant included here since won't be in ANSI C */
#define M_PI       3.1415926535897931160E0  /*Hex  2^ 1 * 1.921FB54442D18 */
#endif

#define E_0 8.854187818E-12	/* epsilon0 +- .000000071E-12 F/m */

/* flags in chkList() in mulDisplay.c (chks direct, local or eval cube lsts) */
#define DIRECT 0
#define LOCAL 1
#define MULTIL 2
#define EVAL 3

/***********************************************************************
 
  Problem, element, configuration, and debug flags.

***********************************************************************/

/* Types of Problems. */
#define FIELD 0                 /* fastlap should only do a field computation.
                                   This is equivalent to computing a RHS.*/
#define GREEN 1                 /* fastlap should solve a linear system with
                                   both right- and left-hand side matrices
                                   as in a Green formulation.*/
#define INDIRECT 2              /* fastlap should solve a linear system with
                                   a left-hand side matrix and a right-hand
                                   side vector as in a single- or double-
                                   layer formulation.*/

/* Shapes of singularities. */
#define POINT 1                 /* Indicates snglrty geometry is a point. */
#define TRIANGLE 3              /* Indicates snglrty geometry is a triangle. */
#define QUADRILAT 4             /* Indicates snglrty geometry is a 
                                   quadrilateral. */

/* Types of singularities. */
#define POINT_SOURCE 1          /* Indicates snglrty is a point source. */
#define CONSTANT_SOURCE 11      /* Indicates snglrty is a constant source. */
#define CONSTANT_DIPOLE 12      /* Indicates snglrty is a constant dipole. */
#define LINEAR_SOURCE 21        /* Indicates snglrty is a linear source. */
#define LINEAR_DIPOLE 22        /* Indicates snglrty is a linear dipole. */

/* Types of iterative methods. */
#define GMRES 0			/* GMRES */

/* Types of preconditioners. */
#define NONE 0
#define OL 2			/* Overlapped preconditioner. Handles both 
				   square and non-square local problems. */
#define SP 3                    /* Forces square local problems by finding 
                                 each sing's' corresponding field point.*/

/* Multipole Configuration. */
#define DNTYPE GRENGD		/* type of downward/eval pass - see above */
#define MULTI ON        	/* ON=> add multipole contribution to P*q */
#define RADINTER ON	        /* ON=> Parent level multis in interlist */
#define NNBRS 2			/* Distance to consider a nearest nbr */
#define ADAPT ON		/* ON=> use adaptive algorithm */
#define OPCNT OFF		/* ON=> dump op-cnt of multiplies & exit */
#define DEFORD 2		/* default expansion order */
#define MAXORDER 6		/* Maximum expansion order */
#define MAXDEP 20		/* maximum partitioning depth */

/* Types of downward/eval passes. */
#define NOLOCL 0	       	/* multipoles evaluated directly, no locals */
#define NOSHFT 1		/* multis to locals w/o local2local shifts */
#define GRENGD 3		/* full Greengard downward pass/eval */

/* Linear System Solution Configuration. */
#define ITRTYP GMRES		/* type of iterative method (GMRES only) */
#define PRECOND OL 		/* Preconditioner, NONE, OL. */
#define PIVOT OFF               /* OFF=> do not pivot in G-J inversion. */
/* (add any new configuration flags to dumpConfig() in mulDisplay.c) */

/* Output Format Configuration. */
#define NOWARN OFF              /* ON=> suppress warning messages on stdout */
#define NOABRT ON               /* ON=> suppress some fatal error traps */
#define CMDDAT ON		/* ON=> dump command line info to output */
#define ITRDAT ON		/* ON=> dump residuals for every iteration */
#define TIMDAT ON		/* ON=> dump time and memory usage numbers */
#define CFGDAT ON		/* ON=> dump configuration to output */
#define MULDAT ON		/* ON=> dump brief multipole setup info */
/* Diagnostic Output. */
#define DISSYN OFF		/* ON=> display synopsis of cubes in lists */
#define DMTCNT OFF		/* ON=> display xform matrix counts by level */
#define PRECNT OFF		/* ON=> display precond mat count */

/* Display of transformation matrices. */
#define DISQ2M OFF		/* ON=> display Q2M matrices when built */
#define DISM2M OFF		/* ON=> display M2M matrices when built */
#define DISM2P OFF		/* ON=> display M2P matrices when built */
#define DISL2P OFF		/* ON=> display L2P matrices when built */
#define DISQ2L OFF		/* ON=> display Q2L matrices when built */
#define DISM2L OFF		/* ON=> display M2L matrices when built */
#define DISL2L OFF		/* ON=> display L2L matrices when built */
#define DALQ2M OFF		/* ON=> display all Q2M matrix build steps */
#define DALM2P OFF		/* ON=> display all M2P matrix build steps */
#define DALL2P OFF		/* ON=> display all L2P matrix build steps */
#define DALQ2L OFF		/* ON=> display all Q2L matrix build steps */

/* Display of other intermediate results. */
#define DUPVEC OFF		/* ON=> display lev 1 upward pass vectors */
#define DISFAC OFF		/* ON=> display factorial fractions in M2L */
#define DPSYSD OFF		/* ON=> display system after direct build */
#define DILIST OFF		/* ON=> display interaction lists */
/* misc debug */
#define CKDLST OFF		/* ON=> check direct list, prnt msg if bad */



