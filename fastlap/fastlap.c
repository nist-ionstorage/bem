/*
* This software is being provided to you, the LICENSEE, by the Massachusetts
* Institute of Technology (M.I.T.) under the following license. By
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
* Copyright 1992 by the Massachusetts Institute of Technology. All rights
* reserved.
*
* THIS SOFTWARE IS PROVIDED "AS IS", AND M.I.T. MAKES NO REPRESENTATIONS OR
* WARRANTIES, EXPRESS OR IMPLIED. By way of example, but not limitation,
* M.I.T. MAKES NO REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS
* FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE OR
* DOCUMENTATION WILL NOT INFRINGE ANY THIRD PARTY PATENTS, COPYRIGHTS,
* TRADEMARKS OR OTHER RIGHTS.
*
* The name of the Massachusetts Institute of Technology or M.I.T. may NOT
* be used in advertising or publicity pertaining to distribution of the
* software. Title to copyright in this software and any associated
* documentation shall at all times remain with M.I.T., and USER agrees to
* preserve same.
*
* Written by: K. Nabors, T. Korsmeyer, and J. White
*
*/

/*
*   This is the top-level of fastlap. It has been constructed to be
*   linkable to FORTRAN or c. Therefore the incoming data is treated as
*   being passed by address, a la FORTRAN. Also, the panel attributes
*   are treated as FORTRAN style arrays and so step one (after some
*   checking of suitability of arguments) is to call 'loadCharge' which
*   takes these arrays and puts them into the singularity structure and
*   to call 'loadFieldpt' which puts the field (or collocation) points
*   into the fieldpoint structure.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#ifndef _TIME_
#include <time.h>
#endif

#include "fastlap.h"
#include "mulLocal.h"
#include "mulMulti.h"
#include "calcp.h"

#define VERTS 4
#define DIMEN 3

/*
*   If you are linking to FORTRAN, the name of this procedure must be
*   fastlap_(..., whereas for C, you must leave it as fastlap(...
*
*   See 'driver.f' for an example of how to call from FORTRAN. See
*   'driver.c' for an example of how to call from c.
*
*   In the makefile provided with the distribution, the above
*   substitution is made for you via sed if you choose to make
*   the version for linking to fortran.
*/

// Function declarations
void mulMatDirect(ssystem *sys);                            // In mulMats.c
void mulMatUp(ssystem *sys);                                // In mulMats.c
void mulMatEval(ssystem *sys);                              // In mulMats.c
void mulMatDown(ssystem *sys);                              // In mulMats.c
void mulDirect(ssystem *sys);                               // In mulDo.c
void mulPrecond(ssystem *sys, int size);                    // In mulDo.c
void spmulPrecond(ssystem *sys, double *work, int size);    // In mulDo.c
void mulUp(ssystem *sys);                                   // In mulDo.c
void mulEval(ssystem *sys);                                 // In mulDo.c
void mulDown(ssystem *sys);                                 // In mulDo.c
void setTranslation(ssystem *sys, int val);                 // In mulSetup.c    

int fastlap(plhsSize,prhsSize,pnumSing,px,pshape,pdtype,plhsType,prhsType,plhsIndex,prhsIndex,plhsVect,prhsVect,pxf,pxnrm,pnumLev,pnumMom,pmaxItr,ptol,pjob, pAreas)
int *plhsSize, *prhsSize, *pnumSing, *pshape, *pdtype, *plhsType, *prhsType, *plhsIndex, *prhsIndex, *pnumLev, *pnumMom, *pmaxItr, *pjob;
double *px, *plhsVect, *prhsVect, *pxf, *pxnrm, *ptol, *pAreas;
{
    int lhsSize, rhsSize, numSing, numLev, autlev, numMom, maxItr;
    snglrty *nq, *snglist, *loadSnglrty();
    fieldpt *nf, *fptlist, *loadFieldpt();
    ssystem *sys, *mulInit();
    double initalltime, ttlsetup, ttlsolve;
    extern size_t memcount;
    extern double prectime, conjtime, dirtime, multime, uptime, downtime;
    extern double structime, dirsetime, mulsetime;
    extern double evaltime;
    int ttliter, i;
    static int callnum = 0;

    if(callnum != 0) {
        restartCalcp();
        /*restartStore(); */ /* Was part of UglyAlloc */
        restartLocal();
        restartMulti();
    }
    callnum++;

    /* Convert passed addresses into ints (this is to maintain FORTRAN
         compatibility).*/
    lhsSize = *plhsSize;
    rhsSize = *prhsSize;
    numSing = *pnumSing;
    numMom = *pnumMom;
    numLev = *pnumLev;
    maxItr = *pmaxItr;

    /* Initialize memory and time counters, etc. */
    memcount = 0;
    prectime = conjtime = dirtime = multime = uptime = downtime = 0.0;
    structime = dirsetime = mulsetime = 0.0;
    evaltime = 0.0;

    autlev = OFF;

    if(numMom < 0 || numMom > MAXORDER) {
        fprintf(stderr,
            "FLE-fastlap: bad expansion order: %d\n", numMom
            );
        exit(0);
    }

    if(numLev < 1 || numLev > MAXDEP) {
#if NOWARN == OFF
        fprintf(stdout,
            "FLW-fastlap: automatic depth ON, selected was: %d\n", numLev);
#endif
        autlev = ON;
    }

#if CFGDAT == ON
    dumpConfig(stdout);
#endif
    starttimer;
    /* Put the passed singularity data into the singularity data structure. */
    snglist = loadSnglrty(numSing,px,pshape,plhsType,prhsType);

    /* Put the passed fieldpoint data into the fieldpoint data structure. */
    fptlist = loadFieldpt(lhsSize, pxf, pdtype, pxnrm);

    /* Get panel geometric attributes. */
    initcalcp(snglist, numMom, pAreas);

    /* Set up cubes with singularities and field points. */
    sys = mulInit(autlev, numLev, numMom, snglist, lhsSize, rhsSize, fptlist);

    if(autlev == ON) numLev = sys->depth;
    stoptimer;
    structime = dtime;

#if CMDDAT == ON
    fprintf(stdout, "\nFastLap ARG LIST SUMMARY:\n");
    fprintf(stdout, "   Number of singularities to process: %d\n", numSing);
    fprintf(stdout, "   Number of field or collocation points: %d\n", lhsSize);
    fprintf(stdout, "   Expansion order: %d\n", numMom);
    fprintf(stdout, "   Number of partitioning levels: %d\n", numLev);
    fprintf(stdout, "   Maximum allowed iterations: %d\n", maxItr);
    fprintf(stdout, "   Convergence tolerance: %g\n", *ptol);
#endif

#if MULDAT == ON
    dumpMulSet(sys, lhsSize, numLev, numMom);
#endif
    fflush(stdout);

    starttimer;
    mulMultiAlloc(MAX(sys->maxq,sys->maxlq), numMom, sys->depth);
    stoptimer;
    initalltime = dtime;        /* save initial allocation time */

    /* To understand the manipulations of sys->p and sys->q, you need to
         know that 'apply' (the top procedure of the multipole algorithm)
         multiplies the vector sys->q by the matrix implied by the singularity
         structure. The result of this multiplication is sys->p. That is:

         {sys->p} = [singularity influence matrix] {sys->q} */

    /* Set up mapping matrices. */
    setTranslation(sys, FALSE);

    matSetup(sys,numSing,snglist);

    /* Compute the field or the RHS for a direct problem. In this case,
     prhsVect points to a vector of singularity strengths and so is ordered
     like the singularities. */
    if(*pjob == FIELD || *pjob == GREEN) {

        /* Use the indices to the RHS vector to load sys->q.*/
        /*cftk this is strictly one strength per sing. at this point.*/
        nq = snglist;
        for(i=0;i<rhsSize;i++) {
            sys->q[nq->index[0]] = prhsVect[i];
            nq = nq->next;
        }

        apply(sys, sys->q, sys->p, lhsSize);
    }
    /* Or else the RHS matrix is the identity and this is an indirect problem.
     In this case, prhsVect points to a vector of field values and so is ordered
     like the fieldpoints. */

    else if(*pjob == INDIRECT) {
        if(rhsSize != lhsSize) {
            fprintf(stderr,
                "FLE-fastlap: Cannot do INDIRECT with unequal length \n right- and left-hand-side vectors. lhsSize = %d and rhsSize = %d\n", lhsSize, rhsSize);
            exit(0);
        }
        nf = fptlist;
        for(i=0;i<rhsSize;i++) {
            sys->p[nf->index] = prhsVect[i];
            nf = nf->next;
        }
    }
    else {
        fprintf(stderr,
            "FLE-fastlap: Asking for undefined fastlap job. %d\n",*pjob);
        exit(0);
    }

    /* If only a field computation is asked for, return to calling routine.*/
    if(*pjob == FIELD) {
        nf = fptlist;
        for(i=0;i<lhsSize;i++) {
            plhsVect[i] = sys->p[nf->index];
            nf = nf->next;
        }
        #if TIMDATFIELD == ON
            ttlsetup = structime + initalltime + dirsetime + mulsetime;
            multime = uptime + downtime + evaltime;
            ttlsolve = dirtime + multime + conjtime;
        
            fprintf(stdout, "\nTIME AND MEMORY USAGE SYNOPSIS\n");
            fprintf(stdout, "   Total time: %g\n", ttlsetup + ttlsolve);
            fprintf(stdout, "       Total setup time: %g\n", ttlsetup);
            fprintf(stdout, "           Data structure setup time: %g\n", structime);
            fprintf(stdout, "           Direct matrix setup time: %g\n", dirsetime);
            fprintf(stdout, "           Multipole matrix setup time: %g\n", mulsetime);
            fprintf(stdout, "           Initial misc. allocation time: %g\n", initalltime);
            fprintf(stdout, "       Total iterative P*q = psi solve time: %g\n", ttlsolve);
            fprintf(stdout, "           P*q product time, direct part: %g\n", dirtime);
            fprintf(stdout, "           Total P*q time, multipole part: %g\n", multime);
            fprintf(stdout, "               Upward pass time: %g\n", uptime);
            fprintf(stdout, "               Downward pass time: %g\n", downtime);
            fprintf(stdout, "               Evaluation pass time: %g\n", evaltime);
            fprintf(stdout, "           Preconditioner solution time: %g\n", prectime);
            fprintf(stdout, "           Iterative loop overhead time: %g\n", conjtime);
        
            fprintf(stdout, "   Total memory allocated: %d kilobytes\n", memcount/1000);
            uallocEfcy(memcount);
            fprintf(stdout, "       Q2M matrix memory allocated: %7.d kilobytes\n",
                memQ2M/1000);
            memcount = memQ2M;
            fprintf(stdout, "       Q2L matrix memory allocated: %7.d kilobytes\n",
                memQ2L/1000);
            memcount += memQ2L;
            fprintf(stdout, "       Q2P matrix memory allocated: %7.d kilobytes\n",
                memQ2P/1000);
            memcount += memQ2P;
            fprintf(stdout, "       L2L matrix memory allocated: %7.d kilobytes\n",
                memL2L/1000);
            memcount += memL2L;
            fprintf(stdout, "       M2M matrix memory allocated: %7.d kilobytes\n",
                memM2M/1000);
            memcount += memM2M;
            fprintf(stdout, "       M2L matrix memory allocated: %7.d kilobytes\n",
                memM2L/1000);
            memcount += memM2L;
            fprintf(stdout, "       M2P matrix memory allocated: %7.d kilobytes\n",
                memM2P/1000);
            memcount += memM2P;
            fprintf(stdout, "       L2P matrix memory allocated: %7.d kilobytes\n",
                memL2P/1000);
            memcount += memL2P;
            fprintf(stdout, "       Miscellaneous mem. allocated:%7.d kilobytes\n",
                memMSC/1000);
            memcount += memMSC;
            fprintf(stdout, "       Total memory (check w/above):%7.d kilobytes\n\n",
                memcount/1000);
        #endif
        return(0);
    }

    /* Note that matrices have been computed. */
    setTranslation(sys, TRUE);
/*cftk see note on use of size down in this routine. */
/*  matSetup(sys,size,snglist);*/
    matSetup(sys,numSing,snglist);

    /* Decide whether to use preconditioning, and which kind. */
#if NOABRT == OFF
    if((sys->depth <= 1) && (PRECOND != NONE)) {
        fprintf(stderr,"FLE-fastlap: You cannot precondition level %d calcs.", sys->depth);
        exit(0);
    }
#endif

    starttimer;
#if PRECOND == OL
    olmulMatPrecond(sys);
#endif
#if PRECOND == SP
    spmulMatPrecond(sys);
#endif
    stoptimer;
    prectime += dtime;

    /* Now we solve the system using the values of the potential and
    normal derivative of the potential found in the structure
    snglrty and picked out by knowing the boundary type. The solution
    is returned as updated values of the potential and
    normal derivative in the structure snglrty. */

    ttliter = ONsolve(sys, snglist, fptlist, lhsSize, maxItr, ptol);

    /*Unload the solution into the arg list vector.*/
    nq = snglist;
    for(i=0;i<lhsSize;i++) {
        plhsVect[i] = sys->q[nq->index[0]];
        nq = nq->next;
    }

#if TIMDAT == ON
    ttlsetup = structime + initalltime + dirsetime + mulsetime;
    multime = uptime + downtime + evaltime;
    ttlsolve = dirtime + multime + prectime + conjtime;

    fprintf(stdout, "\nTIME AND MEMORY USAGE SYNOPSIS\n");
#endif

#ifdef OTHER
    if(TIMDAT == ON) {
        fprintf(stdout,
            "FLW-fastlap: compilation with OTHER flag gives incorrect times\n");
    }
#endif

#if TIMDAT == ON
    fprintf(stdout, "   Total time: %g\n", ttlsetup + ttlsolve);
    fprintf(stdout, "       Total setup time: %g\n", ttlsetup);
    fprintf(stdout, "           Data structure setup time: %g\n", structime);
    fprintf(stdout, "           Direct matrix setup time: %g\n", dirsetime);
    fprintf(stdout, "           Multipole matrix setup time: %g\n", mulsetime);
    fprintf(stdout, "           Initial misc. allocation time: %g\n", initalltime);
    fprintf(stdout, "       Total iterative P*q = psi solve time: %g\n", ttlsolve);
    fprintf(stdout, "           P*q product time, direct part: %g\n", dirtime);
    fprintf(stdout, "           Total P*q time, multipole part: %g\n", multime);
    fprintf(stdout, "               Upward pass time: %g\n", uptime);
    fprintf(stdout, "               Downward pass time: %g\n", downtime);
    fprintf(stdout, "               Evaluation pass time: %g\n", evaltime);
    fprintf(stdout, "           Preconditioner solution time: %g\n", prectime);
    fprintf(stdout, "           Iterative loop overhead time: %g\n", conjtime);

    fprintf(stdout, "   Total memory allocated: %zd kilobytes\n", memcount/1000);
    /* uallocEfcy(memcount);*/
    fprintf(stdout, "       Q2M matrix memory allocated: %7.zd kilobytes\n",
        memQ2M/1000);
    memcount = memQ2M;
    fprintf(stdout, "       Q2L matrix memory allocated: %7.zd kilobytes\n",
        memQ2L/1000);
    memcount += memQ2L;
    fprintf(stdout, "       Q2P matrix memory allocated: %7.zd kilobytes\n",
        memQ2P/1000);
    memcount += memQ2P;
    fprintf(stdout, "       L2L matrix memory allocated: %7.zd kilobytes\n",
        memL2L/1000);
    memcount += memL2L;
    fprintf(stdout, "       M2M matrix memory allocated: %7.zd kilobytes\n",
        memM2M/1000);
    memcount += memM2M;
    fprintf(stdout, "       M2L matrix memory allocated: %7.zd kilobytes\n",
        memM2L/1000);
    memcount += memM2L;
    fprintf(stdout, "       M2P matrix memory allocated: %7.zd kilobytes\n",
        memM2P/1000);
    memcount += memM2P;
    fprintf(stdout, "       L2P matrix memory allocated: %7.zd kilobytes\n",
        memL2P/1000);
    memcount += memL2P;
    fprintf(stdout, "       Miscellaneous mem. allocated:%7.zd kilobytes\n",
        memMSC/1000);
    memcount += memMSC;
    fprintf(stdout, "       Total memory (check w/above):%7.zd kilobytes\n\n",
        memcount/1000);
#endif
    return(ttliter);
}

/*
 * snglrty returns a list of snglrty structs derived from passed data:
 * shape, vertices, and type.
 */
snglrty *loadSnglrty(int numSing, double *px, int *pshape, int *plhsType, int *prhsType)
{
    int i,j,k;
    snglrty *snglist, *nq=NULL;

    /* Initialize snglrty list head pointer. */
    snglist = NULL;

    for(i=0;i<numSing;i++) {

        /* Load structs from passed vectors. */

        if(pshape[i] == POINT) {

            /* Setup points in snglrty structs. */

            /* Allocate snglrty struct to fill in. */
            if(snglist == NULL) {
    CALLOC(snglist, 1, snglrty, ON, AMSC);
    nq = snglist;
            }
            else {
    CALLOC(nq->next, 1, snglrty, ON, AMSC);
    nq = nq->next;
            }

            /* Fill first vertex with the point coordinates. */
            for(j=0;j<1;j++) {
    for(k=0;k<3;k++) {
        (nq->corner[j])[k] = px[(i*VERTS*DIMEN)+(j*DIMEN)+k];
    }
            }

            /* Fill in descriptors. */
            nq->shape = pshape[i];
            nq->lhsType = plhsType[i];
            nq->rhsType = prhsType[i];
        }

        else if(pshape[i] == QUADRILAT) {

            /* Setup quads in snglrty structs. */

            /* Allocate snglrty struct to fill in. */
            if(snglist == NULL) {
    CALLOC(snglist, 1, snglrty, ON, AMSC);
    nq = snglist;
            }
            else {
    CALLOC(nq->next, 1, snglrty, ON, AMSC);
    nq = nq->next;
            }

            /* Fill in corners. */
            for(j=0;j<4;j++) {
    for(k=0;k<3;k++) {
        (nq->corner[j])[k] = px[(i*VERTS*DIMEN)+(j*DIMEN)+k];
    }
            }

            /* Fill in descriptors. */
            nq->shape = pshape[i];
            nq->lhsType = plhsType[i];
            nq->rhsType = prhsType[i];
        }

        else if(pshape[i] == TRIANGLE) {

            /* Setup tris in snglrty structs. */

            /* Allocate snglrty struct to fill in. */
            if(snglist == NULL) {
    CALLOC(snglist, 1, snglrty, ON, AMSC);
    nq = snglist;
            }
            else {
    CALLOC(nq->next, 1, snglrty, ON, AMSC);
    nq = nq->next;
            }

            /* Fill in corners. */
            for(j=0;j<3;j++) {
    for(k=0;k<3;k++) {
        (nq->corner[j])[k] = px[(i*VERTS*DIMEN)+(j*DIMEN)+k];
    }
            }

            /* Fill in descriptors. */
            nq->shape = pshape[i];
            nq->lhsType = plhsType[i];
            nq->rhsType = prhsType[i];
        }
    }
    return(snglist);
}

/*
 * fieldpt returns list of fieldpoint structs derived from passed data:
 * 3-space coordinates.
 */
fieldpt *loadFieldpt(int size, double *px, int *dtype, double *xnrm)
{
    int i;
    fieldpt *fptlist, *nf=NULL;

    /* initialize fieldpoint list head pointer */
    fptlist = NULL;

    for(i=0;i<size;i++) {
        /* load structs from passed vectors. */
        /* allocate snglrty struct to fill in */
        if(fptlist == NULL) {
            CALLOC(fptlist, 1, fieldpt, ON, AMSC);
    nf = fptlist;
        }
        else {
            CALLOC(nf->next, 1, fieldpt, ON, AMSC);
            nf = nf->next;
        }

        /* Fill in coordinates */
        /* cftk this is not x[3] because the panel centroid is done this
             way and I want to duplicate that. Maybe at some point I will
             see why this was considered a good idea.*/
        nf->x = px[(i*DIMEN)];
        nf->y = px[(i*DIMEN+1)];
        nf->z = px[(i*DIMEN+2)];
        /* jt: modifications for derivatives */
        nf->deriv = dtype[i];
        if ( nf->deriv == TRUE ) {   /* We want the normal derivative at the fieldpoint. */
            nf->nrm[0] = xnrm[i*DIMEN];
            nf->nrm[1] = xnrm[i*DIMEN+1];
            nf->nrm[2] = xnrm[i*DIMEN+2];
        }
    }
    return(fptlist);
}

/*
 * matSetup sets up the mapping matrices for the cube data structure.
 * Note, the TranslationComputed flags in each cube should have
 * already been determined.
 */
void matSetup(ssystem *sys, int size, snglrty *snglist)
{
    extern double dirsetime, mulsetime;

    starttimer;
    mulMatDirect(sys);      /* Compute the direct part matrices. */
    stoptimer;
    dirsetime += dtime;

#if DPSYSD == ON
    dissys(sys);
#endif

#if CKDLST == ON
/*  chkList(sys, DIRECT);*/
#endif
/*cftk this routine dumps coefficient info to stdout. It is really only
appropriate to constant strength planar panels, maybe we should junk it.
    dumpnums(ON, size); */          /* save num/type of pot. coeff calculations */

    starttimer;
    mulMatUp(sys);      /* Compute the upward pass matrices. */

#if DNTYPE == NOSHFT
    mulMatDown(sys);        /* find matrices for no L2L shift dwnwd pass */
#endif

#if DNTYPE == GRENGD
    mulMatDown(sys);        /* find matrices for full Greengard dnwd pass*/
#endif

    mulMatEval(sys);        /* set up matrices for evaluation pass */

    stoptimer;
    mulsetime += dtime;     /* save multipole matrix setup time */

/*cftk see call to this above
    dumpnums(OFF, size);*/  /* dump num/type of pot. coeff calculations */

#if CKDLST == ON
    chkList(sys, DIRECT);
    chkList(sys, EVAL);
    chkLowLev(sys, DIRECT);
    dumpList(sys,DIRECT);
    dumpList(sys,EVAL);
    dumpList(sys,LOCAL);
    dumpList(sys,MULTIL);
    exit(0);
#endif

#if DISSYN == ON
    dumpSynop(sys);
#endif

#if DMTCNT == ON
    dumpMatBldCnts(sys);
#endif

}

/*
 * ONsolve solves the linear system.
 */

int ONsolve(ssystem *sys, snglrty *snglist, fieldpt *fptlist, int size, int maxiter, double *tol)
{
    int i, iter = 0;
    double *p, *r, *ap, *z;
    double **bp, **bap;

    /* Allocate space for vectors , r=residual and p=projection, ap = Ap.
         We use vector[1:size+1] in gmres.*/
    CALLOC(r, size+1, double, ON, AMSC);
    CALLOC(z, size+1, double, ON, AMSC);

    /* allocate for accumulated basis vectors.*/
    CALLOC(bp, maxiter+1, double*, ON, AMSC);
    CALLOC(bap, maxiter+1, double*, ON, AMSC);

    /* P is the "pseudo-snglrty" for multipole. Ap is the "pseudo-potential". */
    ASSERT(sys->q != NULL);
    ASSERT(sys->p != NULL);
    p = sys->q;
    ap = sys->p;

    /* Set up the initial residue vector and snglrty guess. */
    for(i=1; i <= size; i++) { r[i] = sys->p[i]; }

    iter = gmres(sys,snglist,fptlist,p,r,ap,z,bp,bap,size,maxiter,tol);
    fflush(stdout);
    return(iter);
}

/*
 * gmres is a preconditioned (possibly) implementation of Saad & Schultz.
 */
int gmres(ssystem *sys, snglrty *snglist, fieldpt *fptlist, double *p, double *r,
          double *ap, double *z, double **bv, double **bh, int size, int maxiter, double *tol)
{
    snglrty *nq;
    fieldpt *nf;
    int iter, i, j;
    double r0, rnorm, apnorm;
    double hi, hip1, length;
    extern double prectime, conjtime;
    double *c=NULL, *s=NULL, *g=NULL, *y=NULL;
#if PRECOND == SP
    double *work=NULL;
#endif
    starttimer;

#if PRECOND == SP
    CALLOC(work, size+1, double, ON, AMSC);
#endif
    CALLOC(c, maxiter+1, double, ON, AMSC);
    CALLOC(s, maxiter+1, double, ON, AMSC);
    CALLOC(g, maxiter+1, double, ON, AMSC);
    CALLOC(y, maxiter+1, double, ON, AMSC);

    /* Set up v^1 and g^0. */
    INNER(rnorm, r, r, size);
    /* Check that r has been loaded with something. */
    if(rnorm == 0.) {
        fprintf(stderr,
            "FLE-gmres: Asking for solution with null RHS.\n"
            );
        exit(0);
    }
    rnorm = sqrt(rnorm);
    r0 = rnorm;
    for(i=1; i <= size; i++) p[i] = r[i] / rnorm;
    for(i=1; i <= maxiter; i++) g[i] = 0.0;         // xxx Yves xxx  here was the 'overflow' error
    g[1] = rnorm;
    stoptimer;
    conjtime += dtime;
    /* Iterate up to maxiter times to drive relative L2 norm (rnorm) less
         than tol. DO AT LEAST 1.*/
#if ITRDAT == ON
    printf("\n\nGMRES RESIDUAL MONITOR:\n");
#endif
    for(iter = 1; (iter <= maxiter) && ((rnorm > *tol) || (iter == 1)); iter++) {
        starttimer;
        /* allocate the back vectors if they haven't been already */
        if(bv[iter] == NULL) {
            CALLOC(bv[iter], size+1, double, ON, AMSC);
            CALLOC(bh[iter], iter+2, double, ON, AMSC);
        }

        /* Save p as the v{iter}. */
        for(i=1; i <= size; i++) bv[iter][i] = p[i];
        stoptimer;
        conjtime += dtime;

        /* Form Av{iter}. */
        /* c2.0 At this point, p is in field-point order, but as the strength
             vector for apply, it needs to be re-ordered in singularity order. */

#if PRECOND == OL
        starttimer;
        for(i=1;i<=size;i++) {
            ap[i] = p[i];
            p[i] = 0;
        }
        /* mulPrecond applies the precondmats to sys->p to get sys->q (the inverse
             of any other piece of the algorithm).  Since ap=sys->p is loaded above
             and is in field point ordering and we have zeroed p=sys->q, we are all
             set to go. On return, p is in singularity order. */
        mulPrecond(sys, size);
        stoptimer;
        prectime += dtime;
#endif
#if PRECOND == SP
        starttimer;
        /* With the SP preconditioner mulPrecond applies the precondmats to sys->p
             to get sys->q (the inverse of any other piece of the algorithm). But it
             operates on ap=sys->p in singularity ordering, so we have to return it to
             that before the call since this vector is in field point ordering.
             On return, p is in singularity order. */
        nf = fptlist;
        for(i=1;i<=size;i++) {
            work[nf->index] = p[nf->index];
            p[nf->index] = 0;
            nf = nf->next;
        }
        spmulPrecond(sys, work, size);
        stoptimer;
        prectime += dtime;
#endif
#if PRECOND == NONE
        nf = fptlist;
        for(i=1;i<=size;i++) {
            ap[i] = p[nf->index];
            nf = nf->next;
        }
        nq = snglist;
        for(i=1;i<=size;i++) {
            p[nq->index[0]] = ap[i];
            nq = nq->next;
        }
#endif

        /* apply does ap = [sys]p, ap is zeroed on entry. */
        apply(sys, p, ap, size);

        starttimer;

        /* Make v^{iter+1} orthogonal to v^{i}, i <= iter. */
        for(j=1; j <= iter; j++) {
            INNER(hi, ap, bv[j], size);
            /* Use modified Gram-Schmidt. */
            for(i=1; i <= size; i++) ap[i] -= hi * bv[j][i];
            bh[iter][j] = hi;
        }

        /* Normalize v^{iter+1}. */
        INNER(apnorm, ap, ap, size);
        apnorm = sqrt(apnorm);
        for(i=1; i <= size; i++) p[i] = ap[i]/apnorm;
        bh[iter][iter+1] = apnorm;

        /* Apply rotations to new h column. */
        for(i=1; i < iter; i++) {
            hi = bh[iter][i];
            hip1 = bh[iter][i+1];
            bh[iter][i] = c[i] * hi - s[i] * hip1;
            bh[iter][i+1] = c[i] * hip1 + s[i] * hi;
        }

        /* Compute new rotations. */
        hi = bh[iter][iter];
        hip1 = bh[iter][iter+1];
        length = sqrt(hi * hi + hip1 * hip1);
        c[iter] = hi/length;
        s[iter] = -hip1/length;

        /* Apply new rotations. */
        bh[iter][iter] = c[iter] * hi - s[iter] * hip1;
        bh[iter][iter+1] = c[iter] * hip1 + s[iter] * hi;
        hi = g[iter];
        g[iter] = c[iter] * hi;
        g[iter+1] = s[iter] * hi;
        /* c2.0 We set the norm to be relative, absolute is another approach.
             Use of the maxnorm could be argued for too. */
        rnorm = ABS(g[iter+1]) / r0;

        stoptimer;
        conjtime += dtime;

#if ITRDAT == ON
        fprintf(stdout, "   Iteration #%d,  ||res|| = %g\n", iter, rnorm);
#endif
        fflush(stdout);
    }
    /* Decrement from the last increment. */
    iter--;

#if ITRDAT == ON
    printf("    Total GMRES iters = %d\n", iter);
#endif
    starttimer;
    /* Compute solution, note, bh is bh[col][row]. */
    for(i=1; i <= iter; i++) y[i] = g[i];
    for(i = iter; i > 0; i--) {
        y[i] /= bh[i][i];
        for(j = i-1; j > 0; j--) {
            y[j] -= bh[i][j]*y[i];
        }
    }

    for(i=1; i <= size; i++) {
        ap[i] = 0.0;
        for(j=1; j <= iter; j++) {
            ap[i] += y[j] * bv[j][i];
        }
    }
    stoptimer;
    conjtime += dtime;

    /* Undo the preconditioning to get the real q. */
#if PRECOND == OL
    starttimer;
    for(i=1; i <= size; i++) { p[i] = 0; }
    mulPrecond(sys, size);
    stoptimer;
    prectime += dtime;
#endif
#if PRECOND == SP
    starttimer;
    /* With the SP preconditioner mulPrecond applies the precondmats to sys->p
         to get sys->q (the inverse of any other piece of the algorithm). But it
         operates on ap=sys->p in singularity ordering, so we have to return it to
         that before the call since this vector is in field point ordering.
         On return, p is in singularity order. */
    nf = fptlist;
    for(i=1;i<=size;i++) {
        work[nf->index] = ap[nf->index];
        nf = nf->next;
        p[i] = 0;
    }
    spmulPrecond(sys, work, size);
    stoptimer;
    prectime += dtime;
#endif
#if PRECOND == NONE
    /* Undo the field point ordering for the ordering of the input. */
    for (nf=fptlist, nq=snglist; (nf!=NULL) && (nq!=NULL); nf = nf->next, nq = nq->next) {
        p[nq->index[0]] = ap[nf->index];
    }
#endif

    if(rnorm > *tol) {
#if NOWARN == OFF
        fprintf(stdout, "FLW-gmres: exiting without converging\n");
#endif
    }
    *tol = rnorm;
    return(iter);
}

/*
apply multiplies the vector q by the influence coefficient matrix
with or without preconditioning to obtain the vector p. It is
NOT assumed that p has been initialized to zero.
*/

void apply(ssystem *sys, double *q, double *p, int size)
{
    extern double dirtime, uptime, downtime, evaltime, prectime;
    int i;
    ASSERT(p == sys->p);
    ASSERT(q == sys->q);

    for(i=1; i <= size; i++) p[i] = 0;

    starttimer;
    mulDirect(sys);
    stoptimer;
    dirtime += dtime;

    starttimer;
    mulUp(sys);
    stoptimer;
    uptime += dtime;

#if DUPVEC == ON
    dumpLevOneUpVecs(sys);
#endif

#if DNTYPE == NOSHFT
    mulDown(sys);       /* do downward pass without local exp shifts */
#endif

#if DNTYPE == GRENGD
    mulDown(sys);                    /* do heirarchical local shift dwnwd pass */
#endif
    stoptimer;
    downtime += dtime;

    starttimer;
#if MULTI == ON
    mulEval(sys);       /* evaluate either locals or multis or both */
#endif
    stoptimer;
    evaltime += dtime;

#if OPCNT == ON
    printops();
    exit(0);
#endif
}
