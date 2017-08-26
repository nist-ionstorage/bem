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
#include "mulStruct.h"

cube *cstack[1024];		/* Stack used in getnbrs routine. */
long memcount;	       	        /* Allocated memory counter */
long memQ2M;			/* Allocated memory counters by function */
long memQ2L;
long memQ2P;
long memL2L;
long memM2M;
long memM2L;
long memM2P;
long memL2P;
long memMSC;
/* 
  global timer and operation count accumulators
*/
double prectime;	       	/* time spent doing back solve for prec */
double conjtime;	        /* time spent doing everything but A*q */
double dirtime;			/* time for direct part of P*q */
double multime;			/* time for multipole part of P*q */
double uptime;			/* time in mulUp(), upward pass */
double downtime;		/* time in mulDown(), downward pass */
double evaltime;		/* time in mulEval(), evaluation pass */
double structime;               /* time in fastlap to load structures and
                                   get panel attributes.*/
double dirsetime;               /* time in matSetup for direct setup */
double mulsetime;               /* time in matSetup for multi setup */

