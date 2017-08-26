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

/* # ***** sort to /src/header
   # ***** */
/* header where rusage and time structs are defined */

#ifdef FOUR
#define NOTOTHER 1
#include <sys/time.h>
#include <sys/resource.h>
struct rusage timestuff;
#endif

#ifdef FIVE
#define NOTOTHER 1
#include <sys/types.h>
#include <sys/param.h>
#include <sys/times.h>
struct tms timestuff;
#endif

/* define macros for time and resident memory usage checks */

static double dtime = 0.0;
static long my_stime, my_utime;

#ifdef NOTOTHER

#ifdef FOUR			/* 4.2,3BSD (tested: Sun4, IBM6000, DEC5000) */
#define starttimer getrusage(RUSAGE_SELF, &timestuff); \
my_stime = timestuff.ru_utime.tv_sec; \
my_utime = timestuff.ru_utime.tv_usec
#define stoptimer getrusage(RUSAGE_SELF, &timestuff); \
dtime = (double)(timestuff.ru_utime.tv_sec - my_stime) \
        + 1.0e-6*(double)(timestuff.ru_utime.tv_usec - my_utime)
#define DUMPRSS			/*  */
#endif /* FOUR */

#ifdef FIVE			/* for System V (tested: HP300) */
#define starttimer times(&timestuff); \
my_utime = timestuff.tms_utime
#define stoptimer times(&timestuff); \
dtime = (timestuff.tms_utime)-my_utime; \
dtime /= HZ
#define DUMPRSS			/*  */
#endif /* FIVE */

#else				/* default - no timers */

#define OTHER			/*  */
#define starttimer		/*  */
#define stoptimer		/*  */
#define DUMPRSS			/*  */

#endif /* NOTOTHER */





