#include <stdlib.h>
#include <stdio.h>

#include "memtracker.h"

#define MAXCNT 5000000

void **mtlist = NULL;
size_t mtcnt;

void mtinit()
{
    mtcnt = 0;
    mtlist = malloc(MAXCNT * sizeof(void *));
}

void mtadd(void *ptr)
{
    if ((mtcnt < MAXCNT) && (mtlist != NULL)) {
        mtlist[mtcnt] = ptr;
        mtcnt++;
    }
}

void mtclear()
{
    size_t j;
    
    if (mtlist != NULL) {
#if TIMDAT
        fprintf(stderr, "memtracker: Clearing %zu allocations...", mtcnt);
#endif
        for (j = 0; j < mtcnt; j++) {
            if (mtlist[j] != NULL) free(mtlist[j]);
        }
        mtcnt = 0;
        free(mtlist);
#if TIMDAT
        fprintf(stderr, "done.\n");
#endif
    }
}
