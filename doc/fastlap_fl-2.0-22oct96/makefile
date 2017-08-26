# Makefile for fastlap Version 2.0
#driverf --   FORTRAN driver to test and play with fastlap
#	     with error checking for sphere problems."
#             Does Green form., single layer, source form.,"
#             and field computations."
#driverc --   C driver to test and play with fastlap."
#             with error checking for sphere problems."
#desing  --   FORTRAN driver to run desingularized sphere."
#             problem in single-layer formulation."
#             Find executable and data in ./desingdir."
#fsp     --   FORTRAN driver to run desingularized free-surface"
#             problem using data from Steve Scorpio (UM)."
#             Find executable and data in ./desingdir."
DIR = .
BIN = .
DBIN = ./desingdir
HDIR = ./src
SDIR = ./src

HDRS = $(HDIR)/mulGlobal.h $(HDIR)/mulStruct.h $(HDIR)/resusage.h

OBJS = $(SDIR)/mulGlobal.o $(SDIR)/mulSetup.o $(SDIR)/mulMats.o $(SDIR)/mulDo.o \
$(SDIR)/mulDisplay.o $(SDIR)/mulMulti.o $(SDIR)/mulLocal.o $(SDIR)/direct.o \
$(SDIR)/calcp.o $(SDIR)/uglyalloc.o 

SRCS = mulGlobal.c mulSetup.c mulMats.c mulDo.c mulDisplay.c mulMulti.c \
mulLocal.c direct.c calcp.c uglyalloc.c  

LIBS = -lm -lc -ldxml
#LIBS = -lm 
#LIBS = /usr/lib/debug/malloc.o -lm

# These compile flags are primarily needed to account for differences in
# the resource usage system call used to get CPU time and memory use data.
# 
# Almost all of the difficulties that users have had in implementing 
# FastLap have come from confusion about which C compiler flags to 
# use.  To date, FastLap has been put up on DEC ALPHA, DEC MIPS, IBM RISC, 
# Sun Sparc, SGI Indigo, SGI 4D, Convex C2, Cray Y.
#
# for most machines running 4.2/3 try:
#CFLAGS = -O -DFOUR
# depending on your optimizer you may need:
#CFLAGS = -O -Olimit 1600 -DFOUR
# debugging will require the g flag, of course:
#CFLAGS =  -g
# for most machines running System V try:
#CFLAGS = -O -DFIVE
# for IBM AIX systems try:
#CC = xlc -O
# for SGI IRIX systems try (at least):
#CC = -cckr
# for HP Snakes you must switch from uglyalloc to malloc and calloc
# because the Snakes are very picky about alignment.  This is done
# by changing the definition of MALCORE and CALCORE to malloc and
# calloc, respectively, in mulGlobal.h.  Having done that, you should
# be able to compile and run without difficulty, although like a lot 
# of codes in our limited PA-RISC experience the highest level of 
# optimization is too agressive to produce an executable which will 
# yield correct results.
#
FFLAGS = -O -w 
CFLAGS = -O -DFOUR -Olimit 1600
#CFLAGS = -g
#FFLAGS = -g -w

default:
	@echo Please specify what to make:
	@echo  "driverf --   FORTRAN driver to test and play with fastlap"
	@echo  "             with error checking for sphere problems."
	@echo  "             Does Green form., single layer, source form.,"
	@echo  "             and field computations."
	@echo  "driverc --   C driver to test and play with fastlap."
	@echo  "             with error checking for sphere problems."
	@echo  "desing  --   FORTRAN driver to run desingularized sphere."
	@echo  "             problem in single-layer formulation."
	@echo  "             Find executable and data in ./desingdir."
	@echo  "fsp     --   FORTRAN driver to run desingularized free-surface"
	@echo  "             problem using data from Steve Scorpio (UM)."
	@echo  "             Find executable and data in ./desingdir."
#	@echo  "extNeum --   C driver for exterior Neumann (sphere) problem"
#	@echo  "             solved by source formulation."
	@echo  "clean   --   Remove object files after a make."

###Primaries:

driverf: $(SDIR)/driverf.o $(SDIR)/geom.o $(SDIR)/fastlapf.o $(OBJS)
	f77 $(FFLAGS) -o $(BIN)/driverf $(SDIR)/driverf.o $(SDIR)/geom.o $(SDIR)/fastlapf.o $(OBJS) $(LIBS)

driverc: $(SDIR)/driverc.o $(SDIR)/fastlap.o $(OBJS)
	cc $(CFLAGS) -o $(BIN)/driverc $(SDIR)/driverc.o $(SDIR)/fastlap.o $(OBJS) $(LIBS)

desing: $(SDIR)/desing.o $(SDIR)/fastlapf.o $(OBJS)
	@echo "  "
	@echo "    You MUST have set PRECOND SP or NONE in"
	@echo "    mulGlobal.h to make a usable desing executable." 
	@echo "  "
	f77 $(FFLAGS) -o $(DBIN)/desing $(SDIR)/desing.o $(SDIR)/fastlapf.o $(OBJS) $(LIBS)

fsp: $(SDIR)/fsp.o $(SDIR)/fastlapf.o $(OBJS)
	@echo "  "
	@echo "    You MUST have set PRECOND SP or NONE in"
	@echo "    mulGlobal.h to make a usable fsp executable." 
	@echo "  "
	f77 $(FFLAGS) -o $(DBIN)/fsp $(SDIR)/fsp.o $(SDIR)/fastlapf.o $(OBJS) $(LIBS)

#extNeum: $(SDIR)/extNeum.o $(SDIR)/fastlap.o $(OBJS)
#	cc $(CFLAGS) -o $(BIN)/extNeum $(SDIR)/extNeum.o $(SDIR)/fastlap.o $(OBJS) $(LIBS)

clean:
	rm $(SDIR)/*.o

###Secondaries:


$(SDIR)/driverf.o: $(SDIR)/driver.f
	cd $(SDIR); f77 $(FFLAGS) -c driver.f -o driverf.o

$(SDIR)/driverc.o: $(SDIR)/driver.c
	cd $(SDIR); cc $(CFLAGS) -c driver.c -o driverc.o

$(SDIR)/desing.o: $(SDIR)/desing.f
	cd $(SDIR); f77 $(FFLAGS) -c desing.f

$(SDIR)/fsp.o: $(SDIR)/fsp.f
	cd $(SDIR); f77 $(FFLAGS) -c fsp.f

#$(SDIR)/extNeum.o: $(SDIR)/extNeum.c
#	cd $(SDIR); cc $(CFLAGS) -c extNeum.c -o extNeum.o

$(SDIR)/geom.o: $(SDIR)/geom.f
	cd $(SDIR); f77 $(FFLAGS) -c geom.f

$(SDIR)/fastlap.o: $(SDIR)/fastlap.c $(HDRS)
	cd $(SDIR); cc $(CFLAGS) -c fastlap.c

$(SDIR)/fastlapf.o: $(SDIR)/fastlapf.c $(HDRS)
	cd $(SDIR); cc $(CFLAGS) -c fastlapf.c

$(SDIR)/mulGlobal.o: $(SDIR)/mulGlobal.c $(HDRS)
	cd $(SDIR); cc $(CFLAGS) -c mulGlobal.c
$(SDIR)/mulSetup.o: $(SDIR)/mulSetup.c $(HDRS)
	cd $(SDIR); cc $(CFLAGS) -c mulSetup.c
$(SDIR)/mulMats.o: $(SDIR)/mulMats.c $(HDRS)
	cd $(SDIR); cc $(CFLAGS) -c mulMats.c
$(SDIR)/mulDo.o: $(SDIR)/mulDo.c $(HDRS)
	cd $(SDIR); cc $(CFLAGS) -c mulDo.c
$(SDIR)/mulDisplay.o: $(SDIR)/mulDisplay.c $(HDRS)
	cd $(SDIR); cc $(CFLAGS) -c mulDisplay.c
$(SDIR)/mulMulti.o: $(SDIR)/mulMulti.c $(HDRS)
	cd $(SDIR); cc $(CFLAGS) -c mulMulti.c
$(SDIR)/mulLocal.o: $(SDIR)/mulLocal.c $(HDRS)
	cd $(SDIR); cc $(CFLAGS) -c mulLocal.c
$(SDIR)/direct.o: $(SDIR)/direct.c $(HDRS)
	cd $(SDIR); cc $(CFLAGS) -c direct.c
$(SDIR)/calcp.o: $(SDIR)/calcp.c $(HDRS)
	cd $(SDIR); cc $(CFLAGS) -c calcp.c
$(SDIR)/uglyalloc.o: $(SDIR)/uglyalloc.c $(HDRS)
	cd $(SDIR); cc $(CFLAGS) -c uglyalloc.c

###Tertiaries:

$(SDIR)/fastlapf.c: $(SDIR)/fastlap.c
# Many Fortran compilers will need an underbar (_) as the last character
# in the name of a routine written in C.  The following sed command will
# insert this. 
	cd $(SDIR); sed s/fastlap\(p/fastlap_\(p/ fastlap.c > fastlapf.c 
# If this is not required by your compiler (IBM under AIX, for instance), 
# use the following line instead.
#	cd $(SDIR); cp fastlap.c fastlapf.c 








