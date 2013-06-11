#
# makefile for Perple_X 6.6.8
#
# To compile the Perple_X programs with this file first edit the compiler
# variables so that they are consistent with the fortran installation on 
# your system. Then type the command
#
#                             make <target>
#
# where target can be one or more of the programs you want compiled, e.g.,
# "make actcor build fluids ctransf frendly meemum pstable pspts psvdraw pssect pt2curv vertex werami"
# will make the most commonly used programs. 
#
#
# NOTE: file conversion from windows to unix may result in a missing carraige
# return/line feed at the end of Perple_X source files. some compilers  
# report this as an error, in which case edit the source file and add a blank line
# at the end of the file.
#
# JADC, Feb 18, 2012    
# 
##################### COMPILER VARIABLES ####################################  
#                   
#                  COMP77 = the name of the local fortran 77 compatible compiler
#COMP77 = g77
COMP77 = gfortran
#                  FFLAGS = the desired compile options
#FFLAGS = -C -O3 -Wpedantic -Wunused 

# pappel: for use with gfortran
# JADC 1/21/13: O2 and O3 cause fp errors in the speciation routine speci2 in gfortran, the optimization
# seems to work if local variables are initialized to zero (-finit-local-zero, even though there are no 
# uninitialized variables). use -ffpe-trap=zero,overflow,underflow to catch fp errors and -fcheck=all 
# for safe run-time code.

FFLAGS = -m64 -O3 -static-libgfortran

# WFM Added 2007Sep05, PAPPEL 2010SEPT08: for 6.6.0
MYOBJ = actcor build fluids ctransf frendly htog meemum pstable pspts psvdraw pssect pt2curv vertex werami unsplt

all:  $(MYOBJ)

clean: 
	rm -f *.o $(MYOBJ)
###################### TARGETS FOR FORTRAN PROGRAMS #########################   
# 
actcor: actcor.o tlib.o 
	$(COMP77) $(FFLAGS) $@.o tlib.o -o $@

build: build.o tlib.o rlib.o flib.o 
	$(COMP77) $(FFLAGS) build.o tlib.o rlib.o flib.o -o $@ 

fluids: fluids.o tlib.o flib.o 
	$(COMP77) $(FFLAGS) $@.o tlib.o flib.o -o $@

ctransf: ctransf.o tlib.o 
	$(COMP77) $(FFLAGS) $@.o tlib.o -o $@

frendly: frendly.o tlib.o rlib.o flib.o olib.o clib.o dlib.o
	$(COMP77) $(FFLAGS) $@.o rlib.o tlib.o flib.o olib.o clib.o dlib.o -o $@

#hptover: hptover.o
#	$(COMP77) $(FFLAGS) $@.o -o $@

htog: htog.o
	$(COMP77) $(FFLAGS) $@.o -o $@

meemum: meemum.o 
	$(COMP77) $(FFLAGS) meemum.o -o $@   

pstable: pstable.o pslib.o pscom.o  tlib.o cont_lib.o 
	$(COMP77) $(FFLAGS) $@.o pslib.o pscom.o  tlib.o cont_lib.o -o $@

pspts: pspts.o pslib.o tlib.o pscom.o 
	$(COMP77) $(FFLAGS) $@.o pslib.o tlib.o pscom.o -o $@

psvdraw: psvdraw.o pslib.o tlib.o pscom.o
	$(COMP77) $(FFLAGS) $@.o pslib.o tlib.o pscom.o -o $@

pssect: psect.o pscom.o pslib.o tlib.o rlib.o flib.o clib.o  dlib.o 
	$(COMP77) $(FFLAGS) psect.o pscom.o pslib.o tlib.o rlib.o flib.o clib.o dlib.o -o $@

pt2curv: pt2curv.o tlib.o
	$(COMP77) $(FFLAGS) $@.o tlib.o -o $@

unsplt: unsplt.o tlib.o rlib.o flib.o clib.o dlib.o
	$(COMP77) $(FFLAGS) $@.o rlib.o tlib.o flib.o clib.o dlib.o -o $@

vertex: vertex.o 
	$(COMP77) $(FFLAGS) vertex.o -o $@

werami: werami.o 
	$(COMP77) $(FFLAGS) werami.o -o $@

# targets missing from '07:
#rk: rk.o flib.o tlib.o
#	$(COMP77) $(FFLAGS) $@.o tlib.o flib.o -o $@
#ps_p_contor: ps_p_contor.o pslib.o
#	$(COMP77) $(FFLAGS) $@.o pslib.o -o $@
#ge0pt: ge0pt.o
#	$(COMP77) $(FFLAGS) $@.o -o $@
#satsurf: satsurf.o flib.o tlib.o
#	$(COMP77) $(FFLAGS) $@.o flib.o tlib.o -o $@
#rewrite,sox,gox,cohscont 

#################################################################################
# Explicitly make objects to override the default compiler choice implemented
# on some machines.

actcor.o: actcor.f
	$(COMP77) $(FFLAGS) -c actcor.f
build.o: build.f
	$(COMP77) $(FFLAGS) -c build.f
fluids.o: fluids.f
	$(COMP77) $(FFLAGS) -c fluids.f
cont_lib.o: cont_lib.f
	$(COMP77) $(FFLAGS) -c cont_lib.f
ctransf.o: ctransf.f
	$(COMP77) $(FFLAGS) -c ctransf.f
frendly.o: frendly.f
	$(COMP77) $(FFLAGS) -c frendly.f
ge0pt.o: ge0pt.f
	$(COMP77) $(FFLAGS) -c ge0pt.f
#hptover.o: hptover.f
#	$(COMP77) $(FFLAGS) -c hptover.f
htog.o: htog.f
	$(COMP77) $(FFLAGS) -c htog.f
psect.o: psect.f
	$(COMP77) $(FFLAGS) -c psect.f
psvdraw.o: psvdraw.f
	$(COMP77) $(FFLAGS) -c psvdraw.f
pscom.o: pscom.f
	$(COMP77) $(FFLAGS) -c pscom.f
pspts.o: pspts.f
	$(COMP77) $(FFLAGS) -c pspts.f
pstable.o: pstable.f
	$(COMP77) $(FFLAGS) -c pstable.f

# NEXT LINE MODIFIED BY pappel (PA@MIN.UNI-KIEL.DE) 2010SEPT08: CHANGED ptcurv.o TO pt2curv.o	

pt2curv.o: pt2curv.f
	$(COMP77) $(FFLAGS) -c pt2curv.f
rk.o: rk.f
	$(COMP77) $(FFLAGS) -c rk.f
meemum.o: meemum.f
	$(COMP77) $(FFLAGS) -c meemum.f
unsplt.o: unsplt.f
	$(COMP77) $(FFLAGS) -c unsplt.f
vertex.o: vertex.f
	$(COMP77) $(FFLAGS) -c vertex.f
werami.o: werami.f
	$(COMP77) $(FFLAGS) -c werami.f
clib.o: clib.f
	$(COMP77) $(FFLAGS) -c clib.f
dlib.o: dlib.f
	$(COMP77) $(FFLAGS) -c dlib.f
flib.o: flib.f
	$(COMP77) $(FFLAGS) -c flib.f
pslib.o: pslib.f
	$(COMP77) $(FFLAGS) -c pslib.f
nlib.o: nlib.f
	$(COMP77) $(FFLAGS) -c nlib.f
rlib.o: rlib.f
	$(COMP77) $(FFLAGS) -c rlib.f
tlib.o: tlib.f
	$(COMP77) $(FFLAGS) -c tlib.f
olib.o: olib.f
	$(COMP77) $(FFLAGS) -c olib.f

# WFM 2007Sep05
.c.o:
	$(CC) $(CFLAGS) $(INCL) -c $<

.f.o:
	$(FC) $(FFLAGS) $(INCL) -c $<
