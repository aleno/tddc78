CFLAGS = -Nmpi

LFLAGS= -lposix4 -lm

all: blurc thresc

clean:
	-$(RM) *.o blur thresc

blurc: ppmio.o gaussw.o blurfilter.o blurmain.o
	icc -o $@ ppmio.o gaussw.o blurfilter.o blurmain.o $(LFLAGS)

thresc: thresmain.o ppmio.o thresfilter.o
	icc -fast -o $@ thresmain.o ppmio.o thresfilter.o $(LFLAGS)

arc:
	tar cf - *.c *.h *.f90 Makefile Makefile.sgi|gzip - > filters.tar.gz
