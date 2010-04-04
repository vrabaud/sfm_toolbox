#
# Unix/Linux Makefile for MATLAB interface to sba
#

MEX=mex
MEXCFLAGS=-I.. -O -fPIC #-g
# WHEN USING LAPACK, CHANGE THE NEXT TWO LINES TO WHERE YOUR COMPILED LAPACK/BLAS & F2C LIBS ARE!
LAPACKBLASLIBS_PATH=/usr/lib
F2CLIBS_PATH=/usr/lib

# I had to specify the absolute path to the libs, otherwise mex linked against their dynamic versions...
INTFACESRCS=sba.c
INTFACESOBJ=sba.mexglx
SBALIB=../libsba.a
# use what follows for 64 bits
LAPACKLIBS=$(LAPACKBLASLIBS_PATH)/liblapack_pic.a $(F2CLIBS_PATH)/libf2c.a

#LAPACKLIBS=$(LAPACKBLASLIBS_PATH)/liblapack.a $(LAPACKBLASLIBS_PATH)/libblas.a $(F2CLIBS_PATH)/libf2c.a
                                 # On systems with a FORTRAN (not f2c'ed) version of LAPACK, libf2c.a is
                                 # not necessary; on others, libf2c.a comes in two parts: libF77.a and libI77.a

LIBS=$(SBALIB) $(LAPACKLIBS) #-ldl

dummy: $(INTFACESOBJ) projac.so

$(INTFACESOBJ): $(INTFACESRCS)
	$(MEX) $(MEXCFLAGS) $(INTFACESRCS) $(LIBS)

projac.so: projac.c
	gcc -Wall -fPIC -O3 -shared -o projac.so projac.c

clean:
	@rm -f $(INTFACESOBJ) projac.so

realclean cleanall: clean

depend:
	makedepend -f Makefile $(INTFACESRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

