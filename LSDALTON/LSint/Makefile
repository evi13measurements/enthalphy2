all: libLSint.a

include ../Makefile.config

INCLUDES := $(INCLUDES) -I../lsutil -I../mm

# Keep the file lists sorted for easier modification.

F90_OBJS = \
	Build_AOBATCH.o\
        BuildBasis.o   \
	dalton-input.o\
	Molecule.o\
	Integral_input.o\
	Integral_interface.o\
	Integral_interface_test.o\
	II_dft_int.o \
	II_dft_ksm.o \
	II_XC_interface.o\
	tco_IntegralDriver.o \
        ThermiteDistribute.o \
        ThermiteDistributeExchange.o \
        ThermiteDistributeExchangeWorker.o \
        ThermiteDriver.o \
        ThermiteOverlapDistribution.o \
        ThermiteIntegrals.o \
        ODbatches.o     \
        ReadMolFile.o  \
	st_Integral_interface.o\
	ls_IntegralInterface.o\
	II_Fragment.o\
	tco_IntegralInterface.o \
	SphCartMatrices.o \
	linsolv_df.o 
F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

# we should really have a separate target for bsm to properly
# track dependencies but we would not like it to be too troublesome.
libLSint.a: $(F_OBJS) $(F90_OBJS) 
#	if grep -q HAVE_BSM ../Makefile.config; then cd bsm && $(MAKE) install prefix=../.. AR="$(AR)"; fi
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(F_OBJS:.o=.F) $(F90_OBJS:.o=.f90)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -- $(F_OBJS:.o=.F) $(F90_OBJS:.o=.f90) > Makefile.depend

include Makefile.depend
######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
#-include Makefile.depend
###################################################################### 
