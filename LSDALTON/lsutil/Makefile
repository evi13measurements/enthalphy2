all: liblsutil.a

include ../Makefile.config

# Keep the file lists sorted for easier modification.

C_OBJS = matop_csr_aux.o

F90_OBJS = \
	dftType.o\
	BasisinfoType.o\
        common_utilities.o \
        file-operations.o \
	lsutilities.o\
	ls_precision.o\
	ls_IO.o\
	lstensor_operations.o\
	par_mod.o\
	lsmpi.o\
	matop_dense.o\
	matop_csr.o\
	lsmatop_dense.o\
	matop_dense_unrest.o\
	mat-operations.o\
	matop_sparse1.o\
	matrix_defop.o\
	matrix_genop.o\
        matrix_utilities.o \
	matrix_module.o\
	MoleculeType.o\
        IntegralOutputType.o\
	IntegralType.o\
	AOType.o\
	ODType.o\
	LSmatrixType.o\
        memory.o \
	Time.o\
	ks-settings.o\
	TYPE-DEF.o\
        fundamental.o

#       AoBatchType.o\

# we should really have a separate target for bsm to properly
# track dependencies but we would not like it to be too troublesome.
liblsutil.a: $(F90_OBJS) $(C_OBJS)
	if grep -q HAVE_BSM ../Makefile.config; then cd bsm && $(MAKE) install prefix=../.. AR="$(AR)"; fi
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(F90_OBJS:.o=.f90) $(C_OBJS:.o=.c)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -- $(F90_OBJS:.o=.f90) $(C_OBJS:.o=.c) > Makefile.depend

include Makefile.depend
######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
#-include Makefile.depend
###################################################################### 

%.o: %.mod
