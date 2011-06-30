all: liblinears.a

include ../Makefile.config

INCLUDES := $(INCLUDES) -I../lsutil -I../LSint -I../geomopt

# Keep the file lists sorted for easier modification.
# PLEASE KEEP THEM SORTED ALPHABETICALLY!!!

F_OBJS = \
        cholesky.o \
	invroutines.o

#ajt Birefring compiles extremely slowly with pgf, so I disable compiling
#ajt it until I have gotten around to simplifying the code.
#ajt Calls to it from linscf.f90 have been commented out.
#ajt0410 Temporarily added again for Radovan's merging with other versions

F90_OBJS = \
        arh_debug.o \
        arh_density.o \
        arh_driver.o \
        average_util.o \
        configuration.o \
	dalton_interface.o \
        dd_utilities.o \
        dd_utilities_unres.o \
	debug.o		\
        decomp.o        \
	densopt.o	\
	dens_subspace.o \
	diag.o		\
	diis.o		\
	dsm.o		\
	dsm_xterm.o     \
        ecdata_module.o \
	fock-eval.o	\
        leastchange.o \
        levelshift.o \
        lsdalton_wrapper.o \
	lowdin.o \
        LSDALTON.o        \
        LSDALTON_RESPONSE.o        \
        LSDALTON_RESPONSE_type.o        \
	lsdalton_modules.o \
        LSlib.o\
        LS_optimizer.o\
	minimize.o \
        optimloc.o \
	print_MO.o\
        prop_contribs.o \
        queue-operations.o \
        queue-operations4.o \
        queue-typedef.o \
        queue-typedef4.o \
        response_driver.o \
        response_prop.o \
	response_prop_type.o\
        rsp_complex.o \
        rsp_cmplx_sym.o\
        rsp_equations.o \
        rsp_precond.o \
        rsp_utilities.o \
        rsp_solver.o    \
        rsp_std_sym.o  \
        rsp-typedef.o  \
        scfopt-typedef.o \
	scratch.o       \
	soeo-loop.o \
	soeo-transform.o \
	soeo-dens.o \
	soeo-matop.o \
	soeo-stability.o \
	soeo-redspace.o \
	soeo-typedef.o \
	soeo-itemroutines.o \
	soeo-debug.o \
	starting_guess.o \
	statistics.o   	\
        trilevel.o 
#	cgop.o
C_OBJS = memory.o

F_CPPS = $(F_OBJS:.o=.i)

pre: $(F_CPPS)

# we should really have a separate target for bsm to properly
# track dependencies but we would not like it to be too troublesome.
liblinears.a: $(F_OBJS) $(F90_OBJS) $(C_OBJS) 
	$(AR) $(ARFLAGS) $@ $?

clean : 
	$(RM) -f *.o *~ *.f *.lst *.i *.L i.*.f v.i.*.f *.a *.mod *.MOD


Makefile.depend depend: $(F_OBJS:.o=.F) $(F90_OBJS:.o=.f90) $(C_OBJS:.o=.c)
	../tools/f90dep.pl $(CPPFLAGS) $(INCLUDES) -M../lsutil -M../LSint -M../geomopt -- $(F_OBJS:.o=.F) $(F90_OBJS:.o=.f90) > Makefile.depend

include Makefile.depend
######################################################################
# read dependencies using sinclude (silent include)
# which continues if Makefile.depend is not readable (incl. doesn't exist!).
# Makefile.depend can be created with 'make depend'.
#-include Makefile.depend
###################################################################### 

# Gao: Fortran compiler creates a file .mod for each module name.
# This confuses the gnumake, which thinks that .mod files are
# Modula2 source code files. We must override this built-in rule
# by adding the following line:
%.o: %.mod
