include ./Makefile.config

# DALTON_LIBS, MODULES and SUBDIRS must be modified at the same time!
DALTON_LIBS = \
-Llinears -llinears     \
-Ldft -ldft		\
-Lgeomopt -lgeomopt     \
-LLSint -lLSint         \
-L$(FMMDIR) -lmm        \
-Lpdpack -lpdpack       \
-Llsutil -llsutil       

MODULES = PD_OBJ LSU_OBJ \
	LIN_OBJ DFT_OBJ GEO_OBJ FMM_OBJ BSM_OBJ LSI_OBJ
SUBDIRS = pdpack lsutil dft LSint geomopt linears $(FMMDIR)

PDLIB = pdpack/pdlib.a

OBJPRG = linears/lsdalton_wrapper.o

#
#     Most common build of Dalton
#
lsdalton.x: $(MODULES)
	@echo "---------------> Linking sequential lsdalton.x ..."
	$(F90) $(FFLAGS) \
	-o $(INSTALLDIR)/lsdalton.x $(OBJPRG) $(DALTON_LIBS) $(LIBS)
#
#       Linux version of the program
#
linux.x: $(MODULES)
	@echo "---------------> Linking sequential lsdalton.x ..."
	$(F90) $(FFLAGS) \
	-o $(INSTALLDIR)/lsdalton.x $(OBJPRG) $(DALTON_LIBS) $(LIBS)
#
#     Linux MPI parallel build (first create sequential build dalton.x, then dalpar.x)
#
#linuxparallel.x: $(MODULES)
#	@echo "---------------> Linking parallel dalpar.x ..."
#	$(F90) $(FFLAGS) \
#	-o $(INSTALLDIR)/lsdalpar.x $(OBJPRG) $(OBJSMEM) \
#	$(OBJSLAVE) $(OBJSAMFI) $(DALTON_LIBS) $(LIBS) \
#	$(MPI_LIB_PATH) $(MPI_LIB) 
##
##	MPI parallel build (first create sequential build dalton.x, then dalpar.x)
##
#parallel.x: dalton.x
#	@echo "---------------> Linking parallel dalpar.x ..."
#	$(F90) $(FFLAGS) \
#	-o $(INSTALLDIR)/lsdalpar.x $(OBJPRG) $(OBJSLAVE) $(OBJSAMFI) \
#	$(DALTON_LIBS) $(LIBS) $(MPI_LIB_PATH) $(MPI_LIB)
##
##	Build with PVMe
##	We will never need MPILIB when using PVM
##       First build master, then slave
##
#dalpvm.x: $(MODULES)
#	$(F90) $(FFLAGS) \
#	-o $(INSTALLDIR)/lsdalpvm.x $(OBJPRG) $(OBJSMEM) \
#	$(OBJSLAVE) $(OBJSAMFI) $(LIBS) \
#	$(DALTON_LIBS) $(PVM_LIB_PATH) \
#	$(PVM_LIBS) $(PVM_INC_PATH)
#	$(F90) $(FFLAGS) -o $(INSTALLDIR)/lsdal_node.x $(OBJNODE) $(OBJSLAVE) \
#	$(OBJSAMFI) $(LIBS) $(DALTON_LIBS) $(PVM_LIB_PATH) \
#	$(PVM_LIBS)  $(PVM_INC_PATH)
##
##	Not tested MPI on Cray yet. Thus no MPILIB, and no OBJSMXM nor
##	OBJSEIS
##
#cray.x: $(MODULES)
#	$(F90) $(FFLAGS) -o $(INSTALLDIR)/lsdalton.x $(OBJPRG) \
#	$(OBJSLAVE) $(OBJSAMFI) $(LIBS) $(DALTON_LIBS)
##
##	This is a proper build for the Cray-T3D
##
#t3d.x: $(MODULES)
#	$(F90) $(FFLAGS) $(MPI_LIB_PATH) $(MPI_LIB) \
#	-o $(INSTALLDIR)/lsdalpar.x $(OBJPRG) $(OBJSLAVE) \
#	$(LIBS) $(OBJSAMFI) $(DALTON_LIBS)
#
#t90.x: $(MODULES)
#	$(F90) $(FFLAGS) $(MPI_LIB_PATH) $(MPI_LIB) \
#	-o $(INSTALLDIR)/lsdalpar.x $(OBJPRG) $(OBJSLAVE) \
#	$(LIBS) $(OBJSAMFI) $(DALTON_LIBS)
#
#int_tester: dft/int_tester.f90 $(MODULES)
#	$(F90) $(F90FLAGS)  -o $@ dft/int_tester.f90 \
#	$(OBJSAMFI) $(OBJSLAVE) $(OBJSMEM) $(DALTON_LIBS) $(LIBS)
#
#
#	Update all dependencies
#
depend  :
	for i in $(SUBDIRS); do ( cd $$i  && touch Makefile.depend && $(MAKE) $@ ); done

#
#	Make it a bit cleaner, remover all .o/.lst/*.f -files
#
clean :
	for i in $(SUBDIRS); do ( cd $$i  && $(MAKE) $@ ); done
	$(RM) -f *~
#
#	We remove the entire source code as well if we do not plan to debug
#
veryclean :
	$(RM) -rf $(SUBDIRS)

PD_OBJ	:
	cd pdpack && $(MAKE) all

LSU_OBJ :
	cd lsutil && $(MAKE) all

DFT_OBJ	: LSU_OBJ
	cd dft && $(MAKE) all

GEO_OBJ : LSU_OBJ LSI_OBJ 
	cd geomopt && $(MAKE) all

LSI_OBJ : LSU_OBJ
	cd LSint && $(MAKE) all

LIN_OBJ : LSU_OBJ LSI_OBJ GEO_OBJ
	cd linears && $(MAKE) all

FMM_OBJ : LSU_OBJ
	cd $(FMMDIR) && $(MAKE) all

BSM_OBJ	:
	if grep -q HAVE_BSM Makefile.config; then cd lsutil/bsm && $(MAKE) AR="$(AR)" install prefix=../..; fi

ftnchek: pre
	ftnchek $(CHEKFLAGS) */*.i > dalton.ftnchek

pre :
	for i in $(SUBDIRS); do ( cd $$i  && $(MAKE) $@ ); done

check : default
	cd test && ./TEST -y linsca
