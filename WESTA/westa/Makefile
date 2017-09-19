include ../Makefile.config

F_OBJS = \
expcd2.o expcid.o wesdum.o wesgp.o wesmom.o \
wesnad.o wesshk.o westa.o westra.o wfovlp.o \
arhpack.o dalton_dummy.o

F_CPPS = $(F_OBJS:.o=.i)

all: $(F_OBJS)

pre: $(F_CPPS)

westa.x:
	cd .. ; $(MAKE) westa.x

ftnchek: pre
	ftnchek $(CHEKFLAGS) *.i > westa.ftnchek

clean : 
	$(RM) -f *.o *~ *.lst *.i *.L i.*.f v.i.*.f *.a

depend:
	../dalton/tools/check_depend $(INCLUDES) > Makefile.depend

include Makefile.depend
