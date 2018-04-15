      subroutine mlcc_wr_ao(h1ao,nnbst,overlap)
#include "implicit.h"
      dimension h1ao(*)
      dimension overlap(*)
      luoain = -1 
      luovlp = -1
      call gpopen(luoain,'MLCC_AOINT','UNKNOWN','SEQUENTIAL',
     &        'FORMATTED',idummy,.false.)
      rewind(luoain)
      write(luoain,*)(h1ao(i),i=1,nnbst)
      call gpclose(luoain,'KEEP')
      call gpopen(luovrlp,'MLCC_OVERLAP','UNKNOWN','SEQUENTIAL',
     &        'FORMATTED',idummy,.false.)
      rewind(luovrlp)
      write(luovrlp,*)(overlap(i),i=1,nnbst)
      call gpclose(luovrlp,'KEEP')
C
      end subroutine mlcc_wr_ao