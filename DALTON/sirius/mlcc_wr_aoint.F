      subroutine mlcc_wr_aoint(h1ao,nnbst)
C     Purpose: Save AO integrals for MLCC code
C     Save to MLCC_AOINT
C
C     Author: Sarai Dery Folkestad
C     Date:   01.03 2017
#include "implicit.h"
      dimension h1ao(*)
C     
      luoain = -1
C
      call gpopen(luoain,'MLCC_AOINT','UNKNOWN','SEQUENTIAL',
     &        'FORMATTED',idummy,.false.)
      rewind(luoain)
C
      write(luoain,*)(h1ao(i),i=1,nnbst)
C
      call gpclose(luoain,'KEEP')
C
      end subroutine mlcc_wr_aoint