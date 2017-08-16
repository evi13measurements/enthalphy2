      subroutine mlcc_wr_aoint(h1ao,nnbst,overlap)
C     Purpose: Save AO integrals for MLCC code
C     Save to MLCC_AOINT
C
C     Author: Sarai Dery Folkestad
C     Date:   01.03 2017
#include "implicit.h"
      dimension h1ao(*)
      dimension overlap(*)
C     
      luoain = -1
      luovlp = -1
C
      call gpopen(luoain,'MLCC_AOINT','UNKNOWN','SEQUENTIAL',
     &        'FORMATTED',idummy,.false.)
      rewind(luoain)
C
      write(luoain,*)(h1ao(i),i=1,nnbst)
C
      call gpclose(luoain,'KEEP')
C
      call gpopen(luovlp,'MLCC_OVERLAP','UNKNOWN','SEQUENTIAL',
     &        'FORMATTED',idummy,.false.)
      rewind(luovlp)
      write(luovlp,*)(overlap(i),i=1,nnbst)
      call gpclose(luovlp,'KEEP')
C
      end subroutine mlcc_wr_aoint