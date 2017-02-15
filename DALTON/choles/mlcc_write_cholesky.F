      subroutine mlcc_rewr_chol(indred,work,lwork)
C
C     Purpose: 
C     Read cholesky vectors and write to new file MLCC_CHOLESKY
C     for use in mlcc module
C
C     Author: Sarai Folkestad
C
C     Subroutine parameters:
C     INDRED: array of reduced indices
C
C     Used from common blocks:
C     nreduc(1): dimension of reduced diagonal 
C     numcho(1): J, number of cholesky vectors
C
#include "implicit.h"
#include "maxorb.h"
#include "maxash.h"
#include "aovec.h"
#include "priunit.h"
#include "iratdef.h"
#include "ccdeco.h"
#include "ccsdsym.h"
#include "choio.h"
C
      dimension work(lwork)
C
      chovec = 1
      koff1 = chovec + nreduc(1)     
      lumlch = -2
C
      call gpopen(lumlch,'MLCC_CHOLESKY','UNKNOWN','SEQUENTIAL',
     &        'UNFORMATTED',idummy,.false.)
C
C             
      do j=1,numcho(1)
C
C
C
        call cc_getcho(work(chovec),1,1,j,indred,work(koff1),nreduc(1))
        write(lumlch)(work(chovec + I - 1),I=1,nreduc(1))
C
      enddo
      call gpclose(lumlch,'KEEP')
C  
      end subroutine mlcc_rewr_chol