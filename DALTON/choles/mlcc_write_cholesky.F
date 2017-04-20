      subroutine mlcc_rewr_chol(indred,work,lwork)
C
C     Purpose: 
C     Read cholesky vectors and write to new file MLCC_CHOLESKY
C     for use in mlcc module
C
C     MLCC_CHOLESKY:
C     nreduc(1),numcho(1)
C     indred  
C     cholesky vectors in AO basis
C  
C     Author: Sarai Folkestad
C
C     Subroutine parameters:
C     INDRED: array of reduced indices
C
C     Used from common blocks:
C     nreduc(1): dimension of reduced diagonal 
C     numcho(1): J, number of cholesky vectors
C     ((1) - because no symmetry)  
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
      double precision work(lwork),indred(nreduc(1))
C
      chovec = 1
      koff1 = chovec + nnbst(1)     
      lumlch = -2
C
      call gpopen(lumlch,'MLCC_CHOLESKY','UNKNOWN','SEQUENTIAL',
     &        'FORMATTED',idummy,.false.)
      rewind(lumlch)
C
C     
C     Write nreduc(1),numcho(1) to file
      write(lumlch,*)numcho(1) 
C
C     Write reduced index array
C     write(lumlch,*)(indred(i),i=1,nreduc(1))
C  
      do j=1,numcho(1)
C
C     Read and write cholesky
C
        call cc_getcho(work(chovec),1,1,j,indred,work(koff1),nreduc(1))
        write(lumlch,*)(work(chovec + I - 1),I=1,nnbst(1))
C
      enddo
      call gpclose(lumlch,'KEEP')
C  
      end subroutine mlcc_rewr_chol
C
