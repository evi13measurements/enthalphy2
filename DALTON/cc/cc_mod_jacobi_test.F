      subroutine cc_mod_jacobi_test(work,lwork)
C
      implicit none 
C
#include "priunit.h"
#include "dummy.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "ccinftap.h"
C
      integer lwork, iside 
      double precision work(lwork)
      character*10 model ! CCSD model variable
C
      character*3 aproxr12
C
      integer kr1am, kr2am, kt1am, kt2am, kend1, lwrk1, iopt
C
      kt1am = 1
      kt2am = kt1am + nt1amx
      kr1am = kt2am + nt2amx
      kr2am = kr1am + nt1amx
      kend1 = kr2am + nt2amx
      lwrk1 = lwork - kend1 + 1
C
C     Read amplitudes from file 
C
      model = 'CCSD      '
C
      iopt = 3
      call cc_rdrsp('R0',idummy,1,iopt,model,work(kt1am),work(kt2am))
C
      write(LUPRI,*) 'IOPT',iopt
C
C      call dzero(work(kt1am), nt1amx)
C      call dzero(work(kt2am),nt2amx)
C
      WRITE(LUPRI,*) 'T1AM'
C
      do i = 1,nrhft
         do a = 1,nvirt
            WRITE(LUPRI,*) nvirt*(i-1)+a,work(kt1am+nvirt*(i-1)+a-1)
         enddo
      enddo
C
      WRITE(LUPRI,*) 'T2AM'
C
      DO I = 1, 15
         WRITE(LUPRI,*) work(kt2am+I-1)
      ENDDO
C
C       call dcopy(nt1amx,t1am,1,work(kr1am),1)
C       call dcopy(nt2amx,t2am,1,work(kr2am),1)
C C
C C      rspim = .true.
C C      call ccrhsn(omeg1,omeg2,t1am,t2am,work(kend1),lwrk1)
C C
C       call dcopy(nt1amx,work(kr1am),1,t1am,1)
C       call dcopy(nt2amx,work(kr2am),1,t2am,1)
C
C      rspim = .false.
C      call ccrhsn(omeg1,omeg2,t1am,t2am,work(kend1),lwrk1)
C
      iside = 1
      call cc_atrr(0.0D0,1,iside,work(kt1am),lwrk1,.false.,
     &             dummy,aproxr12)
C
C      call dcopy(nt1amx,work(kr1am),1,t1am,1)
C      call dcopy(nt2amx,work(kr2am),1,t2am,1)
C
      WRITE(LUPRI,*) 'TRF(SINGLES)'
C
      do i = 1, nrhft
         do a = 1, nvirt
            WRITE(LUPRI,*) nvirt*(i-1)+a,work(kt1am+nvirt*(i-1)+a-1)
         enddo
      enddo
C
      WRITE(LUPRI,*) 'TRF(DOUBLES)'
C
      DO I = 1, 15
         WRITE(LUPRI,*) I,work(kt2am+I-1)
      ENDDO
C
      return
      end