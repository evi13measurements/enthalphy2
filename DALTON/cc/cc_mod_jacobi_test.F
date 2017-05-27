      subroutine cc_mod_jacobi_test(omeg1,omeg2,t1am,t2am,work,lwork)
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
      double precision t1am(*)
      double precision t2am(*)
C
      double precision omeg1(*)
      double precision omeg2(*)
C
      integer lwork, iside 
      double precision work(lwork)
C
      character*3 aproxr12
C
      integer kt1am, kt2am, kend1, lwrk1
C
      kt1am = 1
      kt2am = kt1am + nt1amx
      kend1 = kt2am + nt2amx
      lwrk1 = lwork - kend1 + 1
C
      rspim = .true.
      call ccrhsn(omeg1,omeg2,t1am,t2am,work,lwork)
      rspim = .false.
C
      WRITE(LUPRI,*) 'T1AM'
C
      do i = 1,5
         do a = 1,5
            WRITE(LUPRI,*) nvirt*(i-1)+a,t1am(nvirt*(i-1)+a)
         enddo
      enddo
C
      WRITE(LUPRI,*) 'T2AM'
C
      DO I = 1, 25
         WRITE(LUPRI,*) t2am(I)
      ENDDO
C
      call dcopy(nt1amx,t1am,1,work(kt1am),1)
      call dcopy(nt2amx,t2am,1,work(kt2am),1)
C
      iside = 1
      call cc_atrr(0.0D0,1,iside,work(kt1am),lwrk1,.false.,
     &             dummy,aproxr12)
C
      call dcopy(nt1amx,work(kt1am),1,t1am,1)
      call dcopy(nt2amx,work(kt2am),1,t2am,1)
C
      WRITE(LUPRI,*) 'TRF(SINGLES)'
C
      do i = 1,5
         do a = 1,5
            WRITE(LUPRI,*) nvirt*(i-1)+a,t1am(nvirt*(i-1)+a)
         enddo
      enddo
C
      WRITE(LUPRI,*) 'TRF(DOUBLES)'
C
      DO I = 1, 25
         WRITE(LUPRI,*) t2am(I)
      ENDDO
C
      return
      end