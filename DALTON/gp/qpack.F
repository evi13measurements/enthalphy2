!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
C
#ifdef UNDEF
===========================================================================
/* Comdeck qpack_log */
 1-Mar-1989 (c) Copyright Hans Joergen Aa. Jensen
 Package for trace back, and maybe later call trees and timing
 statistics.
 24-Aug-1998 hjaaj: string operations take too much time, have modified QENTER
  (see "Chj1"), have disabled LQTREE counting because not used in this version
  (see "Chj2")

 Modules:
 QDUMP(LUPRI)
 QENTER(IDENT)
 QEXIT(IDENT)
 QINIT
 QTRACE(LUPRI)

 Common blocks:
 QSTACK
 QSTACC (for character variables)

 Ideas for new modules:
 QTREE(LUPRI)
 QSTAT(LUPRI)
 QDUMP(LUPRI) call QTRACE and dump QTREE

 This could maybe be implemented with:
 REAL QTIME
 COMMON /QSTACK/ QTIME(3,MXTREE),
                 LTIME(MXQLVL), LVLQTR(MXTREE),
                 LEVELQ,        LQTREE
 COMMON /QSTACC/ QNAMES(MXQLVL),QTREE(MXTREE)
 in QENTER:
   LEVELQ = LEVELQ + 1
   LQTREE = LQTREE + 1
   IF (LQTREE .LE. MXTREE) THEN
      LVLQTR(LQTREE) = LEVELQ
      QTREE (LQTREE) = IDENT
      QTIME (1,LQTREE) = system time
      QTIME (2,LQTREE) = cpu time
      QTIME (3,LQTREE) = wall time
   END IF
   LTIME(LEVELQ)  = LQTREE
 in QEXIT:
   JTREE  = LVLQTR(LEVELQ)
   IF (JTREE .LE. MXTREE)
  *   QTIME(i,JTREE) = i_time - QTIME(i,JTREE)
   LEVELQ = LEVELQ - 1
 in QINIT:
   LEVELQ = 0
   LQTREE = 0
 with this information, child time can be subtracted from module time:
      TIMCHL = 0.0D0
      DO 100 I = JTREE+1,MIN(LQTREE,MXTREE)
         IF (LVLQTR(I) .LE. LVLQTR(JTREE)) GO TO 200
         TIMCHL = TIMCHL + QTIME(I)
  100 CONTINUE
  200 CONTINUE
      WRITE (LUPRI,'()')
     *   QTREE(JTREE),LVLQTR(JTREE),
     *   QTIME(JTREE),TIMCHL,QTIME(JTREE)-TIMCHL
 a summary can also be made, with accumulated values for each called
 module.
===========================================================================
#endif
C  /* Deck qpackage */
      SUBROUTINE QENTER(IDENT)
C
C  1-Mar-1989 Hans Joergen Aa. Jensen
C  l.r. 980824-hjaaj: removed slow blank fill of QNAMES
C    because of many system routine calls (marked "chj1"),
C    disabled LQTREE (marked "Chj2"), not used in this version
C
C  Package for tracebacks and statistics.
C
      CHARACTER*(*) IDENT
C
#include "priunit.h"
C
      COMMON /QSTACK/ LEVELQ, LQTREE
      PARAMETER (MXQLVL = 100)
      CHARACTER*12  QNAMES(MXQLVL)
      COMMON /QSTACC/ QNAMES
C
Chj1  CHARACTER*8  IDENT8, BLANK8
Chj1  PARAMETER (BLANK8 = '        ')
C
C     Initialize common blocks in BLOCK DATA QBDINIT
C
      EXTERNAL QBDINIT
C
      LEVELQ = LEVELQ + 1
!     print *, IDENT,' entry, LEVELQ =',LEVELQ,
!    &(' ',QNAMES(i),i=1,levelq-1)
Chj2  LQTREE = LQTREE + 1
      IF (LEVELQ .LE. MXQLVL) THEN
Chj1     IDENT8 = IDENT//BLANK8
Chj1     QNAMES(LEVELQ) = IDENT8
Chj1 980824-hjaaj: these string operations take too long time!
         QNAMES(LEVELQ) = IDENT
      END IF
      RETURN
C
      ENTRY QEXIT(IDENT)
Chj1  IDENT8 = IDENT//BLANK8
Chj1  IF (QNAMES(LEVELQ) .EQ. IDENT8) THEN
      ICHECK = MIN(12,LEN(IDENT))
      IF (QNAMES(LEVELQ) .EQ. IDENT(1:ICHECK)) THEN
         LEVELQ = LEVELQ - 1
      ELSE
         WRITE (LUPRI,'(//A/A,I6,2(/2A),//A)')
     *      ' QEXIT error, exit ID does not match entry ID',
     *      ' ----- stack pt =',LEVELQ,
     *      ' ----- entry ID =',QNAMES(LEVELQ),
     *      ' ----- exit  ID =',IDENT(1:ICHECK),
     *      ' Stack dump follows:'
         CALL QDUMP(LUPRI)
         CALL QUIT('QEXIT error, exit ID does not match entry ID')
      END IF
!     print *, IDENT,' exit,  LEVELQ =',LEVELQ,
!    &(' ',QNAMES(i),i=1,levelq)
      RETURN
C
      ENTRY QINIT
C     reinitalize Q package:
      LEVELQ = 0
      LQTREE = 0
      DO I = 1,MXQLVL
         QNAMES(I) = '-noname-'
      END DO
      RETURN
C
      END
C  /* Deck qbdinit */
      BLOCK DATA QBDINIT
      COMMON /QSTACK/ LEVELQ, LQTREE
      PARAMETER (MXQLVL = 100)
      CHARACTER*12  QNAMES(MXQLVL)
      COMMON /QSTACC/ QNAMES
C
      DATA LEVELQ, LQTREE /0, 0/
      DATA QNAMES /MXQLVL*'-noname-'/
      END
C  /* Deck qdump */
      SUBROUTINE QDUMP(LUPRI)
C
      COMMON /QSTACK/ LEVELQ, LQTREE
      PARAMETER (MXQLVL = 100)
      CHARACTER*12  QNAMES(MXQLVL)
      COMMON /QSTACC/ QNAMES
C
      IF (LEVELQ .GT. 0) THEN
         WRITE (LUPRI,'(//A)')
     &      ' QDUMP calls QTRACE for dump of internal trace stack'
         CALL QTRACE(LUPRI)
      END IF
      RETURN
      END
C  /* Deck qtrace */
      SUBROUTINE QTRACE(LUPRI)
C
C  1-Mar-1989 Hans Joergen Aa. Jensen
C
C  Print traceback from Q stack
C
      COMMON /QSTACK/ LEVELQ, LQTREE
      PARAMETER (MXQLVL = 100)
      CHARACTER*12  QNAMES(MXQLVL)
      COMMON /QSTACC/ QNAMES
C
      IF (LEVELQ .EQ. 0) GO TO 9000
      WRITE (LUPRI,'(//A)') ' QTRACE dump of internal trace stack'
      WRITE (LUPRI,'(3(/A))')
     *   ' ========================',
     *   '      level    module',
     *   ' ========================'
      DO 100 I = LEVELQ,1,-1
         IF (I .LE. MXQLVL) THEN
            WRITE (LUPRI,'(I11,4X,A)') I,QNAMES(I)
         ELSE
            WRITE (LUPRI,'(I11,4X,A)') I,'not recorded'
         END IF
  100 CONTINUE
      WRITE (LUPRI,'(A/)') ' ========================'
C920522-hjaaj: do not call TRACE because it will abort on some computers
C9000 CALL our_own_traceback
 9000 CONTINUE
      RETURN
      END
