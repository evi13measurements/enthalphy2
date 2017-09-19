      SUBROUTINE PRORB
C
      CALL QUIT('Dummy PRORB called')
      STOP 'STOP - Dummy PRORBO called'
      END
      SUBROUTINE PUNMO
C
      CALL QUIT('Dummy PUNMO called')
      STOP 'STOP - Dummy PUNMO called'
      END
      SUBROUTINE MM_PROPS
C
      CALL QUIT('Dummy MM_PROPS called')
      STOP 'STOP - Dummy MM_PROPS called'
      END
      SUBROUTINE DOMEP
C
      CALL QUIT('Dummy DOMEP called')
      STOP 'STOP - Dummy DOMEP called'
      END
      SUBROUTINE WRIPRO
C
      CALL QUIT('Dummy WRIPRO called')
      STOP 'STOP - Dummy WRIPRIO called'
      END
      SUBROUTINE FORMSUP
C
      CALL QUIT('Dummy FORMSUP called')
      STOP 'STOP - Dummy FORMSUP called'
      END
C  /* Deck GPIO_INI_DALTON */
      SUBROUTINE GPIO_INI_DALTON
C
C 1-Mar-2000 K.Ruud
C 10-Apr-2000 HJAaJ: changed initialization to -9xxx from 0,
C                    because unit 0 is stderr on many systems.
C                    The number -9xxx is used for easier debugging.
C
#include <mxcent.h>
#include <eribuf.h>
#include <priunit.h>
#include <inftap.h>
#include <ccinftap.h>
#include <r12int.h>
C
C     priunit.h:
C
      LUERR  = -8000
      LUW4   = -8004
      LUCMD  = -8005
      LUPRI  = -8006
      LUPOT  = -8009
      NINFO  = 0
      NWARN  = 0
      IPRERR = 0
C
C     inftap.h:
C
      LUCME  = -9003
      LUMOL  = -9004
      LUPROP = -9005
      LUSOL  = -9006
      LUINTA = -9007
      LUONEL = -9008
      LUSUPM = -9009
      LUTLM  = -9010
      LUDA1  = -9011
      LUITMP = -9012
      LU2DER = -9013
      LUDASP = -9014
      LURDR  = -9016
      LURDI  = -9017
      LUGDR  = -9018
      LUGDI  = -9019
      LUGDT  = -9020
      LURDT  = -9021
      LUDFCK = -9022
      LUSFDA = -9023
      LUFDC  = -9024
      LUWLK  = -9025
      LUPAO  = -9026
      LUPAS  = -9027
      LUNR1  = -9028
      LUNR3  = -9029
      LUNR5  = -9030
      LUINTR = -9032
      LUMOLDEN=-9033
      LUPMOM = -9034
      LUMOM  = -9035
      LUEIND = -9036
      LUENUC = -9037
      LUESITE= -9038
      LUEOBAR= -9039
      LUVDWSE= -9040
      LUENSA = -9041
      LUQM3E = -9042
      LUQM3P = -9043
      LUOSCR = -9044
      LUMMOV = -9045
      LUOVER = -9046
      LUNDPF = -9047
      LUSRINT =-9048
      LUNMPF = -9049
C
C     /R12INT/ (WK/UniKA/26-11-2002).
C
      DO I = 1, 5
         LUR12(I) = -9100 - MXCOOR - I
      ENDDO
C
C     /ERITAP/
C
      DO I = 0, MXCOOR
         LUAORC(I) = -10000 - I
      ENDDO
C
C     /RSPTAP/
C
      LUAHSO = -9201
      LUCRV1 = -9202
      LUCRV2 = -9203
      LUXYVE = -9204
      LUCRVE = -9205
      LURSP3 = -9206
      LURSP4 = -9207
      LURSP5 = -9208
      LUMHSO = -9209
      LURSP  = -9210
C
C     /SIRTAP/
C
      LUINTM = -9301
      LUIT1  = -9302
      LUIT2  = -9303
      LUIT3  = -9304
      LUIT5  = -9305
      LUINF  = -9306
      LUH2AC = -9307
      LUSIFC = -9308
C     next two for NEWTRA
      LUORDA = -9309
      LUMINT = -9310
C
      FNSOL  = 'AOSOLINT'
      ABARDR = 'ABACUS.RD'
      ABARDI = 'ABACUS.RDI'
      ABAGDR = 'ABACUS.GD'
      ABAGDI = 'ABACUS.GDI'
      ABAGDT = 'ABACUS.GDT'
      ABARDT = 'ABACUS.RDT'
      ABADFK = 'ABACUS.DFK'
      ABASF  = 'ABACUS.SF'
      ABATLM = 'ABACUS.TLM'
      ABAWLK = 'DALTON.WLK'
      ABAIRC = 'DALTON.IRC'
      ABATRJ = 'DALTON.TRJ'
      ABANR1 = 'ABAENR.RST'
      ABANR3 = 'ABAENR.BVC'
      ABANR5 = 'ABAENR.SVC'
      FNINTM = 'MOTWOINT'
      FNSUPM = 'AOSUPINT'
      FNONEL = 'AOONEINT'
      FNSIFC = 'SIRIFC'
      LBSIFC = 'SIR IPH '
C
      RETURN
      END
      subroutine trace()
         print *, ' westa - dummy dalton routine trace called'
      end
      subroutine calc_mep()
         STOP ' westa - dummy dalton routine calc_mep called'
      end
      subroutine sir_intopen()
         print *, ' westa - dummy dalton routine sir_intopen called'
      end
      subroutine sir_prtinp()
         STOP ' westa - dummy dalton routine sir_prtinp called'
      end
      subroutine tractl()
         STOP ' westa - dummy dalton routine tractl called'
      end
      subroutine dpotrf()
         STOP ' westa - dummy dalton routine dpotrf called'
      end
C -- end of dalton_dummy.F --
