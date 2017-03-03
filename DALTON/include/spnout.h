! -- FILE: spnout.h --
      LOGICAL DOSD, DODSO, DOFC, DOSDFC, DOPSO, DOSELE, ANISON,
     &        FCFIN, SPNISO, NCSPNI,
CPFP Added SOSOSC
     &        SOS, SOSSPN, SOSOCC, SOSOCS, SOSOSC
Cend-PFP
      COMMON /SPNOUT/ ABUND,                                             ! real*8
     &        ISPPRI, ISOTPS(MXCENT), ISINGUL,                           ! integer
     &        NSTATS, NSTATT, NSTATI, NSTATF, NITRST, NUCSPI,
     &        DOSD, DODSO, DOFC, DOSDFC, DOPSO, DOSELE, ANISON,          ! logical
     &        FCFIN, SPNISO, NCSPNI(MXCENT),
     &        SOS, SOSSPN, SOSOCC, SOSOCS, SOSOSC, SINGUL
! -- end of spnout.h --
