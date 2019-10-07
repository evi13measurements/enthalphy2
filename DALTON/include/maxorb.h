!
!     File: maxorb.h
!
!     MXSHEL = maximum number of shells (insert shell definition here).
!              (if modified: also  change MXSHEL for __CVERSION__ in infpar.h)
!     MXPRIM = maximum number of primitives.
!     MXCORB = maximum number of orbitals (possibly contracted).
!     MAXOCC = maximum number of occupied orbitals
!
!     IF you change any of these parameters you should rebuild with "make".
!
      INTEGER    MXSHEL, MXPRIM, MXCORB, MXORBT, MAXOCC
      PARAMETER (MXSHEL = 1500, MXPRIM = 15000 )
!PFP  usefull for Local-auto
      PARAMETER (MXCORB = 5000, MXORBT = MXCORB*(MXCORB + 1)/2 )
       PARAMETER (MAXOCC = 1500 )
!       PARAMETER (MXCORB = 10000, MXORBT = MXCORB*(MXCORB + 1)/2 )
!       PARAMETER (MAXOCC = 2500 )
!end-PFP

!     MXCORB_CC = max number of orbitals in CC modules
!     (normally less than MXCORB because of a lot of static allocations
!      in the CC module for address pointers)

      INTEGER    MXCORB_CC
!PFP usefull for Local-auto
       PARAMETER (MXCORB_CC = 600 )
!       PARAMETER (MXCORB_CC = 1200 )
!end-PFP
! -- end of maxorb.h --
