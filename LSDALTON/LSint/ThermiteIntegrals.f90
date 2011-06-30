!> @file 
!> Contains routines to evaluate integrals of primitive hermite gaussian 
!> Thermite integral module
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
MODULE Thermite_integrals
use typedef
use thermite_OD

!****** INTEGRAND
TYPE Integrand
Character(len=80)    :: Operator
TYPE(OverlapPointer) :: P,Q
Integer              :: nAngmom
Integer              :: nPrimitives
Integer              :: maxContracted
Integer              :: minAngmom
Integer              :: maxAngmom
Integer              :: startAngmom
Integer              :: endAngmom
!Dimensions according to nPrimitives
!Integer,pointer      :: iprimP(:)
!Integer,pointer      :: iprimQ(:)
Real(realk),pointer  :: distance(:,:)
Real(realk),pointer  :: squaredDistance(:)
Real(realk),pointer  :: exponents(:)
Real(realk),pointer  :: reducedExponents(:)
Real(realk),pointer  :: integralPrefactor(:)
!One dimension turning into six: nOrbA,nOrbB,nOrbC,nOrbD,nPassP,nPassQ)
Real(realk),pointer  :: theta(:) 
!One dimension turning into seven: nOrbA,nOrbB,nOrbC,nOrbD,nPassP,nPassQ,iDeriv)
Real(realk),pointer  :: integral(:) 
Real(realk)          :: ORIGO(3)
LOGICAL              :: samePQ
LOGICAL              :: kinetic
END TYPE Integrand

CONTAINS
!> \brief wrapper to recurrence relation for hermite integrals Eq. (9.9.18-9.9.20) in the book \f[  R^{n}_{t+1,u,v}(p,\textbf{R}_{PC}) = t R^{n+1}_{t-1,u,v}(p,\textbf{R}_{PC}) + X_{pC} R^{n+1}_{t,u,v}(p,\textbf{R}_{PC})  \f]
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param WTUV the output of the recurrence 
!> \param WJ000 Eq. (9.9.14) in the book \f$ R^{n}_{000}(p,\textbf{R}_{PC}) = (-2p)^{n}F_{n}(pR_{PC}^{2}) \f$
!> \param sharedTUV a TUVitem which contains tabulated boys function and TUVindexing
!> \param Rpq distance between overlap distribution P and Q (X,Y and Z distance) 
!> \param jmin minimum 
!> \param jmax maximum angular momentum ( N where t+u+v =< N )
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param ntuv number of TUV components
!> \param orderAP logical determining the ordering of WTUV and WJ000 [(ntuv,nPrim) or (nPrim,ntuv)]
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE wtuvRecurrence(WTUV,WJ000,sharedTUV,Rpq,jmin,jmax,&
     &                    nPrim,ntuv,orderAP,lupri,iprint)
implicit none
Integer,intent(in)       :: jmin,jmax,nPrim,ntuv,lupri,iprint
TYPE(TUVitem),intent(in) :: SharedTUV
Real(realk),intent(inout)   :: WTUV(:)
Real(realk),intent(in)   :: WJ000(:),Rpq(:,:)
Logical,intent(in)       :: orderAP
IF (orderAP) THEN
  CALL wtuvRecurrenceAP(WTUV,WJ000,sharedTUV,Rpq,jmin,jmax,&
     &                  nPrim,ntuv,lupri,iprint)
ELSE
  CALL wtuvRecurrencePA(WTUV,WJ000,sharedTUV,Rpq,jmin,jmax,&
     &                  nPrim,ntuv,lupri,iprint)
ENDIF
END SUBROUTINE wtuvRecurrence

!> \brief Implemented recurrence relation for hermite integrals Eq. (9.9.18-9.9.20) in the book assuming ntuv,nPrim ordering
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param WTUV the output of the recurrence 
!> \param WJ000 Eq. (9.9.14) in the book \f$ R^{n}_{000}(p,\textbf{R}_{PC}) = (-2p)^{n}F_{n}(pR_{PC}^{2}) \f$
!> \param sharedTUV a TUVitem which contains tabulated boys function and TUVindexing
!> \param Rpq distance between overlap distribution P and Q (X,Y and Z distance) 
!> \param jmin minimum 
!> \param jmax maximum angular momentum ( N where t+u+v =< N )
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param ntuv number of TUV components
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE wtuvRecurrenceAP(WTUV,WJ000,sharedTUV,Rpq,jmin,jmax,&
     &                      nPrim,ntuv,lupri,iprint)
use memory_handling
implicit none
Integer,intent(in):: jmin,jmax,nPrim,ntuv,lupri,iprint
TYPE(TUVitem),intent(in) :: SharedTUV
Real(realk),intent(in) :: WJ000(0:jmax,nPrim),Rpq(3,nPrim)
Real(realk),intent(inout) :: WTUV(ntuv,nPrim)
!
INTEGER,PARAMETER :: MAXJ = 50
Real(realk) :: Xtemp,Ytemp,Ztemp
Integer     :: n,j,t,u,v,ntuvfull,ituv,dir
Real(realk),pointer :: CUR(:,:),OLD(:,:),TEMP(:,:)
Real(realk) :: W100,W010,W001,W200,W110,W101,W020,W011,W002,WJ,WJ1,WJ2,WJ3
integer :: tm1,um1,vm1,tm2,um2,vm2,m1,ituvm,ituvm2,ioffp,k

IF (jmax.EQ.0) THEN
  CALL DCOPY(nPrim,WJ000,1,WTUV,1)
ELSE IF (jmax.EQ.1) THEN
  IF (jmin.EQ.1) THEN
    DO n=1,nPrim
      WJ1 = WJ000(1,n)
      WTUV(1,n) = Rpq(1,n)*WJ1
      WTUV(2,n) = Rpq(2,n)*WJ1
      WTUV(3,n) = Rpq(3,n)*WJ1
    ENDDO
  ELSE
    DO n=1,nPrim
      WTUV(1,n) = WJ000(0,n)
      WJ1 = WJ000(1,n)
      WTUV(2,n) = Rpq(1,n)*WJ1
      WTUV(3,n) = Rpq(2,n)*WJ1
      WTUV(4,n) = Rpq(3,n)*WJ1
    ENDDO
  ENDIF
ELSE IF (jmax.EQ.2) THEN
   IF (jmin.EQ.2) THEN
    DO n=1,nPrim
      Xtemp = Rpq(1,n)
      Ytemp = Rpq(2,n)
      Ztemp = Rpq(3,n)
      WJ1 = WJ000(1,n)
      WJ2 = WJ000(2,n)
      WTUV(1,n) = WJ1 + Xtemp*Xtemp*WJ2
      WTUV(2,n) =       Xtemp*Ytemp*WJ2
      WTUV(3,n) =       Xtemp*Ztemp*WJ2
      WTUV(4,n) = WJ1 + Ytemp*Ytemp*WJ2
      WTUV(5,n) =       Ytemp*Ztemp*WJ2
      WTUV(6,n) = WJ1 + Ztemp*Ztemp*WJ2
    ENDDO
   ELSE IF (jmin.EQ.1) THEN
    DO n=1,nPrim
      Xtemp = Rpq(1,n)
      Ytemp = Rpq(2,n)
      Ztemp = Rpq(3,n)
      WJ1 = WJ000(1,n)
      WJ2 = WJ000(2,n)
      WTUV(1,n) = Xtemp*WJ1
      WTUV(2,n) = Ytemp*WJ1
      WTUV(3,n) = Ztemp*WJ1
      WTUV(4,n) = WJ1 + Xtemp*Xtemp*WJ2
      WTUV(5,n) =       Xtemp*Ytemp*WJ2
      WTUV(6,n) =       Xtemp*Ztemp*WJ2
      WTUV(7,n) = WJ1 + Ytemp*Ytemp*WJ2
      WTUV(8,n) =       Ytemp*Ztemp*WJ2
      WTUV(9,n) = WJ1 + Ztemp*Ztemp*WJ2
    ENDDO
   ELSE
    DO n=1,nPrim
      Xtemp = Rpq(1,n)
      Ytemp = Rpq(2,n)
      Ztemp = Rpq(3,n)
      WJ  = WJ000(0,n)
      WJ1 = WJ000(1,n)
      WJ2 = WJ000(2,n)
      WTUV(1,n)  = WJ
      WTUV(2,n)  = Xtemp*WJ1
      WTUV(3,n)  = Ytemp*WJ1
      WTUV(4,n)  = Ztemp*WJ1
      WTUV(5,n)  = WJ1 + Xtemp*Xtemp*WJ2
      WTUV(6,n)  =       Xtemp*Ytemp*WJ2
      WTUV(7,n)  =       Xtemp*Ztemp*WJ2
      WTUV(8,n)  = WJ1 + Ytemp*Ytemp*WJ2
      WTUV(9,n)  =       Ytemp*Ztemp*WJ2
      WTUV(10,n) = WJ1 + Ztemp*Ztemp*WJ2
    ENDDO
   ENDIF
ELSE  ! J > 2
  ntuvfull = (jmax+1)*(jmax+2)*(jmax+3)/6
  ioffp    = jmin*(jmin+1)*(jmin+2)/6+1
  call mem_alloc(CUR,ntuvfull,nPrim)
  call mem_alloc(OLD,ntuvfull,nPrim)

  DO j=jmax-3,0,-1
    TEMP => CUR
    CUR  => OLD
    OLD  => TEMP
    DO n=1,nPrim
      Xtemp = Rpq(1,n)
      Ytemp = Rpq(2,n)
      Ztemp = Rpq(3,n)
      WJ   = WJ000(j,n)
      WJ1  = WJ000(j+1,n)
      WJ2  = WJ000(j+2,n)
      WJ3  = WJ000(j+3,n)
      W100 = Xtemp*WJ2
      W010 = Ytemp*WJ2
      W001 = Ztemp*WJ2
      W200 = WJ2 + Xtemp*Xtemp*WJ3
      W110 =       Xtemp*Ytemp*WJ3
      W101 =       Xtemp*Ztemp*WJ3
      W020 = WJ2 + Ytemp*Ytemp*WJ3
      W011 =       Ytemp*Ztemp*WJ3
      W002 = WJ2 + Ztemp*Ztemp*WJ3
      CUR(1,n)  = WJ                           !000
      CUR(2,n)  = Xtemp*WJ1                    !100
      CUR(3,n)  = Ytemp*WJ1                    !010
      CUR(4,n)  = Ztemp*WJ1                    !001
      CUR(5,n)  = WJ1 + Xtemp*W100             !200
      CUR(6,n)  =       Ytemp*W100             !110
      CUR(7,n)  =       Ztemp*W100             !101
      CUR(8,n)  = WJ1 + Ytemp*W010             !020
      CUR(9,n)  =       Ztemp*W010             !011
      CUR(10,n) = WJ1 + Ztemp*W001             !002
      CUR(11,n) = 2.0d0*W100 + Xtemp*W200      !300
      CUR(12,n) =       W010 + Xtemp*W110      !210
      CUR(13,n) =       W001 + Xtemp*W101      !201
      CUR(14,n) =       W100 + Ytemp*W110      !120
      CUR(15,n) =              Xtemp*W011      !111
      CUR(16,n) =       W100 + Ztemp*W101      !102
      CUR(17,n) = 2.0d0*W010 + Ytemp*W020      !030
      CUR(18,n) =       W001 + Ytemp*W011      !021
      CUR(19,n) =       W010 + Ztemp*W011      !012
      CUR(20,n) = 2.0d0*W001 + Ztemp*W002      !003
    ENDDO
    ituv = 20
    DO k=4,jmax-j
      DO t=k,0,-1
        DO u=k-t,0,-1
          v=k-t-u
          IF ((t.ge.u).AND.(t.ge.v)) THEN
            tm1 = t-1
            um1 = u
            vm1 = v
            tm2 = t-2
            um2 = u
            vm2 = v
            m1  = tm1
            dir = 1
          ELSEIF (u.ge.v) THEN
            tm1 = t
            um1 = u-1
            vm1 = v
            tm2 = t
            um2 = u-2
            vm2 = v
            m1  = um1
            dir = 2
          ELSE
            tm1 = t
            um1 = u
            vm1 = v-1
            tm2 = t
            um2 = u
            vm2 = v-2
            m1  = vm1
            dir = 3
          ENDIF
          ituv   = ituv + 1
          ituvm  = sharedTUV%tuvIndex(tm1,um1,vm1)
          ituvm2 = sharedTUV%tuvIndex(tm2,um2,vm2)
          DO n=1,nPrim
            CUR(ituv,n) = m1*OLD(ituvm2,n) + Rpq(dir,n)*OLD(ituvm,n)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DO n=1,nPrim
    CALL DCOPY(ntuv,CUR(ioffp,n),1,WTUV(1,n),1)
  ENDDO
  call mem_dealloc(CUR)
  call mem_dealloc(OLD)
ENDIF

END SUBROUTINE wtuvRecurrenceAP

!> \brief Implemented recurrence relation for hermite integrals Eq. (9.9.18-9.9.20) in the book assuming nPrim,ntuv ordering
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param WTUV the output of the recurrence 
!> \param WJ000 Eq. (9.9.14) in the book \f$ R^{n}_{000}(p,\textbf{R}_{PC}) = (-2p)^{n}F_{n}(pR_{PC}^{2}) \f$
!> \param sharedTUV a TUVitem which contains tabulated boys function and TUVindexing
!> \param Rpq distance between overlap distribution P and Q (X,Y and Z distance) 
!> \param jmin minimum 
!> \param jmax maximum angular momentum ( N where t+u+v =< N )
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param ntuv number of TUV components
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE wtuvRecurrencePA(WTUV,WJ000,sharedTUV,Rpq,jmin,jmax,&
     &                      nPrim,ntuv,lupri,iprint)
use memory_handling
implicit none
Integer,intent(in)           :: jmin,jmax,nPrim,ntuv,lupri,iprint
TYPE(TUVitem),intent(in)     :: SharedTUV
Real(realk),intent(in)       :: WJ000(0:jmax,nPrim),Rpq(nPrim,3)
Real(realk),intent(inout)    :: WTUV(nPrim,ntuv)
!
INTEGER,PARAMETER :: MAXJ = 50
Real(realk) :: Xtemp,Ytemp,Ztemp
Integer     :: n,j,t,u,v,ntuvfull,ituv,dir
Real(realk),pointer :: CUR(:,:),OLD(:,:),TEMP(:,:)
Real(realk) :: W100,W010,W001,W200,W110,W101,W020,W011,W002,WJ,WJ1,WJ2,WJ3
integer :: tm1,um1,vm1,tm2,um2,vm2,m1,ituvm,ituvm2,ioffp
integer :: idp1,idir,k

IF (jmax.EQ.0) THEN
  CALL DCOPY(nPrim,WJ000,1,WTUV,1)
ELSE IF (jmax.EQ.1) THEN
  IF (jmin.EQ.1) THEN
    DO idir=1,3
      DO n=1,nPrim
        WTUV(n,idir) = Rpq(n,idir)*WJ000(1,n)
      ENDDO
    ENDDO
  ELSE
    DO n=1,nPrim
      WTUV(n,1) = WJ000(0,n)
    ENDDO
    DO idir=1,3
      idp1 = idir+1
      DO n=1,nPrim
        WTUV(n,idp1) = Rpq(n,idir)*WJ000(1,n)
      ENDDO
    ENDDO
  ENDIF
ELSE IF (jmax.EQ.2) THEN
   IF (jmin.EQ.2) THEN
    DO n=1,nPrim
      Xtemp = Rpq(n,1)
      Ytemp = Rpq(n,2)
      Ztemp = Rpq(n,3)
      WJ1 = WJ000(1,n)
      WJ2 = WJ000(2,n)
      WTUV(n,1) = WJ1 + Xtemp*Xtemp*WJ2
      WTUV(n,2) =       Xtemp*Ytemp*WJ2
      WTUV(n,3) =       Xtemp*Ztemp*WJ2
      WTUV(n,4) = WJ1 + Ytemp*Ytemp*WJ2
      WTUV(n,5) =       Ytemp*Ztemp*WJ2
      WTUV(n,6) = WJ1 + Ztemp*Ztemp*WJ2
    ENDDO
   ELSE IF (jmin.EQ.1) THEN
    DO n=1,nPrim
      Xtemp = Rpq(n,1)
      Ytemp = Rpq(n,2)
      Ztemp = Rpq(n,3)
      WJ1 = WJ000(1,n)
      WJ2 = WJ000(2,n)
      WTUV(n,1) = Xtemp*WJ1
      WTUV(n,2) = Ytemp*WJ1
      WTUV(n,3) = Ztemp*WJ1
      WTUV(n,4) = WJ1 + Xtemp*Xtemp*WJ2
      WTUV(n,5) =       Xtemp*Ytemp*WJ2
      WTUV(n,6) =       Xtemp*Ztemp*WJ2
      WTUV(n,7) = WJ1 + Ytemp*Ytemp*WJ2
      WTUV(n,8) =       Ytemp*Ztemp*WJ2
      WTUV(n,9) = WJ1 + Ztemp*Ztemp*WJ2
    ENDDO
   ELSE
    DO n=1,nPrim
      Xtemp = Rpq(n,1)
      Ytemp = Rpq(n,2)
      Ztemp = Rpq(n,3)
      WJ  = WJ000(0,n)
      WJ1 = WJ000(1,n)
      WJ2 = WJ000(2,n)
      WTUV(n,1)  = WJ
      WTUV(n,2)  = Xtemp*WJ1
      WTUV(n,3)  = Ytemp*WJ1
      WTUV(n,4)  = Ztemp*WJ1
      WTUV(n,5)  = WJ1 + Xtemp*Xtemp*WJ2
      WTUV(n,6)  =       Xtemp*Ytemp*WJ2
      WTUV(n,7)  =       Xtemp*Ztemp*WJ2
      WTUV(n,8)  = WJ1 + Ytemp*Ytemp*WJ2
      WTUV(n,9)  =       Ytemp*Ztemp*WJ2
      WTUV(n,10) = WJ1 + Ztemp*Ztemp*WJ2
    ENDDO
   ENDIF
ELSE  ! J > 2
  ntuvfull = (jmax+1)*(jmax+2)*(jmax+3)/6
  ioffp    = jmin*(jmin+1)*(jmin+2)/6+1
  call mem_alloc(CUR,nPrim,ntuvfull)
  call mem_alloc(OLD,nPrim,ntuvfull)

  DO j=jmax-3,0,-1
    TEMP => CUR
    CUR  => OLD
    OLD  => TEMP
    DO n=1,nPrim
      Xtemp = Rpq(n,1)
      Ytemp = Rpq(n,2)
      Ztemp = Rpq(n,3)
      WJ   = WJ000(j,n)
      WJ1  = WJ000(j+1,n)
      WJ2  = WJ000(j+2,n)
      WJ3  = WJ000(j+3,n)
      W100 = Xtemp*WJ2
      W010 = Ytemp*WJ2
      W001 = Ztemp*WJ2
      W200 = WJ2 + Xtemp*Xtemp*WJ3
      W110 =       Xtemp*Ytemp*WJ3
      W101 =       Xtemp*Ztemp*WJ3
      W020 = WJ2 + Ytemp*Ytemp*WJ3
      W011 =       Ytemp*Ztemp*WJ3
      W002 = WJ2 + Ztemp*Ztemp*WJ3
      CUR(n,1)  = WJ                           !000
      CUR(n,2)  = Xtemp*WJ1                    !100
      CUR(n,3)  = Ytemp*WJ1                    !010
      CUR(n,4)  = Ztemp*WJ1                    !001
      CUR(n,5)  = WJ1 + Xtemp*W100             !200
      CUR(n,6)  =       Ytemp*W100             !110
      CUR(n,7)  =       Ztemp*W100             !101
      CUR(n,8)  = WJ1 + Ytemp*W010             !020
      CUR(n,9)  =       Ztemp*W010             !011
      CUR(n,10) = WJ1 + Ztemp*W001             !002
      CUR(n,11) = 2.0d0*W100 + Xtemp*W200      !300
      CUR(n,12) =       W010 + Xtemp*W110      !210
      CUR(n,13) =       W001 + Xtemp*W101      !201
      CUR(n,14) =       W100 + Ytemp*W110      !120
      CUR(n,15) =              Xtemp*W011      !111
      CUR(n,16) =       W100 + Ztemp*W101      !102
      CUR(n,17) = 2.0d0*W010 + Ytemp*W020      !030
      CUR(n,18) =       W001 + Ytemp*W011      !021
      CUR(n,19) =       W010 + Ztemp*W011      !012
      CUR(n,20) = 2.0d0*W001 + Ztemp*W002      !003
    ENDDO
    ituv = 20
    DO k=4,jmax-j
      DO t=k,0,-1
        DO u=k-t,0,-1
          v=k-t-u
          IF ((t.ge.u).AND.(t.ge.v)) THEN
            tm1 = t-1
            um1 = u
            vm1 = v
            tm2 = t-2
            um2 = u
            vm2 = v
            m1  = tm1
            dir = 1
          ELSEIF (u.ge.v) THEN
            tm1 = t
            um1 = u-1
            vm1 = v
            tm2 = t
            um2 = u-2
            vm2 = v
            m1  = um1
            dir = 2
          ELSE
            tm1 = t
            um1 = u
            vm1 = v-1
            tm2 = t
            um2 = u
            vm2 = v-2
            m1  = vm1
            dir = 3
          ENDIF
          ituv   = ituv + 1
          ituvm  = sharedTUV%tuvIndex(tm1,um1,vm1)
          ituvm2 = sharedTUV%tuvIndex(tm2,um2,vm2)
          DO n=1,nPrim
            CUR(n,ituv) = m1*OLD(n,ituvm2) + Rpq(n,dir)*OLD(n,ituvm)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CALL DCOPY(ntuv*nPrim,CUR(1,ioffp),1,WTUV,1)
  call mem_dealloc(CUR)
  call mem_dealloc(OLD)
ENDIF

END SUBROUTINE wtuvRecurrencePA

!!$SUBROUTINE RECURRENCE(CUR,OLD,JVAL,RJ000,Xpq,Ypq,Zpq,JMAX,&
!!$    &                  NPrim,NTUV,TUVindex,zeroX,zeroY,zeroZ)
!!$ IMPLICIT NONE
!!$ real(realk)             :: CUR(nPrim,NTUV)
!!$ real(realk)             :: OLD(nPrim,NTUV)
!!$ real(realk)             :: RJ000(nPrim,0:Jmax)
!!$ INTEGER                 :: J,K,T,U,V,TUV,JVAL,JMAX,NTUV
!!$ INTEGER                 :: nPrim
!!$ INTEGER,PARAMETER       :: MAXJ = 50
!!$ INTEGER                 :: TUVindex(0:MAXJ,0:MAXJ,0:MAXJ)
!!$ real(realk)             :: Xpq(nPrim),Ypq(nPrim),Zpq(nPrim)
!!$ INTEGER                 :: MAXT,MAXU,MAXV,M1T,M2T,M1U,M2U
!!$ INTEGER                 :: M1V,M2V,IUMAX,I
!!$ INTEGER                 :: zeroX,zeroY,zeroZ,M,MP1
!!$ LOGICAL                 :: NOTSAMEX,NOTSAMEY,NOTSAMEZ
!!$ Real(realk)             :: TMIN1,UMIN1,VMIN1
!!$
!!$ NOTSAMEX=zeroX .EQ. 1
!!$ NOTSAMEY=zeroY .EQ. 1
!!$ NOTSAMEZ=zeroZ .EQ. 1
!!$
!!$ IF (JVAL .EQ. 1) THEN
!!$    !     JVAL = 1
!!$    !     ========
!!$    CALL DCOPY(NPrim,RJ000(1,JMAX-1),1,CUR(1,1),1)
!!$    !LOOP CAN BE PREVENTED IF SAME CENTER
!!$    IF(NOTSAMEX) THEN
!!$       !      DO I = 1, Nprim
!!$       !         CUR(I,2) = Xpq(I)*RJ000(I,JMAX)
!!$       !      ENDDO
!!$
!!$       M = MOD(NPrim,4)
!!$       IF(M .EQ. 0) THEN
!!$          DO I = 1,NPrim,4
!!$             CUR(I,2)   = Xpq(I)*RJ000(I,JMAX)
!!$             CUR(I+1,2) = Xpq(I+1)*RJ000(I+1,JMAX)
!!$             CUR(I+2,2) = Xpq(I+2)*RJ000(I+2,JMAX)
!!$             CUR(I+3,2) = Xpq(I+3)*RJ000(I+3,JMAX) 
!!$          ENDDO
!!$       ELSE
!!$          DO I = 1,M
!!$             CUR(I,2) = Xpq(I)*RJ000(I,JMAX) 
!!$          ENDDO
!!$          IF(NPrim .GT. 4)THEN
!!$             MP1=M+1
!!$             DO I = MP1,NPrim,4
!!$                CUR(I,2)   = Xpq(I)*RJ000(I,JMAX) 
!!$                CUR(I+1,2) = Xpq(I+1)*RJ000(I+1,JMAX) 
!!$                CUR(I+2,2) = Xpq(I+2)*RJ000(I+2,JMAX) 
!!$                CUR(I+3,2) = Xpq(I+3)*RJ000(I+3,JMAX) 
!!$             ENDDO
!!$          ENDIF
!!$       ENDIF
!!$    ENDIF
!!$    IF(NOTSAMEY) THEN
!!$       !      DO I = 1, Nprim
!!$       !         CUR(I,3) = Ypq(I)*RJ000(I,JMAX)
!!$       !      ENDDO
!!$       M = MOD(NPrim,4)
!!$       IF(M .EQ. 0) THEN
!!$          DO I = 1,NPrim,4
!!$             CUR(I,3)   = Ypq(I)*RJ000(I,JMAX)
!!$             CUR(I+1,3) = Ypq(I+1)*RJ000(I+1,JMAX)
!!$             CUR(I+2,3) = Ypq(I+2)*RJ000(I+2,JMAX)
!!$             CUR(I+3,3) = Ypq(I+3)*RJ000(I+3,JMAX) 
!!$          ENDDO
!!$       ELSE
!!$          DO I = 1,M
!!$             CUR(I,3) = Ypq(I)*RJ000(I,JMAX) 
!!$          ENDDO
!!$          IF(NPrim .GT. 4)THEN
!!$             MP1=M+1
!!$             DO I = MP1,NPrim,4
!!$                CUR(I,3)   = Ypq(I)*RJ000(I,JMAX) 
!!$                CUR(I+1,3) = Ypq(I+1)*RJ000(I+1,JMAX) 
!!$                CUR(I+2,3) = Ypq(I+2)*RJ000(I+2,JMAX) 
!!$                CUR(I+3,3) = Ypq(I+3)*RJ000(I+3,JMAX) 
!!$             ENDDO
!!$          ENDIF
!!$       ENDIF
!!$    ENDIF
!!$    IF(NOTSAMEZ) THEN
!!$       !      DO I = 1, Nprim
!!$       !         CUR(I,4) = Zpq(I)*RJ000(I,JMAX)
!!$       !      ENDDO
!!$       M = MOD(NPrim,4)
!!$       IF(M .EQ. 0) THEN
!!$          DO I = 1,NPrim,4
!!$             CUR(I,4)   = Zpq(I)*RJ000(I,JMAX)
!!$             CUR(I+1,4) = Zpq(I+1)*RJ000(I+1,JMAX)
!!$             CUR(I+2,4) = Zpq(I+2)*RJ000(I+2,JMAX)
!!$             CUR(I+3,4) = Zpq(I+3)*RJ000(I+3,JMAX) 
!!$          ENDDO
!!$       ELSE
!!$          DO I = 1,M
!!$             CUR(I,4) = Zpq(I)*RJ000(I,JMAX) 
!!$          ENDDO
!!$          IF(NPrim .GT. 4)THEN
!!$             MP1=M+1
!!$             DO I = MP1,NPrim,4
!!$                CUR(I,4)   = Zpq(I)*RJ000(I,JMAX) 
!!$                CUR(I+1,4) = Zpq(I+1)*RJ000(I+1,JMAX) 
!!$                CUR(I+2,4) = Zpq(I+2)*RJ000(I+2,JMAX) 
!!$                CUR(I+3,4) = Zpq(I+3)*RJ000(I+3,JMAX) 
!!$             ENDDO
!!$          ENDIF
!!$       ENDIF
!!$    ENDIF
!!$ ELSE
!!$
!!$    !     JVAL > 1
!!$    !     ========
!!$
!!$    MAXT   = JMAX
!!$    MAXU   = JMAX
!!$    MAXV   = JMAX
!!$
!!$    !        R(0,0,0)
!!$
!!$    CALL DCOPY(Nprim,RJ000(1,JMAX-JVAL),1,CUR,1)
!!$
!!$    !        R(T,0,0)
!!$
!!$    IF (NOTSAMEX) THEN !NOT SAME CENTER
!!$       DO I = 1, NPrim
!!$          CUR(I,2) = Xpq(I)*OLD(I,1)
!!$       ENDDO
!!$       DO T = 2, MIN(MAXT,JVAL)
!!$          !         TMIN1 = DFLOAT(T - 1)
!!$          TMIN1 = T - 1
!!$          TUV   = TUVINDEX(T  ,0,0)
!!$          M1T   = TUVINDEX(T-1,0,0)
!!$          M2T   = TUVINDEX(T-2,0,0)
!!$          DO I = 1, Nprim
!!$             CUR(I,TUV) = Xpq(I)*OLD(I,M1T) + TMIN1*OLD(I,M2T)
!!$          ENDDO
!!$       ENDDO
!!$    ELSE
!!$       DO T = 2, MIN(MAXT,JVAL), 2
!!$          !        TMIN1 = DFLOAT(T - 1)
!!$          TMIN1 = T - 1
!!$          TUV   = TUVINDEX(T,0,0)
!!$          M2T   = TUVINDEX(T-2,0,0)
!!$
!!$          !         CALL DCOPY(Nprim,OLD(1,M2T),1,CUR(1,TUV),1)
!!$          !         CALL DSCAL(Nprim,TMIN1,CUR(1,TUV),1) 
!!$          !         DO I = 1, NPrim
!!$          !            CUR(I,TUV) = TMIN1*OLD(I,M2T)
!!$          !         ENDDO
!!$
!!$          M = MOD(NPrim,4)
!!$          IF(M .EQ. 0) THEN
!!$             DO I = 1,NPrim,4
!!$                CUR(I,TUV)   = TMIN1*OLD(I,M2T)   
!!$                CUR(I+1,TUV) = TMIN1*OLD(I+1,M2T)   
!!$                CUR(I+2,TUV) = TMIN1*OLD(I+2,M2T)   
!!$                CUR(I+3,TUV) = TMIN1*OLD(I+3,M2T)   
!!$             ENDDO
!!$          ELSE
!!$             DO I = 1,M
!!$                CUR(I,TUV) = TMIN1*OLD(I,M2T)
!!$             ENDDO
!!$             IF(NPrim .GT. 4)THEN
!!$                MP1=M+1
!!$                DO I = MP1,NPrim,4
!!$                   CUR(I,TUV)   = TMIN1*OLD(I,M2T)   
!!$                   CUR(I+1,TUV) = TMIN1*OLD(I+1,M2T)   
!!$                   CUR(I+2,TUV) = TMIN1*OLD(I+2,M2T)   
!!$                   CUR(I+3,TUV) = TMIN1*OLD(I+3,M2T)   
!!$                ENDDO
!!$             ENDIF
!!$          ENDIF
!!$
!!$       ENDDO
!!$    END IF
!!$
!!$    !        R(T,U,0)
!!$
!!$    IF (NOTSAMEY) THEN !NOT SAME CENTER
!!$       DO T = 0, MIN(MAXT,JVAL - 1), zeroX
!!$          TUV = TUVINDEX(T,1,0)
!!$          M1U = TUVINDEX(T,0,0)
!!$          DO I = 1, NPrim
!!$             CUR(I,TUV) = Ypq(I)*OLD(I,M1U)
!!$          ENDDO
!!$       ENDDO
!!$       DO U = 2, MIN(MAXU,JVAL)
!!$          !        UMIN1  = DFLOAT(U - 1)
!!$          UMIN1  = U - 1
!!$          DO T = 0, MIN(MAXT,JVAL - U), zeroX
!!$             TUV = TUVINDEX(T,U  ,0)
!!$             M1U = TUVINDEX(T,U-1,0)
!!$             M2U = TUVINDEX(T,U-2,0)
!!$             DO I = 1, Nprim
!!$                CUR(I,TUV) = Ypq(I)*OLD(I,M1U) + UMIN1*OLD(I,M2U)
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$    ELSE
!!$       DO U = 2, MIN(MAXU,JVAL), 2
!!$          !        UMIN1  = DFLOAT(U - 1)
!!$          UMIN1  = U - 1
!!$          DO T = 0, MIN(MAXT,JVAL - U), zeroX
!!$             TUV = TUVINDEX(T,U  ,0)
!!$             M2U = TUVINDEX(T,U-2,0)
!!$             !            CALL DCOPY(Nprim,OLD(1,M2U),1,CUR(1,TUV),1)
!!$             !            CALL DSCAL(nPrim,UMIN1,CUR(1,TUV),1) 
!!$             !            DO I = 1, NPrim
!!$             !               CUR(I,TUV) = UMIN1*OLD(I,M2U)
!!$             !            ENDDO
!!$
!!$             M = MOD(NPrim,4)
!!$             IF(M .EQ. 0) THEN
!!$                DO I = 1,NPrim,4
!!$                   CUR(I,TUV)   = UMIN1*OLD(I,M2U)   
!!$                   CUR(I+1,TUV) = UMIN1*OLD(I+1,M2U)   
!!$                   CUR(I+2,TUV) = UMIN1*OLD(I+2,M2U)   
!!$                   CUR(I+3,TUV) = UMIN1*OLD(I+3,M2U)   
!!$                ENDDO
!!$             ELSE
!!$                DO I = 1,M
!!$                   CUR(I,TUV) = UMIN1*OLD(I,M2U)
!!$                ENDDO
!!$                IF(NPrim .GT. 4)THEN
!!$                   MP1=M+1
!!$                   DO I = MP1,NPrim,4
!!$                      CUR(I,TUV)   = UMIN1*OLD(I,M2U)   
!!$                      CUR(I+1,TUV) = UMIN1*OLD(I+1,M2U)   
!!$                      CUR(I+2,TUV) = UMIN1*OLD(I+2,M2U)   
!!$                      CUR(I+3,TUV) = UMIN1*OLD(I+3,M2U)   
!!$                   ENDDO
!!$                ENDIF
!!$             ENDIF
!!$
!!$          ENDDO
!!$       ENDDO
!!$    END IF
!!$
!!$    !        R(T,U,V)
!!$
!!$    IF (NOTSAMEZ) THEN !IF NOT SAME CENTER
!!$       IUMAX  = JVAL - 1
!!$       DO U = 0, MIN(MAXU,IUMAX), zeroY
!!$          DO T = 0, MIN(MAXT,IUMAX - U), zeroX
!!$             TUV = TUVINDEX(T,U,1)
!!$             M1V = TUVINDEX(T,U,0)
!!$             DO I = 1, NPrim
!!$                CUR(I,TUV) = Zpq(I)*OLD(I,M1V)
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$       DO V = 2, MIN(MAXV,JVAL)
!!$          !        VMIN1  = DFLOAT(V - 1)
!!$          VMIN1  = V - 1
!!$          IUMAX  = JVAL - V
!!$          DO U = 0, MIN(MAXU,IUMAX), zeroY
!!$             DO T = 0, MIN(MAXT,IUMAX - U), zeroX
!!$                TUV = TUVINDEX(T,U,V  )
!!$                M1V = TUVINDEX(T,U,V-1)
!!$                M2V = TUVINDEX(T,U,V-2)
!!$                DO I = 1, NPrim
!!$                   CUR(I,TUV) = Zpq(I)*OLD(I,M1V)+VMIN1*OLD(I,M2V)
!!$                ENDDO
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$    ELSE
!!$       DO V = 2, MIN(MAXV,JVAL), 2
!!$          !        VMIN1  = DFLOAT(V - 1)
!!$          VMIN1  = V - 1
!!$          IUMAX  = JVAL - V
!!$          DO U = 0, MIN(MAXU,IUMAX), zeroY
!!$             DO T = 0, MIN(MAXT,IUMAX - U), zeroX
!!$                TUV = TUVINDEX(T,U,V  )
!!$                M2V = TUVINDEX(T,U,V-2)
!!$                !               CALL DCOPY(Nprim,OLD(1,M2V),1,CUR(1,TUV),1)
!!$                !               CALL DSCAL(nPrim,VMIN1,CUR(1,TUV),1) 
!!$                !               DO I = 1, Nprim
!!$                !                  CUR(I,TUV) = VMIN1*OLD(I,M2V)
!!$                !               ENDDO
!!$
!!$
!!$                M = MOD(NPrim,4)
!!$                IF(M .EQ. 0) THEN
!!$                   DO I = 1,NPrim,4
!!$                      CUR(I,TUV)   = VMIN1*OLD(I,M2V)   
!!$                      CUR(I+1,TUV) = VMIN1*OLD(I+1,M2V)   
!!$                      CUR(I+2,TUV) = VMIN1*OLD(I+2,M2V)   
!!$                      CUR(I+3,TUV) = VMIN1*OLD(I+3,M2V)   
!!$                   ENDDO
!!$                ELSE
!!$                   DO I = 1,M
!!$                      CUR(I,TUV) = VMIN1*OLD(I,M2V)
!!$                   ENDDO
!!$                   IF(NPrim .GT. 4)THEN
!!$                      MP1=M+1
!!$                      DO I = MP1,NPrim,4
!!$                         CUR(I,TUV)   = VMIN1*OLD(I,M2V)   
!!$                         CUR(I+1,TUV) = VMIN1*OLD(I+1,M2V)   
!!$                         CUR(I+2,TUV) = VMIN1*OLD(I+2,M2V)   
!!$                         CUR(I+3,TUV) = VMIN1*OLD(I+3,M2V)   
!!$                      ENDDO
!!$                   ENDIF
!!$                ENDIF
!!$
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$    END IF
!!$ END IF
!!$
!!$END SUBROUTINE RECURRENCE

!> \brief wrapper to Overlap hermite integral
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param SJ000 \f$ S^{0:jmax}_{000} \f$
!> \param PQ Integrand (info about the overlap distributions P and Q)
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE buildSJ000(SJ000,PQ,nPrim,Jmax,LUPRI,IPRINT)
 IMPLICIT NONE
 TYPE(Integrand),intent(in)  :: PQ
 real(realk),intent(inout)   :: SJ000(0:Jmax,nPrim)
 INTEGER,intent(in)          :: nPrim,LUPRI,IPRINT,Jmax
! INTEGER  :: J,I
 CALL buildSJ000_OP(SJ000,Nprim,Jmax,PQ%reducedExponents&
      &,PQ%exponents,PQ%squaredDistance)

! IF(IPRINT .GT. 25) THEN
!    DO I=1,NPrim
!       DO J=0,Jmax
!          WRITE(LUPRI,'(2X,A6,I4,A1,I2,A2,F16.9)')'SJ000(',I,',',J,')=',SJ000(J,I)
!       ENDDO
!    ENDDO
! ENDIF
END SUBROUTINE buildSJ000

!> \brief Overlap hermite integral \f[ S^{j}_{000} = \left(-2 \alpha \right)^{j} \left( \frac{\gamma}{2\pi} \right)^{\frac{3}{2}} \exp(-\gamma R_{pq}^{2}) \f]
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param SJ000 \f$ S^{0:jmax}_{000} \f$
!> \param PQ Integrand (info about the overlap distributions P and Q)
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE buildSJ000_OP(SJ000,nPrim,Jmax,Alpha,P,R2)
use memory_handling
 IMPLICIT NONE
 real(realk),intent(inout) :: SJ000(0:Jmax,nPrim)
 INTEGER,intent(in)        :: nPrim,Jmax
 REAL(REALK),intent(in)    :: ALPHA(nPrim),P(nPrim),R2(nprim)
 !
 Real(realk), parameter :: PI=3.14159265358979323846D0
 Real(realk), parameter :: OneHalf=1.5d0, Two = 2.0d0
 REAL(REALK),pointer    :: TEMP(:),TEMP2(:)
 INTEGER                :: J,I
 call mem_alloc(TEMP,nprim)
 call mem_alloc(TEMP2,nprim)
 DO I=1,NPrim
    TEMP(I)=Alpha(I)*R2(I)
 ENDDO
 DO I=1,NPrim
    TEMP2(I)=EXP(-TEMP(I))
 ENDDO
 DO I=1,NPrim
    DO J=0,Jmax
       SJ000(J,I)=((-Two*Alpha(I))**J)&
            &*((PI/P(I))**OneHalf)*TEMP2(I)
    ENDDO
 ENDDO
 call mem_dealloc(TEMP)
 call mem_dealloc(TEMP2)
END SUBROUTINE buildSJ000_OP

!> \brief wrapper to 2 electron Coulomb hermite integral \f$ R^{n}_{000}(\alpha,\textbf{R}_{PQ}) =  \left(-2 \alpha \right)^{n}  F_{n}(\alpha R_{PQ}^{2}) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param RJ000 \f$ R^{0:jmax}_{000} \f$
!> \param PQ Integrand (info about the overlap distributions P and Q)
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param integral containing tabulated boys function \f$ F_{n} \f$ 
!> \param HIGH logical to switch between high accuracy and low accuracy versions
SUBROUTINE buildRJ000(RJ000,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,HIGH)
 IMPLICIT NONE
 INTEGER,intent(in)         :: nPrim,LUPRI,IPRINT,Jmax
 REAL(REALK),intent(inout)  :: RJ000(0:Jmax,nPrim)
 TYPE(Integrand),intent(in) :: PQ
 TYPE(integralitem),intent(in) :: integral
 LOGICAL,intent(in)  :: HIGH

 IF (HIGH) THEN
    CALL buildRJ000_OP_HA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
         &PQ%reducedExponents,PQ%squaredDistance,&
         &PQ%integralPrefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)
 ELSE
    CALL buildRJ000_OP_LA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
         &PQ%reducedExponents,PQ%squaredDistance,&
         &PQ%integralPrefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)
 ENDIF

END SUBROUTINE buildRJ000

!> \brief wrapper to 1 electron hermite Coulomb integral \f$ R^{n}_{000}(p,\textbf{R}_{PC}) =  \left(-2 p \right)^{n}  F_{n}(p R_{PC}^{2}) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param RJ000 \f$ R^{0:jmax}_{000} \f$
!> \param PQ Integrand (info about the overlap distributions P and Q)
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param integral containing tabulated boys function \f$ F_{n} \f$ 
!> \param HIGH logical to switch between high accuracy and low accuracy versions
SUBROUTINE buildNuclearRJ000(RJ000,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,HIGH)
 IMPLICIT NONE
 INTEGER,intent(in)        :: nPrim,LUPRI,IPRINT,Jmax
 REAL(REALK),intent(inout) :: RJ000(0:Jmax,nPrim)
 TYPE(Integrand),intent(in):: PQ
 TYPE(integralitem),intent(in) :: integral
 LOGICAL,intent(in) :: HIGH

 IF (HIGH) THEN
    CALL buildRJ000_OP_HA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
         &PQ%Exponents,PQ%squaredDistance,&
         &PQ%integralPrefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)
 ELSE
    CALL buildRJ000_OP_LA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
         &PQ%Exponents,PQ%squaredDistance,&
         &PQ%integralPrefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2)
 ENDIF

END SUBROUTINE buildNuclearRJ000

!> \brief Low accuracy version of electron hermite Coulomb integral \f$ R^{n}_{000}(\alpha,\textbf{R}_{PQ}) =  \left(-2 \alpha \right)^{n}  F_{n}(\alpha R_{PQ}^{2}) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine, B. Jansik  \endlatexonly
!> \date 2009-02-05
!> \param RJ000 \f$ R^{0:jmax}_{000} \f$
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param alpha reduced exponents
!> \param squared distance
!> \param Prefactor (Coulomb: \f$ \frac{2 \pi^{\frac{5}{2}}}{pq\sqrt{p+q}} \f$, Nuclear att: \f$ \frac{2 \pi}{p}\f$)
!> \param TABFJW tabulated boys function \f$ F_{n} \f$ 
!> \param nTABFJW1 dimension 1 of TABFJW
!> \param nTABFJW2 dimension 2 of TABFJW
SUBROUTINE buildRJ000_OP_LA(RJ000,nPrim,Jmax,LUPRI,IPRINT,alpha,R2,Prefactor,TABFJW,nTABFJW1,nTABFJW2)
 IMPLICIT NONE
 INTEGER,intent(in)        :: nPrim,Jmax,Lupri,Iprint,nTABFJW1,nTABFJW2
 REAL(REALK),intent(inout) :: RJ000(0:Jmax,nPrim)
 REAL(REALK),intent(in)    :: alpha(nPrim),R2(nprim),Prefactor(nprim)
 REAL(REALK),intent(in)    :: TABFJW(0:nTABFJW1,0:nTABFJW2)
!
 INTEGER         :: I
 REAL(REALK)     :: D2JP36,WVAL!,WVALS(Nprim,3),WVALU(Nprim)
 REAL(REALK),PARAMETER :: HALF =0.5d0,D1=1.d0,D2 = 2.D0, D4 = 4.D0, D10=10.d0,D100=100d0
 Real(realk),parameter :: D12 = 12.D0, TENTH = 0.01D0
 REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6.D0, COEF4 = D1/24.D0
 REAL(REALK), PARAMETER :: COEF5 = - D1/120.D0, COEF6 = D1/720.D0
 Integer :: IPNT,J
 Real(realk) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA
 REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092D0
 REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686D0
 REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909D0
 REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346D0
 Real(realk), parameter :: PI=3.14159265358979323846D0
 REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730D00
 REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
 REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
 Real(realk) :: R!W2,W3,W4,W5,W6
 REAL(REALK), PARAMETER :: SMALL = 1.D-15!0.000001D0

 D2JP36 = 2*JMAX + 36

 DO I = 1, Nprim
    WVAL = alpha(I)*R2(I)
    !  0 < WVAL < 0.000001
    IF (WVAL .LT. SMALL) THEN         
       RJ000(0,I) = D1
       DO J=1,JMAX
          RJ000(J,I)= D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT
       ENDDO
       !  0 < WVAL < 12 
    ELSE IF (WVAL .LT. D12) THEN
       IPNT = NINT(D100*WVAL)
       WDIFF = WVAL - TENTH*IPNT
       !    W2    = WDIFF*WDIFF
       !    W3    = W2*WDIFF
       !    W4    = W3*WDIFF
       !    W5    = W4*WDIFF
       !    W6    = W5*WDIFF

       !    W2    = W2*COEF2
       !    W3    = W3*COEF3
       !    W4    = W4*COEF4
       !    W5    = W5*COEF5
       !    W6    = W6*COEF6

       DO J=0,JMAX
          R = TABFJW(J,IPNT)
          R = R -TABFJW(J+1,IPNT)*WDIFF
          !    R = R + TABFJW(J+2,IPNT)*W2
          !    R = R + TABFJW(J+3,IPNT)*W3
          !    R = R + TABFJW(J+4,IPNT)*W4
          !    R = R + TABFJW(J+5,IPNT)*W5
          !    R = R + TABFJW(J+6,IPNT)*W6
          RJ000(J,I) = R
       ENDDO

       !  12 < WVAL <= (2J+36) 
    ELSE IF (WVAL.LE.D2JP36) THEN
       REXPW = HALF*EXP(-WVAL)
       RWVAL = D1/WVAL
       GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
       RJ000(0,I) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
       DO J=1,JMAX
          RJ000(J,I) = RWVAL*((J - HALF)*RJ000(J-1,I)-REXPW)
       ENDDO
       !  (2J+36) < WVAL 
    ELSE
       RWVAL = PID4/WVAL
       RJ000(0,I) = SQRT(RWVAL)
       RWVAL = RWVAL*PID4I
       DO J = 1, JMAX
          RJ000(J,I) = RWVAL*(J - HALF)*RJ000(J-1,I)
       ENDDO
    END IF
 ENDDO

 ! Scaling
 DO I=1,nPrim
    PREF = Prefactor(I)
    RJ000(0,I) = PREF*RJ000(0,I)
    IF (jmax.GT.0) THEN
       D2MALPHA = -2*alpha(I)
       DO j=1,jmax
          PREF = PREF*D2MALPHA
          RJ000(J,I) = PREF*RJ000(J,I)
       ENDDO
    ENDIF
 ENDDO

END SUBROUTINE buildRJ000_OP_LA

!> \brief High accuracy version of electron hermite Coulomb integral \f$ R^{n}_{000}(\alpha,\textbf{R}_{PQ}) =  \left(-2 \alpha \right)^{n}  F_{n}(\alpha R_{PQ}^{2}) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine, B. Jansik  \endlatexonly
!> \date 2009-02-05
!> \param RJ000 \f$ R^{0:jmax}_{000} \f$
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param alpha reduced exponents
!> \param squared distance
!> \param Prefactor (Coulomb: \f$ \frac{2 \pi^{\frac{5}{2}}}{pq\sqrt{p+q}} \f$, Nuclear att: \f$ \frac{2 \pi}{p}\f$)
!> \param TABFJW tabulated boys function \f$ F_{n} \f$ 
!> \param nTABFJW1 dimension 1 of TABFJW
!> \param nTABFJW2 dimension 2 of TABFJW
SUBROUTINE buildRJ000_OP_HA(RJ000,nPrim,Jmax,LUPRI,IPRINT,alpha,R2,Prefactor,TABFJW,nTABFJW1,nTABFJW2)
 IMPLICIT NONE
 INTEGER,intent(in)         :: nPrim,Jmax,Lupri,Iprint,nTABFJW1,nTABFJW2
 REAL(REALK),intent(inout)  :: RJ000(0:Jmax,nPrim)
 REAL(REALK),intent(in)     :: alpha(nPrim),R2(nprim),Prefactor(nprim)
 REAL(REALK),intent(in)     :: TABFJW(0:nTABFJW1,0:nTABFJW2)
!
 INTEGER         :: I
 REAL(REALK)     :: D2JP36,WVAL!,WVALS(Nprim,3),WVALU(Nprim)
 REAL(REALK),PARAMETER :: HALF =0.5d0,D1=1.d0,D2 = 2.D0, D4 = 4.D0, D10=10.d0,D100=100d0
 Real(realk),parameter :: D12 = 12.D0, TENTH = 0.01D0
 REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6.D0, COEF4 = D1/24.D0
 REAL(REALK), PARAMETER :: COEF5 = - D1/120.D0, COEF6 = D1/720.D0
 Integer :: IPNT,J
 Real(realk) :: WDIFF,RWVAL,REXPW,GVAL,PREF,D2MALPHA
 REAL(REALK), PARAMETER :: GFAC0 =  D2*0.4999489092D0
 REAL(REALK), PARAMETER :: GFAC1 = -D2*0.2473631686D0
 REAL(REALK), PARAMETER :: GFAC2 =  D2*0.321180909D0
 REAL(REALK), PARAMETER :: GFAC3 = -D2*0.3811559346D0
 Real(realk), parameter :: PI=3.14159265358979323846D0
 REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730D00
 REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
 REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI
 Real(realk) :: W2,W3,R!W4,W5,W6,
! LOGICAL :: maxJgt0
 REAL(REALK), PARAMETER :: SMALL = 1.D-15!0.000001D0

 D2JP36 = 2*JMAX + 36

 DO I = 1, Nprim
    WVAL = alpha(I)*R2(I)
    !  0 < WVAL < 0.000001
    IF (WVAL .LT. SMALL) THEN         
       RJ000(0,I) = D1
       DO J=1,JMAX
          RJ000(J,I)= D1/(2*J + 1) !THE BOYS FUNCTION FOR ZERO ARGUMENT
       ENDDO
       !  0 < WVAL < 12 
    ELSE IF (WVAL .LT. D12) THEN
       IPNT = NINT(D100*WVAL)
       WDIFF = WVAL - TENTH*IPNT
       W2    = WDIFF*WDIFF
       W3    = W2*WDIFF
       !    W4    = W3*WDIFF
       !    W5    = W4*WDIFF
       !    W6    = W5*WDIFF

       W2    = W2*COEF2
       W3    = W3*COEF3
       !    W4    = W4*COEF4
       !    W5    = W5*COEF5
       !    W6    = W6*COEF6

       DO J=0,JMAX
          R = TABFJW(J,IPNT)
          R = R -TABFJW(J+1,IPNT)*WDIFF
          R = R + TABFJW(J+2,IPNT)*W2
          R = R + TABFJW(J+3,IPNT)*W3
          !    R = R + TABFJW(J+4,IPNT)*W4
          !    R = R + TABFJW(J+5,IPNT)*W5
          !    R = R + TABFJW(J+6,IPNT)*W6
          RJ000(J,I) = R
       ENDDO

       !  12 < WVAL <= (2J+36) 
    ELSE IF (WVAL.LE.D2JP36) THEN
       REXPW = HALF*EXP(-WVAL)
       RWVAL = D1/WVAL
       GVAL  = GFAC0 + RWVAL*(GFAC1 + RWVAL*(GFAC2 + RWVAL*GFAC3))
       RJ000(0,I) = SQRPIH*SQRT(RWVAL) - REXPW*GVAL*RWVAL
       DO J=1,JMAX
          RJ000(J,I) = RWVAL*((J - HALF)*RJ000(J-1,I)-REXPW)
       ENDDO
       !  (2J+36) < WVAL 
    ELSE
       RWVAL = PID4/WVAL
       RJ000(0,I) = SQRT(RWVAL)
       RWVAL = RWVAL*PID4I
       DO J = 1, JMAX
          RJ000(J,I) = RWVAL*(J - HALF)*RJ000(J-1,I)
       ENDDO
    END IF
 ENDDO

 ! Scaling
 DO I=1,nPrim
    PREF = Prefactor(I)
    RJ000(0,I) = PREF*RJ000(0,I)
    IF (jmax.GT.0) THEN
       D2MALPHA = -2*alpha(I)
       DO j=1,jmax
          PREF = PREF*D2MALPHA
          RJ000(J,I) = PREF*RJ000(J,I)
       ENDDO
    ENDIF
 ENDDO

END SUBROUTINE buildRJ000_OP_HA

!> \brief tabulation of the boys function or incomplete gamma function \f$ F_{n} \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param LUPRI logical unit number for output printing
!> \param JMX tabulate to this order of angular momentum
!> \param sharedTUV stores the tabulated values
!>
!> ***** Tabulation of incomplete gamma function *****
!>   For J = JMX a power series expansion is used, see for
!>   example Eq.(39) given by V. Saunders in "Computational
!>   Techniques in Quantum Chemistry and Molecular Physics",
!>   Reidel 1975.  For J < JMX the values are calculated
!>   using downward recursion in J.
!>
SUBROUTINE gammaTabulation(LUPRI,JMX,sharedTUV)
 use TYPEDEF
 IMPLICIT NONE
 TYPE(TUVitem),intent(inout) :: sharedTUV
 INTEGER,intent(in)       :: JMX,LUPRI
!
 INTEGER           :: MAXJ0,IADR,IPOINT,IORDER,JADR,JMAX,J
 INTEGER,PARAMETER :: MXQN=13
 INTEGER,PARAMETER :: MAXJ = 4*(MXQN - 1) + 2
 REAL(REALK)       :: DENOM,D2MAX1,R2MAX1,TERM,SUM,REXPW,WVAL,D2WAL
 REAL(REALK), PARAMETER :: HALF = 0.5D0,  TEN6 = 1.0D6
 REAL(REALK), PARAMETER :: D1 = 1.D0, D10 = 10.D0
 REAL(REALK), PARAMETER :: D2 = 2.D0, D4 = 4.D0, D12 = 12.D0, TENTH = 0.01D0
 REAL(REALK), PARAMETER :: COEF2 = HALF,  COEF3 = - D1/6.D0, COEF4 = D1/24.D0
 REAL(REALK), PARAMETER :: COEF5 = - D1/120.D0, COEF6 = D1/720.D0

 REAL(REALK), PARAMETER :: PI    = 3.14159265358979323846D00
 REAL(REALK), PARAMETER :: SQRTPI = 1.77245385090551602730D00
 REAL(REALK), PARAMETER :: R2PI52 = 5.91496717279561287782D00
 REAL(REALK), PARAMETER :: SQRPIH = SQRTPI/D2
 REAL(REALK), PARAMETER :: PID4 = PI/D4, PID4I = D4/PI

 REAL(REALK), PARAMETER :: GFAC30 =  .4999489092D0 
 REAL(REALK), PARAMETER :: GFAC31 = -.2473631686D0
 REAL(REALK), PARAMETER :: GFAC32 =  .321180909D0
 REAL(REALK), PARAMETER :: GFAC33 = -.3811559346D0
 REAL(REALK), PARAMETER :: GFAC20 = .4998436875D0
 REAL(REALK), PARAMETER :: GFAC21 = -.24249438D0
 REAL(REALK), PARAMETER :: GFAC22 =  .24642845D0
 REAL(REALK), PARAMETER :: GFAC10 =  .499093162D0
 REAL(REALK), PARAMETER :: GFAC11 = -.2152832D0
 REAL(REALK), PARAMETER :: GFAC00 =  .490D0

 IF (JMX .GT. MAXJ) THEN
    WRITE (LUPRI,'(//A,I5,A,I3)')&
         &      ' GAMTAB ERROR: JMX =',JMX,', which is greater than',MAXJ
    CALL LSQUIT('GAMTAB ERROR: JMX greater than limit.',lupri)
 END IF
 JMAX = JMX + 6
 MAXJ0 = JMAX
 !
 !     WVAL = 0.0
 !
 IADR = 1
 DENOM = D1
 DO J = 0,JMAX
    SHAREDTUV%TABFJW(J,0) = D1/DENOM
    IADR = IADR + 1201
    DENOM = DENOM + D2
 ENDDO
 !
 !     WVAL = 0.1, 0.2, 0.3,... 12.0
 !
 IADR = IADR - 1201
 !D2MAX1 = DFLOAT(2*JMAX + 1)
 D2MAX1 = 2*JMAX + 1
 R2MAX1 = D1/D2MAX1
 DO IPOINT = 1,1200
    !  WVAL = TENTH*DFLOAT(IPOINT)
    WVAL = TENTH*IPOINT
    D2WAL = WVAL + WVAL
    IADR = IADR + 1
    TERM = R2MAX1
    SUM = TERM
    DENOM = D2MAX1
    DO IORDER = 2,200
       DENOM = DENOM + D2
       TERM = TERM*D2WAL/DENOM
       SUM = SUM + TERM
       IF (TERM .LE. 1.0D-15) EXIT
    ENDDO
    REXPW = EXP(-WVAL)
    SHAREDTUV%TABFJW(JMAX,IPOINT) = REXPW*SUM
    DENOM = D2MAX1
    JADR = IADR
    DO J = 1,JMAX
       DENOM = DENOM - D2
       SHAREDTUV%TABFJW(JMAX-J,IPOINT) = (SHAREDTUV%TABFJW(JMAX-J+1,IPOINT)*D2WAL + REXPW)/DENOM
       JADR = JADR - 1201
    ENDDO
 ENDDO
END SUBROUTINE gammaTabulation

!> \brief initialization routine to init the TUVitem and 
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param sharedTUV contain tabulated boys function and tuvindexing
!> \param Integral contain temporary arrays and sharedTUV
!> \param Input contain info about the requested integral 
!> \param OD_LHS contain the left hand side overlapdistribution 
!> \param OD_RHS contain the right hand side overlapdistribution 
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE initTUVitem(sharedTUV,Integral,Input,OD_LHS,OD_RHS,LUPRI,IPRINT)
 implicit none
 TYPE(TUVitem),intent(inout),target :: sharedTUV
 TYPE(Integralitem),intent(inout)   :: Integral
 TYPE(IntegralInput),intent(in)     :: Input
 TYPE(ODITEM),intent(in)            :: OD_LHS,OD_RHS
 Integer,intent(in)                 :: LUPRI,IPRINT
 !
 LOGICAL              :: TABULATE_BOYS
 Integer,parameter :: MAXAOJ=12,MAXDER=2
 Integer,parameter :: MAXJ=4*MAXAOJ+MAXDER
 !Integer,parameter :: MAXAOJ=12,MAXDER=2
 !Integer,parameter :: MAXJ=12
 Integer,parameter :: NODES=1201,ORDER=7
 Integer,parameter :: NTABFJ000=NODES*(MAXJ+ORDER)
 Integer,parameter :: NTUVMAX = (MAXJ+1)*(MAXJ+2)*(MAXJ+3)/6
 !
 TABULATE_BOYS = .FALSE.
 IF(INPUT%Operator(1:7) .EQ. 'Coulomb') TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator(1:6).EQ.'Nucrep') TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator(1:4) .EQ. 'Erfc') TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator(1:4) .EQ. 'Erf ') TABULATE_BOYS = .TRUE.

 IF(TABULATE_BOYS)THEN
    !Pretabulate Gammafunction/Boysfunction 
    NULLIFY(sharedTUV%TABFJW)  !MXQN=13     MAXJ = 4*(MXQN - 1) + 2 
    ALLOCATE(SharedTUV%TABFJW(0:MAXJ+6,0:1200))     !121*(MAXJ + 7) = 6897
    SharedTUV%nTABFJW1=MAXJ+6
    SharedTUV%nTABFJW2=1200
    CALL GAMMATABULATION(lupri,MAXJ,sharedTUV)    !Jmax is set to 10 is that ok?
 ENDIF
 NULLIFY(SharedTUV%TUVindex)
 NULLIFY(SHAREDTUV%Tindex)
 NULLIFY(SHAREDTUV%Uindex)
 NULLIFY(SHAREDTUV%Vindex)
 ALLOCATE(SharedTUV%TUVindex(0:MAXJ,0:MAXJ,0:MAXJ))
 ALLOCATE(SHAREDTUV%Tindex(NTUVMAX))
 ALLOCATE(SHAREDTUV%Uindex(NTUVMAX))
 ALLOCATE(SHAREDTUV%Vindex(NTUVMAX))

 CALL integralTUVindex(SharedTUV,MAXJ,LUPRI,IPRINT)

 CALL BuildPrecalculatedSphmat(sharedTUV,OD_LHS,OD_RHS,LUPRI,IPRINT)

 integral%TUV => sharedTUV

END SUBROUTINE initTUVitem

!> \brief build spherical transformation matrices and attach to sharedTUV 
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param sharedTUV to store the matrices
!> \param OD_LHS contain the left hand side overlapdistribution 
!> \param OD_RHS contain the right hand side overlapdistribution 
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE buildPrecalculatedSphmat(sharedTUV,OD_LHS,OD_RHS,LUPRI,IPRINT)
 IMPLICIT NONE
 TYPE(ODITEM),intent(in)     :: OD_LHS,OD_RHS
 TYPE(TUVitem),intent(inout) :: SharedTUV
 INTEGER,intent(in)          :: LUPRI,IPRINT
 !
 INTEGER  :: I,MAXANGMOM

 MAXANGMOM = 0
 DO I=1,OD_LHS%nbatches
    MAXANGMOM = MAX(MAXANGMOM,OD_LHS%BATCH(I)%AO(1)%p%maxangmom)
    MAXANGMOM = MAX(MAXANGMOM,OD_LHS%BATCH(I)%AO(2)%p%maxangmom)
 ENDDO
 DO I=1,OD_RHS%nbatches
    MAXANGMOM = MAX(MAXANGMOM,OD_RHS%BATCH(I)%AO(1)%p%maxangmom)
    MAXANGMOM = MAX(MAXANGMOM,OD_RHS%BATCH(I)%AO(2)%p%maxangmom)
 ENDDO

 SharedTUV%nSPHMAT = MAXANGMOM + 1
 CALL setPrecalculatedSphmat(SharedTUV%SPH_MAT,SharedTUV%nSPHMAT,LUPRI,IPRINT)

END SUBROUTINE buildPrecalculatedSphmat

!> \brief set the TUVindex,Tindex,Uindex and Vindex  
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param sharedTUV to store the info
!> \param MAXJ the maximum angular momentum for this calculation
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE integralTUVindex(sharedTUV,MAXJ,LUPRI,IPRINT)
 implicit none
 TYPE(TUVitem),intent(inout)  :: sharedTUV
 Integer,intent(in)        :: LUPRI,IPRINT,MAXJ
!
 Integer :: TUV,T,U,V,J
!
 TUV=0
 DO J = 0, MAXJ
    DO T = J,0,-1       
       DO U = J-T,0,-1
          TUV=TUV+1
          V=J-T-U
          sharedTUV%TUVindex(T,U,V)=TUV
          sharedTUV%Tindex(TUV)=T
          sharedTUV%Uindex(TUV)=U
          sharedTUV%Vindex(TUV)=V
       ENDDO
    ENDDO
 ENDDO
END SUBROUTINE integralTUVindex

!> \brief deallocate memory for the TUVitem
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param sharedTUV to store the info
!> \param Input contain info about the requested integral 
SUBROUTINE freeTUVitem(sharedTUV,Input)
 implicit none
 TYPE(TUVItem),intent(inout)  :: sharedTUV
 TYPE(IntegralInput),intent(in) :: Input
!
 LOGICAL              :: TABULATE_BOYS
!
 TABULATE_BOYS = .FALSE.
 IF(INPUT%Operator(1:7) .EQ. 'Coulomb') TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator(1:6).EQ.'Nucrep') TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator(1:4) .EQ. 'Erfc') TABULATE_BOYS = .TRUE.
 IF(INPUT%Operator(1:4) .EQ. 'Erf ') TABULATE_BOYS = .TRUE.

 IF(TABULATE_BOYS) DEALLOCATE(sharedTUV%TABFJW)
 DEALLOCATE(SharedTUV%TUVindex)
 DEALLOCATE(SHAREDTUV%Tindex)
 DEALLOCATE(SHAREDTUV%Uindex)
 DEALLOCATE(SHAREDTUV%Vindex)
 CALL freePrecalculatedSphmat(SharedTUV%SPH_MAT,SharedTUV%nSPHMAT)
END SUBROUTINE freeTUVitem

!> \brief multipole moment hermite integrals 
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param INTEGRAL store temporary arrays
!> \param Input contain info about the requested integral 
!> \param PQ Integrand containing info about overlapdistribution P and Q
!> \param nPrim number of primitive functions
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param CARTORDER cartesian order of the multipole moments
!> \param INPUTORIGO origin of the multipole moments defined in input 
!>
!> The final MOM(T,U,V) integrals are arranged as follows:
!> S(000)
!> S(100) S(010) S(001)
!> S(200) S(110) S(101) S(020) S(011) S(002)
!> S(300) S(210) S(201) S(120) S(111) S(102) S(030) S(021) S(012) S(003)
!>
SUBROUTINE getMomtuv(INTEGRAL,INPUT,PQ,nprim,LUPRI,IPRINT,CARTORDER,INPUTORIGO)
 implicit none
 TYPE(Integrand),intent(in)         :: PQ
 TYPE(Integralitem),intent(inout)   :: integral
 TYPE(IntegralInput),intent(in)     :: Input
 LOGICAL,intent(in)                 :: INPUTORIGO
 INTEGER,intent(in)                 :: LUPRI,nprim,IPRINT,CARTORDER
!
 INTEGER                 :: J,T,U,V
 INTEGER                 :: ntuv,I
 INTEGER                 :: Jmax,Jstart,e,f,g,efg,nefg
 INTEGER                 :: Xorder,Yorder,Zorder,C
 real(realk)             :: DISTANCE(nprim,3)
 Real(realk), parameter  :: OneHalf=1.5d0
 Real(realk), parameter  :: PI=3.14159265358979323846D0

 Xorder=CARTORDER
 Yorder=CARTORDER
 Zorder=CARTORDER

! NPrim=PQ%nPrimitives
 JMAX=PQ%endAngmom!+CARTORDER
 Jstart=PQ%startAngmom !always zero
 ntuv=(JMAX+1)*(JMAX+2)*(JMAX+3)/6
 nEFG=(CARTORDER+1)*(CARTORDER+2)*(CARTORDER+3)/6
 INTEGRAL%nTUV=ntuv
 INTEGRAL%nEFG=nefg

 IF (JMAX+CARTORDER .EQ. 0) THEN
    !Special case JMAX = 0 A SIMPEL OVERLAP MATRIX WITH ONLY S ORBITALS
    !WRITE(LUPRI,*)'SPECIAL CASE OF CARTISIAN MULTIPOLE MOMENTS - OVERLAP'
    DO I=1,NPrim
       Integral%Rtuv(I)=(PI/PQ%P%p%exponents(I))**OneHalf
    ENDDO
 ELSE
    IF(CARTORDER .EQ. 0) THEN
!!$      !SIMPEL OVERLAP
!!$      NULLIFY(INTEGRAL%Wtuv)
!!$      ALLOCATE(INTEGRAL%Wtuv(Nprim,nTUV))
!!$      DO J = Jstart, PQ%endAngmom
!!$         DO T = J,0,-1
!!$            DO U = J-T,0,-1
!!$               V=J-T-U
!!$               DO I=1,nPrim
!!$                  INTEGRAL%Wtuv(I,INTEGRAL%TUVINDEX(T,U,V))=INTTEMP2(I,INTEGRAL%TUVindex(T,U,V))
!!$               ENDDO
!!$            ENDDO
!!$         ENDDO
!!$      ENDDO
    ELSE
       !      IF(INPUT%orderAngPrim)THEN
       !         CALL QUIT('multipole moments should be run with orderAngPrim = .FALSE.')
       !      ELSE
       IF(INPUTORIGO)THEN
          DO I = 1,Nprim
             DISTANCE(I,1) = PQ%P%p%center(1,I)-INPUT%ORIGO(1)
             DISTANCE(I,2) = PQ%P%p%center(2,I)-INPUT%ORIGO(2)
             DISTANCE(I,3) = PQ%P%p%center(3,I)-INPUT%ORIGO(3)
             !        WRITE(LUPRI,'(2X,A,3F16.9)')'DISTANCE =',DISTANCE(I,1),DISTANCE(I,2),DISTANCE(I,3)
          ENDDO
       ELSE
          DO I = 1,Nprim
             DISTANCE(I,1) = PQ%P%p%center(1,I)-PQ%P%p%ODcenter(1)
             DISTANCE(I,2) = PQ%P%p%center(2,I)-PQ%P%p%ODcenter(2)
             DISTANCE(I,3) = PQ%P%p%center(3,I)-PQ%P%p%ODcenter(3)
             !        WRITE(LUPRI,'(2X,A,3F16.9)')'DISTANCE =',DISTANCE(I,1),DISTANCE(I,2),DISTANCE(I,3)
          ENDDO
       ENDIF
       !      CALL LS_DZERO(INTEGRAL%Rtuv,nTUV*nPrim*nEFG)
       CALL MultipoleRecurrence(INTEGRAL%Rtuv,JMAX,CARTORDER,&
            & distance(:,1),distance(:,2),distance(:,3),&
            & Xorder,Yorder,Zorder,NPrim,nTUV,nEFG,&
            & PQ%P%p%Exponents,lupri)
       !   ENDIF
    ENDIF
 ENDIF
 !    Print section
 !    =============
 IF (IPRINT .GE. 10) THEN
    CALL LSHEADER(LUPRI,'Output from ERITUV')
    WRITE (LUPRI,'(2X,A13,I10)') 'MAX angmom ', JMAX
    WRITE (LUPRI,'(2X,A13,I10)') 'NPrim  ', NPrim
    WRITE (LUPRI,'(2X,A13,I10)') 'NTUV  ', nTUV
    WRITE (LUPRI,'(2X,A13,I10)') 'NEFG  ', nEFG
    IF (IPRINT .GE. 20) THEN
       CALL LSHEADER(LUPRI,'Hermite integrals M(t,u,v)')
       efg=0
       DO C=0,CARTORDER
          DO e = C,0,-1
             DO f = C-e,0,-1
                g=C-e-f
                efg=efg+1
                WRITE(LUPRI,*)'MULTIPOLE MOMENT LEVEL = (',e,',',f,',',g,')   EFG=',efg
                DO J = 0, JMAX
                   DO T = J,0,-1
                      DO U = J-T,0,-1
                         V=J-T-U
                         WRITE(LUPRI,*)'MULTIPOLE MOMENT LEVEL = (',T,',',U,',',V,')   TUV=',INTEGRAL%TUV%TUVINDEX(T,U,V)
                         WRITE (LUPRI,'(2X,A7,I1,A1,I1,A1,I1,A1,2X,5F12.8/,(12X,5F12.8))')&
                              & 'CARMOM(',T,',',U,',',V,')', &
                              &(INTEGRAL%Rtuv(INTEGRAL%TUV%TUVINDEX(T,U,V)+(efg-1)*nTUV+(I-1)*nTUV*nEFG),I=1,NPrim)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    END IF
 END IF
 Integral%nPrim=nPrim

END SUBROUTINE getMomtuv

!> \brief multipole moment hermite integral recurrence relations \f$ M^{e+1}_{t} = t M^{e}_{t-1}+X_{PC}M^{e}_{t}+\frac{1}{2p}M^{e}_{t+1} \f$ 
!> \author \latexonly T. Kj{\ae}rgaard \endlatexonly
!> \date 2009-02-05
!> \param OUTPUT of the recurrence relations
!> \param Jstart start angular momentum
!> \param JMAX maximum angular momentum
!> \param CARTORDER cartesian order of the multipole moments
!> \param Xpq distance in X axis
!> \param Ypq distance in Y axis
!> \param Zpq distance in Z axis
!> \param Xorder multipole moment expansion order in X
!> \param Yorder multipole moment expansion order in Y
!> \param Zorder multipole moment expansion order in Z
!> \param nPrim number of primitive functions
!> \param nTUV number of TUV components
!> \param nEFG number of cartesian moment components
!> \param TUVindex index array 
!> \param EXPONENTS
!> \param LUPRI logical unit number for output printing
!>
!> Se section 9.5.3 in the Book (page 356)
!> \f[ M^{efg}_{tuv} = M^{e}_{t}  M^{f}_{u}  M^{g}_{v} \f]
!> \f[ M^{e+1}_{t} = t M^{e}_{t-1}+X_{PC}M^{e}_{t}+\frac{1}{2p}M^{e}_{t+1} \f] 
!> starting with \f$ M^{0}_{000} = \sqrt{\frac{\pi}{p}} \f$
!>
SUBROUTINE multipoleRecurrence(OUTPUT,JMAX,CARTORDER,&
    & Xpq,Ypq,Zpq,Xorder,Yorder,Zorder,NPrim,nTUV,nEFG,&
    & EXPONENTS,lupri)
 IMPLICIT NONE
 INTEGER,PARAMETER        :: MAXJ = 50
 INTEGER,intent(in)       :: Xorder,Yorder,Zorder,nTUV,nEFG,Jmax,nprim
 !INTEGER,intent(in)       :: TUVindex(0:MAXJ,0:MAXJ,0:MAXJ)
INTEGER,intent(in)       :: CARTORDER
 real(realk),intent(inout):: OUTPUT(NTUV*nPrim*nEFG)
 real(realk),intent(in)   :: Xpq(nPrim),Ypq(nPrim),Zpq(nPrim),EXPONENTS(nPrim)
!
 INTEGER                 :: J,T,U,V,TUV
 real(realk)               :: P(nPrim),M000(nPrim)
 INTEGER                 :: I,e,f,g,efg
 INTEGER                 :: lupri,C
 real(realk)             :: MX(nprim,0:Xorder,0:Xorder)
 real(realk)             :: MY(nprim,0:Yorder,0:Yorder)
 real(realk)             :: MZ(nprim,0:Zorder,0:Zorder)
 Real(realk), parameter  :: PI=3.14159265358979323846D0
 !#IFDEF VAR_MKL
 !real(realk)             :: TEMP(nPrim)
 !#ENDIF

 !CALL LSHEADER(LUPRI,'Multipole recurrence')
 DO I = 1, Nprim
    P(I)=1/(2.d0*EXPONENTS(I))   
 ENDDO
 !#IFDEF VAR_MKL
 !  DO I=1,NPrim
 !     TEMP(I)=PI/EXPONENTS(I)
 !  ENDDO
 !call vdsqrt(nprim,TEMP,M000)
 !#ELSE
 DO I=1,NPrim
    M000(I)=sqrt(PI/EXPONENTS(I))
 ENDDO
 !#ENDIF

 CALL LS_DZERO(MX,nPrim*(Xorder+1)*(Xorder+1))
 IF(Xorder .GE. 1)THEN
    IF(Xorder .GT. 1)THEN
       CALL STANDARDLOOPX(M000,MX,Xpq,Xorder,nPrim,P,JMAX,lupri)
    ELSEIF(Xorder .EQ. 1)THEN
       CALL ORDER_EQ_ONE_LOOP(M000,MX,Xpq,Xorder,nPrim,JMAX,lupri)
    ENDIF
 ELSE
    CALL DCOPY(nPrim,M000,1,MX(:,0,0),1)
 ENDIF

 CALL LS_DZERO(MY,nPrim*(Yorder+1)*(Yorder+1))
 IF(Yorder .GE. 1)THEN
    IF(Yorder .GT. 1)THEN
       CALL STANDARDLOOPX(M000,MY,Ypq,Yorder,nPrim,P,JMAX,lupri)
    ELSEIF(Yorder .EQ. 1)THEN
       CALL ORDER_EQ_ONE_LOOP(M000,MY,Ypq,Yorder,nPrim,JMAX,lupri)
    ENDIF
 ELSE
    CALL DCOPY(nPrim,M000,1,MY(:,0,0),1)
 ENDIF

 CALL LS_DZERO(MZ,nPrim*(Zorder+1)*(Zorder+1))
 IF(Zorder .GE. 1)THEN
    IF(Zorder .GT. 1)THEN
       CALL STANDARDLOOPX(M000,MZ,Zpq,Zorder,nPrim,P,JMAX,lupri)
    ELSEIF(Zorder .EQ. 1)THEN
       CALL ORDER_EQ_ONE_LOOP(M000,MZ,Zpq,Zorder,nPrim,JMAX,lupri)
    ENDIF
 ELSE
    CALL DCOPY(nPrim,M000,1,MZ(:,0,0),1)
 ENDIF
 efg=0
 !Nefg=(CARTORDER+1)*(CARTORDER+2)*(CARTORDER+3)/6
 DO C=0,CARTORDER
    DO e = C,0,-1
       DO f = C-e,0,-1
          g=C-e-f
          efg=efg+1
          TUV=0
          DO J = 0, JMAX
             DO T = J,0,-1
                DO U = J-T,0,-1
                   V=J-T-U
                   TUV=TUV+1
                   IF(T .LE. e)THEN
                      IF(U .LE. f)THEN
                         IF(V .LE. g)THEN
                            DO I=1,nPrim
                               OUTPUT(TUV+(efg-1)*nTUV+(I-1)*nTUV*nEFG)=MX(I,e,T)*MY(I,f,U)*MZ(I,g,V)
                            ENDDO
                         ELSE
                            DO I=1,nPrim
                               OUTPUT(TUV+(efg-1)*nTUV+(I-1)*nTUV*nEFG)=0.D0
                            ENDDO
                         ENDIF
                      ELSE
                         DO I=1,nPrim
                            OUTPUT(TUV+(efg-1)*nTUV+(I-1)*nTUV*nEFG)=0.D0
                         ENDDO
                      ENDIF
                   ELSE
                      DO I=1,nPrim
                         OUTPUT(TUV+(efg-1)*nTUV+(I-1)*nTUV*nEFG)=0.D0
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
 ENDDO

END SUBROUTINE multipoleRecurrence

!> \brief single direction multipole moment hermite integral recurrence relations \f$ M^{e+1}_{t} = t M^{e}_{t-1}+X_{PC}M^{e}_{t}+\frac{1}{2p}M^{e}_{t+1} \f$ 
!> \author \latexonly T. Kj{\ae}rgaard \endlatexonly
!> \date 2009-02-05
!> \param INPUT \f$ M^{0}_{0}\f$
!> \param M output \f$ M^{e}_{t}\f$ 
!> \param Xpq distance in given direction
!> \param Xorder expansion order in given direction
!> \param nPrim number of primitive functions
!> \param P exponents
!> \param LIMIT the maximum angular moment for the recurrence
!> \param LUPRI logical unit number for output printing
!> 
!> written for X loop (therefore Xpq and Xorder) but general for any direction
!>
SUBROUTINE StandardLoopX(INPUT,M,Xpq,Xorder,nPrim,P,LIMIT,lupri)
 IMPLICIT NONE
 INTEGER,intent(in)       :: nPrim,LIMIT,Xorder,lupri
 real(realk),intent(in)   :: Xpq(nPrim),P(nPrim),INPUT(nPrim)
 real(realk),intent(inout):: M(nPrim,0:Xorder,0:Xorder)
!
 INTEGER   :: C,T,I

 CALL DCOPY(nPrim,INPUT,1,M(1,0,0),1)

 !M(prim,e,t)
 DO C=1,Xorder-1
    IF(C .EQ. 1 )THEN
       DO I = 1, Nprim
          M(I,1,0)=Xpq(I)*INPUT(I) 
       ENDDO
       DO I = 1, Nprim
          M(I,1,1)=INPUT(I)              
       ENDDO
    ELSEIF(C .EQ. 2)THEN
       DO I = 1, Nprim
          M(I,2,0)=Xpq(I)*M(I,1,0)+P(I)*M(I,1,1)     
       ENDDO
       DO I = 1, Nprim
          M(I,2,1)=M(I,1,0)+Xpq(I)*M(I,1,1)          
       ENDDO
       DO I = 1, Nprim
          M(I,2,2)=2*M(I,1,1)                      
       ENDDO
    ELSE
       DO I = 1, Nprim
          M(I,C,0)= Xpq(I)*M(I,C-1,0)+P(I)*M(I,C-1,1)     
       ENDDO
       DO T = 1,C-2
          DO I = 1, Nprim
             M(I,C,T)=T*M(I,C-1,T-1)+Xpq(I)*M(I,C-1,T)+P(I)*M(I,C-1,T+1)
          ENDDO
       ENDDO
       !        T=C-1
       DO I = 1, Nprim
          M(I,C,C-1)=(C-1)*M(I,C-1,C-2)+Xpq(I)*M(I,C-1,C-1)
       ENDDO
       !        T=C
       DO I = 1, Nprim
          M(I,C,C)=C*M(I,C-1,C-1)
       ENDDO
    ENDIF
 ENDDO
 C=Xorder
 IF(C .EQ. 2)THEN
    DO I = 1, Nprim
       M(I,2,0)=Xpq(I)*M(I,1,0)+P(I)*M(I,1,1)     
    ENDDO
    IF(LIMIT .GT. 0)THEN
       DO I = 1, Nprim
          M(I,2,1)=M(I,1,0)+Xpq(I)*M(I,1,1)          
       ENDDO
       IF(LIMIT .GT. 1)THEN
          DO I = 1, Nprim
             M(I,2,2)=2*M(I,1,1)                      
          ENDDO
       ENDIF
    ENDIF
 ELSE
    DO I = 1, Nprim
       M(I,C,0)= Xpq(I)*M(I,C-1,0)+P(I)*M(I,C-1,1)     
    ENDDO
    DO T = 1,MIN(LIMIT,C-2)
       DO I = 1, Nprim
          M(I,C,T)=T*M(I,C-1,T-1)+Xpq(I)*M(I,C-1,T)+P(I)*M(I,C-1,T+1)
       ENDDO
    ENDDO
    IF(LIMIT .GT. C-2)THEN
       !        T=C-1
       DO I = 1, Nprim
          M(I,C,C-1)=(C-1)*M(I,C-1,C-2)+Xpq(I)*M(I,C-1,C-1)
       ENDDO
       IF(LIMIT .GT. C-1)THEN
          !        T=C
          DO I = 1, Nprim
             M(I,C,C)=C*M(I,C-1,C-1)
          ENDDO
       ENDIF
    ENDIF
 ENDIF

END SUBROUTINE StandardLoopX

!> \brief special multipole moment recurrence relations for expansion order equal 1
!> \author \latexonly T. Kj{\ae}rgaard  \endlatexonly
!> \date 2009-02-05
!> \param INPUT \f$ M^{0}_{0}\f$
!> \param M output \f$ M^{e}_{t}\f$ 
!> \param Xpq distance in given direction
!> \param Xorder expansion order in given direction
!> \param nPrim number of primitive functions
!> \param LUPRI logical unit number for output printing
!> 
!> written for X loop (therefore Xpq and Xorder) but general for any direction
!>
SUBROUTINE order_eq_one_loop(INPUT,M,Xpq,Xorder,nPrim,LIMIT,lupri)
 IMPLICIT NONE
 INTEGER,intent(in)     :: nPrim,Xorder,lupri,LIMIT
 real(realk),intent(inout) :: M(nPrim,0:Xorder,0:Xorder)
 real(realk),intent(in) :: Xpq(nPrim),INPUT(nPrim)
!
 INTEGER                :: I
 CALL DCOPY(nPrim,INPUT,1,M(1,0,0),1)
 DO I = 1, Nprim
    M(I,1,0)=Xpq(I)*INPUT(I)  
 ENDDO
 IF(LIMIT .GT. 0)THEN
    DO I = 1, Nprim
       M(I,1,1)=INPUT(I)                        
    ENDDO
 ENDIF
END SUBROUTINE order_eq_one_loop

!> \brief wrapper which copy integralprefactors and reducedexponents and then call secondary wrapper to 2 electron complementary error function modified Coulomb hermite integral
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param RJ000 \f$ \tilde{R}^{0:jmax}_{000} \f$
!> \param PQ Integrand (info about the overlap distributions P and Q)
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param integral containing tabulated boys function \f$ F_{n} \f$ 
!> \param Omega argument in the complementary error function erfc\f$(\omega r_{12})\f$
!> \param HIGH logical to switch between high accuracy and low accuracy versions
SUBROUTINE buildErfRJ000(RJ000,PQ,nPrim,Jmax,LUPRI,IPRINT,integral,OMEGA,HIGH)
 IMPLICIT NONE
 INTEGER,intent(in)           :: nPrim,LUPRI,IPRINT,Jmax
 REAL(REALK),intent(inout)    :: RJ000(0:Jmax,nPrim)
 REAL(REALK),intent(in)       :: OMEGA
 TYPE(Integrand),intent(in)   :: PQ
 TYPE(integralitem),intent(in):: integral
 LOGICAL,intent(in)           :: HIGH
!
 REAL(REALK) :: Prefactor(nPrim),alpha(nPrim)

 CALL DCOPY(nPrim,PQ%integralPrefactor,1,Prefactor,1)
 CALL DCOPY(nPrim,PQ%reducedExponents,1,alpha,1)
 !COPY BECAUSE BUILD_ERF2_RJ000 changes the values of prefactor and ALPHA
 CALL buildErf2RJ000(OMEGA,RJ000,nPrim,Jmax,LUPRI,IPRINT,&
      &alpha,PQ%squaredDistance,&
      &Prefactor,INTEGRAL%TUV%TABFJW,INTEGRAL%TUV%nTABFJW1,INTEGRAL%TUV%nTABFJW2,HIGH)

END SUBROUTINE buildErfRJ000

!> \brief secondary wrapper to 2 electron complementary error function modified Coulomb hermite integral \f$ \tilde{R}^{n}_{000}(\alpha,\textbf{R}_{PQ}) =  \left(-2 \alpha \frac{\omega^{2}}{\alpha + \omega} \right)^{n}  \sqrt{\frac{\omega^{2}}{\alpha + \omega}} F_{n}(\alpha \frac{\omega^{2}}{\alpha + \omega} R_{PQ}^{2}) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param OMEGA argument in the complementary error function erfc\f$(\omega r_{12})\f$
!> \param RJ000 \f$ \tilde{R}^{0:jmax}_{000} \f$
!> \param nPrim number of primitive ( nPrim for P * nPrim for Q ) 
!> \param Jmax maximum angular momentum ( N where t+u+v =< N )
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
!> \param alpha reduced exponent (to be modified to \f$ \frac{\omega^{2}}{\alpha + \omega} \f$)
!> \param R2 squared distance
!> \param Prefactor \f$ \frac{2 \pi^{\frac{5}{2}}}{pq\sqrt{p+q}} \f$ to be modified to \f$ \sqrt{\frac{\omega^{2}}{\alpha + \omega}} \frac{2 \pi^{\frac{5}{2}}}{pq\sqrt{p+q}} \f$
!> \param TABFJW tabulated boys function \f$ F_{n} \f$ 
!> \param nTABFJW1 dimension 1 of TABFJW
!> \param nTABFJW2 dimension 2 of TABFJW
!> \param HIGH logical to switch between high accuracy and low accuracy versions
SUBROUTINE buildErf2RJ000(OMEGA,RJ000,nPrim,Jmax,LUPRI,IPRINT,alpha,R2,Prefactor,TABFJW,nTABFJW1,nTABFJW2,HIGH)
 IMPLICIT NONE
 INTEGER,intent(in)         :: nPrim,Jmax,Lupri,Iprint,nTABFJW1,nTABFJW2
 REAL(REALK),intent(inout)  :: RJ000(0:Jmax,nPrim),alpha(nPrim),Prefactor(nprim)
 REAL(REALK),intent(in)     :: R2(nprim),TABFJW(0:nTABFJW1,0:nTABFJW2),OMEGA
 LOGICAL,intent(in)         :: HIGH
!
 Real(realk) :: OMEGA2,BETAAT,SQBETA
 INTEGER     :: I

 OMEGA2=OMEGA*OMEGA
 DO I = 1, NPrim
    BETAAT = OMEGA2 / (ALPHA(I) + OMEGA2)  ! beta = omega^2/(alpha+omega^2) 
    SQBETA = SQRT(BETAAT)
    Prefactor(I) = Prefactor(I)*SQBETA
    ALPHA(I) = ALPHA(I)*BETAAT         
 ENDDO

 IF (HIGH) THEN
    !HIGH ACCURACY
    CALL buildRJ000_OP_HA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
         &alpha,R2,Prefactor,TABFJW,nTABFJW1,nTABFJW2)
 ELSE
    !LOW ACCURACY
    CALL buildRJ000_OP_LA(RJ000,nPrim,Jmax,LUPRI,IPRINT,&
         &alpha,R2,Prefactor,TABFJW,nTABFJW1,nTABFJW2)
 ENDIF

END SUBROUTINE buildErf2RJ000

!> \brief Distribute the primitive hermite integrals stored using tuvPQ into OUTPUT(tuvQ,tuvP,nPrimP,nPrimQ) using tuvP and tuvQ
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param tuvTUV output (ntuvQ,ntuvP,nPrimPQ) = (ntuvQ,ntuvP,nPrimP,nPrimQ)
!> \param WTUV input (ntuvPQ,nPrimPQ) 
!> \param TUVindex indexing for the tuv components
!> \param nPrim number of primitive orbitals = nPrimP*nPrimQ
!> \param startP start angular momentum for the P overlap distribution
!> \param endP end angular momentum for the P overlap distribution
!> \param startQ start angular momentum for the Q overlap distribution
!> \param endQ end angular momentum for the Q overlap distribution
!> \param ioffPQ ofset for the tuvPQ index
!> \param nTUVEFGPQ number of tuv components and derivative efg components
!> \param ntuvPQ  number of tuv components for PQ 
!> \param ntuvP number of tuv components for P overlap distribution
!> \param ntuvQ  number of tuv components for Q overlap distribution
!> \param ideriv derivative index
!> \param lupri logical unit number for output printing
SUBROUTINE distributeHermiteQ1(tuvTUV,WTUV,TUVindex,nPrim,startP,endP,startQ,endQ,&
    &                         ioffPQ,nTUVEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv,lupri)
 implicit none
 Integer,intent(in) :: nPrim,startP,endP,startQ,endQ,ioffPQ!,ioffP,ioffQ
 Integer,intent(in) :: ntuvEFGPQ,ntuvPQ,ntuvP,ntuvQ,ideriv
 Integer,pointer    :: TUVindex(:,:,:)
 Real(realk),intent(inout) :: tuvTUV(ntuvQ,ntuvP,nPrim)
 Real(realk),intent(in)    :: WTUV(ntuvEFGPQ,nPrim)
 !
 Integer     :: TUVPQindex(nTUVP*nTUVQ)
 Integer     :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ,ituvPQ
 Integer     :: iPrimPQ,iOFF,lupri,iituvP,iituvQ,inTUVPQ

 iPrimPQ=1
 iOFF=(ideriv-1)*ntuvPQ-ioffPQ
 ituvP = 0
 inTUVPQ = 0
 DO jP = startP,endP
    DO tP=jP,0,-1
       DO uP=jP-tP,0,-1
          vP=jP-tP-uP
          ituvP = ituvP+1
          ituvQ = 0
          DO jQ = startQ,endQ
             DO tQ=jQ,0,-1
                DO uQ=jQ-tQ,0,-1
                   vQ=jQ-tQ-uQ
                   ituvQ=ituvQ+1
                   ituvPQ=TUVindex(tP+tQ,uP+uQ,vP+vQ)+iOFF
                   inTUVPQ = inTUVPQ+1
                   TUVPQindex(inTUVPQ)=ituvPQ 
                   tuvTUV(ituvQ,ituvP,iPrimPQ) = WTUV(ituvPQ,iPrimPQ)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
 ENDDO
 DO iPrimPQ=2,nPrim
    inTUVPQ = 0
    DO iituvP = 1,ituvP
       DO iituvQ = 1,ituvQ
          inTUVPQ = inTUVPQ+1
          ituvPQ = TUVPQindex(inTUVPQ)
          tuvTUV(iituvQ,iituvP,iPrimPQ) = WTUV(ituvPQ,iPrimPQ)
       ENDDO
    ENDDO
 ENDDO

END SUBROUTINE distributeHermiteQ1

!> \brief Print the Distribute primitive hermite integrals stored using (tuvQ,tuvP,nPrimPQ)
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param tuvTUV to be printed
!> \param nPrim number of primitive orbitals = nPrimP*nPrimQ
!> \param startP start angular momentum for the P overlap distribution
!> \param endP end angular momentum for the P overlap distribution
!> \param startQ start angular momentum for the Q overlap distribution
!> \param endQ end angular momentum for the Q overlap distribution
!> \param ntuvP number of tuv components for P overlap distribution
!> \param ntuvQ  number of tuv components for Q overlap distribution
!> \param lupri logical unit number for output printing
SUBROUTINE printHermitePQ(tuvTUV,nPrim,startP,endP,startQ,endQ,&
    &                   ntuvP,ntuvQ,lupri)
 implicit none
 Integer,intent(in)     :: nPrim,startP,endP,startQ,endQ
 Real(realk),intent(in) :: tuvTUV(ntuvQ,ntuvP,nPrim)
 Integer,intent(in)     :: ntuvP,ntuvQ,lupri
 !
 Integer     :: jP,jQ,tP,uP,vP,tQ,uQ,vQ,ituvP,ituvQ
 Integer     :: iPrimPQ

 ituvP = 0
 WRITE(LUPRI,'(3X,A)') '***************************************************************'
 WRITE(LUPRI,'(3X,A)') '***                         HermitePQ'
 WRITE(LUPRI,'(3X,A)') '***************************************************************'
 DO jP = startP,endP
    DO tP=jP,0,-1
       DO uP=jP-tP,0,-1
          vP=jP-tP-uP
          ituvP = ituvP+1
          ituvQ = 0
          DO jQ = startQ,endQ
             DO tQ=jQ,0,-1
                DO uQ=jQ-tQ,0,-1
                   vQ=jQ-tQ-uQ
                   ituvQ=ituvQ+1
                   WRITE(LUPRI,'(5X,A,I1,A,I1,A,I1,A,I1,A,I1,A,I1,A)') &
                        &            'W(',tP,',',uP,',',vP,'|',tQ,',',uQ,',',vQ,') ='

                   WRITE(LUPRI,'(5X,6F10.4/,(5X,6F10.4))') &
                        &           (tuvTUV(ituvQ,ituvP,iPrimPQ),iPrimPQ=1,nPrim)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
 ENDDO

END SUBROUTINE printHermitePQ

!> \brief Contract Ecoefficients with Integrals \f$ OUT(ijk,nOrb,nPrim) = \sum_{ntuv} Ecoeffs(ijk,ntuv,nPrim)*IN(ntuv,nOrb,nPrim) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param IntegralIN input integrals to be contracted with Ecoeffs Ordering(ntuv,nOrb,nPrim)
!> \param IntegralOUT output of contraction with Ecoeffs Ordering(ntuv,nOrb,nPrim)
!> \param Ecoeffs cartesian or hermite Ecoefficients
!> \param startE start index in the Ecoefficient 
!> \param nOrb either nPrimP*ntuvP (since the P overlap has not been worked on yet) or nContQ*nCompOrbQ (since the Q overlap have been contracted)
!> \param nAng number of tuv components
!> \param ijk number of cartesian components or number of spherical components
!> \param nPrim number of primitive functions for the overlap being worked on
SUBROUTINE contractEcoeff1(IntegralIN,IntegralOUT,Ecoeffs,startE,nOrb,nAng,ijk,nPrim)
 implicit none           
 Integer,intent(in)     :: nOrb,nAng,ijk,nPrim,startE
 Real(realk),intent(in) :: IntegralIN(nAng,nOrb,nPrim)
 Real(realk),intent(inout) :: IntegralOUT(ijk,nOrb,nPrim)
 Real(realk),pointer    :: Ecoeffs(:)
 !
 Integer      :: iPrimP,iE
 Real(realk)  :: D1=1.0d0,D0=0.0d0

 IF (ijk.EQ.1) THEN
    CALL contractEcoeff1ss(IntegralIN,IntegralOUT,Ecoeffs,nOrb,nPrim)
 ELSE
    iE = startE
    DO iPrimP = 1,nPrim
       CALL DGEMM('N','N',ijk,nOrb,nAng,D1,Ecoeffs(iE),ijk,&
            &          IntegralIN(1,1,iPrimP),nAng,D0,IntegralOUT(1,1,iPrimP),ijk)
       iE = iE+nAng*ijk
    ENDDO
 ENDIF

END SUBROUTINE contractEcoeff1

!> \brief special case of Contract Ecoefficients with Integrals when ijk=1 \f$ OUT(1,nOrb,nPrim) = \sum_{ntuv} Ecoeffs(1,1,nPrim)*IN(1,nOrb,nPrim) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param IntegralIN input integrals to be contracted with Ecoeffs Ordering(nOrb,nPrim)
!> \param IntegralOUT output of contraction with Ecoeffs Ordering(nOrb,nPrim)
!> \param Ecoeffs cartesian or hermite Ecoefficients
!> \param nOrb either nPrimP*ntuvP (since the P overlap has not been worked on yet) or nContQ*nCompOrbQ (since the Q overlap have been contracted)
!> \param nPrim number of primitive functions for the overlap being worked on
SUBROUTINE contractEcoeff1ss(IntegralIN,IntegralOUT,Ecoeffs,nOrb,nPrim)
 implicit none
 Integer,intent(in)     :: nOrb,nPrim
 Real(realk),intent(in),dimension(nOrb,nPrim) :: IntegralIN
 Real(realk),intent(inout),dimension(nOrb,nPrim) :: IntegralOUT
 Real(realk),intent(in),dimension(nPrim)      :: Ecoeffs
 !
 Integer     :: iPrim,iOrb
 real(realk)  :: E
! Real(realk),parameter :: D1=1.0d0,D0=0.0d0,E

 DO iPrim = 1,nPrim
    E = Ecoeffs(iPrim)
    DO iOrb=1,nOrb
       IntegralOUT(iOrb,iPrim) = E*IntegralIN(iOrb,iPrim)
    ENDDO
 ENDDO

END SUBROUTINE contractEcoeff1ss

!> \brief Contract contraction coeffficients AP ordering \f$ OUT(dim,nCont,nPasses) = \sum_{nPrim} IN(dim,nPrim,nPasses) CC(nPrim,nCont) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param PrimInt to be contracted with contraction coefficients (dim,nPrim,nPasses)
!> \param ContInt output after contraction (dim,nCont,nPasses)
!> \param CC empty input - just allocated memory - to be used for contactionmatrix  
!> \param P overlap which is contracted
!> \param nPrim number of primitive functions
!> \param nCont number of contracted functions
!> \param nPasses number of passes
!> \param nDim dimension of the unaffected dimension
!> \param nC1 number of contracted functions on center A (for P) or C (for Q)
!> \param nC2 number of contracted functions on center B (for P) or D (for Q)
!> \param iA1 angular moment or center A (for P) or C (for Q) (used for family basis set)
!> \param iA2 angular moment or center B (for P) or D (for Q) (used for family basis set)
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE contractBasis_AP(PrimInt,ContInt,CC,P,nPrim,nCont,nPasses,nDim,&
    &                     nC1,nC2,iA1,iA2,LUPRI,IPRINT)
 implicit none
 Integer,intent(in) :: nPrim,nCont,nPasses,nDim,nC1,nC2,iA1,iA2,LUPRI,IPRINT
 Real(realk),intent(in),dimension(nDim,nPrim,nPasses) :: PrimInt
 Real(realk),intent(inout),dimension(nDim,nCont,nPasses) :: ContInt
 Real(realk),intent(in),dimension(nPrim,nC1*nC2)      :: CC
 TYPE(Overlap),intent(in) :: P
 !
 Integer     :: iPass,iPrim,iCont,iDim
 Real(realk),parameter :: D1=1.0d0,D0=0.0d0
 Real(realk) :: TMP

 CALL ConstructContraction_AP(CC,P,nPrim,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
IF(nCont .LE. 5)THEN
   DO iPass=1,nPasses
      iPrim = 1
      DO iCont = 1,nCont
         TMP=CC(iPrim,iCont)
         DO iDim = 1,nDim
            ContInt(iDim,iCont,iPass) = PrimInt(iDim,iPrim,iPass)*TMP
         ENDDO
      ENDDO
   ENDDO
   DO iPass=1,nPasses
      DO iCont = 1,nCont
         DO iPrim =2,nPrim
            TMP=CC(iPrim,iCont)
            DO iDim = 1,nDim
               ContInt(iDim,iCont,iPass) = ContInt(iDim,iCont,iPass) + PrimInt(iDim,iPrim,iPass)*TMP
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ELSE
   DO iPass=1,nPasses
      CALL DGEMM('N','N',nDim,nCont,nPrim,D1,PrimInt(1,1,iPass),nDim,&
           &     CC,nPrim,D0,ContInt(1,1,iPass),nDim)
   ENDDO
ENDIF

END SUBROUTINE contractBasis_AP

!> \brief Contract contraction coeffficients PA ordering \f$ OUT(nCont,nPasses,dim) = \sum_{nPrim} CC(nCont,nPrim)*IN(nPrim,nPasses,dim) \f$
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param PrimInt to be contracted with contraction coefficients (nPrim,nPasses,dim)
!> \param ContInt output after contraction (nCont,nPasses,dim)
!> \param CC empty input - just allocated memory - to be used for contactionmatrix  
!> \param P overlap which is contracted
!> \param nPrim number of primitive functions
!> \param nCont number of contracted functions
!> \param nPasses number of passes
!> \param nDim dimension of the unaffected dimension
!> \param nC1 number of contracted functions on center A (for P) or C (for Q)
!> \param nC2 number of contracted functions on center B (for P) or D (for Q)
!> \param iA1 angular moment or center A (for P) or C (for Q) (used for family basis set)
!> \param iA2 angular moment or center B (for P) or D (for Q) (used for family basis set)
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE contractBasis_PA(PrimInt,ContInt,CC,P,nPrim,nCont,nPasses,nDim,&
    &                     nC1,nC2,iA1,iA2,LUPRI,IPRINT)
 implicit none
 Integer,intent(in) :: nPrim,nCont,nPasses,nDim,nC1,nC2,iA1,iA2,LUPRI,IPRINT
 !Real(realk),dimension(nDim,nPrim,nPasses) :: PrimInt
 !Real(realk),dimension(nDim,nCont,nPasses) :: ContInt
 Real(realk),intent(in),dimension(nPrim,nPasses*nDim) :: PrimInt!TEMP
 Real(realk),intent(inout),dimension(nCont,nPasses*nDim) :: ContInt!TEMP
 !Real(realk),dimension(nPrim,nC1,nC2)      :: CC
 Real(realk),intent(in),dimension(nC1*nC2,nPrim)      :: CC
 TYPE(Overlap),intent(in) :: P
 !
 Integer     :: iPass,iCont,iPrim
 Real(realk),parameter :: D1=1.0d0,D0=0.0d0
 Real(realk) :: TMP

 CALL ConstructContraction_PA(CC,P,nPrim,nC1,nC2,iA1,iA2,LUPRI,IPRINT)
 IF(nCont .LE. 5)THEN
    DO iPass=1,nPasses*nDim
       DO iCont = 1,nCont
          TMP = CC(iCont,1)*PrimInt(1,iPass)
          DO iPrim = 2,nPrim
             TMP = TMP+CC(iCont,iPrim)*PrimInt(iPrim,iPass)
          ENDDO
          ContInt(iCont,iPass) = TMP
       ENDDO
    ENDDO
ELSE
 CALL DGEMM('N','N',nCont,nDim*nPasses,nPrim,D1,CC,nCont,&
      &     PrimInt,nPrim,D0,ContInt,nCont)
 !DO iPass=1,nPasses
 !  CALL DGEMM('N','N',nCont,nDim*nPasses,nPrim,D1,CC,nCont,&
 !       &     PrimInt(:,iPass,:),nPrim,D0,ContInt(:,iPass,:),nCont)
 !ENDDO
ENDIF
END SUBROUTINE contractBasis_PA

!> \brief Print after basisset contraction
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param Contracted multiarrray to be printet
!> \param nAng first dimension
!> \param nOrb second dimension
!> \param nPrim third dimension
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE printContracted(Contracted,nAng,nOrb,nPrim,LUPRI,IPRINT)
 Real(realk),intent(in) :: Contracted(nAng,nOrb,nPrim)
 Integer,intent(in)     :: nAng,nOrb,nPrim
 !
 Integer :: iAng,iOrb,iPrim
 IF (IPRINT.GT.50) THEN
    WRITE(LUPRI,'(5X,A)') '**********************************************************'
    WRITE(LUPRI,'(5X,A)') '***                   Contracted'
    WRITE(LUPRI,'(5X,A)') '**********************************************************'
    DO iAng=1,nAng
       DO iOrb=1,nOrb
          WRITE(LUPRI,'(5X,A,I3,A,I3)') 'iAngP =',iAng,' iOrb =',iOrb
          WRITE(LUPRI,'(17X,5F10.4)') (Contracted(iAng,iOrb,iPrim),iPrim=1,nPrim)
       ENDDO
    ENDDO
 ENDIF
END SUBROUTINE printContracted

!> \brief wrapper to add calculated contracted integrals to temporary array order AP
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param TUVQ store the integrals until the Q overlap is fully contracted
!> \param OldTUVQ completed contracted integrals
!> \param totOrbQ total number of Orbitals on overlap distribution Q
!> \param startOrb start orbital to place in TUVQ
!> \param endOrb end orbital to place in TUVQ
!> \param nTUVP size of 2. dimension number of tuv components
!> \param nPrimP number of primitive orbitals on overlap distribution P
!> \param nCompQ number of orbital components on overlap distribution Q
!> \param nContQ number of contracted orbitals on overlap distribution Q
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE addToTUVQ_AP(TUVQ,OldTUVQ,totOrbQ,startOrb,endOrb,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
 implicit none
 Integer,intent(in)     :: totOrbQ,startOrb,endOrb,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT
 Real(realk),intent(inout) :: TUVQ(nPrimP,nTUVP,totOrbQ)
 Real(realk),intent(in) :: OldTUVQ(nCompQ,nTUVP,nPrimP,nContQ)

 CALL addToTUVQ1_AP(TUVQ(:,:,startOrb:endOrb),OldTUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)

END SUBROUTINE addToTUVQ_AP

!> \brief routine to add calculated contracted integrals to temporary array order AP
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param TUVQ store the integrals until the Q overlap is fully contracted
!> \param OldTUVQ completed contracted integrals
!> \param nTUVP size of 2. dimension number of tuv components
!> \param nPrimP number of primitive orbitals on overlap distribution P
!> \param nCompQ number of orbital components on overlap distribution Q
!> \param nContQ number of contracted orbitals on overlap distribution Q
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE addToTUVQ1_AP(TUVQ,OldTUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
 implicit none
 Integer,intent(in)     :: nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT
 Real(realk),intent(inout) :: TUVQ(nPrimP,nTUVP,nCompQ,nContQ)
 Real(realk),intent(in) :: OldTUVQ(nCompQ,nTUVP,nPrimP,nContQ)
 !
 Integer :: iTUVP,iPrimP,iCompQ,iContQ
 !
 DO iContQ=1,nContQ
    DO iCompQ=1,nCompQ
       DO iTUVP=1,nTUVP
          DO iPrimP=1,nPrimP
             TUVQ(iPrimP,iTUVP,iCompQ,iContQ) = OldTUVQ(iCompQ,iTUVP,iPrimP,iContQ)
          ENDDO
       ENDDO
    ENDDO
 ENDDO

 CALL printTUVQ(TUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
END SUBROUTINE addToTUVQ1_AP

!> \brief wrapper to add calculated contracted integrals to temporary array order PA
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param TUVQ store the integrals until the Q overlap is fully contracted
!> \param OldTUVQ completed contracted integrals
!> \param totOrbQ total number of Orbitals on overlap distribution Q
!> \param startOrb start orbital to place in TUVQ
!> \param endOrb end orbital to place in TUVQ
!> \param nTUVP size of 2. dimension number of tuv components
!> \param nPrimP number of primitive orbitals on overlap distribution P
!> \param nCompQ number of orbital components on overlap distribution Q
!> \param nContQ number of contracted orbitals on overlap distribution Q
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE addToTUVQ_PA(TUVQ,OldTUVQ,totOrbQ,startOrb,endOrb,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
 implicit none
 Integer,intent(in)     :: totOrbQ,startOrb,endOrb,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT
 Real(realk),intent(inout) :: TUVQ(nPrimP,nTUVP,totOrbQ)
 Real(realk),intent(in) :: OldTUVQ(nContQ,nPrimP,nTUVP,nCompQ)

 CALL addToTUVQ1_PA(TUVQ(:,:,startOrb:endOrb),OldTUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)

END SUBROUTINE addToTUVQ_PA

!> \brief routine to add calculated contracted integrals to temporary array order PA
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param TUVQ store the integrals until the Q overlap is fully contracted
!> \param OldTUVQ completed contracted integrals
!> \param nTUVP size of 2. dimension number of tuv components
!> \param nPrimP number of primitive orbitals on overlap distribution P
!> \param nCompQ number of orbital components on overlap distribution Q
!> \param nContQ number of contracted orbitals on overlap distribution Q
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE addToTUVQ1_PA(TUVQ,OldTUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
 implicit none
 Real(realk),intent(inout) :: TUVQ(nPrimP,nTUVP,nCompQ,nContQ)
 Real(realk),intent(in) :: OldTUVQ(nContQ,nPrimP,nTUVP,nCompQ)
 Integer,intent(in)     :: nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT
 !
 Integer :: iTUVP,iPrimP,iCompQ,iContQ
 !
 DO iContQ=1,nContQ
    DO iCompQ=1,nCompQ
       DO iTUVP=1,nTUVP
          DO iPrimP=1,nPrimP
             TUVQ(iPrimP,iTUVP,iCompQ,iContQ) = OldTUVQ(iContQ,iPrimP,iTUVP,iCompQ)
          ENDDO
       ENDDO
    ENDDO
 ENDDO

END SUBROUTINE addToTUVQ1_PA

!> \brief print routine for TUVQ the temporary array to store integrals fully Q contracted
!> \author \latexonly T. Kj{\ae}rgaard, S. Reine  \endlatexonly
!> \date 2009-02-05
!> \param TUVQ the multi dimension array to be printet 
!> \param nTUVP size of 2. dimension
!> \param nPrimP size of 1. dimension
!> \param nCompQ size of 3. dimension
!> \param nContQ size of 4. dimension
!> \param LUPRI logical unit number for output printing
!> \param IPRINT integer containing printlevel
SUBROUTINE printTUVQ(TUVQ,nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT)
 implicit none
 Real(realk),intent(in) :: TUVQ(nPrimP,nTUVP,nCompQ,nContQ)
 Integer,intent(in)     :: nTUVP,nPrimP,nCompQ,nContQ,LUPRI,IPRINT
 !
 Integer :: iTUVP,iPrimP,iCompQ,iContQ
 !
 IF (IPRINT.GT.10) THEN
    CALL LSHEADER(LUPRI,'TUVQ')
    DO iContQ=1,nContQ
       DO iCompQ=1,nCompQ
          DO iTUVP=1,nTUVP
             WRITE(LUPRI,'(5X,A,I3,A,I3,A,I3)') 'iTUVP =',iTUVP,' iCompQ =',iCompQ,' iContQ =',iContQ
             WRITE(LUPRI,'(5X,5F10.4)')  (TUVQ(iPrimP,iTUVP,iCompQ,iContQ), iPrimP=1,nPrimP)
          ENDDO
       ENDDO
    ENDDO
 ENDIF
END SUBROUTINE printTUVQ

END MODULE Thermite_integrals
