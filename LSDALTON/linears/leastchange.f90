module leastchange_module
use decompMod
use TYPEDEF
use matrix_operations

contains

subroutine leastchange_to_oao_basis(decomp,A)
implicit none 
type(decompItem),intent(inout) :: decomp
type(Matrix) :: A, wrk

   call MAT_INIT(wrk,A%nrow,A%ncol)
   call mat_mul(decomp%U_inv,A,'n','n',1.d0,0.d0,wrk)
   call mat_mul(wrk,decomp%U_inv,'n','t',1.d0,0.d0,A)
   call mat_free(wrk)

end subroutine leastchange_to_oao_basis


function leastchange_propint(decomp,ls,lupri,luerr)
implicit none
type(decompItem),intent(inout) :: decomp
integer, parameter  :: nderiv=2, nMAT=10
TYPE(lsitem) :: ls
Real(realk),pointer :: leastchange_propint(:,:,:)
type(Matrix) :: propint(nMAT)
INTEGER               :: i, n,lupri,luerr
 
     n   = ls%setting%BASIS(1)%p%REGULAR%nbast

     DO i=1,nMAT
        call mat_init(propint(I),n,n)
     ENDDO

     call II_get_carmom(lupri,luerr,ls%setting,propint,nMAT,nderiv)

     call mat_free(propint(1))
     call mat_free(propint(6))
     call mat_free(propint(7))
     call mat_free(propint(9))

     call leastchange_to_oao_basis(decomp,propint(2))
     call leastchange_to_oao_basis(decomp,propint(3))
     call leastchange_to_oao_basis(decomp,propint(4))
     call leastchange_to_oao_basis(decomp,propint(5))
     call leastchange_to_oao_basis(decomp,propint(8))
     call leastchange_to_oao_basis(decomp,propint(10))

     allocate(leastchange_propint(n,n,6))

     do i=1,4
      call mat_to_full(propint(i+1),1d0,leastchange_propint(:,:,i))
      call mat_free(propint(i+1))
     enddo
     
     call mat_to_full(propint(8),1d0,leastchange_propint(:,:,5))
     call mat_free(propint(8))
     
     call mat_to_full(propint(10),1d0,leastchange_propint(:,:,6))
     call mat_free(propint(10))

     
  return
end function leastchange_propint


subroutine leastchange_blocktransf(T,CMO,n,nocc)
implicit none
integer     :: n
real(realk) :: T(n,n), CMO(n,n), optwrk, dummy
integer     :: nocc, nvirt,LWORK,INFO
real(realk),allocatable :: U(:,:), S(:), VT(:,:),WORK(:), F(:,:),CB(:,:)
character*27 :: msg

 nvirt = n - nocc

!optimal memory
 call dgesvd('A','A',nocc,nocc,dummy,nocc,dummy,dummy,nocc,dummy,nocc,optwrk,-1,INFO)
 LWORK=nint(optwrk)

 allocate(S(nocc),U(nocc,nocc),VT(nocc,nocc),CB(nocc,nocc),WORK(LWORK))

 CB = CMO(1:nocc,1:nocc)

 call dgesvd('A','A',nocc,nocc,CB,nocc,S,U,nocc,VT,nocc,WORK,LWORK,INFO)

 deallocate(CB,S,WORK)

 if (INFO.ne.0) then
   write(msg,'(A,I3)') 'dgesvd failed with INFO=',INFO
   call lsquit(msg,-1)
 endif

 
 allocate(F(nocc,nocc))

 call dgemm('T','T',nocc,nocc,nocc,1d0,VT,nocc,U,nocc,0d0,F,nocc)

 deallocate(U,VT)

 call dgemm('N','N',nocc,nocc,nocc,1d0,CMO(1,1),n,F,nocc,0d0,T(1,1),n)
 call dgemm('N','N',nvirt,nocc,nocc,1d0,CMO(nocc+1,1),n,F,nocc,0d0,T(nocc+1,1),n)

 deallocate(F)

!optimal memory
 call dgesvd('A','A',nvirt,nvirt,dummy,nvirt,dummy,dummy,nvirt,dummy,nvirt,optwrk,-1,INFO)
 LWORK=nint(optwrk)
 
 allocate(S(nvirt),U(nvirt,nvirt),VT(nvirt,nvirt),CB(nvirt,nvirt),WORK(LWORK))

 CB = CMO(nocc+1:n,nocc+1:n)

 call dgesvd('A','A',nvirt,nvirt,CB,nvirt,S,U,nvirt,VT,nvirt,WORK,LWORK,INFO)

 deallocate(CB,S,WORK)

 if (INFO.ne.0) then
   write(msg,'(A,I3)') 'dgesvd failed with INFO=',INFO
   call lsquit(msg,-1)
 endif

 
 allocate(F(nvirt,nvirt))

 call dgemm('T','T',nvirt,nvirt,nvirt,1d0,VT,nvirt,U,nvirt,0d0,F,nvirt)
 
 deallocate(U,VT)

 call dgemm('N','N',nocc,nvirt,nvirt,1d0,CMO(1,nocc+1),n,F,nvirt,0d0,T(1,nocc+1),n)
 call dgemm('N','N',nvirt,nvirt,nvirt,1d0,CMO(nocc+1,nocc+1),n,F,nvirt,0d0,T(nocc+1,nocc+1),n)

 deallocate(F)

end subroutine leastchange_blocktransf

subroutine leastchange_orbspread(spread,uindex,n,&
                     &T,DIPX,DIPY,DIPZ,SECX,SECY,SECZ)
implicit none
integer :: n, uindex(2*n)
real(realk) :: spread(n),T(n,n),DIPX(n,n),DIPY(n,n),DIPZ(n,n)
real(realk) :: SECX(n,n),SECY(n,n),SECZ(n,n)
real(realk),allocatable :: tmp(:), tmpT(:,:)
real(realk) :: sec,dip
integer     :: i, itmp
real(realk), external :: ddot


 call leastchange_iuindex(uindex(n+1),uindex,n)

 allocate(tmp(n),tmpT(n,n))

 call leastchange_rowreorder(tmpT,T,uindex(n+1),n)

 T = tmpT

 deallocate(tmpT)


 

 do i=1,n
     call dsymv('u',n,1d0,SECX,n,T(1,i),1,0d0,tmp,1)
     sec = ddot(n,T(1,i),1,tmp,1)
     call dsymv('u',n,1d0,DIPX,n,T(1,i),1,0d0,tmp,1)
     dip = ddot(n,T(1,i),1,tmp,1)
     spread(i) = sec - dip*dip


     call dsymv('u',n,1d0,SECY,n,T(1,i),1,0d0,tmp,1)
     sec = ddot(n,T(1,i),1,tmp,1)
     call dsymv('u',n,1d0,DIPY,n,T(1,i),1,0d0,tmp,1)
     dip = ddot(n,T(1,i),1,tmp,1)
     spread(i) = spread(i) + sec - dip*dip
     

     call dsymv('u',n,1d0,SECZ,n,T(1,i),1,0d0,tmp,1)
     sec = ddot(n,T(1,i),1,tmp,1)
     call dsymv('u',n,1d0,DIPZ,n,T(1,i),1,0d0,tmp,1)
     dip = ddot(n,T(1,i),1,tmp,1)
     spread(i) = spread(i) + sec - dip*dip

     spread(i) = sqrt(spread(i))
 enddo

 deallocate(tmp)
end subroutine leastchange_orbspread

function leastchange_idmin(n,vec)
integer :: leastchange_idmin, n, i
real(realk) :: vec(n)
  
   leastchange_idmin=1
   do i=2,n
       if ((vec(leastchange_idmin)-vec(i))>1d-6) leastchange_idmin=i
   enddo

return
end function leastchange_idmin


subroutine leastchange_rowreorder(B,A,uindex,n)
implicit none
integer :: n, i
integer :: uindex(n)
real(realk) :: A(n,n), B(n,n)

   do i=1,n
      call dcopy(n,A(uindex(i),1),n,B(i,1),n)
   enddo

end subroutine leastchange_rowreorder

subroutine leastchange_iuindex(iuindex,uindex,n)
integer :: n, i
integer :: iuindex(n),uindex(n)
  
   do i=1,n
      iuindex(uindex(i))=i
   enddo

end subroutine leastchange_iuindex   


subroutine leastchange_rowsort_diag(uindex,CMO,nocc,n)
implicit none
integer :: n,nocc, uindex(n)
real(realk) :: CMO(n,n)
integer, external :: idamax
integer :: i,j,tmp

 do i=1,nocc
       j=idamax(n-i+1,CMO(i,i),1)
       j=j+i-1
       if (j.gt.i) then
        call dswap(n, CMO(j,1),n,CMO(i,1),n)   
        tmp=uindex(i); uindex(i)=uindex(j); uindex(j)=tmp
       endif
 enddo
 
end subroutine leastchange_rowsort_diag

subroutine leastchange_rowsort_occup(uindex,CMO,nocc,n)
implicit none
integer :: n,nocc, uindex(n)
real(realk) :: CMO(n,n)
real(realk),allocatable :: occnum(:)
integer, external     :: idamax
real(realk), external :: ddot
integer :: i,j,tmp
real(realk) :: rtmp


 allocate(occnum(n))
 do i=1,n
   occnum(i) = ddot(nocc,CMO(i,1),n,CMO(i,1),n)
 enddo
 
 do i=1,nocc
       j=idamax(n-i+1,occnum(i),1)
       j=j+i-1
       if (j.gt.i) then
        call dswap(n, CMO(j,1),n,CMO(i,1),n)   
        tmp=uindex(i); uindex(i)=uindex(j); uindex(j)=tmp
        rtmp=occnum(i); occnum(i)=occnum(j); occnum(j)=rtmp
       endif
 enddo

 deallocate(occnum)
 
end subroutine leastchange_rowsort_occup

end module leastchange_module

subroutine leastchangeOrbspreadStandalone(mx,ls,CMO,lupri,luerr)
use TYPEDEF
use matrix_operations
implicit none
real(realk)              :: mx
TYPE(lsitem)             :: ls
type(Matrix)             :: CMO
integer                  :: lupri, luerr

integer                  :: n, i
integer, parameter       :: nderiv=2, nMAT=10
type(Matrix)             :: propint(nMAT), PROPT
real(realk)              :: sec,dip
real(realk), allocatable :: T(:,:), PROPTf(:,:)
real(realk), allocatable :: spread(:)
real(realk), external    :: ddot
integer, external        :: idamax


  n=CMO%nrow

  DO i=1,nMAT
     call mat_init(propint(I),n,n)
  ENDDO

     call II_get_carmom(lupri,luerr,ls%setting,propint,nMAT,nderiv)

! Free what we do not need
     call mat_free(propint(1))
     call mat_free(propint(6))
     call mat_free(propint(7))
     call mat_free(propint(9))

! SEC = SECX + SECY + SECZ
     call mat_daxpy(1d0,propint(8),propint(5))
     call mat_free(propint(8))
     call mat_daxpy(1d0,propint(10),propint(5))
     call mat_free(propint(10))

! SEC, prepare T and PROPT=SEC*T 
     call mat_init(PROPT,n,n)
     call mat_mul(propint(5),CMO,'n','n',1d0,0d0,PROPT)
     call mat_free(propint(5))

     allocate(PROPTf(n,n))
     call mat_to_full(PROPT,1d0,PROPTf)
     call mat_free(PROPT)

     allocate(T(n,n))
     call mat_to_full(CMO,1d0,T)

     allocate(spread(n))

 do i=1,n
     sec = ddot(n,T(1,i),1,PROPTf(1,i),1)
     spread(i) = sec
 enddo


!DIPX
     call mat_init(PROPT,n,n)
     call mat_mul(propint(2),CMO,'n','n',1d0,0d0,PROPT)
     call mat_free(propint(2))

     call mat_to_full(PROPT,1d0,PROPTf)
     call mat_free(PROPT)

 do i=1,n
     dip = ddot(n,T(1,i),1,PROPTf(1,i),1)
     spread(i) = spread(i) - dip*dip
 enddo

!DIPY
     call mat_init(PROPT,n,n)
     call mat_mul(propint(3),CMO,'n','n',1d0,0d0,PROPT)
     call mat_free(propint(3))

     call mat_to_full(PROPT,1d0,PROPTf)
     call mat_free(PROPT)

 do i=1,n
     dip = ddot(n,T(1,i),1,PROPTf(1,i),1)
     spread(i) = spread(i) - dip*dip
 enddo

!DIPZ
     call mat_init(PROPT,n,n)
     call mat_mul(propint(4),CMO,'n','n',1d0,0d0,PROPT)
     call mat_free(propint(4))

     call mat_to_full(PROPT,1d0,PROPTf)
     call mat_free(PROPT)

 do i=1,n
     dip = ddot(n,T(1,i),1,PROPTf(1,i),1)
     spread(i) = spread(i) - dip*dip
     
     spread(i) = sqrt(spread(i))
 enddo

 deallocate(T,PROPTf)

 mx = spread(idamax(n,spread,1))

 deallocate(spread)

end subroutine leastchangeOrbspreadStandalone



subroutine leastchange_lcv(decomp,CMO,nocc,ls)
use leastchange_module
use decompMod
implicit none
type(decompItem) :: decomp
type(Matrix) :: CMO
type(lsitem) :: ls
type(Matrix) :: tmp
integer :: n, nocc, i
integer, allocatable :: uindex(:)
real(realk),allocatable :: CMOf(:,:), Tf(:,:)

  !convert CMO to orthonormal basis
  n = CMO%nrow
  call mat_init(tmp,n,n)
  call mat_mul(decomp%U,CMO,'n','n',1d0,0d0,tmp)
  allocate (CMOf(n,n))
  call mat_to_full(tmp,1d0,CMOf)
  call mat_free(tmp)


  !sorting index
  allocate(uindex(2*n))
  do i=1,n
    uindex(i)=i
  enddo

  !occupation sorting
  call leastchange_rowsort_occup(uindex,CMOf,nocc,n)

  allocate (Tf(n,n))

  !lcv basis
  call leastchange_blocktransf(Tf,CMOf,n,nocc)

  !inverse sorting index
  call leastchange_iuindex(uindex(n+1),uindex,n)

  !reorder to original
  call leastchange_rowreorder(CMOf,Tf,uindex(n+1),n)


  !back to gcscf basis
  deallocate(Tf)

  call mat_init(tmp,n,n)
  call mat_set_from_full(CMOf,1d0,tmp)

  deallocate(CMOf)
 
  call mat_mul(decomp%U_inv,tmp,'n','n',1d0,0d0,CMO)
  
  call mat_free(tmp)


end subroutine leastchange_lcv

subroutine leastchange_lcm(decomp,CMO,nocc,ls)
use leastchange_module
use decompMod
implicit none
type(decompItem),intent(in) :: decomp
type(Matrix) :: CMO
type(lsitem) :: ls
type(Matrix) :: tmp
integer :: n, nocc
real(realk),allocatable :: CMOf(:,:), Tf(:,:)

  !convert CMO to orthonormal basis
  n = CMO%nrow
  call mat_init(tmp,n,n)
  call mat_mul(decomp%U,CMO,'n','n',1d0,0d0,tmp)
  allocate (CMOf(n,n))
  call mat_to_full(tmp,1d0,CMOf)
  call mat_free(tmp)

  allocate (Tf(n,n))

  !lcm basis
  call leastchange_blocktransf(Tf,CMOf,n,nocc)

  !back to gcscf basis

  call mat_init(tmp,n,n)
  call mat_set_from_full(Tf,1d0,tmp)

  deallocate(CMOf,Tf)
 
  call mat_mul(decomp%U_inv,tmp,'n','n',1d0,0d0,CMO)
  
  call mat_free(tmp)


end subroutine leastchange_lcm
