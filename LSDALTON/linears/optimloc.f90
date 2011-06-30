module orbspread_module
use typedef
use matrix_module
use matrix_operations

type orbspread_data
 integer :: norb
 integer :: m
 real(realk), pointer :: spread2(:)
 type(Matrix) :: R(3)
 type(Matrix) :: Q
 type(Matrix), pointer :: G
 type(Matrix) :: propint(10)
 type(Matrix) :: tmpM(4)
end type orbspread_data

contains

subroutine orbspread_free(orbspread_input)
implicit none
type(orbspread_data) :: orbspread_input
integer :: i

     deallocate(orbspread_input%spread2) 
 
     call mat_free(orbspread_input%Q)

     do i=1,3
      call mat_free(orbspread_input%R(i))

      call mat_free(orbspread_input%propint(i+1))
     enddo

     do i=1,4
      call mat_free(orbspread_input%tmpM(i))
     enddo

end subroutine orbspread_free


subroutine orbspread_init(orbspread_input,norb,ls)
implicit none
type(orbspread_data), intent(inout) :: orbspread_input
integer             , intent(in) :: norb
TYPE(lsitem) :: ls

integer, parameter  :: nderiv=2, nMAT=10
integer i,n

! init and compute propint property integrals
! propint(2:4) -> DIPX,DIPY,DIPZ
! propint(5)   -> SECX+SECY+SECZ 

     n   = ls%setting%BASIS(1)%p%REGULAR%nbast

     DO i=1,nMAT
        call mat_init(orbspread_input%propint(I),n,n)
     ENDDO

     call II_get_carmom(6,6,ls%setting,orbspread_input%propint,nMAT,nderiv)

     call mat_free(orbspread_input%propint(1))
     call mat_free(orbspread_input%propint(6))
     call mat_free(orbspread_input%propint(7))
     call mat_free(orbspread_input%propint(9))

! propint(5)   -> SECX+SECY+SECZ 
     call mat_daxpy(1d0,orbspread_input%propint(8),orbspread_input%propint(5))
     call mat_free(orbspread_input%propint(8))
     call mat_daxpy(1d0,orbspread_input%propint(10),orbspread_input%propint(5))
     call mat_free(orbspread_input%propint(10))


! init R Q and tmpM matrices
     call mat_init(orbspread_input%Q,norb,norb)
     do i=1,3
     call mat_init(orbspread_input%R(i),norb,norb)
     enddo
     do i=1,4
     call mat_init(orbspread_input%tmpM(i),norb,norb)
     enddo

! set norb
     orbspread_input%norb =  norb

! set m power
     orbspread_input%norb =  2

! allocate spread2

     allocate(orbspread_input%spread2(norb))

end subroutine orbspread_init


subroutine orbspread_update(orbspread_input,CMO)
implicit none
type(orbspread_data), intent(inout) :: orbspread_input
type(Matrix), intent(in) :: CMO

type(Matrix) :: tmp
integer :: nbas,norb,i
real(realk), allocatable :: tmpv(:)

     nbas=CMO%nrow
     norb=CMO%ncol

     call mat_init(tmp,nbas,norb)

!   R(1:3)
    do i=1,3
     call mat_mul(orbspread_input%propint(i+1),CMO,'n','n',1d0,0d0,tmp)
     call mat_mul(CMO,tmp,'t','n',1d0,0d0,orbspread_input%R(i))
    enddo

!   Q
     call mat_mul(orbspread_input%propint(5),CMO,'n','n',1d0,0d0,tmp)
     call mat_mul(CMO,tmp,'t','n',1d0,0d0,orbspread_input%Q)

     call mat_free(tmp)

!   spread2
    allocate(tmpv(norb))

    call mat_extract_diagonal(orbspread_input%spread2,orbspread_input%Q)

    do i=1,3
    call mat_extract_diagonal(tmpv,orbspread_input%R(i))
    call daxpy(norb,-1d0,tmpv**2,1,orbspread_input%spread2,1)
    enddo
   
    deallocate(tmpv)

end subroutine orbspread_update

subroutine orbspread_localize(CMO,ls)
implicit none
type(Matrix) , intent(inout ):: CMO
TYPE(lsitem) , intent(in)    :: ls

type(orbspread_data) :: orbspread_input
type(Matrix) :: X, G
integer :: norb, i



  norb=CMO%ncol

  call mat_init(X,norb,norb)
  call mat_init(G,norb,norb)

  call orbspread_init(orbspread_input,norb,ls)

  call orbspread_update(orbspread_input,CMO)

  call orbspread_gradx(G,norb,orbspread_input)

  do i=1,100
  
    if((sqrt(mat_sqnorm2(G))/norb).le.1e-2) exit

    !X -> solver(orbspread_hesslin, G)

    call orbspread_updatecmo(CMO,X)

    call orbspread_update(orbspread_input,CMO)

    call orbspread_gradx(G,norb,orbspread_input)

  enddo

  call orbspread_free(orbspread_input)
  
end subroutine orbspread_localize

subroutine orbspread_updatecmo(CMO,X)
use matrix_util
implicit none
type(Matrix), intent(inout) :: CMO
type(Matrix), intent(in)    :: X

integer    ::  norb,nbas
type(Matrix) :: expX, tmp

    norb=CMO%ncol
    nbas=CMO%nrow

    call mat_init(expX,norb,norb)
 
    call matrix_exponential(X,expX,1d-8)

    call mat_init(tmp,nbas,norb)

    call mat_mul(CMO,expX,'n','n',1d0,0d0,tmp)

    call mat_free(expX)

    call mat_assign(CMO,tmp)

    call mat_free(tmp)

end subroutine orbspread_updatecmo


subroutine orbspread_gradx(G,norb,inp)
implicit none
Type(Matrix) , target, intent(inout) :: G
integer                              :: norb
type(orbspread_data), intent(inout)  :: inp

real(realk)  :: diagR(norb,3), tmp(norb)
integer      :: x, m

  m=inp%m

  do x=1, 3
  call mat_extract_diagonal(diagR(:,x),inp%R(x))
  enddo

  tmp =  (inp%spread2**(m-1))
  call mat_dmul(tmp,inp%Q,'n',-2d0*m,0d0,G)

  do x=1,3
  tmp = diagR(:,x)*(inp%spread2**(m-1))
  call mat_dmul(tmp,inp%R(x),'n',4d0*m,1d0,G)
  enddo

  call mat_trans(G,inp%tmpM(1))
  call mat_daxpy(-1d0,inp%tmpM(1),G)

  inp%G => G

end subroutine  orbspread_gradx
  
subroutine orbspread_precond(Pout,mu,norb,inp)
type(Matrix), intent(inout) :: Pout
real(realk), intent(in)   :: mu
integer, intent(in)       :: norb
type(orbspread_data), intent(in) :: inp

real(realk)  :: P(norb,norb)
real(realk)  :: diagQ(norb), diagR(norb,3), spm2(norb),spm1(norb)
real(realk)  :: Qkl, Rklx, elm, tmpk, tmpl
integer      :: x,k,l,m

  m=inp%m


  call mat_extract_diagonal(diagQ,inp%Q)

  do x=1, 3
  call mat_extract_diagonal(diagR(:,x),inp%R(x))
  enddo

 spm1 = inp%spread2**(m-1)
 spm2 = inp%spread2**(m-2)
 

 do k=1,norb
 do l=k,norb
     call mat_get_elm(inp%Q,k,l,Qkl) 

     P(k,l) = 4*m*(m-1)*(spm2(k) +  spm2(l))*Qkl*Qkl

     P(k,l) = P(k,l) + 2*m*(diagQ(k)-diagQ(l))*(spm1(l) - spm1(k));

     tmpk = 0d0; tmpl= 0d0

     do x=1,3
     call mat_get_elm(inp%R(x),k,l,Rklx) 
     P(k,l) = P(k,l) + 4*m*spm1(k)*(diagR(k,x)*diagR(k,x)-diagR(k,x)*diagR(l,x));
     P(k,l) = P(k,l) + 4*m*spm1(l)*(diagR(l,x)*diagR(l,x)-diagR(k,x)*diagR(l,x));

      P(k,l) = P(k,l) - 8*m*(spm1(k) + spm1(l))*Rklx*Rklx
      P(k,l) = P(k,l) -16*m*(m-1)*Qkl*(spm2(k)*diagR(k,x) + spm2(l)*diagR(l,x))*Rklx

      tmpk = tmpk + diagR(k,x)*Rklx;
      tmpl = tmpl + diagR(l,x)*Rklx;

     enddo

     P(k,l) = P(k,l) + 16*m*(m-1)*(spm2(k)*tmpk*tmpk + spm2(l)*tmpl*tmpl)


     P(k,l)= P(k,l)-mu

     P(l,k)=P(k,l)
    
 enddo
 enddo


 call mat_set_from_full(P,1d0,Pout)


end subroutine orbspread_precond



end module orbspread_module
!> \brief hessian linear transformation for power m orbital spread locality measure
!> \author B. Jansik
!> \date 2010
!> \param Hv linear transformation of H with trial vector V
!> \param V, trial vector
!> \param mu, level shift (H-muI)
!> \param norb, size of matrices (same as number of orbitals involved in localization)
!> \param orbspread_input structure holding orbspread related data
!!> \param m power
!!> \param spread2 vector of squared orbital spreads 
!!> \param R(3) matrix array holding XDIPLEN,YDIPLEN,ZDIPLEN components
!!> \param Q matrix holding sum XXSECMOM+YYSECMOM+ZZSECMOM
!!> \param tmpM(4) matrix array used as workspace
!> all parameters must be allocated prior to subroutine call
!> all type(matrix) arguments are of norb,norb size
subroutine orbspread_hesslin(Hv,V,mu,norb,orbspread_input)!m,spread2,R,Q,tmpM)
use orbspread_module
implicit none

Type(Matrix), intent(inout) :: Hv
Type(Matrix), intent(in)  :: V
real(realk), intent(in)   :: mu
integer, intent(in)       :: norb
type(orbspread_data), intent(in), target :: orbspread_input

integer       :: m
Type(Matrix)  :: R(3), RV(3)
Type(Matrix), pointer  ::  Q, G, QV ,tmpM
real(realk), pointer   :: spread2(:)
real(realk)  :: diagQV(norb), diagR(norb,3), diagRV(norb,3), tmp(norb)
integer      :: x,y,i

 !pointer assignments
 !stupid f90, no support for pointer arrays!
  do i=1,3
   call mat_clone(R(i),orbspread_input%R(i))
   call mat_clone(RV(i),orbspread_input%tmpM(i))
  enddo
  Q       => orbspread_input%Q
  QV      => orbspread_input%tmpM(4)
  tmpM    => orbspread_input%tmpM(4) !it's ok, QV will not be needed at that point
  G       => orbspread_input%G

  spread2 => orbspread_input%spread2

  m       =  orbspread_input%m


  !job
  call mat_mul(Q,V,'n','n',1d0,0d0,QV)

  call mat_extract_diagonal(diagQV,QV)


  do x=1, 3
  call mat_mul(R(x),V,'n','n',1d0,0d0,RV(x))
  call mat_extract_diagonal(diagRV(:,x),RV(x))
  call mat_extract_diagonal(diagR(:,x),R(x))
  enddo


  tmp = diagQV*(spread2**(m-2))
  call mat_dmul(tmp,Q,'n',-4d0*m*(m-1),0d0,Hv)


  tmp =  (spread2**(m-1))
  call mat_dmul(tmp,QV,'t',-2d0*m,1d0,Hv)

  call mat_dmul(tmp,QV,'n',-2d0*m,1d0,Hv)

  do x=1, 3
  tmp = diagR(:,x)*diagQV*(spread2**(m-2))
  call mat_dmul(tmp,R(x),'n',8d0*m*(m-1),1d0,Hv)

  tmp = diagR(:,x)*diagRV(:,x)*(spread2**(m-2))
  call mat_dmul(tmp,Q,'n',8d0*m*(m-1),1d0,Hv)
 
     do y=1, 3
     tmp = diagR(:,x)*diagR(:,y)*diagRV(:,x)*(spread2**(m-2))
     call mat_dmul(tmp,R(y),'n',-16d0*m*(m-1),1d0,Hv)
     enddo

  tmp = diagRV(:,x)*(spread2**(m-1))
  call mat_dmul(tmp,R(x),'n',8d0*m,1d0,Hv)

  tmp = diagR(:,x)*(spread2**(m-1))
  call mat_dmul(tmp,RV(x),'t',4d0*m,1d0,Hv)

  call mat_dmul(tmp,RV(x),'n',4d0*m,1d0,Hv)

  enddo

  call mat_mul(V,G,'n','n',0.5d0,1d0,Hv)

  call mat_trans(Hv,tmpM)
  call mat_daxpy(-1d0,tmpM,Hv)

  call mat_scal(2d0,Hv)

  if (mu.ne.0d0) call mat_daxpy(-mu,V,Hv)

end subroutine orbspread_hesslin

!> \brief unitest for orbspread_hesslin() subroutine
!> \author B. Jansik
!> \date 2010
!> \param passed 
!> passed shold be set to .true. otherwise orbspread_hesslin()
!> must be considered broken.
subroutine orbspread_hesslin_unitest(passed)
use orbspread_module
implicit none
logical, intent(out)   :: passed
integer,parameter      :: m=3,norb=4

type(orbspread_data)   :: inp

type(Matrix) :: Hv,V,T, G
real(realk),parameter :: mu=-1.32d0
integer :: i

!allocations
 call mat_init(Hv,norb,norb)
 call mat_init(V,norb,norb)
 call mat_init(inp%Q,norb,norb)
 call mat_init(T,norb,norb)
 call mat_init(G,norb,norb)
 do i=1,3
  call mat_init(inp%R(i),norb,norb)
 enddo
 do i=1,4
  call mat_init(inp%tmpM(i),norb,norb)
 enddo

 allocate(inp%spread2(norb))

 inp%m = m

!initializations
V%elms=(/  0.00000,   0.10165,  -0.20779,   0.80376,&
&         -0.10165,   0.00000,   0.43667,   0.67566,&
&          0.20779,  -0.43667,   0.00000,  -0.01383,&
&         -0.80376,  -0.67566,   0.01383,   0.00000/)


inp%spread2=(/0.16223,  0.14987,  0.39164,  0.61443/)


inp%R(1)%elms=(/  0.734535,  0.584218,  0.701755,  0.064086,&
&             0.853723,  0.796857,  0.305236,  0.681963,&
&             0.625375,  0.441937,  0.741858,  0.748711,&
&             0.766346,  0.463687,  0.927562,  0.019981/)

inp%R(2)%elms=(/  0.557369,  0.435741,  0.736332,  0.455053,&
&             0.027335,  0.545388,  0.923669,  0.803679,&
&             0.984998,  0.065604,  0.471317,  0.816302,&
&             0.662640,  0.481126,  0.147809,  0.717182/)

inp%R(3)%elms=(/  0.6461293,  0.6788147,  0.6815679,  0.1965552,&
&             0.2743679,  0.7154769,  0.3921278,  0.2958226,&
&             0.6670042,  0.3111366,  0.8163864,  0.0074360,&
&             0.7020067,  0.8659539,  0.7503724,  0.0757745/)

inp%Q%elms=(/0.047692,  0.495820,  0.258083,  0.845362,&
&        0.423167,  0.148130,  0.953802,  0.640069,&
&        0.588172,  0.639629,  0.395089,  0.084721,&
&        0.244598,  0.104311,  0.949004,  0.040549/)

T%elms=(/       0.00000,   -6.50228,   25.47338,   23.90821,&
&               6.50228,    0.00000,    6.68831,   14.07902,&
&             -25.47338,   -6.68831,    0.00000,  -10.29937,&
&             -23.90821,  -14.07902,   10.29937,    0.00000/)


!test
 call orbspread_gradx(G,norb,inp)

 call orbspread_hesslin(Hv,V,mu,norb,inp)


 call mat_daxpy(-1d0,Hv,T)

 passed = mat_sqnorm2(T).le.1d-9

!deallocations
 call mat_free(Hv)
 call mat_free(V)
 call mat_free(inp%Q)
 call mat_free(T)
 call mat_free(G)
 do i=1,3
  call mat_free(inp%R(i))
 enddo
 do i=1,4
  call mat_free(inp%tmpM(i))
 enddo

 deallocate(inp%spread2)

end subroutine orbspread_hesslin_unitest


!> \brief unitest for orbspread_gradx() subroutine
!> \author B. Jansik
!> \date 2010
!> \param passed 
!> passed shold be set to .true. otherwise orbspread_gradx()
!> must be considered broken.
subroutine orbspread_gradx_unitest(passed)
use orbspread_module
implicit none
logical, intent(out)   :: passed
integer,parameter      :: m=3,norb=4

type(orbspread_data)   :: inp

type(Matrix) :: G,T
integer :: i

!allocations
 call mat_init(G,norb,norb)
 call mat_init(inp%Q,norb,norb)
 call mat_init(T,norb,norb)
 do i=1,3
  call mat_init(inp%R(i),norb,norb)
 enddo
 do i=1,4
  call mat_init(inp%tmpM(i),norb,norb)
 enddo

 allocate(inp%spread2(norb))

 inp%m = m

inp%spread2=(/0.16223,  0.14987,  0.39164,  0.61443/)


inp%R(1)%elms=(/  0.734535,  0.584218,  0.701755,  0.064086,&
&             0.853723,  0.796857,  0.305236,  0.681963,&
&             0.625375,  0.441937,  0.741858,  0.748711,&
&             0.766346,  0.463687,  0.927562,  0.019981/)

inp%R(2)%elms=(/  0.557369,  0.435741,  0.736332,  0.455053,&
&             0.027335,  0.545388,  0.923669,  0.803679,&
&             0.984998,  0.065604,  0.471317,  0.816302,&
&             0.662640,  0.481126,  0.147809,  0.717182/)

inp%R(3)%elms=(/  0.6461293,  0.6788147,  0.6815679,  0.1965552,&
&             0.2743679,  0.7154769,  0.3921278,  0.2958226,&
&             0.6670042,  0.3111366,  0.8163864,  0.0074360,&
&             0.7020067,  0.8659539,  0.7503724,  0.0757745/)

inp%Q%elms=(/0.047692,  0.495820,  0.258083,  0.845362,&
&        0.423167,  0.148130,  0.953802,  0.640069,&
&        0.588172,  0.639629,  0.395089,  0.084721,&
&        0.244598,  0.104311,  0.949004,  0.040549/)

T%elms=(/    0.00000,   0.06159,   2.02191,  -0.76216,&
&           -0.06159,   0.00000,   0.85115,   1.00137,&
&           -2.02191,  -0.85115,   0.00000,   0.88169,&
&            0.76216,  -1.00137,  -0.88169,   0.00000/)


!test

 call orbspread_gradx(G,norb,inp)

 call mat_daxpy(-1d0,G,T)

 passed = mat_sqnorm2(T).le.1d-9

!deallocations
 call mat_free(G)
 call mat_free(inp%Q)
 call mat_free(T)
 do i=1,3
  call mat_free(inp%R(i))
 enddo
 do i=1,4
  call mat_free(inp%tmpM(i))
 enddo

 deallocate(inp%spread2)

end subroutine orbspread_gradx_unitest

subroutine orbspread_precond_unitest(passed)
use orbspread_module
implicit none
logical, intent(out)   :: passed
integer,parameter      :: m=3,norb=4
real(realk),parameter :: mu=0d0
type(orbspread_data)   :: inp

type(Matrix) :: P,T
integer :: i

!allocations
 call mat_init(P,norb,norb)
 call mat_init(inp%Q,norb,norb)
 call mat_init(T,norb,norb)
 do i=1,3
  call mat_init(inp%R(i),norb,norb)
 enddo
 do i=1,4
  call mat_init(inp%tmpM(i),norb,norb)
 enddo

 allocate(inp%spread2(norb))

 inp%m = m

inp%spread2=(/0.16223,  0.14987,  0.39164,  0.61443/)


inp%R(1)%elms=(/   0.22097, 0.75410, 0.46831, 0.36425,&
&                  0.75410, 0.22471, 0.49422, 0.48281,&
&                  0.46831, 0.49422, 0.88116, 0.41692,&
&                  0.36425, 0.48281, 0.41692, 0.60770/)

inp%R(2)%elms=(/    0.27603,  0.79741,  0.46847,  0.15278,&
&                   0.79741,  0.74890,  0.59191,  0.94024,&
&                   0.46847,  0.59191,  0.35630,  0.60522,&
&                   0.15278,  0.94024,  0.60522,  0.95474/)

inp%R(3)%elms=(/   0.83232,  0.48356,  0.60922,  0.40515,&
&                  0.48356,  0.45738,  0.33015,  0.28682,&
&                  0.60922,  0.33015,  0.28448,  0.56644,&
&                  0.40515,  0.28682,  0.56644,  0.26361/)

inp%Q%elms=(/  0.14532,  0.47040,  0.67501,  0.47688,&
&              0.47040,  0.60246,  0.63225,  0.75969,&
&              0.67501,  0.63225,  0.90917,  0.55272,&
&              0.47688,  0.75969,  0.55272,  0.16080/)

T%elms=(/16.26034, 11.36277,  5.86034,   4.28595,&
&        11.36277,  6.87555,  6.51157,  42.06788,&
&         5.86034,  6.51157, 13.85434,  29.38276,&
&         4.28595, 42.06788, 29.38276, 165.78081/)

!test

 call orbspread_precond(P,mu,norb,inp)

 call mat_daxpy(-1d0,P,T)

 passed = mat_sqnorm2(T).le.1d-9

!deallocations
 call mat_free(P)
 call mat_free(inp%Q)
 call mat_free(T)
 do i=1,3
  call mat_free(inp%R(i))
 enddo
 do i=1,4
  call mat_free(inp%tmpM(i))
 enddo

 deallocate(inp%spread2)

end subroutine orbspread_precond_unitest

