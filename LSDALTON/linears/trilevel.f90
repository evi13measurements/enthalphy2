!> @file 
!> Contains the trilevel and atoms starting guess in addition to the construction of the GCbasis se PCCP 2009, 11, 5805-5813 
!> Trilevel and Atoms module
!> \author Branislav Jansik documented by \latexonly T. Kj{\ae}rgaard\endlatexonly
!> \date 2010-03-03
module trilevel_module
use typedef
use Matrix_Operations

type trilevel_atominfo
     integer :: LUPRI
     integer :: ND,NA ! number of disticnt atoms, number of atoms
     integer, pointer :: NATOM(:)!length ND : FOR A GIVEN UNIQUE ATOM THIS POINTS TO A ATOM OF SAME TYPE
     integer, pointer :: UATOMTYPE(:) !length ND : FOR A GIVEN UNIQUE ATOM THIS IS THE ATOMTYPE IN FULL 
end type trilevel_atominfo

contains

!> \brief diagonalize a angular moment block of the atomic h1 matrix 
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param angular moment
!> \param basis_size size of basis
!> \param bCMO MO coefficients
!> \param F full Fock matrix 
!> \param S full Overlap matrix 
!> \param nbast full dimensions 
subroutine trilevel_diag_per_ang(ang,basis_size,bCMO,F,S,nbast)
implicit none

integer, intent(in)     :: ang, basis_size(:), nbast
real(realk), intent(in) :: F(nbast,nbast), S(nbast,nbast)
real(realk),target      :: bCMO( basis_size(ang+1), basis_size(ang+1))
!
real(realk), pointer    :: bF(:,:), bS(:,:), eig(:), wrk(:)
integer                 :: i, j, k, istart, info, nb, lwrk
integer,     pointer    :: indexlist(:) 
 
 nb =  basis_size(ang+1)

 indexlist => trilevel_indexlist(ang,basis_size)

!create subblock matrix
 allocate(bS(nb,nb))
 bF => bCMO
 !allocate(bCMO(basis_size(ang+1),basis_size(ang+1)))

 do i=1, nb
  do j=1, nb
       bF(i,j)=F(indexlist(i),indexlist(j))
       bS(i,j)=S(indexlist(i),indexlist(j))
  enddo
 enddo

! diagonalize
!querry optimal work size

 if (nb.gt.1) then
  allocate(eig(nb))
  call dsygv(1,'V','U',nb,bF,nb,bS,nb,eig,eig,-1,info)
  
!run diagonalization
  lwrk = eig(1)
  allocate(wrk(lwrk))
  call dsygv(1,'V','U',nb,bF,nb,bS,nb,eig,wrk,lwrk,info)
 
  deallocate(eig,wrk)
 else
  bCMO(1,1)=1d0
 endif


!cleanup
deallocate(indexlist,bS)

end subroutine trilevel_diag_per_ang

!> \brief 
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param opt Contains info about SCF optimization
!> \param D density matrix
!> \param H1 one electron contribution to the fock matrix
!> \param F fock matrix
!> \param Etotal energy 
!> \param newlupri logical unit number for output
!> \param newluerr logical unit number for error output
!> \param ls lsitem structure containing integral info
   SUBROUTINE trilevel_get_fock(opt,D,H1,F,Etotal,newlupri,newluerr,setting)
   ! ===================================================================
   ! di_get_fock obtains total fock matrix and corresponding energy.
   ! WE have to go through the interface to dalton before the fock
   ! evaluator learns how to handle arbitrary-type arrays.
   ! ===================================================================
      use opttype
      IMPLICIT NONE
      type(optItem),intent(in)   :: opt
      TYPE(Matrix), INTENT(IN)   :: D 
      TYPE(Matrix),intent(inout) :: F
      type(lssetting),intent(inout) :: setting
      real(realk), INTENT(OUT) :: Etotal
      TYPE(Matrix),intent(in)  :: H1
      integer,intent(in) :: newlupri,newluerr
      real(realk)   :: edfty, edfty_a, edfty_b
      integer nbast
      logical  :: Dsym

!     Two-electron part: G(D)
      IF(.not.opt%cfg_unres) THEN !closed shell case
!        Coulomb and exchange
         Dsym = .TRUE.!symmetric Density matrix
         call II_get_Fock_mat(newlupri,newluerr,setting,D,Dsym,F,1)
         Etotal = fockenergy_f(opt%cfg_unres,F,D,H1)
!        Exchange-correlation
         !IF(cfg_isdft) THEN
         if (opt%calctype == opt%dftcalc) then
            nbast = D%nrow
            call II_get_xc_fock_mat(newlupri,newluerr,setting,nbast,D,F,Edfty)
            Etotal = Etotal + Edfty
         ENDIF
      ELSE               !unrestricted open shell case
         Dsym = .TRUE. !symmetric Density matrix
         call II_get_Fock_mat(newlupri,newluerr,setting,D,Dsym,F,1)
         Etotal = fockenergy_f(opt%cfg_unres,F,D,H1)
         if (opt%calctype == opt%dftcalc) then
            nbast = D%nrow
            call II_get_xc_fock_mat(newlupri,newluerr,setting,nbast,D,F,Edfty)
            Etotal = Etotal + Edfty
         ENDIF
      ENDIF
!     Add one-electron part: F(D) = h + G(D)
      call mat_daxpy(1d0,H1,F)

   contains      
      double precision function fockenergy_F(unres,F,D,H1)
         !E = Tr(h + F)D + hnuc ! No factor since molecular orbitals
         implicit none
         logical, intent(in) :: unres
         TYPE(matrix), intent(in) :: F,D,H1
         double precision :: hD, FD, fac
         double precision, external :: DF_E_ROB_CORR
         integer :: ndim

         ndim=F%nrow
         fac=2d0
         if(unres) fac=1d0
         !Get the one-electron hamiltonian
         !Tr(FD)
         !Tr(hD)
         hD =       mat_dotproduct(D,H1)
         FD = 0.5d0*mat_dotproduct(D,F)
         !E(HF) = Tr(h +F)D + hnuc
         fockenergy_F = (hd + FD)*fac !Stinne: Remove call to boxed DF contribution (obsolete), 08-06-2010
      end function fockenergy_F

   END SUBROUTINE trilevel_get_fock

!> \brief the grand canonical SCF loop
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param opt optItem containing info about scf optimization
!> \param D density matrix
!> \param CMO MO coefficients
!> \param H1 one electron contribution to the fock matrix
!> \param F fock matrix
!> \param S Overlap matrix
!> \param ai trilevel_atominfo structure containing info about the unique atoms
!> \param ls lsitem structure containing integral,molecule,basis info
!> \param iatom unique atom index
!> \param newlupri logical unit number for output
!> \param newluerr logical unit number for error output
SUBROUTINE trilevel_gcscfloop(opt,D,CMO,H1,F,S,ai,setting,molecule,basis,iatom,newlupri,newluerr)
  use av_utilities
  use matrix_util
  use opttype
  USE scf_stats
  use linsca_diis 

   implicit none
   type(optItem),intent(in) :: opt
   type(trilevel_atominfo),intent(in) :: ai
   integer,intent(in) :: iatom
   type(lssetting),intent(inout) :: setting
   type(moleculeinfo),intent(in) :: molecule
   type(basisinfo),intent(in)    :: basis
   TYPE(Matrix),intent(inout)   :: D, F, S, H1, CMO
   integer,intent(in) :: newlupri,newluerr
!
   TYPE(Matrix)            :: grad
   TYPE(util_HistoryStore) :: queue
   real(realk) :: E, gradnrm
   integer     :: iteration
   logical :: energy_converged,firsttime
   integer :: itype,nbast,iAO
   type(avItem) :: av
   type(moleculeinfo),target :: atomicmolecule

   itype = ai%UATOMTYPE(iatom)
   nbast = molecule%atom(ai%NATOM(iatom))%nContOrbREG

   call av_set_default_config(av)
   av%lupri = opt%lupri
   av%CFG_averaging = av%CFG_AVG_DIIS
   av%trilevel_gcscf = .true.

   if (opt%calctype == opt%dftcalc) then
      av%cfg_set_type = av%CFG_THR_dft
      av%diis_history_size  = av%cfg_settings(av%CFG_SET_type)%max_history_size
   endif
   call queue_init(av,queue)

   call scf_stats_init(opt)

   call mat_init(grad,H1%nrow,H1%nrow)

!
! GCSCF iterations
!
   CALL build_atomicmolecule(molecule,atomicmolecule,ai%NATOM(iatom),opt%lupri)
   do iAO=1,4
      setting%molecule(iAO)%p => atomicmolecule
      setting%fragment(iAO)%p => atomicmolecule
   enddo

   DO iteration = 1, opt%cfg_max_linscf_iterations

      call trilevel_get_fock(opt, D, H1, F, E, &
           &newlupri,newluerr,setting)

      call get_AO_gradient(F, D, S, grad)

      gradnrm = sqrt(mat_sqnorm2(grad))

      call scf_stats_update(iteration,gradnrm,E,opt)

      IF(gradnrm < 7.5d-4) then 
            EXIT
      ENDIF

      CALL add_to_queue(av, F, D, S, E, grad, queue) 

      call diis(av,queue,D,F)

      CALL trilevel_get_density_blocked(D,CMO,F,S,basis,itype,nbast)
   END DO

   print *, "gcscf loop done"
  call queue_free(av,queue)
  CALL mat_free(grad)
  call scf_stats_shutdown
  call free_moleculeinfo(atomicmolecule)

END SUBROUTINE trilevel_gcscfloop

!> \brief 
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param D density matrix
!> \param CMO MO coefficients 
!> \param F fock matrix
!> \param S Overlaps matrix
!> \param ls lsitem containing info about the integral evaluation and basis,..
!> \param itype atom type of the current atom of interest
!> \param nbast number of basis functions for this atom
subroutine trilevel_get_density_blocked(D,CMO,F,S,basis,itype,nbast)
implicit none
type(Matrix)           ::  F,D,S,CMO
type(basisinfo)        :: basis
integer                ::  nAngmom, kmult,ipos, ang, nbast,nb,itype
real(realk),allocatable::  bCMO(:,:),bD(:,:),occ(:), Fmat(:), Smat(:)
integer, pointer       ::  perm(:), iperm(:), basis_size(:)


  nAngmom = BASIS%REGULAR%ATOMTYPE(itype)%nAngmom

  basis_size=>trilevel_set_basis_size(BASIS,itype) !size of basis for a given angmom

  !build list of permutation to obtain (1S,2S,3S,2Px,3Px,2Py,3Py,2Pz,3Pz,3Pz,...)
  perm=>trilevel_blockdiagonal_permutation(nAngmom,nbast,basis_size) 
  !sorting
  iperm=>trilevel_isort(perm,nbast)

  allocate(occ(nbast))
  occ = 0d0;
  CMO%elms=0d0;

  kmult=1; ipos =1
  do ang=0,nAngmom-1  

     nb=basis_size(ang+1)
     if (nb.ne.0) then

     allocate(bCMO(nb,nb))
     !diagonal angular block => bCMO
     call trilevel_diag_per_ang(ang,basis_size,bCMO,F%elms,S%elms,nbast)
     !set occupation numbers based on element table in ecdata
     call trilevel_set_occ(occ(ipos:nbast),ang,BASIS%REGULAR%ATOMTYPE(itype)%Charge)
     !Build full MO coefficient matrix from the bCMO blocks
     call trilevel_mdiag(CMO%elms,ipos,nbast,bCMO,nb,kmult)

     deallocate(bCMO)

     endif
     kmult = kmult + 2
  enddo
  !reorder the MO coefficient matrix
  call trilevel_reorder2d(CMO%elms,nbast,iperm)
  !build Density matrix eq. 7 from article
  call trilevel_density_from_orbs(D%elms,occ,CMO%elms,nbast)
  
  deallocate(basis_size,perm,iperm,occ)
  nullify(basis_size,perm,iperm)

end subroutine trilevel_get_density_blocked

!> \brief build trilevel_blockdiagonal_permutation  
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param nAngmom number of angular moments
!> \param nbast number of basis functions
!> \param basis_size size of basis
function trilevel_blockdiagonal_permutation(nAngmom,nbast,basis_size)
implicit none
integer, pointer :: trilevel_blockdiagonal_permutation(:)
integer          :: nAngmom, nbast
integer          :: basis_size(nAngmom)
integer          :: i,j,k,l,istart,iend,ioff

  allocate(trilevel_blockdiagonal_permutation(nbast))

  k=1;l=1
  istart=1
  do i=1,nAngmom
   iend = istart  + k*basis_size(i) -1
   do ioff=0,(k -1)
   do j=(istart+ioff),iend,k
      trilevel_blockdiagonal_permutation(l) = j 
      l=l+1
   enddo
   enddo
   istart= iend+1
   k=k+2
  enddo

return
end function trilevel_blockdiagonal_permutation

!> \brief sorting routine
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param arr  vector of integers to be sorted 
!> \param n length of vector
function trilevel_isort(arr,n)
implicit none
integer, pointer :: trilevel_isort(:)
integer          :: n
integer          :: arr(n),a,b
integer          :: i,j

 allocate(trilevel_isort(n))

 do i=1,n
  trilevel_isort(i)=i
 enddo

 do j=2,n
  a=arr(trilevel_isort(j)); b=trilevel_isort(j)
  do i=j-1,1,-1
    if (arr(trilevel_isort(i)).le.a) goto 10
!   arr(i+1)=arr(i)
    trilevel_isort(i+1)=trilevel_isort(i)
  enddo
  i=0
10 continue 
! arr(i+1)=a
  trilevel_isort(i+1)=b
 enddo

 return
end function trilevel_isort

!> \brief reorder routine
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param A matrix to be reordered
!> \param n dimension
!> \param perm reorder according to this vector 
subroutine trilevel_reorder2d(A,n,perm)
implicit none
integer     :: n,i
real(realk) :: A(n,n)
real(realk) :: tmp(n,n)
integer     :: perm(n)


  do i=1,n
     tmp(:,i) = A(:,perm(i))
  enddo

  do i=1,n
     A(i,:)   = tmp(perm(i),:)
  enddo


end subroutine trilevel_reorder2d

!> \brief calculates the size of the basisset
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param ls the lsitem containing the basis info
function  trilevel_set_basis_size(basis,itype)
implicit none
integer, pointer :: trilevel_set_basis_size(:)
Type(basisinfo)  :: basis
integer          :: itype
!
integer                          :: i

 allocate(trilevel_set_basis_size(BASIS%REGULAR%ATOMTYPE(itype)%nAngmom))

 do i=1, BASIS%REGULAR%ATOMTYPE(itype)%nAngmom
    trilevel_set_basis_size(i)= BASIS%REGULAR%ATOMTYPE(itype)%SHELL(i)%norb
 enddo

 return
end function trilevel_set_basis_size

!> \brief Eq. 7 from PCCP 2009, 11, 5805-5813
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param D density matrix
!> \param occ occupation numbers
!> \param CMO MO coefficients
!> \param n dimension of density matrix
subroutine trilevel_density_from_orbs(D,occ,CMO,n)
implicit none
integer     :: n, i, j
real(realk) :: occ(n), CMO(n,n), D(n,n)
real(realk) :: tmp(n,n)

 do j=1, n
      tmp(:,j)=occ(j)*CMO(:,j)
 enddo

 call dgemm('n','t',n,n,n,1d0,tmp,n,CMO,n,0d0,D,n)


end subroutine trilevel_density_from_orbs

!> \brief 
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param D output density matrix  
!> \param istart start index
!> \param nbast dimension of D
!> \param bD 
!> \param n dimension of bD
!> \param kmult
subroutine trilevel_mdiag(D,istart,nbast,bD,n,kmult)
implicit none
integer     :: istart,nbast,n,kmult,i
real(realk) :: bD(n,n),D(nbast,nbast)

 do i=1, kmult
  D(istart:istart+n-1,istart:istart+n-1)=bD
  istart=istart+n
 enddo

end subroutine trilevel_mdiag


!> \brief set occupation from elementtabel (from ecdata)
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param occ occupation number 
!> \param ang angular moment
!> \param charge charge of the atom of interest 
subroutine trilevel_set_occ(occ,ang,charge)
use EcData
implicit none
real(realk)   :: occ(:)
integer       :: ang, charge


 select case(ang)
 case(0); occ(1:ElementTable(charge)%nocc_s)=&
         &ElementTable(charge)%occ_s(1:ElementTable(charge)%nocc_s);
 case(1); occ(1:ElementTable(charge)%nocc_p)=&
         &ElementTable(charge)%occ_p(1:ElementTable(charge)%nocc_p);
 case(2); occ(1:ElementTable(charge)%nocc_d)=&
         &ElementTable(charge)%occ_d(1:ElementTable(charge)%nocc_d);
 case(3); occ(1:ElementTable(charge)%nocc_f)=&
         &ElementTable(charge)%occ_f(1:ElementTable(charge)%nocc_f);

 case default
 end select
end subroutine trilevel_set_occ

!> \brief 
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param 
!> \param 
!> \param 
!> \param 
!> \param 
function trilevel_indexlist(ang,basis_size)
implicit none
integer, pointer :: trilevel_indexlist(:)
integer          :: ang, basis_size(ang+1)
integer          :: nb, i,j, k, istart
 
 nb =  basis_size(ang+1)
 k=1;
 istart=1
 do i=1,ang
   istart = istart + k*basis_size(i)
   k = k+2;
 enddo

!create list of indexes of elements with same angular and lowest magnetic number (i.e. px1,px2,px3)
 allocate(trilevel_indexlist(nb))

 j = 1
 do i=istart, istart -1 + (k*nb), k
   trilevel_indexlist(j)=i  
   j=j+1
 enddo

return
end function trilevel_indexlist

!> \brief change the input basis in the ls%setting to the grand canonical basis Eq. 8
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param ls lsitem structure containing info about the integral eval and basis
!> \param CMO MO coefficients
!> \param ai trilevel_atominfo containing info about the unique atoms
!> \param iatom the unique atom index
subroutine trilevel_convert_ao2gcao(ls,CMO,nbast,ai,iatom)
use BUILDAOBATCH
implicit none
TYPE(lsitem) :: ls!,atomic_ls
type(trilevel_atominfo) :: ai
real(realk)         :: CMO(nbast,nbast)
real(realk),allocatable :: bCMO(:,:)
integer             :: nAngmom, nbast,ang,nb,i,j,itype,icharge
integer             :: nprim, norb,iatom,nrow,ncol
integer, pointer    :: basis_size(:), indexlist(:)
real(realk), pointer:: contrCoeffs(:)
real(realk), allocatable :: CCtmp(:,:)

  itype = ai%UATOMTYPE(iatom) !type of distinct atom in full
  nAngmom = ls%input%BASIS%REGULAR%ATOMTYPE(itype)%nAngmom

  basis_size=>trilevel_set_basis_size(ls%setting%BASIS(1)%p,itype)

 do ang=0, nAngmom-1

   nb=basis_size(ang+1)
   if (nb.eq.0) cycle

   indexlist => trilevel_indexlist(ang,basis_size)

   allocate(bCMO(nb,nb))

   do i=1, nb
    do j=1, nb
       bCMO(i,j)=CMO(indexlist(i),indexlist(j))
    enddo
   enddo

   nprim      = ls%input%BASIS%REGULAR%&
                &ATOMTYPE(itype)%SHELL(ang+1)%nprim
   norb       = ls%input%BASIS%REGULAR%&
                &ATOMTYPE(itype)%SHELL(ang+1)%norb
   !TEST
   nrow      = ls%input%BASIS%REGULAR%&
                &ATOMTYPE(itype)%SHELL(ang+1)%segment(1)%nrow
   ncol       = ls%input%BASIS%REGULAR%&
                &ATOMTYPE(itype)%SHELL(ang+1)%segment(1)%ncol
   IF((nrow.NE.nprim).OR.(ncol.NE.norb))THEN
      WRITE(6,*)'nprim',nprim
      WRITE(6,*)'norb',norb
      WRITE(6,*)'nrow',nrow
      WRITE(6,*)'ncol',ncol
      CALL LSQUIT('Error in trilevel_convert_ao2gcao: The basis is &
           &either segmented as it should not be or the basis is somehow &
           &corrupted. This is probably because you try to run a &
           &DALTON.INP file with the .RUN WAVE keyword, but use a &
           &LSDALTON_ONLY PRECOMPILER FLAG. You need to have .RUN LINSCA &
           &WHEN YOU USE LSDALTON_ONLY PRECOMPILER FLAG. TK',-1)
   ENDIF
   contrCoeffs=>ls%input%BASIS%REGULAR%&
                &ATOMTYPE(itype)%SHELL(ang+1)%segment(1)%elms
   allocate(CCtmp(nprim,norb))

   CCtmp = reshape(contrCoeffs(1:nprim*norb), (/ nprim, norb /))

   CCtmp = matmul(CCtmp,bCMO)
   
   contrCoeffs = reshape(CCtmp, (/ nprim*norb /))
   
   deallocate(bCMO,CCtmp,indexlist)
 enddo 

 deallocate(basis_size)

end subroutine trilevel_convert_ao2gcao

!!$!> \brief copy ls to minls 
!!$!> \author Branislav Jansik
!!$!> \date 2010-03-03
!!$!> \param ls lsitem structure containing full integral,molecule,basis info
!!$!> \param minls lsitem structure containing valens integral,molecule,basis info
!!$SUBROUTINE trilevel_alloc_sync_ls(minls,ls)
!!$use typedef
!!$implicit none
!!$TYPE(LSITEM) :: minls,ls
!!$  call trilevel_alloc_sync_daltoninput(minls%input,ls%input)
!!$  call II_init_setting(minls%setting)
!!$  call II_set_default_setting(minls%setting,minls%input)
!!$  minls%lupri = ls%lupri
!!$  minls%luerr = ls%luerr
!!$END SUBROUTINE trilevel_alloc_sync_ls
!!$
!!$!> \brief 
!!$!> \author Branislav Jansik
!!$!> \date 2010-03-03
!!$!> \param 
!!$!> \param 
!!$!> \param 
!!$!> \param 
!!$!> \param 
!!$SUBROUTINE trilevel_alloc_sync_daltoninput(NDALTON,DALTON)
!!$use typedef
!!$IMPLICIT NONE
!!$TYPE(daltoninput) :: DALTON
!!$TYPE(daltoninput) :: NDALTON
!!$!TYPE(daltonitem),pointer :: NDALTON(:)
!!$
!!$! STRUCTURE
!!$  NDALTON = DALTON
!!$  !ALL STRUCTURES IN THIS STRUCTURE LIKE, THE DALTONITEM IS COPIED BY THIS
!!$  !LINE, ALTHOUGH ALL STRUCTURES WHICH INCLUDES POINTERS ARE NOT COPIED CORRECTLY
!!$  !WE THEREFORE NEED TO ALLOCATE ALL THE POINTERS AND THEN COPY THE CONTENT OF
!!$  ! THESE POINTERS 
!!$
!!$! THE MOLECULE  (CONTAINS THE POINTER MOLECULEINFO%ATOM)
!!$  NULLIFY(NDALTON%MOLECULE)
!!$  ALLOCATE(NDALTON%MOLECULE)
!!$  NULLIFY(NDALTON%BASIS)
!!$  ALLOCATE(NDALTON%BASIS)
!!$  NDALTON%MOLECULE = DALTON%MOLECULE
!!$  NDALTON%BASIS = DALTON%BASIS
!!$  !AND NOW COPY THE MOLECULEINFO%ATOM
!!$  NULLIFY(NDALTON%MOLECULE%ATOM)
!!$  ALLOCATE(NDALTON%MOLECULE%ATOM(DALTON%MOLECULE%nAtoms))
!!$  NDALTON%MOLECULE%ATOM = DALTON%MOLECULE%ATOM
!!$  
!!$  !THE SAME IS DONE FOR THE BASISSETINFO 
!!$  CALL trilevel_ALLOC_SYNC_BASISSETINFO(NDALTON%BASIS%REGULAR,DALTON%BASIS%REGULAR)
!!$  !CALL trilevel_ALLOC_SYNC_BASISSETINFO(NDALTON%BASIS%HUCKEL,DALTON%BASIS%HUCKEL)
!!$  CALL trilevel_ALLOC_SYNC_BASISSETINFO(NDALTON%BASIS%AUXILIARY,DALTON%BASIS%AUXILIARY)
!!$
!!$  !AND THE BLOCK STRUCTURES - BUT THIS IS ONLY USED IN INTEGRAL EVALUATION
!!$!  NULLIFY(NDALTON%LHSblock%blocks)
!!$!  ALLOCATE(NDALTON%LHSblock%blocks(DALTON%LHSblock%numBlocks))
!!$!  NDALTON%LHSblock%blocks = DALTON%LHSblock%blocks
!!$!  NULLIFY(NDALTON%RHSblock%blocks)
!!$!  ALLOCATE(NDALTON%RHSblock%blocks(DALTON%RHSblock%numBlocks))
!!$!  NDALTON%RHSblock%blocks = DALTON%RHSblock%blocks
!!$
!!$END SUBROUTINE trilevel_ALLOC_SYNC_DALTONINPUT

SUBROUTINE trilevel_alloc_sync_BASISSETINFO(NBASISINFO,BASISINFO)
IMPLICIT NONE
TYPE(BASISSETINFO) :: BASISINFO,NBASISINFO
INTEGER            :: I,J,K,L,nrow,nsize,icharge,maxcharge,ncol

  IF (BASISINFO%natomtypes.eq.0) THEN
     NBASISINFO%natomtypes=0
     return
  ENDIF
  NBASISINFO%natomtypes = BASISINFO%natomtypes
  NULLIFY(NBASISINFO%ATOMTYPE)
  ALLOCATE(NBASISINFO%ATOMTYPE(BASISINFO%natomtypes))
  NBASISINFO%ATOMTYPE = BASISINFO%ATOMTYPE  
  maxcharge = 0
  NBASISINFO%nAtomtypes = BASISINFO%nAtomtypes
  DO J=1,BASISINFO%nAtomtypes
     icharge = BASISINFO%ATOMTYPE(J)%charge
     maxcharge = MAX(maxcharge,icharge)
     !NO need to allocate SHELL
     NBASISINFO%ATOMTYPE(J)%nAngmom = &
          &BASISINFO%ATOMTYPE(J)%nAngmom
     DO K=1,BASISINFO%ATOMTYPE(J)%nAngmom
        !NO need to allocate segments
        NBASISINFO%ATOMTYPE(J)%SHELL(K)%nsegments = &
             &BASISINFO%ATOMTYPE(J)%SHELL(K)%nsegments
        DO L=1,BASISINFO%ATOMTYPE(J)%SHELL(K)%nsegments
           nrow=BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%nrow
           ncol=BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%ncol
           nsize=nrow*ncol
           NULLIFY(NBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms)
           NULLIFY(NBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents)
           ALLOCATE(NBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms(nSIZE))
           ALLOCATE(NBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents(nrow))
           NBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms(:) = &
                & BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%elms(:)
           NBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents(:) = &
                & BASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%Exponents(:)
           NBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%nrow = nrow
           NBASISINFO%ATOMTYPE(J)%SHELL(K)%segment(L)%ncol = ncol 
        ENDDO
     ENDDO
  ENDDO
  NBASISINFO%labelindex = BASISINFO%labelindex
  IF(BASISINFO%labelindex.EQ.0)THEN
     NULLIFY(NBASISINFO%chargeindex)
     ALLOCATE(NBASISINFO%chargeindex(BASISINFO%nchargeindex))
     NBASISINFO%nchargeindex = BASISINFO%nchargeindex
     NBASISINFO%chargeindex(1:BASISINFO%nchargeindex) = &
          &BASISINFO%chargeindex(1:BASISINFO%nchargeindex)
  ENDIF

END SUBROUTINE trilevel_ALLOC_SYNC_BASISSETINFO

!> \brief loop over all atoms and construct the grand canonical basis 
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param opt optItem containing info about scf optimization
!> \param ls lsitem structure containing integral,molecule,basis info
!> \param ai trilevel_atominfo structure containing info about the unique atoms
SUBROUTINE trilevel_gcbasis(opt,ls,ai,LUPRI,LUERR)
use READMOLEFILE
use BUILDBASISSET
use opttype
IMPLICIT NONE
type(optItem),intent(in) :: opt
INTEGER             :: I,LUPRI,LUERR,IPRINT
TYPE(lsitem),intent(in) :: ls
type(trilevel_atominfo) :: ai
TYPE(lsitem),pointer :: atomic_ls
Type(Matrix)         :: F, H1, S , D, CMO
integer              :: nbast, len,iAO,itype
type(moleculeinfo),target :: atomicmolecule
TYPE(lssetting)           :: atomicSetting

do i=1, ai%ND
   call io_free(ls%setting%IO)
   call io_init(ls%setting%IO)
    itype = ai%UATOMTYPE(i)
   !print statements
   len = len_trim(ls%input%BASIS%REGULAR%ATOMTYPE(itype)%NAME)
   write (lupri,'(X,A,X,A,X,A,X,I3)') 'Level 1 atomic calculation on', &
  & ls%input%BASIS%REGULAR%ATOMTYPE(itype)%NAME(1:len), 'Charge', &
  & ls%input%BASIS%REGULAR%ATOMTYPE(itype)%Charge
   write (*,'(X,A,X,A,X,A,X,I3)') 'Level 1 atomic calculation on', &
  & ls%input%BASIS%REGULAR%ATOMTYPE(itype)%NAME(1:len), 'Charge', &
  & ls%input%BASIS%REGULAR%ATOMTYPE(itype)%Charge
 
   nbast = ls%input%molecule%atom(ai%NATOM(i))%nContOrbREG
   CALL mat_init(F,nbast,nbast)
   CALL mat_init(H1,nbast,nbast)
   CALL mat_init(S,nbast,nbast)
   CALL mat_init(D,nbast,nbast)
   CALL mat_init(CMO,nbast,nbast)
 
   !we change the setting to point to the atom of interest 
   CALL build_atomicmolecule(ls%input%molecule,atomicmolecule,ai%NATOM(i),lupri)
   CALL II_init_setting(atomicSetting)
   call II_set_default_setting(atomicSetting,ls%input)
   do iAO=1,4
      atomicSetting%molecule(iAO)%p => atomicmolecule
      atomicSetting%fragment(iAO)%p => atomicmolecule
   enddo
   atomicSetting%numNodes = 1
   atomicSetting%numFragments = 1
   atomicSetting%scheme%fragment = .FALSE.
   atomicSetting%scheme%densfit = .FALSE.
   atomicSetting%scheme%grdone = 0 !the DFT grid should be made for all unique atoms
 
   !build atomic overlap
   CALL II_get_overlap(lupri,luerr,atomicSetting,S)
   !build atomic h1
   CALL II_get_h1(lupri,luerr,atomicSetting,H1)
 
   !build Density matrix eq. 7 from article
   CALL trilevel_get_density_blocked(D,CMO,H1,S,ls%input%basis,itype,nbast)
   !use Density as start guess for a SCF convergence for this atom
   call trilevel_gcscfloop(opt,D,CMO,H1,F,S,ai,atomicSetting,ls%input%molecule,ls%input%basis,&
      &                    i,lupri,luerr)
  !use the MO coefficients from the gcscf loop to create the gcao basis (Eq. 8)
  !WARNING: this changes the input basis in ls to the grand canonical basis
  call trilevel_convert_ao2gcao(ls,CMO%elms,nbast,ai,i)

  CALL II_free_setting(atomicSetting)
  call free_moleculeinfo(atomicmolecule)
  CALL mat_free(F)
  CALL mat_free(H1)
  CALL mat_free(S)
  CALL mat_free(D)
  CALL mat_free(CMO)

enddo

end subroutine trilevel_gcbasis

!> \brief initialise the trilevel_atominfo containing info about the unique atoms
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param ai trilevel_atominfo containing info about the unique atoms to be built
!> \param ls lsitem structure containing the entire molecule and basis
!> \param lupri output logical unit number
subroutine trilevel_atominfo_init(ai,ls,LUPRI)
use READMOLEFILE
use BUILDBASISSET
IMPLICIT NONE
TYPE(trilevel_atominfo) :: ai !Atomic Info
TYPE(lsitem) :: ls
INTEGER             :: LUPRI,IPRINT,I

ai%LUPRI = LUPRI
ai%ND = ls%input%BASIS%REGULAR%nAtomtypes
ai%NA = ls%input%MOLECULE%Natoms
ALLOCATE(ai%NATOM(ai%ND),ai%UATOMTYPE(ai%ND))

IPRINT=ls%input%dalton%BASPRINT
CALL BUILD_DISTINCT_ATOMS(LUPRI,ai%ND,ai%NA,ai%NATOM,ai%UATOMTYPE,ls,IPRINT)

end subroutine trilevel_atominfo_init

!> \brief deallocation routines for freeing the trilevel_atominfo containing info about the unique atoms
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param ai trilevel_atominfo containing info about the unique atoms to be built
subroutine trilevel_atominfo_free(ai)
  use daltoninfo
IMPLICIT NONE
TYPE(trilevel_atominfo) :: ai !Atomic Info
integer                 :: I

DEALLOCATE(ai%NATOM,ai%UATOMTYPE) 

end subroutine trilevel_atominfo_free

!> \brief build the valens basis as a subset of full, ls contains on input the full basis and as output the valens basis
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param ls lsitem structure containing integral,molecule,basis info
!> \param ai trilevel_atominfo structure containing info about the unique atoms
subroutine trilevel_full2valence(VBASIS,ls,ai)
use READMOLEFILE
use BUILDBASISSET
use BUILDAOBATCH
use EcData
IMPLICIT NONE
INTEGER             :: I,LUPRI,IPRINT
TYPE(lsitem) :: ls
TYPE(BASISINFO) :: VBASIS
TYPE(trilevel_atominfo) :: ai
!TYPE(lsitem),pointer :: atomic_ls(:)
integer             :: nbast,itype,nAngmom,charge,ang,nb,nocc,icharge,jatom
integer, pointer    :: basis_size(:)

VBASIS = ls%input%basis
CALL trilevel_ALLOC_SYNC_BASISSETINFO(VBASIS%REGULAR,ls%input%basis%REGULAR)
CALL trilevel_ALLOC_SYNC_BASISSETINFO(VBASIS%AUXILIARY,ls%input%basis%AUXILIARY)

do i=1, ai%ND
   itype = ai%UATOMTYPE(i)
   !points to function which calculates the size of the basis
   basis_size=>trilevel_set_basis_size(ls%setting%BASIS(1)%p,itype)
   jatom = ai%NATOM(i) !an atom in the full input molecule
   nAngmom = ls%input%BASIS%REGULAR%ATOMTYPE(itype)%nAngmom
   charge  = ls%input%BASIS%REGULAR%ATOMTYPE(itype)%Charge
   
   VBASIS%REGULAR%ATOMTYPE(itype)%ToTnorb = &
        & ElementTable(charge)%nocc_s + ElementTable(charge)%nocc_p &
        &+ElementTable(charge)%nocc_d + ElementTable(charge)%nocc_f
   
   if (ElementTable(charge)%nocc_s .gt. 0 ) VBASIS%REGULAR%ATOMTYPE(itype)%nAngmom = 1
   if (ElementTable(charge)%nocc_p .gt. 0 ) VBASIS%REGULAR%ATOMTYPE(itype)%nAngmom = 2
   if (ElementTable(charge)%nocc_d .gt. 0 ) VBASIS%REGULAR%ATOMTYPE(itype)%nAngmom = 3
   if (ElementTable(charge)%nocc_f .gt. 0 ) VBASIS%REGULAR%ATOMTYPE(itype)%nAngmom = 4
   
   do ang=0, nAngmom-1
      
      nb=basis_size(ang+1)
      if (nb.eq.0) cycle
      
      select case(ang)
      case(0); nocc=ElementTable(charge)%nocc_s
      case(1); nocc=ElementTable(charge)%nocc_p/3
      case(2); nocc=ElementTable(charge)%nocc_d/5
      case(3); nocc=ElementTable(charge)%nocc_f/7
      case default; nocc=0;
      end select
      VBASIS%REGULAR%ATOMTYPE(itype)%SHELL(ang+1)%norb = nocc
      VBASIS%REGULAR%ATOMTYPE(itype)%SHELL(ang+1)%segment(1)%ncol = nocc
      if (nocc.eq.0) then
         VBASIS%REGULAR%ATOMTYPE(itype)%SHELL(ang+1)%nprim = 0
         VBASIS%REGULAR%ATOMTYPE(itype)%SHELL(ang+1)%segment(1)%nrow = 0
      endif
   enddo
   deallocate(basis_size)
   nullify(basis_size)
enddo

call determine_nbast(ls%input%MOLECULE,VBASIS%REGULAR,&
      &ls%setting%scheme%DoSpherical,ls%setting%scheme%uncont)

END SUBROUTINE trilevel_full2valence

!> \brief build trilevel atomic density as a diagonal matrix with occupation numbers on the diagonal
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param D density matrix
!> \param ls lsitem structure containing basis info
subroutine trilevel_ATOMS_density(D,ls,REGBASIS)
  use BUILDAOBATCH
implicit none
TYPE(lsitem) :: ls
TYPE(basissetinfo) :: REGbasis
Type(Matrix) :: D
real(realk),allocatable  :: occ(:), tmp(:,:)
integer :: nAtoms, nAngmom, norb, ipos, charge,icharge
integer :: itype, ang, i, nbast, kmult

 nbast = REGBASIS%nbast
 allocate(occ(nbast))

 occ = 0d0;

 nAtoms = ls%input%MOLECULE%nAtoms

 ipos = 1
 do i=1, nAtoms
    IF(REGBASIS%labelindex .EQ.0)THEN
       icharge = INT(ls%input%MOLECULE%ATOM(i)%charge) 
       itype = REGBASIS%chargeindex(icharge)
    ELSE
       itype = ls%input%MOLECULE%ATOM(i)%IDtype(1)
    ENDIF
    nAngmom = REGBASIS%ATOMTYPE(itype)%nAngmom
    charge  = REGBASIS%ATOMTYPE(itype)%Charge

    kmult = 1
    do ang = 0,nAngmom-1
       norb = REGBASIS%ATOMTYPE(itype)%SHELL(ang+1)%norb
       call trilevel_set_occ(occ(ipos:nbast),ang,charge)  
       ipos = ipos + (norb*kmult)
       kmult = kmult + 2
    enddo
 enddo

 allocate(tmp(nbast,nbast))
 tmp = 0d0;

 do i=1,nbast
  tmp(i,i) = occ(i)
 enddo

 call mat_set_from_full(tmp,1d0,D,unres3=.true.)

 deallocate(occ,tmp)
end subroutine trilevel_ATOMS_density

!> \brief create a list to convert valence to full
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param list 
!> \param vlist 
!> \param len
!> \param ls lsitem structure containing full integral,molecule,basis info
!> \param vbasis valence basis
subroutine trilevel_setlist_valence2full(list,vlist,len,ls,vbasis)
use BUILDAOBATCH
implicit none
TYPE(lsitem), intent(in)      :: ls
type(basissetinfo)            :: vbasis
integer, pointer              :: vlist(:,:), list(:,:)
integer, intent(out)          :: len
!
integer :: nAtoms, nAngmom, vnAngmom, norb, ipos, charge
integer :: itype, ang, i,j,k, kmult, icharge
integer :: istart, iend, vistart, viend, vnorb

 nAtoms = ls%input%MOLECULE%nAtoms

 len = 0
 do i=1, nAtoms
    IF(vbasis%labelindex .EQ.0)THEN
       icharge = INT(ls%input%MOLECULE%ATOM(i)%charge) 
       itype = ls%input%BASIS%REGULAR%chargeindex(icharge)
    ELSE
       itype = ls%input%MOLECULE%ATOM(i)%IDtype(1)
    ENDIF

    nAngmom = vbasis%ATOMTYPE(itype)%nAngmom
    
    len = len + nAngmom
 enddo

 allocate(list(len,2),vlist(len,2))

 istart = 1; vistart=1; k=1
 do i=1, nAtoms
    IF(vbasis%labelindex .EQ.0)THEN
       icharge = INT(ls%input%MOLECULE%ATOM(i)%charge) 
       itype = ls%input%BASIS%REGULAR%chargeindex(icharge)
    ELSE
       itype = ls%input%MOLECULE%ATOM(i)%IDtype(1)
    ENDIF

      vnAngmom = vbasis%ATOMTYPE(itype)%nAngmom
       nAngmom = ls%input%BASIS%REGULAR%ATOMTYPE(itype)%nAngmom

      kmult = 1
      do ang = 0,vnAngmom-1
         norb = ls%input%BASIS%REGULAR%ATOMTYPE(itype)%SHELL(ang+1)%norb
        vnorb = vbasis%ATOMTYPE(itype)%SHELL(ang+1)%norb

         iend =  istart  + (vnorb*kmult) -1
        viend = vistart  + (vnorb*kmult) -1

         list(k,1)=  istart;   list(k,2)=  iend
        vlist(k,1)= vistart;  vlist(k,2)= viend

         k = k + 1
         istart =  istart + (norb*kmult)
        vistart = viend + 1
         kmult = kmult + 2
      enddo

      do ang=vnAngmom, nAngmom-1
         norb =  ls%input%BASIS%REGULAR%ATOMTYPE(itype)%SHELL(ang+1)%norb
         istart = istart + (norb*kmult)
         kmult = kmult +2
      enddo

 enddo

end subroutine trilevel_setlist_valence2full

!> \brief convert the Density matrix from the valence basis to full basis
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param D full density matrix
!> \param Dval valence density matrix
!> \param list
!> \param vlist
!> \param len
subroutine trilevel_density_valence2full(D,Dval,list,vlist,len)
implicit none
Type(Matrix) :: D,Dval
real(realk),allocatable  :: tD(:,:), tDval(:,:)
integer :: i,j,k,len
integer, pointer :: vlist(:,:), list(:,:)

 allocate(tD(D%ncol,D%nrow),tDval(Dval%ncol,Dval%nrow))

 call mat_to_full(Dval,1d0,tDval)
 tD = 0d0

 do i=1,len
 if (vlist(i,1).gt.vlist(i,2)) cycle
  do j=1,len
    if (vlist(j,1).gt.vlist(j,2)) cycle

    tD(list(i,1):list(i,2),list(j,1):list(j,2)) = &
   & tDval(vlist(i,1):vlist(i,2),vlist(j,1):vlist(j,2))

  enddo
 enddo

 call mat_set_from_full(tD,1d0,D,'Dmat')

 deallocate(tD,tDval)
end subroutine trilevel_density_valence2full

!> \brief convert the MO coefficients from the valence basis to full basis
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param CMO MO coefficients in full basis
!> \param CMO MO coefficients in valence basis
!> \param list
!> \param vlist
!> \param len
subroutine trilevel_cmo_valence2full(CMO,vCMO,list,vlist,len)
implicit none
Type(Matrix) :: CMO,vCMO
real(realk),allocatable  :: tCmo(:,:), tvCMO(:,:)
integer :: i,j,k,  len
integer, pointer :: vlist(:,:), list(:,:)

 allocate(tCMO(CMO%ncol,CMO%nrow),tvCMO(vCMO%ncol,vCMO%nrow))

 call mat_to_full(vCMO,1d0,tvCMO)
 tCMO = 0d0

 do i=1,len
 if (vlist(i,1).gt.vlist(i,2)) cycle

    tCMO(list(i,1):list(i,2),1:vCMO%ncol) = &
   & tvCMO(vlist(i,1):vlist(i,2),1:vCMO%ncol)
 enddo

 k = vCMO%ncol + 1;
 do i=1, len-1
 if (vlist(i,1).gt.vlist(i,2)) cycle
   do j=list(i,2)+1 , list(i+1,1)-1
      tCMO(j,k) = 1d0
      k = k + 1
   enddo  
 enddo

 do j=list(len,2)+1, CMO%ncol
      tCMO(j,k) = 1d0
      k = k + 1
 enddo

 call mat_set_from_full(tCMO,1d0,CMO)

 deallocate(tCMO,tvCMO)
end subroutine trilevel_cmo_valence2full

!> \brief read density from disk
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param 
!> \param 
!> \param 
!> \param 
!> \param 
subroutine trilevel_readdens(D,lupri)
use files
implicit none
Type(Matrix) :: D
integer      :: lu, lupri, IDUM,LDUM
data lu /-1/

call lsopen(lu,'vdens.restart','OLD','UNFORMATTED')
rewind lu
call mat_read_from_disk(lu,D)
call lsclose(lu,'KEEP')

WRITE(LUPRI,*)
WRITE(LUPRI,*) '*** RESTART FROM DENSITY ON DISK - READ FROM vdens.restart  ***'
WRITE(LUPRI,*)
WRITE(*,*)
WRITE(*,*) '*** RESTART FROM DENSITY ON DISK - READ FROM vdens.restart  ***'
WRITE(*,*)

end subroutine trilevel_readdens

end module trilevel_module

!> \brief main wrapper routine to obtain the trilevel_basis
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param opt optItem containing info about scf optimization
!> \param ls lsitem structure containing integral,molecule,basis info
SUBROUTINE trilevel_basis(opt,ls)
use trilevel_module
use typedef
use ecdata
use opttype
implicit none
type(optItem), intent(inout) :: opt
type(optItem)             :: gcopt
TYPE(lsitem) :: ls
TYPE(trilevel_atominfo) :: ai
integer                 :: CFG_averaging_sav, matrix_sav
logical                 :: unres_sav

  opt%optlevel = 1

  matrix_sav  = matrix_type
  matrix_type = mtype_dense

  gcopt = opt
  gcopt%cfg_unres = .false.
 
  !initialise element tables : list of s,p,d,f orbital occupation
  call ecdata_init

  !initialise atominfo : find # unique atoms, and map to atoms in molecule
  call trilevel_atominfo_init(ai,ls,gcopt%LUPRI)

  !main driver to build the grand canonical basis
  call trilevel_gcbasis(gcopt,ls,ai,gcopt%LUPRI,gcopt%LUERR)

  !free atominfo
  call trilevel_atominfo_free(ai)
  matrix_type = matrix_sav

  !CFG_averaging = CFG_averaging_sav
  !we have written screening matrices which should not 
  !be used in the rest of the calculations so we reset the IO  
  call io_free(ls%setting%IO)
  call io_init(ls%setting%IO)
END SUBROUTINE trilevel_basis

!> \brief the atoms start guess routine. The gcao basis have been created by a call to trilevel_basis and now a start D is made
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param decomp Contains matrices from OAO decomposition of overlap matrix
!> \param D density matrix 
!> \param H1 one electron contribution to the fock matrix
!> \param ls lsitem structure containing all info about integrals,basis,..
SUBROUTINE atoms_start(config,D,H1,ls)
use configuration
use trilevel_module
use typedef
use matrix_util
use lsdalton_fock_module
use dal_interface
use initial_guess
implicit none
type(ConfigItem),intent(in) :: config
TYPE(lsitem),intent(inout) :: ls
Type(Matrix),target        :: H1
Type(Matrix),intent(inout) :: D
!
Type(Matrix) :: Dval,F,Cmo
TYPE(trilevel_atominfo) :: ai
real(realk)         :: E
real(realk),allocatable :: eival(:) 
integer      :: nbast
logical      :: dalink

  !a diagonal matrix with occupation numbers on the diagonal
  call trilevel_ATOMS_density(D,ls,ls%input%basis%regular)

  nbast = ls%input%BASIS%REGULAR%nbast

  CALL mat_init(F,nbast,nbast)
  CALL mat_init(Cmo,nbast,nbast)

  !We cannot use DaLink in the 0'th iteration - this gives a diagonal Fock matrix
  ! => bad starting guess. If DaLink is requested, turn it off and then back on after
  ! the 0'th iteration. /Stinne, Thomas, Brano 19/11-2009
  dalink = .false.
  if (ls%setting%scheme%DALINK) then
     ls%setting%scheme%DALINK = .FALSE.
     dalink = .true.
  endif
  ls%setting%scheme%DFTELS = 1.d0 !the density is not idempotent so it would give a wrong number of electrons
  ! Iteration 0 : The density matrix is not idempotent; a diagonalization gives a proper 
  ! idempotent density matrix  
  call di_get_fock_LSDALTON(D,H1,F,E,config%decomp%lupri,config%decomp%luerr,ls)
  write(*,*) ' Iteration 0 energy:', E
  write(config%decomp%lupri,*) ' Iteration 0 energy:', E

  if (config%decomp%cfg_unres) then; allocate(eival(2*nbast)); else
    allocate(eival(nbast))
  endif
  call mat_diag_f(F,config%decomp%S,eival,Cmo)

  !Asymetrizing starting guess if .ASYM is in input
  ! 21.04.2010 C. Nygaard
  !Only works if HOMO and LUMO are of different symmetry
  if (config%decomp%cfg_unres .and. config%opt%cfg_asym) then
    call asymmetrize_starting_guess (Cmo, config%decomp)
  endif

  call density_from_orbs(config%decomp%cfg_unres,config%decomp%nocc,config%decomp%nocca,config%decomp%noccb,Cmo,D)
  deallocate(eival)

  !Turn DaLink back on, if requested:
  if (dalink) then
     ls%setting%scheme%DALINK = .true.
  endif 
  ls%setting%scheme%DFTELS = ls%input%dalton%DFTELS 

  CALL mat_free(F)
  CALL mat_free(Cmo)

END SUBROUTINE ATOMS_START

!> \brief Main Driver routine for the trilevel start guess
!> \author Branislav Jansik
!> \date 2010-03-03
!> \param D density matrix
!> \param ls structure containing info about molecule and integral evaluation
!> \param config structure containing basicly all info
SUBROUTINE trilevel_start(D,ls,config)
use configuration
use files
use trilevel_module
use typedef
use ecdata
use lsdalton_fock_module
use dal_interface
use ks_settings
use direct_dens_util
use daltoninfo
implicit none
TYPE(lsitem),target :: ls
TYPE(BASISINFO),target :: vbasis !valens basis 
TYPE(trilevel_atominfo) :: ai
Type(Matrix) :: Dval,F,D,Cmo, fCmo
Type(Matrix),target :: H1,S
Type(lsint_fock_data_type) :: lsint_fock_data_sav
real(realk)         :: E, mx
real(realk)         :: maxelm_save, maxstep_save,thr_save
real(realk),allocatable :: eival(:) 
integer      :: nbast, len, nocc
integer, pointer :: vlist(:,:), list(:,:)
integer :: densL2_lun, idum, ldum,iAO
integer(8) :: fperm, vperm
logical :: restart_from_dens, no_rhdiis, dalink, vdens_exists
type(ConfigItem) :: config

  !initialise the IO so that old screening matrices are not used 
  call io_free(ls%setting%IO)
  call io_init(ls%setting%IO)

  !build a minls used in the valens calculation
!  call trilevel_alloc_sync_ls(minls,ls)

  !build a trilevel_atominfo structure containing info for each distinct atom
  call trilevel_atominfo_init(ai,ls,config%LUPRI)

  !loops over all distinct atoms and build valensbasis
  call trilevel_full2valence(vbasis,ls,ai)
  !set setting basis to the valens basis

  do iAO=1,4
     ls%setting%basis(iAO)%p => vbasis
  enddo

#ifdef HAVE_BSM
 if (config%opt%CFG_prefer_BSM) then
    call bsm_lib_getpermutation(fperm)
    call ls_initbsm(VBASIS%REGULAR,ls%input)
    call bsm_lib_getpermutation(vperm)
 endif
#endif

  write(config%lupri,'(1X,A)') 'Level 2 molecular calculation'
  write(   * ,'(1X,A)') 'Level 2 molecular calculation'

  nbast = VBASIS%REGULAR%nbast !nbast for valence basis
  call mat_init(Dval,nbast,nbast)

  !select starting from atoms or from file
  INQUIRE(file='vdens.restart',EXIST=vdens_exists) 
  restart_from_dens = config%diag%cfg_restart.AND.vdens_exists
  if (restart_from_dens) then
     !.RESTART specified in input and if a density 
     !have been written to file
     call trilevel_readdens(Dval,config%lupri)
  else
     !default option
     call trilevel_ATOMS_density(Dval,ls,VBASIS%REGULAR)
  endif

  CALL mat_init(F,nbast,nbast)
  CALL mat_init(H1,nbast,nbast)
  CALL mat_init(S,nbast,nbast)
  CALL mat_init(Cmo,nbast,nbast)

  !Get overlap and H1
  CALL II_get_overlap(config%decomp%lupri,config%decomp%luerr,ls%setting,S)
  CALL II_get_h1(config%decomp%lupri,config%decomp%luerr,ls%setting,H1)

  config%decomp%S => S

  !save L3 inputs
  lsint_fock_data_sav = lsint_fock_data
  lsint_fock_data%ls  => ls
  lsint_fock_data%H1  => H1
  lsint_fock_data%lupri = config%lupri
  lsint_fock_data%luerr = config%decomp%luerr
  !maxstep_save = cfg_max_step    
  !maxelm_save = cfg_max_element
  !thr_save = cfg_convergence_threshold !We don't need to converge hard on level 2
  config%opt%set_convergence_threshold = config%opt%cfg_convergence_threshold*config%opt%cfg_level2_convfactor

  !We cannot use DaLink in the 0'th iteration - this gives a diagonal Fock matrix
  ! => bad starting guess. If DaLink is requested, turn it off and then back on after
  ! the 0'th iteration. /Stinne, Thomas, Brano 19/11-2009
  dalink = .false.
  if (ls%setting%scheme%DALINK) then
     ls%setting%scheme%DALINK = .FALSE.
     dalink = .true.
  endif
  ! Iteration 0
  ls%setting%scheme%DFTELS = 1.d0 !the density is not idempotent so it would give a wrong number of electrons
  if (.not.restart_from_dens) then
   call di_get_fock_LSDALTON(Dval,H1,F,E,config%lupri,config%decomp%luerr,ls)
   write(*,*) ' Iteration 0 energy:', E
   write(config%lupri,*) ' Iteration 0 energy:', E
   allocate(eival(nbast))
   call mat_diag_f(F,S,eival,Cmo)
   deallocate(eival)
   call density_from_orbs(config%decomp%cfg_unres,config%decomp%nocc,config%decomp%nocca,config%decomp%noccb,Cmo,Dval)
  endif
  ls%setting%scheme%DFTELS = ls%input%dalton%DFTELS
  !Turn DaLink back on, if requested:
  if (dalink) then
     ls%setting%scheme%DALINK = .true.
  endif

  !initialize incremental scheme
  if (config%opt%cfg_incremental) call ks_init_incremental_fock(nbast)

  !Does optimization method require overlap decomposition?
  no_rhdiis =(config%opt%cfg_density_method == config%opt%cfg_f2d_direct_dens .or. &
            & config%opt%cfg_density_method == config%opt%cfg_f2d_arh .or. &
            & config%decomp%cfg_check_converged_solution .or. &
            & config%decomp%cfg_rsp_nexcit > 0 ) 

  !Do decomposition
  if (no_rhdiis .or. config%decomp%cfg_lcv) then   
     call decomp_init(nbast,config%decomp)
     call dd_mat_eigenvalues_to_aux(config%decomp%cfg_unres,S)
     config%decomp%lcv_basis = .false.
     call decomposition(config%decomp)
  endif

  !Prepare for L2 calculation in localized CMO basis.
  !Do localization of initial orbitals and do decomposition
  !to get local orthonormal CMO basis
  if (config%decomp%cfg_lcv.and.(.not.restart_from_dens).and. no_rhdiis) then
     nocc = config%decomp%nocc
     call leastchange_lcv(config%decomp,Cmo,nocc,ls)

     call leastchangeOrbspreadStandalone(mx,ls,Cmo,config%decomp%lupri,config%decomp%luerr)
     write(*,*) 'Orbspread standalone: ', mx

     config%decomp%lcv_basis = .true.
     call mat_init(config%decomp%lcv_CMO,nbast,nbast)
     call mat_assign(config%decomp%lcv_CMO,Cmo)

     call save_decomposition(config%decomp)
     call decomposition(config%decomp)
  endif

  if (config%av%CFG_averaging == config%av%CFG_AVG_van_lenthe) then !FIXME: put this somewhere else!
     call mat_init(config%av%Fprev,nbast,nbast)
     call mat_init(config%av%Dprev,nbast,nbast)
  endif

  ! L2 minimization
  config%opt%optlevel = 2
  call scfloop(H1,F,Dval,S,E,ls,config)
  write(   * ,'(1X,A)') 'L2 minimization done!'
  write(config%lupri,'(1X,A)') 'L2 minimization done!'
  config%opt%optlevel = 3
 
   if (config%av%CFG_averaging == config%av%CFG_AVG_van_lenthe) then
      call mat_free(config%av%Fprev)
      call mat_free(config%av%Dprev)
   endif

  !Dump level 2 density matrix to disk
  densL2_lun = -1
  call lsopen(densL2_lun,'densL2.restart','NEW','UNFORMATTED')
  call mat_write_to_disk(densL2_lun,Dval)
  call lsclose(densL2_lun,'KEEP')

  !release data related to incremental fock
  if (config%opt%cfg_incremental) call ks_free_incremental_fock
  !restore L3 inputs
  lsint_fock_data = lsint_fock_data_sav
  config%solver%set_max_step    = config%solver%cfg_max_step
  config%solver%set_max_element = config%solver%cfg_max_element
  config%solver%set_local       = .false.
  config%opt%set_convergence_threshold = config%opt%cfg_convergence_threshold
  if (config%av%cfg_averaging == config%av%cfg_avg_van_lenthe) then
     config%av%vanlentheCounter = 0
     config%diag%CFG_lshift = config%diag%cfg_lshift_vanlenthe
     config%av%usequeue = .false.
   endif
  !cfg_max_step = maxstep_save
  !cfg_max_element = maxelm_save
  !cfg_dd_local = .false.
  !cfg_convergence_threshold = thr_save

  ! move dens.restart file to vdens.restart
#ifdef SYS_AIX
  ! ToDo: rename does not work. Combination aix/xlf90_r on huge. /SR 2010-03-30
  write(config%LUPRI,'(1X,A)') 'Warning: Renaming of dens.restart to vdens.restart not working for AIX/XLF90'
#else
  call rename('dens.restart','vdens.restart')
#endif

  ! get Valence Cmo
  allocate(eival(nbast))
  call mat_diag_f(F,S,eival,Cmo)
  deallocate(eival)


  CALL mat_free(H1)
  CALL mat_free(F)
  CALL mat_free(S)
  !restore overlap decomposition
  if (config%decomp%cfg_lcv.and.(.not.restart_from_dens).and. no_rhdiis) &
  & call restore_decomposition(config%decomp)

  ! localization of Valence CMO, by least change algorithm
  if (config%decomp%cfg_lcv) then
      nocc = config%decomp%nocc
      call leastchange_lcv(config%decomp,Cmo,nocc,ls)

      call leastchangeOrbspreadStandalone(mx,ls,Cmo,config%decomp%lupri,config%decomp%luerr)
      write(*,*) 'Orbspread standalone: ', mx
  endif
  ! shutdown of decomposition and dd
  if (no_rhdiis .or. config%decomp%cfg_lcv) then
     call decomp_shutdown(config%decomp)
     !call dd_shutdown(config%decomp%cfg_unres)
  endif

  ! create integer list for conversion to full basis
  call trilevel_setlist_valence2full(list,vlist,len,ls,vbasis%regular)

  ! free mindalton structure
!  call dalton_free(minls%input)
  ! initialize decomp%lcv_CMO, if BSM, we need to temporarily switch
  ! from valence BSM permutation to full BSM permutation 
#ifdef HAVE_BSM
 if (config%opt%CFG_prefer_BSM) call bsm_lib_setpermutation(fperm)
#endif

  nbast = ls%input%BASIS%REGULAR%nbast !full basis nbast 
  call mat_init(config%decomp%lcv_CMO,nbast,nbast)

#ifdef HAVE_BSM
 if (config%opt%CFG_prefer_BSM) call bsm_lib_setpermutation(vperm)
#endif
  ! convert valence CMO to full cmo
  call trilevel_cmo_valence2full(config%decomp%lcv_Cmo,Cmo,list,vlist,len)
  config%decomp%lcv_basis = .true.

  call mat_free(Cmo)
  ! convert valence Dval to full D
  call trilevel_density_valence2full(D,Dval,list,vlist,len)
  CALL mat_free(Dval)
  deallocate(list,vlist)
  !we restore the default settings, which sets the basis to full basis
  !and set the GRDONE=0 so that the dft grid is calculated with the full 
  !basis on the next kohn-sham matrix build
  call II_set_default_setting(ls%setting,ls%input)
  call determine_nbast(ls%input%MOLECULE,ls%input%BASIS%REGULAR,&
      &ls%setting%scheme%DoSpherical,ls%setting%scheme%uncont)
  call io_free(ls%setting%IO)
  call io_init(ls%setting%IO)! we change basis, can nolonger use the mat on disk

  call leastchangeOrbspreadStandalone(mx,ls,config%decomp%lcv_Cmo,config%decomp%lupri,config%decomp%luerr)
  write(*,*) 'Orbspread standalone full CMO: ', mx

  ! release valence BSM permutation and set the full
  ! instead
#ifdef HAVE_BSM
 if (config%opt%CFG_prefer_BSM) then
  call bsm_lib_free(vperm)
  call bsm_lib_setpermutation(fperm)
 endif
#endif
 call trilevel_atominfo_free(ai)
 call free_basissetinfo(vbasis%regular)
 call free_basissetinfo(vbasis%auxiliary)

END SUBROUTINE trilevel_start

