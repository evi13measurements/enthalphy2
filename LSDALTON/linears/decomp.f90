!> @file
!> Contains module for decomposition of overlap matrix.

!> \brief Decomposition of overlap matrix and routines for transforming to OAO basis.
!> \author S. Host
!> \date 2007
module decompMod
use matrix_module

!> \brief Contains settings for decomposition and OAO decomposed overlap.
!> \author S. Host
!> \date 2007
type DecompItem
      !> Logical unit number for LSDALTON.OUT
      integer     :: lupri
      !> Logical unit number for LSDALTON.ERR
      integer     :: luerr
   !CONFIGURATIONS:
   !===============
      !> Do Cholesky decomposition
      logical     :: cholesky_decomp
      !> Do Lowdin decomposition by diagonalization
      logical     :: lowdin_diagonalize
      !> Do iterative Lowdin decomposition
      logical     :: lowdin_iterative
      !> Do iterative Lowdin decomposition in quadruple precision
      logical     :: lowdin_qiterative
      !> Construct Fock/KS matrix directly in OAO without transforming to AO basis
      logical     :: cfg_do_in_oao
      !> Use Least-Change Valence basis 
      logical     :: cfg_lcv
      !> Use Least-Change Molecular basis
      logical     :: cfg_lcm
      !> Use Least-Change Valence basis, brute force!
      logical     :: lcv_basis
      !> Use Grand-Canonical basis. Default for ATOMS and TRILEVEL
      logical     :: cfg_gcbasis
      !> True if unrestricted calculation
      logical     :: cfg_unres
      !> Number of occupied orbitals (if restricted)
      integer     :: nocc
      !> Has number of alpha electrons been specified explicitly in input?
      logical     :: alpha_specified
      !> Number of occupied alpha orbitals (if unrestricted)
      integer     :: nocca
      !> Has number of beta electrons been specified explicitly in input?
      logical     :: beta_specified
      !> Number of occupied beta orbitals (if unrestricted)
      integer     :: noccb
      !> I *think* this means number of unpaired electrons - why is it called spin????
      integer     :: spin 
      !> Number of active electrons
      integer     :: nactive

      !> Maximum number of iterations for iterative HOMO-LUMO gap
      integer     :: cfg_homolumo_maxit
      !> If HL gap is not required for Hessian eival/response calcs, we can just print a warning if not converged
      logical     :: cfg_hlgap_needed
      !> If HL gap is not required and did not converge, don't try to print it
      logical     :: cfg_hlgap_converged
      !> HL gap converged in this number of iterations
      integer     :: cfg_hlgap_nit_conv
      !> Maximum number of iterations for iterative Hessian eigenvalue
      integer     :: cfg_check_maxit
      !> True if lowest Hessian eigenvalue should be calculated
      logical     :: cfg_check_converged_solution
      !> Can be set > 1 if more than one lowest Hessian eigenvalues are wanted
      integer     :: cfg_hessian_nvec

      !> Number of excitation energies requested
      integer     :: cfg_rsp_nexcit
      !> Should MOs be used to construct starting guess for excitation energies?
      logical     :: cfg_rsp_mostart
      !> Do we want number of start vectors for exc. energies to be different (i.e. larger) than no of excitation energies? 
      logical     :: cfg_startvectors
      !> Number of startvectors, if cfg_startvectors = .true.
      integer     :: cfg_no_of_startvectors

      !> Do not use preconditioning when solving for Hessian eigenvalue
      logical     :: cfg_noprec
      !> Relative convergence threshold used when doing PCG preconditioning
      real(realk) :: cfg_micro_thresh  
   !SETTINGS:
   !=========
      !> Set to true when we use PCG for preconditioning (i.e. do not add 2nd order contribution to linear transformation)
      logical     :: pcg_preconditioning
   !Matrices:
   !=========
      !> Least-Change Valence molecular orbitals
      type(Matrix) :: lcv_CMO
      !> Decomposed overlap matrix, S = U*Ut
      type(Matrix) :: U
      !> Inverse of decomposed overlap matrix, S = U*Ut
      type(Matrix) :: U_inv
      !> Used if decomposed overlap matrix should be save
      type(Matrix) :: U_sav
      !> Used if decomposed overlap matrix should be save
      type(Matrix) :: U_inv_sav
      !> Density matrix in OAO basis DU = U*D*Ut
      type(Matrix) :: DU
      !> Fock/KS matrix in OAO basis FU = U_inv_t*F*U_inv
      type(Matrix) :: FU
      !> Fock/KS matrix in OAO basis projected on occupied space FUP = DU*FU*DU
      type(Matrix) :: FUP
      !> Projector onto virtual space QU = I - DU
      type(Matrix) :: QU
      !> Fock/KS matrix in OAO basis projected on virtual space FUQ = QU*FU*QU
      type(Matrix) :: FUQ
      !> Pointer to overlap matrix
      type(matrix),pointer :: S
   !DEBUG VARIABLES:
   !================
      logical      :: debug_rsp_linsca
   !INFO VARIABLES:
   !===============
      !> Print info from orbital energy / Hessian eigenvalue solvers
      logical      :: info_stability
      !> Print reduced space info from orbital energy / Hessian eigenvalue solvers
      logical      :: info_stability_redspace
      !> Print info from normalization routines, orbital energy / Hessian eigenvalue solvers
      logical      :: info_dd_normalize
      logical      :: info_rsp_precond
end type DecompItem

!radovan: this module is used in the response code
!         which is used also in DIRAC
!         this module depends on many things that DIRAC
!         does not have - i have to skip all the following
!         to allow compilation in DIRAC
#ifndef PRG_DIRAC
contains

   !> \brief Default setting for OAO decomposition.
   !> \author S. Host
   !> \date March 2010
   subroutine decomp_set_default_config(decomp)
   implicit none
      !> Contains settings for decomposition and OAO decomposed overlap
      type(decompItem) :: decomp
   
      !CONFIGURATIONS:
      !===============
         decomp%cholesky_decomp      = .false.
         decomp%lowdin_diagonalize   = .true.
         decomp%lowdin_iterative     = .false.
         decomp%lowdin_qiterative    = .false.
   
         decomp%cfg_do_in_oao        = .false.
         decomp%cfg_lcv              = .false.
         decomp%cfg_lcm              = .false.
         decomp%lcv_basis            = .false.
         decomp%cfg_gcbasis          = .false.
   
         decomp%cfg_unres            = .false.
         decomp%alpha_specified      = .false.
         decomp%beta_specified       = .false.
         decomp%spin                 = 0
   
         decomp%cfg_homolumo_maxit   = 5000
         decomp%cfg_hlgap_needed     = .false.
         decomp%cfg_hlgap_converged  = .false.
         decomp%cfg_check_maxit      = 40
         decomp%cfg_check_converged_solution = .false.
         decomp%cfg_hessian_nvec     = 1
   
         decomp%cfg_rsp_nexcit       = 0
         decomp%cfg_rsp_mostart      = .false.
         decomp%cfg_startvectors     = .false.
   
         decomp%cfg_noprec           = .false.
   
         decomp%cfg_micro_thresh     = 1.0d-2
      !SETTINGS:
      !=========
         decomp%pcg_preconditioning = .false. 
      !DEBUG VARIABLES:
      !================
         decomp%debug_rsp_linsca     = .false.
      !INFO VARIABLES:
      !===============
         decomp%info_stability         = .false.
         decomp%info_stability_redspace = .false.
         decomp%info_dd_normalize    = .false.
         decomp%info_rsp_precond     = .false.
      !DATA:
      !=====
         nullify(decomp%S)
   end subroutine decomp_set_default_config

   !> \brief Save decomposed overlap in U_sav, U_inv_sav
   !> \author B. Jansik
   !> \date 2008
   subroutine save_decomposition(decomp)
   use matrix_operations
   implicit none
       !> Contains settings for decomposition and OAO decomposed overlap
        type(DecompItem)    :: decomp
    
        call mat_init(decomp%U_sav,decomp%U%nrow,decomp%U%ncol)
        call mat_init(decomp%U_inv_sav,decomp%U_inv%nrow,decomp%U_inv%ncol)
                       
        call mat_assign(decomp%U_sav,decomp%U) 
        call mat_assign(decomp%U_inv_sav,decomp%U_inv)
     
   end subroutine save_decomposition
   
   !> \brief Restore decomposed overlap from U_sav, U_inv_sav
   !> \author B. Jansik
   !> \date 2008
   subroutine restore_decomposition(decomp)
   use matrix_operations
   implicit none
       !> Contains settings for decomposition and OAO decomposed overlap
        type(DecompItem)    :: decomp
   
        call mat_assign(decomp%U,decomp%U_sav) 
        call mat_assign(decomp%U_inv, decomp%U_inv_sav)
   
        call mat_free(decomp%U_sav)
        call mat_free(decomp%U_inv_sav)
     
   end subroutine restore_decomposition

   !> \brief Initialize OAO matrices.
   !> \author S. Host
   !> \date March 2010
   subroutine decomp_init(ndim,decomp)
   use matrix_operations
     implicit none
     !> Number of basis functions
     integer, intent(in) :: ndim
     !> Contains settings for decomposition and OAO decomposed overlap
     type(DecompItem)    :: decomp

     call mat_init(decomp%U,ndim,ndim)
     call mat_init(decomp%U_inv,ndim,ndim)
     call mat_init(decomp%DU,ndim,ndim)
     call mat_init(decomp%FU,ndim,ndim)
     call mat_init(decomp%FUP,ndim,ndim)
     call mat_init(decomp%FUQ,ndim,ndim)
     call mat_init(decomp%QU,ndim,ndim)

   end subroutine decomp_init

   !> \brief Free OAO matrices.
   !> \author S. Host
   !> \date March 2010
   subroutine decomp_shutdown(decomp)
   use matrix_operations
   implicit none
     !> Contains settings for decomposition and OAO decomposed overlap
     type(DecompItem)    :: decomp

     call mat_free(decomp%U)
     call mat_free(decomp%U_inv)
     call mat_free(decomp%DU)
     call mat_free(decomp%FU)
     call mat_free(decomp%FUP)
     call mat_free(decomp%FUQ)
     call mat_free(decomp%QU)
     nullify(decomp%S)
   end subroutine decomp_shutdown

   !> \brief Orthonormal decomposition of overlap matrix.
   !> \author S. Host, B. jansik
   !> \date 2005
   !>
   !> Get cholesky decomposition of S = Ut*U \n
   !>    or Lowdin decomposition of S = S^1/2 * S^1/2 \n
   !> For Lowdin, options are
   !> - Lowdin by diagonalization
   !> - Lowdin by iterative scheme
   !> - Lowdin by iterative scheme in quadruple precision
   !>
   subroutine decomposition(decomp)
   use matrix_operations
   use lowdin_module
   use lstiming
   implicit none
         !> Contains settings for decomposition and OAO decomposed overlap
         type(DecompItem)                          :: decomp
         type(Matrix)                              :: S_sav
         real(realk),dimension(:), allocatable     :: work1
         real(realk),dimension(:,:), allocatable   :: U_full, U_inv_full 
         integer, dimension(:),allocatable         :: IPVT
         real(realk)                               :: RCOND, dummy(2), tmstart, tmend
         integer                                   :: IERR, i, j, fulldim, ludecomp, ndim

   call lstimer('START ',tmstart,tmend,decomp%lupri)

   ndim = decomp%S%nrow
   if (decomp%cfg_unres) then 
      fulldim = 2*ndim
   else
      fulldim = ndim 
   endif

   if (decomp%cfg_unres) then
      if (.not. decomp%cholesky_decomp) then
         write(decomp%lupri,*) 'FALLBACK: open shell currently works only with cholesky decomposition'
         write(decomp%lupri,*) '- Ask Dr. Jansik to fix it for Lowdin!'
         decomp%cholesky_decomp = .true.
      endif
   endif

   if (decomp%lcv_basis) then
      !compute overlap S in lcv basis, save S in zeta_A basis
      call mat_init(S_sav,ndim,ndim)
      call mat_assign(S_sav,decomp%S)
      call mat_mul(decomp%lcv_CMO,decomp%S,'t','n',1d0,0d0,decomp%U)
      call mat_mul(decomp%U,decomp%lcv_CMO,'n','n',1d0,0d0,decomp%S)
   endif

   if (decomp%cholesky_decomp .or. ndim < 3) then
      !print *, "Preparing to do cholesky decomposition..."
      write (decomp%lupri,*) "Preparing to do cholesky decomposition..."
      allocate(U_full(fulldim,fulldim)) 
      allocate(U_inv_full(fulldim,fulldim))
      allocate(work1(fulldim))
      allocate(IPVT(fulldim))

        !write(lupri,*) 'Overlap S:'
        !call mat_print(decomp%S,1,decomp%S%nrow,1,decomp%S%nrow,lupri)

      U_full = 0.0d0
      call mat_to_full(decomp%S,1.0d0,U_full)

        !write (lupri,*) 'Overlap S, full:'
        !call OUTPUT(U_full, 1, fulldim, 1, fulldim, fulldim, fulldim, 1, lupri)

      !Set lower half of U = 0:
      do i = 1, fulldim
         do j = 1, i-1
            U_full(i,j) = 0.0d0
         enddo
      enddo

      call dchdc(U_full,fulldim,fulldim,work1,0,0,IERR)

      !Invert U and Ut:
      IPVT = 0 ; RCOND = 0.0d0  

      U_inv_full = U_full
      call DGECO(U_inv_full,fulldim,fulldim,IPVT,RCOND,work1)
      call DGEDI(U_inv_full,fulldim,fulldim,IPVT,dummy,work1,01)

      IPVT = 0 ; RCOND = 0.0d0  

      !Convert framework:
      call mat_set_from_full(U_full,1.0d0,decomp%U)
      call mat_set_from_full(U_inv_full,1.0d0,decomp%U_inv)
      
      deallocate(U_full) 
      deallocate(U_inv_full) 
      deallocate(work1)
      deallocate(IPVT)
      call lstimer('CHOLDC',tmstart,tmend,decomp%lupri)

   else if (decomp%lowdin_diagonalize) then
      !print *, "Preparing to do S^1/2 decomposition..."
      write (decomp%lupri,*) "Preparing to do S^1/2 decomposition..."

      if(decomp%cfg_unres) call lsquit('lowdin it not tested for unres',decomp%lupri)

      call lowdin_diagonalize(decomp%lupri,decomp%cfg_unres,decomp%S,decomp%U,decomp%U_inv)
      call lstimer('LWDIAG',tmstart,tmend,decomp%lupri)
   else if (decomp%lowdin_iterative) then
      !print *, "Preparing to do iterative S^1/2 decomposition..."
      write (decomp%lupri,*) "Preparing to do iterative S^1/2 decomposition..."
      call lowdin_schulz(decomp%S,decomp%U,decomp%U_inv,decomp%lupri)

      call lstimer('LWITER',tmstart,tmend,decomp%lupri)

   else if (decomp%lowdin_qiterative) then
      !print *, "Preparing to do iterative S^1/2 decomposition..."
      write (decomp%lupri,*) "Preparing to do iterative S^1/2 decomposition in quadruple precision..."

! Quadruple precision not supported by gfortran-4.1.1
#ifndef GFORTRAN
      call lowdin_qschulz(decomp%cfg_unres,decomp%lupri,decomp%S,decomp%U,decomp%U_inv)
#else
      call lsquit('Qschultz does not support quadruple precision',decomp%lupri)
#endif

      call lstimer('LWQITER',tmstart,tmend,decomp%lupri)
   else 
     WRITE(decomp%LUPRI,'(/A)') &
     &     'No type of decomposition specified '
          CALL lsQUIT(' No type of decomposition specified ',decomp%lupri)
   endif

    if (decomp%lcv_basis) then
      !compute  U = lc_CMO*(lc_CMO'*S*lc_CMO)^-0.5
      !        iU = U'*S
     
      call mat_assign(decomp%S,S_sav)
      call mat_assign(S_sav,decomp%U_inv) 

      call mat_mul(decomp%lcv_CMO,S_sav,'n','n',1d0,0d0,decomp%U_inv)
      call mat_mul(decomp%U_inv,decomp%S,'t','n',1d0,0d0,decomp%U)

      call mat_free(S_sav)
      call mat_free(decomp%lcv_CMO)
    endif

   end subroutine decomposition 

   !> \brief Transform gradient/residual/Fock matrix to OAO basis: F_OAO = U_inv_t * F_AO * U_inv
   !> \author S. Host
   !> \date 2005
   subroutine res_to_oao_basis(decomp, res, resU)
   use matrix_operations
   implicit none
        !> Contains settings for decomposition and OAO decomposed overlap
        type(decompItem),intent(in)  :: decomp
        !> Residual in AO basis
        type(Matrix), intent(in)     :: res
        !> Residual in OAO basis (output)
        type(Matrix), intent(inout)  :: resU
        type(Matrix)                 :: wrk

   call MAT_INIT(wrk,res%nrow,res%ncol)

   call mat_mul(decomp%U_inv,res,'t','n',1.d0,0.d0,wrk)
   call mat_mul(wrk,decomp%U_inv,'n','n',1.d0,0.d0,resU)

   call mat_free(wrk)
   end subroutine res_to_oao_basis

   !> \brief Transform gradient/residual/Fock matrix from OAO basis: F_AO = Ut * F_OAO * U
   !> \author S. Host
   !> \date 2005
   subroutine res_from_oao_basis(decomp, resU, res)
   use matrix_operations
   implicit none
        !> Contains settings for decomposition and OAO decomposed overlap
        type(decompItem),intent(in)  :: decomp
        !> Residual in OAO basis
        type(Matrix), intent(in)     :: resU
        !> Residual in AO basis  (output)
        type(Matrix), intent(inout)  :: res
        type(Matrix)                 :: wrk

   call MAT_INIT(wrk,res%nrow,res%ncol)

   call mat_mul(decomp%U,resU,'t','n',1.d0,0.d0,wrk)
   call mat_mul(wrk,decomp%U,'n','n',1.d0,0.d0,res)

   call mat_free(wrk)
   end subroutine res_from_oao_basis

   !> \brief Transform trial vector/density matrix to OAO basis: X_OAO = U * X_AO * Ut
   !> \author S. Host
   !> \date 2005
   subroutine x_to_oao_basis(decomp, x, xU)  !Used for arh
   use matrix_operations
   implicit none
        !> Contains settings for decomposition and OAO decomposed overlap
        type(decompItem),intent(in)  :: decomp
        !> Trial vector/density matrix in AO basis
        type(Matrix), intent(in)     :: x
        !> Trial vector/density matrix in OAO basis (output)
        type(Matrix), intent(inout)  :: xU 
        type(Matrix)                 :: wrk
                                                                                                                     
   call MAT_INIT(wrk,xU%nrow,xU%ncol)
   call mat_mul(decomp%U,x,'n','n',1.d0,0.d0,wrk)
   call mat_mul(wrk,decomp%U,'n','t',1.d0,0.d0,xU)
   call mat_free(wrk)
   end subroutine x_to_oao_basis

   !> \brief Transform trial vector/density matrix from OAO basis: X_AO = U_inv * X_OAO * U_inv_t
   !> \author S. Host
   !> \date 2005
   subroutine x_from_oao_basis(decomp, xU, x)
   use matrix_operations
   implicit none
        !> Contains settings for decomposition and OAO decomposed overlap
        type(decompItem),intent(in)  :: decomp
        !> Trial vector/density matrix in OAO basis
        type(Matrix), intent(in)     :: xU
        !> Trial vector/density matrix in AO basis (output)
        type(Matrix), intent(inout)  :: x
        type(Matrix)                 :: wrk

   call MAT_INIT(wrk,xU%nrow,xU%ncol)
   call mat_mul(decomp%U_inv,xU,'n','n',1.d0,0.d0,wrk)
   call mat_mul(wrk,decomp%U_inv,'n','t',1.d0,0.d0,x)
   call mat_free(wrk)
   end subroutine x_from_oao_basis

   !> \brief Construct OAO matrices FU, DU, QU, FUQ, FUQ defined in type DecompItem above. 
   !> \author S. Host
   !> \date 2005
   subroutine get_oao_transformed_matrices(decomp,F,D)
   use matrix_operations
   implicit none
        !> Contains settings for decomposition and OAO decomposed overlap
        type(decompItem),intent(inout)          :: decomp
        !> Fock/KS matrix in AO basis
        type(Matrix), intent(in)                :: F
        !> Density matrix in AO basis
        type(Matrix), intent(in)                :: D
        type(Matrix)                            :: wrk
        integer                                 :: ndim

      ndim = decomp%S%nrow
      call mat_init(wrk,ndim,ndim)

      if (decomp%cfg_do_in_oao) then
         decomp%DU = D
         decomp%FU = F
      else
         call mat_mul(decomp%U,D,'n','n',1.d0,0.d0,wrk)
         call mat_mul(wrk,decomp%U,'n','t',1.d0,0.d0,decomp%DU)
         !write(lupri,*) 'DU:'
         !call mat_print(DU,1,FU%nrow,1,FU%ncol,lupri)

         call mat_mul(decomp%U_inv,F,'t','n',1.d0,0.d0,wrk)
         call mat_mul(wrk,decomp%U_inv,'n','n',1.d0,0.d0,decomp%FU)
      endif

       call mat_add_identity(1.0d0, -1.0d0, decomp%DU, decomp%QU)
      !write(lupri,*) 'FU:'
      !call mat_print(FU,1,FU%nrow,1,FU%ncol,lupri)

       call mat_mul(decomp%DU,decomp%FU,'n','n',1.d0,0.d0,wrk)
       call mat_mul(wrk,decomp%DU,'n','n',1.d0,0.d0,decomp%FUP)
      !write(lupri,*) 'FUP:'
      !call mat_print(FUP,1,FUP%nrow,1,FUP%ncol,lupri)

       call mat_mul(decomp%QU,decomp%FU,'n','n',1.d0,0.d0,wrk)
       call mat_mul(wrk,decomp%QU,'n','n',1.d0,0.d0,decomp%FUQ)
      !write(lupri,*) 'FUQ:'
      !call mat_print(FUQ,1,FUQ%nrow,1,FUQ%ncol,lupri)

   call mat_free(wrk)
   end subroutine get_oao_transformed_matrices

   !> \brief Project out redundancies in OAO basis, P(X) = DU*X*QU + QU*X*DU.
   !> \author S. Host
   !> \date 2005
   subroutine project_oao_basis(decomp, X, symm, X_proj)
   use matrix_operations
   implicit none
       !> Contains settings for decomposition and OAO decomposed overlap
       type(decompItem),intent(in) :: decomp       
       !> Matrix to be projected (OAO basis)
       TYPE(matrix), intent(in)    :: X
       !> Symmetry of X, 0 = nonsymmetric, 1 = symmetric, 2 = antisymmetric
       integer, intent(in)         :: symm
       !> Projected matrix in OAO basis (output)
       TYPE(matrix), intent(inout) :: X_proj
       TYPE(matrix)                :: scr, PXQ, QXP

   call MAT_INIT(PXQ,X%nrow,X%ncol)
   call MAT_INIT(scr,X%nrow,X%ncol)

   if (symm == 1 .or. symm == 2) then
      call MAT_MUL(decomp%DU,X,'n','n',1.0d0,0.0d0,scr)
      call MAT_MUL(scr,decomp%QU,'n','n',1.0d0,0.0d0,PXQ)

      X_proj = PXQ

      call mat_trans(PXQ, scr)
       if (symm == 2) then 
         call mat_scal(-1.0d0,scr)
      endif
      call mat_daxpy(1.0d0,scr,X_proj)
   else
      call MAT_INIT(QXP,X%nrow,X%ncol)

      call MAT_MUL(decomp%DU,X,'n','n',1.0d0,0.0d0,scr)
      call MAT_MUL(scr,decomp%QU,'n','n',1.0d0,0.0d0,PXQ)
      call MAT_MUL(decomp%QU,X,'n','n',1.0d0,0.0d0,scr)
      call MAT_MUL(scr,decomp%DU,'n','n',1.0d0,0.0d0,QXP)

      call MAT_ADD(1.0d0,PXQ,1.0d0,QXP,X_proj)

      call MAT_FREE(QXP)
   endif

   call MAT_FREE(PXQ)
   call MAT_FREE(scr)
   end subroutine project_oao_basis
#endif /* ifndef PRG_DIRAC */

end module decompMod

!radovan: this module is used in the response code
!         which is used also in DIRAC
!         this module depends on many things that DIRAC
!         does not have - i have to skip all the following
!         to allow compilation in DIRAC
#ifndef PRG_DIRAC
!> \brief Lowdin decomposition by diagonalization.
!> \author B. Jansik
!> \date 2006
subroutine lowdin_diagonalize(lupri,unres,S,U,U_inv)
  use matrix_operations
  use lowdin_module
  implicit none
  !> Logical unit number for output
  integer, intent(in) :: lupri
  !> True if unrestricted calculation
  logical, intent(in) :: unres
  !> Overlap matrix
  type(matrix), intent(in) :: S
  !> Decomposition of overlap matrix, S = U * Ut
  type(matrix), intent(in) :: U
  !> Inverse of decomposition of overlap matrix
  type(matrix), intent(in) :: U_inv
  real(realk), allocatable :: S_full(:,:), Ut_full(:,:), Ut_inv_full(:,:)
  integer dim

  dim = S%ncol

  ! We decompose the alpha part of the matrix which we expect to be spin
  ! independent.
  if(unres) then
     allocate(S_full(dim,2*dim))
     call mat_to_full2(S, 1d0, S_full)
  else
     allocate(S_full(dim,dim))
     call mat_to_full(S,1.0d0,S_full)
  end if

  allocate(Ut_inv_full(dim,dim), Ut_full(dim,dim))
  call lowdin_diag(S%nrow,S_full,Ut_full,Ut_inv_full,lupri)
  call mat_set_from_full(Ut_inv_full, 1d0, U_inv)
  call mat_set_from_full(Ut_full,     1d0, U)
  deallocate(S_full, Ut_inv_full, Ut_full)
end subroutine lowdin_diagonalize
#endif /* ifndef PRG_DIRAC */
