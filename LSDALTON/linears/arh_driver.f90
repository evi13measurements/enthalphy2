!> @file 
!> Contains Augmented Roothaan-Hall / Direct Density optimization drivers.

!> \brief Augmented RH/Direct dens driver.
!> \author S. Host
!> \date 2007
!>
!> ARH: S. Host, B. Jansik, J. Olsen et al. PCCP 10, 5344 (2008)   \n
!>      S. Host, J. Olsen, B. Jansik et al. JCP 129, 124106 (2008) \n
!>                                                                 \n
!> Direct density: P. Salek, S. Host, L. Thogersen et al. JCP 126, 114110 (2007)
!>
MODULE ARHmodule
   use direct_dens_util
   use files
   use queue_ops
   use arhDensity

contains

   !> \brief This routine calls the appropriate ARH solver 
   !> \author S. Host
   !> \date 2007
   !>
   !> Calculate the density that minimizes the Eepsilon = 2Tr(F*D(X)) + ARH terms. \n
   !> D ,the start guess, is updated during the optimization
   !> the optimized density is returned in D. .\n
   !> This routine is a bit obsolete (made more sense when it was possible to do
   !> more than one Newton iteration). Could be merged with driver.
   !>
   subroutine arh_get_density(arh,decomp,F,fifoqueue,Dnew,SCF_iteration)
   use matrix_util
    implicit none
    !> Contains solver info (ARH/TrFD)
    type(solverItem),intent(inout) :: arh
    !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
    type(decompItem),intent(inout) :: decomp
    !> Current Fock/KS matrix 
    type(Matrix), intent(in)       :: F
    !> Contains Fock/KS and density matrices from previous SCF iterations
    type(modFIFO),intent(inout)    :: fifoqueue
    !> Input: Current density matrix. Output: New density matrix.
    type(matrix), intent(inout)    :: Dnew 
    !> Current SCF iteration
    integer, intent(in)            :: SCF_iteration
    type(Matrix)                   :: x,wrk,wrk2,hes
    real(realk)                    :: gradnorm, xnorm
    integer                        :: ndim, ndens, i, hesdim
    real(realk)                    :: t1, t2
!Print variables
      real(realk), allocatable    :: weights(:)
      real(realk)                 :: actual_change_norm, expanded_change_norm
      type(Matrix), pointer       :: Dpointer, Fpointer, D0, Di, dummyp, Fi, F0
      type(debugItem)             :: debug
      type(DDitem)                :: DD

    ndim = decomp%S%nrow
    call mat_init(x,ndim,ndim)
    call mat_init(wrk,ndim,ndim)
    call mat_init(wrk2,ndim,ndim)

    !write (arh%lupri,*) 'Incoming density, arh_get_density, OAO basis:'
    !call MAT_PRINT(D, 1, D%nrow, 1, D%ncol, arh%LUPRI)

    if (arh%set_arhterms) then
       ndens = fifoqueue%offset
       if (ndens > 0) then 
          allocate(arh%fifometric(ndens,ndens))
          allocate(arh%inv_fifometric(ndens,ndens))
          allocate(arh%fifoM(ndens,ndens))
          call fifo_inverse_metric(arh,fifoqueue)
          call arh_get_M(arh,fifoqueue)
       endif
    endif
 
     !write(arh%lupri,*) 'F, AO:'
     !call mat_print(F,1,ndim,1,ndim,arh%lupri)
     !write(arh%lupri,*) 'D, AO:'
     !call mat_print(Dnew,1,ndim,1,ndim,arh%lupri)

    if (arh%step_accepted) call get_oao_transformed_matrices(decomp,F,Dnew)

    !FIXME: repair debug stuff for new input structure
    if (arh%debug_dd) then
       arh%scfit = SCF_iteration
       call dd_debug_homolumo(decomp,arh%diag_hlgap)
       if (SCF_iteration > arh%cfg_nits_debug) then
          call DD_init(decomp,DD)
          call DD_homolumo_and_heseigen(DD,decomp,debug,.false.,0,fifoqueue=fifoqueue)
          arh%iter_hlgap = debug%iter_hlgap
          arh%heseival   = debug%heseival
          arh%arheival   = debug%arheival
          call DD_shutdown(decomp,DD)
       endif
    else if (arh%debug_dd_homolumo) then
       call dd_debug_homolumo(decomp,debug%diag_hlgap)
    endif

    if (arh%debug_hessian) then
       if (decomp%cfg_unres) call lsquit('Debug routine not tested for unrestricted',decomp%lupri)
       hesdim = ndim*(ndim+1)/2 - ndim
       call mat_init(hes,hesdim,hesdim)
       call debug_get_hessian(arh,decomp,fifoqueue,hes)
       call util_diag(arh%lupri,hes,.false.,0.25d0,'Lowest Hessian eigenvalue:')
       call mat_free(hes)
    endif

    !write(arh%lupri,*) 'FU:'
    !call mat_print(decomp%FU,1,ndim,1,ndim,arh%lupri)
    !write(arh%lupri,*) 'DU:'
    !call mat_print(decomp%DU,1,ndim,1,ndim,arh%lupri)

    call get_OAO_gradient(decomp%FU, decomp%DU, wrk) !wrk = gradient
    call mat_scal(0.25d0,wrk) !To match linear transformation, also divided by 4!

    gradnorm = sqrt(mat_sqnorm2(wrk)) 
    write (arh%lupri,*) 'OAO gradnorm', gradnorm

    if (arh%cfg_nodamp) then
       call arh_PCG(arh, decomp, wrk, x, 0.0d0, antisymmetric,fifoqueue)
    else if (arh%cfg_arh_truncate .or. arh%cfg_arh_crop) then
       call arh_crop_solver(decomp, arh, debug, wrk, antisymmetric, x, fifoqueue)
    else
       call lsquit('Something wrong, no ARH solver chosen',decomp%lupri)
    endif
     !write(arh%lupri,*) 'New X **'
     !call mat_print(x,1,grad%nrow,1,grad%nrow,arh%lupri)

     !xnorm = SQRT(mat_sqnorm2(x))
     !write(arh%lupri,*) 'X norm orthogonal basis',xnorm
     if (arh%xnorm < 0.2d0 .and. SCF_iteration > 4) then
        arh%set_local = .true.
        if (arh%cfg_2nd_order_local .and. .not. arh%set_do_2nd_order) then
           write(arh%lupri,*) 'Local region - switching to second order optimization!'
           arh%set_do_2nd_order = .true.
           arh%cfg_arh_truncate = .false.
        endif
     endif

     !Now we construct the new density from X:
     CALL LSTIMER('START ',t1,t2,decomp%lupri)
     call oao_density_param(x,decomp%DU,wrk) !wrk = new D(X)
     call oao_purify(wrk,wrk2) !wrk2 = new purified density in oao basis
     CALL LSTIMER('NEW D ',t1,t2,arh%lupri)

     if (arh%set_arhterms) then
        !Determine part of the change in density which can be expanded in the subspace:
        arh%D_para_D_tot = 0.0d0 !So we don't get undefined output if ndens = 0
        if (ndens > 0) then
           allocate(weights(ndens))

           call arh_get_weights(arh,wrk2,x,fifoqueue,weights)

           call mat_zero(x) !x is now part of the change in density which can be expanded in the subspace
           call get_from_modFIFO(fifoqueue, 0, dummyp, D0)
           do i = 1, fifoqueue%offset
              call get_from_modFIFO(fifoqueue, i, dummyp, Di) 
              call MAT_ADD(1.0d0, Di, -1.0d0, D0, wrk)
              call mat_daxpy(weights(i),wrk,x)
           enddo

           !Actual change in density is Dchol - DU
           !Part of change which can be expanded is res
           call MAT_ADD(1.0d0, wrk2, -1.0d0, D0, wrk)  
           actual_change_norm = SQRT(mat_sqnorm2(wrk))
           expanded_change_norm = SQRT(mat_sqnorm2(x))

           arh%D_para_D_tot = expanded_change_norm/actual_change_norm
           deallocate(weights)
        endif
     endif

     !Stinne 24/1-07: Why tranform to AO basis? We already calculated the corresponding Fock matrix...
     !30/1-07: Because otherwise a wrong density is written in dens.restart!!
     if (decomp%cfg_do_in_oao) then
        call mat_assign(Dnew,wrk2)
     else
        call x_from_oao_basis(decomp,wrk2, Dnew) 
     endif

    if (associated(arh%fifometric)) then
       deallocate(arh%fifometric)
       nullify(arh%fifometric)
    endif
    if (associated(arh%inv_fifometric)) then
       deallocate(arh%inv_fifometric)
       nullify(arh%inv_fifometric)
    endif
    if (associated(arh%fifoM)) then
       deallocate(arh%fifoM)
       nullify(arh%fifoM)
    endif
    call mat_free(x)
    call mat_free(wrk)
    call mat_free(wrk2)
   end subroutine arh_get_density

!> \brief Reduced space solver used for solving Augmented Roothaan-Hall equations
!> \author S. Host
!> \date 2007
!>
!> The solver is based on the Conjugate Residual OPtimal vectors (CROP) scheme: \n
!>  M. Ziolkowski, V. Weijo, P. Jorgensen et al. JCP 128, 204105 \n
!> There are two options: 
!> - Run with the full subspace and keep vectors on disk (cfg_arh_crop=.true.) - this is invoked with .ARH FULL under
!> the *LINSCA section. 
!> - Run with a truncated number of vectors in memory (cfg_arh_truncate=.true.) 
!> (this is standard when using .ARH). Number of vectors to be kept may be set with .MICROVECS. The minimum (and default) is 2, corresponding
!> to two previous trial- and sigma-vectors. Plus the current ones, this means that 3 trial- and sigmas in total are stored. 
!>
!FIXME: Make it an option whether vectors should be kept on disk or in memory.
   subroutine arh_crop_solver(decomp, arh, debug, Grad, symm, x, fifoqueue)
   use direct_dens_util
   use levelshift_mod
      implicit none
      !> Contains decomposed overlap matrix (Löwdin, Cholesky or other)
      type(decompItem),intent(in)    :: decomp
      !> Contains solver info (ARH/TrFD)
      type(solverItem),intent(inout) :: arh
      !> If requested, some debug output is put in here to be printed later
      type(debugItem),intent(inout)  :: debug
      !> SCF gradient in orthonormal AO (OAO) basis
      type(Matrix), intent(in)    :: Grad
      !> If symm = 1, trial vector is symmetric. If symm = 2, trial vector is antisymmetric.
      integer, intent(in)         :: symm
      !> Contains Fock/KS and density matrices from previous SCF iterations
      TYPE(modFIFO),intent(inout) :: fifoqueue
      !> Output. The X that minimizes E = TrFD(X) + ARH terms(X)
      type(Matrix),intent(inout)  :: x
      type(Matrix)                :: scrmat, res, resP, xsave
      integer                     :: i, j, k, l, matdim, max_it, redspacedim_save, xsave_lu, ndens, dampdim, it
      integer                     :: lub, lusigma
      real(realk)                 :: err, errsave, t1, t2, mumax, thresh
      logical                     :: fileexists, done, do_LS
      real(realk)                 :: mu
!Truncate subspace:
      TYPE(modFIFO),target        :: vectorsubspace
      integer                     :: maxvec
!New levelshift
      type(Matrix)                :: xF, sigmaF
!Sparse1 matrices:
      real(realk)                 :: cutoff
!Queue on disk:
!      integer                     :: queue_lu
!      logical                     :: exppoint
      type(lshiftItem)             :: lshift

   if (arh%set_arhterms) then
      ndens = fifoqueue%offset
      write(arh%lupri,*) 'Number of densities in queue:', ndens
   endif
   done = .false.

   !The SCF energy from the previous iteration is found in the queue:
   matdim = Grad%nrow
   max_it = 200  !Stinne's first guess as to how many iterations are needed

   arh%set_optxelm = .false.
   call set_levelshift_config(arh,lshift)

   !Threshold for convergence:
   thresh = arh%cfg_micro_thresh*sqrt(mat_sqnorm2(Grad))
   if (matrix_type == mtype_sparse1) then
      call mat_inquire_cutoff(cutoff)
      if (thresh < 10.0d0*cutoff) then
      thresh = 10.0d0*thresh
      write(arh%lupri,*) 'Sparse matrices: Convergence threshold reset to:', thresh 
      endif
   endif

   if (arh%cfg_arh_truncate) then
      maxvec = arh%cfg_arh_microvecs
   else
      maxvec = max_it
   endif
   if (arh%cfg_arh_newdamp) then
      call mat_init(xF,matdim,matdim)
      call mat_init(sigmaF,matdim,matdim)
      allocate(arh%Ared(maxvec+1,maxvec+1), arh%Gred(maxvec+1))
      allocate(arh%Sred(maxvec+1,maxvec+1), arh%CROPmat(maxvec,maxvec))
      dampdim = maxvec+1
   else
      allocate(arh%Ared(maxvec,maxvec), arh%Gred(maxvec))
      allocate(arh%Sred(maxvec,maxvec), arh%CROPmat(maxvec,maxvec))
      dampdim = maxvec
   endif
      
   call mat_init(scrmat,matdim,matdim)
   call mat_init(res,matdim,matdim)
   call mat_init(resP,matdim,matdim)

   mu = 0.0d0 !Stinne test 1/7-09: Try to remove calc of HL gap, it takes too long...
   if (arh%cfg_fixed_shift) mu = -arh%cfg_fixed_shift_param

   CALL LSTIMER('START ',t1,t2,arh%lupri)
   if (arh%info_lineq) then 
      if (arh%cfg_arh_truncate) then
         WRITE(arh%LUPRI, "('CROP scheme, truncated subspace (CORE) - max. no. of vectors is', i3)") arh%cfg_arh_microvecs
      else if (arh%cfg_arh_crop) then
         WRITE(arh%LUPRI, "('CROP scheme with full subspace (DISK)')")
      endif
      if (.not. arh%cfg_arh_crop_safe) then
         WRITE(arh%LUPRI, "('WARNING: Double preconditioning has been removed!')")
      endif
   endif
   arh%trustradius_decreased = .false.

   !The maximum number of rejections is currently set to six
   ! - if there are that many rejections, the calculation is definitely
   ! unhealthy and should be stopped. This is almost always caused by lack of
   ! integral accuracy.
   if (arh%Nrejections > 6) then
     WRITE(arh%LUPRI, "('Too many rejections - probably related to lack of integral accuracy!')")
     CALL lsQUIT('Too many rejections - probably related to lack of integral accuracy!',decomp%lupri)
   endif

   arh%Ared = 0.0d0 ; arh%Gred = 0.0d0
   arh%Sred = 0.0d0 ; arh%CROPmat = 0.0d0

   if (arh%cfg_arh_truncate) then
      call modfifo_init(vectorsubspace,maxvec+1,matdim,arh%cfg_arh_disk_micro) 
   else
      lusigma = -1 ; lub = -1
      CALL LSOPEN(lusigma,'sigmavecs','unknown','UNFORMATTED')
      CALL LSOPEN(lub,'bvecs','unknown','UNFORMATTED')
   endif
   !gradnorm = SQRT(mat_sqnorm2(Grad))

   !Save first trial b vector = gradient:
   !call mat_copy(-1.0d0, Grad, scrmat)
   !call project_oao_basis(decomp, scrmat, symm, x)
   call project_oao_basis(decomp, Grad, symm, x)
   call mat_scal(-1.0d0,x)
   !First linear transformation:
   call arh_lintrans(arh,decomp,x,symm,0.0d0,scrmat,fifoqueue)

   if (arh%cfg_arh_truncate) then
      call add_to_modFIFO(vectorsubspace, x, scrmat, 0.0d0, .false.)
   else
      call mat_write_to_disk(lub,x)
      call mat_write_to_disk(lusigma,scrmat)
   endif
   arh%Gred(1) = mat_dotproduct(x,Grad)
   arh%Sred(1,1) = mat_dotproduct(x,x)
   arh%Ared(1,1) = mat_dotproduct(x,scrmat)

   !We need to call levelshift here, otherwise mu is wrong when constructing 1st residual
   call main_levelshift(lshift,arh%Ared,arh%Sred,arh%Gred,maxvec,lub,1,matdim,mu,vectorsubspace=vectorsubspace)
   write(arh%lupri,*) 'First mu:', mu 
   call mat_add(-1.d0,Grad,-1.d0, scrmat, res) !res = 1st residual
   call mat_daxpy(mu,x,res)

   !Preconditioning:
   call arh_precond(arh,decomp,res,symm,mu,resP)

   arh%CROPmat(1,1) = mat_dotproduct(res,resP)
   if (arh%info_crop) then
      write(arh%lupri,*) 'First elements in reduced space:'
      write(arh%lupri,*) 'Sred(1,1) =', arh%Sred(1,1) 
      write(arh%lupri,*) 'Ared(1,1) =', arh%Ared(1,1) 
      write(arh%lupri,*) 'CROPmat(1,1) =', arh%CROPmat(1,1) 
   endif

   !if (arh%debug_diag_redspace) then
   !   write(arh%lupri,*) 'First element in reduced space:', arh%Ared(2,2)
   !   debug%final_redspace_eival = arh%Ared(2,2)
   !endif

   j = 0 ; l = 0 ; k = 5 ; i = 0
   do i = 1, max_it
      if (i == max_it) then
         if (arh%set_optxelm) then
            if (arh%INFO_LINEQ) write (arh%lupri,*) 'MAXIT reached: Optimization of xmax failed, I will use old x'
            call mat_init(xsave,matdim,matdim)
            rewind(xsave_lu)
            call mat_read_from_disk(xsave_lu,xsave)
            call mat_max_elm(xsave, arh%maxelm)
            arh%final_redspacedim = redspacedim_save
            mu = mumax
            x = xsave
            call mat_free(xsave)
            exit
         else
            WRITE(arh%LUPRI,'(/A)') &
            &     'ARH: Linear equations not converged'
            CALL lsQUIT('ARH: Linear equations not converged',decomp%lupri)
         endif
      endif
      if (i > 1) call main_levelshift(lshift,arh%Ared,arh%Sred,arh%Gred,dampdim,lub,&
                             & i,matdim,mu,xF,sigmaF,vectorsubspace)
      !call arh_crop_x_and_res(fifoqueue,Grad,symm,mu,b_current,scrmat,resP,res)
      call arh_crop_x_and_res(arh,decomp,fifoqueue,Grad,symm,mu,x,scrmat,resP,res)
      !res = resP Stinne bugfix February 2010
      !x = b_current
      !Convergence?
      err = SQRT(mat_sqnorm2(res))
      it = i !"A do-variable within a DO body shall not appear in a variable definition context" pointed out by Thomas
      call arh_test_convergence(arh,err,errsave,x,mu,mumax,thresh,it,j,k,xsave_lu,redspacedim_save,maxvec,do_LS,done)
      lshift%optxelm = arh%set_optxelm
      if (do_LS) call main_levelshift(lshift,arh%Ared,arh%Sred,arh%Gred,dampdim,lub,i,matdim,mu,xF,sigmaF,vectorsubspace)
      if (done) exit
      !call arh_crop_intermed_sub(vectorsubspace,i,symm,Grad,mu,resP,b_current,scrmat,xF,sigmaF)
      call arh_crop_intermed_sub(arh,decomp,lub,lusigma,vectorsubspace,i,symm,Grad,mu,res,x,scrmat,xF,sigmaF) !resP -> res, Stinne bugfix February 2010
      !Save current on disc:
      if (arh%cfg_arh_truncate) then
         call add_to_modFIFO(vectorsubspace,x,scrmat,0.0d0,.false.)
      else
         call mat_write_to_disk(lub,x)
         call mat_write_to_disk(lusigma,scrmat)
      endif
      !call arh_crop_setup_redsp(vectorsubspace,i,symm,Grad,mu,b_current,scrmat,resP,xF,sigmaF)
      call arh_crop_setup_redsp(arh,decomp,lub,lusigma,vectorsubspace,i,symm,Grad,mu,x,scrmat,resP,res,xF,sigmaF)
   enddo

   arh%current_mu = mu

   !if (arh%debug_diag_redspace .and. .not. arh%cfg_arh_crop) then
   !   write (arh%lupri,'("Final lowest redspace eigenvalue: ",F12.6, "      ")') debug%final_redspace_eival
   !endif

   CALL LSTIMER('CROP solver',t1,t2,arh%lupri)

   if (arh%cfg_arh_crop) call mat_scal(-1.0d0,x) !HACK!!!!

   call arh_get_TR_denom(arh, Grad, x, scrmat, decomp%cfg_unres, mu)

   if (arh%cfg_arh_truncate) then
      call modfifo_free(vectorsubspace)
   else
      CALL LSCLOSE(lusigma,'DELETE')
      CALL LSCLOSE(lub,'DELETE')
   endif

   INQUIRE(file='xsave',EXIST=fileexists)
   if (fileexists) call LSCLOSE(xsave_lu,'DELETE')

   if (arh%cfg_arh_newdamp) then
      call mat_free(xF)
      call mat_free(sigmaF)
   endif
   deallocate(arh%Ared, arh%Gred)
   deallocate(arh%Sred, arh%CROPmat)
   call mat_free(scrmat)
   call mat_free(res)
   call mat_free(resP)
   end subroutine arh_crop_solver

!> \brief Set the parameters that must be passed to level shift module
!> \author S. Host
!> \date March 2010
   subroutine set_levelshift_config(arh,lshift)
   use levelshift_mod
   implicit none
      !> Contains solver info (ARH/TrFD)
      type(solverItem),intent(in) :: arh
      !> Contains info used for level shifting
      type(lshiftItem),intent(inout) :: lshift

      lshift%lupri              = arh%lupri
      lshift%fixed_shift        = arh%cfg_fixed_shift       
      lshift%fixed_shift_param  = arh%cfg_fixed_shift_param 
      lshift%arh_crop           = arh%cfg_arh_crop          
      lshift%arh_truncate       = arh%cfg_arh_truncate      
      lshift%arh_newdamp        = arh%cfg_arh_newdamp       
      lshift%min_lshift         = arh%cfg_min_lshift        
      lshift%max_element        = arh%set_max_element
      lshift%max_step           = arh%set_max_step
      lshift%optxelm            = arh%set_optxelm
      lshift%info_levelshift    = arh%info_levelshift
   end subroutine set_levelshift_config

! ------------ Debug section! ---------------!

!> \brief Print debug info.
!> \author S. Host
!> \date 2007
   subroutine arh_debug_print(arh)
   use arhDensity
   implicit none
        !> Contains solver info (ARH/TrFD)
        type(solverItem),intent(in) :: arh

   if (.not. arh%step_accepted) then
      write (arh%lupri, "(i5, F36.6, F12.6, '                ///')") &
           & arh%scfit, arh%final_redspace_eival, arh%arheival
   else          
      write (arh%lupri, "(i5, F12.6, F12.6, F12.6, F12.6, F12.6, '    ///')") &
           & arh%scfit, arh%diag_hlgap, arh%iter_hlgap, arh%final_redspace_eival, &
           & arh%arheival, arh%heseival
   endif
   end subroutine arh_debug_print

END MODULE ARHmodule
