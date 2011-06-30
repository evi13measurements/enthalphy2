module soeo_loop

use files
use soeo_transform
use soeo_dens
use soeo_stability
use soeo_redspace
use matrix_util
use Fock_evaluator
use Matrix_Operations
use soeo_typedef
use soeo_itemroutines
use soeo_debug

implicit none

!Contains:
!  soeoloop

contains

!> \brief Makes Second Order Ensemble Optimization
!> \author C. Nygaard
!> \date 2010
!> \param Fao The Fock matrix in AO basis
!> \param H1 The one-electron part of the Fock matrix
!> \param S The overlap matrix
!> \param Dao The density matrix in the AO basis
!
!> Ordinary RH-optimization (arh or similar) should be done first
!>  FC = SCe
!> such that a coefficient matrix C is available
!> together with an AO-density-matrix and an overlap matrix
!=======================================================================
subroutine soeoloop (Fao, H1, S, Dao, decomp)

implicit none

!Needed matrices
type(decompItem), intent(in) :: decomp
type(soeoItem)               :: soeo
type(matrix), intent(in)     :: H1, S, Dao, Fao
type(matrix)                 :: K, deltatheta
type(matrix)                 :: Dmo_diff
type(matrix)                 :: grad_m, grad_v
real(realk), allocatable     :: occs(:), orbE(:)
real(realk)                  :: ratio, Dmo_diffnorm, gradnorm
real(realk)                  :: thresh, trust, mu !
integer                      :: maxiter           ! Those are for the red-space
integer                      :: Nelec, Nocc, Nact, Nvirt
integer                      :: i, j, m, l, iter
logical                      :: stable, loopdone, tmplog, aufbau, grandcan
type(matrix)                 :: tmpmat, tmpvec

!Other
real(realk)                  :: tstart, tloop1, tloop2, tend, dt
real(realk)                  :: tmp, val1, val2, x
real(realk), allocatable     :: hessian(:,:), gradient(:), step(:), smallhessian(:,:)
real(realk), allocatable     :: hessian2(:,:), gradient2(:)
integer                      :: tmpint, tmpdim

        write (decomp%lupri, *) "SOEOLOOP STARTED!"
        write (decomp%lupri, *) "============================================="

call cpu_time (tstart)

!Initializations:
!-----------------------------------------------------------------------
call soeoItem_init (Dao, Fao, S, decomp, soeo)
print *, 'soeoItem_init done'


!DEBUGGING PRINTS

write (decomp%lupri, *) 'Old Dao ='
call mat_print (Dao, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'Old Fao ='
call mat_print (Fao, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'Old S ='
call mat_print (S, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'My Dao ='
call mat_print (soeo%Dao, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'My Fao ='
call mat_print (soeo%Fao, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'My Dmo ='
call mat_print (soeo%Dmo, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'My Fmo ='
call mat_print (soeo%Fmo, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'My C ='
call mat_print (soeo%C, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)

!END DEBUGGING PRINTS


if (decomp%cfg_unres) then
  Nelec = 2*decomp%nocc + decomp%nactive
else
  Nelec = decomp%nocc
endif
call mat_init (K        , soeo%Nbast, soeo%Nbast)
call mat_init (Dmo_diff , soeo%Nbast, soeo%Nbast)
call mat_init (grad_m   , soeo%Nbast, soeo%Nbast)
allocate (occs(soeo%Nbast))
if (decomp%cfg_unres) then
  allocate (orbE(2*soeo%Nbast))
else
  allocate (orbE(soeo%Nbast))
endif

print *, 'initializations done'
!-----------------------------------------------------------------------


!Printing active space:
!-----------------------------------------------------------------------
write (decomp%lupri, *)
if (decomp%cfg_unres) then
  write (decomp%lupri, *) "Number of electrons:", Nelec
else
  write (decomp%lupri, *) "Number of electron pairs:", Nelec
endif
write (decomp%lupri, *) "Number of basis-functions:", soeo%Nbast
write (decomp%lupri, *) "Size of occupied space:", soeo%Nocc
write (decomp%lupri, *) "Size of active space:", soeo%Nact
write (decomp%lupri, *) "Size of virtual space:", soeo%Nbast - soeo%Nocc - soeo%Nact
write (decomp%lupri, *)
!-----------------------------------------------------------------------

!Testing if we have an aufbau solution
! - and setting the new active space
!Performing linesearch along active occupation numbers if not aufbau
!-----------------------------------------------------------------------
write (decomp%lupri, *) "Testing if we have an aufbau solution"

!FOR DEBUGGING

grandcan = .true.
Nocc = soeo%Nocc ; Nact = soeo%Nact ; Nvirt = soeo%Nbast - soeo%Nocc - soeo%Nact
if (grandcan) then
  call soeo_checkorbs (decomp%cfg_unres, soeo%Dmo, soeo%Fmo, Nocc, Nact, Nvirt, aufbau)
  call soeoItem_update_space (Nocc, Nact, decomp, soeo)
else
  call soeo_isit_aufbau (decomp%cfg_unres, soeo%Dmo, soeo%Fmo, Nocc, Nact, Nvirt, aufbau)
  call soeoItem_update_space (Nocc, Nact, decomp, soeo)
endif

!!We have H-atom, Nbast = 5 with 1sa = 0, 1sb = 1, we want 1sa = x, 1sb = 1
!! (0 <= x <= 1)
!! and we want the occupations not to chagne from this
!
!call mat_zero (soeo%Dmo)
!x = 0.0d0
!call mat_create_ab_elms (1, 1, (/ x, 1.0d0 /), soeo%Dmo)
!call soeo_oldtheta_init (soeo%Dmo, soeo%oldtheta, decomp)
!call soeoItem_update_space (0, 0, decomp, soeo)
!
!call soeo_getnew_Dao (soeo%Dmo, soeo%C, soeo%S, K, soeo%Dao)
!call FCK_get_Fock (soeo%Dao, soeo%Fao, soeo%Etotal)
!
!write (decomp%lupri, *) 'x =', x


write (decomp%lupri, *)
if (decomp%cfg_unres) then
  write (decomp%lupri, *) "Number of electrons:", Nelec
else
  write (decomp%lupri, *) "Number of electron pairs:", Nelec
endif
write (decomp%lupri, *) "Number of basis-functions:", soeo%Nbast
write (decomp%lupri, *) "Size of occupied space:", soeo%Nocc
write (decomp%lupri, *) "Size of active space:", soeo%Nact
write (decomp%lupri, *) "Size of virtual space:", soeo%Nbast - soeo%Nocc - soeo%Nact
write (decomp%lupri, *)

call mat_init (deltatheta, soeo%Nact, 1)
call mat_init (grad_v, soeo%Nact, 1)
!aufbau = .false.

if (aufbau) then
  write (decomp%lupri, *) "We have an aufbau solution"
else
  write (decomp%lupri, *) "We do not have an aufbau solution"
  write (decomp%lupri, *) "   - performing linesearch"

  if (decomp%cfg_unres) then
    allocate (step(2))
  else
    allocate (step(1))
  endif
  tmp = 3.141592653589793 / 4.0d0
  do i=soeo%Nocc+1,soeo%Nocc+soeo%Nact
    call mat_get_ab_elms (soeo%oldtheta, i, 1, step)
    step = step + tmp
    call mat_create_ab_elms (i, 1, step, soeo%oldtheta)
  enddo
  deallocate (step)

!  call soeo_linesearch (soeo, decomp)
endif
!END DEBUGGING
!-----------------------------------------------------------------------

call mat_zero (K) ; call mat_zero (deltatheta)
loopdone = .false.
iter = 0
mu = 0.0d0

if (.not. aufbau) then
  !Start the soeo-loop
  !---------------------------------------------------------------------
  loopdone = .true.
  do
    call cpu_time (tloop1)
    iter = iter + 1
    write (decomp%lupri, *)
    write (decomp%lupri, *) "soeo-iter =", iter

    !Get the new matrices:
    !-------------------------------------------------------------------
    call soeoItem_stepforth (K, deltatheta, decomp, soeo)
    !-------------------------------------------------------------------

!DEBUGGING
write (decomp%lupri, *) "norm(K,deltatheta), trust =", soeo_norm (K, deltatheta), soeo%trust
!END DEBUGGING

    !Test for convergence:
    !-------------------------------------------------------------------
    call soeo_get_gradient (soeo, grad_m, grad_v, decomp)
    gradnorm = soeo_norm (grad_m, grad_v)
    write (decomp%lupri, *) "gradnorm =", gradnorm

!    call mat_add (1.0d0, soeo%old_Dmo, -1.0d0, soeo%Dmo, Dmo_diff)
!    Dmo_diffnorm = mat_sqnorm2 (Dmo_diff)
!    Dmo_diffnorm = sqrt(Dmo_diffnorm)
!    write (decomp%lupri, *) "Dmo_diffnorm =", Dmo_diffnorm

    write (decomp%lupri, *) "Etotal and dE:", soeo%Etotal, soeo%dE
  
    if (gradnorm < soeo%macrothresh .and. iter > 1) then
      exit
    endif
    if (iter > soeo%macromaxiter) then
      print *, "Maximum soeo-iterations reached - sorry"
      exit
    endif
    !-------------------------------------------------------------------

    !Stats
    !-------------------------------------------------------------------
    if (iter == 1) then
      write (decomp%lupri, '(a103)') " OOO  iter    Etotal        &
                                          & dEpred        &
                                          & dE            &
                                          & ratio         &
                                          & mu            &
                                          & gradnorm      "
      write (decomp%lupri, '(a103)') " OOO------------------------&
                                          &---------------&
                                          &---------------&
                                          &---------------&
                                          &---------------&
                                          &---------------"
      print ('(a100)'), "   iter    Etotal        &
                                          & dEpred        &
                                          & dE            &
                                          & ratio         &
                                          & mu            &
                                          & gradnorm      "
      print ('(a100)'), " ------------------------&
                                          &---------------&
                                          &---------------&
                                          &---------------&
                                          &---------------&
                                          &---------------"
    endif
!    if (soeo%dE > 1.0d-7 .or. soeo%dE < -1.0d-7) then
    if (soeo%dE /= 0.0d0) then
      ratio = soeo%dEpred/soeo%dE
    else
      ! ratio is not changed from last iteration (because last iteration was rejected)
    endif
    write (decomp%lupri, '(" OOO ", i4, 6f15.8)') iter, soeo%Etotal, soeo%dEpred, soeo%dE, ratio, mu, gradnorm
    print ('("  ", i4, 6f15.8)'), iter, soeo%Etotal, soeo%dEpred, soeo%dE, ratio, mu, gradnorm
    !-------------------------------------------------------------------

    !Update the trust-radius
    !-------------------------------------------------------------------
    if (iter > 1) then
      if (ratio > 0.75d0) then
        !increase trust radius
        soeo%trust = 1.2d0 * soeo%trust
      elseif (ratio < 0.25d0) then
        !decrease trust radius
        soeo%trust = 0.7d0 * soeo%trust
        if (ratio < 0.0d0) then
          !rejection
          if (soeo%dE > 0.0d0) then
            print *, 'Rejection: dE > 0.0d0'
            write (decomp%lupri, *) ' OOO Rejection: dE > 0.0d0'
          elseif (soeo%dEpred > 0.0d0) then
            print *, 'Rejection: dEpred > 0.0d0'
            write (decomp%lupri, *) ' OOO Rejection: dEpred > 0.0d0'
          endif
          call soeoItem_stepback (soeo)
        elseif (soeo%dE > 0.0d0) then
          !rejection
          print *, 'Rejection: dE > 0.0d0'
          write (decomp%lupri, *) ' OOO Rejection: dE > 0.0d0'
          call soeoItem_stepback (soeo)
        endif
      endif
    endif
    !-------------------------------------------------------------------

    !Solve the soeo-equations in reduced space
    !-------------------------------------------------------------------
    call mat_zero (K)

!FOR DEBUGGING
!tmpint = (soeo%Nbast * (soeo%Nbast + 1) / 2 - soeo%Nbast) + soeo%Nact
!if (decomp%cfg_unres) then
!  tmpint = 2 * tmpint
!endif
!allocate (hessian(tmpint,tmpint), gradient(tmpint))
!allocate (hessian2(tmpint,tmpint), gradient2(tmpint))
!if (decomp%cfg_unres) then
!  tmpdim = tmpint/2
!else
!  tmpdim = tmpint
!endif
!tmpdim = tmpdim - soeo%Nact
!call mat_init (tmpvec, tmpdim, 1)
!
!call soeo_get_full_hessian (soeo, hessian, decomp)
!write (decomp%lupri, *) 'hessian in iteration', iter
!call output (hessian, 1, tmpint, 1, tmpint, &
!           & tmpint, tmpint, 1, decomp%lupri)
!
!call mat_to_vec ('a', grad_m, tmpvec)
!if (decomp%cfg_unres) then
!  allocate (step(2))
!else
!  allocate (step(1))
!endif
!do i=1,tmpdim
!  call mat_get_ab_elms (tmpvec, i, 1, step)
!  gradient(i) = step(1)
!  if (decomp%cfg_unres) then
!    gradient(i+tmpint/2) = step(2)
!  endif
!enddo
!do i=1,soeo%Nact
!  call mat_get_ab_elms (grad_v, i, 1, step)
!  gradient(tmpdim+i) = step(1)
!  if (decomp%cfg_unres) then
!    gradient(tmpint/2+tmpdim+i) = step(2)
!  endif
!enddo
!write (decomp%lupri, *) 'gradient in iteration', iter
!call output (gradient, 1, tmpint, 1, 1, tmpint, 1, 1, decomp%lupri)
!deallocate (step)
!!Finite difference Hessian
!call soeo_finite_difference (soeo, gradient2, hessian2, decomp)
!write (decomp%lupri, *) 'hessian from finite difference in iteration', iter
!call output (hessian2, 1, tmpint, 1, tmpint, &
!           & tmpint, tmpint, 1, decomp%lupri)
!write (decomp%lupri, *) 'and the gradient'
!call output (gradient2, 1, tmpint, 1, 1, tmpint, 1, 1, decomp%lupri)
!
!write (decomp%lupri, *) 'difference between hessians:'
!hessian = hessian-hessian2
!call output (hessian, 1, tmpint, 1, tmpint, &
!           & tmpint, tmpint, 1, decomp%lupri)
!tmp = 0.0d0
!do i=1,tmpint
!  do j=1,tmpint
!    tmp = tmp + hessian(i,j)*hessian(i,j)
!  enddo
!enddo
!write (decomp%lupri, *) 'diffnorm hessian', tmp
!write (decomp%lupri, *) 'difference between gradients:'
!gradient = gradient-gradient2
!call output (gradient, 1, tmpint, 1, 1, tmpint, 1, 1, decomp%lupri)
!tmp=0.0d0
!do i=1,tmpint
!  tmp = tmp+gradient(i)
!enddo
!write (decomp%lupri, *) 'diffnorm gradient', tmp
!
!call mat_free (tmpvec)
!
!!allocate (smallhessian(tmpint/2,tmpint/2))
!!do i=1,tmpint/2
!!  do j=1,tmpint/2
!!    smallhessian(i,j) = hessian(i,j) + hessian(tmpint/2+i,tmpint/2+j)
!!    smallhessian(i,j) = 0.5d0 * smallhessian(i,j)
!!  enddo
!!enddo
!!write (decomp%lupri, *) 'average alpha - beta hessian in iteration', iter
!!call output (smallhessian, 1, tmpint/2, 1, tmpint/2, &
!!           & tmpint/2, tmpint/2, 1, decomp%lupri)
!!tmplog = .true.
!!do i=1,tmpint/2
!!  do j=1,tmpint/2
!!    if (smallhessian(i,j)-smallhessian(j,i)>1.0d-5 .or.&
!!         & smallhessian(j,i)-smallhessian(i,j)>1.0d-5) then
!!      tmplog = .false.
!!    endif
!!  enddo
!!enddo
!!if (tmplog) then
!!  write (decomp%lupri, *) 'average hessian is symmetric'
!!else
!!  write (decomp%lupri, *) 'average hessian is NOT symmetric'
!!endif
!!
!!!printing different parts of the hessian
!!write (decomp%lupri, *) 'tt-hessian, alpha'
!!call output (hessian, 11, 13, 11, 13, 26, 26, 1, decomp%lupri)
!!write (decomp%lupri, *) 'tt-hessian, beta'
!!call output (hessian, 24, 26, 24, 26, 26, 26, 1, decomp%lupri)
!
!deallocate (hessian, gradient, hessian2, gradient2)
!!deallocate (smallhessian)
!if (iter == 1) then
!!  stop 'VELOCIRAPTOR'
!endif
!!END DEBUGGING

    call soeo_solver (soeo, orbE, K, deltatheta, mu, decomp)
    !-------------------------------------------------------------------

    !Finding the predicted energy change
    !-------------------------------------------------------------------
    call soeo_get_Epred (soeo, K, deltatheta, decomp)
    !-------------------------------------------------------------------

    !Maybe some averaging here?

    !Energy?
    !-------------------------------------------------------------------
    call cpu_time (tloop2)
    dt = tloop2 - tloop1

    write (decomp%lupri, '("Iteration ", i3, " done in ", f15.10, " seconds")') &
                  & iter, dt
!    print '("Iteration ", i3, " done in ", f15.10, " seconds")', iter, dt

  enddo
  !---------------------------------------------------------------------
else
  write (decomp%lupri, *) "Your system is stable with the previously calculated&
                   & occupation numbers - "
  write (decomp%lupri, *) "No need for SOEO-calculation."
endif

call fck_get_fock (soeo%Dao, soeo%Fao, soeo%Etotal)
write (decomp%lupri, *) "Final energy (fck_get_fock)=", soeo%Etotal

write (decomp%lupri, *) 'Final Fmo ='
call mat_print (soeo%Fmo, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'Final C ='
call mat_print (soeo%C, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
write (decomp%lupri, *) 'Final gradient ='
call mat_print (grad_m, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)
call mat_print (grad_v, 1, soeo%Nact, 1, 1, decomp%lupri)
write (decomp%lupri, *) 'Final Dmo ='
call mat_print (soeo%Dmo, 1, soeo%Nbast, 1, soeo%Nbast, decomp%lupri)

do i=1,soeo%Nbast
  occs(i) = soeo_diag_elma (soeo%Dmo, i)
enddo
write (decomp%lupri, *) "Occupation numbers:", occs
write (decomp%lupri, *) "sum:", sum(occs)
if (decomp%cfg_unres) then
  do i=1,soeo%Nbast
    occs(i) = soeo_diag_elmb (soeo%Dmo, i)
  enddo
  write (decomp%lupri, *) "Occupation numbers in beta orbitals:", occs
  write (decomp%lupri, *) "sum beta:", sum(occs)
  write (decomp%lupri, *) "Total number of electrons should be:", Nelec
else
  write (decomp%lupri, *) "Number of electron pairs should be:", Nelec
endif



!write (decomp%lupri, *) "Final MO-density matrix:"
!call mat_print (Dmo, 1, Nbast, 1, Nbast, decomp%lupri)
!write (decomp%lupri, *) "Final MO_Fock matrix:"
!call mat_print (Fmo, 1, Nbast, 1, Nbast, decomp%lupri)

!Deallocations:
!-----------------------------------------------------------------------
call mat_free (K)
call mat_free (Dmo_diff)
call mat_free (grad_m)
deallocate (occs, orbE)

call mat_free (deltatheta)
call mat_free (grad_v)

call soeoItem_shutdown (soeo)

!-----------------------------------------------------------------------

        print *, "SOEOLOOP FINISHED!"
        write (decomp%lupri, *) "SOEOLOOP FINISHED!"
        write (decomp%lupri, *) "Number of soeo-iterations =", iter
        write (decomp%lupri, *) "============================================="

call cpu_time (tend)
dt = tend - tstart
write (decomp%lupri, '("Time spend in soeoloop:", f15.10)') dt

end subroutine soeoloop
!=======================================================================

end module soeo_loop
