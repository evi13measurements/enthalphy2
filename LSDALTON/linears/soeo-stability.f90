module soeo_stability

use soeo_matop
use soeo_transform
use soeo_dens
use soeo_redspace
use soeo_typedef
use Fock_evaluator

implicit none

!Contains:
!  soeo_isit_aufbau
!  soeo_checkorbs
!  soeo_linesearch

contains

!> \brief Tests from Dmo and Fmo if we have an aufbau solution
!> \author C. Nygaard
!> \date 2010-06-24
!> \param decomp Information
!> \param Dmo The MO density matrix
!> \param Fmo The MO Fock matrix (contains orbital energies
!> \param aufbau Logical, true if we have an aufbau solution
!=======================================================================
subroutine soeo_isit_aufbau (unres, Dmo, Fmo, Nocc, Nact, Nvirt, aufbau)

implicit none

logical, intent(in)      :: unres
type(matrix), intent(in) :: Dmo, Fmo
integer, intent(inout)   :: Nocc, Nact, Nvirt
logical, intent(out)     :: aufbau

integer                  :: Nbast, i
integer, allocatable     :: homo(:), lumo(:)
real(realk), allocatable :: ehomo(:), elumo(:), tmp1(:), tmp2(:)

Nbast = Dmo%nrow

if (unres) then
  allocate (tmp1(2), tmp2(2), homo(2), ehomo(2), lumo(2), elumo(2))
else
  allocate (tmp1(1), tmp2(1), homo(1), ehomo(1), lumo(1), elumo(1))
endif

!Set a minimum value of ehomo and a maximum value of elumo
homo = 0 ; ehomo = -1.0d5
lumo = 0 ; elumo = 1.0d5

do i=1,Nbast
  call mat_get_ab_elms (Dmo, i, i, tmp1)
  call mat_get_ab_elms (Fmo, i, i, tmp2)
  !For the alpha part
  if (tmp1(1)-1.0d0 < 1.0d-5 .and. tmp1(1)-1.0d0 > -1.0d-5) then !MO i is occupied
    if (tmp2(1) >= ehomo(1)) then
      ehomo(1) = tmp2(1)
      homo(1) = i
    endif
  elseif (tmp1(1) < 1.0d-5 .and. tmp1(1) > -1.0d-5) then !MO i is virtual
    if (tmp2(1) < elumo(1)) then
      elumo(1) = tmp2(1)
      lumo(1) = i
    endif
  else
    call lsquit ('Wrong density matrix in soeo_isit_aufbau',-1)
  endif
  !...and if unrestricted also for the beta part
  if (unres) then
    if (tmp1(2)-1.0d0 < 1.0d-5 .and. tmp1(2)-1.0d0 > -1.0d-5) then !MO i is occupied
      if (tmp2(2) >= ehomo(2)) then
        ehomo(2) = tmp2(2)
        homo(2) = i
      endif
    elseif (tmp1(2) < 1.0d-5 .and. tmp1(2) > -1.0d-5) then !MO i is virtual
      if (tmp2(2) < elumo(2)) then
        elumo(2) = tmp2(2)
        lumo(2) = i
      endif
    else
      call lsquit ('Wrong density matrix in soeo_isit_aufbau',-1)
    endif
  endif
enddo

aufbau = .true.
!check if ehomo < elumo
if (ehomo(1) > elumo(1)) then
  aufbau = .false.
endif
if (unres) then
  if (ehomo(2) > elumo(2)) then
    aufbau = .false.
  endif
endif
!set the new active space
print *, 'homo =', homo
print *, 'lumo =', lumo
if (unres) then
  if (lumo(1) < lumo(2)) then
    Nocc = lumo(1) - 1
  else
    Nocc = lumo(2) - 1
  endif
  if (homo(1) > homo(2)) then
    Nvirt = Nbast - homo(1)
  else
    Nvirt = Nbast - homo(2)
  endif
  Nact = Nbast - Nocc - Nvirt
else
  Nocc = lumo(1) - 1
  Nact = homo(1) - lumo(1) + 1
  Nvirt = Nbast - homo(1)
endif

deallocate (tmp1, tmp2, homo, ehomo, lumo, elumo)

end subroutine soeo_isit_aufbau
!=======================================================================

!> \brief Tests from Dmo and Fmo if we have a gcaufbau solution
!> \author C. Nygaard
!> \date 2010-06-24
!> \param decomp Information
!> \param Dmo The MO density matrix
!> \param Fmo The MO Fock matrix (contains orbital energies)
!> \param aufbau Logical, true if we have a gcaufbau solution
!=======================================================================
subroutine soeo_checkorbs (unres, Dmo, Fmo, Nocc, Nact, Nvirt, aufbau)

implicit none

logical, intent(in)      :: unres
type(matrix), intent(in) :: Dmo, Fmo
integer, intent(inout)   :: Nocc, Nact, Nvirt
logical, intent(out)     :: aufbau

integer                  :: Nbast, i
integer, allocatable     :: occplus(:), occminus(:), virtplus(:), virtminus(:)
real(realk), allocatable :: tmp1(:), tmp2(:)

Nbast = Dmo%nrow

if (unres) then
  allocate (occplus(2), occminus(2), virtplus(2), virtminus(2))
  allocate (tmp1(2), tmp2(2))
else
  allocate (occplus(1), occminus(1), virtplus(1), virtminus(1))
  allocate (tmp1(1), tmp2(1))
endif

!Set a minimum value of ehopo and a maximum value of eluno
occplus = 0 ; occminus = 0
virtplus = 0 ; virtminus = 0

do i=1,Nbast
  call mat_get_ab_elms (Dmo, i, i, tmp1)
  call mat_get_ab_elms (Fmo, i, i, tmp2)

  !For the alpha part
  if (tmp1(1)-1.0d0 < 1.0d-5 .and. tmp1(1)-1.0d0 > -1.0d-5) then !MO i is occupied
    if (tmp2(1) > 0.0d0) then
      occplus(1) = occplus(1) + 1
    else
      occminus(1) = occminus(1) + 1
    endif
  elseif (tmp1(1) < 1.0d-5 .and. tmp1(1) > -1.0d-5) then !MO i is virtual
    if (tmp2(1) < 0.0d0) then
      virtminus(1) = virtminus(1) + 1
    else
      virtplus(1) = virtplus(1) + 1
    endif
  else
    call lsquit ('Wrong density matrix in soeo_isit_aufbau',-1)
  endif
  !...and if unrestricted also for the beta part
  if (unres) then
    if (tmp1(2)-1.0d0 < 1.0d-5 .and. tmp1(2)-1.0d0 > -1.0d-5) then !MO i is occupied
      if (tmp2(2) > 0.0d0) then
        occplus(2) = occplus(2) + 1
      else
        occminus(2) = occminus(2) + 1
      endif
    elseif (tmp1(2) < 1.0d-5 .and. tmp1(2) > -1.0d-5) then !MO i is virtual
      if (tmp2(2) < 0.0d0) then
        virtminus(2) = virtminus(2) + 1
      else
        virtplus(2) = virtplus(2) + 1
      endif
    else
      call lsquit ('Wrong density matrix in soeo_isit_aufbau',-1)
    endif
  endif
enddo

aufbau = .true.
!Check if there are any positive occ or negative virt
if (occplus(1) /= 0 .or. virtminus(1) /= 0) then
  aufbau = .false.
endif
if (unres) then
  if (occplus(2) /= 0 .or. virtminus(2) /= 0) then
    aufbau = .false.
  endif
endif

!set the new active space
if (unres) then
  if (occminus(1) < occminus(2)) then
    Nocc = occminus(1)
  else
    Nocc = occminus(2)
  endif
  if (virtplus(1) < virtplus(2)) then
    Nvirt = virtplus(1)
  else
    Nvirt = virtplus(2)
  endif
  Nact = Nbast - Nocc - Nvirt
else
  Nocc = occminus(1)
  Nvirt = virtplus(1)
  Nact = Nbast - Nocc - Nvirt
endif
if (Nact /= 0) then
  aufbau = .false.
endif

!print *, 'Nocc, Nact, Nvirt =', Nocc, Nact, Nvirt
!print *, 'occminus =', occminus
!print *, 'occplus =', occplus
!print *, 'virtminus =', virtminus
!print *, 'virtplus =', virtplus

deallocate (tmp1, tmp2, occplus, occminus, virtplus, virtminus)

end subroutine soeo_checkorbs
!=======================================================================

!> \brief Performs line-search in the direction of negative eigenvalues
!> \author C. Nygaard
!> \date 2010
!> \param Nocc Number of occupied MO orbitals
!> \param Nact Number of active MO orbitals
!> \param S The overlap matrix
!> \param C The MO orbital coefficients
!> \param Fmo The Fock matrix in MO basis
!> \param Dmo The density matrix in MO basis
!> \param oldtheta The occupation angles to be changed
!> \param decomp Information
!=======================================================================
subroutine soeo_linesearch (soeo, decomp)

implicit none

!I/O:
type(soeoItem), intent(inout) :: soeo
type(decompItem), intent(in)  :: decomp
!Other:
real(realk), parameter        :: pi = 3.141592653589793
integer                       :: i, j, Nbast, iter, minenergy
type(matrix)                  :: step, theta
real(realk)                   :: alpha(3), energy(3)
real(realk), allocatable      :: tmp(:)
real(realk)                   :: thresh, a, num, denom
real(realk)                   :: alphaleft, energyleft
real(realk)                   :: alpharight, energyright
real(realk)                   :: newalpha, newenergy
real(realk)                   :: minalpha
type(matrix)                  :: Dmotheta, Fao, Dao
                                 !We need new Fao and Dao to find the
                                 ! energy for every Dmotheta

write (decomp%lupri, *) 'linesearch started'

if (decomp%cfg_unres) then
  allocate (tmp(2))
else
  allocate (tmp(1))
endif

thresh = 1.0d-8

call mat_init (Dmotheta, soeo%Nbast, soeo%Nbast)
call mat_init (Fao     , soeo%Nbast, soeo%Nbast)
call mat_init (Dao     , soeo%Nbast, soeo%Nbast)
call mat_init (step    , soeo%Nbast, 1)
call mat_init (theta   , soeo%Nbast, 1)

call mat_zero (Dmotheta)

!Creating the step direction
!-----------------------------------------------------------------------
call mat_zero (step)
do i=soeo%Nocc+1,soeo%Nocc+soeo%Nact
!  call mat_get_ab_elms (soeo%Dmo, i, i, tmp)
!  if (tmp(1)-1.0d0 < 1.0d-5 .and. tmp(1)-1.0d0 > -1.0d-5) then !occupied
!    tmp(1) = -1.0d0
!  elseif (tmp(1) < 1.0d-5 .and. tmp(1) > -1.0d-5) then !virtual
!    tmp(1) = 1.0d0
!  else
!    call lsquit ('Error in soeo_linesearch_2: Orbital is neither occ or virt!!',-1)
!  endif
!  if (decomp%cfg_unres) then
!    if (tmp(2)-1.0d0 < 1.0d-5 .and. tmp(2)-1.0d0 > -1.0d-5) then !occupied
!      tmp(2) = -1.0d0
!    elseif (tmp(2) < 1.0d-5 .and. tmp(2) > -1.0d-5) then !virtual
!      tmp(2) = 1.0d0
!    else
!      call lsquit ('Error in soeo_linesearch_2: Orbital is neither occ or virt!!',-1)
!    endif
!  endif
  tmp = 1.0d0
  call mat_create_ab_elms (i, 1, tmp, step)
enddo
!
!write (decomp%lupri, *) 'step ='
!call mat_print (step, 1, soeo%Nbast, 1, 1, decomp%lupri)
!
!-----------------------------------------------------------------------

!The first 3 guesses for alpha
!-----------------------------------------------------------------------
alpha(1) = 0.0d0 ; alpha(2) = pi/4.0d0 ; alpha(3) = pi/2.0d0

do i=1,3
  call mat_add (1.0d0, soeo%oldtheta, alpha(i), step, theta)
  call soeo_getnew_Dmo (theta, Dmotheta, decomp)
  call util_MO_to_AO_2 (soeo%S, soeo%C, Dmotheta, Dao, .false.)
  call FCK_get_fock (Dao, Fao, energy(i))
enddo

!print *, 'alpha =', alpha
!print *, 'energy =', energy

!-----------------------------------------------------------------------

iter = 0
do
  iter = iter + 1

  num = (energy(2)-energy(3))*(alpha(1)-alpha(3))
  num = num - (energy(1)-energy(3))*(alpha(2)-alpha(3))

  denom = (alpha(2)-alpha(1))*(alpha(2)-alpha(3))*(alpha(1)-alpha(3))

  a = num/denom

!print *, 'a =', a

  !The new alpha
  !Is the stationary point a maximum or a minimum?
  if (a > 0.0d0) then !minimum, newalpha = alpha(minimum)
    minalpha = - (energy(2)-energy(3))
    minalpha = minalpha / (2.0d0*a*(alpha(2)-alpha(3)))
    minalpha = minalpha + 0.5*(alpha(2)+alpha(3))

!print *, 'minalpha =', minalpha

    !The new alpha-vector
    if (minalpha < alpha(1)) then !newalpha --> alphaleft
  
      if (minalpha < 0.0d0) then
  
        alphaleft = 0.0d0
        call soeo_getnew_Dmo (soeo%oldtheta, Dmotheta, decomp)
        call util_MO_to_AO_2 (soeo%S, soeo%C, Dmotheta, Dao, .false.)
        call FCK_get_fock (Dao, Fao, energyleft)

        if (alpha(1) > 1.0d-5) then
          alpharight = alpha(1)
          energyright = energy(1)
        else
          alpharight = alpha(2)
          energyright = energy(2)
        endif

        newalpha = 0.5d0 * alpharight
        call mat_add (1.0d0, soeo%oldtheta, newalpha, step, theta)
        call soeo_getnew_Dmo (theta, Dmotheta, decomp)
        call util_MO_to_AO_2 (soeo%S, soeo%C, Dmotheta, Dao, .false.)
        call FCK_get_fock (Dao, Fao, newenergy)

      else
  
        alphaleft = minalpha
        call mat_add (1.0d0, soeo%oldtheta, alphaleft, step, theta)
        call soeo_getnew_Dmo (theta, Dmotheta, decomp)
        call util_MO_to_AO_2 (soeo%S, soeo%C, Dmotheta, Dao, .false.)
        call FCK_get_fock (Dao, Fao, energyleft)

        alpharight = alpha(2) ; energyright = energy(2)
        newalpha = alpha(1)   ; newenergy = energy(1)
  
      endif
  
    elseif (minalpha > alpha(3)) then !newalpha --> alpharight
  
      if (minalpha > pi/2.0d0) then
  
        alpharight = pi/2.0d0
        call mat_add (1.0d0, soeo%oldtheta, alpharight, step, theta)
        call soeo_getnew_Dmo (theta, Dmotheta, decomp)
        call util_MO_to_AO_2 (soeo%S, soeo%C, Dmotheta, Dao, .false.)
        call FCK_get_fock (Dao, Fao, energyright)

        if (alpha(3)-pi/2.0d0 < -1.0d-5) then
          alphaleft = alpha(3)
          energyleft = energy(3)
        else
          alphaleft = alpha(2)
          energyleft = energy(2)
        endif
        newalpha = 0.5d0 * (alphaleft + alpharight)
        call mat_add (1.0d0, soeo%oldtheta, newalpha, step, theta)
        call soeo_getnew_Dmo (theta, Dmotheta, decomp)
        call util_MO_to_AO_2 (soeo%S, soeo%C, Dmotheta, Dao, .false.)
        call FCK_get_fock (Dao, Fao, newenergy)
  
      else
  
        alpharight = minalpha
        call mat_add (1.0d0, soeo%oldtheta, alpharight, step, theta)
        call soeo_getnew_Dmo (theta, Dmotheta, decomp)
        call util_MO_to_AO_2 (soeo%S, soeo%C, Dmotheta, Dao, .false.)
        call FCK_get_fock (Dao, Fao, energyright)
  
        alphaleft = alpha(2) ; energyleft = energy(2)
        newalpha = alpha(3)  ; newenergy = energy(3)
  
      endif
  
    else !newalpha is between alpha(1) and alpha(3)
 
      newalpha = minalpha 
      call mat_add (1.0d0, soeo%oldtheta, newalpha, step, theta)
      call soeo_getnew_Dmo (theta, Dmotheta, decomp)
      call util_MO_to_AO_2 (soeo%S, soeo%C, Dmotheta, Dao, .false.)
      call FCK_get_fock (Dao, Fao, newenergy)
  
      do i=2,3
        if (alpha(i) > newalpha) then
          alpharight = alpha(i)  ; energyright = energy(i)
          alphaleft = alpha(i-1) ; energyleft = energy(i-1)
          exit
        endif
      enddo
  
    endif
  
  else !maximum, newalpha = halfway towards the lowest end of the curve
    minenergy = 1
    if (energy(3) < energy(minenergy)) then
      minenergy = 3
    endif

    newalpha = 0.5d0 * (alpha(2) + alpha(minenergy))
    call mat_add (1.0d0, soeo%oldtheta, newalpha, step, theta)
    call soeo_getnew_Dmo (theta, Dmotheta, decomp)
    call util_MO_to_AO_2 (soeo%S, soeo%C, Dmotheta, Dao, .false.)
    call FCK_get_fock (Dao, Fao, newenergy)

    if (minenergy == 1) then

      !alphaleft = 0.0d0
      alphaleft = 0.5d0 * alpha(1)
      call mat_add (1.0d0, soeo%oldtheta, alphaleft, step, theta)
      call soeo_getnew_Dmo (soeo%oldtheta, Dmotheta, decomp)
      call util_MO_to_AO_2 (soeo%S, soeo%C, Dmotheta, Dao, .false.)
      call FCK_get_fock (Dao, Fao, energyleft)

      alpharight = alpha(2)
      energyright = energy(2)

    elseif (minenergy == 3) then

      alphaleft = alpha(2)
      energyleft = energy(2)

      !alpharight = pi/2.0d0
      alpharight = 2.0d0*alpha(3)
      if (alpharight > pi/2.0d0) then
        alpharight = pi/2.0d0
      endif
      call mat_add (1.0d0, soeo%oldtheta, alpharight, step, theta)
      call soeo_getnew_Dmo (theta, Dmotheta, decomp)
      call util_MO_to_AO_2 (soeo%S, soeo%C, Dmotheta, Dao, .false.)
      call FCK_get_fock (Dao, Fao, energyright)

    endif

  endif

  !Convergence?

!print *, 'difference alpha(2) and newalpha =', newalpha-alpha(2)

  if (newalpha-alpha(2) < thresh .and. alpha(2)-newalpha < thresh) then
    exit
  endif

  alpha(1) = alphaleft  ; energy(1) = energyleft
  alpha(3) = alpharight ; energy(3) = energyright
  alpha(2) = newalpha   ; energy(2) = newenergy

!print *, 'alpha =', alpha
!print *, 'energy =', energy

enddo


!newalpha = pi/4.0d0

call mat_daxpy (newalpha, step, soeo%oldtheta)

        write (decomp%lupri, *) "Line_Search gives alpha =", newalpha

call mat_free (Dmotheta)
call mat_free (Fao)
call mat_free (Dao)
call mat_free (step)
call mat_free (theta)

deallocate (tmp)

end subroutine soeo_linesearch
!=======================================================================

end module soeo_stability
