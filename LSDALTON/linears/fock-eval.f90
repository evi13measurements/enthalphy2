MODULE Fock_evaluator
  use matrix_operations

! This module evaluates the Fock/Kohn-Sham matrix using chosen
! algorithm, matrix representation, etc.
!
   PUBLIC

CONTAINS

   SUBROUTINE FCK_get_fock(d, F, Etotal)
      use dal_interface
      use lsdalton_fock_module
      IMPLICIT NONE
      TYPE(Matrix), intent(in)    :: D
      type(matrix), intent(inout) :: F  !output 
      real(realk), INTENT(OUT) :: Etotal

      !cfg_nfock = cfg_nfock + 1

      CALL di_get_fock_LSDALTON(D,lsint_fock_data%H1,F,Etotal, &
           & lsint_fock_data%lupri,lsint_fock_data%luerr,  &
           & lsint_fock_data%ls)
   END SUBROUTINE FCK_get_fock

   SUBROUTINE FCK_get_fock_LSINT(D,h1,F,Etotal,LSint,newlupri,newluerr,ls)
      use dal_interface
      use TYPEDEF
      IMPLICIT NONE
      TYPE(Matrix), intent(in)    :: D, h1
      type(matrix), intent(inout) :: F  !output 
      real(realk), INTENT(OUT) :: Etotal
      LOGICAL                  :: LSint
      INTEGER                  :: newlupri,newluerr
      type(lsitem) :: ls

      !cfg_nfock = cfg_nfock + 1
      CALL di_get_fock_LSDALTON(D,h1,F,Etotal,newlupri,newluerr,ls)
    END SUBROUTINE FCK_get_fock_LSINT

   subroutine fck_scale_virt_fock(S,H1,F,D)
     !to create a N-1 potential for the virtual orbitals
     !allowing for more 'occupied-like' virtual orbitals
     !good when configuration shifts are needed!!!
     !** F_mod = F + (n-1/n - 1)Q^T (F - H1) Q
     type(Matrix), intent(in) :: S,H1,D
     type(Matrix), intent(inout) :: F
     type(matrix) :: Q,wrk1,wrk2
     real(realk) :: konst
     integer :: ndim
     
     ndim = S%nrow
     call mat_init(Q,ndim,ndim)
     call mat_init(wrk1,ndim,ndim)
     call mat_init(wrk2,ndim,ndim)

     !** wrk2 = F - h1
     call mat_add(1.d0,F,-1.d0,H1,wrk2)
     !** Q = 1 - DS
     call mat_mul(D,S,'n','n',1.d0,0.d0,wrk1)
     call mat_add_identity(1.d0,-1.d0,wrk1,Q)
     !** wrk1 = (F - h1) Q = wrk2 Q
     call mat_mul(wrk2,Q,'n','n',1.d0,0.d0,wrk1)
     !** konst = -1/n
     konst = -1.d0/(2.d0*cfg_nocc)
     !** F_mod = F + konst Q^T (F - H1) Q = F + konst Q^T wrk1
     call mat_mul(Q,wrk1,'T','n',konst,1.d0,F)

     call mat_free(Q)
     call mat_free(wrk1)
     call mat_free(wrk2)
   end subroutine fck_scale_virt_fock

   ! If the virtual orbitals are scaled in the calculation
   ! then this routine unscales them such that the "actual" 
   ! Fock matrix is outputted.
   subroutine fck_unscale_virt(H1,S,D,F,nocc)
     use dal_interface
     ! F = F_mod - k Q^T G(D) Q = F_mod - k/(1+k)Q^T (F_mod - h1) Q  where
     ! k = (n-1)/n - 1 = -1/n
     implicit none
     type(matrix), intent(in)    :: H1,S,D
     type(Matrix),intent(inout)  :: F
     !> Number of occupied orbitals
     integer, intent(in) :: nocc
     type(Matrix) :: Q,wrk1,wrk2
     real(realk) :: konst,Etotal
     integer :: ndim

     if (.false.) then
     ndim = S%nrow
     call mat_init(Q,ndim,ndim)
     call mat_init(wrk1,ndim,ndim)
     call mat_init(wrk2,ndim,ndim)
     !wrk2 = F_mod - h1
     call mat_add(1.d0,F,-1.d0,H1,wrk2)
     !Q = 1 - DS
     call mat_mul(D,S,'n','n',1.d0,0.d0,wrk1)
     call mat_add_identity(1.d0,-1.0d0,wrk1,Q)
     !wrk1 = (F_mod - h1) Q = wrk2 Q
     call mat_mul(wrk2,Q,'n','n',1.d0,0.d0,wrk1)
     !F = F_mod + konst Q^T (F_mod - h1) Q = F_mod + konst Q^T wrk1
     konst = 1.d0/(nocc*2.d0 - 1.d0)
     call mat_mul(Q,wrk1,'t','n',konst,1.d0,F)
     call mat_free(Q)
     call mat_free(wrk1)
     call mat_free(wrk2)
     else
       !cfg_scale_virt = .false.
       CALL lsquit('replaced di_get_fock(d,F,Etotal) with this quit statement',-1)
     endif
   END SUBROUTINE fck_unscale_virt
  
END MODULE Fock_evaluator
