module cholesky_integrals_class
!
!  Cholesky integrals class
!  Written by Sarai D. Folkestad and Eirik F. Kjønstad, 21 Apr 2017
!
   use input_output
   use mlcc_types
!
   type cholesky_integrals
!
      integer(i15) :: n_J  = 0 ! Number of Cholesky vectors 
      integer(i15) :: n_ao = 0 ! Number of atomic orbitals
!
   contains
!	
      procedure :: init => init_cholesky_integrals
      procedure :: read_ia => read_ia_cholesky_integrals
!
   end type cholesky_integrals
!
!  Private variables and routine
!
      private :: read_and_transform_cholesky_vectors
!
contains
!
   subroutine init_cholesky_integrals(chol,mo_coef,n_occ,n_vir)
!
      implicit none
!
      class(cholesky_integrals) :: chol
!
      integer(i15) :: n_occ,n_vir 
      real(dp), dimension(n_occ+n_vir,n_occ+n_vir) :: mo_coef
!
      write(unit_output,*) 'In init_cholesky_integrals'
      call flshfo(unit_output)
!
!     Read the Cholesky vectors from file (in the AO basis),
!     transform them to the MO basis, and save to file
!
      call read_and_transform_cholesky_vectors(chol,mo_coef,n_occ,n_vir)
!
   end subroutine init_cholesky_integrals
!
   subroutine read_ia_cholesky_integrals(chol,L_ia_J,n_ov)
!
!     Read IA Cholesky Integrals 
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, 25 Apr 2017
!
!     Reads the Cholesky vector L_ia^J from file
!     and places it in the incoming vector L_ia_J
!
      implicit none 
!
      class(cholesky_integrals) :: chol 
!
      integer(i15) :: n_ov ! n_occ * n_vir
      real(dp), dimension(n_ov,chol%n_J) :: L_ia_J
!
   end subroutine read_ia_cholesky_integrals
!
   subroutine read_and_transform_cholesky_vectors(chol,mo_coef,n_occ,n_vir)
!
!     Read and Transform Cholesky Vectors
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, 20 Apr 2017
!
!     Reads the AO Cholesky vectors from file, transforms the vectors 
!     to the MO basis, and saves the MO vectors to file
!
      use input_output
      use workspace
      use mlcc_oo_utilities
!
      implicit none
!
      class(cholesky_integrals) :: chol
!
      integer(i15) :: n_occ,n_vir
!
      real(dp), dimension(n_occ+n_vir,n_occ+n_vir) :: mo_coef
!
      integer(i15) :: n_mo
!
      integer(i15) :: unit_chol_ao    = -1
      integer(i15) :: unit_chol_mo_ij = -1
      integer(i15) :: unit_chol_mo_ia = -1
      integer(i15) :: unit_chol_mo_ab = -1
!
      integer(i15) :: n_ao_sq_packed      = 0
!
      real(dp), dimension(:,:), allocatable :: chol_ao
      real(dp), dimension(:,:), allocatable :: chol_ao_sq
      real(dp), dimension(:,:), allocatable :: chol_mo_sq
!
      real(dp), dimension(:,:), allocatable :: X
!
      integer(i15) :: i,j,a,b,k 
!
!     Open Dalton file MLCC_CHOLESKY (see mlcc_write_cholesky.F)
! 
      call generate_unit_identifier(unit_chol_ao)
      open(unit=unit_chol_ao,file='mlcc_cholesky',status='old',form='formatted')
      rewind(unit_chol_ao)
!
!     Read the number of Cholesky vectors (n_J) and 
!     the number of atomic orbitals (n_ao)
!
      read(unit_chol_ao,*) chol%n_ao, chol%n_J
!
!     Open files for MO Cholesky vectors 
! 
      call generate_unit_identifier(unit_chol_mo_ij)
      call generate_unit_identifier(unit_chol_mo_ia)
      call generate_unit_identifier(unit_chol_mo_ab)
!
      open(unit_chol_mo_ij,file='cholesky_ij',status='new',form='unformatted')
      rewind(unit_chol_mo_ij)
!
      open(unit_chol_mo_ia,file='cholesky_ia',status='new',form='unformatted')
      rewind(unit_chol_mo_ia)
!
      open(unit_chol_mo_ab,file='cholesky_ab',status='new',form='unformatted')
      rewind(unit_chol_mo_ab)
!
!     Allocate packed and unpacked Cholesky AO, and 
!     unpacked Cholesky MO vectors
!
      n_ao_sq_packed = packed_size(chol%n_ao)
      n_mo           = n_occ + n_vir
!
      call allocator(chol_ao,n_ao_sq_packed,1)
      call allocator(chol_ao_sq,chol%n_ao,chol%n_ao) 
      call allocator(chol_mo_sq,n_mo,n_mo)
!
      chol_ao    = zero
      chol_ao_sq = zero
      chol_mo_sq = zero
!
!     Allocate an intermediate, X
!
      call allocator(X,chol%n_ao,n_mo)
!
      X = zero
!
!     Loop over the number of Cholesky vectors,
!     reading them one by one 
!
      do j = 1,chol%n_J
!
!        Read Cholesky AO vector
!
         read(unit_chol_ao,*) (chol_ao(i,1),i=1,n_ao_sq_packed)
!
!        Unpack/square up AO vector 
!
         call squareup(chol_ao,chol_ao_sq,chol%n_ao)
!
!        Transform the AO vectors to form the Cholesky MO vectors
!
         call dgemm('N','N',     &
                     chol%n_ao,  &
                     n_mo,       &
                     chol%n_ao,  &
                     one,        &
                     chol_ao_sq, &
                     chol%n_ao,  &
                     mo_coef,    &
                     chol%n_ao,  &
                     zero,       &
                     X,          &
                     chol%n_ao)
!
         call dgemm('T','N',     &
                     n_mo,       &
                     n_mo,       &
                     chol%n_ao,  &
                     one,        &
                     mo_coef,    &
                     chol%n_ao,  &
                     X,          &
                     chol%n_ao,  &
                     zero,       &
                     chol_mo_sq, &
                     n_mo)
!
!        Write the MO vectors to files in blocks
!
         write(unit_chol_mo_ij) ((chol_mo_sq(i,j),i=1,n_occ),k=1,n_occ)
         write(unit_chol_mo_ia) ((chol_mo_sq(i,a),i=1,n_occ),a=n_occ+1,n_mo)
         write(unit_chol_mo_ab) ((chol_mo_sq(a,b),a=n_occ+1,n_mo),b=n_occ+1,n_mo)
!
      enddo
!
!     Close files 
!
      close(unit_chol_mo_ij)
      close(unit_chol_mo_ia)
      close(unit_chol_mo_ab)
!  
   end subroutine read_and_transform_cholesky_vectors
!
end module cholesky_integrals_class
