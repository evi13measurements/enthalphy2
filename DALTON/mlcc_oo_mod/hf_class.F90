module hf_class
!
!  Date:    April 2017
!  Authors: Sarai D. Folkestad and Eirik F. Kjønstad
   use workspace
   use mlcc_types
   use input_output
   use mlcc_oo_utilities
!
   implicit none
!
   type hartree_fock
!
      integer(i15) :: n_occ
      integer(i15) :: n_vir
      integer(i15) :: n_mo 
      integer(i15) :: n_J
      integer(i15) :: n_ao
!
      real(dp), dimension(:,:), allocatable :: mo_coef
      real(dp), dimension(:,:), allocatable :: fock_diagonal
!
      real(dp), dimension(:,:), allocatable :: fock_matrix_ij
      real(dp), dimension(:,:), allocatable :: fock_matrix_ia
      real(dp), dimension(:,:), allocatable :: fock_matrix_ai
      real(dp), dimension(:,:), allocatable :: fock_matrix_ab
!
      real(dp) :: nuclear_potential
      real(dp) :: scf_energy
!
   contains
!
!     Public routines
!
      procedure   :: init => init_hartree_fock
!
!     Private routines
!
      procedure  :: read_cholesky_ij => read_cholesky_ij_hartree_fock
      procedure  :: read_cholesky_ia => read_cholesky_ia_hartree_fock
      procedure  :: read_cholesky_ai => read_cholesky_ai_hartree_fock
      procedure  :: read_cholesky_ab => read_cholesky_ab_hartree_fock
!
      procedure :: read_hf_info              => read_hf_info_hartree_fock
      procedure :: read_transform_cholesky   => read_transform_cholesky_hartree_fock 
      procedure :: allocate_fock_matrix      => allocate_fock_matrix_hartree_fock
!
   end type hartree_fock
!
contains 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!                                          !!!!
!!!!  Public class subroutines and functions  !!!!
!!!!                                          !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine init_hartree_fock(wfn)
!
!  Initialization of Hartree-Fock object
!
!  Calls read_hf_info - reads MLCC_HF_INFO
!  Initializes Cholesky vectors  
!  Allocates Fock matrix and sets it to 0. Fock matrix is constructed in derived types.
!  
      implicit none
!
      class(hartree_fock) :: wfn
!
      write(unit_output,*) 'In init_hartree_fock'
      call flshfo(unit_output)
!
!     Initializing HF variables by reading MLCC_HF_INFO
!
      write(unit_output,*) 'Read and HF info...'
      call flshfo(unit_output)
      call wfn % read_hf_info        
!
!     Initializing integral and Cholesky variables
!     
      write(unit_output,*) 'Read and transform Cholesky vectors...'
      call flshfo(unit_output)      
      call wfn % read_transform_cholesky
!
!     Allocate Fock matrix and set to 0
!
      write(unit_output,*) 'Allocate Fock matrix blocks...'
      call flshfo(unit_output)
      call wfn % allocate_fock_matrix
      

   end subroutine init_hartree_fock
!
   subroutine read_hf_info_hartree_fock(wfn)
!
!  Date:    April 2017
!  Authors: Sarai D. Folkestad and Eirik F. Kjønstad
!
!  Reads MLCC_HF_INFO and initializes HF variables n_occ, n_vir, n_mo, orbital_coef, fock_diagonal
!
!  MLCC_HF_INFO is written in the mlcc_write_sirifc subroutine 
!  and called from wr_sirifc subroutine in siropt module.
!  
      use workspace
!
      implicit none
!
      class(hartree_fock) :: wfn
!
      integer(i15)        :: unit_identifier_hf = -1 
      integer(i15)        :: n_lambda
      integer(i15)        :: i,j
!
!     Open mlcc_hf_info
!     ---------------------
!
      call generate_unit_identifier(unit_identifier_hf)
      open(unit=unit_identifier_hf,file='mlcc_hf_info',status='old',form='formatted')
      rewind(unit_identifier_hf)
!
!     Read mlcc_hf_info
!     ---------------------
!  
      read(unit_identifier_hf,*) wfn % n_mo, wfn % n_occ, n_lambda, wfn % nuclear_potential, wfn % scf_energy
!
!     Setting n_vir
!
      wfn % n_vir = (wfn % n_mo) - (wfn % n_occ)

!      
!     Allocate space for Fock diagonal and coefficients.
!
      call allocator(wfn % fock_diagonal, wfn % n_mo,1)
      wfn % fock_diagonal = zero
!
      call allocator(wfn % mo_coef,n_lambda,1)
      wfn % mo_coef = zero
!
!     Read in Fock diagonal and coefficients
!
      read(unit_identifier_hf,*) (wfn % fock_diagonal(i,1),i=1,wfn % n_mo)
      read(unit_identifier_hf,*) (wfn % mo_coef(i,1),i=1,n_lambda)   
!
!     Done with file
!     ------------------
!    
      close(unit_identifier_hf)
!
   end subroutine read_hf_info_hartree_fock
!
  subroutine read_transform_cholesky_hartree_fock(wfn)
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
      class(hartree_fock) :: wfn
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
      read(unit_chol_ao,*) wfn%n_ao, wfn%n_J
!
!     Open files for MO Cholesky vectors 
! 
      call generate_unit_identifier(unit_chol_mo_ij)
      call generate_unit_identifier(unit_chol_mo_ia)
      call generate_unit_identifier(unit_chol_mo_ab)
!
      open(unit_chol_mo_ij,file='cholesky_ij',status='unknown',form='unformatted')
      rewind(unit_chol_mo_ij)
!
      open(unit_chol_mo_ia,file='cholesky_ia',status='unknown',form='unformatted')
      rewind(unit_chol_mo_ia)
!
      open(unit_chol_mo_ab,file='cholesky_ab',status='unknown',form='unformatted')
      rewind(unit_chol_mo_ab)
!
!     Allocate packed and unpacked Cholesky AO, and 
!     unpacked Cholesky MO vectors
!
      n_ao_sq_packed = packed_size(wfn%n_ao)
!
      call allocator(chol_ao,n_ao_sq_packed,1)
      call allocator(chol_ao_sq,wfn%n_ao,wfn%n_ao) 
      call allocator(chol_mo_sq,wfn%n_mo,wfn%n_mo)
!
      chol_ao    = zero
      chol_ao_sq = zero
      chol_mo_sq = zero
!
!     Allocate an intermediate, X
!
      call allocator(X,wfn%n_ao,wfn%n_mo)
!
      X = zero
!
!     Loop over the number of Cholesky vectors,
!     reading them one by one 
!
      do j = 1,wfn%n_J
!
!        Read Cholesky AO vector
!
         read(unit_chol_ao,*) (chol_ao(i,1),i=1,n_ao_sq_packed)
!
!        Unpack/square up AO vector 
!
         call squareup(chol_ao,chol_ao_sq,wfn%n_ao)
!
!        Transform the AO vectors to form the Cholesky MO vectors
!
         call dgemm('N','N',      &
                     wfn%n_ao,    &
                     wfn%n_mo,    &
                     wfn%n_ao,    &
                     one,         &
                     chol_ao_sq,  &
                     wfn%n_ao,    &
                     wfn%mo_coef, &
                     wfn%n_ao,    &
                     zero,        &
                     X,           &
                     wfn%n_ao)
!
         call dgemm('T','N',      &
                     wfn%n_mo,    &
                     wfn%n_mo,    &
                     wfn%n_ao,    &
                     one,         &
                     wfn%mo_coef, &
                     wfn%n_ao,    &
                     X,           &
                     wfn%n_ao,    &
                     zero,        &
                     chol_mo_sq,  &
                     wfn%n_mo)
!
!        Write the MO vectors to files in blocks
!
         write(unit_chol_mo_ij) ((chol_mo_sq(i,j),i=1,wfn%n_occ),k=1,wfn%n_occ)
         write(unit_chol_mo_ia) ((chol_mo_sq(i,a),i=1,wfn%n_occ),a=wfn%n_occ+1,wfn%n_mo)
         write(unit_chol_mo_ab) ((chol_mo_sq(a,b),a=wfn%n_occ+1,wfn%n_mo),b=wfn%n_occ+1,wfn%n_mo)
!
      enddo
!
!     Close files 
!
      close(unit_chol_mo_ij)
      close(unit_chol_mo_ia)
      close(unit_chol_mo_ab)
!  
   end subroutine read_transform_cholesky_hartree_fock
!
   subroutine allocate_fock_matrix_hartree_fock(wfn)
!
      use workspace
!
      implicit none
!  
      class(hartree_fock) :: wfn   
!
      call allocator(wfn % fock_matrix_ij, wfn % n_occ, wfn % n_occ)
      call allocator(wfn % fock_matrix_ia, wfn % n_occ, wfn % n_vir)
      call allocator(wfn % fock_matrix_ai, wfn % n_vir, wfn % n_occ)
      call allocator(wfn % fock_matrix_ab, wfn % n_vir, wfn % n_vir)

!
      wfn % fock_matrix_ij = zero
      wfn % fock_matrix_ia = zero
      wfn % fock_matrix_ai = zero
      wfn % fock_matrix_ab = zero
!
   end subroutine allocate_fock_matrix_hartree_fock
!
!
   subroutine read_cholesky_ij_hartree_fock(wfn,L_ij_J)
!
!     Read mo ij cholesky vectors from file and places them in L_ij_J array
!
      implicit none
!
      class(hartree_fock)      :: wfn
      real(dp), dimension(:,:) :: L_ij_J
!
      integer(i15) :: unit_chol_mo_ij =-1
      integer(i15) :: i=0,j=0
!
!     Prepare for reading. Generate unit idientifier, open file and rewind.
!
      call generate_unit_identifier(unit_chol_mo_ij)
      
      open(unit=unit_chol_mo_ij,file='cholesky_ij',status='unknown',form='unformatted')
      rewind(unit_chol_mo_ij)
!
!     Read vectors
!
      do j = 1,wfn % n_J
         read(unit_chol_mo_ij) (L_ij_J(i,j), i=1,(wfn % n_occ)*(wfn % n_occ))
      enddo
!
!     Close file
!
      close(unit_chol_mo_ij)    
!   
   end subroutine read_cholesky_ij_hartree_fock
!
   subroutine read_cholesky_ia_hartree_fock(wfn,L_ia_J)
!
!     Read mo ia cholesky vectors from file and places them in L_ia_J array
!
      implicit none
!
      class(hartree_fock)      :: wfn
      real(dp), dimension(:,:) :: L_ia_J
!
      integer(i15) :: unit_chol_mo_ia =-1
      integer(i15) :: i=0,j=0
!
!     Prepare for reading. Generate unit idientifier, open and rewind file.
!
      call generate_unit_identifier(unit_chol_mo_ia)
      
      open(unit=unit_chol_mo_ia,file='cholesky_ia',status='unknown',form='unformatted')
      rewind(unit_chol_mo_ia)
!
!     Read vectors
!
      do j = 1,wfn % n_J
         read(unit_chol_mo_ia) (L_ia_J(i,j), i=1,(wfn % n_occ)*(wfn % n_vir))
      enddo
!
!     Close file
!
      close(unit_chol_mo_ia)    
!   
   end subroutine read_cholesky_ia_hartree_fock
!
!
   subroutine read_cholesky_ai_hartree_fock(wfn,L_ai_J)
!
!     Read mo ai cholesky vectors from file and places them in L_ai_J array
!
      implicit none
!
      class(hartree_fock)      :: wfn
      real(dp), dimension(:,:) :: L_ai_J
!
      integer(i15) :: unit_chol_mo_ai =-1
      integer(i15) :: i=0,j=0
!
!     Prepare for reading. Generate unit idientifier, open and rewind file.
!
      call generate_unit_identifier(unit_chol_mo_ai)
      
      open(unit=unit_chol_mo_ai,file='cholesky_ai',status='unknown',form='unformatted')
      rewind(unit_chol_mo_ai)
!
      do j = 1,wfn % n_J
         read(unit_chol_mo_ai) (L_ai_J(i,j), i=1,(wfn % n_occ)*(wfn % n_vir))
      enddo
!
!     Close file.
!
      close(unit_chol_mo_ai)    
!
   end subroutine read_cholesky_ai_hartree_fock
!   
!
   subroutine read_cholesky_ab_hartree_fock(wfn,L_ab_J,start,end,ab_dim,reorder)
!
!     Read mo ab cholesky vectors, with batching if needed. 
!
!        If reorder = .true. L_ba_J is returned with batching over a,
!                      otherwize L_ab_J is returned with batching over b.
!
!
      implicit none
!
      class(hartree_fock)      :: wfn
      real(dp), dimension(:,:) :: L_ab_J
!
      integer(i15) :: unit_chol_mo_ab=-1
      integer(i15) :: ab_dim
      integer(i15) :: a=0,b=0,j=0,i=0
      integer(i15) :: start,end
      integer(i15) :: batch_length=0
      integer(i15) :: throw_away_index
      real(dp)     :: throw_away
!
      logical :: reorder
!
!     Prepare for reading. Generate unit identifier, open and rewind file
!  
      call generate_unit_identifier(unit_chol_mo_ab)
      open(unit=unit_chol_mo_ab,file='cholesky_ab',status='unknown',form='unformatted')
      rewind(unit_chol_mo_ab)
!
!     Calculating batch length
!
      batch_length = end-start+1
!
      if (.not. reorder) then
!
!  
         if (start .ne. 1) then
!  
!           Calculate index of last element to throw away
!  
            throw_away_index=index_two(wfn%n_vir,start-1,wfn%n_vir)
!  
!           Throw away all elements from 1 to throw_away_index, then read from batch start.
!  
            do j = 1,wfn%n_J
              read(unit_chol_mo_ab) (throw_away,i=1,throw_away_index),(L_ab_J(a,j),a=1,ab_dim)
            enddo
!
         else
!  
!           Read from start
!  
            do j = 1,wfn%n_J
              read(unit_chol_mo_ab)(L_ab_J(a,j),a=1,ab_dim)
            enddo
!
         endif
!
      else ! Reorder L_ab_J is L_ba_J
!
         throw_away_index = index_two(wfn%n_vir,start-1,wfn%n_vir)
!
!        Reading vectors
!
         do j = 1,wfn%n_J
!
            if (start .eq. 1) then 
               read(unit_chol_mo_ab)((L_ab_J(index_two(b,a,wfn%n_vir),j),b=1,wfn%n_vir),a=1,batch_length)
            else
               read(unit_chol_mo_ab)(throw_away,i=1,throw_away_index), &
                                    ((L_ab_J(index_two(b,a,wfn%n_vir),j),b=1,wfn%n_vir),a=1,batch_length)
            endif
!
         enddo
!     
      endif  ! Reorder
! 
!     Close file
!    
      close(unit_chol_mo_ab)
!
   end subroutine read_cholesky_ab_hartree_fock
!
end module hf_class
