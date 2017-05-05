module hf_class
!
!
!
!                      Hartree-Fock (HF) class module                                 
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017         
!                                                                           
!
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
!  Use general tools
!
   use types
   use utils
   use workspace
   use input_output
!
   implicit none
!
!  ::::::::::::::::::::::::::::::::::::
!  -::- Definition of the HF class -::-
!  ::::::::::::::::::::::::::::::::::::
!   
   type :: hartree_fock
!
!     Model name 
!
      character(len=7) :: name = 'HF     '
!
!     Orbital information attributes
!
      integer(i15) :: n_o  ! Number of occupied orbitals
      integer(i15) :: n_v  ! Number of virtual orbitals
      integer(i15) :: n_ao ! Number of atomic orbitals (AOs)
      integer(i15) :: n_mo ! Number of molecular orbitals (MOs)
      integer(i15) :: n_J  ! Number of Cholesky vectors
!
      real(dp), dimension(:,:), allocatable :: mo_coef ! MO coefficient matrix
!
!     Fock matrix attributes
!
      real(dp), dimension(:,:), allocatable :: fock_diagonal  ! diagonal vector
!
!     Energy attributes
!
      real(dp) :: energy            ! Same as scf_energy for HF class, different for descendants
!
      real(dp) :: nuclear_potential ! Nuclear potential energy term
      real(dp) :: scf_energy        ! The Hartree-Fock (HF/SCF) energy
!
   contains
!
!     Initialization and driver routines
!
      procedure :: init => init_hartree_fock
      procedure :: drv  => drv_hartree_fock
!
!     Routines to read MO Cholesky vectors from file
!
      procedure, non_overridable :: read_cholesky_ij => read_cholesky_ij_hartree_fock ! occ-occ
      procedure, non_overridable :: read_cholesky_ia => read_cholesky_ia_hartree_fock ! occ-vir
      procedure, non_overridable :: read_cholesky_ai => read_cholesky_ai_hartree_fock ! vir-occ
      procedure, non_overridable :: read_cholesky_ab => read_cholesky_ab_hartree_fock ! vir-vir
!
!     Routines needed to initialize HF     
!
!     read_hf_info            : sets attributes from file (n_o,n_v,scf_energy,...)
!     read_transform_cholesky : reads AO Cholesky vectors, transforms to MO basis, and
!                               saves the MO vectors to file
!
      procedure, non_overridable :: read_hf_info            => read_hf_info_hartree_fock
      procedure, non_overridable :: read_transform_cholesky => read_transform_cholesky_hartree_fock 
!
   end type hartree_fock
!
contains 
!
!  ::::::::::::::::::::::::::::::::::::::::::::
!  -::- Initialization and driver routines -::-
!  ::::::::::::::::::::::::::::::::::::::::::::
! 
   subroutine init_hartree_fock(wf)
!
!     Initialization of Hartree-Fock object
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Performs the following tasks:
!
!     1. Sets HF orbital and energy information by reading from file (read_hf_info)
!     2. Transforms AO Cholesky vectors to MO basis and saves to file (read_transform_cholesky)
!     3. Allocates the Fock matrix and sets it to zero (note: the matrix is constructed in 
!        the descendant classes) 
!
      implicit none
!
      class(hartree_fock) :: wf
!
!     Initialize HF attributes
!
      call wf%read_hf_info        
!
!     Initialize Cholesky vectors
!     
      call wf%read_transform_cholesky
!
   end subroutine init_hartree_fock
!
!
   subroutine drv_hartree_fock(wf)
!
!     HF Driver
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Lets the user know there is no driver for Hartree-Fock and exits
!     the program if called. The module reads Hartree-Fock information 
!     from files and contains no independent solver.
!
      implicit none 
!
      class(hartree_fock) :: wf
!
      write(unit_output,*) 'ERROR: There is no driver for the Hartree-Fock class'
      call exit
!
   end subroutine drv_hartree_fock
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::-
!  :::::::::::::::::::::::::::::::::::::::::
!
   subroutine read_hf_info_hartree_fock(wf)
!
!     Read HF Info
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Reads the file mlcc_hf_info and sets the following HF attributes: 
!     n_o, n_v, n_mo, orbital_coef, and fock_diagonal
!
!     The file mlcc_hf_info is written in the mlcc_write_sirifc 
!     subroutine, which is called from the wr_sirifc subroutine in
!     the siropt module of the DALTON suite
!  
      implicit none
!
      class(hartree_fock) :: wf
!
      integer(i15) :: unit_identifier_hf = -1 ! Unit identifier for mlcc_hf_info file
!
      integer(i15) :: n_lambda = 0 ! n_ao * n_mo, read but discarded
      integer(i15) :: i = 0, j = 0
!
!     Open the file mlcc_hf_info
!
      call generate_unit_identifier(unit_identifier_hf)
      open(unit=unit_identifier_hf, file='mlcc_hf_info', status='old', form='formatted')
      rewind(unit_identifier_hf)
!
!     Read mlcc_hf_info into HF variables
!  
      read(unit_identifier_hf,*) wf%n_mo, wf%n_o, n_lambda, &
                                  wf%nuclear_potential, wf%scf_energy
!
!     Set the energy equal to the read SCF energy
!
      wf%energy = wf%scf_energy
!
!     Calculate the number of virtuals
!
      wf%n_v = wf%n_mo - wf%n_o
!      
!     Allocate the Fock diagonal and the MO coefficients
!
      call allocator(wf%fock_diagonal, wf%n_mo, 1)
      wf%fock_diagonal = zero
!
      call allocator(wf%mo_coef, n_lambda, 1)
      wf%mo_coef = zero
!
!     Read in the Fock diagonal and MO coefficients
!
      read(unit_identifier_hf,*) (wf%fock_diagonal(i,1), i = 1, wf%n_mo)
      read(unit_identifier_hf,*) (wf%mo_coef(i,1), i = 1, n_lambda) 
!
!     Close the mlcc_hf_info file
!    
      close(unit_identifier_hf)
!
   end subroutine read_hf_info_hartree_fock
!
!
   subroutine read_transform_cholesky_hartree_fock(wf)
!
!     Read and Transform Cholesky
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, 20 Apr 2017
!
!     Reads the AO Cholesky vectors from file, transforms the vectors 
!     to the MO basis, and saves the MO vectors to file
!
      implicit none
!
      class(hartree_fock) :: wf
!
      integer(i15) :: unit_chol_ao    = -1 ! Unit identifier for mlcc_cholesky file
      integer(i15) :: unit_chol_mo_ij = -1 ! cholesky_ij file
      integer(i15) :: unit_chol_mo_ia = -1 ! cholesky_ia file
      integer(i15) :: unit_chol_mo_ab = -1 ! cholesky_ab file
!
      integer(i15) :: n_ao_sq_packed = 0 ! Packed dimensionality of (n_ao x n_ao) matrix
!
      real(dp), dimension(:,:), allocatable :: chol_ao    ! Packed AO Cholesky vector
      real(dp), dimension(:,:), allocatable :: chol_ao_sq ! Unpacked AO Cholesky vector
      real(dp), dimension(:,:), allocatable :: chol_mo_sq ! Unpacked MO Cholesky vector
!
      real(dp), dimension(:,:), allocatable :: X ! An intermediate matrix
!
      integer(i15) :: i,j,a,b,k 
!
!     Open Dalton file mlcc_cholesky (see mlcc_write_cholesky.F)
! 
      call generate_unit_identifier(unit_chol_ao)
      open(unit=unit_chol_ao, file='mlcc_cholesky', status='old', form='formatted')
      rewind(unit_chol_ao)
!
!     Read the number of Cholesky vectors (n_J) and 
!     the number of atomic orbitals (n_ao)
!
      read(unit_chol_ao,*) wf%n_ao, wf%n_J
!
!     Open files for MO Cholesky vectors 
! 
      call generate_unit_identifier(unit_chol_mo_ij)
      call generate_unit_identifier(unit_chol_mo_ia)
      call generate_unit_identifier(unit_chol_mo_ab)
!
      open(unit_chol_mo_ij, file='cholesky_ij', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ij)
!
      open(unit_chol_mo_ia, file='cholesky_ia', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ia)
!
      open(unit_chol_mo_ab, file='cholesky_ab', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ab)
!
!     Allocate packed and unpacked Cholesky AO, and 
!     unpacked Cholesky MO vectors
!
      n_ao_sq_packed = packed_size(wf%n_ao)
!
      call allocator(chol_ao, n_ao_sq_packed, 1)
      call allocator(chol_ao_sq, wf%n_ao, wf%n_ao) 
      call allocator(chol_mo_sq, wf%n_mo, wf%n_mo)
!
      chol_ao    = zero
      chol_ao_sq = zero
      chol_mo_sq = zero
!
!     Allocate an intermediate, X
!
      call allocator(X, wf%n_ao, wf%n_mo)
!
      X = zero
!
!     Loop over the number of Cholesky vectors,
!     reading them one by one 
!
      do j = 1, wf%n_J
!
!        Read Cholesky AO vector
!
         read(unit_chol_ao,*) (chol_ao(i,1), i = 1, n_ao_sq_packed)
!
!        Unpack/square up AO vector 
!
         call squareup(chol_ao, chol_ao_sq, wf%n_ao)
!
!        Transform the AO vectors to form the Cholesky MO vectors
!
         call dgemm('N','N',     &
                     wf%n_ao,    &
                     wf%n_mo,    &
                     wf%n_ao,    &
                     one,        &
                     chol_ao_sq, &
                     wf%n_ao,    &
                     wf%mo_coef, &
                     wf%n_ao,    &
                     zero,       &
                     X,          &
                     wf%n_ao)
!
         call dgemm('T','N',     &
                     wf%n_mo,    &
                     wf%n_mo,    &
                     wf%n_ao,    &
                     one,        &
                     wf%mo_coef, &
                     wf%n_ao,    &
                     X,          &
                     wf%n_ao,    &
                     zero,       &
                     chol_mo_sq, &
                     wf % n_mo)
!
!        Write the MO vectors to files in blocks
!
         write(unit_chol_mo_ij) ((chol_mo_sq(i,k), i = 1, wf%n_o), k = 1, wf%n_o)
         write(unit_chol_mo_ia) ((chol_mo_sq(i,a), i = 1, wf%n_o), a = wf%n_o + 1, wf%n_mo)
         write(unit_chol_mo_ab) ((chol_mo_sq(a,b), a = wf%n_o + 1, wf%n_mo), b = wf%n_o + 1, wf%n_mo)
!
      enddo
!
!     Close files 
!
      close(unit_chol_ao)
      close(unit_chol_mo_ij)
      close(unit_chol_mo_ia)
      close(unit_chol_mo_ab)
!  
   end subroutine read_transform_cholesky_hartree_fock
!
!
   subroutine read_cholesky_ij_hartree_fock(wf,L_ij_J)
!
!     Read Cholesky IJ vectors
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Reads the MO Cholesky IJ (occ-occ) vectors from file and 
!     places them in the incoming L_ij_J matrix
!
      implicit none
!
      class(hartree_fock) :: wf
!
      real(dp), dimension((wf%n_o)**2, wf%n_J) :: L_ij_J ! L_ij^J
!
      integer(i15) :: unit_chol_mo_ij = -1 ! Unit identifier for cholesky_ij file 
      integer(i15) :: i = 0, j = 0
!
!     Prepare for reading: generate unit idientifier, open file, and rewind
!
      call generate_unit_identifier(unit_chol_mo_ij)
      open(unit=unit_chol_mo_ij, file='cholesky_ij', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ij)
!
!     Read the Cholesky vectors into the L_ij_J matrix
!
      do j = 1, wf%n_J
         read(unit_chol_mo_ij) (L_ij_J(i,j), i = 1, (wf%n_o)**2)
      enddo
!
!     Close file
!
      close(unit_chol_mo_ij)    
!   
   end subroutine read_cholesky_ij_hartree_fock
!
!
   subroutine read_cholesky_ia_hartree_fock(wf,L_ia_J)
!
!     Read Cholesky IA 
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Reads the MO Cholesky IA (occ-vir) vectors from file and
!     places them in the incoming L_ia_J matrix
!
      implicit none
!
      class(hartree_fock) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v), wf%n_J) :: L_ia_J ! L_ia^J
!
      integer(i15) :: unit_chol_mo_ia = -1 ! Unit identifier for cholesky_ia file
      integer(i15) :: i = 0, j = 0
!
!     Prepare for reading: generate unit idientifier, open, and rewind file
!
      call generate_unit_identifier(unit_chol_mo_ia)
      open(unit=unit_chol_mo_ia, file='cholesky_ia', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ia)
!
!     Read Cholesky vectors into the L_ia_J matrix
!
      do j = 1, wf%n_J
         read(unit_chol_mo_ia) (L_ia_J(i,j), i = 1, (wf%n_o)*(wf%n_v))
      enddo
!
!     Close file
!
      close(unit_chol_mo_ia)    
!   
   end subroutine read_cholesky_ia_hartree_fock
!
!
   subroutine read_cholesky_ai_hartree_fock(wf, L_ai_J)
!
!     Read Cholesky AI 
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Reads the MO Cholesky AI (vir-occ) vectors from file and
!     places them in the incoming L_ai_J matrix
!
      implicit none
!
      class(hartree_fock) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), wf%n_J) :: L_ai_J ! L_ai^J
!
      real(dp), dimension(:,:), allocatable :: L_ia_J       
!
      integer(i15) :: i = 0, j = 0, a = 0, ia = 0, ai = 0
!
!     Allocation
!
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      L_ia_J = zero
!     
!     Get Cholesky IA vector 
!
      call wf%read_cholesky_ia(L_ia_J)
!
!     Reorder and save in AI vector 
!      
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
!           Needed indices
!
            ai = index_two(a, i, wf%n_v)
            ia = index_two(i, a, wf%n_o)
!
            do j = 1, wf%n_J
!
               L_ai_J(ai, j) = L_ia_J(ia, j)
!
            enddo
!
         enddo
      enddo
!
!     Deallocate temporary vector 
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)   
!
   end subroutine read_cholesky_ai_hartree_fock
!   
!
   subroutine read_cholesky_ab_hartree_fock(wf, L_ab_J, first, last, ab_dim, reorder)
!
!     Read Cholesky AB 
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Reads the MO Cholesky AB (vir-vir) vectors from file and
!     places them in the incoming L_ab_J matrix, with batching 
!     if necessary
!
!     If reorder = .true.,  L_ba_J is returned with batching over a
!     If reorder = .false., L_ab_J is returned with batching over b
!
      implicit none
!
      class(hartree_fock) :: wf
!
      integer(i15), intent(in) :: first   ! First index (can differ from 1 when batching)
      integer(i15), intent(in) :: last    ! Last index  (can differ from n_v when batching)
      integer(i15), intent(in) :: ab_dim  ! Dimension of ab index (not n_v^2 when batching)      
      logical, intent(in)      :: reorder ! See description above
!
      real(dp), dimension(ab_dim, wf%n_J) :: L_ab_J ! L_ab^J
!
      integer(i15) :: unit_chol_mo_ab = -1 ! Unit identifier for cholesky_ab file
!
      integer(i15) :: a = 0, b = 0, j = 0, i = 0
      integer(i15) :: batch_length = 0
!
      integer(i15) :: throw_away_index = 0
      real(dp)     :: throw_away = 0
!
!     Prepare for reading: generate unit identifier, open, and rewind file
!  
      call generate_unit_identifier(unit_chol_mo_ab)
      open(unit=unit_chol_mo_ab, file='cholesky_ab', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ab)
!
!     Calculating batch length
!
      batch_length = last - first + 1
!
      if (.not. reorder) then
!
         if (first .ne. 1) then
!  
!           Calculate index of last element to throw away
!  
            throw_away_index = index_two(wf%n_v, first - 1, wf%n_v)
!  
!           Throw away all elements from 1 to throw_away_index, then read from batch start
!  
            do j = 1, wf%n_J
!
              read(unit_chol_mo_ab) (throw_away, i = 1, throw_away_index),&
                                    (L_ab_J(a,j), a = 1, ab_dim)
!
            enddo
!
         else
!  
!           Read from the start of each entry
!  
            do j = 1, wf%n_J
!
              read(unit_chol_mo_ab) (L_ab_J(a,j), a = 1, ab_dim)
!
            enddo
!
         endif
!
      else ! Reorder L_ab_J is L_ba_J
!
         throw_away_index = index_two(wf%n_v, first - 1, wf%n_v)
!
!        Reading vectors
!
         do j = 1, wf%n_J
!
            if (first .eq. 1) then
! 
               read(unit_chol_mo_ab) ((L_ab_J(index_two(b, a, wf%n_v), j), b = 1, wf%n_v), a = 1, batch_length)
!
            else
!
               read(unit_chol_mo_ab) (throw_away, i = 1, throw_away_index),&
                                     ((L_ab_J(index_two(b, a, wf%n_v), j), b = 1, wf%n_v), a = 1, batch_length)
!
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
!
end module hf_class
