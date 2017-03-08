!Interface routines to new MLCC code
!Requires a call to MLCC_input, mlcc_driver, mlcc_HF_reader and mlcc_ao_reader
!Dalton input is a mess, where is the best way to call the routine SIR_NEWINP?
!
!
subroutine mlcc_interface_input(lunit,word,lupri)
!
!  Read in input from  DALTON.INP, starting with *MLCC until next * and store in input string
!  Pass on input to mlcc_input in mlcc_module
!
!  Rolf H. Myhre 
!  December 2016
!
   use mlcc_input_data
   use mlcc_workspace
!
!
   implicit none
!
   integer, intent(in)              :: lunit !input unit
   integer, intent(in)              :: lupri !general output unit
!
   character(len=7), intent(inout)  :: word
   character(len=24)                :: buffer
!
   integer                          :: r_off, length, last, actual
   integer                          :: i_count,i, r_len
!
   r_off=1
!
   if (word(1:7) .eq. '*MLCC  ') then
!
      do
!
!        Read next keyword      
         read (lunit,'(a)') buffer
!
!        Skip comments
         do while (buffer(1:1) .eq. '!' .or. buffer(1:1) .eq. '#' )
            read (lunit,'(a)') buffer
         end do
!            
!        Stop when hitting the next *
         if (buffer(1:1) .eq.'*') then
            word = buffer(1:7)
            exit
         end if
!
!        
         select case (trim(buffer))
!        
            case('.MLACTIVE')
!        
               mlcc_active=.true.
!        
            case('.PRINT')
!        
               read (lunit,*) print_mlcc
!        
            case default
               write(lupri,*) 'Keyword ', trim(buffer), ' not recognized in mlcc_input'
               stop
!        
         end select
!
!
      end do
!
   else 
      write(lupri,*) 'mlcc_input called without *MLCC'
      call flshfo(lupri)
      call quit('mlcc_input, something is wrong')
   end if
!
   if(print_mlcc .ge. 3) then
!
      write(lupri,*)
      write(lupri,*) 'Output from mlcc3_input'
      write(lupri,*) 'print_mlcc', print_mlcc
      write(lupri,*) 'mlcc_active', mlcc_active
      write(lupri,*)
!
   end if
!
end subroutine mlcc_interface_input
!
!
subroutine mlcc_interface_drv(work,lwork,lupri)
!
!  Interface routine for the MLCC driver. 
!
!  Rolf H. Myhre 
!  December 2016
!
   use mlcc_drive, only: mlcc_drv
!
   implicit none
!
   integer, intent(in)                    :: lupri !general output unit
   integer, intent(in)                    :: lwork !free space in work
   real*8, intent(in), dimension(lwork)   :: work !work statis array
!
   call mlcc_drv(work,lwork,lupri)
!
end subroutine mlcc_interface_drv
!
subroutine mlcc_iajb(vec)
!
!  MLCC iajb MO integral calculator
!  Authors Sarai Folkestad, Eirik Kjønstad, January 2017
!
!  Purpose: read in atomic orbitals and calculate iajb MO integrals
!
   use mlcc_data
   use mlcc_types
!
   implicit none
!
   real(dp), dimension(n_t2am_pack,1) :: vec
!
   vec = zero
!
end subroutine mlcc_iajb
!
subroutine mlcc_get_cholesky
!
!  MLCC Cholesky vector reader, and transformator
!  Authors:  Eirik Kjønstad and Sarai Folkestad, January 2017
!
!  Purpose: Read Cholesky vectors in AO basis, transform to MO basis and write to files CHOLESKY_IJ,
!           CHOLESKY_IA, CHOLESKY_AI and CHOLESKY_AB according to PQ-type.
!
   use mlcc_data
   use mlcc_workspace
   use mlcc_utilities
!
   implicit none
!
   integer       :: ludiag,lucho,lucho_ij,lucho_ia,lucho_ab,lumlch
   integer       :: i,j,k,idum,a,b
   integer       :: n_twoel_diag
   integer, dimension(:,:), pointer       :: index_reduced  => null()
   real(dp), dimension(:,:), pointer      :: cho_diag       => null()  
   real(dp), dimension(:,:), pointer      :: cho_ao         => null()
   integer       :: lencho
   real(dp), dimension(:,:), pointer      :: cholesky_ao,cholesky_ao_sq    => null()   
   real(dp), dimension(:,:), pointer      :: cholesky_mo_sq,cholesky_mo,X  => null()
!
   n_twoel_diag = packed_size(n_basis)
!
! Open MLCC_CHOLESKY - See mlcc_write_cholesky.F
!  
   lumlch = -1
   call gpopen(lumlch,'MLCC_CHOLESKY','OLD','SEQUENTIAL','FORMATTED',idum,.false.)
   rewind(lumlch)
!
! Read number (n_J) and length (n_reduced) of Cholesky vectors
!
   read(lumlch,*)n_reduced,n_J
!
! Read reduced index array
!
   call allocator_int(index_reduced,n_reduced,1)
   index_reduced=0
!
   read(lumlch,*)(index_reduced(i,1),i=1,n_reduced)
!
!========================================
!
! Preparation of files for cholesky in MO basis:
! ij-type:
  lucho_ij = -1
  call gpopen(lucho_ij,'CHOLESKY_IJ','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
  rewind(lucho_ij)
! ai-type:
  lucho_ia = -1
  call gpopen(lucho_ia,'CHOLESKY_IA','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
  rewind(lucho_ia)
! ab-type:
  lucho_ab = -1
  call gpopen(lucho_ab,'CHOLESKY_AB','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
  rewind(lucho_ab)
!
!
!
! Allocation
! Cholesky AO array, intermediate X and cholesky MO array
!
   call allocator(cholesky_ao,n_twoel_diag,1)
   call allocator(cholesky_ao_sq,n_basis,n_basis)
   call allocator(X,n_basis,n_orbitals)
   call allocator(cholesky_mo_sq,n_orbitals,n_orbitals)
!
   cholesky_ao = zero
   cholesky_ao_sq = zero
   cholesky_mo_sq = zero
   X = zero
   do j=1,n_J
!
! Read Cholesky AO:
!
      read(lumlch,*) (cholesky_ao(i,1),i=1,n_twoel_diag)
! Square up of cholesky_ao:
!
      call squareup(cholesky_ao,cholesky_ao_sq,n_basis)
!
! Transform to MO
!  
      call dgemm('N','N',n_basis,n_orbitals,n_basis,one,cholesky_ao_sq,n_basis,orb_coefficients,n_basis,zero,X,n_basis)
      call dgemm('T','N',n_orbitals,n_orbitals,n_basis,one,orb_coefficients,n_basis,X,n_basis,zero,cholesky_mo_sq,n_orbitals)
!
! Packing MO Cholesky and save to file
!
      write(lucho_ij)((cholesky_mo_sq(i,k),i=1,n_occ),k=1,n_occ)
      write(lucho_ia)((cholesky_mo_sq(i,a),i=1,n_occ),a=n_occ+1,n_orbitals)
      write(lucho_ab)((cholesky_mo_sq(a,b),a=n_occ+1,n_orbitals),b=n_occ+1,n_orbitals)
!
  enddo
!
!  Close all files
!
   call gpclose(lumlch,'KEEP')
   call gpclose(lucho_ij,'KEEP')
   call gpclose(lucho_ia,'KEEP')
   call gpclose(lucho_ab,'KEEP')
!
!  Deallocation
!
   call deallocator(cholesky_ao,n_twoel_diag,1)
   call deallocator(cholesky_ao,n_basis,n_basis)
   call deallocator_int(index_reduced,n_reduced,1)
   call deallocator(cholesky_mo_sq,n_orbitals,n_orbitals)
   call deallocator(X,n_basis,n_orbitals)
!
!================================
end subroutine mlcc_get_cholesky
!
subroutine hf_reader

   use mlcc_data
   use mlcc_workspace

!
!  Hartree-Fock reader routine
!  Authors Henrik Koch, Rolf H. Myhre, Sarai Folkestad, Eirik Kjønstad
!  January 2017
!
!  Purpose: read in data from LUSIFC file
!
      use mlcc_data
      use mlcc_workspace
!
      implicit none
!
!
      integer  :: lusifc = -1 
      integer  :: idummy = 0
      integer  :: n_symmetries, n_basis_sym, n_orbitals_sym
      integer  :: i,j
!      
!     
!     Open Sirius Fock file
!     ---------------------
!
      call gpopen(lusifc,'SIRIFC','OLD',' ','UNFORMATTED',idummy,'.FALSE.')
      rewind(lusifc)
!
!     Read in various stuff from Sirius Fock file. Things depending on symmetry is mostly
!     discarded at the end of the subroutine as we do not use symmetry. Information in 
!     file should be in Cholesky orbital format. If Cholesky orbitals has been generated,
!     SIRIFC will contain the data in the Cholesky basis
!
      call mollab('TRCCINT ',lusifc,luprint)
!      
      read(lusifc) n_symmetries, n_orbitals, n_basis, n_lambda, n_occ, &
      &            n_orbitals_sym, n_basis_sym, nuclear_potential, scf_energy
!
!      
      if (n_symmetries /= 1) call quit('error in mlcc_mod_init: not implemented with symmetry')
!      
!      
!     Calculate number of virtuals and amplitudes
!
      n_vir          = n_orbitals - n_occ
      n_ov           = n_occ*n_vir
      n_oo           = n_occ*n_occ
      n_vv           = n_vir*n_vir
      n_oo_packed    = n_occ*(n_occ+1)/2
      n_oov          = n_vir*n_occ*n_occ
      n_basis_2_pack = n_basis*(n_basis+1)/2
      n_ooo          = n_occ*n_occ*n_occ
      n_vv_packed    = n_vir*(n_vir+1)/2
      n_ov_ov_packed = n_ov*(n_ov+1)/2
      n_ovv          = n_occ*n_vir*n_vir
!      
!     Allocate space for Fock diagonal and coefficients.
!
      call allocator(fock_diagonal,n_orbitals,1)
      fock_diagonal = zero
!      
      call allocator(orb_coefficients,n_lambda,1)
      orb_coefficients = zero
!
!     Read in Fock diagonal and coefficients
!
      read(lusifc) (fock_diagonal(i,1),i=1,n_orbitals)
      read(lusifc) (orb_coefficients(i,1),i=1,n_lambda) 
!
!     Done with file
!
      call gpclose(lusifc,'KEEP')
!

   end subroutine hf_reader
!
   subroutine mlcc_fock()
!  Purpose: Construct Fock matrix in MO basis
!
!  One-electron part obtained by reading one-electron integrals from MLCC_AOINT.
!  MLCC_AOINT is written to file from the mlcc_wr_aoint routine in sirius/mlcc_wr_aoint.F
!
!  Two-electron contributions calculated using Cholesky vectors
      use mlcc_data
      use mlcc_workspace
      use mlcc_utilities
      use mlcc_cholesky
!
      implicit none
      real(dp), dimension(:,:), pointer      :: fock_ao  => null()
      real(dp), dimension(:,:), pointer      :: ao_int  => null()
      real(dp), dimension(:,:), pointer      :: X => null()
      integer                                :: luaoin = -1
      integer                                :: idummy,i,j,k,a,b,ij,kk,ik,kj,ii,jj,ji,ai,ib,bi,ia,aj,ja,ab
      real(dp), dimension(:,:), pointer      :: g_ij_kl => null()
      real(dp), dimension(:,:), pointer      :: g_ab_ij => null()
      real(dp), dimension(:,:), pointer      :: g_ai_jb => null()
      real(dp), dimension(:,:), pointer      :: g_ia_jk => null()
      real(dp), dimension(:,:), pointer      :: g_ai_jk => null()
      real(dp), dimension(:,:), pointer      :: L_ij_J => null()
      real(dp), dimension(:,:), pointer      :: L_ia_J => null()
      real(dp), dimension(:,:), pointer      :: L_ai_J => null()
      real(dp), dimension(:,:), pointer      :: L_ab_J => null()
      integer                                :: available, required,max_batch_length,n_batch,batch_start
      integer                                :: batch_end,batch_length,g_off
      integer                                :: b_batch = 0
      logical                                :: debug = .true.
!
   call allocator(mo_fock_mat,n_orbitals,n_orbitals)
!
!
!
!!! ONE-ELECTRON CONTRIBUTION !!!
!
!
!  Allocate for one-electron ao integrals
!
      call allocator(ao_int,n_basis_2_pack,1)
!
!  Read one-electron ao integrals
!
      call gpopen(luaoin,'MLCC_AOINT','UNKNOWN','SEQUENTIAL','FORMATTED',idummy,.false.)
      rewind(luaoin)
!
      read(luaoin,*)(ao_int(i,1),i=1,n_basis_2_pack)
!
      call gpclose(luaoin,'KEEP')
!
!  Allocate ao fock matrix, add one-electron contributions
!
      call allocator(fock_ao,n_basis,n_basis)
      call squareup(ao_int,fock_ao,n_basis)   
!
!  Deallocation of ao integrals
!   
      call deallocator(ao_int,n_basis_2_pack,1) 
!
!  Transform to one-electron part to mo and add to mo_fock_mat 
!
      call allocator(X,n_basis,n_orbitals)
!
      call dgemm('N','N',n_basis,n_orbitals,n_basis,one,fock_ao,n_basis,orb_coefficients,n_basis,zero,X,n_basis)
      call dgemm('T','N',n_orbitals,n_orbitals,n_basis,one,orb_coefficients,n_basis,X,n_basis,zero,mo_fock_mat,n_orbitals)

      call deallocator(X,n_basis,n_orbitals)
      call deallocator(fock_ao,n_basis,n_basis)
!
!
!!! TWO-ELECTRON CONTRIBUTION !!!
!
!
!
!!  occupied-occupied block: F_ij = h_ij + sum_k (2*g_ijkk - g_ikkj) !!
!  
!
!  Allocation for L_ij_J_pack
! 
   call allocator(L_ij_J,n_oo,n_J)
   call allocator(g_ij_kl,n_oo,n_oo)
   L_ij_J=zero
   g_ij_kl=zero
!
!  Read cholesky vectors
!
   call get_cholesky_ij(L_ij_J)
!
!  g_ij_kl = sum_J(L_ij_J*L_J_kl)
!
   call dgemm('N','T',n_oo,n_oo,n_J &
      ,one,L_ij_J,n_oo,L_ij_J,n_oo &
      ,zero,g_ij_kl,n_oo)
!
!  Add two-electron contributions to occupied-occupied block
!
   do i=1,n_occ
      do j=1,n_occ
         ij=index_two(i,j,n_occ)
         do k = 1,n_occ
            kk=index_two(k,k,n_occ)
            ik=index_two(i,k,n_occ)
            kj=index_two(k,j,n_occ)
            mo_fock_mat(i,j)=mo_fock_mat(i,j)+two*g_ij_kl(ij,kk)-g_ij_kl(ik,kj)
         enddo
      enddo
   enddo
!
!     Deallocate g_ij_kl
! 
      call deallocator(g_ij_kl,n_oo,n_oo)
!
!!    Occupied-vacant blocks F_ai=F_ia=0 because this Fock matrix satisfies the HF equations !! OBS T1!
!
!
!     Allocation for g_ia_jk 
!
     call allocator(L_ia_J,n_ov,n_J)
     call allocator(g_ia_jk,n_ov,n_oo)
!
!     Reading Cholesky vector L_ia_J
!
      call get_cholesky_ia(L_ia_J)
!
!     g_ia_jk
!
      call dgemm('N','T',n_ov,n_oo,n_J,one,L_ia_J,n_ov,L_ij_J,n_oo,zero,g_ia_jk,n_ov)
      call deallocator(L_ia_J,n_ov,n_J)
!
!     Allocation for g_ai_jk 
!
      call allocator(L_ai_J,n_ov,n_J)
      call allocator(g_ai_jk,n_ov,n_oo)
!
!     Reading Cholesky vector L_ai_J
!
      call get_cholesky_ai(L_ai_J)
!
!     g_ai_jk
!
      call dgemm('N','T',n_ov,n_oo,n_J,one,L_ai_J,n_ov,L_ij_J,n_oo,zero,g_ai_jk,n_ov)
      call deallocator(L_ai_J,n_ov,n_J)
!
!     Adding terms to Fock matrix
!
      do i=1,n_occ
         do a=1,n_vir
            do j=1,n_occ
!
!              Needed indices
!
               ia=index_two(i,a,n_occ)
               ja=index_two(j,a,n_occ)
!
               ai=index_two(a,i,n_vir)
               aj=index_two(a,j,n_vir)
!
               jj=index_two(j,j,n_occ)
               ji=index_two(j,i,n_occ)
               ij=index_two(i,j,n_occ)
!
               mo_fock_mat(i,a+n_occ)=mo_fock_mat(i,a+n_occ)+two*g_ia_jk(ia,jj)-g_ia_jk(ja,ij)
               mo_fock_mat(a+n_occ,i)=mo_fock_mat(a+n_occ,i)+two*g_ai_jk(ai,jj)-g_ai_jk(aj,ji)
            enddo
         enddo
      enddo
      call deallocator(g_ia_jk,n_ov,n_oo)
      call deallocator(g_ai_jk,n_ov,n_oo)
!
!!    Vacant-vacant block F_ab = h_ab + sum_k (2*g_abkk - g_akkb) !!
!
      call allocator(g_ab_ij,n_vv,n_oo)
      g_ab_ij=zero
!
!     Batching over a
!
!
!     Setup of variables needed for batching
!
      available = get_available()
      required = n_vir*n_vir*n_J*4
!
      call n_one_batch(required,available,max_batch_length,n_batch,n_vir)
!
      batch_start=1
      batch_end=0
      batch_length=0
!
!     Start batchig loop
!
      do b_batch = 1,n_batch

!
!        Get batch limits  and  length of batch
!
         call one_batch_limits(batch_start,batch_end,b_batch,max_batch_length,n_vir)
         batch_length=batch_end-batch_start+1
!
!        Allocation of L_ab_J
!
         call allocator(L_ab_J,n_vir*batch_length,n_J)
         L_ab_J=zero
!
!        Read Cholesky vectors
!
         call get_cholesky_ab(L_ab_J,batch_start,batch_end,n_vir*batch_length,.false.)
!
!        g_ab_ij=sum_J L_ab_J* L_ij_J
!
         g_off = index_two(1,batch_start,n_vir)
!
         call dgemm('N','T',n_vir*batch_length,n_oo,n_J &
            ,one,L_ab_J,n_vir*batch_length,L_ij_J,n_oo &
            ,one,g_ab_ij(g_off,1),n_vv)
!
!        Deallocation of L_ab_J
!
         call deallocator(L_ab_J,batch_length*n_vir,n_J)
!
      enddo ! batching done
!
!     Deallocation of L_ij_J
!
      call deallocator(L_ij_J,n_oo,n_J)
!
!     Allocate for g_ai_jb
!
      call allocator(g_ai_jb,n_ov,n_ov)
      call allocator(L_ai_J,n_ov,n_J)
      call allocator(L_ia_J,n_ov,n_J)
!
!
!     Reading Cholesky vector L_ia_J and L_ai_J
!
      call get_cholesky_ia(L_ia_J)
      call get_cholesky_ai(L_ai_J)
      call dgemm('N','T',n_ov,n_ov,n_J,one,L_ai_J,n_ov,L_ia_J,n_ov,zero,g_ai_jb,n_ov)
!
!     Deallocate L_ia_J
!
     call deallocator(L_ia_J,n_ov,n_J)
     call deallocator(L_ai_J,n_ov,n_J)
!
!     Calculation of two-electron terms for virtual-virtual blocks
!
      do a = 1,n_vir
         do b = 1,n_vir
            ab=index_two(a,b,n_vir)
            do i = 1,n_occ
               ii=index_two(i,i,n_occ)
               ai=index_two(a,i,n_vir)
               bi=index_two(b,i,n_vir)
               ia=index_two(i,a,n_occ)
               ib=index_two(i,b,n_occ)
               mo_fock_mat(n_occ+a,n_occ+b)=mo_fock_mat(n_occ+a,n_occ+b)+two*g_ab_ij(ab,ii)-g_ai_jb(ai,ib)
            enddo
         enddo 
      enddo
     call deallocator(g_ab_ij,n_vv,n_oo)
     call deallocator(g_ai_jb,n_ov,n_ov)
!
!    Clean-up of Fock matrix
!
     call mlcc_cleanup(mo_fock_mat,n_orbitals,n_orbitals)
!
!     Prints
!
      if (debug) then
         do i=1,n_occ
            write(luprint,*)(mo_fock_mat(i,j),j=1,n_occ)
         enddo
         do i=1,n_vir
            write(luprint,*)(mo_fock_mat(j,i+n_occ),j=1,n_occ)
         enddo
         do i=1,n_vir
            write(luprint,*)(mo_fock_mat(i+n_occ,j),j=1,n_occ)
         enddo
         do i=1,5
            write(luprint,*)(mo_fock_mat(i+n_occ,j+n_occ),j=1,5) 
         enddo
      endif
!
!     Save the blocks of the Fock matrix in memory (ij,ia,ai,ab)
!
      F_i_j = zero
      F_i_a = zero
      F_a_i = zero
      F_a_b = zero
!
      do i = 1,n_occ
         do j = 1,n_occ
            F_i_j(i,j) = mo_fock_mat(i,j)
         enddo
      enddo
!
      do i = 1,n_occ
         do a = 1,n_vir
            F_i_a(i,a) = mo_fock_mat(i,n_occ+a)
            F_a_i(a,i) = mo_fock_mat(n_occ+a,i)
         enddo
      enddo
!
      do a = 1,n_vir
         do b = 1,n_vir
            F_a_b(a,b) = mo_fock_mat(n_occ+a,n_occ+b)
         enddo
      enddo
!

   call deallocator(mo_fock_mat,n_orbitals,n_orbitals)
!
   end subroutine mlcc_fock

