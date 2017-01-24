module mlcc_mod_init
!
!  MLCC initialiasion routines
!  Authors Henrik Koch, Rolf H. Myhre, Sarai Folkestad, Eirik Kjønstad
!  January 2017
!
   use mlcc_mod_work
   use mlcc_data
!
contains
   subroutine hf_reader(lupri)
!
!  Hartree-Fock reader routine
!  Authors Henrik Koch, Rolf H. Myhre, Sarai Folkestad, Eirik Kjønstad
!  January 2017
!
!  Purpose: read in data from LUSIFC file
!
      implicit none
!
      integer, intent(in)                    :: lupri !general output unit
!
      integer  :: lusifc = -1 
      integer  :: idummy = 0
      integer  :: n_symmetries, n_basis_sym, n_orbitals_sym
      integer  :: i
!      
!
!     
!     Open Sirius Fock file
!     ---------------------
!
      call gpopen(lusifc,'SIRIFC','OLD',' ','UNFORMATTED',idummy,'.FALSE.')
      rewind(lusifc)
!
!
!     Read in various stuff from Sirius Fock file. Things depending on symmetry is mostly
!     discarded at the end of the subroutine as we do not use symmetry. Information in 
!     file should be in Cholesky orbital format. If Cholesky orbitals has been generated,
!     SIRIFC will contain the data in the Cholesky basis
!
      call mollab('TRCCINT ',lusifc,lupri)
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
      n_t1am         = n_vir*n_occ
      n_t2am         = n_t1am*n_t1am
      n_t2am_pack    = n_t1am*(n_t1am+1)/2
!      
!     Allocate space for Fock diagonal and coefficients.
!
      call allocator(fock_diagonal,n_orbitals,1)
!      
      call allocator(orb_coefficients,n_lambda,1)
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
end module mlcc_mod_init