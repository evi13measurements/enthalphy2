      subroutine mlcc_rewr_sirifc(work,lwork)
#include "implicit.h"
      integer  :: lu_sir   = -1 
      integer  :: lu_mlcc  = -1 
      integer  :: idummy   = 0
      integer  :: n_symmetries, n_basis_sym, n_orbitals_sym, n_lambda
      integer  :: n_occ,n_orbitals_sym, n_basis_sym 
      integer. :: nuclear_potential, scf_energy          
      integer  :: i,j
C     
C     Open Sirius Fock file
C     ---------------------
C
      lu_sir = -1
      lu_mlcc  = -1
      call gpopen(lu_sir,'SIRIFC','OLD',' ','UNFORMATTED',idummy,'.FALSE.')
      rewind(lu_sir)
C
C     Read in various stuff from Sirius Fock file. Things depending on symmetry is mostly
C     discarded at the end of the subroutine as we do not use symmetry. Information in 
C     file should be in Cholesky orbital format. If Cholesky orbitals has been generated,
C     SIRIFC will contain the data in the Cholesky basis
C
      call mollab('TRCCINT ',lu_sir,lupri)
C      
      read(lu_sir) n_symmetries, n_orbitals, n_basis, n_lambda, n_occ, &
      &            n_orbitals_sym, n_basis_sym, &
      &            nuclear_potential, scf_energy
C
C      
      if (n_symmetries /= 1) call quit('error in mlcc_mod_init: not implemented with symmetry')

      call gpclose(lu_sir,'KEEP')
!
      end subroutine mlcc_rewr_sirifc
