      subroutine mlcc_rewr_sirifc(work,lwork)
C
C     Writing info for mlcc calculation to MLCC_HF_INFO
C  
C     Date:    April 2017
C     Authors: Henrik Koch, Rolf H. Myhre, 
C              Sarai D. Folkestad, Eirik F. Kjønstad
C
      implicit none 
C
#include "priunit.h"
C
      integer lu_sir 
      integer lu_mlcc 
      integer idummy 
      integer n_symmetries, n_basis_sym, n_orbitals_sym 
      integer n_lambda
      integer n_occ, n_orbitals, n_basis        
      integer i,j
      integer lwork
      integer fock_diagonal, orbital_coeffs 
      integer lend 
      integer lfree

      double precision nuclear_potential, scf_energy  
      double precision, dimension(lwork) :: work
C     
C     Open Sirius Fock file
C     ---------------------
C
      call gpopen(lu_sir,'SIRIFC','OLD',' ','UNFORMATTED',
     &             idummy,'.FALSE.')
      rewind(lu_sir)
C
C     Read in various stuff from Sirius Fock file. Things depending on symmetry is mostly
C     discarded at the end of the subroutine as we do not use symmetry. Information in 
C     file should be in Cholesky orbital format. If Cholesky orbitals has been generated,
C     SIRIFC will contain the data in the Cholesky basis
C
      call mollab('TRCCINT ',lu_sir,lupri)
C      
      read(lu_sir) n_symmetries, n_orbitals, n_basis, n_lambda, n_occ,
     &            n_orbitals_sym, n_basis_sym, 
     &            nuclear_potential, scf_energy
C
C      
      if (n_symmetries .ne. 1)
     &  call quit('error in mlcc: not implemented with symmetry')
C     
C     Memmory handling
C
      fock_diagonal  = 1
      orbital_coeffs = fock_diagonal  + n_orbitals
      lend           = orbital_coeffs + n_lambda
      lfree          = lwork - lend + 1 
C
      if (lfree .lt. 0) then
         write(lupri,*)'Not enough memory in mlcc_write_sirifc'
      endif
C
C     Read fock diagonal and orbital coefficients
C
      read(lu_sir) (work(fock_diagonal+i-1),i=1,n_orbitals)
      read(lu_sir) (work(orbital_coeffs+i-1),i=1,n_lambda)
C
      call gpclose(lu_sir,'KEEP')
C
C     Write to mlcc_hf_info
C
      call  gpopen(lu_mlcc,'MLCC_HF_INFO','UNKNOWN','SEQUENTIAL',
     &        'FORMATTED',idummy,.false.)
      rewind(lu_mlcc)
C
      write(lu_mlcc,*) n_orbitals, n_occ, n_lambda,
     &                 nuclear_potential, scf_energy
C
      write(lu_mlcc,*) (work(fock_diagonal+i-1),i=1,n_orbitals)
C
      write(lu_mlcc,*) (work(orbital_coeffs+i-1),i=1,n_lambda)
C
      call gpclose(lu_mlcc,'KEEP')
C
      end subroutine mlcc_rewr_sirifc
