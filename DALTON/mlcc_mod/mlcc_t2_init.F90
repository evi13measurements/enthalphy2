module mlcc_t2_init
contains
   subroutine t2_init
   !
   ! Purpose: Calculate initial guess for t2 amplitudes
   !
   !  t2(aibj)=g_iajb/(e_a+e_b-e_j-e_i)
   !
   !  Author: Sarai Folkestad
   !
      use mlcc_data
      use mlcc_utilities
      use mlcc_workspace
   !
      implicit none
   !
      real(dp),dimension(:,:),pointer  :: cholesky_ia,g_iajb   => null()
      integer                          :: lucho_ia,idum
      integer                          :: i,j,a,b,ai,bj,aibj
   !
   !
   !
   !
   ! Allocations
   !
      call allocator(cholesky_ia,n_ov,n_J)
      call allocator(g_iajb,n_ov,n_ov)
      cholesky_ia=zero
      g_iajb=zero
   !
   !  IO: read cholesky vectors from CHOLESKY_IA
   !
      lucho_ia=-1
      call gpopen(lucho_ia,'CHOLESKY_IA','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
      rewind(lucho_ia)
   !
   !
   !
      do j=1,n_J
         read(lucho_ia)(cholesky_ia(i,j),i=1,n_ov)
      enddo
   !
      call gpclose(lucho_ia,'KEEP')
   !
   !  Two electron integrals g_iajb = L^J_ia L^J_jb
   !
      call dgemm('N','T',n_ov,n_ov,n_J,one,cholesky_ia,n_ov,cholesky_ia,n_ov,zero,g_iajb,n_ov) 
   !
   ! t2 amplitude guess
   !
      call packin(t2am,g_iajb,n_ov)
   !
   do a=1,n_vir
      do b=1,n_vir
         do i=1,n_occ
            do j=1,n_occ
               ai = index_two(a,i,n_vir)
               bj = index_two(b,j,n_vir)
               if (ai .le. bj) then
                  aibj = index_packed(ai,bj)
                  t2am(aibj,1) = -t2am(aibj,1)/(fock_diagonal(n_occ+a,1)+fock_diagonal(n_occ+b,1) &
                                                            -fock_diagonal(i,1)-fock_diagonal(j,1))
               endif
            enddo
         enddo
      enddo
   enddo
   !
   ! MP2 energy. OBS! Only works if orbitals are canonical
   !
   call mp2_energy(t2am,g_iajb)
   !
   ! Deallocations
   !
   call deallocator(cholesky_ia,n_ov,n_J)
   call deallocator(g_iajb,n_ov,n_ov)
   !   
   end subroutine t2_init
   !
   subroutine mp2_energy(t2,g_iajb)
   !
   !  Purpose: Calculate MP2 energy.  Only works if orbitals are canonical!
   !  E^(2)=sum_iajb t^(1)_aibj*L_iajb
   !
   !  Author: Sarai Folkestad
   !
      use mlcc_data
      use mlcc_utilities
      use mlcc_workspace
   !
      implicit none
      real(dp),dimension(:,:)                :: t2,g_iajb
      real(dp),dimension(:,:), pointer       :: L_iajb   => null()
      integer                                :: i,j,a,b,ai,aj,bi,bj,aibj,ajbi
      real(dp)                               :: E
   !
   !
   ! Allocation
   !
      call allocator(L_iajb,n_ov,n_ov)
      L_iajb=zero
   !
   ! Calculate correction to energy
   !
      E=zero
      do i=1,n_occ
         do j=1,n_occ
            do a=1,n_vir
               do b=1,n_vir
                  ai=index_two(a,i,n_vir)
                  aj=index_two(a,j,n_vir)
                  bi=index_two(b,i,n_vir)
                  bj=index_two(b,j,n_vir)
                  L_iajb(ai,bj)=two*g_iajb(ai,bj)-g_iajb(aj,bi)
                  E=E+L_iajb(ai,bj)*t2(index_packed(ai,bj),1)
               enddo
            enddo
         enddo
      enddo
   !
   ! Printouts
   !
   write(ml_lupri,*)'MLCC MP2 ENERGY - Only correct for canonical orbitals'
   write(ml_lupri,*)'MP2 energy correction', E
   write(ml_lupri,*)'MP2 energy',E+scf_energy
   !
   ! Deallocation
   !
      call deallocator(L_iajb,n_ov,n_ov)
   !
   end subroutine mp2_energy
end module mlcc_t2_init
   