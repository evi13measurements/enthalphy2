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
      use mlcc_cholesky
   !
      implicit none
   !
      real(dp),dimension(:,:),pointer  :: cholesky_ia,g_iajb   => null()
      integer                          :: lucho_ia,idum
      integer                          :: i,j,a,b,ai,bj,aibj,ia,jb
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
      call read_cholesky_ia(cholesky_ia)
   !  Two electron integrals g_iajb = L^J_ia L^J_jb
   !
      call dgemm('N','T',n_ov,n_ov,n_J &
         ,one,cholesky_ia,n_ov,cholesky_ia,n_ov & 
         ,zero,g_iajb,n_ov) 
   !
   ! t2 amplitude guess
   !
   !
   do a=1,n_vir
      do b=1,n_vir
         do i=1,n_occ
            do j=1,n_occ
               ai = index_two(a,i,n_vir)
               bj = index_two(b,j,n_vir)
               ia = index_two(i,a,n_occ)
               jb = index_two(j,b,n_occ)
               if (ai .le. bj) then
                  aibj = index_packed(ai,bj)
                  t2am(aibj,1) = -g_iajb(ia,jb)/(fock_diagonal(n_occ+a,1)+fock_diagonal(n_occ+b,1) &
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
      integer                                :: i,j,a,b,ia,ja,ib,jb,aibj,ajbi,bj,ai
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
!                 Calculating all needed indices
                  ia=index_two(i,a,n_occ)
                  ja=index_two(j,a,n_occ)
                  ib=index_two(i,b,n_occ)
                  jb=index_two(j,b,n_occ)
                  ai=index_two(a,i,n_vir)
                  bj=index_two(b,j,n_vir)
                  aibj=index_packed(ai,bj)
!                 Constructing L_iajb
                  L_iajb(ia,jb)=two*g_iajb(ia,jb)-g_iajb(ja,ib)
!                 Calculating energy contribution
                  E=E+L_iajb(ia,jb)*t2(aibj,1)
               enddo
            enddo
         enddo
      enddo
   !
   ! Printouts
   !
   write(luprint,*)'MLCC MP2 ENERGY - Only correct for canonical orbitals'
   write(luprint,*)'MP2 energy correction', E
   write(luprint,*)'MP2 energy',E+scf_energy
   !
   ! Deallocation
   !
      call deallocator(L_iajb,n_ov,n_ov)
   !
   end subroutine mp2_energy
end module mlcc_t2_init
   