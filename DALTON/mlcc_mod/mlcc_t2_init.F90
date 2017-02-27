module mlcc_t2_init
   contains
   subroutine t2_init(t2)
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
      real(dp),dimension(:,:)          :: t2
      real(dp),dimension(:,:),pointer  :: cholesky_ia,g_iajb
      integer                          :: lucho_ia,idum,n_aibj
      integer                          :: i,j,a,b,ai,bj,aibj
   !
   !
      n_aibj = n_vir*n_occ
   !
      call allocator(cholesky_ia,n_aibj,n_J)
      call allocator(g_iajb,n_aibj,n_aibj)
      cholesky_ia=zero
      g_iajb=zero
   !
   !  IO: read from CHOLESKY_IA
   !
      lucho_ia=-1
      call gpopen(lucho_ia,'CHOLESKY_IA','UNKNOWN','SEQUENTIAL','UNFORMATTED',idum,.false.)
      rewind(lucho_ia)
   !
      do j=1,n_J
         read(lucho_ia)(cholesky_ia(i,j),i=1,n_aibj)
      enddo
   !
      call gpclose(lucho_ia,'KEEP')
   !
   !
   !  g_iajb=L^J_ia L^J_jb
   !
      call dgemm('N','T',n_aibj,n_aibj,n_J,one,cholesky_ia,n_aibj,cholesky_ia,n_aibj,zero,g_iajb,n_aibj) 
   !
   !
      call packin(t2,g_iajb,n_aibj)
   !
   do a=1,n_vir
      do b=1,n_vir
         do i=1,n_occ
            do j=1,n_occ
               ai=index_t1(i,a)
               bj=index_t1(j,b)
               aibj=index_t2(ai,bj)
               t2(aibj,1)=t2(aibj,1)/(fock_diagonal(n_occ+a,1)+fock_diagonal(n_occ+b,1) &
                  -fock_diagonal(i,1)-fock_diagonal(j,1))
            enddo
         enddo
      enddo
   enddo
   call mp2_energy(t2,g_iajb)
   !
      call deallocator(cholesky_ia,n_aibj,n_J)
      call deallocator(g_iajb,n_aibj,n_aibj)
   !   
   end subroutine t2_init
   !
   subroutine mp2_energy(t2,g_iajb)
   !  Only works if orbitals are canonical!!!
   !  E^(2)=sum_iajb t^(1)_aibj*L_iajb
      use mlcc_data
      use mlcc_utilities
      use mlcc_workspace
   !
      implicit none
      real(dp),dimension(:,:)                :: t2,g_iajb
      real(dp),dimension(:,:), pointer       :: L_iajb
      integer                                :: i,j,a,b,nai,naj,nbi,nbj,naibj,najbi,n_occvirr
      real(dp)                               :: E
   !
      n_occvirr=n_occ*n_vir
      call allocator(L_iajb,n_occvirr,n_occvirr)
      L_iajb=zero
   !
      E=zero
      do i=1,n_occ
         do j=1,n_occ
            do a=1,n_vir
               do b=1,n_vir
                  nai=index_t1(i,a)
                  naj=index_t1(j,a)
                  nbi=index_t1(i,b)
                  nbj=index_t1(j,b)
                  L_iajb(nai,nbj)=2*g_iajb(nai,nbj)-g_iajb(naj,nbi)
                  E=E+L_iajb(nai,nbj)*t2(index_t2(nai,nbj),1)
               enddo
            enddo
         enddo
      enddo
   !
   write(*,*)'MP2 energy correction', E
   write(*,*)'MP2 energy',E+scf_energy
   
      call deallocator(L_iajb,n_occ*n_vir,n_occ*n_vir)
   !
   end subroutine mp2_energy
end module mlcc_t2_init
   