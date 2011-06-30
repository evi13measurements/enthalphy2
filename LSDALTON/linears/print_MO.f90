module print_moorb_grid_mod
  use matrix_operations
  use BUILDAOBATCH
  use typedef
contains
  subroutine print_moorb(ls,CMO,filename,nocc,onlyocc)
    implicit none
    type(matrix)                   :: CMO
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer                        :: nocc
    logical                        :: onlyocc
!
    integer                        :: I, J, nATOMS,nbast,nX,nY,nZ
    real(realk)                    :: deltax, deltay, deltaz,cm(3)
    logical,allocatable            :: atom_inside_grid(:)
    real(realk), allocatable       :: ATOMXYZ(:,:),newXYZ(:,:)

    nATOMS = ls%setting%MOLECULE(1)%p%natoms
    allocate(atom_inside_grid(nATOMS))
    allocate(ATOMXYZ(3,nATOMS))
    allocate(newXYZ(3,nATOMS))
    call calculate_cm(cm,ls) !center of charge
    ! Loop over all atoms (natoms) and store atomic coordinates in ATOMXYZ
    DO I=1,nATOMS !Atom coordinates        
       DO J=1,3                                                         
          ATOMXYZ(J,I) = ls%setting%molecule(1)%p%atom(I)%CENTER(J)
          newXYZ(J,I) = ATOMXYZ(J,I)-cm(J)
       ENDDO
    ENDDO
    ! Default: Gridbox encapsulates all atoms and the distances between gridpoints are given by:
    ! deltax = deltay = deltaz = 0.3 (a.u.)
!    deltax = 0.3d0; deltay=deltax; deltaz=deltax
    nX = 100
    nY = 100
    nZ = 100
    atom_inside_grid = .true.
    call calculate_MOORB(ls,nATOMS,atom_inside_grid, ATOMXYZ,&
         & natoms, nX,nY,nZ,CMO,filename,nocc,onlyocc)
    deallocate(atom_inside_grid, ATOMXYZ)
  end subroutine print_moorb

  subroutine calculate_cm(cm,ls)
    implicit none
    type(lsitem)                   :: ls
    real(realk)                    :: cm(3)
    !
    real(realk) :: fullcharge,X,Y,Z,charge
    integer :: iatom

    CM=0.d0
    fullcharge = 0.d0
    do iatom = 1,ls%setting%MOLECULE(1)%p%natoms
       charge = ls%setting%MOLECULE(1)%p%ATOM(Iatom)%charge
       fullcharge = fullcharge + charge
       X = ls%setting%MOLECULE(1)%p%ATOM(Iatom)%CENTER(1)
       Y = ls%setting%MOLECULE(1)%p%ATOM(Iatom)%CENTER(2)
       Z = ls%setting%MOLECULE(1)%p%ATOM(Iatom)%CENTER(3)
       CM(1) = CM(1)+charge*X 
       CM(2) = CM(2)+charge*Y 
       CM(3) = CM(3)+charge*Z 
    enddo
    CM(1)=CM(1)/fullcharge
    CM(2)=CM(2)/fullcharge
    CM(3)=CM(3)/fullcharge
  end subroutine calculate_cm

  subroutine calculate_moorb(ls,nATOMS_grid,atom_inside_grid, &
       & ATOMXYZ, natoms, nX,nY,nZ, CMO,filename,nocc,onlyocc)
    !*********************************************************
    ! Driver routine for calculating the molecular orbitals
    ! in specified points r in space 
    ! \phi(r)_{m} = \sum_{n} C_{mn} X(r)_{n}
    ! where X(r)_[n] are atomic orbitals.
    !*********************************************************
    implicit none
    type(matrix)                   :: CMO
    type(lsitem)                   :: ls
    Character*(*)                  :: filename
    integer,intent(in)             :: natoms, nATOMS_grid,nX,nY,nZ,nocc 
    logical, intent(in)            :: atom_inside_grid(nATOMS),onlyocc
    real(realk), intent(in)        :: ATOMXYZ(3,nATOMS)
    !Local
    integer                        :: I,J,P,Q,Xg,Yg,Zg,nGRIDPOINTS, gridnr,iFilename,iFilename2
    integer                        :: orbnr, nORBITALS,iorb,iunit,ig
    real(realk)                    :: X,Y,Z,X1,Y1,Z1,Xgrid, Ygrid, Xc, Yc, Zc,deltax, deltay, deltaz 
    real(realk), allocatable       :: GRIDCOOR(:,:), GAO(:),moorb(:,:)
    Character(len=80)              :: filename2
    Character(len=1)  :: STRING1
    Character(len=2)  :: STRING2
    Character(len=3)  :: STRING3
    Character(len=4)  :: STRING4
    Character(len=5)  :: STRING5
    Character(len=6)  :: STRING6

    nORBITALS = CMO%nrow
    allocate(GAO(nORBITALS))
    ! Determine grid-box. It extends from
    ! X1/Y1/Z1 to Xn/Yn/Zn in the x-/y-/z-directions,
    ! and the number of gridpoints in these directions are nX/nY/nZ.
    ! Thus, the total number of gridpoints is nGRIDPOINT=nX*nY*nZ.
    call DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
         &ATOMXYZ,natoms,nGRIDPOINTS, atom_inside_grid)
    ! Center of grid.
    Xc = X1 + 0.5d0*deltax*(nX-1)
    Yc = Y1 + 0.5d0*deltay*(nY-1)
    Zc = Z1 + 0.5d0*deltaz*(nZ-1)
    allocate(GRIDCOOR(3,nGRIDPOINTS))  !xyz coordinates in grid
    allocate(moorb(nGRIDPOINTS,CMO%nrow)) !MO orbitals
    moorb = 0
    ! gridnr: Counter for grid points, orbnr: Counter for orbitals
    gridnr=0; 
    !Loop over xyz-grid points to determine 3D box of grid points
    DO Xg = 1,nX ! Xg: grid-coordinate number in x direction
       Xgrid = X1 + deltax*(Xg-1) !x grid-coordinate
       DO Yg = 1,nY ! Yg: grid-coordinate number in y direction
          Ygrid = Y1 + deltay*(Yg-1) !y grid-coordinate
          DO Zg = 1,nZ ! Zg: grid-coordinate number in z direction
             gridnr = gridnr+1
             GRIDCOOR(1,gridnr) = Xgrid
             GRIDCOOR(2,gridnr) = Ygrid
             GRIDCOOR(3,gridnr) = Z1 + deltaz*(Zg-1) !z grid-coordinate
             GAO = 0
             orbnr=1
             DO I=1,natoms
                ! X-, Y-, and Z-distances from gridpoint to atom I
                X = GRIDCOOR(1,gridnr) - ATOMXYZ(1,I)
                Y = GRIDCOOR(2,gridnr) - ATOMXYZ(2,I)
                Z = GRIDCOOR(3,gridnr) - ATOMXYZ(3,I)
                call determine_orbitals(I,GAO,orbnr, nORBITALS, X,Y,Z,ls)
             ENDDO
             J = 0
             do iorb = 1,CMO%ncol
                DO P= 1, nORBITALS
                   J = J+1   
                   moorb(gridnr,Iorb) = moorb(gridnr,iorb) + CMO%elms(p+(iorb-1)*CMO%ncol)*GAO(p)
                ENDDO
             enddo
          ENDDO
       ENDDO
    ENDDO
    iFilename2 = 0
    do i=1,80
       if(filename(i:i).EQ.'.')iFilename2 = i
    enddo
    if(iFilename2.NE.0)then
       filename2(1:iFilename2) = FILENAME(1:iFilename2)  
       iFilename = iFilename2+1
       Filename2(iFilename:iFilename+4) = 'cube'
       iFilename = iFilename + 5          
       do i=iFilename,80
          Filename2(i:i) = ' '
       enddo
    else
       iFilename = LEN(filename)
       filename2(1:iFilename) = FILENAME(1:iFilename)  
       Filename2(iFilename:iFilename+4) = 'cube'
       iFilename = iFilename + 5
       do i=iFilename,80
          Filename2(i:i) = ' '
       enddo
    endif
    IUNIT = -1
    print*,'filname2:',trim(filename2)
    CALL LSOPEN(IUNIT,trim(filename2),'UNKNOWN','FORMATTED')
    call write_cubefile(IUNIT,gridnr,moorb,nGRIDPOINTS,CMO%ncol,nocc,onlyocc,natoms,&
         &deltax,deltay,deltaz,nX,nY,nZ,X1,Y1,Z1,ls)
    CALL LSCLOSE(IUNIT,'KEEP')
    deallocate(GRIDCOOR,moorb,GAO)

  end subroutine calculate_moorb

  subroutine write_cubefile(IUNIT,gridnr,moorb,nGRIDPOINTS,nbast,nocc,onlyocc,natoms,&
       &deltax,deltay,deltaz,nX,nY,nZ,X1,Y1,Z1,ls)
    implicit none
    type(lsitem)                   :: ls
    integer,intent(in)             :: IUNIT,natoms, nX,nY,nZ,nocc,nbast,nGRIDPOINTS 
    logical, intent(in)            :: onlyocc
    real(realk), intent(in)        :: X1,Y1,Z1,deltax,deltay,deltaz
    real(realk),intent(in)         :: moorb(nGRIDPOINTS,nbast)
    !
    real(realk)                    :: charge,X,Y,Z
    integer                        :: I,nmo,gridnr,iorb,igr,Xg,Yg,Zg
    Character(len=80)              :: filename2

    WRITE(IUNIT,*)'Dalton Cube File  '
    IF(onlyocc)then
       WRITE(IUNIT,*)'Occupied Molecular Orbitals'
    ELSE
       WRITE(IUNIT,*)'Molecular Orbitals'
    ENDIF
    WRITE(IUNIT,'(I5,F12.6,F12.6,F12.6)')-natoms,X1,Y1,Z1
    WRITE(IUNIT,'(I5,F12.6,F12.6,F12.6)')nX,deltax,0.d0,0.d0
    WRITE(IUNIT,'(I5,F12.6,F12.6,F12.6)')nY,0.d0,deltay,0.d0
    WRITE(IUNIT,'(I5,F12.6,F12.6,F12.6)')nZ,0.d0,0.d0,deltaz
    do i=1,natoms
       charge = ls%setting%molecule(1)%p%ATOM(i)%charge
       X = ls%setting%molecule(1)%p%ATOM(i)%center(1)
       Y = ls%setting%molecule(1)%p%ATOM(i)%center(2)
       Z = ls%setting%molecule(1)%p%ATOM(i)%center(3)
       WRITE(IUNIT,'(I5,F12.6,F12.6,F12.6,F12.6)')INT(charge),charge,X,Y,Z
    enddo
    if(onlyocc)then
       nmo = nocc
    else
       nmo = nbast
    endif
    WRITE(IUNIT,'(10I5)')nmo,(I,I =1,nmo)   
    gridnr=1 
    DO Xg = 1,nX 
       DO Yg = 1,nY 
          DO Zg = 1,nZ 
           gridnr = gridnr+1
           write(iunit,'(6E13.5)') (moorb(gridnr,Iorb),iorb=1,nmo)
       ENDDO
    ENDDO
    enddo
  end subroutine write_cubefile

  subroutine DETERMINE_GRIDBOX(X1,nX,Y1,nY,Z1,nZ,deltax,deltay,deltaz,&
       &ATOMXYZ,natoms,nGRIDPOINTS, atom_inside_grid)
    implicit none
    integer, intent(in)          :: nX,nY,nZ,natoms
    real(realk), intent(in)      :: ATOMXYZ(3,natoms)
    logical, intent(in)          :: atom_inside_grid(natoms)
    integer, intent(out)         :: nGRIDPOINTS
    real(realk), intent(out)     :: deltax,deltay,deltaz,X1,Y1,Z1
    real(realk)                  :: Xn, Yn, Zn,distX,distY,distZ,buffer
    integer                      :: I
    ! Minimum and maximum values in the gridbox. 
    X1 = HUGE(1.0)
    Y1 = X1; Z1 = X1; Xn = - X1; Yn = - X1; Zn = - X1 
    do I = 1,natoms
       if(atom_inside_grid(I)) then
          if(ATOMXYZ(1,I)>Xn) Xn = ATOMXYZ(1,I)
          if(ATOMXYZ(2,I)>Yn) Yn = ATOMXYZ(2,I)
          if(ATOMXYZ(3,I)>Zn) Zn = ATOMXYZ(3,I)
          if(ATOMXYZ(1,I)<X1) X1 = ATOMXYZ(1,I)
          if(ATOMXYZ(2,I)<Y1) Y1 = ATOMXYZ(2,I)
          if(ATOMXYZ(3,I)<Z1) Z1 = ATOMXYZ(3,I)
       endif
    enddo

    ! Add/subtract arbitrary "3.0 a.u." to ensure that we get the electron density
    ! in a box which is larger than then box containing the atoms of interest.
    buffer = 3.d0
    X1 = X1 - buffer; Y1 = Y1 - buffer; Z1 = Z1 - buffer
    Xn = Xn + buffer; Yn = Yn + buffer; Zn = Zn + buffer
    distX = Xn-X1
    distY = Yn-Y1
    distZ = Zn-Z1
    deltax = distX/nX
    deltay = distY/nY
    deltaz = distZ/nZ
    ! Number of gridpoints from: X1 to Xn; Y1 to Yn; Z1 to Zn.
    nGRIDPOINTS = nX*nY*nZ

  end subroutine DETERMINE_GRIDBOX

  subroutine determine_orbitals(I,GAO,orbnr, nORBITALS, X, Y, Z, ls)
    ! Determine value of all atomix orbitals on atom I at point (X,Y,Z)
    ! measured relative to atom center.
    ! The values are stored in the GAO vector.
    implicit none
    type(lsitem)                   :: ls
    integer, intent(in)         :: I, nORBITALS
    integer, intent(inout)       :: orbnr
    real(realk), intent(in)     :: X,Y,Z
    real(realk), intent(inout)  :: GAO(nORBITALS)
    ! Local
    integer                     :: basindex,type,set,nAngmom,J,K,L,M,nsegments,nEXPONENTS,nORBJ,P,Q,R,icharge,itype
    real(realk)                 :: coeff,exponent,R2,contracted_gaussians,GAX,GAY,GAZ,GAXX,GAYY,CINT
    integer, allocatable        :: AVALUE(:), BVALUE(:),CVALUE(:)
    real, allocatable           :: fxyz(:)

    ! total distance from gridpoint to atom I, squared
    R2 = X**2+Y**2+Z**2
    ! basindex = 1 or 2 (2 for auxillary basis sets)
    ! set: integer for describing which basis set is used for atom I.
    ! type: integer for describing the atomtype of atom I.
    IF(LS%setting%basis(1)%p%regular%labelindex .EQ.0)THEN
       icharge = INT(ls%setting%MOLECULE(1)%p%ATOM(i)%charge) 
       itype = ls%setting%basis(1)%p%regular%chargeindex(icharge)
    ELSE
       itype = ls%setting%MOLECULE(1)%p%ATOM(i)%IDtype(1)
    ENDIF
    nAngmom = LS%setting%BASIS(1)%p%REGULAR%ATOMTYPE(itype)%nAngmom
    DO J =1,nAngmom ! Loop over angular moments
       nsegments=LS%setting%BASIS(1)%p%REGULAR%ATOMTYPE(itype)%SHELL(J)%nsegments
       DO K = 1,nsegments ! Loop over segments in basis set
          nEXPONENTS = LS%setting%BASIS(1)%p%REGULAR%ATOMTYPE(itype)%SHELL(J)%segment(K)%nrow
          nORBJ = LS%setting%BASIS(1)%p%REGULAR%ATOMTYPE(itype)%SHELL(J)%segment(K)%ncol
          DO L = 1,nORBJ 
             DO M = 1,nEXPONENTS 
                coeff = LS%setting%basis(1)%p%REGULAR%ATOMTYPE(itype)%SHELL(J)&
                     &%segment(K)%elms(M+(L-1)*nEXPONENTS)
                exponent = LS%setting%basis(1)%p%REGULAR%ATOMTYPE(itype)%SHELL(J)&
                     &%segment(K)%exponents(M)
                GAO(orbnr)=GAO(orbnr)+coeff*EXP(-exponent*R2)
             ENDDO
             IF (J.EQ.1) THEN !1 S-orbital
                ! Orbital is already constructed and stored in GAO(orbnr).
                orbnr=orbnr+1
             ELSEIF (J.EQ.2) THEN ! 3 P-orbitals
                ! Multiply orbital by x, y, or z for px, py, and pz orbitals.
                contracted_gaussians = GAO(orbnr)
                GAO(orbnr)=X*contracted_gaussians
                orbnr=orbnr+1
                GAO(orbnr)=Y*contracted_gaussians
                orbnr=orbnr+1
                GAO(orbnr)=Z*contracted_gaussians
                orbnr=orbnr+1
             ELSEIF (J.EQ.3) THEN ! 5 D spherical orbitals
                GAX  = X*GAO(orbnr)
                GAY  = Y*GAO(orbnr)
                GAZ  = Z*GAO(orbnr)
                GAXX = X*GAX
                GAYY = Y*GAY

                GAO(orbnr) = Y*GAX
                orbnr=orbnr+1

                GAO(orbnr) = Y*GAZ
                orbnr=orbnr+1

                GAO(orbnr) = -0.288675134594813*(GAXX +GAYY)&
                     &+ 0.577350269189626*Z*GAZ
                orbnr=orbnr+1

                GAO(orbnr) = X*GAZ
                orbnr=orbnr+1

                GAO(orbnr) = 0.5d0*GAXX - 0.5d0*GAYY
                orbnr=orbnr+1
             ELSEIF (J.EQ.4) THEN ! 7 F-orbitals
                ! Building the cartesian f-orbitals, fxyz
                ALLOCATE(AVALUE(7*7), BVALUE(7*7), CVALUE(7*7))
                ALLOCATE(fxyz(10))
                ! Constructing the xyz-power of the 10 cartesian f-orbitals
                ! X**AVALUE * Y**BVALUE * Z*CVALUE
                P=0
                DO R = 1,4
                   DO Q = 1,R
                      P=P+1
                      AVALUE(P)=4-R
                      BVALUE(P)=R-Q
                      CVALUE(P)=Q-1
                   ENDDO
                ENDDO
                DO Q = 1, 10
                   fxyz(Q) =&
                        & (X**AVALUE(Q))*(Y**BVALUE(Q))*(Z**CVALUE(Q))*GAO(orbnr)
                ENDDO
                ! The 7 spheric f-orbitals are found as
                ! linear combinations of the "xyz"-f orbitals:
                ! (F1, f2, ..., f7) = (fxyz1, fxyz2, ..., fxyz10) * TM (dimension 10x7)
                ! The transformation matrix TM contains many zeros.
                ! Therefore, the transformations for the 7 spheric 
                ! f-orbitals are written out explicitly.
                GAO(orbnr) = 0.612372435695794*fxyz(2) - 0.204124145231932*fxyz(7)
                orbnr=orbnr+1
                GAO(orbnr) = fxyz(5)
                orbnr=orbnr+1
                GAO(orbnr) =  -0.158113883008419*(fxyz(2)+fxyz(7)) + 0.632455532033676*fxyz(9)
                orbnr=orbnr+1
                GAO(orbnr) =  -0.387298334620742*(fxyz(3)+fxyz(8)) + 0.258198889747161*fxyz(10)
                orbnr=orbnr+1
                GAO(orbnr) = -0.158113883008419*(fxyz(1)+fxyz(4)) + 0.632455532033676*fxyz(6)
                orbnr=orbnr+1
                GAO(orbnr) = 0.5d0*fxyz(3) - 0.5d0*fxyz(8)
                orbnr=orbnr+1
                GAO(orbnr) = 0.204124145231932*fxyz(1) - 0.612372435695794*fxyz(4)
                orbnr=orbnr+1
                DEALLOCATE(AVALUE, BVALUE, CVALUE, fxyz)
             ELSE
                CALL LSQUIT('determine_orbitals is not implemented for g and higher order orbitals',-1)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  end subroutine determine_orbitals

end module print_moorb_grid_mod


