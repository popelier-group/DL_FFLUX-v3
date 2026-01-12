!File that rotates the multipole moments from the ALF to the global frame for use in multipolar Ewald
!Author BCBS - 12/2018

Subroutine fflux_rotate_moments(nstep,CartPoles,SphericalPoles,arrSize)

  Use kinds_f90
  Use omp_lib
  Use setup_module
  Use config_module, Only : natms,nlast,xxx,yyy,zzz,ixyz,ltg,cell
  Use comms_module, Only: idnode
  Use fflux_module
  !Use mpoles_module, Only : mplgfr,mprotx,mproty,mprotz
  Use mpoles_module

  Implicit None

  Integer, Intent( In    )       :: nstep,arrSize
  Real( Kind = wp), Intent( InOut ) :: CartPoles(mximpl,mxatms),SphericalPoles(arrSize,natms)
  Integer                        :: i,j,k,l,a,b,c,d,limit,mlast
  Logical                        :: Bool
  Integer, Allocatable           :: ixyz0(:)
  Real( Kind = wp )              :: Q,T_det,temp
  Real( Kind = wp ), Allocatable :: C_Mat(:,:),InvC_Mat(:,:),dipole(:), &
                                    & Quadrupole(:,:),Temp_Quadrupole(:,:)

  Allocate(C_Mat(3,3))
  Allocate(InvC_Mat(3,3))
  Allocate(dipole(3))
  Allocate(Quadrupole(3,3))
  Allocate(Temp_Quadrupole(3,3))

  C_Mat = 0.0_wp
  InvC_Mat = 0.0_wp
  dipole = 0.0_wp
  temp = 0.0_wp
  !CartPoles = 0.0_wp

 If (PoleOrder .EQ. 1) Then
   limit = 1
 Else If (PoleOrder .EQ. 2) Then
   limit = 4
 Else If (PoleOrder .EQ. 3) Then
   limit = 9
 Else If (PoleOrder .EQ. 4) Then
   limit = 16
 Else If (PoleOrder .EQ. 5) Then
   limit = 25
 End If

 Call Tensor_Conversion(limit,CartPoles,SphericalPoles)
  
  !Rotate all the non-monopole moments from ALF to global frame
  Do i=1,natms

    !Form the C matrix for atom, this contains the unit vectors of the ALF.
    !C_Mat = Form_C_Matrix(i)
    !This is significantly better but still technically needless copying.
    C_Mat(1,:) = lxx(i,:)
    C_Mat(2,:) = lyy(i,:)
    C_Mat(3,:) = lzz(i,:)

    !Transpose the C matrix to get inverse (it is unitary).
    Do j=1,3
      Do k=1,3

        InvC_Mat(j,k) = C_Mat(k,j)

      End Do
    End Do

    !Rotate the moments for this atom and put into mplgfr array.
    Do j=2,4

      dipole(j-1) = ( CartPoles(j,i) * Conv_B2A )

    End Do

    !Add monopole to array.
    mplgfr(1,i) = CartPoles(1,i)

    !If ((ldipole .EQV. .TRUE.) .AND. (nstep .GE. nsdipole) .AND. (modulo(nstep,isdipole) .EQ. 0)) Then
    !  Write(ndipoledt, "(i10,3f20.10)") ltg(i), xxx(i),yyy(i),zzz(i)
    !  Write(ndipoledt, "(3f20.10)") mplgfr(1,i), cell(1)
    !  Write(ndipoledt, "(3f20.10)") dipole(1),dipole(2),dipole(3)
    !End If

    !Rotate dipole and add to array.
    Do j=1,3

      mplgfr(j+1,i) = ( (InvC_Mat(j,1) * dipole(1)) + (InvC_Mat(j,2) * dipole(2)) + (InvC_Mat(j,3) * dipole(3)) )

    End Do

    !These writes are not MPI safe right now.
    If ((ldipole .EQV. .TRUE.) .AND. (nstep .GE. nsdipole) .AND. (modulo(nstep,isdipole) .EQ. 0)) Then
      Write(ndipoledt, "(i10,3f20.10)") ltg(i), xxx(i),yyy(i),zzz(i)
      Write(ndipoledt, "(3f20.10)") mplgfr(1,i), cell(1)
      Write(ndipoledt, "(3f20.10)") mplgfr(2,i),mplgfr(3,i),mplgfr(4,i)
    End If

    !Rotate quadropole moments.
    If (PoleOrder .GE. 3) Then

      !New rotation method
      Do j=1,6
        temp = 0.0_wp
        Do a=1,3
          Do b=1,3

            temp = temp + InvC_Mat(Quad_Ind(j,1),a)*InvC_Mat(Quad_Ind(j,2),b)*CartPoles(Quad_Perm(a,b)+4,i)

          End Do
        End Do
        mplgfr(j+4,i) = ( temp * Conv_B2A * Conv_B2A )!Don't even need this copying, can go straight in.
      End Do 

      !Account for degeneracy in DL_POLY unidimensional multipole treatment.
      mplgfr(6,i) = mplgfr(6,i) * 2.0_wp
      mplgfr(7,i) = mplgfr(7,i) * 2.0_wp
      mplgfr(9,i) = mplgfr(9,i) * 2.0_wp

    End If

    If (PoleOrder .GE. 4) Then

      Do j=1,10
        temp = 0.0_wp
        Do a=1,3
          Do b=1,3
            Do c=1,3
                                                                                                                     
              temp = temp + InvC_Mat(Oct_Ind(j,1),a)*InvC_Mat(Oct_Ind(j,2),b)*InvC_Mat(Oct_Ind(j,3),c)*CartPoles(Oct_Perm(a,b,c)+10,i)
                                                                                                                     
            End Do
          End Do
        End Do
        mplgfr(j+10,i) = ( temp * Conv_B2A * Conv_B2A * Conv_B2A )!Don't even need this copying, can go straight in.
      End Do

      !Account for degeneracy in DL_POLY unidimensional multipole treatment.
      mplgfr(12,i) = mplgfr(12,i) * 3.0_wp
      mplgfr(13,i) = mplgfr(13,i) * 3.0_wp
      mplgfr(14,i) = mplgfr(14,i) * 3.0_wp
      mplgfr(15,i) = mplgfr(15,i) * 6.0_wp
      mplgfr(16,i) = mplgfr(16,i) * 3.0_wp
      mplgfr(18,i) = mplgfr(18,i) * 3.0_wp
      mplgfr(19,i) = mplgfr(19,i) * 3.0_wp 

    End If

    If (PoleOrder .GE. 5) Then

      Do j=1,15
        temp = 0.0_wp
        Do a=1,3
          Do b=1,3
            Do c=1,3
              Do d=1,3

                temp = temp + InvC_Mat(Hex_Ind(j,1),a)*InvC_Mat(Hex_Ind(j,2),b)*InvC_Mat(Hex_Ind(j,3),c)*InvC_Mat(Hex_Ind(j,4),d)*CartPoles(Hex_Perm(a,b,c,d)+20,i)

              End Do
            End Do
          End Do
        End Do
        mplgfr(j+20,i) = ( temp * Conv_B2A * Conv_B2A * Conv_B2A * Conv_B2A)
      End Do

      !Account for degeneracy in DL_POLY unidimensional multipole treatment.
      mplgfr(22,i) = mplgfr(22,i) * 4.0_wp
      mplgfr(23,i) = mplgfr(23,i) * 4.0_wp
      mplgfr(24,i) = mplgfr(24,i) * 6.0_wp
      mplgfr(25,i) = mplgfr(25,i) * 12.0_wp
      mplgfr(26,i) = mplgfr(26,i) * 6.0_wp
      mplgfr(27,i) = mplgfr(27,i) * 4.0_wp
      mplgfr(28,i) = mplgfr(28,i) * 12.0_wp 
      mplgfr(29,i) = mplgfr(29,i) * 12.0_wp
      mplgfr(30,i) = mplgfr(30,i) * 4.0_wp
      mplgfr(32,i) = mplgfr(32,i) * 4.0_wp
      mplgfr(33,i) = mplgfr(33,i) * 6.0_wp
      mplgfr(34,i) = mplgfr(34,i) * 4.0_wp

    End If

    !Generate infinitesimal rotations of multipoles - needed for torques
    Call infinitesimal_rotation(i)

  End Do

! Unwrap the multipole and infinitesimal rotation arrays.
! Note: this code is taken from the mpoles_rotmat_set_halo

  Allocate(ixyz0(1:mxatms))
  ixyz0(1:nlast) = ixyz(1:nlast)
  mlast=natms

  Call mpoles_rotmat_export(-1,mlast,ixyz0,CartPoles)
  Call mpoles_rotmat_export( 1,mlast,ixyz0,CartPoles)

  Call mpoles_rotmat_export(-2,mlast,ixyz0,CartPoles)
  Call mpoles_rotmat_export( 2,mlast,ixyz0,CartPoles)

  Call mpoles_rotmat_export(-3,mlast,ixyz0,CartPoles)
  Call mpoles_rotmat_export( 3,mlast,ixyz0,CartPoles)
 
  Deallocate(ixyz0)
  Deallocate(C_Mat)
  Deallocate(InvC_Mat)
  Deallocate(dipole)
  Deallocate(Quadrupole)
  Deallocate(Temp_Quadrupole)

End subroutine
