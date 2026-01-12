!A subroutine that converts multipole moments from spherical tensors into cartesian tensors.
!See A. J. Stone, theory of intermolecular forces, Appendix E for reference.
!Author - BCBS 03/2019
Subroutine Tensor_Conversion(limit,CartPoles,SphericalPoleArray)

  Use Kinds_f90
  Use setup_module
  Use fflux_module,  Only : PoleOrder,rt5,rt6,rt10,rt15,rt35,rt70, &
                          & rt2_3,rt3_8,rt5_8,rt5_12,rt1_24,maxFeatures
  Use config_module, Only : natms,nlast

  Integer,         Intent( In    ) :: limit
  Integer                          :: i, j
  Real( kind=wp ), Intent( In    ) :: SphericalPoleArray(limit,natms)
  Real( Kind=wp )                  :: CartPoles(mximpl,mxatms)
  Real( Kind=wp ), Allocatable     :: dQ_df(:)
  Real( Kind=wp )                  :: Q

  !Allocate( SphericalPoleArray(limit,natms) )
  Allocate( dQ_df(maxFeatures) )
 
  !SphericalPoleArray = 0.0_wp
  dQ_df = 0.0_wp 

  !Predict all the moments in spherical tensor form.
 !!$OMP  PARALLEL DEFAULT(NONE) &
 !!$OMP& PRIVATE(i,j,Q,dQ_df) SHARED(natms,limit,SphericalPoleArray)
  !!$OMP DO
!THIS PRED CAN BE MOVED TO FFLUX EWALD PROVIDED ROTATE GOES AFTER DERIV
!CONVERSION.
  !Do i=1, natms
  !  Do j=1, limit
!
!      Call fflux_predict_value(i, j+1, Q, dQ_df(:))
!      SphericalPoleArray(j,i) = Q
!
!    End Do
 ! End Do
  !!$OMP END DO
 !!$OMP END PARALLEL


  !Monopole and dipole need no conversion.
  !Simply put monopole into array and reorder dipole components.
  Do i=1, natms

    CartPoles(1,i) = SphericalPoleArray(1,i)
    CartPoles(2,i) = SphericalPoleArray(3,i)
    CartPoles(3,i) = SphericalPoleArray(4,i)
    CartPoles(4,i) = SphericalPoleArray(2,i)

  End Do

 !Convert the quadropole moments using explicit formulae.
 !The order of qudropole components is: xx,xy,xz,yy,yz,zz
  If (PoleOrder .GE. 3) Then

    Do i=1,natms

      CartPoles(5,i)  = (1.0_wp/3.0_wp) * ( 0.5_wp * ( (rt3*SphericalPoleArray(8,i)) - SphericalPoleArray(5,i) ))
      CartPoles(6,i)  = (1.0_wp/3.0_wp) * ( 0.5_wp * rt3 * SphericalPoleArray(9,i) )
      CartPoles(7,i)  = (1.0_wp/3.0_wp) * ( 0.5_wp * rt3 * SphericalPoleArray(6,i) )
      CartPoles(8,i)  = (1.0_wp/3.0_wp) * (-0.5_wp * ( (rt3*SphericalPoleArray(8,i)) + SphericalPoleArray(5,i) ))
      CartPoles(9,i)  = (1.0_wp/3.0_wp) * ( 0.5_wp * rt3 * SphericalPoleArray(7,i) )
      CartPoles(10,i) = (1.0_wp/3.0_wp) * SphericalPoleArray(5,i)

    End Do

  End If

! Convert octupole moments - order is:
! xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz
  If (PoleOrder .GE. 4) Then

    Do i=1,natms

      CartPoles(11,i) = (1.0_wp/15.0_wp)*( (rt5_8 *SphericalPoleArray(15,i)) - (rt3_8 *SphericalPoleArray(11,i)) ) !xxx
      CartPoles(12,i) = (1.0_wp/15.0_wp)*( (rt5_8 *SphericalPoleArray(16,i)) - (rt1_24*SphericalPoleArray(12,i)) ) !xxy
      CartPoles(13,i) = (1.0_wp/15.0_wp)*( (rt5_12*SphericalPoleArray(13,i)) - (0.5_wp*SphericalPoleArray(10,i)) ) !xxz
      CartPoles(14,i) =-(1.0_wp/15.0_wp)*( (rt5_8 *SphericalPoleArray(15,i)) + (rt1_24*SphericalPoleArray(11,i)) ) !xyy
      CartPoles(15,i) = (1.0_wp/15.0_wp)*(  rt5_12*SphericalPoleArray(14,i))                                       !xyz
      CartPoles(16,i) = (1.0_wp/15.0_wp)*(  rt2_3 *SphericalPoleArray(11,i))                                       !xzz
      CartPoles(17,i) =-(1.0_wp/15.0_wp)*( (rt5_8 *SphericalPoleArray(16,i)) + (rt3_8 *SphericalPoleArray(12,i)) ) !yyy
      CartPoles(18,i) =-(1.0_wp/15.0_wp)*( (rt5_12*SphericalPoleArray(13,i)) + (0.5_wp*SphericalPoleArray(10,i)) ) !yyz
      CartPoles(19,i) = (1.0_wp/15.0_wp)*(  rt2_3 *SphericalPoleArray(12,i))                                       !yzz
      CartPoles(20,i) = (1.0_wp/15.0_wp)*SphericalPoleArray(10,i)                                                  !zzz

    End Do

  End If

! Convert the hexadecapole moments - order is:
! xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz
  If (PoleOrder .GE. 5) Then

    Do i=1, natms

      CartPoles(21,i) = ( ((3.0_wp/8.0_wp)*SphericalPoleArray(17,i)) - (0.25_wp*rt5*SphericalPoleArray(20,i)) + ((rt35*SphericalPoleArray(24,i))/8.0_wp) )
      CartPoles(22,i) = ( 0.125_wp*((rt35*SphericalPoleArray(25,i)) - (rt5*SphericalPoleArray(21,i))) )
      CartPoles(23,i) = ( 0.0625_wp*((rt70*SphericalPoleArray(22,i)) - (3.0_wp*rt10*SphericalPoleArray(18,i))) )
      CartPoles(24,i) = ( (0.125_wp*SphericalPoleArray(17,i)) - (0.125_wp*rt35*SphericalPoleArray(24,i)) )
      CartPoles(25,i) = ( 0.0625_wp*((rt70*SphericalPoleArray(23,i)) - (rt10*SphericalPoleArray(19,i))) )
      CartPoles(26,i) = ( 0.5_wp*((0.5_wp*rt5*SphericalPoleArray(20,i)) - SphericalPoleArray(17,i)) )
      CartPoles(27,i) = (-0.125_wp*((rt5*SphericalPoleArray(21,i)) + (rt35*SphericalPoleArray(25,i))) )
      CartPoles(28,i) = (-0.0625_wp*((rt10*SphericalPoleArray(18,i)) + (rt70*SphericalPoleArray(22,i))) )
      CartPoles(29,i) = ( 0.25_wp*rt5*SphericalPoleArray(21,i) )
      CartPoles(30,i) = ( rt5_8*SphericalPoleArray(18,i) )
      CartPoles(31,i) = ( ((3.0_wp/8.0_wp)*SphericalPoleArray(17,i)) + (0.25_wp*rt5*SphericalPoleArray(20,i)) + ((rt35*SphericalPoleArray(24,i))/8.0_wp) )
      CartPoles(32,i) = (-0.0625_wp*((3.0_wp*rt10*SphericalPoleArray(19,i)) + (rt70*SphericalPoleArray(23,i))) )
      CartPoles(33,i) = (-0.5_wp*((0.5_wp*rt5*SphericalPoleArray(20,i)) + SphericalPoleArray(17,i)) )
      CartPoles(34,i) = ( rt5_8*SphericalPoleArray(19,i) )
      CartPoles(35,i) = SphericalPoleArray(17,i)

      Do j=21,35
        CartPoles(j,i) = (1.0_wp/105.0_wp) * CartPoles(j,i)
      End Do

    End Do

  End If

!  Deallocate( SphericalPoleArray )
  Deallocate( dQ_df )

End Subroutine
