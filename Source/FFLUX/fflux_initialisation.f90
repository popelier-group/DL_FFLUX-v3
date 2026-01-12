Subroutine fflux_initialisation

  Use kinds_f90
  Use omp_lib
  Use setup_module,  Only : nffluxdt, mxatms, pi, mxompl
  Use config_module, Only : natms,nlast
!  Use read_field,    Only : ntpmls ! number of different models
  Use site_module,   Only : numsit, ntpmls ! number of atoms in each model
  Use fflux_module,  Only : rt3,rt5,rt6,rt10,rt15,rt35,rt70, &
                          & rt2_3,rt3_8,rt5_8,rt5_12,rt1_24, &
                          & lxx,lyy,lzz,dCxx,dCyy,dCzz,pwd,conv_B2A,features, &
                          & maxFeatures,MultMomNames, &
                          & features_atom_order,dIQA_df, &
                          & conv_Ha_2_10Jmol, conv_Fau, conv_MM, &
                          & Quad_Ind, Oct_Ind, Hex_Ind, &
                          & Quad_Perm,Oct_Perm,Hex_Perm,ewaldFlag,dipole_deriv,&
                          & max_atms_model,Uni_Ind, &
                          & pForce, pEner, pxxx, pyyy, pzzz, Force_Max_Thresh, &
                          & Force_RMS_Thresh, Ener_Max_Thresh, &
                          & Dis_Max_Thresh, Dis_RMS_Thresh, CONVERGED, &
                          & fflux_convergence_criteria, &
                          & fflux_check_convergence, &
                          & set_Force_Max_Thresh, set_Force_RMS_Thresh, &
                          & set_Ener_Max_Thresh, set_Dis_Max_Thresh, &
                          & set_Dis_RMS_Thresh,kernel_names
  Use ewald_module,  Only : kprod

  Implicit None

  Character( Len=225 ) :: CurrentDir, QQQ

  Allocate( dipole_deriv(4) )
  Allocate( Uni_Ind(0:mxompl,0:mxompl,0:mxompl) )
  Allocate( kprod(0:mxompl+1,0:mxompl+1,0:mxompl+1) )
  Allocate( kernel_names(1:3) )

  ! get the current working directory, pwd
  Call getcwd(CurrentDir) 
  pwd = trim(adjustl(CurrentDir))

  ! assign the names of the multipoles

  MultMomNames = (/ "q00 ",                                                        &
                  & "q10 ","q11c","q11s",                                          &
                  & "q20 ","q21c","q21s","q22c","q22s",                            &
                  & "q30 ","q31c","q31s","q32c","q32s","q33c","q33s",              &
                  & "q40 ","q41c","q41s","q42c","q42s","q43c","q43s","q44c","q44s" /)

! List of implemented kernels. Can be extended as needed.
! Currently RBF and RBF-cyclic do the same thing.
  kernel_names(1) = 'rbf'
  kernel_names(2) = 'rbf-cyclic'
  kernel_names(3) = 'periodic'

  Call fflux_read_models

  ! explicit square root values needed in fflux_mm_explicittensor.f90

  rt3   = 1.7320508075689_wp
  rt5   = 2.2360679774998_wp
  rt6   = 2.4494897427832_wp
  rt10  = 3.1622776601684_wp
  rt15  = 3.8729833462074_wp
  rt35  = 5.9160797830996_wp
  rt70  = 8.3666002653408_wp
  rt2_3 = 0.81649658092773_wp
  rt3_8 = 0.61237243569579_wp
  rt5_8 = 0.79056941504209_wp
  rt5_12= 0.64549722436790_wp
  rt1_24= 0.20412414523193_wp

  ! Define the array to convert indices to MM names
  conv_MM(1 ,:) = (/ 0 , 0 /)

  conv_MM(2 ,:) = (/ 1 , 0 /)
  conv_MM(3 ,:) = (/ 1 , 1 /)
  conv_MM(4 ,:) = (/ 1 , 2 /)

  conv_MM(5 ,:) = (/ 2 , 0 /)
  conv_MM(6 ,:) = (/ 2 , 1 /)
  conv_MM(7 ,:) = (/ 2 , 2 /)
  conv_MM(8 ,:) = (/ 2 , 3 /)
  conv_MM(9 ,:) = (/ 2 , 4 /)

  conv_MM(10,:) = (/ 3 , 0 /)
  conv_MM(11,:) = (/ 3 , 1 /)
  conv_MM(12,:) = (/ 3 , 2 /)
  conv_MM(13,:) = (/ 3 , 3 /)
  conv_MM(14,:) = (/ 3 , 4 /)
  conv_MM(15,:) = (/ 3 , 5 /)
  conv_MM(16,:) = (/ 3 , 6 /)

  conv_MM(17,:) = (/ 4 , 0 /)
  conv_MM(18,:) = (/ 4 , 1 /)
  conv_MM(19,:) = (/ 4 , 2 /)
  conv_MM(20,:) = (/ 4 , 3 /)
  conv_MM(21,:) = (/ 4 , 4 /)
  conv_MM(22,:) = (/ 4 , 5 /)
  conv_MM(23,:) = (/ 4 , 6 /)
  conv_MM(24,:) = (/ 4 , 7 /)
  conv_MM(25,:) = (/ 4 , 8 /)

  Uni_Ind(0,0,0) = 1
  If (mxompl >= 1) Then
    Uni_Ind(1,0,0) = 2
    Uni_Ind(0,1,0) = 3
    Uni_Ind(0,0,1) = 4
  End If
  If (mxompl >= 2) Then
    Uni_Ind(2,0,0) = 5
    Uni_Ind(1,1,0) = 6
    Uni_Ind(1,0,1) = 7
    Uni_Ind(0,2,0) = 8
    Uni_Ind(0,1,1) = 9
    Uni_Ind(0,0,2) = 10
  End If
  If (mxompl >= 3) Then
    Uni_Ind(3,0,0) = 11
    Uni_Ind(2,1,0) = 12
    Uni_Ind(2,0,1) = 13
    Uni_Ind(1,2,0) = 14
    Uni_Ind(1,1,1) = 15
    Uni_Ind(1,0,2) = 16
    Uni_Ind(0,3,0) = 17
    Uni_Ind(0,2,1) = 18
    Uni_Ind(0,1,2) = 19
    Uni_Ind(0,0,3) = 20 
  End If
  If (mxompl == 4) Then
    Uni_Ind(4,0,0) = 21
    Uni_Ind(3,1,0) = 22
    Uni_Ind(3,0,1) = 23
    Uni_Ind(2,2,0) = 24
    Uni_Ind(2,1,1) = 25
    Uni_Ind(2,0,2) = 26
    Uni_Ind(1,3,0) = 27
    Uni_Ind(1,2,1) = 28
    Uni_Ind(1,1,2) = 29
    Uni_Ind(1,0,3) = 30
    Uni_Ind(0,4,0) = 31
    Uni_Ind(0,3,1) = 32
    Uni_Ind(0,2,2) = 33
    Uni_Ind(0,1,3) = 34
    Uni_Ind(0,0,4) = 35
  End If

! Define the array that stores the cartesian multipole indices - BCBS 04/2019

! Quadrupole
! Ind arrays for a rank n tensor convert the n indices into a single index.
! This is needed for rotating the quadrupole matrix in unidimensional form.
  Quad_Ind(1,:) = (/ 1, 1 /)
  Quad_Ind(2,:) = (/ 1, 2 /)
  Quad_Ind(3,:) = (/ 1, 3 /)
  Quad_Ind(4,:) = (/ 2, 2 /)
  Quad_Ind(5,:) = (/ 2, 3 /)
  Quad_Ind(6,:) = (/ 3, 3 /)

! Note: permutation arrays encode symmetry of tensors.
  Quad_Perm(1,:) = (/ 1, 2, 3 /)
  Quad_Perm(2,:) = (/ 2, 4, 5 /)
  Quad_Perm(3,:) = (/ 3, 5, 6 /)

! Octupole
  Oct_Ind(1,:)  = (/ 1, 1, 1 /)
  Oct_Ind(2,:)  = (/ 1, 1, 2 /)
  Oct_Ind(3,:)  = (/ 1, 1, 3 /)
  Oct_Ind(4,:)  = (/ 1, 2, 2 /)
  Oct_Ind(5,:)  = (/ 1, 2, 3 /)
  Oct_Ind(6,:)  = (/ 1, 3, 3 /)
  Oct_Ind(7,:)  = (/ 2, 2, 2 /)
  Oct_Ind(8,:)  = (/ 2, 2, 3 /)
  Oct_Ind(9,:)  = (/ 2, 3, 3 /)
  Oct_Ind(10,:) = (/ 3, 3, 3 /)

  Oct_Perm(1,1,:) = (/ 1, 2, 3  /)
  Oct_Perm(1,2,:) = (/ 2, 4, 5  /)
  Oct_Perm(1,3,:) = (/ 3, 5, 6  /)
  Oct_Perm(2,1,:) = (/ 2, 4, 5  /)
  Oct_Perm(2,2,:) = (/ 4, 7, 8  /)
  Oct_Perm(2,3,:) = (/ 5, 8, 9  /)
  Oct_Perm(3,1,:) = (/ 3, 5, 6  /)
  Oct_Perm(3,2,:) = (/ 5, 8, 9  /)
  Oct_Perm(3,3,:) = (/ 6, 9, 10 /)

! Hexadecapole
  Hex_Ind(1,:)  = (/ 1, 1, 1, 1 /)
  Hex_Ind(2,:)  = (/ 1, 1, 1, 2 /)
  Hex_Ind(3,:)  = (/ 1, 1, 1, 3 /)
  Hex_Ind(4,:)  = (/ 1, 1, 2, 2 /)
  Hex_Ind(5,:)  = (/ 1, 1, 2, 3 /)
  Hex_Ind(6,:)  = (/ 1, 1, 3, 3 /)
  Hex_Ind(7,:)  = (/ 1, 2, 2, 2 /)
  Hex_Ind(8,:)  = (/ 1, 2, 2, 3 /)
  Hex_Ind(9,:)  = (/ 1, 2, 3, 3 /)
  Hex_Ind(10,:) = (/ 1, 3, 3, 3 /)
  Hex_Ind(11,:) = (/ 2, 2, 2, 2 /)
  Hex_Ind(12,:) = (/ 2, 2, 2, 3 /)
  Hex_Ind(13,:) = (/ 2, 2, 3, 3 /)
  Hex_Ind(14,:) = (/ 2, 3, 3, 3 /)
  Hex_Ind(15,:) = (/ 3, 3, 3, 3 /)

  Hex_Perm(1,1,1,:) = (/ 1, 2, 3   /)
  Hex_Perm(1,1,2,:) = (/ 2, 4, 5   /)
  Hex_Perm(1,1,3,:) = (/ 3, 5, 6   /)
  Hex_Perm(1,2,1,:) = (/ 2, 4, 5   /)
  Hex_Perm(1,2,2,:) = (/ 4, 7, 8   /)
  Hex_Perm(1,2,3,:) = (/ 5, 8, 9   /)
  Hex_Perm(1,3,1,:) = (/ 3, 5, 6   /)
  Hex_Perm(1,3,2,:) = (/ 5, 8, 9   /)
  Hex_Perm(1,3,3,:) = (/ 6, 9, 10  /)
  Hex_Perm(2,1,1,:) = (/ 2, 4, 5   /)
  Hex_Perm(2,1,2,:) = (/ 4, 7, 8   /)
  Hex_Perm(2,1,3,:) = (/ 5, 8, 9   /)
  Hex_Perm(2,2,1,:) = (/ 4, 7, 8   /)
  Hex_Perm(2,2,2,:) = (/ 7, 11, 12 /)
  Hex_Perm(2,2,3,:) = (/ 8, 12, 13 /)
  Hex_Perm(2,3,1,:) = (/ 5, 8, 9   /)
  Hex_Perm(2,3,2,:) = (/ 8, 12, 13 /)
  Hex_Perm(2,3,3,:) = (/ 9, 13, 14 /)
  Hex_Perm(3,1,1,:) = (/ 3, 5, 6   /)
  Hex_Perm(3,1,2,:) = (/ 5, 8, 9   /)
  Hex_Perm(3,1,3,:) = (/ 6, 9, 10  /)
  Hex_Perm(3,2,1,:) = (/ 5, 8, 9   /)
  Hex_Perm(3,2,2,:) = (/ 8, 12, 13 /)
  Hex_Perm(3,2,3,:) = (/ 9, 13, 14 /)
  Hex_Perm(3,3,1,:) = (/ 6, 9, 10  /)
  Hex_Perm(3,3,2,:) = (/ 9, 13, 14 /)
  Hex_Perm(3,3,3,:) = (/ 10, 14, 15/)
  
  !Array for dipole derivative spherical to Cart indices.
  dipole_deriv(1) = 1
  dipole_deriv(2) = 4
  dipole_deriv(3) = 2
  dipole_deriv(4) = 3

  ! allocate the atom order in the features array

  Allocate( features_atom_order(natms,natms-1) )
  Allocate( features(mxatms,maxFeatures) )

  ! allocate the local axis orientation
  
  Allocate( lxx(mxatms,3) )
  Allocate( lyy(mxatms,3) )
  Allocate( lzz(mxatms,3) )
!Really these should be mxatms*3*3*3
  Allocate( dCxx(mxatms,3,3,3) )
  Allocate( dCyy(mxatms,3,3,3) )
  Allocate( dCzz(mxatms,3,3,3) )

  ! allocate the features and the derivatives

  !                  f^A    Omega,  A,  i
  !If (ewaldFlag .EQV. .FALSE.) Then
  !  Allocate( df_da(nFeatures,mxatms,mxatms,3) ) 
  !End If
  Allocate( dIQA_df(maxFeatures)  ) 

  ! allocate convergence check arrays and threshholds

  Allocate( pForce(natms) )
  Allocate( pEner(natms) )
  
  Allocate( pxxx(natms) )
  Allocate( pyyy(natms) )
  Allocate( pzzz(natms) )
  
  CONVERGED = .FALSE.
  
  pForce(:) = 0.0
  pEner(:) = 0.0
  
  pxxx(:) = 0.0
  pyyy(:) = 0.0
  pzzz(:) = 0.0

  Ener_Max_Thresh  = -1.0_wp
  Force_Max_Thresh = -1.0_wp
  Force_RMS_Thresh = -1.0_wp
  Dis_Max_Thresh   = -1.0_wp
  Dis_RMS_Thresh   = -1.0_wp

  If ( fflux_check_convergence .EQV. .TRUE.) Then
    If ( ( fflux_convergence_criteria .LT. 0.0_wp ) .AND. &
       & ( Force_Max_Thresh .LT. 0.0_wp )           .AND. &
       & ( Force_RMS_Thresh .LT. 0.0_wp )           .AND. &
       & ( Ener_Max_Thresh  .LT. 0.0_wp )           .AND. &
       & ( Dis_Max_Thresh   .LT. 0.0_wp )           .AND. &
       & ( Dis_RMS_Thresh   .LT. 0.0_wp )         ) Then
  
      ! Default Convergence Criteria
      fflux_convergence_criteria = 4
    End If
   
  !  convergence Criteria Taken From:
  !  http://www.psicode.org/psi4manual/master/optking.html  
  
    Select Case (fflux_convergence_criteria)
      Case (1) ! NWCHEM_LOOSE
        Force_Max_Thresh = 0.004500
        Force_RMS_Thresh = 0.003000
        Dis_Max_Thresh   = 0.005400
        Dis_RMS_Thresh   = 0.003600
      Case (2) ! GAU_LOOSE
        Force_Max_Thresh = 0.002500
        Force_RMS_Thresh = 0.001700
        Dis_Max_Thresh   = 0.010000
        Dis_RMS_Thresh   = 0.006700
      Case (3) ! TURBOMOLE
        Ener_Max_Thresh  = 0.000001
        Force_Max_Thresh = 0.001000
        Force_RMS_Thresh = 0.000500
        Dis_Max_Thresh   = 0.001000
        Dis_RMS_Thresh   = 0.000500
      Case (4) ! GAU
        Force_Max_Thresh = 0.000450
        Force_RMS_Thresh = 0.000300
        Dis_Max_Thresh   = 0.001800
        Dis_RMS_Thresh   = 0.001200
      Case (5) ! CFOUR
        Force_RMS_Thresh = 0.000100
      Case (6) ! QCHEM
        Ener_Max_Thresh  = 0.000001
        Force_Max_Thresh = 0.000300
        Dis_Max_Thresh   = 0.001200
      Case (7) ! MOLPRO
        Ener_Max_Thresh  = 0.000001
        Force_Max_Thresh = 0.000300
        Dis_Max_Thresh   = 0.000300
      Case (8) ! INTERFRAG_TIGHT
        Ener_Max_Thresh  = 0.000001
        Force_Max_Thresh = 0.000015
        Force_RMS_Thresh = 0.000010
        Dis_Max_Thresh   = 0.000600
        Dis_RMS_Thresh   = 0.000400
      Case (9) ! GAU_TIGHT
        Force_Max_Thresh = 0.000015
        Force_RMS_Thresh = 0.000010
        Dis_Max_Thresh   = 0.000060
        Dis_RMS_Thresh   = 0.000040
      Case (10) ! GAU_VERYTIGHT
        Force_Max_Thresh = 0.000002
        Force_RMS_Thresh = 0.000001
        Dis_Max_Thresh   = 0.000006
        Dis_RMS_Thresh   = 0.000004
    End Select
  
    If ( set_Ener_Max_Thresh .GT. 0.0_wp ) Then
      ener_Max_Thresh = set_Ener_Max_Thresh
    End If
  
    If ( set_Force_Max_Thresh .GT. 0.0_wp ) Then
      force_Max_Thresh = set_Force_Max_Thresh
    end If
    If ( set_Force_RMS_Thresh .GT. 0.0_wp ) Then
      force_RMS_Thresh = set_Force_Max_Thresh
    End If
  
    If ( set_Dis_Max_Thresh .GT. 0.0_wp ) Then
      dis_Max_Thresh = set_Dis_Max_Thresh
    End If
    If ( set_Dis_RMS_Thresh .GT. 0.0_wp ) Then
      dis_RMS_Thresh = set_Dis_Max_Thresh
    End If
  End If

  ! Bohr2Ang/Ang2Bohr Conversion Factor
  conv_B2A = 0.52917721067_wp
  !conv_B2A = 1.0_wp/1.88971616463_wp
  ! Hartrees to 10.J/mol
  conv_Ha_2_10Jmol = 262549.9638285d0 
  ! Force conversion term
  conv_Fau = conv_Ha_2_10Jmol / conv_B2A 

  write(nffluxdt,'(a)') '# step, E_IQA / Ha, E_vdW / kJ mol-1, E_coul / kJ mol-1, E_kin / kJ mol-1'

End Subroutine fflux_initialisation
