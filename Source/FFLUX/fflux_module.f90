Module fflux_module

!==============================================================================
! FFlux Version 1.0
! 
! Author : Joseph CR Thacker
! Date   : 2017      
!
! All equation numbers given in FFlux files (unless otherwise stated) are 
! from the paper:
!
! "Mills, MJL & Popelier, PLA, J. Chem. Theory. Comput., 2014, 10, 3840-3856"
! 
! and its Supplementary Material.
!
! Want to find an equation? (for example, equation 18). Try in terminal:
!
! $ grep -ir "[[18]]" *.f90
!
!==============================================================================
! There are plenty of code optimisiations that can be done within these
! functions - many values are calculated multiple times. - JCRT
!==============================================================================

  Use kinds_f90  
  Use setup_module,           Only : nffluxdt

  Implicit none

! Constants
  Real( Kind = wp )               :: conv_B2A, conv_Ha_2_10Jmol, conv_Fau

! File Management
  Character(Len=:),  Allocatable :: pwd
  Character(Len=4)               :: MultMomNames(25)  
  
! Kriging Variables
  Integer,           Allocatable :: alf(:,:), globalAlf(:,:), features_atom_order(:,:), &
                                  & atm_label_alias(:),non_alf_atms(:,:),dipole_deriv(:), &
                                  & mpi_list(:)

  Integer                        :: totAtms,maxFeatures,max_atms_model
  
  Real( Kind = wp ), Allocatable :: features(:,:),df_da(:,:,:,:),dIQA_df(:)

  Character(Len=64), Allocatable :: model_names(:)
  Character(Len=64), Allocatable :: kernel_names(:)

! Derived data type for Kriging models - BCBS 5/20
  Type Model

    Character(Len=128) :: SystemName,AtomName
    Integer                       :: nTrain,nFeatures,nPredPerAtm,natms_model,nModels
    Integer,          Allocatable :: ALF(:),non_ALF(:),kernels(:)
    Real( Kind=wp),   Allocatable :: trained_value(:,:,:),weights_krig(:,:),mu(:), &
                                   & Theta(:,:),p_krig(:,:)
    !Procedure(KernelType), Pointer :: Kernel - Placeholder until KernelType
    !implemented.

  End Type Model

! Array of Model types - Size will be specified in KRIGING.txt
  Type ( Model ), Target,  Allocatable :: Models(:)

! Derived data type for pointers to the Model type - 5/20 BCBS
! Needed to construct an array of pointers to Model types
  Type ModelPtr

    Type ( Model ), Pointer :: ptr

  End Type ModelPtr

! Array of pointers that keeps track of which model describes a given atom.
  Type ( ModelPtr ), Allocatable :: model_list(:)

! Local Cartesian Axis Orientation (defined by the ALF)
!ER : these variables are essential to allocate memory, and therefore avoid some corrupted values we observed before.
  Real( Kind = wp ), Allocatable, Save :: lxx(:,:), lyy(:,:), lzz(:,:),& ! ALF Rotation Matrix split into 3. Called C in the paper
                                        & dCxx(:,:,:,:), dCyy(:,:,:,:), dCzz(:,:,:,:) 
  Real( Kind = wp )                    :: y_vec(3), sigma_fflux

! Variables for fflux_mm_explicittensor.f90
  Integer, Save :: conv_MM(25,2) 
  Real( Kind = wp )               :: R_
  Real( Kind = wp )               :: rax,ray,raz, &
                                        & rbx,rby,rbz, &
                                        & Cxx,Cyx,Czx, &
                                        & Cxy,Cyy,Czy, &
                                        & Cxz,Cyz,Czz
  Real( Kind = wp )               :: rt3, rt5, rt6, rt10, rt15, rt35, rt70, &
                                   & rt2_3,rt3_8,rt5_8,rt5_12,rt1_24
  Real( Kind = wp ), Allocatable  :: TAB(:), dT_dq(:,:)

! Ewald variables  -BCBS 12/18
  Logical                         :: ewaldFlag,clustFlag
  Integer                         :: poleOrder
! Variable cut of radii variables - BCBS 5/20
  Real                            :: cutDipole,cutQuad
! Electrostatic convergence test variables - BCBS 8/20
  Logical                         :: lEleConv
! FFLUX print option variables - BCBS 03/2019
  Logical                         :: lfflux,ldipole,lfflux_energy,lfflux_force
  Integer                         :: nsfflux,isfflux,nsdipole,isdipole !Note DL_POLY notation ns = starting timestep, is = timestep interval.

! Variables to store indices for rotation of Cartesian multipole moments - BCBS 04/2019
  Integer, Save                   :: Quad_Ind(6,2),Oct_Ind(10,3),Hex_Ind(15,4), &
                                   & Quad_Perm(3,3),Oct_Perm(3,3,3),Hex_Perm(3,3,3,3)
  Integer, Allocatable            :: Uni_Ind(:,:,:)

! Variables for normalisation
  Real( Kind = wp ), Allocatable  :: training_min_value(:,:), training_max_value(:,:)
  Logical                         :: fflux_normalise_training = .FALSE.

! Variables for fflux_converged.f90
  Real( Kind = wp )               :: Force_Max_Thresh, Force_RMS_Thresh, &
                                   & Ener_Max_Thresh, &
                                   & Dis_Max_Thresh, Dis_RMS_Thresh
  Real( Kind = wp ), Allocatable  :: pForce(:), pEner(:), pxxx(:), pyyy(:), pzzz(:)

! Check convergence (set in read_control.f90)
  Integer                         :: fflux_convergence_criteria = -1
  Logical                         :: fflux_check_convergence = .FALSE.

  Logical                         :: CONVERGED

  Real( Kind = wp )               :: set_Force_Max_Thresh = -1.0_wp, &
                                   & set_Force_RMS_Thresh = -1.0_wp, &
                                   & set_Ener_Max_Thresh  = -1.0_wp, &
                                   & set_Dis_Max_Thresh   = -1.0_wp, &
                                   & set_Dis_RMS_Thresh   = -1.0_wp

  !$OMP THREADPRIVATE(sigma_fflux,y_vec)

  !=========================================  
  !=========================================  
  ! FFlux Functions
  !=========================================  
  !=========================================  

  Contains

  !=========================================  
  ! sign_j                            [[12]]
  !========================================= 

  Real( Kind = wp ) Function sign_j(fdiff)

    Implicit None

    Real( Kind = wp ), Intent( In ) :: fdiff

    If (fdiff .LE. 0.0_wp) Then
      sign_j = 1.0d0
    Else
      sign_j = (-1.0d0)
    End If
 
  End Function

  !=========================================  
  ! Kronecker delta {delta_{ik}}
  !========================================= 

  Real( Kind = wp ) Function kronecker(i,k)

    Implicit None

    Integer, Intent( In ) :: i, k

    If (i .EQ. k) Then
      kronecker = 1.0d0
    Else
      kronecker = 0.0d0
    End If
 
  End Function

  !=========================================  
  ! Delta_{Omega,Ax} Delta_{Omega,Axy} Delta_{Omega,An}    [[33b]]
  !=========================================  

  Function delta_Ax(iatm, Omega)

    Use kinds_f90

    Implicit None

    Integer, Intent( In ) :: iatm, Omega

    Real( Kind = wp ) :: delta_Ax

    If      (Omega .EQ. iatm       ) Then
      delta_Ax = (-1.0_wp)
    Else If (Omega .EQ. alf(iatm,2)) Then
      delta_Ax = (1.0_wp)
    Else
      delta_Ax = 0.0_wp
    End If

  End Function

  Function delta_Axy(iatm, Omega)

    Use kinds_f90

    Implicit None

    Integer, Intent( In ) :: iatm, Omega
    Real( Kind = wp ) :: delta_Axy


    If      (Omega .EQ. iatm       ) Then
      delta_Axy = (-1.0_wp)
    Else If (Omega .EQ. alf(iatm,3)) Then
      delta_Axy = (1.0_wp)
    Else
      delta_Axy = 0.0_wp
    End If

  End Function
  
  Function delta_An(iatm, n, Omega)

    Use kinds_f90

    Implicit None

    Integer, Intent( In ) :: iatm, n, Omega
    Real( Kind = wp ) :: delta_An

    If      (Omega .EQ. iatm       ) Then
      delta_An = (-1.0_wp)
    Else If (Omega .EQ. n) Then
      delta_An = (1.0_wp)
    Else
      delta_An = 0.0_wp
    End If

  End Function

  !=========================================  
  ! C_1                              [[25]]
  !=========================================  
 !Changed by BCBS - 02/2019 for optimisation purposes. Function now outputs a 
 !vector rather than scalar, reduces unnecessary work.
  Function C_1( atm )

    Use kinds_f90
    Use config_module, Only : imcon,cell,xxx,yyy,zzz

    Implicit None
    
    Integer,           Intent( In ) :: atm
    Real( Kind = wp )               :: C_1(3),reciprocalFeature,xdiff,ydiff,zdiff

    reciprocalFeature = 1.0_wp/features(atm,1)

    xdiff = xxx(alf(atm,2)) - xxx(atm)
    ydiff = yyy(alf(atm,2)) - yyy(atm)
    zdiff = zzz(alf(atm,2)) - zzz(atm)

    Call images_s(imcon, cell, xdiff, ydiff, zdiff)

    C_1(1) = ( xdiff * reciprocalFeature )
    C_1(2) = ( ydiff * reciprocalFeature )
    C_1(3) = ( zdiff * reciprocalFeature )

  End Function

  !=========================================  
  ! C_2                              [[27]]
  !=========================================  
 !Changed by BCBS - 02/2019 for optimisation purposes. Function now outputs a
 !vector rather than scalar, reduces unnecessary work.
  Function C_2( atm )
 
    Use kinds_f90
    Use config_module, Only : imcon,cell,xxx,yyy,zzz
 
    Implicit None

    Integer,           Intent( In ) :: atm

    Real( Kind = wp ) :: xdiff1, ydiff1, zdiff1, xdiff2, ydiff2, zdiff2, reciprocal, C_2(3)

    xdiff1 = xxx(alf(atm,2)) - xxx(atm)
    ydiff1 = yyy(alf(atm,2)) - yyy(atm)
    zdiff1 = zzz(alf(atm,2)) - zzz(atm)

    Call images_s(imcon, cell, xdiff1, ydiff1, zdiff1)

    xdiff2 = xxx(alf(atm,3)) - xxx(atm)
    ydiff2 = yyy(alf(atm,3)) - yyy(atm)
    zdiff2 = zzz(alf(atm,3)) - zzz(atm)

    Call images_s(imcon, cell, xdiff2, ydiff2, zdiff2)

    sigma_fflux = - ( (xdiff1 * xdiff2) + (ydiff1 * ydiff2) + (zdiff1 * zdiff2) ) &
                & / ( (xdiff1 * xdiff1) + (ydiff1 * ydiff1) + (zdiff1 * zdiff1) )
          
    y_vec(1) = (sigma_fflux * xdiff1) + xdiff2
    y_vec(2) = (sigma_fflux * ydiff1) + ydiff2
    y_vec(3) = (sigma_fflux * zdiff1) + zdiff2

    reciprocal = 1.0_wp/(dsqrt( (y_vec(1) * y_vec(1)) + (y_vec(2) * y_vec(2)) + (y_vec(3) * y_vec(3)) ))

    C_2(:) = y_vec(:) * reciprocal

  End Function

  !=========================================  
  ! C_3                              [[30]]
  !=========================================  

  Function C_3(atm, C1_Vec, C2_Vec)

    Implicit None

    Integer,           Intent( In ) :: atm
    Real( Kind = wp ), Intent( In ) :: C1_Vec(3),C2_Vec(3)
    Real( Kind = wp ), Dimension(3) :: C_3

    C_3(1) = ( C1_Vec(2) * C2_Vec(3) ) - ( C1_Vec(3) * C2_Vec(2) )
    C_3(2) = ( C1_Vec(3) * C2_Vec(1) ) - ( C1_Vec(1) * C2_Vec(3) )
    C_3(3) = ( C1_Vec(1) * C2_Vec(2) ) - ( C1_Vec(2) * C2_Vec(1) )

  End Function
  
  !=========================================  
  ! [[B15]] - my rederived form [[B15]] equiv.
  !=========================================  
  !Altered by BCBS 03/2019 to account for PBC's

  Function dRaxdotRaxy_da( iatm, Omega, xdiff1,ydiff1,zdiff1,xdiff2,ydiff2,zdiff2 )
    Use kinds_f90
    Use config_module, Only : imcon,cell

    Implicit None

    Integer,           Intent( In ) :: iatm, Omega
    Real( Kind = wp ), Intent( In ) :: xdiff1,ydiff1,zdiff1,xdiff2,ydiff2,zdiff2

    Real( Kind = wp ) :: dRaxdotRaxy_da(3)

    If ( Omega .EQ. iatm ) Then

      dRaxdotRaxy_da(1) = -1.0_wp*(xdiff1 + xdiff2)
      dRaxdotRaxy_da(2) = -1.0_wp*(ydiff1 + ydiff2)
      dRaxdotRaxy_da(3) = -1.0_wp*(zdiff1 + zdiff2)

    Else If ( Omega .EQ. alf(iatm,2) ) Then

      dRaxdotRaxy_da(1) = xdiff2
      dRaxdotRaxy_da(2) = ydiff2
      dRaxdotRaxy_da(3) = zdiff2

    Else If ( Omega .EQ. alf(iatm,3) ) Then

      dRaxdotRaxy_da(1) = xdiff1
      dRaxdotRaxy_da(2) = ydiff1
      dRaxdotRaxy_da(3) = zdiff1

    Else

      dRaxdotRaxy_da(1) = 0.0_wp 
      dRaxdotRaxy_da(2) = 0.0_wp 
      dRaxdotRaxy_da(3) = 0.0_wp 

    End If 

  End Function
  
  !=========================================  
  ! [[B15]] - my rederived form [[B15]] equiv.
  !=========================================  
  !Altered by BCBS 03/2019 to account for PBC's

  Function dRaxdotRax_da( iatm, Omega, xdiff1,ydiff1,zdiff1 ) ! d(Rax . Rax) / d(alpha)

    Use kinds_f90

    Implicit None

    Integer,           Intent( In ) :: iatm, Omega
    Real( Kind = wp ), Intent( In ) :: xdiff1,ydiff1,zdiff1

    Real( Kind = wp ) :: dRaxdotRax_da(3)

    If ( Omega .EQ. iatm ) Then

      dRaxdotRax_da(1) = -2.0_wp * xdiff1
      dRaxdotRax_da(2) = -2.0_wp * ydiff1
      dRaxdotRax_da(3) = -2.0_wp * zdiff1

    Else If ( Omega .EQ. alf(iatm,2) ) Then

      dRaxdotRax_da(1) = 2.0_wp * xdiff1
      dRaxdotRax_da(2) = 2.0_wp * ydiff1
      dRaxdotRax_da(3) = 2.0_wp * zdiff1

    Else

      dRaxdotRax_da(1) = 0.0_wp
      dRaxdotRax_da(2) = 0.0_wp
      dRaxdotRax_da(3) = 0.0_wp

    End If 

  End Function

  !=========================================  
  ! d(sigma)/d(alpha_i)              [[C13]]
  !=========================================

  Function dsigma_da(iatm, Omega)

    Use kinds_f90
    Use config_module, Only : imcon,cell,xxx,yyy,zzz

    Implicit None

    Integer,           Intent( In ) :: iatm, Omega
    
    Real( Kind = wp ) :: dsigma_da(3), RaxdotRax, RaxdotRaxy
    Real( Kind = wp ) :: xdiff1, ydiff1, zdiff1, xdiff2, ydiff2, zdiff2

    xdiff1 = xxx(alf(iatm,2)) - xxx(iatm)
    ydiff1 = yyy(alf(iatm,2)) - yyy(iatm)
    zdiff1 = zzz(alf(iatm,2)) - zzz(iatm)

    Call images_s(imcon, cell, xdiff1, ydiff1, zdiff1)

    xdiff2 = xxx(alf(iatm,3)) - xxx(iatm)
    ydiff2 = yyy(alf(iatm,3)) - yyy(iatm)
    zdiff2 = zzz(alf(iatm,3)) - zzz(iatm)

    Call images_s(imcon, cell, xdiff2, ydiff2, zdiff2)

    RaxdotRax  = (xdiff1*xdiff1) + (ydiff1*ydiff1) + (zdiff1*zdiff1) 
    RaxdotRaxy = (xdiff1*xdiff2) + (ydiff1*ydiff2) + (zdiff1*zdiff2) 

    dsigma_da(:) = ( RaxdotRaxy * dRaxdotRax_da(iatm, Omega, xdiff1,ydiff1,zdiff1) / (RaxdotRax * RaxdotRax)) &
                 & - ( dRaxdotRaxy_da(iatm, Omega, xdiff1,ydiff1,zdiff1,xdiff2,ydiff2,zdiff2) / RaxdotRax ) 

  End Function


  !=========================================  
  ! d(y_j)/d(alpha_i)                [[C11]]
  !=========================================

  Function dyj_da(iatm, Omega, diff)

    Use kinds_f90

    Implicit None

    Integer,           Intent( In ) :: iatm, Omega
    Real( Kind = wp ), Intent( In ) :: diff(:)
                              ! i k
    Real( Kind = wp ) :: dyj_da(3,3), dsigma(3), delta_OmAx, delta_OmAxy
    Integer :: i, k

    dsigma(:) = dsigma_da(iatm, Omega) 
    delta_OmAx = delta_Ax(iatm,Omega)
    delta_OmAxy= delta_Axy(iatm,Omega)
  
    Do i = 1,3
      Do k = 1,3
        dyj_da(i,k) = (dsigma(i)*diff(k)) + (sigma_fflux*kronecker(i,k)*delta_OmAx) + (kronecker(i,k)*delta_OmAxy)
      End Do
    End Do

  End Function

  !=========================================  
  ! d(1/sqrt(y.y))/d(alpha_i)        [[C12]]
  !========================================= 

  Function EqC12(iatm, Omega, diff )
  
    Use kinds_f90

    Implicit None

    Integer,           Intent( In ) :: iatm, Omega
    Real( Kind = wp ), Intent( In ) :: diff(:)

    Integer :: j
    Real( Kind = wp ) :: EqC12(3), dyj(3,3), sqrtydoty, temp

    EqC12(:) = 0.0_wp
    dyj(:,:) = dyj_da(iatm, Omega, diff)

    Do j = 1,3

      EqC12(1) = EqC12(1) + (y_vec(j) * dyj(1,j)) 
      EqC12(2) = EqC12(2) + (y_vec(j) * dyj(2,j))
      EqC12(3) = EqC12(3) + (y_vec(j) * dyj(3,j))

    End Do

    sqrtydoty = dsqrt( (y_vec(1)*y_vec(1)) + (y_vec(2)*y_vec(2)) + (y_vec(3)*y_vec(3)) ) 

    temp = (-1.0_wp) / (sqrtydoty*sqrtydoty*sqrtydoty)  

    EqC12(:) = EqC12(:) * temp
   
  End Function 

  !=========================================  
  ! d(R)/d(alpha_i) for all i         [[18]]
  !=========================================  

  Function dR_da(iatm, jatm, cFeat, Omega)
    Use kinds_f90
    Use config_module, Only : xxx,yyy,zzz,imcon,cell

    Implicit None

    Integer,           Intent( In ) :: iatm, jatm, cFeat, Omega

    Integer :: indx
    Real( Kind =  wp ) :: dR_da(3),xdiff,ydiff,zdiff,reci_feat

    If ( Omega .EQ. iatm ) Then
      
      xdiff = xxx(Omega) - xxx(jatm)
      ydiff = yyy(Omega) - yyy(jatm)
      zdiff = zzz(Omega) - zzz(jatm)

      Call images_s(imcon, cell, xdiff, ydiff, zdiff)

      reci_feat = 1.0_wp/features(iatm,cFeat)

      dR_da(1) = xdiff * reci_feat
      dR_da(2) = ydiff * reci_feat
      dR_da(3) = zdiff * reci_feat

    Else If ( Omega .EQ. jatm ) Then 

       xdiff = xxx(Omega) - xxx(iatm)                      
       ydiff = yyy(Omega) - yyy(iatm)
       zdiff = zzz(Omega) - zzz(iatm)
 
       Call images_s(imcon, cell, xdiff, ydiff, zdiff)
 
       reci_feat = 1.0_wp/features(iatm,cFeat)
 
       dR_da(1) = xdiff * reci_feat
       dR_da(2) = ydiff * reci_feat
       dR_da(3) = zdiff * reci_feat

    Else

      dR_da(1) = 0.0_wp
      dR_da(2) = 0.0_wp
      dR_da(3) = 0.0_wp

    End If

  End Function

  !=========================================  
  ! d(C_1k)/d(alpha_i)      [[C6]] or [[26]]
  !=========================================  

  ! Consider rewriting this equation as calculated (can use if statements for the two different Omega values)

  Function dC1k_da(iatm, Omega, C1_Vec)

    Use kinds_f90
    Use config_module, Only : imcon,cell,xxx,yyy,zzz

    Implicit None
    
    Integer,           Intent( In ) :: iatm, Omega
    Real( Kind = wp ), Intent( In ) :: C1_Vec(3)
  
    Real( Kind = wp ) :: dC1k_da(3,3)!, C1_Vec(3) ! (i,k)
    Real( Kind = wp ) :: inv_RAx, inv_RAx3, diff(3), delta_OmAx
    Integer :: i, k

    inv_RAx = 1.0_wp / features(iatm,1)
    inv_RAx3 = inv_RAx * inv_RAx * inv_RAx

    diff(1) = xxx(alf(iatm,2)) - xxx(iatm)
    diff(2) = yyy(alf(iatm,2)) - yyy(iatm)
    diff(3) = zzz(alf(iatm,2)) - zzz(iatm)

    Call images_s(imcon, cell, diff(1), diff(2), diff(3))

    delta_OmAx = delta_Ax(iatm,Omega) 
    
    Do i = 1,3
      Do k = 1,3
        
        dC1k_da(i,k) = (inv_RAx*delta_OmAx*kronecker(i,k)) - (diff(i)*diff(k)*delta_OmAx*inv_RAx3)

      End Do
    End Do 

  End Function

  !=========================================  
  ! d(C_2k)/d(alpha_i)               [[C10]]
  !========================================= 

  Function dC2k_da(iatm, Omega)
 
    Use kinds_f90
    Use config_module, Only : imcon,cell,xxx,yyy,zzz
 
    Implicit None

    Integer,           Intent( In ) :: iatm, Omega

    Real( Kind = wp ) :: dC2k_da(3,3), dyj(3,3), invsqrt_ydoty, C12(3),diff(3),tmp(3)
    Integer :: k

    tmp = C_2(iatm)

    invsqrt_ydoty = 1.0_wp / dsqrt( y_vec(1)*y_vec(1)+y_vec(2)*y_vec(2)+y_vec(3)*y_vec(3)  )
  
    diff(1) = xxx(alf(iatm,2)) - xxx(iatm) 
    diff(2) = yyy(alf(iatm,2)) - yyy(iatm)
    diff(3) = zzz(alf(iatm,2)) - zzz(iatm) 

    Call images_s(imcon, cell, diff(1), diff(2), diff(3))

    C12(:) = EqC12(iatm,Omega,diff)
    dyj(:,:) = dyj_da(iatm, Omega,diff)

    Do k = 1,3
      dC2k_da(1,k) = (dyj(1,k)*invsqrt_ydoty) + (y_vec(k)*C12(1)) 
      dC2k_da(2,k) = (dyj(2,k)*invsqrt_ydoty) + (y_vec(k)*C12(2)) 
      dC2k_da(3,k) = (dyj(3,k)*invsqrt_ydoty) + (y_vec(k)*C12(3))
    End Do
    
  End Function  
 
  !=========================================  
  ! d(C_3k)/d(alpha_i)               [[C18]]
  !========================================= 

  Function dC3k_da(iatm, Omega, C1_Vec,C2_Vec,dC1k,dC2k)

    Use kinds_f90

    Implicit None

    Integer,           Intent( In ) :: iatm, Omega
    Real( Kind = wp ), Intent( In ) :: C1_Vec(3),c2_Vec(3),dC1k(3,3),dC2k(3,3)

    Real( Kind = wp ) :: dC3k_da(3,3)
    Integer :: i

    Do i = 1,3

      ! k = 1
      dC3k_da(i,1) =   ( (C1_Vec(2)*dC2k(i,3)) + (C2_Vec(3)*dC1k(i,2)) ) &
                   & - ( (C1_Vec(3)*dC2k(i,2)) + (C2_Vec(2)*dC1k(i,3)) )

      ! k = 2
      dC3k_da(i,2) =   ( (C1_Vec(3)*dC2k(i,1)) + (C2_Vec(1)*dC1k(i,3)) ) &
                   & - ( (C1_Vec(1)*dC2k(i,3)) + (C2_Vec(3)*dC1k(i,1)) )

      ! k = 3
      dC3k_da(i,3) =   ( (C1_Vec(1)*dC2k(i,2)) + (C2_Vec(2)*dC1k(i,1)) ) &
                   & - ( (C1_Vec(2)*dC2k(i,1)) + (C2_Vec(1)*dC1k(i,2)) )

    End Do
     
  End Function

  !=========================================  
  ! d(zeta_j)/d(alpha_i)              [[24]]
  !=========================================  

  Function dZj_da(iatm, jatm, Omega, j)

    Use kinds_f90
    Use config_module, Only : xxx,yyy,zzz,imcon,cell

    Implicit None

    Integer, Intent( In )           :: iatm, jatm, Omega, j

    Real( Kind = wp ) :: dZj_da(3), dCjk(3,3), Cjk(3), zeta(3), diff(3), C1_Vec(3), C2_Vec(3), C3_Vec(3),dC1k(3,3),dC2k(3,3)
    Integer :: i,k

    diff(1) = xxx(jatm) - xxx(iatm)
    diff(2) = yyy(jatm) - yyy(iatm)
    diff(3) = zzz(jatm) - zzz(iatm)

    Call images_s(imcon, cell, diff(1), diff(2), diff(3))

    !C1_Vec = C_1(iatm)
    C2_Vec = C_2(iatm)
    !C3_Vec = C_3(iatm,C1_Vec,C2_Vec)

    C1_Vec = lxx(iatm,:)
!    C2_Vec = lyy(iatm,:)
    C3_Vec = lzz(iatm,:)

    Select Case (j)

      Case(1)

        Cjk = C1_Vec
        dCjk(:,:) = dC1k_da(iatm,Omega,C1_Vec)

      Case(2)

        Cjk = C2_Vec
        dCjk(:,:) = dC2k_da(iatm,Omega)

      Case(3)

        Cjk = C3_Vec
        dC1k(:,:) = dC1k_da(iatm,Omega,C1_Vec)
        dC2k(:,:) = dC2k_da(iatm,Omega)
        dCjk(:,:) = dC3k_da(iatm,Omega,C1_Vec,C2_Vec,dC1k,dC2k)

    End Select

    dZj_da(:) = 0.0_wp

    Do i = 1,3
      Do k = 1,3
    
        dZj_da(i) = dZj_da(i) + dCjk(i,k)*diff(k) + Cjk(k)*kronecker(i,k)*delta_An(iatm,jatm,Omega)

      End Do
    End Do

    zeta(1) = (C1_Vec(1)*diff(1)) + (C1_Vec(2)*diff(2)) + (C1_Vec(3)*diff(3))
    zeta(2) = (C2_Vec(1)*diff(1)) + (C2_Vec(2)*diff(2)) + (C2_Vec(3)*diff(3))
    zeta(3) = (C3_Vec(1)*diff(1)) + (C3_Vec(2)*diff(2)) + (C3_Vec(3)*diff(3))

  End Function



  !=========================================  
  ! d(chi)/d(alpha_i) for all i       [[20]] using: [[B12]] &  [[B13]]
  !=========================================  

  Function dChi_da( iatm, Omega )
    Use kinds_f90
    Use config_module, Only : xxx,yyy,zzz,imcon,cell

    Implicit None

    Integer,           Intent( In ) :: iatm, Omega

    Real( Kind = wp ) :: dChi_da(3), &
                       & invsin_Chi, &
                       & RxdotRxy, R_Ax, R_Axy, &
                       & xdiff1, ydiff1, zdiff1, &
                       & xdiff2, ydiff2, zdiff2, temp(3)

    R_Ax  = features(iatm,1)
    R_Axy = features(iatm,2)

    xdiff1 = xxx(alf(iatm,2)) - xxx(iatm)
    ydiff1 = yyy(alf(iatm,2)) - yyy(iatm)
    zdiff1 = zzz(alf(iatm,2)) - zzz(iatm)

    Call images_s(imcon, cell, xdiff1, ydiff1, zdiff1)

    xdiff2 = xxx(alf(iatm,3)) - xxx(iatm)
    ydiff2 = yyy(alf(iatm,3)) - yyy(iatm)
    zdiff2 = zzz(alf(iatm,3)) - zzz(iatm)

    Call images_s(imcon, cell, xdiff2, ydiff2, zdiff2)

    RxdotRxy = (xdiff1*xdiff2) + (ydiff1*ydiff2) + (zdiff1*zdiff2)

    invsin_Chi = (-1.0_wp) / dsin( features(iatm,3) )

    dChi_da(:) = ( ((-1.0_wp) * dR_da(iatm,alf(iatm,2),1,Omega) * RxdotRxy) / (R_Ax * R_Ax * R_Axy)) &
            & + ( dRaxdotRaxy_da(iatm,Omega,xdiff1,ydiff1,zdiff1,xdiff2,ydiff2,zdiff2) / (R_Ax * R_Axy)) &
            & + ( ((-1.0_wp) * dR_da(iatm,alf(iatm,3),2,Omega) * RxdotRxy) / (R_Ax * R_Axy * R_Axy ) )                                                                                                                                  
    dChi_da(:) = dChi_da(:) * invsin_Chi

    temp(:) =  dRaxdotRaxy_da(iatm,Omega,xdiff1,ydiff1,zdiff1,xdiff2,ydiff2,zdiff2)

  End Function 

  !=========================================  
  ! d(Theta_{An})/d(alpha_i)          [[21]]
  !=========================================  

  Function dTheta_da(iatm, jatm, cFeat, zeta3, Omega)

    Use kinds_f90

    Integer,           Intent( In ) :: iatm, jatm, cFeat, Omega
    Real( Kind = wp ), Intent( In ) :: zeta3

    Real( Kind = wp ) :: dTheta_da(3), dZ3(3), invSinTheta, R_An, inv_RAn, dInvR(3)
    Integer :: cFeat_R

    cFeat_R = cFeat-1

    R_An = features(iatm,cFeat_R)

    inv_RAn = 1.0_wp / R_An

    invSinTheta = (-1.0_wp) / dsin(features(iatm,cFeat)) 

    dInvR(:) = (-1.0_wp) * inv_RAn * inv_RAn * dR_da(iatm,jatm,cFeat_R,Omega) ! [[B14]]

    dZ3(:) = dZj_da(iatm,jatm,Omega,3)

    dTheta_da(1) = invSinTheta * ( (zeta3*dInvR(1)) + (inv_RAn*dZ3(1)) )
    dTheta_da(2) = invSinTheta * ( (zeta3*dInvR(2)) + (inv_RAn*dZ3(2)) )
    dTheta_da(3) = invSinTheta * ( (zeta3*dInvR(3)) + (inv_RAn*dZ3(3)) )

  End Function

  !=========================================  
  ! d(Phi_{An})/d(alpha_i)            [[22]]
  !=========================================  

  Function dPhi_da(iatm, jatm, zeta1, zeta2, cFeat, Omega)

    Use kinds_f90

    Implicit None

    Integer,           Intent( In ) :: iatm, jatm, cFeat, Omega
    Real( Kind = wp ), Intent( In ) :: zeta1, zeta2
    
    Real( Kind = wp ) :: cos2phi, dZ1(3), dZ2(3), dPhi_da(3)

    cos2phi = dcos(features(iatm,cFeat))

    cos2phi = cos2phi * cos2phi

    dZ1(:) = dZj_da(iatm,jatm,Omega,1)
    dZ2(:) = dZj_da(iatm,jatm,Omega,2)

    
    dPhi_da(1) = cos2phi * ( ( ((-1.0_wp) * zeta2 * dZ1(1)) / (zeta1*zeta1)) + (dZ2(1)/zeta1) )
    dPhi_da(2) = cos2phi * ( ( ((-1.0_wp) * zeta2 * dZ1(2)) / (zeta1*zeta1)) + (dZ2(2)/zeta1) )
    dPhi_da(3) = cos2phi * ( ( ((-1.0_wp) * zeta2 * dZ1(3)) / (zeta1*zeta1)) + (dZ2(3)/zeta1) )
    
  End Function
  
  !=========================================  
  ! Form C Matrix 
  !========================================= 
!ER  this function serves the purpose of allocation of memory for C matrix. 
  Function Form_C_Matrix(iatm)

    Use kinds_f90 

    Implicit None

    Integer, Intent( In ) :: iatm
    Real ( Kind = wp ) :: Form_C_Matrix(3,3)

    Form_C_Matrix(1,:) = lxx(iatm,:)
    Form_C_Matrix(2,:) = lyy(iatm,:) 
    Form_C_Matrix(3,:) = lzz(iatm,:) 

  End Function

  !=========================================  
  ! Form dC Matrix (the derivative of C)
  !=========================================  
!ER  this function serves the purpose of allocation of memory for dC matrix.
  Function Form_dC_Matrix(iatm,Omega)

    Use kinds_f90 

    Implicit None

    Integer, Intent( In ) :: iatm, Omega
    Real ( Kind = wp ) :: Form_dC_Matrix(3,3,3)

    Form_dC_Matrix(1,:,:) = dCxx(iatm,Omega,:,:)
    Form_dC_Matrix(2,:,:) = dCyy(iatm,Omega,:,:) 
    Form_dC_Matrix(3,:,:) = dCzz(iatm,Omega,:,:) 

  End Function
  
  !=========================================  
  ! [[B15]] - my rederived form [[B15]] equiv.
  !=========================================  
  !Altered by BCBS 03/2019 to account for PBC's
  Function dRdotR_da( iatm, jatm, Omega, xxx, yyy, zzz ) ! d(Rax . Rax) / d(alpha)

    Use kinds_f90
    Use config_module, Only : imcon,cell

    Implicit None

    Integer, Intent( In ) :: iatm, jatm, Omega
    Real( kind = wp ), Intent( In ) :: xxx(*), yyy(*), zzz(*)
    Real( Kind = wp ) :: dRdotR_da(3)

    If ( Omega .EQ. iatm ) Then

      dRdotR_da(1) = ( xxx(iatm) - xxx(jatm) )
      dRdotR_da(2) = ( yyy(iatm) - yyy(jatm) )
      dRdotR_da(3) = ( zzz(iatm) - zzz(jatm) )

      Call images_s(imcon,cell,dRdotR_da(1),dRdotR_da(2),dRdotR_da(3))

      dRdotR_da(1) = 2.0_wp * dRdotR_da(1)
      dRdotR_da(2) = 2.0_wp * dRdotR_da(2)
      dRdotR_da(3) = 2.0_wp * dRdotR_da(3)

      !dRdotR_da(1) = 2.0_wp * ( xxx(iatm) - xxx(jatm) )
      !dRdotR_da(2) = 2.0_wp * ( yyy(iatm) - yyy(jatm) )
      !dRdotR_da(3) = 2.0_wp * ( zzz(iatm) - zzz(jatm) )

    Else If ( Omega .EQ. jatm ) Then

      dRdotR_da(1) = ( xxx(jatm) - xxx(iatm) )
      dRdotR_da(2) = ( yyy(jatm) - yyy(iatm) )
      dRdotR_da(3) = ( zzz(jatm) - zzz(iatm) )

      Call images_s(imcon,cell,dRdotR_da(1),dRdotR_da(2),dRdotR_da(3))
                                                                    
      dRdotR_da(1) = 2.0_wp * dRdotR_da(1)
      dRdotR_da(2) = 2.0_wp * dRdotR_da(2)
      dRdotR_da(3) = 2.0_wp * dRdotR_da(3)

      !dRdotR_da(1) = 2.0_wp * ( xxx(jatm) - xxx(iatm) )
      !dRdotR_da(2) = 2.0_wp * ( yyy(jatm) - yyy(iatm) )
      !dRdotR_da(3) = 2.0_wp * ( zzz(jatm) - zzz(iatm) )

    Else

      dRdotR_da(1) = 0.0_wp
      dRdotR_da(2) = 0.0_wp
      dRdotR_da(3) = 0.0_wp

    End If 

  End Function
  
  !=========================================  
  ! [[39a]] & [[39b]]
  !========================================= 
! Rewritten as 3 explicit functions for optimisation purposes - BCBS.
  Function dwjR_da( iatm, xdiff, ydiff, zdiff, dRi_da, Omega, j, xxx, yyy, zzz )

    Use kinds_f90
    Use config_module, Only : ltg

    Implicit None

    Integer, Intent ( In ) :: iatm, Omega, j
    Real ( Kind = wp ), Intent( In ) :: dRi_da, xdiff, ydiff, zdiff, xxx(*), yyy(*), zzz(*)
    Real ( Kind = wp ) :: dwjR_da(3), dCjk(3,3),C(3)!, C(3,3)
    Integer            :: i,alias
    
    alias = atm_label_alias(ltg(Omega))

    If (j .EQ. 1) Then
      dCjk = dCxx(iatm,alias,:,:)
      C = lxx(iatm,:)
    Else If (j .EQ. 2) Then
      dCjk = dCyy(iatm,alias,:,:)
      C = lyy(iatm,:)
    Else
      dCjk = dCzz(iatm,alias,:,:)
      C = lzz(iatm,:)
    End If

    Do i=1,3
      dwjR_da(i) = (C(i) * dRi_da) + (dCjk(i,1) * xdiff) + (dCjk(i,2) * ydiff) + (dCjk(i,3) * zdiff)
    End Do

  End Function

  Function dwjR_da_1( iatm, xdiff, ydiff, zdiff, dRi_da, Omega, xxx, yyy, zzz )

    Use kinds_f90
    Use config_module, Only : ltg

    Implicit None

    Integer, Intent ( In ) :: iatm, Omega
    Real ( Kind = wp ), Intent( In ) :: dRi_da, xdiff, ydiff, zdiff, xxx(*), yyy(*), zzz(*)
    Real ( Kind = wp ) :: dwjR_da_1(3), dCjk(3,3),C(3)!, C(3,3)
    Integer            :: i,alias

    !alias = atm_label_alias(ltg(Omega))
    If (Omega .EQ. iatm) Then
      alias = 1
    Else If (Omega .EQ. alf(iatm,2)) Then
      alias = 2
    Else If (Omega .EQ. alf(iatm,3)) Then
      alias = 3
    Else
      alias = 100
    End If
    If (alias <= 3) Then
      dwjR_da_1(:) = (lxx(iatm,:) * dRi_da) + (dCxx(iatm,alias,:,1) * xdiff) + (dCxx(iatm,alias,:,2) * ydiff) + (dCxx(iatm,alias,:,3) * zdiff)
    Else
      dwjR_da_1(:) = (lxx(iatm,:) * dRi_da)
    End If

  End Function

  Function dwjR_da_2( iatm, xdiff, ydiff, zdiff, dRi_da, Omega, xxx, yyy, zzz )

    Use kinds_f90
    Use config_module, Only : ltg

    Implicit None

    Integer, Intent ( In ) :: iatm, Omega
    Real ( Kind = wp ), Intent( In ) :: dRi_da, xdiff, ydiff, zdiff, xxx(*), yyy(*), zzz(*)
    Real ( Kind = wp ) :: dwjR_da_2(3), dCjk(3,3),C(3)!, C(3,3)
    Integer            :: i,alias

    !alias = atm_label_alias(ltg(Omega))
    If (Omega .EQ. iatm) Then
      alias = 1
    Else If (Omega .EQ. alf(iatm,2)) Then
      alias = 2
    Else If (Omega .EQ. alf(iatm,3)) Then
      alias = 3
    Else
      alias = 100
    End If

    If (alias <= 3) Then
      Do i=1,3
        dwjR_da_2(i) = (lyy(iatm,i) * dRi_da) + (dCyy(iatm,alias,i,1) * xdiff) + (dCyy(iatm,alias,i,2) * ydiff) + (dCyy(iatm,alias,i,3) * zdiff)
      End Do
    Else
      Do i=1,3
        dwjR_da_2(i) = (lyy(iatm,i) * dRi_da)
      End Do
    End If

  End Function

  Function dwjR_da_3( iatm, xdiff, ydiff, zdiff, dRi_da, Omega, xxx, yyy, zzz )

    Use kinds_f90
    Use config_module, Only : ltg

    Implicit None

    Integer, Intent ( In ) :: iatm, Omega
    Real ( Kind = wp ), Intent( In ) :: dRi_da, xdiff, ydiff, zdiff, xxx(*), yyy(*), zzz(*)
    Real ( Kind = wp ) :: dwjR_da_3(3), dCjk(3,3),C(3)!, C(3,3)
    Integer            :: i,alias

    !alias = atm_label_alias(ltg(Omega))
    If (Omega .EQ. iatm) Then
      alias = 1
    Else If (Omega .EQ. alf(iatm,2)) Then
      alias = 2
    Else If (Omega .EQ. alf(iatm,3)) Then
      alias = 3
    Else
      alias = 100
    End If
    If (alias <= 3) Then
      Do i=1,3
        dwjR_da_3(i) = (lzz(iatm,i) * dRi_da) + (dCzz(iatm,alias,i,1) * xdiff) + (dCzz(iatm,alias,i,2) * ydiff) + (dCzz(iatm,alias,i,3) * zdiff)
      End Do
    Else
      Do i=1,3
        dwjR_da_3(i) = (lzz(iatm,i) * dRi_da)
      End Do
    End If

  End Function
  
  !=========================================  
  ! [[43]]
  !========================================= 

  Function dwAjwBk_da( iatm, jatm, j, k, Omega, xxx, yyy, zzz )

    Use kinds_f90
    Use config_module, Only : ltg

    Implicit None

    Integer, Intent( In ) :: iatm, jatm, j, k, Omega 
    Real( Kind = wp ), Intent( In ) :: xxx(*), yyy(*), zzz(*) 
    Real( Kind = wp ) :: dCAjk(3,3), dCBjk(3,3), dwAjwBk_da(3), CA(3), CB(3)
    Integer           :: i,alias
 
 !   alias = atm_label_alias(ltg(Omega))

! Made explicit by BCBS for optimisation purposes.
    If ((Omega .EQ. iatm) .OR. (Omega .EQ. alf(iatm,2)) .OR. (Omega .EQ. alf(iatm,3))) Then

    If (Omega .EQ. iatm) Then
      alias = 1
    Else If (Omega .EQ. alf(iatm,2)) Then
      alias = 2
    Else If (Omega .EQ. alf(iatm,3)) Then
      alias = 3
    End If

    If ((j .EQ. 1) .AND. (k .EQ. 1)) Then

      Do i=1,1!3
        dwAjwBk_da(:) =   ( lxx(jatm,1) * dCxx(iatm,alias,:,1) ) &
                      !& + ( lxx(iatm,1) * dCxx(jatm,alias,:,1) ) &
                      & + ( lxx(jatm,2) * dCxx(iatm,alias,:,2) ) &
                      !& + ( lxx(iatm,2) * dCxx(jatm,alias,:,2) ) &
                      & + ( lxx(jatm,3) * dCxx(iatm,alias,:,3) )
                      !& + ( lxx(iatm,3) * dCxx(jatm,alias,:,3) )
      End Do

    Else If ((j .EQ. 1) .AND. (k .EQ. 2)) Then

      Do i=1,1!3
        dwAjwBk_da(:) =   ( lyy(jatm,1) * dCxx(iatm,alias,:,1) ) &
                      !& + ( lxx(iatm,1) * dCyy(jatm,alias,:,1) ) &
                      & + ( lyy(jatm,2) * dCxx(iatm,alias,:,2) ) &
                      !& + ( lxx(iatm,2) * dCyy(jatm,alias,:,2) ) &
                      & + ( lyy(jatm,3) * dCxx(iatm,alias,:,3) )
                      !& + ( lxx(iatm,3) * dCyy(jatm,alias,:,3) )
      End Do

    Else If ((j .EQ. 1) .AND. (k .EQ. 3)) Then

      Do i=1,1!3
        dwAjwBk_da(:) =   ( lzz(jatm,1) * dCxx(iatm,alias,:,1) ) &
                      !& + ( lxx(iatm,1) * dCzz(jatm,alias,:,1) ) &
                      & + ( lzz(jatm,2) * dCxx(iatm,alias,:,2) ) &
                      !& + ( lxx(iatm,2) * dCzz(jatm,alias,:,2) ) &
                      & + ( lzz(jatm,3) * dCxx(iatm,alias,:,3) )
                      !& + ( lxx(iatm,3) * dCzz(jatm,alias,:,3) )
      End Do

    Else If ((j .EQ. 2) .AND. (k .EQ. 1)) Then

      Do i=1,1!3
        dwAjwBk_da(:) =   ( lxx(jatm,1) * dCyy(iatm,alias,:,1) ) &
                      !& + ( lyy(iatm,1) * dCxx(jatm,alias,:,1) ) &
                      & + ( lxx(jatm,2) * dCyy(iatm,alias,:,2) ) &
                      !& + ( lyy(iatm,2) * dCxx(jatm,alias,:,2) ) &
                      & + ( lxx(jatm,3) * dCyy(iatm,alias,:,3) )
                      !& + ( lyy(iatm,3) * dCxx(jatm,alias,:,3) )
      End Do

    Else If ((j .EQ. 2) .AND. (k .EQ. 2)) Then

      Do i=1,1!3
        dwAjwBk_da(:) =   ( lyy(jatm,1) * dCyy(iatm,alias,:,1) ) &
                      !& + ( lyy(iatm,1) * dCyy(jatm,alias,:,1) ) &
                      & + ( lyy(jatm,2) * dCyy(iatm,alias,:,2) ) &
                      !& + ( lyy(iatm,2) * dCyy(jatm,alias,:,2) ) &
                      & + ( lyy(jatm,3) * dCyy(iatm,alias,:,3) )
                      !& + ( lyy(iatm,3) * dCyy(jatm,alias,:,3) )
      End Do

    Else If ((j .EQ. 2) .AND. (k .EQ. 3)) Then

      Do i=1,1!3
        dwAjwBk_da(:) =   ( lzz(jatm,1) * dCyy(iatm,alias,:,1) ) &
                      !& + ( lyy(iatm,1) * dCzz(jatm,alias,:,1) ) &
                      & + ( lzz(jatm,2) * dCyy(iatm,alias,:,2) ) &
                      !& + ( lyy(iatm,2) * dCzz(jatm,alias,:,2) ) &
                      & + ( lzz(jatm,3) * dCyy(iatm,alias,:,3) )
                      !& + ( lyy(iatm,3) * dCzz(jatm,alias,:,3) )
      End Do

    Else If ((j .EQ. 3) .AND. (k .EQ. 1)) Then

      Do i=1,1!3
        dwAjwBk_da(:) =   ( lxx(jatm,1) * dCzz(iatm,alias,:,1) ) &
                      !& + ( lzz(iatm,1) * dCxx(jatm,alias,:,1) ) &
                      & + ( lxx(jatm,2) * dCzz(iatm,alias,:,2) ) &
                      !& + ( lzz(iatm,2) * dCxx(jatm,alias,:,2) ) &
                      & + ( lxx(jatm,3) * dCzz(iatm,alias,:,3) )
                      !& + ( lzz(iatm,3) * dCxx(jatm,alias,:,3) )
      End Do

    Else If ((j .EQ. 3) .AND. (k .EQ. 2)) Then

      Do i=1,1!3
        dwAjwBk_da(:) =   ( lyy(jatm,1) * dCzz(iatm,alias,:,1) ) &
                      !& + ( lzz(iatm,1) * dCyy(jatm,alias,:,1) ) &
                      & + ( lyy(jatm,2) * dCzz(iatm,alias,:,2) ) &
                      !& + ( lzz(iatm,2) * dCyy(jatm,alias,:,2) ) &
                      & + ( lyy(jatm,3) * dCzz(iatm,alias,:,3) )
                      !& + ( lzz(iatm,3) * dCyy(jatm,alias,:,3) )
      End Do

    Else If ((j .EQ. 3) .AND. (k .EQ. 3)) Then

      Do i=1,1!3
        dwAjwBk_da(:) =   ( lzz(jatm,1) * dCzz(iatm,alias,:,1) ) &
                      !& + ( lzz(iatm,1) * dCzz(jatm,alias,:,1) ) &
                      & + ( lzz(jatm,2) * dCzz(iatm,alias,:,2) ) &
                      !& + ( lzz(iatm,2) * dCzz(jatm,alias,:,2) ) &
                      & + ( lzz(jatm,3) * dCzz(iatm,alias,:,3) )
                      !& + ( lzz(iatm,3) * dCzz(jatm,alias,:,3) )
      End Do

    End If

    Else If ((Omega .EQ. jatm) .OR. (Omega .EQ. alf(jatm,2)) .OR. (Omega .EQ. alf(jatm,3))) Then

    If (Omega .EQ. jatm) Then
      alias = 1
    Else If (Omega .EQ. alf(jatm,2)) Then
      alias = 2
    Else If (Omega .EQ. alf(jatm,3)) Then
      alias = 3
    End If

    If ((j .EQ. 1) .AND. (k .EQ. 1)) Then
                                                                   
      Do i=1,1!3
        !dwAjwBk_da(:) =   ( lxx(jatm,1) * dCxx(iatm,alias,:,1) ) &
        dwAjwBk_da(:) = ( lxx(iatm,1) * dCxx(jatm,alias,:,1) ) &
                      !& + ( lxx(jatm,2) * dCxx(iatm,alias,:,2) ) &
                      & + ( lxx(iatm,2) * dCxx(jatm,alias,:,2) ) &
                      !& + ( lxx(jatm,3) * dCxx(iatm,alias,:,3) ) &
                      & + ( lxx(iatm,3) * dCxx(jatm,alias,:,3) )
      End Do
                                                                   
    Else If ((j .EQ. 1) .AND. (k .EQ. 2)) Then
                                                                   
      Do i=1,1!3
        !dwAjwBk_da(:) =   ( lyy(jatm,1) * dCxx(iatm,alias,:,1) ) &
        dwAjwBk_da(:) = ( lxx(iatm,1) * dCyy(jatm,alias,:,1) ) &
                      !& + ( lyy(jatm,2) * dCxx(iatm,alias,:,2) ) &
                      & + ( lxx(iatm,2) * dCyy(jatm,alias,:,2) ) &
                      !& + ( lyy(jatm,3) * dCxx(iatm,alias,:,3) ) &
                      & + ( lxx(iatm,3) * dCyy(jatm,alias,:,3) )
      End Do
                                                                   
    Else If ((j .EQ. 1) .AND. (k .EQ. 3)) Then
                                                                   
      Do i=1,1!3
        !dwAjwBk_da(:) =   ( lzz(jatm,1) * dCxx(iatm,alias,:,1) ) &
        dwAjwBk_da(:) = ( lxx(iatm,1) * dCzz(jatm,alias,:,1) ) &
                      !& + ( lzz(jatm,2) * dCxx(iatm,alias,:,2) ) &
                      & + ( lxx(iatm,2) * dCzz(jatm,alias,:,2) ) &
                      !& + ( lzz(jatm,3) * dCxx(iatm,alias,:,3) ) &
                      & + ( lxx(iatm,3) * dCzz(jatm,alias,:,3) )
      End Do
                                                                   
    Else If ((j .EQ. 2) .AND. (k .EQ. 1)) Then
                                                                   
      Do i=1,1!3
        !dwAjwBk_da(:) =   ( lxx(jatm,1) * dCyy(iatm,alias,:,1) ) &
        dwAjwBk_da(:) =  ( lyy(iatm,1) * dCxx(jatm,alias,:,1) ) &
                      !& + ( lxx(jatm,2) * dCyy(iatm,alias,:,2) ) &
                      & + ( lyy(iatm,2) * dCxx(jatm,alias,:,2) ) &
                      !& + ( lxx(jatm,3) * dCyy(iatm,alias,:,3) ) &
                      & + ( lyy(iatm,3) * dCxx(jatm,alias,:,3) )
      End Do
                                                                   
    Else If ((j .EQ. 2) .AND. (k .EQ. 2)) Then
                                                                   
      Do i=1,1!3
        !dwAjwBk_da(:) =   ( lyy(jatm,1) * dCyy(iatm,alias,:,1) ) &
        dwAjwBk_da(:) = ( lyy(iatm,1) * dCyy(jatm,alias,:,1) ) &
                      !& + ( lyy(jatm,2) * dCyy(iatm,alias,:,2) ) &
                      & + ( lyy(iatm,2) * dCyy(jatm,alias,:,2) ) &
                      !& + ( lyy(jatm,3) * dCyy(iatm,alias,:,3) ) &
                      & + ( lyy(iatm,3) * dCyy(jatm,alias,:,3) )
      End Do
                                                                   
    Else If ((j .EQ. 2) .AND. (k .EQ. 3)) Then
                                                                   
      Do i=1,1!3
        !dwAjwBk_da(:) =   ( lzz(jatm,1) * dCyy(iatm,alias,:,1) ) &
        dwAjwBk_da(:) = ( lyy(iatm,1) * dCzz(jatm,alias,:,1) ) &
                      !& + ( lzz(jatm,2) * dCyy(iatm,alias,:,2) ) &
                      & + ( lyy(iatm,2) * dCzz(jatm,alias,:,2) ) &
                      !& + ( lzz(jatm,3) * dCyy(iatm,alias,:,3) ) &
                      & + ( lyy(iatm,3) * dCzz(jatm,alias,:,3) )
      End Do
                                                                   
    Else If ((j .EQ. 3) .AND. (k .EQ. 1)) Then
                                                                   
      Do i=1,1!3
        !dwAjwBk_da(:) =   ( lxx(jatm,1) * dCzz(iatm,alias,:,1) ) &
        dwAjwBk_da(:) = ( lzz(iatm,1) * dCxx(jatm,alias,:,1) ) &
                      !& + ( lxx(jatm,2) * dCzz(iatm,alias,:,2) ) &
                      & + ( lzz(iatm,2) * dCxx(jatm,alias,:,2) ) &
                      !& + ( lxx(jatm,3) * dCzz(iatm,alias,:,3) ) &
                      & + ( lzz(iatm,3) * dCxx(jatm,alias,:,3) )
      End Do
                                                                   
    Else If ((j .EQ. 3) .AND. (k .EQ. 2)) Then
                                                                   
      Do i=1,1!3
        !dwAjwBk_da(:) =   ( lyy(jatm,1) * dCzz(iatm,alias,:,1) ) &
        dwAjwBk_da(:) = ( lzz(iatm,1) * dCyy(jatm,alias,:,1) ) &
                      !& + ( lyy(jatm,2) * dCzz(iatm,alias,:,2) ) &
                      & + ( lzz(iatm,2) * dCyy(jatm,alias,:,2) ) &
                      !& + ( lyy(jatm,3) * dCzz(iatm,alias,:,3) ) &
                      & + ( lzz(iatm,3) * dCyy(jatm,alias,:,3) )
      End Do
                                                                   
    Else If ((j .EQ. 3) .AND. (k .EQ. 3)) Then
                                                                   
      Do i=1,1!3
        !dwAjwBk_da(:) =   ( lzz(jatm,1) * dCzz(iatm,alias,:,1) ) &
        dwAjwBk_da(:) = ( lzz(iatm,1) * dCzz(jatm,alias,:,1) ) &
                      !& + ( lzz(jatm,2) * dCzz(iatm,alias,:,2) ) &
                      & + ( lzz(iatm,2) * dCzz(jatm,alias,:,2) ) &
                      !& + ( lzz(jatm,3) * dCzz(iatm,alias,:,3) ) &
                      & + ( lzz(iatm,3) * dCzz(jatm,alias,:,3) )
      End Do
                                                                   
    End If
    End If


  End Function

!Previous function has been split into 9 functions for optimisation purposes
!BCBS.
  Function dwAjwBk_xx_da( iatm, jatm, Omega, xxx, yyy, zzz )

    Use kinds_f90

    Implicit None

    Integer, Intent( In ) :: iatm, jatm, Omega 
    Real( Kind = wp ), Intent( In ) :: xxx(*), yyy(*), zzz(*) 
    Real( Kind = wp ) :: dwAjwBk_xx_da(3)

      If (Omega .EQ. iatm) Then
        dwAjwBk_xx_da(:) =   ( lxx(jatm,1) * dCxx(iatm,1,:,1) ) &
                      & + ( lxx(jatm,2) * dCxx(iatm,1,:,2) ) &
                      & + ( lxx(jatm,3) * dCxx(iatm,1,:,3) )
      Else If (Omega .EQ. alf(iatm,2)) Then
        dwAjwBk_xx_da(:) =   ( lxx(jatm,1) * dCxx(iatm,2,:,1) ) &
                      & + ( lxx(jatm,2) * dCxx(iatm,2,:,2) ) &
                      & + ( lxx(jatm,3) * dCxx(iatm,2,:,3) )
      Else If (Omega .EQ. alf(iatm,3)) Then
        dwAjwBk_xx_da(:) =   ( lxx(jatm,1) * dCxx(iatm,3,:,1) ) &
                      & + ( lxx(jatm,2) * dCxx(iatm,3,:,2) ) &
                      & + ( lxx(jatm,3) * dCxx(iatm,3,:,3) )
       Else If (Omega .EQ. jatm) Then
         dwAjwBk_xx_da(:) = ( lxx(iatm,1) * dCxx(jatm,1,:,1) ) &
                       & + ( lxx(iatm,2) * dCxx(jatm,1,:,2) ) &
                       & + ( lxx(iatm,3) * dCxx(jatm,1,:,3) )
       Else If (Omega .EQ. alf(jatm,2)) Then
         dwAjwBk_xx_da(:) = ( lxx(iatm,1) * dCxx(jatm,2,:,1) ) &
                       & + ( lxx(iatm,2) * dCxx(jatm,2,:,2) ) &
                       & + ( lxx(iatm,3) * dCxx(jatm,2,:,3) )
       Else If (Omega .EQ. alf(jatm,3)) Then
         dwAjwBk_xx_da(:) = ( lxx(iatm,1) * dCxx(jatm,3,:,1) ) &
                       & + ( lxx(iatm,2) * dCxx(jatm,3,:,2) ) &
                       & + ( lxx(iatm,3) * dCxx(jatm,3,:,3) )
       End If

  End Function

  Function dwAjwBk_yx_da( iatm, jatm, Omega, xxx, yyy, zzz )

    Use kinds_f90

    Implicit None

    Integer, Intent( In ) :: iatm, jatm, Omega 
    Real( Kind = wp ), Intent( In ) :: xxx(*), yyy(*), zzz(*) 
    Real( Kind = wp ) :: dwAjwBk_yx_da(3)

      If (Omega .EQ. iatm) Then
        dwAjwBk_yx_da(:) =   ( lxx(jatm,1) * dCyy(iatm,1,:,1) ) &
                      & + ( lxx(jatm,2) * dCyy(iatm,1,:,2) ) &
                      & + ( lxx(jatm,3) * dCyy(iatm,1,:,3) )
      Else If (Omega .EQ. alf(iatm,2)) Then
        dwAjwBk_yx_da(:) =   ( lxx(jatm,1) * dCyy(iatm,2,:,1) ) &
                      & + ( lxx(jatm,2) * dCyy(iatm,2,:,2) ) &
                      & + ( lxx(jatm,3) * dCyy(iatm,2,:,3) )
      Else If (Omega .EQ. alf(iatm,3)) Then
        dwAjwBk_yx_da(:) =   ( lxx(jatm,1) * dCyy(iatm,3,:,1) ) &
                      & + ( lxx(jatm,2) * dCyy(iatm,3,:,2) ) &
                      & + ( lxx(jatm,3) * dCyy(iatm,3,:,3) )
       Else If (Omega .EQ. jatm) Then
         dwAjwBk_yx_da(:) =  ( lyy(iatm,1) * dCxx(jatm,1,:,1) ) &   
                       & + ( lyy(iatm,2) * dCxx(jatm,1,:,2) ) &
                       & + ( lyy(iatm,3) * dCxx(jatm,1,:,3) )
       Else If (Omega .EQ. alf(jatm,2)) Then
         dwAjwBk_yx_da(:) =  ( lyy(iatm,1) * dCxx(jatm,2,:,1) ) &   
                       & + ( lyy(iatm,2) * dCxx(jatm,2,:,2) ) &
                       & + ( lyy(iatm,3) * dCxx(jatm,2,:,3) )
       Else If (Omega .EQ. alf(jatm,3)) Then
         dwAjwBk_yx_da(:) =  ( lyy(iatm,1) * dCxx(jatm,3,:,1) ) &   
                       & + ( lyy(iatm,2) * dCxx(jatm,3,:,2) ) &
                       & + ( lyy(iatm,3) * dCxx(jatm,3,:,3) )
       End If

  End Function

  Function dwAjwBk_zx_da( iatm, jatm, Omega, xxx, yyy, zzz )

    Use kinds_f90

    Implicit None

    Integer, Intent( In ) :: iatm, jatm, Omega 
    Real( Kind = wp ), Intent( In ) :: xxx(*), yyy(*), zzz(*) 
    Real( Kind = wp ) :: dwAjwBk_zx_da(3)

      If (Omega .EQ. iatm) Then
        dwAjwBk_zx_da(:) =   ( lxx(jatm,1) * dCzz(iatm,1,:,1) ) &
                      & + ( lxx(jatm,2) * dCzz(iatm,1,:,2) ) &
                      & + ( lxx(jatm,3) * dCzz(iatm,1,:,3) )
      Else If (Omega .EQ. alf(iatm,2)) Then
        dwAjwBk_zx_da(:) =   ( lxx(jatm,1) * dCzz(iatm,2,:,1) ) &
                      & + ( lxx(jatm,2) * dCzz(iatm,2,:,2) ) &
                      & + ( lxx(jatm,3) * dCzz(iatm,2,:,3) )
      Else If (Omega .EQ. alf(iatm,3)) Then
        dwAjwBk_zx_da(:) =   ( lxx(jatm,1) * dCzz(iatm,3,:,1) ) &
                      & + ( lxx(jatm,2) * dCzz(iatm,3,:,2) ) &
                      & + ( lxx(jatm,3) * dCzz(iatm,3,:,3) )
       Else If (Omega .EQ. jatm) Then
         dwAjwBk_zx_da(:) = ( lzz(iatm,1) * dCxx(jatm,1,:,1) ) &
                       & + ( lzz(iatm,2) * dCxx(jatm,1,:,2) ) &
                       & + ( lzz(iatm,3) * dCxx(jatm,1,:,3) )
       Else If (Omega .EQ. alf(jatm,2)) Then
         dwAjwBk_zx_da(:) = ( lzz(iatm,1) * dCxx(jatm,2,:,1) ) &
                       & + ( lzz(iatm,2) * dCxx(jatm,2,:,2) ) &
                       & + ( lzz(iatm,3) * dCxx(jatm,2,:,3) )
       Else If (Omega .EQ. alf(jatm,3)) Then
         dwAjwBk_zx_da(:) = ( lzz(iatm,1) * dCxx(jatm,3,:,1) ) &
                       & + ( lzz(iatm,2) * dCxx(jatm,3,:,2) ) &
                       & + ( lzz(iatm,3) * dCxx(jatm,3,:,3) )
       End If

  End Function

  Function dwAjwBk_xy_da( iatm, jatm, Omega, xxx, yyy, zzz )

    Use kinds_f90

    Implicit None

    Integer, Intent( In ) :: iatm, jatm, Omega 
    Real( Kind = wp ), Intent( In ) :: xxx(*), yyy(*), zzz(*) 
    Real( Kind = wp ) :: dwAjwBk_xy_da(3)

      If (Omega .EQ. iatm) Then
        dwAjwBk_xy_da(:) = ( lyy(jatm,1) * dCxx(iatm,1,:,1) ) &
                      & + ( lyy(jatm,2) * dCxx(iatm,1,:,2) ) &
                      & + ( lyy(jatm,3) * dCxx(iatm,1,:,3) )
      Else If (Omega .EQ. alf(iatm,2)) Then
        dwAjwBk_xy_da(:) = ( lyy(jatm,1) * dCxx(iatm,2,:,1) ) &
                      & + ( lyy(jatm,2) * dCxx(iatm,2,:,2) ) &
                      & + ( lyy(jatm,3) * dCxx(iatm,2,:,3) )
      Else If (Omega .EQ. alf(iatm,3)) Then
        dwAjwBk_xy_da(:) = ( lyy(jatm,1) * dCxx(iatm,3,:,1) ) &
                      & + ( lyy(jatm,2) * dCxx(iatm,3,:,2) ) &
                      & + ( lyy(jatm,3) * dCxx(iatm,3,:,3) )
      Else If (Omega .EQ. jatm) Then
         dwAjwBk_xy_da(:) = ( lxx(iatm,1) * dCyy(jatm,1,:,1) ) &
                       & + ( lxx(iatm,2) * dCyy(jatm,1,:,2) ) &
                       & + ( lxx(iatm,3) * dCyy(jatm,1,:,3) )
      Else If (Omega .EQ. alf(jatm,2)) Then
         dwAjwBk_xy_da(:) = ( lxx(iatm,1) * dCyy(jatm,2,:,1) ) &
                       & + ( lxx(iatm,2) * dCyy(jatm,2,:,2) ) &
                       & + ( lxx(iatm,3) * dCyy(jatm,2,:,3) )
      Else If (Omega .EQ. alf(jatm,3)) Then
         dwAjwBk_xy_da(:) = ( lxx(iatm,1) * dCyy(jatm,3,:,1) ) &
                       & + ( lxx(iatm,2) * dCyy(jatm,3,:,2) ) &
                       & + ( lxx(iatm,3) * dCyy(jatm,3,:,3) )
      End If

  End Function

  Function dwAjwBk_yy_da( iatm, jatm, Omega, xxx, yyy, zzz )

    Use kinds_f90

    Implicit None

    Integer, Intent( In ) :: iatm, jatm, Omega 
    Real( Kind = wp ), Intent( In ) :: xxx(*), yyy(*), zzz(*) 
    Real( Kind = wp ) :: dwAjwBk_yy_da(3)

      If (Omega .EQ. iatm) Then
        dwAjwBk_yy_da(:) = ( lyy(jatm,1) * dCyy(iatm,1,:,1) ) &
                      & + ( lyy(jatm,2) * dCyy(iatm,1,:,2) ) &
                      & + ( lyy(jatm,3) * dCyy(iatm,1,:,3) )
      Else If (Omega .EQ. alf(iatm,2)) Then
        dwAjwBk_yy_da(:) = ( lyy(jatm,1) * dCyy(iatm,2,:,1) ) &
                      & + ( lyy(jatm,2) * dCyy(iatm,2,:,2) ) &
                      & + ( lyy(jatm,3) * dCyy(iatm,2,:,3) )
      Else If (Omega .EQ. alf(iatm,3)) Then
        dwAjwBk_yy_da(:) = ( lyy(jatm,1) * dCyy(iatm,3,:,1) ) &
                      & + ( lyy(jatm,2) * dCyy(iatm,3,:,2) ) &
                      & + ( lyy(jatm,3) * dCyy(iatm,3,:,3) )
      Else If (Omega .EQ. jatm) Then
        dwAjwBk_yy_da(:) = ( lyy(iatm,1) * dCyy(jatm,1,:,1) ) &
                      & + ( lyy(iatm,2) * dCyy(jatm,1,:,2) ) &
                      & + ( lyy(iatm,3) * dCyy(jatm,1,:,3) )
      Else If (Omega .EQ. alf(jatm,2)) Then
        dwAjwBk_yy_da(:) = ( lyy(iatm,1) * dCyy(jatm,2,:,1) ) &
                      & + ( lyy(iatm,2) * dCyy(jatm,2,:,2) ) &
                      & + ( lyy(iatm,3) * dCyy(jatm,2,:,3) )
      Else If (Omega .EQ. alf(jatm,3)) Then
        dwAjwBk_yy_da(:) = ( lyy(iatm,1) * dCyy(jatm,3,:,1) ) &
                      & + ( lyy(iatm,2) * dCyy(jatm,3,:,2) ) &
                      & + ( lyy(iatm,3) * dCyy(jatm,3,:,3) )
      End If

  End Function

  Function dwAjwBk_zy_da( iatm, jatm, Omega, xxx, yyy, zzz )

    Use kinds_f90

    Implicit None

    Integer, Intent( In ) :: iatm, jatm, Omega 
    Real( Kind = wp ), Intent( In ) :: xxx(*), yyy(*), zzz(*) 
    Real( Kind = wp ) :: dwAjwBk_zy_da(3)

      If (Omega .EQ. iatm) Then
        dwAjwBk_zy_da(:) = ( lyy(jatm,1) * dCzz(iatm,1,:,1) ) &
                      & + ( lyy(jatm,2) * dCzz(iatm,1,:,2) ) &
                      & + ( lyy(jatm,3) * dCzz(iatm,1,:,3) )
      Else If (Omega .EQ. alf(iatm,2)) Then
        dwAjwBk_zy_da(:) = ( lyy(jatm,1) * dCzz(iatm,2,:,1) ) &
                      & + ( lyy(jatm,2) * dCzz(iatm,2,:,2) ) &
                      & + ( lyy(jatm,3) * dCzz(iatm,2,:,3) )
      Else If (Omega .EQ. alf(iatm,3)) Then
        dwAjwBk_zy_da(:) = ( lyy(jatm,1) * dCzz(iatm,3,:,1) ) &
                      & + ( lyy(jatm,2) * dCzz(iatm,3,:,2) ) &
                      & + ( lyy(jatm,3) * dCzz(iatm,3,:,3) )
      Else If (Omega .EQ. jatm) Then
        dwAjwBk_zy_da(:) = ( lzz(iatm,1) * dCyy(jatm,1,:,1) ) &
                      & + ( lzz(iatm,2) * dCyy(jatm,1,:,2) ) &
                      & + ( lzz(iatm,3) * dCyy(jatm,1,:,3) )
      Else If (Omega .EQ. alf(jatm,2)) Then
        dwAjwBk_zy_da(:) = ( lzz(iatm,1) * dCyy(jatm,2,:,1) ) &
                      & + ( lzz(iatm,2) * dCyy(jatm,2,:,2) ) &
                      & + ( lzz(iatm,3) * dCyy(jatm,2,:,3) )
      Else If (Omega .EQ. alf(jatm,3)) Then
        dwAjwBk_zy_da(:) = ( lzz(iatm,1) * dCyy(jatm,3,:,1) ) &
                      & + ( lzz(iatm,2) * dCyy(jatm,3,:,2) ) &
                      & + ( lzz(iatm,3) * dCyy(jatm,3,:,3) )
      End If

  End Function

  Function dwAjwBk_xz_da( iatm, jatm, Omega, xxx, yyy, zzz )

    Use kinds_f90

    Implicit None

    Integer, Intent( In ) :: iatm, jatm, Omega 
    Real( Kind = wp ), Intent( In ) :: xxx(*), yyy(*), zzz(*) 
    Real( Kind = wp ) :: dwAjwBk_xz_da(3)

      If (Omega .EQ. iatm) Then
        dwAjwBk_xz_da(:) =   ( lzz(jatm,1) * dCxx(iatm,1,:,1) ) &
                      & + ( lzz(jatm,2) * dCxx(iatm,1,:,2) ) &
                      & + ( lzz(jatm,3) * dCxx(iatm,1,:,3) )
      Else If (Omega .EQ. alf(iatm,2)) Then
        dwAjwBk_xz_da(:) =   ( lzz(jatm,1) * dCxx(iatm,2,:,1) ) &
                      & + ( lzz(jatm,2) * dCxx(iatm,2,:,2) ) &
                      & + ( lzz(jatm,3) * dCxx(iatm,2,:,3) )
      Else If (Omega .EQ. alf(iatm,3)) Then
        dwAjwBk_xz_da(:) =   ( lzz(jatm,1) * dCxx(iatm,3,:,1) ) &
                      & + ( lzz(jatm,2) * dCxx(iatm,3,:,2) ) &
                      & + ( lzz(jatm,3) * dCxx(iatm,3,:,3) )
      Else If (Omega .EQ. jatm) Then
        dwAjwBk_xz_da(:) = ( lxx(iatm,1) * dCzz(jatm,1,:,1) ) &
                      & + ( lxx(iatm,2) * dCzz(jatm,1,:,2) ) &
                      & + ( lxx(iatm,3) * dCzz(jatm,1,:,3) )
      Else If (Omega .EQ. alf(jatm,2)) Then
        dwAjwBk_xz_da(:) = ( lxx(iatm,1) * dCzz(jatm,2,:,1) ) &
                      & + ( lxx(iatm,2) * dCzz(jatm,2,:,2) ) &
                      & + ( lxx(iatm,3) * dCzz(jatm,2,:,3) )
      Else If (Omega .EQ. alf(jatm,3)) Then
        dwAjwBk_xz_da(:) = ( lxx(iatm,1) * dCzz(jatm,3,:,1) ) &
                      & + ( lxx(iatm,2) * dCzz(jatm,3,:,2) ) &
                      & + ( lxx(iatm,3) * dCzz(jatm,3,:,3) )
      End If

  End Function

  Function dwAjwBk_yz_da( iatm, jatm, Omega, xxx, yyy, zzz )

    Use kinds_f90

    Implicit None

    Integer, Intent( In ) :: iatm, jatm, Omega 
    Real( Kind = wp ), Intent( In ) :: xxx(*), yyy(*), zzz(*) 
    Real( Kind = wp ) :: dwAjwBk_yz_da(3)

      If (Omega .EQ. iatm) Then
        dwAjwBk_yz_da(:) =   ( lzz(jatm,1) * dCyy(iatm,1,:,1) ) &
                      & + ( lzz(jatm,2) * dCyy(iatm,1,:,2) ) &
                      & + ( lzz(jatm,3) * dCyy(iatm,1,:,3) )
      Else If (Omega .EQ. alf(iatm,2)) Then
        dwAjwBk_yz_da(:) =   ( lzz(jatm,1) * dCyy(iatm,2,:,1) ) &
                      & + ( lzz(jatm,2) * dCyy(iatm,2,:,2) ) &
                      & + ( lzz(jatm,3) * dCyy(iatm,2,:,3) )
      Else If (Omega .EQ. alf(iatm,3)) Then
        dwAjwBk_yz_da(:) =   ( lzz(jatm,1) * dCyy(iatm,3,:,1) ) &
                      & + ( lzz(jatm,2) * dCyy(iatm,3,:,2) ) &
                      & + ( lzz(jatm,3) * dCyy(iatm,3,:,3) )
      Else If (Omega .EQ. jatm) Then
        dwAjwBk_yz_da(:) = ( lyy(iatm,1) * dCzz(jatm,1,:,1) ) &
                      & + ( lyy(iatm,2) * dCzz(jatm,1,:,2) ) &
                      & + ( lyy(iatm,3) * dCzz(jatm,1,:,3) )
      Else If (Omega .EQ. alf(jatm,2)) Then
        dwAjwBk_yz_da(:) = ( lyy(iatm,1) * dCzz(jatm,2,:,1) ) &
                      & + ( lyy(iatm,2) * dCzz(jatm,2,:,2) ) &
                      & + ( lyy(iatm,3) * dCzz(jatm,2,:,3) )
      Else If (Omega .EQ. alf(jatm,3)) Then
        dwAjwBk_yz_da(:) = ( lyy(iatm,1) * dCzz(jatm,3,:,1) ) &
                      & + ( lyy(iatm,2) * dCzz(jatm,3,:,2) ) &
                      & + ( lyy(iatm,3) * dCzz(jatm,3,:,3) )
      End If

  End Function

  Function dwAjwBk_zz_da( iatm, jatm, Omega, xxx, yyy, zzz )

    Use kinds_f90

    Implicit None

    Integer, Intent( In ) :: iatm, jatm, Omega 
    Real( Kind = wp ), Intent( In ) :: xxx(*), yyy(*), zzz(*) 
    Real( Kind = wp ) :: dwAjwBk_zz_da(3)

      If (Omega .EQ. iatm) Then
        dwAjwBk_zz_da(:) =   ( lzz(jatm,1) * dCzz(iatm,1,:,1) ) &
                      & + ( lzz(jatm,2) * dCzz(iatm,1,:,2) ) &
                      & + ( lzz(jatm,3) * dCzz(iatm,1,:,3) )
      Else If (Omega .EQ. alf(iatm,2)) Then
        dwAjwBk_zz_da(:) =   ( lzz(jatm,1) * dCzz(iatm,2,:,1) ) &
                      & + ( lzz(jatm,2) * dCzz(iatm,2,:,2) ) &
                      & + ( lzz(jatm,3) * dCzz(iatm,2,:,3) )
      Else If (Omega .EQ. alf(iatm,3)) Then
        dwAjwBk_zz_da(:) =   ( lzz(jatm,1) * dCzz(iatm,3,:,1) ) &
                      & + ( lzz(jatm,2) * dCzz(iatm,3,:,2) ) &
                      & + ( lzz(jatm,3) * dCzz(iatm,3,:,3) )
      Else If (Omega .EQ. jatm) Then
        dwAjwBk_zz_da(:) = ( lzz(iatm,1) * dCzz(jatm,1,:,1) ) &
                      & + ( lzz(iatm,2) * dCzz(jatm,1,:,2) ) &
                      & + ( lzz(iatm,3) * dCzz(jatm,1,:,3) )
      Else If (Omega .EQ. alf(jatm,2)) Then
        dwAjwBk_zz_da(:) = ( lzz(iatm,1) * dCzz(jatm,2,:,1) ) &
                      & + ( lzz(iatm,2) * dCzz(jatm,2,:,2) ) &
                      & + ( lzz(iatm,3) * dCzz(jatm,2,:,3) )
      Else If (Omega .EQ. alf(jatm,3)) Then
        dwAjwBk_zz_da(:) = ( lzz(iatm,1) * dCzz(jatm,3,:,1) ) &
                      & + ( lzz(iatm,2) * dCzz(jatm,3,:,2) ) &
                      & + ( lzz(iatm,3) * dCzz(jatm,3,:,3) )
      End If
 
  End Function

  !=========================================  
  ! Form Partial C Matrix 
  !========================================= 
  Function Form_Partial_C_Matrix(iatm,idx)

    Use kinds_f90 

    Implicit None

    Integer, Intent( In ) :: iatm,idx
    Real ( Kind = wp ) :: Form_Partial_C_Matrix(3)

    If (idx == 1) Then 
      Form_Partial_C_Matrix(:) = lxx(iatm,:)
    End If
    If (idx == 2) Then 
      Form_Partial_C_Matrix(:) = lyy(iatm,:)
    End If
    If (idx == 3) Then 
      Form_Partial_C_Matrix(:) = lzz(iatm,:)
    End If

  End Function

  !=============================================  
  ! Form Partial dC Matrix (the derivative of C)
  !=============================================  
  Function Form_Partial_dC_Matrix(iatm,Omega,idx)

    Use kinds_f90 

    Implicit None

    Integer, Intent( In ) :: iatm, Omega, idx
    Real ( Kind = wp ) :: Form_Partial_dC_Matrix(3,3),dC1(1:3,1:3),dC2(1:3,1:3)

    If (idx == 1 ) Then
      Form_Partial_dC_Matrix(:,:) = dCxx(iatm,Omega,:,:)
      !Form_Partial_dC_Matrix = dC1k_da(iatm,Omega,lxx(iatm,:))
    End If
    If (idx == 2 ) Then
      Form_Partial_dC_Matrix(:,:) = dCyy(iatm,Omega,:,:)
      !Form_Partial_dC_Matrix = dC2k_da(iatm,Omega)
    End If
    If (idx == 3 ) Then
      Form_Partial_dC_Matrix(:,:) = dCzz(iatm,Omega,:,:)
      !dC1 = dC1k_da(iatm,Omega,lxx(iatm,:))
      !dC2 = dC2k_da(iatm,Omega)
      !Form_Partial_dC_Matrix = dC3k_da(iatm,Omega,lxx(iatm,:),lyy(iatm,:),dC1,dC2)
    End If

  End Function

End Module
