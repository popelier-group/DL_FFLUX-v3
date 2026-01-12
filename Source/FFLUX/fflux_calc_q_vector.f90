
!=========================================  
! calculate the vector "q" for MultipoleMoments 
!=========================================  

Subroutine calc_q_vector(iatm,jatm,xxx,yyy,zzz,R_,rax,ray,raz,rbx,rby,rbz,Cxx,Cyx,Czx,Cxy,Cyy,Czy,Cxz,Cyz,Czz) 

  Use kinds_f90
  Use config_module, Only : imcon,cell
  Use fflux_module,  Only : Form_C_Matrix
!  Use fflux_module, Only : Form_C_Matrix, & 
!                           R_, &
!                           rax,ray,raz, &
!                           rbx,rby,rbz, &
!                           Cxx,Cyx,Czx, &
!                           Cxy,Cyy,Czy, &
!                           Cxz,Cyz,Czz

  Implicit None

  Integer,           Intent( In ) :: iatm, jatm
  Real( Kind = wp ), Intent( In ) :: xxx(*), yyy(*), zzz(*)
  Real( Kind = wp ), Intent( Out ) :: R_, rax, ray, raz, &
                                        & rbx, rby, rbz, &
                                        & Cxx, Cyx, Czx, &
                                        & Cxy, Cyy, Czy, &
                                        & Cxz, Cyz, Czz

  Real( Kind = wp ) :: RAB(3), RAB_UNITVEC(3), C_i(3,3), C_j(3,3), q(0:15) 

  RAB(1) = xxx(jatm) - xxx(iatm)
  RAB(2) = yyy(jatm) - yyy(iatm)
  RAB(3) = zzz(jatm) - zzz(iatm)

! Added by BCBS 03/2019 to account for PBC's
  Call images_s(imcon,cell,RAB(1),RAB(2),RAB(3))

  C_i = Form_C_Matrix(iatm)
  C_j = Form_C_Matrix(jatm)

  R_ = dsqrt( dot_product(RAB,RAB) )

  RAB_UNITVEC(:) = RAB(:) / R_

  rax = dot_product( RAB_UNITVEC, C_i(1,:) )
  ray = dot_product( RAB_UNITVEC, C_i(2,:) )
  raz = dot_product( RAB_UNITVEC, C_i(3,:) )

  rbx = dot_product( RAB_UNITVEC, C_j(1,:) )
  rby = dot_product( RAB_UNITVEC, C_j(2,:) )
  rbz = dot_product( RAB_UNITVEC, C_j(3,:) )

  Cxx = dot_product( C_i(1,:), C_j(1,:) )
  Cyx = dot_product( C_i(2,:), C_j(1,:) )
  Czx = dot_product( C_i(3,:), C_j(1,:) )
                                
  Cxy = dot_product( C_i(1,:), C_j(2,:) )
  Cyy = dot_product( C_i(2,:), C_j(2,:) )
  Czy = dot_product( C_i(3,:), C_j(2,:) )
                                
  Cxz = dot_product( C_i(1,:), C_j(3,:) )
  Cyz = dot_product( C_i(2,:), C_j(3,:) )
  Czz = dot_product( C_i(3,:), C_j(3,:) )

End Subroutine

Subroutine calc_q_derivs(iatm,jatm,Omega,xxx,yyy,zzz,dq_da)

  Use kinds_f90
  Use omp_lib
  Use setup_module
  Use config_module, Only : natms,imcon,cell,ltg
  Use fflux_module, Only : dRdotR_da, dwjR_da_1,dwjR_da_2,dwjR_da_3,dwAjwBk_da,&
                         & atm_label_alias,alf,lxx,lyy,lzz,dwAjwBk_xx_da,dwAjwBk_yx_da,&
                         & dwAjwBk_zx_da,dwAjwBk_xy_da,dwAjwBk_yy_da,dwAjwBk_zy_da,&
                         & dwAjwBk_xz_da,dwAjwBk_yz_da,dwAjwBk_zz_da

  Implicit None

  Integer, Intent ( In ) :: iatm, jatm, Omega
  Real( Kind = wp ), Intent( In ) :: xxx(*), yyy(*), zzz(*)
  Real( Kind = wp ), Intent( InOut ) :: dq_da(0:15,natms,3)
  Real(Kind = wp ) :: dRi_da,xdiff,ydiff,zdiff !dq_da(0:15,natms,3)
  Integer :: x, y, z,alias
  Character( Len=200 ) :: printChar
  
  If ( Omega .EQ. iatm) Then
    dRi_da = -1.0_wp
  Else If ( Omega .EQ. jatm ) Then
    dRi_da = 1.0_wp
  Else 
    dRi_da = 0.0_wp
  End If 
  
  xdiff = xxx(jatm) - xxx(iatm)
  ydiff = yyy(jatm) - yyy(iatm)
  zdiff = zzz(jatm) - zzz(iatm)

! Added by BCBS 03/2019 to account for PBC's
  Call images_s(imcon,cell,xdiff,ydiff,zdiff)  

!  dq_da(:,:,:) = 0.0_wp
  
  x = 1; y = 2; z = 3;

! Define the q derivatives
  dq_da(0,Omega,:) = dRdotR_da(iatm,jatm,Omega,xxx,yyy,zzz)                     ! R   . R
 If (mximpl .GT. 0) Then
  If ((Omega .EQ. iatm) .OR. (Omega .EQ. alf(iatm,2)) .OR. (Omega .EQ. alf(iatm,3))) Then
    dq_da(1,Omega,:) = dwjR_da_1(iatm,xdiff,ydiff,zdiff,dRi_da,Omega,xxx,yyy,zzz) ! x_A . R
    dq_da(2,Omega,:) = dwjR_da_2(iatm,xdiff,ydiff,zdiff,dRi_da,Omega,xxx,yyy,zzz) ! y_A . R
    dq_da(3,Omega,:) = dwjR_da_3(iatm,xdiff,ydiff,zdiff,dRi_da,Omega,xxx,yyy,zzz) ! z_A . R
  Else
    dq_da(1,Omega,:) = lxx(iatm,:)*dRi_da
    dq_da(2,Omega,:) = lyy(iatm,:)*dRi_da
    dq_da(3,Omega,:) = lzz(iatm,:)*dRi_da
  End If
  If ((Omega .EQ. jatm) .OR. (Omega .EQ. alf(jatm,2)) .OR. (Omega .EQ. alf(jatm,3))) Then
    dq_da(4,Omega,:) = dwjR_da_1(jatm,xdiff,ydiff,zdiff,dRi_da,Omega,xxx,yyy,zzz) ! x_B . R
    dq_da(5,Omega,:) = dwjR_da_2(jatm,xdiff,ydiff,zdiff,dRi_da,Omega,xxx,yyy,zzz) ! y_B . R
    dq_da(6,Omega,:) = dwjR_da_3(jatm,xdiff,ydiff,zdiff,dRi_da,Omega,xxx,yyy,zzz) ! z_B . R
  Else
    dq_da(4,Omega,:) = lxx(jatm,:)*dRi_da
    dq_da(5,Omega,:) = lyy(jatm,:)*dRi_da
    dq_da(6,Omega,:) = lzz(jatm,:)*dRi_da
  End If

  !alias = atm_label_alias(ltg(Omega))
  !If avoids redundant work where terms are 0 by definition.
  !iatm and jatm are in different models by definition => Omega can't belong to
  !ALF of both.
  !If (alias <= 3) Then
    !dq_da(7,Omega,:) = dwAjwBk_da(iatm,jatm,x,x,Omega,xxx,yyy,zzz)                ! x_A . x_B
    !dq_da(8,Omega,:) = dwAjwBk_da(iatm,jatm,y,x,Omega,xxx,yyy,zzz)                ! y_A . x_B
    !dq_da(9,Omega,:) = dwAjwBk_da(iatm,jatm,z,x,Omega,xxx,yyy,zzz)                ! z_A . x_B
    dq_da(7,Omega,:) = dwAjwBk_xx_da(iatm,jatm,Omega,xxx,yyy,zzz)                ! x_A . x_B
    dq_da(8,Omega,:) = dwAjwBk_yx_da(iatm,jatm,Omega,xxx,yyy,zzz)                ! y_A . x_B
    dq_da(9,Omega,:) = dwAjwBk_zx_da(iatm,jatm,Omega,xxx,yyy,zzz)                ! z_A . x_B
                                                                                           
    !dq_da(10,Omega,:) = dwAjwBk_da(iatm,jatm,x,y,Omega,xxx,yyy,zzz)                ! x_A . y_B
    dq_da(10,Omega,:) = dwAjwBk_xy_da(iatm,jatm,Omega,xxx,yyy,zzz)                ! x_A . y_B
    !dq_da(11,Omega,:) = dwAjwBk_da(iatm,jatm,y,y,Omega,xxx,yyy,zzz)                ! y_A . y_B
    dq_da(11,Omega,:) = dwAjwBk_yy_da(iatm,jatm,Omega,xxx,yyy,zzz)                ! y_A . y_B
    !dq_da(12,Omega,:) = dwAjwBk_da(iatm,jatm,z,y,Omega,xxx,yyy,zzz)                ! z_A . y_B
    dq_da(12,Omega,:) = dwAjwBk_zy_da(iatm,jatm,Omega,xxx,yyy,zzz)                ! z_A . y_B
                                                                                           
    !dq_da(13,Omega,:) = dwAjwBk_da(iatm,jatm,x,z,Omega,xxx,yyy,zzz)                ! x_A . z_B
    !dq_da(14,Omega,:) = dwAjwBk_da(iatm,jatm,y,z,Omega,xxx,yyy,zzz)                ! y_A . z_B
    !dq_da(15,Omega,:) = dwAjwBk_da(iatm,jatm,z,z,Omega,xxx,yyy,zzz)                ! z_A . z_B
    dq_da(13,Omega,:) = dwAjwBk_xz_da(iatm,jatm,Omega,xxx,yyy,zzz)                ! x_A . z_B  
    dq_da(14,Omega,:) = dwAjwBk_yz_da(iatm,jatm,Omega,xxx,yyy,zzz)                ! y_A . z_B
    dq_da(15,Omega,:) = dwAjwBk_zz_da(iatm,jatm,Omega,xxx,yyy,zzz)                ! z_A . z_B

  !Else
  !  dq_da(7,Omega,:) = 0.0_wp
  !  dq_da(8,Omega,:) = 0.0_wp
  !  dq_da(9,Omega,:) = 0.0_wp
  !                     
  !  dq_da(10,Omega,:) =0.0_wp
  !  dq_da(11,Omega,:) =0.0_wp
  !  dq_da(12,Omega,:) =0.0_wp
  !                     
  !  dq_da(13,Omega,:) =0.0_wp
  !  dq_da(14,Omega,:) =0.0_wp
  !  dq_da(15,Omega,:) =0.0_wp
  !End If
 End If

End Subroutine 
