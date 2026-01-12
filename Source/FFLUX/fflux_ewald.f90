Subroutine fflux_ewald(rcut,alpha,epsq,engcpe_rc,vircpe_rc,engcpe_rl,vircpe_rl,engcpe_ex,vircpe_ex,stress,keyens,megfrz,nsteql,nstfce,nstep,rlnk)
 
  Use kinds_f90
  Use config_module, Only : natms,chge,xxx,yyy,zzz,list,imcon,cell, &
                          & fxx,fyy,fzz,ltg,nlast,lsi,lsa
  Use setup_module
  Use comms_module!,  Only : mxnode,gsum
  Use fflux_module,  Only : maxFeatures,poleOrder,alf,lEleConv,model_list, &
                          & non_alf_atms,max_atms_model,atm_label_alias, &
                          & conv_B2A,dipole_deriv,mpi_list,rt5,rt6,rt10, &
                          & rt15,rt35,rt70,rt2_3,rt3_8,rt5_8,rt5_12,rt1_24
  Use ewald_module
  Use mpoles_module, Only : mplgfr
  Use site_module,   Only : chgsit
  Use omp_lib

  Implicit None

  Real( Kind = wp ),                        Intent( In    ) :: alpha,epsq,rlnk
  Real( Kind = wp ),                        Intent( InOut ) :: engcpe_rc,vircpe_rc,engcpe_rl,vircpe_rl,engcpe_ex,vircpe_ex,rcut
  Real( Kind = wp ), Dimension(1:9),        Intent( InOut ) :: stress

  Integer, Intent( In    )       :: keyens,megfrz,nsteql,nstfce,nstep
  Integer                        :: i,j,k,l,limit,temp,indx,counter,fail,local_index,mpi_counter,chunk,ind,arrSize
  Real( Kind = wp )              :: Q,engacc,viracc,diff,talp2,alpsqrpi,dQ_dx,dQ_dy,dQ_dz,tmp(3), &
                                  & dq20(3),dq21c(3),dq21s(3),dq22c(3),dq22s(3),B2A_sq,dq30(3), &
                                  & dq31c(3),dq31s(3),dq32c(3),dq32s(3),dq33c(3),dq33s(3),dq40(3), &
                                  & dq41c(3),dq41s(3),dq42c(3),dq42s(3),dq43c(3),dq43s(3),dq44c(3),&
                                  & dq44s(3),strs1,strs2,strs3,strs5,strs6,strs9,xt,yt,zt,xtmp,ytmp,ztmp
  Real( Kind = wp ), Allocatable :: dQ_df(:),xxt(:),yyt(:),zzt(:),rrt(:),dQ_da(:,:,:,:),CartPoles(:,:), &
                                  & fxc(:),fyc(:),fzc(:),dQG_da(:,:,:,:),efx(:),efy(:),efz(:), &
                                  & SphericalPoles(:,:)


! Real space Ewald variables - saved becuase they are calculate once only.
  Logical,                  Save :: newjob = .TRUE.
  Real( Kind = wp ),        Save :: drewd,rdrewd,alp2,co1,co2,co3,co4,co5,exclcoef, &
                                    twzz,twtwz,fozz
  Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer

  Allocate(dQ_df(maxFeatures))
  If (mximpl .EQ. 0) Then
    Allocate(dQ_da(max_atms_model,nlast,1,3))
    Allocate(dQG_da(max_atms_model,nlast,1,3))
  Else
    Allocate(dQ_da(max_atms_model,nlast,mximpl,3))
    Allocate(dQG_da(max_atms_model,nlast,mximpl,3))
  End If
  Allocate(xxt(1:mxlist))
  Allocate(yyt(1:mxlist))
  Allocate(zzt(1:mxlist))
  Allocate(rrt(1:mxlist))
  Allocate(CartPoles(mximpl,mxatms))
  Allocate(efx(1:natms))
  Allocate(efy(1:natms))
  Allocate(efz(1:natms))

  mpi_counter = 0
  Do i=1,Size(mpi_list)
    If (mpi_list(i) .EQ. 0) Then
      exit
    Else
      mpi_counter = mpi_counter + 1
    End If
  End Do
! FFLUX force arrays for forces that need to be communicated between domains.
! Mpi only and in the special case of a molecule that is split across domains.
  Allocate(fxc(mpi_counter))
  Allocate(fyc(mpi_counter))
  Allocate(fzc(mpi_counter))

  Allocate(alpha_pow(0:mxompl+mxompl))

  If (PoleOrder .EQ. 1) Then
    Allocate(SphericalPoles(1,natms))
    arrSize = 1
  Else If (PoleOrder .EQ. 2) Then
    Allocate(SphericalPoles(4,natms))
    arrSize = 4
  Else If (PoleOrder .EQ. 3) Then
    Allocate(SphericalPoles(9,natms))
    arrSize = 9
  Else If (PoleOrder .EQ. 4) Then
    Allocate(SphericalPoles(16,natms))
    arrSize = 16
  Else
    Allocate(SphericalPoles(25,natms))
    arrSize = 25
  End If

  If (mxompl .GT. 0) Then
    alpha_pow(0) = 1.0_wp
    alpha_pow(1) = alpha
    Do i=2,mxompl+mxompl
      alpha_pow(i) = alpha_pow(i-1)*alpha
    End Do
  End If

  fxc = 0.0_wp; fyc = 0.0_wp; fzc = 0.0_wp
! FFLUX force accumulators - needed for stress tensor.
  efx = 0.0_wp; efy = 0.0_wp; efz = 0.0_wp

  chge = 0.0_wp

  B2A_sq = Conv_B2A*Conv_B2A

! Get derivatives of Q wrt global positions
! Derivatives are wrt position of i
! -> Deriv is wrt i i.e. don't need 1 to nlast as we don't need forces on halo
! atoms. However, we need Q \exists 1 to nlast. 
! Note still double predicting charges - change this.
! Note derivs are of local frame moments.
  Call Dlp2Atm_pos()
  
! Currently L'=0 has redundant work (MPI) because of this.
! Need to change self + reci sum MPI to get rid of it (see L'=1,2)
  If ((mxnode > 1) .AND. (poleOrder .EQ. 1)) Then
    limit = nlast
  Else
    limit = natms
  End If

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP& SHARED(alf,dQ_da,dipole_deriv,conv_B2A,SphericalPoles,natms,chge) &
!$OMP& SHARED(model_list,poleOrder,ltg,non_alf_atms,nlast,lsi,lsa,limit,mximpl) &
!$OMP& SHARED(B2A_sq,rt5_8,rt5_12,rt3_8,rt1_24,rt2_3,rt5,rt35,rt10,rt70) &
!$OMP& PRIVATE(i,j,k,l,Q,dQ_dx,dQ_dy,dQ_dz,dq20,dq21c,dq21s,dq22c,dq22s) &
!$OMP& PRIVATE(dq30,dq31c,dq31s,dq32c,dq32s,dq33c,dq33s,dq40,dq41c,dq41s)&
!$OMP& PRIVATE(dq42c,dq42s,dq43c,dq43s,dq44s,dq44c)

  If (mximpl .EQ. 0) Then
    Do i=1,limit
      !Loop over alf atoms
      Do j=1,3
        Call fflux_derivs(alf(i,j),i,2,Q,dQ_dx,dQ_dy,dQ_dz)
        dQ_da(j,i,1,1) = dQ_dx/conv_B2A
        dQ_da(j,i,1,2) = dQ_dy/conv_B2A
        dQ_da(j,i,1,3) = dQ_dz/conv_B2A
        If (j .EQ. 1) chge(i) = Q
      End Do
      !Loop over non_alf_atms
      Do j=1,model_list(ltg(i))%ptr%natms_model-3
        l = local_index(non_alf_atms(ltg(i),j), nlast,lsi,lsa)
        Call fflux_derivs(l,i,2,Q,dQ_dx,dQ_dy,dQ_dz)
        dQ_da(j+3,i,1,1) = dQ_dx/conv_B2A
        dQ_da(j+3,i,1,2) = dQ_dy/conv_B2A
        dQ_da(j+3,i,1,3) = dQ_dz/conv_B2A
      End Do
    End Do
  Else If (poleOrder .GE. 2) Then
    !$OMP DO 
    Do i=1,limit
      !Loop over moments.
      Do k=1,4
        !Loop over alf atoms
        Do j=1,3
          Call fflux_derivs(alf(i,j),i,k+1,Q,dQ_dx,dQ_dy,dQ_dz)
          dQ_da(j,i,dipole_deriv(k),1) = dQ_dx/conv_B2A
          dQ_da(j,i,dipole_deriv(k),2) = dQ_dy/conv_B2A
          dQ_da(j,i,dipole_deriv(k),3) = dQ_dz/conv_B2A
          If ((j .EQ. 1) .AND. (i <= natms)) SphericalPoles(k,i) = Q
        End Do
        !Loop over non_alf_atms
        Do j=1,model_list(ltg(i))%ptr%natms_model-3
          l = local_index(non_alf_atms(ltg(i),j), nlast,lsi,lsa)
          Call fflux_derivs(l,i,k+1,Q,dQ_dx,dQ_dy,dQ_dz)
          dQ_da(j+3,i,dipole_deriv(k),1) = dQ_dx/conv_B2A
          dQ_da(j+3,i,dipole_deriv(k),2) = dQ_dy/conv_B2A
          dQ_da(j+3,i,dipole_deriv(k),3) = dQ_dz/conv_B2A

        End Do
      End Do
    End Do
    !$OMP END DO
  End If
  !Local quadrupole moment derivatives.
  If (PoleOrder .GE. 3) Then
   !$OMP DO
    Do i=1,limit
      !Loop over alf atoms.
      Do j=1,3

        !Derivatives of local spherical moments.
        Call fflux_derivs(alf(i,j),i,6,Q,dQ_dx,dQ_dy,dQ_dz)
        dq20(1)  = dQ_dx; dq20(2)  = dQ_dy; dq20(3)  = dQ_dz
        If ((j .EQ. 1) .AND. (i <= natms)) SphericalPoles(5,i) = Q

        Call fflux_derivs(alf(i,j),i,7,Q,dQ_dx,dQ_dy,dQ_dz)
        dq21c(1) = dQ_dx; dq21c(2) = dQ_dy; dq21c(3) = dQ_dz
        If ((j .EQ. 1) .AND. (i <= natms)) SphericalPoles(6,i) = Q

        Call fflux_derivs(alf(i,j),i,8,Q,dQ_dx,dQ_dy,dQ_dz)
        dq21s(1) = dQ_dx; dq21s(2) = dQ_dy; dq21s(3) = dQ_dz
        If ((j .EQ. 1) .AND. (i <= natms)) SphericalPoles(7,i) = Q

        Call fflux_derivs(alf(i,j),i,9,Q,dQ_dx,dQ_dy,dQ_dz)
        dq22c(1) = dQ_dx; dq22c(2) = dQ_dy; dq22c(3) = dQ_dz
        If ((j .EQ. 1) .AND. (i <= natms)) SphericalPoles(8,i) = Q

        Call fflux_derivs(alf(i,j),i,10,Q,dQ_dx,dQ_dy,dQ_dz)
        dq22s(1) = dQ_dx; dq22s(2) = dQ_dy; dq22s(3) = dQ_dz
        If ((j .EQ. 1) .AND. (i <= natms)) SphericalPoles(9,i) = Q

        !Convert to local Cartesian
        dQ_da(j,i,5,:) = (1.0_wp/3.0_wp) * ( 0.5_wp * ( (rt3*dq22c(:)) - dq20(:) ))*Conv_B2A
        dQ_da(j,i,6,:) = (1.0_wp/3.0_wp) * ( 0.5_wp * rt3 * dq22s(:) )*Conv_B2A
        dQ_da(j,i,7,:) = (1.0_wp/3.0_wp) * ( 0.5_wp * rt3 * dq21c(:) )*Conv_B2A
        dQ_da(j,i,8,:) = (1.0_wp/3.0_wp) * (-0.5_wp * ( (rt3*dq22c(:)) + dq20(:) ))*Conv_B2a
        dQ_da(j,i,9,:) = (1.0_wp/3.0_wp) * ( 0.5_wp * rt3 * dq21s(:) )*Conv_B2A
        dQ_da(j,i,10,:)= (1.0_wp/3.0_wp) * dq20(:) * Conv_B2A

      End Do
      !Loop over non-alf atoms.
      Do j=1,model_list(ltg(i))%ptr%natms_model-3

        l = local_index(non_alf_atms(ltg(i),j), nlast,lsi,lsa)
        !Derivatives of local spherical moments.
        Call fflux_derivs(l,i,6,Q,dQ_dx,dQ_dy,dQ_dz)
        dq20(1)  = dQ_dx; dq20(2)  = dQ_dy; dq20(3)  = dQ_dz
                                                                                             
        Call fflux_derivs(l,i,7,Q,dQ_dx,dQ_dy,dQ_dz)
        dq21c(1) = dQ_dx; dq21c(2) = dQ_dy; dq21c(3) = dQ_dz
                                                                                             
        Call fflux_derivs(l,i,8,Q,dQ_dx,dQ_dy,dQ_dz)
        dq21s(1) = dQ_dx; dq21s(2) = dQ_dy; dq21s(3) = dQ_dz
                                                                                             
        Call fflux_derivs(l,i,9,Q,dQ_dx,dQ_dy,dQ_dz)
        dq22c(1) = dQ_dx; dq22c(2) = dQ_dy; dq22c(3) = dQ_dz
                                                                                             
        Call fflux_derivs(l,i,10,Q,dQ_dx,dQ_dy,dQ_dz)
        dq22s(1) = dQ_dx; dq22s(2) = dQ_dy; dq22s(3) = dQ_dz
                                                                                             
        !Convert to local Cartesian
        dQ_da(j+3,i,5,:) = (1.0_wp/3.0_wp) * ( 0.5_wp * ( (rt3*dq22c(:)) - dq20(:) ))*Conv_B2A
        dQ_da(j+3,i,6,:) = (1.0_wp/3.0_wp) * ( 0.5_wp * rt3 * dq22s(:) )*Conv_B2A
        dQ_da(j+3,i,7,:) = (1.0_wp/3.0_wp) * ( 0.5_wp * rt3 * dq21c(:) )*Conv_B2A
        dQ_da(j+3,i,8,:) = (1.0_wp/3.0_wp) * (-0.5_wp * ( (rt3*dq22c(:)) + dq20(:) ))*Conv_B2a
        dQ_da(j+3,i,9,:) = (1.0_wp/3.0_wp) * ( 0.5_wp * rt3 * dq21s(:) )*Conv_B2A
        dQ_da(j+3,i,10,:)= (1.0_wp/3.0_wp) * dq20(:) * Conv_B2A

      End Do
    End Do
   !$OMP END DO
  End If
  !Octupole moments
  If (PoleOrder .GE. 4) Then
   !$OMP DO
    Do i=1,limit
      !Loop over ALF
      Do j=1,3

        !Derivatives of local spherical moments.
        Call fflux_derivs(alf(i,j),i,11,Q,dQ_dx,dQ_dy,dQ_dz)
        dq30(1) = dQ_dx; dq30(2) = dQ_dy; dq30(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(10,i) = Q

        Call fflux_derivs(alf(i,j),i,12,Q,dQ_dx,dQ_dy,dQ_dz)
        dq31c(1) = dQ_dx; dq31c(2) = dQ_dy; dq31c(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(11,i) = Q

        Call fflux_derivs(alf(i,j),i,13,Q,dQ_dx,dQ_dy,dQ_dz)
        dq31s(1) = dQ_dx; dq31s(2) = dQ_dy; dq31s(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(12,i) = Q

        Call fflux_derivs(alf(i,j),i,14,Q,dQ_dx,dQ_dy,dQ_dz)
        dq32c(1) = dQ_dx; dq32c(2) = dQ_dy; dq32c(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(13,i) = Q

        Call fflux_derivs(alf(i,j),i,15,Q,dQ_dx,dQ_dy,dQ_dz)
        dq32s(1) = dQ_dx; dq32s(2) = dQ_dy; dq32s(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(14,i) = Q

        Call fflux_derivs(alf(i,j),i,16,Q,dQ_dx,dQ_dy,dQ_dz)
        dq33c(1) = dQ_dx; dq33c(2) = dQ_dy; dq33c(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(15,i) = Q

        Call fflux_derivs(alf(i,j),i,17,Q,dQ_dx,dQ_dy,dQ_dz)
        dq33s(1) = dQ_dx; dq33s(2) = dQ_dy; dq33s(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(16,i) = Q

        dQ_da(j,i,11,:) = (1.0_wp/15.0_wp)*((rt5_8 *dq33c(:)) - (rt3_8 *dq31c(:)))*B2A_sq !xxx
        dQ_da(j,i,12,:) = (1.0_wp/15.0_wp)*((rt5_8 *dq33s(:)) - (rt1_24*dq31s(:)))*B2A_sq !xxy
        dQ_da(j,i,13,:) = (1.0_wp/15.0_wp)*((rt5_12*dq32c(:)) - (0.5_wp*dq30(:)))*B2A_sq  !xxz
        dQ_da(j,i,14,:) =-(1.0_wp/15.0_wp)*((rt5_8 *dq33c(:)) + (rt1_24*dq31c(:)))*B2A_sq !xyy
        dQ_da(j,i,15,:) = (1.0_wp/15.0_wp)*( rt5_12*dq32s(:))*B2A_sq                      !xyz
        dQ_da(j,i,16,:) = (1.0_wp/15.0_wp)*( rt2_3 *dq31c(:))*B2A_sq                      !xzz
        dQ_da(j,i,17,:) =-(1.0_wp/15.0_wp)*((rt5_8 *dq33s(:)) + (rt3_8 *dq31s(:)))*B2A_sq !yyy
        dQ_da(j,i,18,:) =-(1.0_wp/15.0_wp)*((rt5_12*dq32c(:)) + (0.5_wp*dq30(:)))*B2A_sq  !yyz
        dQ_da(j,i,19,:) = (1.0_wp/15.0_wp)*( rt2_3 *dq31s(:))*B2A_sq                      !yzz
        dQ_da(j,i,20,:) = (1.0_wp/15.0_wp)*dq30(:)*B2A_sq                                 !zzz

      End Do
      !Loop over non-ALF atoms
      Do j=1,model_list(ltg(i))%ptr%natms_model-3

        l = local_index(non_alf_atms(ltg(i),j), nlast,lsi,lsa)

        !Derivatives of local spherical moments.
        Call fflux_derivs(l,i,11,Q,dQ_dx,dQ_dy,dQ_dz)
        dq30(1) = dQ_dx; dq30(2) = dQ_dy; dq30(3) = dQ_dz

        Call fflux_derivs(l,i,12,Q,dQ_dx,dQ_dy,dQ_dz)
        dq31c(1) = dQ_dx; dq31c(2) = dQ_dy; dq31c(3) = dQ_dz

        Call fflux_derivs(l,i,13,Q,dQ_dx,dQ_dy,dQ_dz)
        dq31s(1) = dQ_dx; dq31s(2) = dQ_dy; dq31s(3) = dQ_dz

        Call fflux_derivs(l,i,14,Q,dQ_dx,dQ_dy,dQ_dz)
        dq32c(1) = dQ_dx; dq32c(2) = dQ_dy; dq32c(3) = dQ_dz

        Call fflux_derivs(l,i,15,Q,dQ_dx,dQ_dy,dQ_dz)
        dq32s(1) = dQ_dx; dq32s(2) = dQ_dy; dq32s(3) = dQ_dz

        Call fflux_derivs(l,i,16,Q,dQ_dx,dQ_dy,dQ_dz)
        dq33c(1) = dQ_dx; dq33c(2) = dQ_dy; dq33c(3) = dQ_dz

        Call fflux_derivs(l,i,17,Q,dQ_dx,dQ_dy,dQ_dz)
        dq33s(1) = dQ_dx; dq33s(2) = dQ_dy; dq33s(3) = dQ_dz

        dQ_da(j+3,i,11,:) = (1.0_wp/15.0_wp)*((rt5_8 *dq33c(:)) - (rt3_8 *dq31c(:)))*B2A_sq !xxx
        dQ_da(j+3,i,12,:) = (1.0_wp/15.0_wp)*((rt5_8 *dq33s(:)) - (rt1_24*dq31s(:)))*B2A_sq !xxy
        dQ_da(j+3,i,13,:) = (1.0_wp/15.0_wp)*((rt5_12*dq32c(:)) - (0.5_wp*dq30(:)))*B2A_sq  !xxz
        dQ_da(j+3,i,14,:) =-(1.0_wp/15.0_wp)*((rt5_8 *dq33c(:)) + (rt1_24*dq31c(:)))*B2A_sq !xyy
        dQ_da(j+3,i,15,:) = (1.0_wp/15.0_wp)*( rt5_12*dq32s(:))*B2A_sq                      !xyz
        dQ_da(j+3,i,16,:) = (1.0_wp/15.0_wp)*( rt2_3 *dq31c(:))*B2A_sq                      !xzz
        dQ_da(j+3,i,17,:) =-(1.0_wp/15.0_wp)*((rt5_8 *dq33s(:)) + (rt3_8 *dq31s(:)))*B2A_sq !yyy
        dQ_da(j+3,i,18,:) =-(1.0_wp/15.0_wp)*((rt5_12*dq32c(:)) + (0.5_wp*dq30(:)))*B2A_sq  !yyz
        dQ_da(j+3,i,19,:) = (1.0_wp/15.0_wp)*( rt2_3 *dq31s(:))*B2A_sq                      !yzz
        dQ_da(j+3,i,20,:) = (1.0_wp/15.0_wp)*dq30(:)*B2A_sq                                 !zzz

      End Do
    End Do
   !$OMP END DO
  End If
  !Hexadecapole moments
  If (PoleOrder .EQ. 5) Then
   !$OMP DO
    Do i=1,limit
      !ALF
      Do j=1,3

        !Derivatives of local spherical moments.
        Call fflux_derivs(alf(i,j),i,18,Q,dQ_dx,dQ_dy,dQ_dz)
        dq40(1) = dQ_dx; dq40(2) = dQ_dy; dq40(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(17,i) = Q

        Call fflux_derivs(alf(i,j),i,19,Q,dQ_dx,dQ_dy,dQ_dz)
        dq41c(1) = dQ_dx; dq41c(2) = dQ_dy; dq41c(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(18,i) = Q

        Call fflux_derivs(alf(i,j),i,20,Q,dQ_dx,dQ_dy,dQ_dz)
        dq41s(1) = dQ_dx; dq41s(2) = dQ_dy; dq41s(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(19,i) = Q

        Call fflux_derivs(alf(i,j),i,21,Q,dQ_dx,dQ_dy,dQ_dz)
        dq42c(1) = dQ_dx; dq42c(2) = dQ_dy; dq42c(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(20,i) = Q

        Call fflux_derivs(alf(i,j),i,22,Q,dQ_dx,dQ_dy,dQ_dz)
        dq42s(1) = dQ_dx; dq42s(2) = dQ_dy; dq42s(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(21,i) = Q

        Call fflux_derivs(alf(i,j),i,23,Q,dQ_dx,dQ_dy,dQ_dz)
        dq43c(1) = dQ_dx; dq43c(2) = dQ_dy; dq43c(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(22,i) = Q

        Call fflux_derivs(alf(i,j),i,24,Q,dQ_dx,dQ_dy,dQ_dz)
        dq43s(1) = dQ_dx; dq43s(2) = dQ_dy; dq43s(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(23,i) = Q

        Call fflux_derivs(alf(i,j),i,25,Q,dQ_dx,dQ_dy,dQ_dz)
        dq44c(1) = dQ_dx; dq44c(2) = dQ_dy; dq44c(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(24,i) = Q

        Call fflux_derivs(alf(i,j),i,26,Q,dQ_dx,dQ_dy,dQ_dz)
        dq44s(1) = dQ_dx; dq44s(2) = dQ_dy; dq44s(3) = dQ_dz
        If (j .EQ. 1) SphericalPoles(25,i) = Q

        dQ_da(j,i,21,:) = ( ((3.0_wp/8.0_wp)*dq40(:)) - (0.25_wp*rt5*dq42c(:)) + ((rt35*dq44c(:))/8.0_wp) ) !xxxx
        dQ_da(j,i,22,:) = ( 0.125_wp*((rt35*dq44s(:)) - (rt5*dq42s(:))) )                                   !xxxy
        dQ_da(j,i,23,:) = ( 0.0625_wp*((rt70*dq43c(:)) - (3.0_wp*rt10*dq41c(:))) )                          !xxxz
        dQ_da(j,i,24,:) = ( (0.125_wp*dq40(:)) - (0.125_wp*rt35*dq44c(:)) )                                 !xxyy
        dQ_da(j,i,25,:) = ( 0.0625_wp*((rt70*dq43s(:)) - (rt10*dq41s(:))) )                                 !xxyz
        dQ_da(j,i,26,:) = ( 0.5_wp*((0.5_wp*rt5*dq42c(:)) - dq40(:)) )                                      !xxzz
        dQ_da(j,i,27,:) = (-0.125_wp*((rt5*dq42s(:)) + (rt35*dq44s(:))) )                                   !xyyy
        dQ_da(j,i,28,:) = (-0.0625_wp*((rt10*dq41c(:)) + (rt70*dq43c(:))) )                                 !xyyz
        dQ_da(j,i,29,:) = ( 0.25_wp*rt5*dq42s(:) )                                                          !xyzz
        dQ_da(j,i,30,:) = ( rt5_8*dq41c(:) )                                                                !xzzz
        dQ_da(j,i,31,:) = ( ((3.0_wp/8.0_wp)*dq40(:)) + (0.25_wp*rt5*dq42c(:)) + ((rt35*dq44c(:))/8.0_wp) ) !yyyy
        dQ_da(j,i,32,:) = (-0.0625_wp*((3.0_wp*rt10*dq41s(:)) + (rt70*dq43s(:))) )                          !yyyz
        dQ_da(j,i,33,:) = (-0.5_wp*((0.5_wp*rt5*dq42c(:)) + dq40(:)) )                                      !yyzz
        dQ_da(j,i,34,:) = ( rt5_8*dq41s(:) )                                                                !yzzz
        dQ_da(j,i,35,:) = dq40(:)                                                                           !zzzz

        dQ_da(j,i,21:35,:) = (1.0_wp/105.0_wp)*B2A_sq*Conv_B2A*dQ_da(j,i,21:35,:)

      End Do
      !NON-ALF
      Do j=1,model_list(ltg(i))%ptr%natms_model-3

        l = local_index(non_alf_atms(ltg(i),j), nlast,lsi,lsa)

        !Derivatives of local spherical moments.
        Call fflux_derivs(l,i,18,Q,dQ_dx,dQ_dy,dQ_dz)
        dq40(1) = dQ_dx; dq40(2) = dQ_dy; dq40(3) = dQ_dz

        Call fflux_derivs(l,i,19,Q,dQ_dx,dQ_dy,dQ_dz)
        dq41c(1) = dQ_dx; dq41c(2) = dQ_dy; dq41c(3) = dQ_dz

        Call fflux_derivs(l,i,20,Q,dQ_dx,dQ_dy,dQ_dz)
        dq41s(1) = dQ_dx; dq41s(2) = dQ_dy; dq41s(3) = dQ_dz

        Call fflux_derivs(l,i,21,Q,dQ_dx,dQ_dy,dQ_dz)
        dq42c(1) = dQ_dx; dq42c(2) = dQ_dy; dq42c(3) = dQ_dz

        Call fflux_derivs(l,i,22,Q,dQ_dx,dQ_dy,dQ_dz)
        dq42s(1) = dQ_dx; dq42s(2) = dQ_dy; dq42s(3) = dQ_dz

        Call fflux_derivs(l,i,23,Q,dQ_dx,dQ_dy,dQ_dz)
        dq43c(1) = dQ_dx; dq43c(2) = dQ_dy; dq43c(3) = dQ_dz

        Call fflux_derivs(l,i,24,Q,dQ_dx,dQ_dy,dQ_dz)
        dq43s(1) = dQ_dx; dq43s(2) = dQ_dy; dq43s(3) = dQ_dz

        Call fflux_derivs(l,i,25,Q,dQ_dx,dQ_dy,dQ_dz)
        dq44c(1) = dQ_dx; dq44c(2) = dQ_dy; dq44c(3) = dQ_dz

        Call fflux_derivs(l,i,26,Q,dQ_dx,dQ_dy,dQ_dz)
        dq44s(1) = dQ_dx; dq44s(2) = dQ_dy; dq44s(3) = dQ_dz

        dQ_da(j+3,i,21,:) = ( ((3.0_wp/8.0_wp)*dq40(:)) - (0.25_wp*rt5*dq42c(:)) + ((rt35*dq44c(:))/8.0_wp) ) !xxxx
        dQ_da(j+3,i,22,:) = ( 0.125_wp*((rt35*dq44s(:)) - (rt5*dq42s(:))) )                                   !xxxy
        dQ_da(j+3,i,23,:) = ( 0.0625_wp*((rt70*dq43c(:)) - (3.0_wp*rt10*dq41c(:))) )                          !xxxz
        dQ_da(j+3,i,24,:) = ( (0.125_wp*dq40(:)) - (0.125_wp*rt35*dq44c(:)) )                                 !xxyy
        dQ_da(j+3,i,25,:) = ( 0.0625_wp*((rt70*dq43s(:)) - (rt10*dq41s(:))) )                                 !xxyz
        dQ_da(j+3,i,26,:) = ( 0.5_wp*((0.5_wp*rt5*dq42c(:)) - dq40(:)) )                                      !xxzz
        dQ_da(j+3,i,27,:) = (-0.125_wp*((rt5*dq42s(:)) + (rt35*dq44s(:))) )                                   !xyyy
        dQ_da(j+3,i,28,:) = (-0.0625_wp*((rt10*dq41c(:)) + (rt70*dq43c(:))) )                                 !xyyz
        dQ_da(j+3,i,29,:) = ( 0.25_wp*rt5*dq42s(:) )                                                          !xyzz
        dQ_da(j+3,i,30,:) = ( rt5_8*dq41c(:) )                                                                !xzzz
        dQ_da(j+3,i,31,:) = ( ((3.0_wp/8.0_wp)*dq40(:)) + (0.25_wp*rt5*dq42c(:)) + ((rt35*dq44c(:))/8.0_wp) ) !yyyy
        dQ_da(j+3,i,32,:) = (-0.0625_wp*((3.0_wp*rt10*dq41s(:)) + (rt70*dq43s(:))) )                          !yyyz
        dQ_da(j+3,i,33,:) = (-0.5_wp*((0.5_wp*rt5*dq42c(:)) + dq40(:)) )                                      !yyzz
        dQ_da(j+3,i,34,:) = ( rt5_8*dq41s(:) )                                                                !yzzz
        dQ_da(j+3,i,35,:) = dq40(:)                                                                           !zzzz

        dQ_da(j+3,i,21:35,:) = (1.0_wp/105.0_wp)*B2A_sq*Conv_B2A*dQ_da(j+3,i,21:35,:)

      End Do

    End Do
   !$OMP END DO
  End If
 !$OMP END PARALLEL

  If (poleOrder .GT. 1) Then
    Call fflux_rotate_moments(nstep,CartPoles,SphericalPoles,arrSize)
  End If

  Call Atm2Dlp_pos()

 !Avoid redundant computation in serial.
  If (mxnode .EQ. 1) Then
    Do i=natms+1,nlast
      dQ_da(:,i,:,:) = dQ_da(:,ltg(i),:,:)
      If (poleOrder == 1) chge(i) = chge(ltg(i))
    End Do
  End If

  If (mxompl .EQ. 0) Then
    limit = 0
  Else If (mxompl .EQ. 1) Then
    limit = 3
  Else If (mxompl .EQ. 2) Then
    limit = 9
  Else If (mxompl .EQ. 3) Then
    limit = 19
  Else If (mxompl .EQ. 4) Then
    limit = 34
  End If

  !If ((mxnode > 1) .AND. (poleOrder .EQ. 1)) Then
  !  l = nlast
  !Else
  !  l = natms
  !End If
  l = natms

  !Now compute and store the global Cartesian multipole moments.
  Do i=1,l
    !Loop over ALF
    Do j=1,3

      dQG_da(j,i,1,:) = dQ_da(j,i,1,:)

      Do k=1,limit

        Call fflux_cmoment_derivs(i,j,k,dQ_dx,dQ_dy,dQ_dz,dQ_da,CartPoles)

        dQG_da(j,i,k+1,1) = dQ_dx
        dQG_da(j,i,k+1,2) = dQ_dy
        dQG_da(j,i,k+1,3) = dQ_dz

      End Do
    End Do
    !Loop over non-ALF
    Do j=1,model_list(ltg(i))%ptr%natms_model-3

      dQG_da(j+3,i,1,:) = dQ_da(j+3,i,1,:)

      Do k=1,limit

        Call fflux_cmoment_derivs(i,j+3,k,dQ_dx,dQ_dy,dQ_dz,dQ_da,CartPoles)

        dQG_da(j+3,i,k+1,1) = dQ_dx
        dQG_da(j+3,i,k+1,2) = dQ_dy
        dQG_da(j+3,i,k+1,3) = dQ_dz

      End Do
    End Do
  End Do

  If (lEleConv) Then
    If (mximpl .EQ. 0) Then
      arrSize = 1
    Else
      arrSize = mximpl
    End If
    Call fflux_electro_converge(rcut,alpha,epsq,engcpe_rc,vircpe_rc,engcpe_rl,vircpe_rl,engcpe_ex,vircpe_ex,stress,keyens,megfrz,nsteql,nstfce,nstep,dQ_da,CartPoles,fxc,fyc,fzc,mpi_counter,dQG_da,efx,efy,efz,arrSize)
  End If

  !Call the relevant Ewald routines
  Call ewald_check(keyens,megfrz,nsteql,nstfce,nstep)

  !Calculate the reciprocal space Ewald sum using the spme method
  If (poleOrder .EQ. 1) Then
    Call ewald_spme_forces(alpha,epsq,engcpe_rc,vircpe_rc,stress,dQ_da,efx,efy,efz,mpi_counter,fxc,fyc,fzc)
  Else If ((poleOrder .GT. 1) .AND. (poleOrder .LE. 3)) Then
    Call ewald_spme_mforces_d(alpha,epsq,engcpe_rc,vircpe_rc,stress,dQG_da,efx,efy,efz,mpi_counter,fxc,fyc,fzc)
  Else If (poleOrder .GE. 4) Then
    Call ewald_spme_mforces(alpha,epsq,engcpe_rc,vircpe_rc,stress,dQG_da,efx,efy,efz)
  End If

!One off calculation of necessary variables for real space Ewald.
  If (newjob) Then
    newjob = .FALSE.

    fail = 0
    Allocate(erc(0:mxgele),fer(0:mxgele), Stat=fail)
    If (fail > 0) Then
      Write(nrite,'(/,1x,a,i0)') 'ewald_real allocation failure, node:', idnode
      Call error(0)
    End If

    !interpolation interval

    drewd = rcut/Real(mxgele-4,wp)

    !Reciprocal of interpolation interval

    rdrewd = 1.0_wp/drewd

    !Generate error function complement tables for ewald sum

    Call erfcgen(rcut,alpha,mxgele,erc,fer)

    If ((poleOrder .EQ. 2) .OR. (poleOrder .EQ. 3)) Then
      !Coefficients for exponential in recurrence relation
     
      talp2 = 2.0_wp*alpha*alpha
      alpsqrpi = 1.0_wp/(alpha*sqrpi)
     
      co1 = talp2*alpsqrpi
      co2 = talp2*co1
      co3 = talp2*co2
      co4 = talp2*co3
      co5 = talp2*co4
     
      alp2 = alpha*alpha
                                                        
      exclcoef = r4pie0*alpha /sqrpi/epsq
                  
      twzz=-2.0_wp*alpha**3 *r4pie0/(3.0_wp*sqrpi*epsq)
      twtwz=4.0_wp*alpha**5 *r4pie0/(5.0_wp*sqrpi*epsq)
      fozz=12.0_wp*alpha**5 *r4pie0/(5.0_wp*sqrpi*epsq)

    End If

  End If

 !Calculate the real space Ewald sum

!!$OMP  PARALLEL DEFAULT(NONE) &
!!$OMP& PRIVATE(i,limit,k,j,xxt,yyt,zzt,rrt,engacc,viracc,mxompl,fxc,fyc,fzc,chunk,mpi_counter) &
!!$OMP& SHARED(natms,list,xxx,yyy,zzz,poleOrder,rcut,alpha,epsq,dQ_da,dQG_da) &
!!$OMP& SHARED(drewd,rdrewd,alp2,co1,co2,co3,co4,co5,exclcoef,twzz,twtwz,fozz,erc,fer) &
!!$OMP& REDUCTION(+:stress,engcpe_rl,vircpe_rl,efx,efy,efz)

! chunk = Int(natms/OMP_GET_NUM_THREADS())

! !$OMP DO SCHEDULE(STATIC,chunk)
  Do i=1,natms
   
    limit=list(0,i)

    Do k=1, limit
      j=list(k,i)
      
      xxt(k)=xxx(i)-xxx(j)
      yyt(k)=yyy(i)-yyy(j)
      zzt(k)=zzz(i)-zzz(j)

    End Do

    Do k=1,limit
      rrt(k)=Sqrt(xxt(k)**2 + yyt(k)**2 + zzt(k)**2)
    End Do
  
    If (poleOrder .EQ. 1) Then
      Call ewald_real_forces(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress,drewd,rdrewd,erc,fer,dQ_da,fxc,fyc,fzc,mpi_counter,efx,efy,efz)
    Else If ((poleOrder .GT. 1) .AND. (poleOrder .LE. 3)) Then
      Call ewald_real_mforces_d(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress,drewd,rdrewd,alp2,co1,co2,co3,co4,co5,exclcoef,twzz,twtwz,fozz,erc,fer,dQG_da,efx,efy,efz,fxc,fyc,fzc,mpi_counter)
    Else
      Call ewald_real_mforces(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress,drewd,rdrewd,erc,fer,dQG_da,efx,efy,efz)
    End If
    engcpe_rl = engcpe_rl + engacc
    vircpe_rl = vircpe_rl + viracc

  End Do
 !!$OMP END DO

!!$OMP END PARALLEL

! Calculate the ewald excluded corrections.
  Do i=1,natms
   
    limit=list(-1,i)-list(0,i)
   
    If (limit > 0) Then
      Do k=1,limit
        j=list(list(0,i)+k,i)
        
        xxt(k)=xxx(i)-xxx(j)
        yyt(k)=yyy(i)-yyy(j)
        zzt(k)=zzz(i)-zzz(j)
      End Do

      Do k=1,limit
        rrt(k)=Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2)
      End Do

      If (poleOrder .EQ. 1) Then
        Call ewald_excl_forces(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress,dQ_da,fxc,fyc,fzc,mpi_counter,efx,efy,efz)
      Else If ((poleOrder .GT. 1) .AND. (poleOrder .LE. 3)) Then
        Call ewald_excl_mforces_d(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress,dQG_da,efx,efy,efz,fxc,fyc,fzc,mpi_counter)
      Else
        Call ewald_excl_mforces(i,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engacc,viracc,stress,dQG_da,efx,efy,efz)
      End If
      engcpe_ex = engcpe_ex + engacc
      vircpe_ex = vircpe_ex + viracc

    End If
  End Do

! Need to do comms if any molecules are split across domains
  If ((mxnode > 1) .AND. (mpi_counter > 0)) Then
    Call fflux_comms(fxc,fyc,fzc,mpi_counter,rlnk,efx,efy,efz)
  End If

!Compute real end excluded sum contributions to the stress tensor.
  strs1 = 0.0_wp; strs2 = 0.0_wp; strs3 = 0.0_wp
  strs5 = 0.0_wp; strs6 = 0.0_wp; strs9 = 0.0_wp
  Do i=1,natms
    !Check if i is the heaviest atom in the ALF, if so no position correction
    !needed
    If ((i .EQ. alf(alf(i,2),2)) .AND. (i .EQ. alf(alf(i,3),2))) Then

      !Standard f.r for i, no correction as i is in main cell by definition
      vircpe_rl = vircpe_rl - ((efx(i)*xxx(i)) + (efy(i)*yyy(i)) + (efz(i)*zzz(i)))

      !Update stress tensor components.
      strs1 = strs1 + xxx(i)*efx(i)
      strs2 = strs2 + xxx(i)*efy(i)
      strs3 = strs3 + xxx(i)*efz(i)
      strs5 = strs5 + yyy(i)*efy(i)
      strs6 = strs6 + yyy(i)*efz(i)
      strs9 = strs9 + zzz(i)*efz(i)

    Else

      xt = xxx(i) - xxx(alf(i,2))
      yt = yyy(i) - yyy(alf(i,2))
      zt = zzz(i) - zzz(alf(i,2))

      xtmp = xxx(alf(i,2))
      ytmp = yyy(alf(i,2))
      ztmp = zzz(alf(i,2))
      Call images_s(imcon, cell, xtmp,ytmp,ztmp)
      xt = xxx(i) - xtmp
      yt = yyy(i) - ytmp
      zt = zzz(i) - ztmp

      Call images_s(imcon, cell, xt,yt,zt)
      xt = xt + xtmp; yt = yt + ytmp; zt = zt + ztmp;
      !Virial contribution
      vircpe_rl = vircpe_rl - ((efx(i)*xt) + (efy(i)*yt) + (efz(i)*zt))

      !Update stress tensor components.
      strs1 = strs1 + xt*efx(i)
      strs2 = strs2 + xt*efy(i)
      strs3 = strs3 + xt*efz(i)
      strs5 = strs5 + yt*efy(i)
      strs6 = strs6 + yt*efz(i)
      strs9 = strs9 + zt*efz(i)

    End If
  End Do

  stress(1) = stress(1) + strs1
  stress(2) = stress(2) + strs2
  stress(3) = stress(3) + strs3
  stress(4) = stress(4) + strs2
  stress(5) = stress(5) + strs5
  stress(6) = stress(6) + strs6
  stress(7) = stress(7) + strs3
  stress(8) = stress(8) + strs6
  stress(9) = stress(9) + strs9

!Add FFLUX forces to global force arrays.
  fxx(1:natms) = fxx(1:natms) + efx(1:natms)
  fyy(1:natms) = fyy(1:natms) + efy(1:natms)
  fzz(1:natms) = fzz(1:natms) + efz(1:natms)

! Need to do comms if any molecules are split across domains
!  If ((mxnode > 1) .AND. (mpi_counter > 0)) Then
!    Call fflux_comms(fxc,fyc,fzc,mpi_counter,rlnk)
!  End If

  !Correct for self interaction
  If (mxnode > 1) Then
    If (mximpl > 0 .and. mxompl <=2) Then
      Call gsum(engsic)
    End If
  End If

  Deallocate(dQ_df)
  Deallocate(dQ_da)
  Deallocate(dQG_da)
  Deallocate(xxt)
  Deallocate(yyt)
  Deallocate(zzt)
  Deallocate(rrt)
  Deallocate(CartPoles)
  Deallocate(fxc)
  Deallocate(fyc)
  Deallocate(fzc)
  Deallocate(alpha_pow)
  Deallocate(efx)
  Deallocate(efy)
  Deallocate(efz)
  Deallocate(SphericalPoles)

End Subroutine fflux_ewald
