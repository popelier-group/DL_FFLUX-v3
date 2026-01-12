Subroutine fflux_calc_forces(engmet,virmet,engcpe,vircpe,xxt,yyt,zzt,stress,rcut,epsq,nstep)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 -- FFLUX subroutine for calculating E_IQA and force terms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use omp_lib
  Use kinds_f90
  Use config_module, Only : natms,ltg,chge,fxx,fyy,fzz,xxx,yyy,zzz,imcon,cell,lsite,nlast,lsa,lsi
  Use setup_module
  Use mpoles_module, Only : mplgfr
  Use fflux_module,  Only : pwd,features,alf,dR_da,non_alf_atms,conv_Fau, &
                          & conv_Ha_2_10Jmol,conv_B2A, dCxx,dCyy,dCzz, &
                          & dC1k_da,dC2k_da,dC3k_da, lxx,lyy,lzz, C_1,C_2,C_3, &
                          & poleOrder,ewaldFlag, globalAlf, dChi_da,dTheta_da,dPhi_da, &
                          & model_list,maxFeatures,max_atms_model,atm_label_alias,clustFlag, &
                          & mpi_list,lfflux_energy,lfflux_force, CONVERGED, &
                          & fflux_check_convergence
  Use comms_module

  Implicit None

  Real( Kind = wp ), Dimension( 1:natms ),  Intent( In    ) :: xxt,yyt,zzt
  Real( Kind = wp ),                        Intent(   Out ) :: engmet,virmet,engcpe,vircpe,rcut,epsq
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
  Integer,                                  Intent( In    ) :: nstep

  Integer                                                   :: i,k,nPred,j,local_index,l1,l2,l3,cFeat,limit,l,idi,counter,ind
  Integer, Allocatable                                      :: local_non_alf(:)
  Character( LEN = 2   )                                    :: AtomStr
  Character( LEN = 225 )                                    :: filename,QQQ
  Real( Kind = wp ), Allocatable                            :: dIQA_df(:),wxx(:),wyy(:),wzz(:),dM_da(:,:,:,:),dq_da(:,:,:),cPoles(:,:),E_atm(:)
  Real( Kind = wp )                                         :: E_IQA,&
                                                               & strs1,strs2,strs3,strs5,strs6,strs9,&
                                                               & fix, qmol, aaa, rcutsq,zeta1,zeta2,zeta3,&
                                                               & xdiff,ydiff,zdiff,xt,yt,zt,strsConv,xtmp,ytmp,ztmp,&
                                                               & dIQA_dx,dIQA_dy,dIQA_dz,Q,dQ_dx,dQ_dy,dQ_dz
  Real( Kind = wp ), Allocatable, Dimension(:)              :: Localfxx,Localfyy,Localfzz,temp
  Logical                                                   :: exists

  Allocate( dIQA_df( maxFeatures ) )
  Allocate( wxx( 3 ) )
  Allocate( wyy( 3 ) )
  Allocate( wzz( 3 ) )
  Allocate( Localfxx( 1:Size(fxx) ) )
  Allocate( Localfyy( 1:Size(fyy) ) )
  Allocate( Localfzz( 1:Size(fzz) ) )
  Allocate( temp(3) )
  Allocate( local_non_alf(Size(non_alf_atms(1,:))) )
  Allocate( dq_da(0:15,natms,3) )
  Allocate( cPoles(mximpl,mxatms) )
  If (fflux_check_convergence) Then
    Allocate(E_atm(1:natms))
  End If

  If (poleOrder .EQ. 1) Then
    limit = 1
  Else If (PoleOrder .EQ. 2) Then
    limit = 4
  Else If (PoleOrder .EQ. 3) Then
    limit = 9
  Else If (PoleOrder .EQ. 4) Then
    limit = 16
  Else
    limit = 25
  End If

  Allocate( dM_da(natms,max_atms_model,limit,3) )

! initialise potential energy and virial
  engmet=0.0_wp
  virmet=0.0_wp

  engcpe=0.0_wp
  vircpe=0.0_wp

  temp=0.0_wp

! convert from dlp to atomic units
  Call Dlp2Atm(rcut,engmet,virmet,engcpe,vircpe)

  !Compute local alf
  !Do i=1,nlast
  !  alf(i,1) = i
  !  alf(i,2) = local_index(globalAlf(ltg(i),2), nlast,lsi,lsa)
  !  alf(i,3) = local_index(globalAlf(ltg(i),3), nlast,lsi,lsa)
  !End Do
 
  !Compute local ALF
  counter = 1
  ! added by MN. Changed the size of the mpi_list from natms to nlast 
  ! in order to be able to store all the atoms that need to be communicated in
  ! the fflux_comms.f90
  Deallocate(mpi_list) 
  Allocate(mpi_list(nlast)) 
    mpi_list = 0
  exists = .FALSE.
  Do i=1,nlast
    alf(i,1) = i
    alf(i,2) = local_index(globalAlf(ltg(i),2), nlast,lsi,lsa)
    alf(i,3) = local_index(globalAlf(ltg(i),3), nlast,lsi,lsa)
    If (i <= natms) Then
      If (alf(i,2) > natms) Then
        Do j=1,Size(mpi_list)
          If (alf(i,2) .EQ. mpi_list(j)) Then
            exists = .TRUE.
            exit
          Else If (mpi_list(j) .EQ. 0) Then
            exit
          End If
        End Do
        If (exists .EQV. .FALSE.) Then
          mpi_list(counter) = alf(i,2)
          counter = counter + 1
        End If
        exists = .FALSE.
      End If

      If (alf(i,3) > natms) Then
        Do j=1,Size(mpi_list)
          If (alf(i,3) .EQ. mpi_list(j)) Then
            exists = .TRUE.
            exit
          Else If (mpi_list(j) .EQ. 0) Then
            exit
          End If
        End Do
        If (exists .EQV. .FALSE.) Then
          mpi_list(counter) = alf(i,3)
          counter = counter + 1
        End If
        exists = .FALSE.
      End If

      !Non_alf
      Do j=1,model_list(ltg(i))%ptr%natms_model-3
        ind = local_index(non_alf_atms(ltg(i),j), nlast,lsi,lsa)
        If (ind > natms) Then
          Do k=1,Size(mpi_list)
            If (ind .EQ. mpi_list(k)) Then
              exists = .TRUE.
              exit
            Else If (mpi_list(k) .EQ. 0) Then
              exit
            End If
          End Do
          If (exists .EQV. .FALSE.) Then
            mpi_list(counter) = ind
            counter = counter + 1
          End If
          exists = .FALSE.
        End If
      End Do

    End If
  End Do
  
! determine rcutsq 
  rcutsq = rcut * rcut

  Localfxx(:) = 0.0_wp
  Localfyy(:) = 0.0_wp
  Localfzz(:) = 0.0_wp
!Calculate features
  Do i = 1, nlast

    !Calculate wrapped cooridnates needed for feature and C matrix calculation.

     wxx(1) = xxx(alf(i,1)) - xxx(alf(i,1))
     wxx(2) = xxx(alf(i,2)) - xxx(alf(i,1))
     wxx(3) = xxx(alf(i,3)) - xxx(alf(i,1))
                                              
     wyy(1) = yyy(alf(i,1)) - yyy(alf(i,1))
     wyy(2) = yyy(alf(i,2)) - yyy(alf(i,1))
     wyy(3) = yyy(alf(i,3)) - yyy(alf(i,1))
                                              
     wzz(1) = zzz(alf(i,1)) - zzz(alf(i,1))
     wzz(2) = zzz(alf(i,2)) - zzz(alf(i,1))
     wzz(3) = zzz(alf(i,3)) - zzz(alf(i,1))

    Call images( imcon, cell, Size(wxx), wxx, wyy, wzz )

    Call fflux_calc_features(i,wxx,wyy,wzz,nstep)

  End Do

! IQA forces
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP& SHARED(natms,lfflux_energy,non_alf_atms,alf,ltg,nlast,lsi,lsa) &
!$OMP& SHARED(Localfxx,Localfyy,Localfzz,model_list,E_atm,fflux_check_convergence) &
!$OMP& PRIVATE(i,j,E_IQA,dIQA_df,dIQA_dx,dIQA_dy,dIQA_dz) &
!$OMP& REDUCTION(+:engmet)

 !$OMP DO
  Do i=1,natms
    !Loop over alf altoms
    Call fflux_predict_value(i,1,E_IQA,dIQA_df)
    engmet = engmet + E_IQA
    If (fflux_check_convergence) Then
      E_atm(i) = E_IQA
    End If

    If (lfflux_energy) Write(nenergydt, '(i10,3f20.10)') i,E_IQA

    Do j=1,3
      Call fflux_derivs(i,alf(i,j), 1,E_IQA,dIQA_dx,dIQA_dy,dIQA_dz)
      Localfxx(i) = Localfxx(i) - dIQA_dx
      Localfyy(i) = Localfyy(i) - dIQA_dy
      Localfzz(i) = Localfzz(i) - dIQA_dz
    End Do
    !Loop over non alf atoms
    Do j=1,model_list(ltg(i))%ptr%natms_model-3
      Call fflux_derivs(i,local_index(non_alf_atms(ltg(i),j), nlast,lsi,lsa), 1,E_IQA,dIQA_dx,dIQA_dy,dIQA_dz)
      Localfxx(i) = Localfxx(i) - dIQA_dx
      Localfyy(i) = Localfyy(i) - dIQA_dy
      Localfzz(i) = Localfzz(i) - dIQA_dz
    End Do
  End Do
 !$OMP END DO
!$OMP END PARALLEL

  If (lfflux_force) Then
    Do i=1,natms
      Write(nforcedt, '(i10,3f20.10,3f20.10,3f20.10)') i,Localfxx(i),Localfyy(i),Localfzz(i)
    End Do
  End If

!old electrostatics method - clusters
 If (clustFlag) Then

  If (mximpl .EQ. 0) Then                    
    Do i=1,natms
      !Loop over alf atoms
      Do j=1,3
        Call fflux_derivs(i,alf(i,j),2,Q,dQ_dx,dQ_dy,dQ_dz)
        If (j .EQ. 1) chge(i) = Q
        dM_da(i,atm_label_alias(ltg(alf(i,j))),1,1) = dQ_dx
        dM_da(i,atm_label_alias(ltg(alf(i,j))),1,2) = dQ_dy
        dM_da(i,atm_label_alias(ltg(alf(i,j))),1,3) = dQ_dz
      End Do
      !Loop over non_alf_atms
      Do j=1,model_list(ltg(i))%ptr%natms_model-3
        l = local_index(non_alf_atms(i,j), nlast,lsi,lsa)
        Call fflux_derivs(i,l,2,Q,dQ_dx,dQ_dy,dQ_dz)
        dM_da(i,atm_label_alias(ltg(l)),1,1) = dQ_dx
        dM_da(i,atm_label_alias(ltg(l)),1,2) = dQ_dy
        dM_da(i,atm_label_alias(ltg(l)),1,3) = dQ_dz
      End Do
    End Do
  Else
    Do i=1,natms
      !Loop over moments.
      Do k=1,limit
        !Loop over alf atoms
        Do j=1,3
          Call fflux_derivs(i,alf(i,j),K+1,Q,dQ_dx,dQ_dy,dQ_dz)
          If (j .EQ. 1) mplgfr(k,i) = Q
          dM_da(i,atm_label_alias(ltg(alf(i,j))),k,1) = dQ_dx
          dM_da(i,atm_label_alias(ltg(alf(i,j))),k,2) = dQ_dy
          dM_da(i,atm_label_alias(ltg(alf(i,j))),k,3) = dQ_dz
        End Do
        !Loop over non_alf_atms
        Do j=1,model_list(ltg(i))%ptr%natms_model-3
          l = local_index(non_alf_atms(i,j), nlast,lsi,lsa)
          Call fflux_derivs(i,l,k+1,Q,dQ_dx,dQ_dy,dQ_dz)
          dM_da(i,atm_label_alias(ltg(l)),k,1) = dQ_dx
          dM_da(i,atm_label_alias(ltg(l)),k,2) = dQ_dy
          dM_da(i,atm_label_alias(ltg(l)),k,3) = dQ_dz
        End Do
      End Do
    End Do
  End If

  Do i=1,natms
    Do j=i,natms
      If (( ANY(j .EQ. alf(i,:)) .OR. ANY(j .EQ. non_alf_atms(i,:))) .EQV. .FALSE.) Then
        Call fflux_calc_mp_forces(i,j,rcut,epsq,engcpe,Localfxx,Localfyy,Localfzz,limit,dM_da,dq_da)
      End If
    End Do
  End Do

 End If

  strsConv = conv_B2A*conv_Fau
  strs1 = 0.0_wp; strs2 = 0.0_wp; strs3 = 0.0_wp
  strs5 = 0.0_wp; strs6 = 0.0_wp; strs9 = 0.0_wp


! Virial calculation corrected for PBC's + MPI
  Do i=1,natms
    !Check if i is the heaviest atom in the ALF, if so no position correction
    !needed
    If ((i .EQ. alf(alf(i,2),2)) .AND. (i .EQ. alf(alf(i,3),2))) Then

      !Standard f.r for i, no correction as i is in main cell by definition
      virmet = virmet - ((Localfxx(i)*xxx(i)) + (Localfyy(i)*yyy(i)) + (Localfzz(i)*zzz(i)))

      !Update stress tensor components.
      strs1 = strs1 + xxx(i)*Localfxx(i)
      strs2 = strs2 + xxx(i)*Localfyy(i)
      strs3 = strs3 + xxx(i)*Localfzz(i)
      strs5 = strs5 + yyy(i)*Localfyy(i)
      strs6 = strs6 + yyy(i)*Localfzz(i)
      strs9 = strs9 + zzz(i)*Localfzz(i)

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
      virmet = virmet - ((Localfxx(i)*xt) + (Localfyy(i)*yt) + (Localfzz(i)*zt))

      !Update stress tensor components.
      strs1 = strs1 + xt*Localfxx(i)
      strs2 = strs2 + xt*Localfyy(i)
      strs3 = strs3 + xt*Localfzz(i)
      strs5 = strs5 + yt*Localfyy(i)
      strs6 = strs6 + yt*Localfzz(i)
      strs9 = strs9 + zt*Localfzz(i)

    End If
  End Do

  fxx(:) = fxx(:) + Localfxx(:)
  fyy(:) = fyy(:) + Localfyy(:)
  fzz(:) = fzz(:) + Localfzz(:)

  stress(1) = stress(1) + strs1*strsConv
  stress(2) = stress(2) + strs2*strsConv
  stress(3) = stress(3) + strs3*strsConv
  stress(4) = stress(4) + strs2*strsConv
  stress(5) = stress(5) + strs5*strsConv
  stress(6) = stress(6) + strs6*strsConv
  stress(7) = stress(7) + strs3*strsConv
  stress(8) = stress(8) + strs6*strsConv
  stress(9) = stress(9) + strs9*strsConv

  If (fflux_check_convergence .EQV. .TRUE.) Then  
    Call fflux_geometry_opt_converged(E_atm, nstep)  ! todo: is virmet
  End If

! convert back from atomic to dlp units
  Call Atm2Dlp(rcut,engmet,virmet,engcpe,vircpe)

  Deallocate( dIQA_df )
  Deallocate( wxx )
  Deallocate( wyy )
  Deallocate( wzz )
  Deallocate( Localfxx )
  Deallocate( Localfyy )
  Deallocate( Localfzz )
  Deallocate( temp )
  Deallocate( local_non_alf )
  Deallocate( dM_da )
  Deallocate( dq_da )
  Deallocate( cPoles)
  If (allocated(E_atm)) Deallocate(E_atm)

End Subroutine fflux_calc_forces
