Subroutine fflux_calc_features(iatm,wxx,wyy,wzz,nstep)

  Use kinds_f90
  Use config_module, Only : natms,xxx,yyy,zzz,imcon,cell,ltg,lsi,lsa,nlast
  Use fflux_module,  Only : lxx,lyy,lzz,features,alf,C_1,C_2,C_3, &
                          & df_da,dR_da,dChi_da,dTheta_da,dPhi_da, &
                          & non_alf_atms,conv_B2A, ewaldFlag, &
                          & dCxx,dCyy,dCzz,dC1k_da,dC2k_da,dC3k_da, &
                          & model_list,atm_label_alias
  Use comms_module

  Implicit None

  Integer, Intent( In ) :: iatm, nstep
  Real( Kind = wp ), Intent( In ) :: wxx(*),wyy(*),wzz(*)
  Integer :: j, jatm, cFeat,local_index,alias,k
  Real( Kind = wp ) :: xdiff,ydiff,zdiff
  Real( Kind = wp ) :: xdiff1,ydiff1,zdiff1
  Real( Kind = wp ) :: xdiff2,ydiff2,zdiff2 
  Real( Kind = wp ) :: temp, temp2, zeta(3),t1(3),t2(3),t3(3)

    cFeat = 1

  !======================================!
  !                  1st feature, R^{Ax} !
  !======================================!

    xdiff1 = wxx(2)     
    ydiff1 = wyy(2)    
    zdiff1 = wzz(2)
    features(iatm,cFeat) = dsqrt( xdiff1*xdiff1 + ydiff1*ydiff1 + zdiff1*zdiff1 )

    cFeat = cFeat + 1

  !======================================!
  !                 2nd feature, R^{Axy} !
  !======================================!

    xdiff2 = wxx(3)     
    ydiff2 = wyy(3)    
    zdiff2 = wzz(3)     
    features(iatm,cFeat) = dsqrt( xdiff2*xdiff2 + ydiff2*ydiff2 + zdiff2*zdiff2  )

    cFeat = cFeat + 1

  !======================================!
  !                 3rd feature, chi^{A} !
  !======================================!

    temp = xdiff1*xdiff2 + ydiff1*ydiff2 + zdiff1*zdiff2
    features(iatm,cFeat) = dacos( temp / (features(iatm,1)*features(iatm,2)) ) 

    !Calculate C matrix - Moved here by BCBS as non_alf features require it.
    lxx(iatm,:) = C_1(iatm)
    lyy(iatm,:) = C_2(iatm)
    lzz(iatm,:) = C_3(iatm,lxx(iatm,:),lyy(iatm,:))

    If (iatm <= nlast) Then
      Do j=1,3
        jatm = alf(iatm,j)
        !If (jatm <= natms) Then
          dCxx(iatm,j,:,:) = dC1k_da(iatm,jatm, lxx(iatm,:))
          dCyy(iatm,j,:,:) = dC2k_da(iatm,jatm)
          dCzz(iatm,j,:,:) = dC3k_da(iatm,jatm, lxx(iatm,:),lyy(iatm,:),dCxx(iatm,j,:,:),dCyy(iatm,j,:,:))
          !alias = atm_label_alias(ltg(jatm))
          !If (alias <= 3) Then
          !  !C matrix of iatm differentiated w.r.t. position of jatm.
          !  dCxx(iatm,alias,:,:) = dC1k_da(iatm,jatm, lxx(iatm,:))
          !  dCyy(iatm,alias,:,:) = dC2k_da(iatm,jatm)
          !  dCzz(iatm,alias,:,:) = dC3k_da(iatm,jatm, lxx(iatm,:),lyy(iatm,:),dCxx(iatm,alias,:,:),dCyy(iatm,alias,:,:))
          !End If
        !End If
      End Do
    End If

    cFeat = cFeat + 1

!    Do jatm = 1, Size( non_alf_atms(ltg(iatm),:) ) 
     Do jatm = 1, model_list(ltg(iatm))%ptr%natms_model-3

        j = local_index(non_alf_atms(ltg(iatm),jatm), nlast,lsi,lsa)
    
        xdiff = xxx(j) - xxx(iatm)
        ydiff = yyy(j) - yyy(iatm)
        zdiff = zzz(j) - zzz(iatm)

        Call images_s(imcon,cell,xdiff,ydiff,zdiff)

    !======================================!
    !                 n-th feature, R^{An} !  
    !======================================!

        features(iatm,cFeat) = dsqrt( xdiff*xdiff + ydiff*ydiff + zdiff*zdiff  )
        
        ! calculate [[18]] 
        
        !df_da(cFeat,iatm,iatm,:) = dR_da( iatm, j, cFeat, iatm, wxx, wyy, wzz )
        !df_da(cFeat,j,iatm,   :) = dR_da( iatm, j, cFeat, j,    wxx, wyy, wzz )
    

        cFeat = cFeat + 1

    !======================================!
    !                           Theta^{An} !
    !======================================!

        zeta(3) = lzz(iatm,1)*(xdiff) + lzz(iatm,2)*(ydiff) + lzz(iatm,3)*(zdiff) ! [[23]]

        features(iatm,cFeat) = dacos( zeta(3) / features(iatm,cFeat-1) )
        
        ! calculate [[21]] 
        ! Check the manipulation of this cFEAT
!        df_da(cFeat,iatm,iatm,       :) = dTheta_da(iatm, j, cFeat, zeta(3), iatm,       wxx,wyy,wzz)
!        df_da(cFeat,j,iatm,          :) = dTheta_da(iatm, j, cFeat, zeta(3), j,          wxx,wyy,wzz) 
!        df_da(cFeat,alf(iatm,2),iatm,:) = dTheta_da(iatm, j, cFeat, zeta(3), alf(iatm,2),wxx,wyy,wzz)
!        df_da(cFeat,alf(iatm,3),iatm,:) = dTheta_da(iatm, j, cFeat, zeta(3), alf(iatm,3),wxx,wyy,wzz)

        cFeat = cFeat + 1

    !======================================!
    !                             phi^{An} !
    !======================================!

        zeta(2) = lyy(iatm,1)*(xdiff) + lyy(iatm,2)*(ydiff) + lyy(iatm,3)*(zdiff) ! [[23]]
        zeta(1) = lxx(iatm,1)*(xdiff) + lxx(iatm,2)*(ydiff) + lxx(iatm,3)*(zdiff) ! [[23]]
        features(iatm,cFeat) = datan2( zeta(2) , zeta(1) ) 
        
        ! calculate [[22]] 
        ! Check the manipulation of this cFEAT due to the theta/phi swap
        
        !df_da(cFeat,iatm,iatm,       :) = dPhi_da(iatm, j, zeta(1), zeta(2), cFeat, iatm,       wxx,wyy,wzz)
        !df_da(cFeat,j,iatm,          :) = dPhi_da(iatm, j, zeta(1), zeta(2), cFeat, j,          wxx,wyy,wzz)
        !df_da(cFeat,alf(iatm,2),iatm,:) = dPhi_da(iatm, j, zeta(1), zeta(2), cFeat, alf(iatm,2),wxx,wyy,wzz)
        !df_da(cFeat,alf(iatm,3),iatm,:) = dPhi_da(iatm, j, zeta(1), zeta(2), cFeat, alf(iatm,3),wxx,wyy,wzz)

        cFeat = cFeat + 1

    End Do

End Subroutine 
