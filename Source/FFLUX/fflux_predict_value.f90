Subroutine fflux_predict_value( iatm, nf, Q_est, dQ_df )

  Use omp_lib
  Use kinds_f90
  Use setup_module,  Only : pi, twopi
  Use fflux_module,  Only : features,sign_j,atm_label_alias, &
                          & alf,model_list,maxFeatures
  Use config_module, Only : ltg,xxx,yyy,zzz,lsi,lsa,nlast
  Use comms_module, Only  : idnode

  Implicit None

! Inputs
  Integer, Intent( In ) :: iatm, nf

! Outputs
  Real( Kind = wp ) :: Q_est, dQ_df(maxFeatures)
 
! Local variables 
  Integer :: j, h, gInd
  Real( Kind = wp ) :: expo,fdiff, dQ_df_temp(maxFeatures)

! Theta values are defined to be 1/(2*l^2) where l is the lengthscale for any kernel.
! Keep this is mind when converting between theta/lengthscale values.


  Q_est = 0.0_wp
  dQ_df(:) = 0.0_wp

  gInd = ltg(iatm)
    Do j = 1, model_list(gInd)%ptr%nTrain

      expo = 0.0_wp

      Do h = 1, model_list(gInd)%ptr%nFeatures

        ! calculate exponent of Gaussian Kernel
        fdiff = model_list(gInd)%ptr%trained_value(nf,j,h) - features(iatm,h)

        !rbf-cyclic
        If (model_list(gInd)%ptr%kernels(h) < 3) Then

          ! Cyclic Feature Correction - Added by MJB 04/2020
          If ( mod(h, 3) .EQ. 0 ) Then
            fdiff = sign(1.0_wp, fdiff)*(Mod(Abs(fdiff)+pi,twopi) -pi)
          End If

          expo = expo + (model_list(gInd)%ptr%Theta(nf,h) * fdiff * fdiff)

          ! calculate the derivatives of dQ_df(h)
          dQ_df_temp(h) = ( model_list(gInd)%ptr%weights_krig(nf,j)*sign_j(fdiff) & 
                        & * (-2.0_wp) * model_list(gInd)%ptr%Theta(nf,h) * dabs(fdiff) )

        !Periodic kernel
        Else If (model_list(gInd)%ptr%kernels(h) == 3)Then

          expo = expo + 4.0_wp*model_list(gInd)%ptr%Theta(nf,h)*(Sin(fdiff/2.0_wp)**2)

          dQ_df_temp(h) = model_list(gInd)%ptr%weights_krig(nf,j)*(-4.0_wp)*model_list(gInd)%ptr%Theta(nf,h)*sign_j(fdiff)*Sin(dabs(fdiff)/2.0_wp)*Cos(dabs(fdiff)/2.0_wp)

        End If

      End Do

      expo = exp((-1.0_wp)*expo)

      ! Accumulate the predicted value

      Q_est = Q_est + (model_list(gInd)%ptr%weights_krig(nf,j) * expo)
   
      ! Multiply the derivatives of the features by the exponential

      dQ_df_temp(:) = dQ_df_temp(:) * expo

      ! Accumulate the dQ_df over nTrainedPoints

      dQ_df(:) = dQ_df(:) + dQ_df_temp(:)

    End Do

  Q_est = Q_est + model_list(gInd)%ptr%mu(nf)


End Subroutine
