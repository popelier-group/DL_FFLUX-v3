! Put external multipoles here

Subroutine fflux_calc_mp_forces(iatm,jatm,rcut,epsq,engcpe,Localfxx,Localfyy,Localfzz,arrSize,dM_da,dq_da)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! in a periodic system using 1/r potential with no truncation or damping
!
! copyright - daresbury laboratory
! author    - t.forester february 1993
! amended   - i.t.todorov november 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use omp_lib
  Use config_module, Only : natms,ltg,list,chge,xxx,yyy,zzz,atmnam,fxx!,fyy,fzz
  Use setup_module 
  Use mpoles_module, Only : mplgfr
  Use fflux_module,  Only : conv_B2A,conv_Ha_2_10Jmol,maxFeatures,conv_MM,model_list, &
                          & df_da,alf,non_alf_atms,poleOrder,ewaldFlag,rt10, &
                          & max_atms_model,atm_label_alias,model_list

  Implicit None

  Integer,                                  Intent( In    ) :: iatm,jatm,arrSize
  Real( Kind = wp ),                        Intent( In    ) :: rcut,epsq
  Real( Kind = wp ),                        Intent( In    ) :: dM_da(natms,max_atms_model,arrSize,3)
  Real( Kind = wp ),                        Intent( InOut ) :: engcpe
  Real( Kind = wp ),                        Intent( InOut ) :: Localfxx(Size(fxx)),Localfyy(Size(fxx)),Localfzz(Size(fxx)),dq_da(0:15,natms,3)

  Integer                                                   :: nPred_A,nPred_B,i,j,k,L,limit,Ind_A,Ind_B

  Real( Kind = wp )                                         :: coul, QA, QB, TAB, dTAB_dq(0:15), atom_coul(natms),dQ_dx,dQ_dy,dQ_dz
  !Variables to be passed into calc_q_vector for use in explicit tensor, moved here by BCBS 01/2019 for purposes of OpenMP parallelisation.
  Real( Kind = wp )                                         :: R_, rax, ray, raz, & 
                                                                 & rbx, rby, rbz, &
                                                                 & Cxx, Cyx, Czx, &
                                                                 & Cxy, Cyy, Czy, &
                                                                 & Cxz, Cyz, Czz, &
                                                                 & Q_prod,QAT_prod,QBT_prod
  Real( Kind = wp ), Allocatable                            :: dTAB_da(:)

! Determine Multipole Moments on atoms
! If ( (nPredPerAtm .GT. 1) .AND. (iatm .NE. jatm) ) Then

   Allocate( dTAB_da( 3 ) )
    
!  ER: Before using higher rank monopoles, make sure that the Kriging.txt fil you are using is set to read the Kriging files need.

  ! MM loop starts here, replacing nPred

! ! nPred is a number defning the kind of multipole you are using, for example nPred=02 is the monopole Q00, nPred=05 is the dipole Q11s (l=1, m=-1), nPred=09 is the quadrupole Q22c (l=2, m=2). You can look into the details fflux_mm_explicittensor.f90 and fflux_predict_value.f90 for a deeper understanding of nPred.

!  To change the rank of the multipole moments you want to consider you just have to change the inequality in the if statement. Of course if you want to save some cpu time, you can change the range of nPred you are interesting in. For example L<=4 doesn't need hexadecapole.

!  Monopole     : nPred=2
!  Dipole       : nPred=3,5
!  Quadrupole   : nPred=6,10
!  Octopole     : nPred=11,17
!  Hexadecapole : nPred=18,26

 If (poleOrder == 1) Then
   limit = 2
 Else If (poleOrder == 2) Then
   limit = 5
 Else If (poleOrder == 3) Then
   limit = 10
 Else If (poleOrder == 4) Then
   limit = 17
 Else If (poleOrder == 5) Then
   limit = 26
 Else 
   Write(nffluxdt, '(a)') "Invalid option for L in CONTROL file"
 End If

 Do nPred_A = 2,limit

   If (mximpl .EQ. 0) Then
     QA = chge(iatm)   
   Else
     QA = mplgfr(nPred_A-1,iatm)
   End If
   Ind_A = atm_label_alias(ltg(iatm))

   Call calc_q_vector( iatm,jatm,xxx,yyy,zzz,R_,rax,ray,raz,rbx,rby,rbz,Cxx,Cyx,Czx,Cxy,Cyy,Czy,Cxz,Cyz,Czz )

   Do nPred_B = 2,limit

     L = conv_MM(npred_A-1,1)+conv_MM(npred_B-1,1)+1
      If (L .LE. poleOrder) Then

        Call ExplicitMultipoleTensor( nPred_A - 1,nPred_B - 1,TAB,dTAB_dq,R_,rax,ray,raz,rbx,rby,rbz,Cxx,Cyx,Czx,Cxy,Cyy,Czy,Cxz,Cyz,Czz )

        If (mximpl .EQ. 0) Then
          QB = chge(jatm)
        Else
          QB = mplgfr(nPred_B-1,jatm)
        End If
        Ind_B = atm_label_alias(ltg(jatm)) 

        Q_prod = QA*QB
        QAT_prod = QA*TAB
        QBT_prod = QB*TAB

        !Loop over ALF of iatm
        Do j=1,3

          i = alf(iatm,j)

          Call calc_q_derivs(iatm,jatm,i,xxx,yyy,zzz,dq_da)
         
          dTAB_da(:) = 0.0_wp

          Do k = 0, 15
            dTAB_da(1) = dTAB_da(1) + ( dTAB_dq(k) * dq_da(k,i,1) ) 
            dTAB_da(2) = dTAB_da(2) + ( dTAB_dq(k) * dq_da(k,i,2) ) 
            dTAB_da(3) = dTAB_da(3) + ( dTAB_dq(k) * dq_da(k,i,3) ) 
          End Do

          Localfxx(i) = Localfxx(i) - ( (dM_da(i,Ind_A,nPred_A-1,1) * QBT_prod) + (Q_prod * dTAB_da(1)) )
          Localfyy(i) = Localfyy(i) - ( (dM_da(i,Ind_A,nPred_A-1,2) * QBT_prod) + (Q_prod * dTAB_da(2)) )
          Localfzz(i) = Localfzz(i) - ( (dM_da(i,Ind_A,nPred_A-1,3) * QBT_prod) + (Q_prod * dTAB_da(3)) )

        End Do

        !Loop over non ALF of iatm
        Do j=1,model_list(iatm)%ptr%natms_model-3
                                                                                                          
          i = non_alf_atms(iatm,j)
                                                                                                          
          Call calc_q_derivs(iatm,jatm,i,xxx,yyy,zzz,dq_da)
                                                                         
          dTAB_da(:) = 0.0_wp
                                                                                                          
          Do k = 0, 15
            dTAB_da(1) = dTAB_da(1) + ( dTAB_dq(k) * dq_da(k,i,1) ) 
            dTAB_da(2) = dTAB_da(2) + ( dTAB_dq(k) * dq_da(k,i,2) ) 
            dTAB_da(3) = dTAB_da(3) + ( dTAB_dq(k) * dq_da(k,i,3) ) 
          End Do
                                                                                                          
          Localfxx(i) = Localfxx(i) - ( (dM_da(i,Ind_A,nPred_A-1,1) * QBT_prod))! + (Q_prod * dTAB_da(1)) )
          Localfyy(i) = Localfyy(i) - ( (dM_da(i,Ind_A,nPred_A-1,2) * QBT_prod))! + (Q_prod * dTAB_da(2)) )
          Localfzz(i) = Localfzz(i) - ( (dM_da(i,Ind_A,nPred_A-1,3) * QBT_prod))! + (Q_prod * dTAB_da(3)) )

        End Do

        !Loop over ALF of jatm
        Do j=1,3

          i = alf(jatm,j)

          Call calc_q_derivs(iatm,jatm,i,xxx,yyy,zzz,dq_da)

          dTAB_da(:) = 0.0_wp
                                                                         
          Do k = 0, 15
            dTAB_da(1) = dTAB_da(1) + ( dTAB_dq(k) * dq_da(k,i,1) ) 
            dTAB_da(2) = dTAB_da(2) + ( dTAB_dq(k) * dq_da(k,i,2) ) 
            dTAB_da(3) = dTAB_da(3) + ( dTAB_dq(k) * dq_da(k,i,3) ) 

          End Do

          Localfxx(i) = Localfxx(i) - ( (Q_prod * dTAB_da(1)) + (QAT_prod * dM_da(i,Ind_B,nPred_B-1,1)) )
          Localfyy(i) = Localfyy(i) - ( (Q_prod * dTAB_da(2)) + (QAT_prod * dM_da(i,Ind_B,nPred_B-1,2)) )
          Localfzz(i) = Localfzz(i) - ( (Q_prod * dTAB_da(3)) + (QAT_prod * dM_da(i,Ind_B,nPred_B-1,3)) )

        End Do

        !Loop over non ALF of jatm
        Do j=1,model_list(jatm)%ptr%natms_model-3
                                                                                                          
          i = non_alf_atms(jatm,j)
 
          Call calc_q_derivs(iatm,jatm,i,xxx,yyy,zzz,dq_da)
 
          dTAB_da(:) = 0.0_wp
 
          Do k = 0, 15
            dTAB_da(1) = dTAB_da(1) + ( dTAB_dq(k) * dq_da(k,i,1) ) 
            dTAB_da(2) = dTAB_da(2) + ( dTAB_dq(k) * dq_da(k,i,2) ) 
            dTAB_da(3) = dTAB_da(3) + ( dTAB_dq(k) * dq_da(k,i,3) ) 
          End Do

          !Localfxx(i) = Localfxx(i) - ( (Q_prod * dTAB_da(1)) + (QAT_prod * dM_da(i,Ind_B,nPred_B-1,1)) )
          !Localfyy(i) = Localfyy(i) - ( (Q_prod * dTAB_da(2)) + (QAT_prod * dM_da(i,Ind_B,nPred_B-1,2)) )
          !Localfzz(i) = Localfzz(i) - ( (Q_prod * dTAB_da(3)) + (QAT_prod * dM_da(i,Ind_B,nPred_B-1,3)) )

          Localfxx(i) = Localfxx(i) - ( (QAT_prod * dM_da(i,Ind_B,nPred_B-1,1)) )
          Localfyy(i) = Localfyy(i) - ( (QAT_prod * dM_da(i,Ind_B,nPred_B-1,2)) )
          Localfzz(i) = Localfzz(i) - ( (QAT_prod * dM_da(i,Ind_B,nPred_B-1,3)) )

        End Do

        coul = QA * TAB * QB
        atom_coul(iatm) = atom_coul(iatm) + ( coul )

!       calculate potential energy
        engcpe = engcpe + coul

      End If
     !End If
    End Do
  End Do

  Deallocate( dTAB_da )
 
! End If

End Subroutine 

