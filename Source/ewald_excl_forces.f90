Subroutine ewald_excl_forces &
           (iatm,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engcpe_ex,vircpe_ex,stress,dQ_da,fxc,fyc,fzc,mpi_counter,ffx,ffy,ffz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! in a periodic system using ewald's method
!
! Note: exclusion correction terms
!       frozen pairs are ignored by default, they are not dealt with here
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module, Only : natms,ltg,list,chge,fxx,fyy,fzz,nlast,lsi,lsa
  Use comms_module,  Only : idnode
  Use fflux_module,  Only : alf,model_list,non_alf_atms, &
                          & max_atms_model,atm_label_alias,mpi_list

  Implicit None

  Integer,                                  Intent( In    ) :: iatm,mpi_counter
  Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
  Real( Kind = wp ),                        Intent(   Out ) :: engcpe_ex,vircpe_ex
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
  Real( Kind = wp ), Dimension(max_atms_model,nlast,1,3), Intent(In) :: dQ_da
  Real( Kind = wp ), Dimension(mpi_counter),Intent( InOut ) :: fxc,fyc,fzc
  Real( Kind = wp ), Dimension(natms),      Intent( InOut ) :: ffx,ffy,ffz

  Real( Kind = wp ), Parameter :: a1 =  0.254829592_wp
  Real( Kind = wp ), Parameter :: a2 = -0.284496736_wp
  Real( Kind = wp ), Parameter :: a3 =  1.421413741_wp
  Real( Kind = wp ), Parameter :: a4 = -1.453152027_wp
  Real( Kind = wp ), Parameter :: a5 =  1.061405429_wp
  Real( Kind = wp ), Parameter :: pp =  0.3275911_wp
  Real( Kind = wp ), Parameter :: rr3  = 1.0_wp/3.0_wp
  Real( Kind = wp ), Parameter :: r10  = 0.1_wp
  Real( Kind = wp ), Parameter :: r42  = 1.0_wp/42.0_wp
  Real( Kind = wp ), Parameter :: r216 = 1.0_wp/216.0_wp

  Integer           :: limit,idi,jatm,m
  Real( Kind = wp ) :: chgea,chgprd,rsq,rrr,alpr,alpr2, &
                       erfr,egamma,exp1,tt,             &
                       fix,fiy,fiz,fx,fy,fz,            &
                       strs1,strs2,strs3,strs5,strs6,strs9

!FFLUX variables for extra force terms - BCBS

  Real( Kind = wp ) :: dQ_dx,dQ_dy,dQ_dz,Q,poti,potj
  Integer           :: local_index,j,k,ind

! 'potential' accumulators
  poti = 0.0_wp
  potj = 0.0_wp

! initialise potential energy and virial

  engcpe_ex=0.0_wp
  vircpe_ex=0.0_wp

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

! global identity of iatm

  idi=ltg(iatm)

! ignore interaction if the charge is zero

  chgea = chge(iatm)

  If (Abs(chgea) > zero_plus) Then

     chgea = chgea*r4pie0/epsq

! load forces

     fix=fxx(iatm)
     fiy=fyy(iatm)
     fiz=fzz(iatm)

! Get list limit

     limit=list(-1,iatm)-list(0,iatm)

! start of primary loop for forces evaluation

     Do m=1,limit

! atomic index and charge
        jatm=list(list(0,iatm)+m,iatm)
        chgprd=chge(jatm)

! interatomic distance

        rrr=rrt(m)

! interaction validity and truncation of potential

        If (Abs(chgprd) > zero_plus .and. rrr < rcut) Then

! charge product

           chgprd=chgprd*chgea

! Squared distance

           rsq=rrr**2

! calculate forces

           alpr =rrr*alpha
           alpr2=alpr*alpr

! calculate error function and derivative

           If (alpr < 1.0e-2_wp) Then
! close particles (core-shell units) - small distances limit

              erfr=2.0_wp*chgprd*(alpha/sqrpi) * &
              (1.0_wp+alpr2*(-rr3+alpr2*(r10+alpr2*(-r42+alpr2*r216))))

              poti = erfr/chge(jatm)
              potj = erfr/chge(iatm)

              egamma=-4.0_wp*chgprd*(alpha**3/sqrpi) * &
              (rr3+alpr2*(-2.0_wp*r10+alpr2*(3.0_wp*r42-4.0_wp*alpr2*r216)))

           Else
! distant particles - traditional

              exp1=Exp(-(alpha*rrr)**2)
              tt  =1.0_wp/(1.0_wp+pp*alpha*rrr)

              erfr=chgprd * &
              (1.0_wp-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)/rrr

              poti = erfr/chge(jatm)
              potj = erfr/chge(iatm)

              egamma=-(erfr-2.0_wp*chgprd*(alpha/sqrpi)*exp1)/rsq

           End If

! calculate forces

           fx = egamma*xxt(m)
           fy = egamma*yyt(m)
           fz = egamma*zzt(m)

           !FFLUX force terms
           !Loop over ALF atoms.
           Do j=1,3

             If (alf(iatm,j) .EQ. iatm) Then
               dQ_dx = dQ_da(j,iatm,1,1)
               dQ_dy = dQ_da(j,iatm,1,2)
               dQ_dz = dQ_da(j,iatm,1,3)

               ffx(iatm) = ffx(iatm) + (dQ_dx*potj)
               ffy(iatm) = ffy(iatm) + (dQ_dy*potj)
               ffz(iatm) = ffz(iatm) + (dQ_dz*potj)
             Else If (alf(iatm,j) <= natms) Then
               dQ_dx = dQ_da(j,iatm,1,1)
               dQ_dy = dQ_da(j,iatm,1,2)
               dQ_dz = dQ_da(j,iatm,1,3)

               ffx(alf(iatm,j)) = ffx(alf(iatm,j)) + (dQ_dx*potj)
               ffy(alf(iatm,j)) = ffy(alf(iatm,j)) + (dQ_dy*potj)
               ffz(alf(iatm,j)) = ffz(alf(iatm,j)) + (dQ_dz*potj)
             Else
               Do k=1,mpi_counter
                 If (alf(iatm,j) .EQ. mpi_list(k)) Then
                   dQ_dx = dQ_da(j,iatm,1,1)
                   dQ_dy = dQ_da(j,iatm,1,2)
                   dQ_dz = dQ_da(j,iatm,1,3)
                                                                      
                   fxc(k) = fxc(k) + (dQ_dx*potj)
                   fyc(k) = fyc(k) + (dQ_dy*potj)
                   fzc(k) = fzc(k) + (dQ_dz*potj)
                 End If
               End Do
             End If
                                                                           
             If (jatm <= natms) Then 
               If (alf(jatm,j) .EQ. iatm) Then
                 dQ_dx = dQ_da(j,jatm,1,1)
                 dQ_dy = dQ_da(j,jatm,1,2)
                 dQ_dz = dQ_da(j,jatm,1,3)

                 ffx(iatm) = ffx(iatm) + (dQ_dx*poti)
                 ffy(iatm) = ffy(iatm) + (dQ_dy*poti)
                 ffz(iatm) = ffz(iatm) + (dQ_dz*poti)
               Else If (alf(jatm,j) <= natms) Then
                 dQ_dx = dQ_da(j,jatm,1,1)
                 dQ_dy = dQ_da(j,jatm,1,2)
                 dQ_dz = dQ_da(j,jatm,1,3)

                 ffx(alf(jatm,j)) = ffx(alf(jatm,j)) + (dQ_dx*poti)
                 ffy(alf(jatm,j)) = ffy(alf(jatm,j)) + (dQ_dy*poti)
                 ffz(alf(jatm,j)) = ffz(alf(jatm,j)) + (dQ_dz*poti)
               Else
                 Do k=1,mpi_counter
                   If (alf(jatm,j) .EQ. mpi_list(k)) Then
                     dQ_dx = dQ_da(j,jatm,1,1)
                     dQ_dy = dQ_da(j,jatm,1,2)
                     dQ_dz = dQ_da(j,jatm,1,3)
                                                                        
                     fxc(k) = fxc(k) + (dQ_dx*poti)
                     fyc(k) = fyc(k) + (dQ_dy*poti)
                     fzc(k) = fzc(k) + (dQ_dz*poti)
                   End If
                 End Do
               End If
             End If
                                                                           
           End Do
           !Loop over non-ALF atoms of iatm
           Do j=1,model_list(ltg(iatm))%ptr%natms_model-3
                                                                           
             ind = local_index(non_alf_atms(idi,j), nlast,lsi,lsa)

             If (ind <= natms) Then                           
               dQ_dx = dQ_da(j+3,iatm,1,1)
               dQ_dy = dQ_da(j+3,iatm,1,2)
               dQ_dz = dQ_da(j+3,iatm,1,3)
                                                                           
               ffx(ind) = ffx(ind) + (dQ_dx*potj)
               ffy(ind) = ffy(ind) + (dQ_dy*potj)
               ffz(ind) = ffz(ind) + (dQ_dz*potj)
             Else 
               Do k=1,mpi_counter
                 If (ind .EQ. mpi_list(k)) Then
                   dQ_dx = dQ_da(j+3,iatm,1,1)
                   dQ_dy = dQ_da(j+3,iatm,1,2)
                   dQ_dz = dQ_da(j+3,iatm,1,3)
                                                                               
                   fxc(k) = fxc(k) + (dQ_dx*potj)
                   fyc(k) = fyc(k) + (dQ_dy*potj)
                   fzc(k) = fzc(k) + (dQ_dz*potj)
                 End If
               End Do
             End If
                                                      
           End Do
           !Loop over non-ALF atoms of jatm
           Do j=1,model_list(ltg(jatm))%ptr%natms_model-3
                                                                           
             ind = local_index(non_alf_atms(ltg(jatm),j), nlast,lsi,lsa)

             If ((jatm <= natms) .AND. (ind <= natms)) Then
               dQ_dx = dQ_da(j+3,jatm,1,1)
               dQ_dy = dQ_da(j+3,jatm,1,2)
               dQ_dz = dQ_da(j+3,jatm,1,3)

               If (ind .EQ. iatm) Then
                 ffx(iatm) = ffx(iatm) + (dQ_dx*poti)
                 ffy(iatm) = ffy(iatm) + (dQ_dy*poti)
                 ffz(iatm) = ffz(iatm) + (dQ_dz*poti)
               Else
                 ffx(ind) = ffx(ind) + (dQ_dx*poti)
                 ffy(ind) = ffy(ind) + (dQ_dy*poti)
                 ffz(ind) = ffz(ind) + (dQ_dz*poti)
               End If
             Else If (jatm <= natms) Then
               Do k=1,mpi_counter
                 If (ind .EQ. mpi_list(k)) Then
                   dQ_dx = dQ_da(j+3,jatm,1,1)
                   dQ_dy = dQ_da(j+3,jatm,1,2)
                   dQ_dz = dQ_da(j+3,jatm,1,3)

                   fxc(k) = fxc(k) + (dQ_dx*poti) 
                   fyc(k) = fyc(k) + (dQ_dy*poti)
                   fzc(k) = fzc(k) + (dQ_dz*poti)
                 End If
               End Do
             End If
                                                                           
           End Do

           fix=fix+fx
           fiy=fiy+fy
           fiz=fiz+fz

           If (jatm <= natms) Then

              fxx(jatm)=fxx(jatm)-fx
              fyy(jatm)=fyy(jatm)-fy
              fzz(jatm)=fzz(jatm)-fz

           End If

           If (jatm <= natms .or. idi < ltg(jatm)) Then

! add potential energy and virial

              engcpe_ex = engcpe_ex - erfr
              vircpe_ex = vircpe_ex - egamma*rsq

! add stress tensor

              strs1 = strs1 + xxt(m)*fx
              strs2 = strs2 + xxt(m)*fy
              strs3 = strs3 + xxt(m)*fz
              strs5 = strs5 + yyt(m)*fy
              strs6 = strs6 + yyt(m)*fz
              strs9 = strs9 + zzt(m)*fz

           End If

        End If

     End Do

! load back forces

     fxx(iatm)=fix
     fyy(iatm)=fiy
     fzz(iatm)=fiz

! complete stress tensor

     stress(1) = stress(1) + strs1
     stress(2) = stress(2) + strs2
     stress(3) = stress(3) + strs3
     stress(4) = stress(4) + strs2
     stress(5) = stress(5) + strs5
     stress(6) = stress(6) + strs6
     stress(7) = stress(7) + strs3
     stress(8) = stress(8) + strs6
     stress(9) = stress(9) + strs9

  End If

End Subroutine ewald_excl_forces
