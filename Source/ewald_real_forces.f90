!Subroutine ewald_real_forces &
!           (iatm,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engcpe_rl,vircpe_rl,stress)

Subroutine ewald_real_forces &
           (iatm,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engcpe_rl,vircpe_rl,stress,drewd,rdrewd,erc,fer,dQ_da,fxc,fyc,fzc,mpi_counter,ffx,ffy,ffz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! in a periodic system using ewald's method
!
! Note: real space terms
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov april 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode
  Use setup_module
  Use config_module, Only : natms,ltg,list,chge,fxx,fyy,fzz,nlast,lsi,lsa, &
                          & xxx,yyy,zzz
  Use fflux_module,  Only : model_list,alf,non_alf_atms,max_atms_model, &
                          & atm_label_alias,mpi_list

  Implicit None

  Integer,                                  Intent( In    ) :: iatm,mpi_counter
  Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
  Real( Kind = wp ),                        Intent(   Out ) :: engcpe_rl,vircpe_rl
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

  Real( Kind = wp ),                         Intent( In    ) :: drewd,rdrewd
  Real( Kind = wp ), Dimension( 0:mxgele ),  Intent( In    ) :: erc,fer
  Real( Kind = wp ), Dimension( max_atms_model,nlast,1,3 ), Intent(In) :: dQ_da
  Real( Kind = wp ), Dimension(mpi_counter), Intent( InOut ) :: fxc,fyc,fzc
  Real( Kind = wp ), Dimension(natms),       Intent( InOut ) :: ffx,ffy,ffz

!  Logical,           Save :: newjob = .true.
!  Real( Kind = wp ), Save :: drewd,rdrewd

  Logical           :: Flag
  Integer           :: fail,m,idi,jatm,k,tempInt,z
  Real( Kind = wp ) :: chgea,chgprd,rrr,ppp,egamma,   &
                       fix,fiy,fiz,fx,fy,fz,          &
                       vk0,vk1,vk2,gk0,gk1,gk2,t1,t2, &
                       strs1,strs2,strs3,strs5,strs6,strs9
  !Real( Kind = wp ), Save :: tmpx(1:7000),tmpy(1:7000),tmpz(1:7000)

! Added for FFLUX force term - BCBS 01/21
  Real( Kind = wp )              :: Q,dQ_dx,dQ_dy,dQ_dz,poti,potj
  Integer                        :: local_index,j,ind

!  Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer

!  If (newjob) Then
!     newjob = .false.
!     fail=0
!     Allocate (erc(0:mxgele),fer(0:mxgele), Stat=fail)
!     If (fail > 0) Then
!        Write(nrite,'(/,1x,a,i0)') 'ewald_real_forces allocation failure, node: ', idnode
!        Call error(0)
!     End If
!
! interpolation interval
!
!     drewd = rcut/Real(mxgele-4,wp)
!
! reciprocal of interpolation interval
!
!     rdrewd = 1.0_wp/drewd
!
! generate error function complement tables for ewald sum
!
!     Call erfcgen(rcut,alpha,mxgele,erc,fer)
!  End If

! initialise potential energy and virial

  engcpe_rl=0.0_wp
  vircpe_rl=0.0_wp

  !If (iatm .EQ. 1) Then
  !  tmpx = 0.0_wp
  !  tmpy = 0.0_wp
  !  tmpz = 0.0_wp
  !End If

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

! start of primary loop for forces evaluation

     Do m=1,list(0,iatm)

! atomic index and charge

        jatm=list(m,iatm)
        chgprd=chge(jatm)

! interatomic distance

        rrr=rrt(m)

! interaction validity and truncation of potential

        If (Abs(chgprd) > zero_plus .and. rrr < rcut) Then

! charge product

           chgprd=chgprd*chgea

! calculate forces

           k   = Int(rrr*rdrewd)
           ppp = rrr*rdrewd - Real(k,wp)

! calculate forces using 3pt interpolation

           gk0 = fer(k) ; If (k == 0) gk0 = gk0*rrr
           gk1 = fer(k+1)
           gk2 = fer(k+2)

           t1 = gk0 + (gk1 - gk0)*ppp
           t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

           egamma = (t1 + (t2-t1)*ppp*0.5_wp)*chgprd

! calculate forces

           fx = egamma*xxt(m)
           fy = egamma*yyt(m)
           fz = egamma*zzt(m)

           fix=fix+fx
           fiy=fiy+fy
           fiz=fiz+fz

           If (jatm <= natms) Then

              fxx(jatm)=fxx(jatm)-fx
              fyy(jatm)=fyy(jatm)-fy
              fzz(jatm)=fzz(jatm)-fz

           End If

! calculate interaction energy using 3-point interpolation

           vk0 = erc(k)
           vk1 = erc(k+1)
           vk2 = erc(k+2)

           t1 = vk0 + (vk1 - vk0)*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           potj = (((t1+(t2-t1)*ppp*0.5_wp)*chgprd)/chge(iatm))
           poti = (((t1+(t2-t1)*ppp*0.5_wp)*chgprd)/chge(jatm))

           !FFLUX forces
           !Loop over ALF atoms
           Do j=1,3
                
             If (alf(iatm,j) .EQ. iatm) Then
               dQ_dx = dQ_da(j,iatm,1,1)
               dQ_dy = dQ_da(j,iatm,1,2)
               dQ_dz = dQ_da(j,iatm,1,3)

               ffx(iatm) = ffx(iatm) - (dQ_dx*potj)
               ffy(iatm) = ffy(iatm) - (dQ_dy*potj)
               ffz(iatm) = ffz(iatm) - (dQ_dz*potj)
             !If ensures forces are only computed for local (non halo) atoms
             Else If (alf(iatm,j) <= natms) Then
               dQ_dx = dQ_da(j,iatm,1,1)
               dQ_dy = dQ_da(j,iatm,1,2)
               dQ_dz = dQ_da(j,iatm,1,3)

               ffx(alf(iatm,j)) = ffx(alf(iatm,j)) - (dQ_dx*potj)
               ffy(alf(iatm,j)) = ffy(alf(iatm,j)) - (dQ_dy*potj)
               ffz(alf(iatm,j)) = ffz(alf(iatm,j)) - (dQ_dz*potj)
             Else
               Do k=1,mpi_counter
                 If (alf(iatm,j) .EQ. mpi_list(k)) Then
                   dQ_dx = dQ_da(j,iatm,1,1)
                   dQ_dy = dQ_da(j,iatm,1,2)
                   dQ_dz = dQ_da(j,iatm,1,3)
                                                                      
                   fxc(k) = fxc(k) - (dQ_dx*potj)
                   fyc(k) = fyc(k) - (dQ_dy*potj)
                   fzc(k) = fzc(k) - (dQ_dz*potj)
                 End If
               End Do
             End If

             !If ensures forces are only computed for local (non halo) atoms.
             !Also avoids double counting.
             If ((jatm <= natms) .AND. (alf(jatm,j) <= natms)) Then

               dQ_dx = dQ_da(j,jatm,1,1)
               dQ_dy = dQ_da(j,jatm,1,2)
               dQ_dz = dQ_da(j,jatm,1,3)

               ffx(alf(jatm,j)) = ffx(alf(jatm,j)) - (dQ_dx*poti)
               ffy(alf(jatm,j)) = ffy(alf(jatm,j)) - (dQ_dy*poti)
               ffz(alf(jatm,j)) = ffz(alf(jatm,j)) - (dQ_dz*poti)
             Else If (jatm <= natms) Then
               Do k=1,mpi_counter
                 If (alf(jatm,j) .EQ. mpi_list(k)) Then
                   dQ_dx = dQ_da(j,jatm,1,1)                                           
                   dQ_dy = dQ_da(j,jatm,1,2)
                   dQ_dz = dQ_da(j,jatm,1,3)
 
                   fxc(k) = fxc(k) - (dQ_dx*poti)
                   fyc(k) = fyc(k) - (dQ_dy*poti)
                   fzc(k) = fzc(k) - (dQ_dz*poti)
                 End If
               End Do
             End If

           End Do
           !Loop over non-ALF atoms of iatm
           Do j=1,model_list(ltg(iatm))%ptr%natms_model-3

             ind = local_index(non_alf_atms(idi,j), nlast,lsi,lsa)

             If (ind <= natms) Then
               dQ_dx = dQ_da(j+3,iatm,1,1)
               dQ_dy = dQ_da(j+3,iatm,1,2)
               dQ_dz = dQ_da(j+3,iatm,1,3)

               ffx(ind) = ffx(ind) - (dQ_dx*potj)
               ffy(ind) = ffy(ind) - (dQ_dy*potj)
               ffz(ind) = ffz(ind) - (dQ_dz*potj)
             Else
             !NEED TO CHECK THIS WORKS FOR NON ALF AS WELL.
               Do k=1,mpi_counter
                 If (ind .EQ. mpi_list(k)) Then
                   dQ_dx = dQ_da(j+3,iatm,1,1)
                   dQ_dy = dQ_da(j+3,iatm,1,2)
                   dQ_dz = dQ_da(j+3,iatm,1,3)
                                                      
                   fxc(k) = fxc(k) - (dQ_dx*potj)
                   fyc(k) = fyc(k) - (dQ_dy*potj)
                   fzc(k) = fzc(k) - (dQ_dz*potj)
                 End If
               End Do
             End If
 
           End Do
           !Loop over non-ALF atoms of jatm
           Do j=1,model_list(ltg(jatm))%ptr%natms_model-3

             ind = local_index(non_alf_atms(ltg(jatm),j), nlast,lsi,lsa)

             !If (jatm <= natms) Then
             If ((jatm <= natms) .AND. (ind <= natms)) Then
               dQ_dx = dQ_da(j+3,jatm,1,1)
               dQ_dy = dQ_da(j+3,jatm,1,2)
               dQ_dz = dQ_da(j+3,jatm,1,3)

               ffx(ind) = ffx(ind) - (dQ_dx*poti)
               ffy(ind) = ffy(ind) - (dQ_dy*poti)
               ffz(ind) = ffz(ind) - (dQ_dz*poti)
             Else If (jatm <= natms) Then
               Do k=1,mpi_counter
                 If (ind .EQ. mpi_list(k)) Then
                   dQ_dx = dQ_da(j+3,jatm,1,1)                          
                   dQ_dy = dQ_da(j+3,jatm,1,2)
                   dQ_dz = dQ_da(j+3,jatm,1,3)
                                                      
                   fxc(k) = fxc(k) - (dQ_dx*poti)
                   fyc(k) = fyc(k) - (dQ_dy*poti)
                   fzc(k) = fzc(k) - (dQ_dz*poti)
                 End If
               End Do
             End If
 
           End Do

           If (jatm <= natms .or. idi < ltg(jatm)) Then

             engcpe_rl = engcpe_rl + (t1 + (t2-t1)*ppp*0.5_wp)*chgprd

! calculate virial

             vircpe_rl = vircpe_rl - egamma*rrr**2

! calculate stress tensor

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

End Subroutine ewald_real_forces
