Subroutine d_ene_trq_mpoles(vircpe_dt,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating the change in energy produced by
! an infinitesimal rotation of multipoles
!
! Reference : Sagui, Pedersen, Darden, J. Chem. Phys. 120, 73 (2004)
!             doi: 10.1063/1.1630791
!
! copyright - daresbury laboratory
! author    - h.a.boateng april 2015
! amended   - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module, Only : imcon,natms,ltg,fxx,fyy,fzz,xxx,yyy,zzz,cell
  Use mpoles_module, Only : mprotm,mptrqx,mptrqy,mptrqz

  Implicit None

  Real( Kind = wp ),                   Intent(   Out ) :: vircpe_dt
  Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress

  Integer           :: idi,j,iatm,jatm

  Real( Kind = wp ) :: rsq,rrr,p1(1:3),p2(1:3),u(1:3),v(1:3),w(1:3),magp1,magp2, &
                       p2u(1:3),p1perp2(1:3),p2perp1(1:3),magp1perp2,magp2perp1, &
                       tx,ty,tz,dedu,dedv,dedw,stx(1:3),sty(1:3),stz(1:3),       &
                       dux_dx,dux_dy,dux_dz,duy_dy,duy_dz,duz_dz,rmag,rmag3,     &
                       fx,fy,fz,fix,fiy,fiz,dedux,deduy,deduz,                   &
                       tmptx,tmpty,tmptz,xdf,ydf,zdf,p3(1:3),magp3,              &
                       strs1,strs2,strs3,strs4,strs5,strs6,strs7,strs8,strs9,tmp(6)

! initialise virial

  vircpe_dt=0.0_wp

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs4=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs7=0.0_wp
  strs8=0.0_wp
  strs9=0.0_wp

  Do iatm = 1, natms

! load forces

     fix=fxx(iatm)
     fiy=fyy(iatm)
     fiz=fzz(iatm)

! global identity of iatm

     idi=ltg(iatm)

     If (mprotm(iatm)%flag == 1) Then

! p1 and p2 define the local frame

        p1 = mprotm(iatm)%p1 !BCBS -> Redefine these as vectors used to form ALF?
        p2 = mprotm(iatm)%p2

        magp1 = Sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
        magp2 = Sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))

        p2u = p2/magp2

        !BCBS addition
        p3(1) = p1(2)*p2(3) - p1(3)*p2(2)
        p3(2) = p1(3)*p2(1) - p1(1)*p2(3)
        p3(3) = p1(1)*p2(2) - p1(2)*p2(1)

        magp3 = Sqrt(p3(1)*p3(1)+p3(2)*p3(2)+p3(3)*p3(3))

! standard basis for coordinate system

        u(1)  = mprotm(iatm)%mtrxa(1)
        u(2)  = mprotm(iatm)%mtrxa(4)
        u(3)  = mprotm(iatm)%mtrxa(7)

        v(1)  = mprotm(iatm)%mtrxa(2)
        v(2)  = mprotm(iatm)%mtrxa(5)
        v(3)  = mprotm(iatm)%mtrxa(8)

        w(1)  = mprotm(iatm)%mtrxa(3)
        w(2)  = mprotm(iatm)%mtrxa(6)
        w(3)  = mprotm(iatm)%mtrxa(9)

! change in energy (E) due to infinitesimal rotation (torque => \tau) dE_{\omega} = -\tau * d\omega
! we omit the negative here and introduce it in the final force computation, i.e., because force is
! the negative of the change in energy with respect to energy, we'll add the magnitude of the force
! computed instead of subtracting the magnitude

        tmptx = mptrqx(iatm)
        tmpty = mptrqy(iatm)
        tmptz = mptrqz(iatm)

        tx = tmptx*u(1) + tmpty*u(2) + tmptz*u(3)
        !ty = (tmptx*p2(1) + tmpty*p2(2) + tmptz*p2(3)) / magp2
        ty = tmptx*v(1) + tmpty*v(2) + tmptz*v(3)
        tz = tmptx*w(1) + tmpty*w(2) + tmptz*w(3)

! component of p2 perpendicular to p1

        p2perp1    = p2 - (u(1)*p2(1) + u(2)*p2(2) + u(3)*p2(3))*u
        magp2perp1 = Sqrt(p2perp1(1)*p2perp1(1)+p2perp1(2)*p2perp1(2)+p2perp1(3)*p2perp1(3))

! component of p1 perpendicular to p2

        p1perp2    = p1 - (p2u(1)*p1(1) + p2u(2)*p1(2) + p2u(3)*p1(3))*p2u
        magp1perp2 = Sqrt(p1perp2(1)*p1perp2(1)+p1perp2(2)*p1perp2(2)+p1perp2(3)*p1perp2(3))

! For p1

! dedu = 0.0 since movement du along the u axis does not rotate the frame
        !Rotation about w axis
        dedv   = tz/magp1
        !Rotation about p2u axis
        !dedw   =-(ty)/magp1perp2
        !Rotation about v axis.
        dedw = -ty/magp1

        stx(1) = dedv*v(1) + dedw*w(1)
        sty(1) = dedv*v(2) + dedw*w(2)
        stz(1) = dedv*v(3) + dedw*w(3)

! for p2

        !Rotation about u axis
        dedw   = (tx)/magp2perp1
        !dedw = tx/(p2(1)*v(1) + p2(2)*v(2) + p2(3)*v(3))

        !Rotation about w axis -> BCBS addition
        dedu   =-(tz)/(p2(1)*v(1) + p2(2)*v(2) + p2(3)*v(3))
        dedv   = (tz)/(p2(1)*u(1) + p2(2)*u(2) + p2(3)*u(3))
        !dedu = 0.0_wp
        !dedv = 0.0_wp

        !Rotation about v axis
        dedw = dedw - (ty)/(p2(1)*u(1) + p2(2)*u(2) + p2(3)*u(3))

        stx(2) = dedu*u(1) + dedv*v(1) + dedw*w(1)
        sty(2) = dedu*u(2) + dedv*v(2) + dedw*w(2)
        stz(2) = dedu*u(3) + dedv*v(3) + dedw*w(3)

! for p3 - BCBS
!        dedu = ty/magp3
!        dedv = -tx/magp3
!
!        stx(3) = dedu*u(1) + dedv*v(1)
!        sty(3) = dedu*u(2) + dedv*v(2)
!        stz(3) = dedu*u(3) + dedv*v(3)

! now compute forces and virial

        fx = 0.0_wp ; fy = 0.0_wp ; fz = 0.0_wp

        Do j = 1, 2

           jatm = mprotm(iatm)%mbnd(j)

           If (jatm > 0) Then
              xdf = xxx(jatm) - xxx(iatm)
              ydf = yyy(jatm) - yyy(iatm)
              zdf = zzz(jatm) - zzz(iatm)

              !Call images_s(imcon,cell,xdf,ydf,zdf)

              !Surely these vecs need to be p1 and p2.
              !If unit vecs then yes... otherwise no but then below derivs are
              !wrong.
              !If (j .EQ. 1) Then
              !  xdf = p1(1)
              !  ydf = p1(2)
              !  zdf = p1(3)
              !Else
              !  xdf = p2(1)
              !  ydf = p2(2)
              !  zdf = p2(3)
              !End If

              !Derivs below are for unit vecs - these vecs aren't...
              !Units don't add up... forces are coming out as N/m in the end...
              !Derivs are wrt xxx(jatm),yyy(jatm),zzz(jatm) i.e. global coords
              !of jatm - BCBS

              !rsq   = xdf*xdf + ydf*ydf + zdf*zdf
              !rrr   = sqrt(rsq)
              !rmag  = 1.0_wp/rrr
              !rmag3 = 1.0_wp/(rsq*rrr)

              !dux_dx = rmag - xdf * xdf * rmag3
              !dux_dy =      - xdf * ydf * rmag3       ! duy_dx = dux_dy
              !dux_dz =      - xdf * zdf * rmag3       ! duz_dx = dux_dz
              !duy_dy = rmag - ydf * ydf * rmag3
              !duy_dz =      - ydf * zdf * rmag3       ! duz_dy = duy_dz
              !duz_dz = rmag - zdf * zdf * rmag3

              dedux  = stx(j)
              deduy  = sty(j)
              deduz  = stz(j)

! Now to find the forces (derivatives of energy with respect to cartesian positions),
! i.e. fx=dedx, fy=dedy, fz=dedz

              !fx = dedux * dux_dx + deduy * dux_dy + deduz * dux_dz
              !fy = dedux * dux_dy + deduy * duy_dy + deduz * duy_dz
              !fz = dedux * dux_dz + deduy * duy_dz + deduz * duz_dz

              !If (j .EQ. 1) Then
              !  fx = dedux 
              !  fy = deduy
              !  fz = deduz

              !  dedux = stx(3)
              !  deduy = sty(3)
              !  deduz = stz(3)

              !  fx = fx - deduy*p2(3) + deduz*p2(2)
              !  fy = fy + dedux*p2(3) - deduz*p2(1)
              !  fz = fz - dedux*p2(2) + deduy*p2(1)
              !Else
              !  fx = dedux 
              !  fy = deduy  
              !  fz = deduz

              !  dedux = stx(3)
              !  deduy = sty(3)
              !  deduz = stz(3)

              !  fx = fx + deduy*p1(3) - deduz*p1(2)
              !  fy = fy - dedux*p1(3) + deduz*p1(1)
              !  fz = fz + dedux*p1(2) - deduy*p1(1)

              !End If

              fx = -1.0_wp*dedux
              fy = -1.0_wp*deduy
              fz = -1.0_wp*deduz

              fix = fix + fx
              fiy = fiy + fy
              fiz = fiz + fz
              tmp(iatm) = tmp(iatm) + fx

              If (jatm <= natms) Then

                 fxx(jatm)=fxx(jatm)-fx
                 tmp(jatm) = tmp(jatm) - fx
                 fyy(jatm)=fyy(jatm)-fy
                 fzz(jatm)=fzz(jatm)-fz

              End If

              If (jatm <= natms .or. idi < ltg(jatm)) Then

! calculate virial

                 vircpe_dt = vircpe_dt - (fx*xdf + fy*ydf + fz*zdf)

! calculate stress tensor

                 strs1 = strs1 + xdf*fx
                 strs2 = strs2 + xdf*fy
                 strs3 = strs3 + xdf*fz
                 strs4 = strs4 + ydf*fx
                 strs5 = strs5 + ydf*fy
                 strs6 = strs6 + ydf*fz
                 strs7 = strs7 + zdf*fx
                 strs8 = strs8 + zdf*fy
                 strs9 = strs9 + zdf*fz

              End If

           End If

        End Do

! load back forces

        fxx(iatm)=fix
        fyy(iatm)=fiy
        fzz(iatm)=fiz

     End If

! complete stress tensor (and symmetrise)

     stress(1) = stress(1) + strs1
     stress(2) = stress(2) + 0.5_wp * (strs2 + strs4)
     stress(3) = stress(3) + 0.5_wp * (strs3 + strs7)
     stress(4) = stress(4) + 0.5_wp * (strs2 + strs4)
     stress(5) = stress(5) + strs5
     stress(6) = stress(6) + 0.5_wp * (strs6 + strs8)
     stress(7) = stress(7) + 0.5_wp * (strs3 + strs7)
     stress(8) = stress(8) + 0.5_wp * (strs6 + strs8)
     stress(9) = stress(9) + strs9

  End Do

  Print *, 'D_ENE_TRQ'
  Do iatm =1,natms
    Print *, iatm, tmp(iatm)
  End Do

End Subroutine d_ene_trq_mpoles
