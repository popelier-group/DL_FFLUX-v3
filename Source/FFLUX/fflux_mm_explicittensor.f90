!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          !
! Subroutine to calculate the analytical form of the orientation tensor    !
! as well as it first and second derivatives with respect to Qi            !
!                                                                          !
! Author: O. Lamarche                                                      !
! Date:   06/04/04                                                         !
!                                                                          !
! Update: J CR Thacker                                                     !
! Subroutine now has explicit statements for all dipole moments.           !
! Wherein the dipole moment currently evaluated is stored in vars:         !
! rAlpha, rBeta, cxA, cxB, ..., etc.                                       !
!                                                                          !
! Qi code:                                                                 !
!                                                                          !
! 0:   R.R                                                                 !
! 1:   Xa.R,  2:   Ya.R,  3:   Za.R                                        !
! 4:   Xb.R,  5:   Yb.R,  6:   Zb.R                                        !
! 7:   Xa.Xb, 8:   Ya.Xb, 9:   Za.Xb                                       !
! 10:  Xa.Yb, 11:  Ya.Yb, 12:  Za.Yb                                       !
! 13:  Xa.Zb, 14:  Ya.Zb, 15:  Za.Zb                                       !
!                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitMultipoleTensor( i_idx,j_idx,T0,T1,R_,rax,ray,raz,rbx,rby,rbz,Cxx,Cyx,Czx,Cxy,Cyy,Czy,Cxz,Cyz,Czz )

   Use omp_lib
   Use kinds_f90
   Use setup_module, Only : nffluxdt
!   use fflux_module, Only : conv_MM, R_, &
!                          & rax,ray,raz, &
!                          & rbx,rby,rbz, &
!                          & Cxx,Cyx,Czx, &
!                          & Cxy,Cyy,Czy, &
!                          & Cxz,Cyz,Czz, &
   Use fflux_module, Only : conv_MM, &
                          & rt3, rt5, rt6, & 
                          & rt10, rt15, rt35, rt70

   implicit none

   integer, intent(in) :: i_idx, j_idx
   Real( Kind = wp ), Intent( In ) :: R_, rax, ray, raz, &
                                        & rbx, rby, rbz, &
                                        & Cxx, Cyx, Czx, &
                                        & Cxy, Cyy, Czy, &
                                        & Cxz, Cyz, Czz
   real(Kind = wp)      :: T0, T1(0:15), &
                        & rAlpha, rBeta, cAB, &
                        & cxA, cyA, czA, &
                        & cxB, cyB, czB, &
                        & R1, R2, R3, R4, &
                        & R5, R6, R7
   integer             :: la, ka, lb, kb, &
                          ia, ak, ib, bk, &
                        & oa, ob
   character(2)        :: dipind
   character(4)        :: idx

   la = conv_MM(i_idx,1)
   ka = conv_MM(i_idx,2)
   lb = conv_MM(j_idx,1)
   kb = conv_MM(j_idx,2)

   ! Intialise the local arrays for the tensor and the first derivatives.
   T0 = 0.0D0
   T1 = 0.0D0

   ! Explicit square root values needed in this subroutine
    rt3 = 0.17320508075689D+01
    rt5 = 0.22360679774998D+01
    rt6 = 0.24494897427832D+01
   rt10 = 0.31622776601684D+01
   rt15 = 0.38729833462074D+01
   rt35 = 0.59160797830996D+01
   rt70 = 0.83666002653408D+01

   ! Calculate:
   ! R^-1 ... R^-7
   R1 = 1.0D0 / (R_)
   R2 = R1 * R1
   R3 = R2 * R1
   R4 = R3 * R1
   R5 = R4 * R1
   R6 = R5 * R1
   R7 = R6 * R1


   ! Initialise the local values { ia, ak, ib, bk }
   ia = la
   ak = ka
   ib = lb
   bk = kb

   ! If one of the multipole moments is a dipole (l = 1),
   ! it is important to convert the k in such a way that:
   !    T{10,lk} or
   !    T{lk,10} or
   !    T{10,10}
   ! are always called but the correct value of 'k' for the dipole
   ! moment is used.

   if ( ia .EQ. 1 ) then ! if (la == dipole )

      ak = 0 ! ensure that labels {l,k = 1,0} for dipole on atom 1.

      ! This is used to ensure the correct indexing of the T1 array for dipoles.
      if     ( ka .EQ. 0 ) then ! if (k = z) !ORIGNAL
         oa = 2
         rAlpha=raz
         !cxA=Cxz; cyA=Cyz; czA=Czz
         cxA=Czx; cyA=Czy; czA=Czz
      elseif ( ka .EQ. 1 ) then ! if (k = x) !ORIGNAL
         oa  = 0
         rAlpha=rax
         !cxA=Cxx; cyA=Cyx; czA=Czx
         cxA=Cxx; cyA=Cxy; czA=Cxz
      elseif ( ka .EQ. 2 ) then ! if (k = y) !ORIGNAL
         oa = 1
         rAlpha=ray
         !cxA=Cxy; cyA=Cyy; czA=Czy
         cxA=Cyx; cyA=Cyy; czA=Cyz
      end if

      !BCBS replaced case select for optimisation (writing was very slow).
      if ( ib .EQ. 1 ) then ! if both the multipole moments are dipoles
          If ((ka .EQ. 0) .AND. (kb .EQ. 0)) Then
              cAB=Czz
          Else If ((ka .EQ. 0) .AND. (kb .EQ. 1)) Then  
              cAB=Czx
          Else If ((ka .EQ. 0) .AND. (kb .EQ. 2)) Then  
              cAB=Czy
          Else If ((ka .EQ. 1) .AND. (kb .EQ. 0)) Then  
              cAB=Cxz
          Else If ((ka .EQ. 1) .AND. (kb .EQ. 1)) Then  
              cAB=Cxx
          Else If ((ka .EQ. 1) .AND. (kb .EQ. 2)) Then  
              cAB=Cxy
          Else If ((ka .EQ. 2) .AND. (kb .EQ. 0)) Then  
              cAB=Cyz
          Else If ((ka .EQ. 2) .AND. (kb .EQ. 1)) Then  
              cAB=Cyx
          Else If ((ka .EQ. 2) .AND. (kb .EQ. 2)) Then  
              cAB=Cyy
          End If
      end if 
   end if

   ! same as previous "if" but for dipole on atom 2
   !BCBS bug fix - indices on C terms weren't reversed.
   if ( ib .EQ. 1 ) then
      bk = 0
      if     ( kb .EQ. 0 ) then
         ob = 2
         rBeta=rbz
         cxB=Cxz; cyB=Cyz; czB=Czz
         !cxB=Czx; cyB=Czy; czB=Czz
      elseif ( kb .EQ. 1 ) then
         ob = 0
         rBeta=rbx
         cxB=Cxx; cyB=Cyx; czB=Czx
         !cxB=Cxx; cyB=Cxy; czB=Cxz
      elseif ( kb .EQ. 2 ) then
         ob = 1
         rBeta=rby
         cxB=Cxy; cyB=Cyy; czB=Czy
         !cxB=Cyx; cyB=Cyy; czB=Cyz
      end if
   end if

   ! write case selector - dipoles always represented '10',
   ! regardless of orientation.
   !write(idx,"(4I1)") ia, ak, ib, bk

!
!  L = 1
!  Charge-charge       00, 00
!
        If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R1
          T1(0)     = 1.0_wp*R3 * (-0.5_wp)
!
!  L = 2
!  Dipole-charge       10, 00
!
        Else If ((ia .EQ. 1) .AND. (ak .EQ. 0) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R2 * rAlpha
          T1(0)     = 1.0_wp*R4 * (-1.5_wp) * rAlpha
          T1(1+oa)  = 1.0_wp*R3
!
!  Charge-dipole       00, 10
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 1) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R2 * (-1.0_wp)*rBeta
          T1(0)     = 1.0_wp*R4 * 1.5_wp * rBeta
          T1(4+ob)  = 1.0_wp*R3 * (-1.0_wp)
!
!  L = 3
!  Quadrupole-charge    20, 00
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 0) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R3 * ((3.0D0/2.0D0)*(raz**(2.0D0))-(1.0D0/2.0D0))
          T1(0)     = 1.0_wp*R5 * ((-15.0D0/4.0D0)*(raz**(2.0D0))+(3.0D0/4.0D0))
          T1(3)     = 1.0_wp*R4 * (3.0D0*raz)
!
!                      21c, 00
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 1) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R3 * (rt3*rax*raz)
          T1(0)     = 1.0_wp*R5 * (rt3*(-5.0D0/2.0D0)*rax*raz)
          T1(1)     = 1.0_wp*R4 * (rt3*raz)
          T1(3)     = 1.0_wp*R4 * (rt3*rax)
!
!                      21s, 00
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 2) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R3 * (rt3*ray*raz)
          T1(0)     = 1.0_wp*R5 * (rt3*(-5.0D0/2.0D0)*ray*raz)
          T1(2)     = 1.0_wp*R4 * (rt3*raz)
          T1(3)     = 1.0_wp*R4 * (rt3*ray)
!
!                      22c, 00
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 3) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R3 * (rt3*(1.0D0/2.0D0)*(rax*rax-ray*ray))
          T1(0)     = 1.0_wp*R5 * (rt3*(-5.0D0/4.0D0)*(rax*rax-ray*ray))
          T1(1)     = 1.0_wp*R4 * (rt3*rax)
          T1(2)     = 1.0_wp*R4 * (-1.0D0*rt3*ray)
!
!                      22s, 00
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 4) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R3 * (rt3*rax*ray)
          T1(0)     = 1.0_wp*R5 * (rt3*(-5.0D0/2.0D0)*rax*ray)
          T1(1)     = 1.0_wp*R4 * (rt3*ray)
          T1(2)     = 1.0_wp*R4 * (rt3*rax)
!
!  Charge-quadrupole   00, 20
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 0)) Then
          !T0        = R3 * ((3.0D0/2.0D0)*(rbz**(2.0D0))-(1.0D0/2.0D0))
          !T1(0)     = R5 * ((-15.0D0/4.0D0)*(rbz**(2.0D0))+(3.0D0/4.0D0))
          !T1(6)     = R4 * (3.0D0*rbz)
          T0        = 1.0_wp*R3 * ((3.0D0/2.0D0)*(rbz**(2.0D0))-(1.0D0/2.0D0))
          T1(0)     = 1.0_wp*R5 * ((-15.0D0/4.0D0)*(rbz**(2.0D0))+(3.0D0/4.0D0))
          T1(6)     = 1.0_wp*R4 * (3.0D0*rbz)
!
!                      00, 21c
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 1)) Then
          !T0        = R3 * (rt3*rbx*rbz)
          !T1(0)     = R5 * (rt3*(-5.0D0/2.0D0)*rbx*rbz)
          !T1(4)     = R4 * (rt3*rbz)
          !T1(6)     = R4 * (rt3*rbx)
          T0        = 1.0_wp*R3 * (rt3*rbx*rbz)
          T1(0)     = 1.0_wp*R5 * (rt3*(-5.0D0/2.0D0)*rbx*rbz)
          T1(4)     = 1.0_wp*R4 * (rt3*rbz) * (1.0_wp)
          T1(6)     = 1.0_wp*R4 * (rt3*rbx) * (1.0_wp)
!
!                      00, 21s
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 2)) Then
          !T0        = R3 * (rt3*rby*rbz)
          !T1(0)     = R5 * (rt3*(-5.0D0/2.0D0)*rby*rbz)
          !T1(5)     = R4 * (rt3*rbz)
          !T1(6)     = R4 * (rt3*rby)
          T0        = 1.0_wp*R3 * (rt3*rby*rbz)
          T1(0)     = 1.0_wp*R5 * (rt3*(-5.0D0/2.0D0)*rby*rbz)
          T1(5)     = 1.0_wp*R4 * (rt3*rbz) * (1.0_wp)
          T1(6)     = 1.0_wp*R4 * (rt3*rby) * (1.0_wp)
!
!                      00, 22c
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 3)) Then
          !T0        = R3 * (rt3*(1.0D0/2.0D0)*(rbx*rbx-rby*rby))
          !T1(0)     = R5 * (rt3*(-5.0D0/4.0D0)*(rbx*rbx-rby*rby))
          !T1(4)     = R4 * (rt3*rbx)
          !T1(5)     = R4 * (-1.0D0*rt3*rby)
          T0        = 1.0_wp*R3 * (rt3*(1.0D0/2.0D0)*(rbx*rbx-rby*rby))
          T1(0)     = 1.0_wp*R5 * (rt3*(-5.0D0/4.0D0)*(rbx*rbx-rby*rby))
          T1(4)     = 1.0_wp*R4 * (rt3*rbx)*(1.0_wp)
          T1(5)     = 1.0_wp*R4 * (rt3*rby)*(-1.0_wp)
!
!                      00, 22s
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 4)) Then
          !T0        = R3 * (rt3*rbx*rby)
          !T1(0)     = R5 * (rt3*(-5.0D0/2.0D0)*rbx*rby)
          !T1(4)     = R4 * (rt3*rby)
          !T1(5)     = R4 * (rt3*rbx)
          T0        = 1.0_wp*R3 * (rt3*rbx*rby)
          T1(0)     = 1.0_wp*R5 * (rt3*(-5.0D0/2.0D0)*rbx*rby)
          T1(4)     = 1.0_wp*R4 * (rt3*rby) * (1.0_wp)
          T1(5)     = 1.0_wp*R4 * (rt3*rbx) * (1.0_wp)
!
!  Dipole-dipole       10, 10
!
        Else If ((ia .EQ. 1) .AND. (ak .EQ. 0) .AND. (ib .EQ. 1) .AND. (bk .EQ. 0)) Then
          T0              = 1.0_wp*R3 * ( (-3.0D0 * rAlpha * rBeta) + cAB )
          T1(0)           = 1.0_wp*R5 * ( ((7.5D0) * rAlpha * rBeta) - (1.5D0 * cAB) )
          T1(1+oa)        = 1.0_wp*R4 * ( -3.0D0 * rBeta )
          T1(4+ob)        = 1.0_wp*R4 * ( -3.0D0 * rAlpha )
          T1(7+(ob*3)+oa) = 1.0_wp*R3 
!
!  L = 4
!  Octopole-charge     30, 00
!
        Else If ((ia .EQ. 3) .AND. (ak .EQ. 0) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R4 * ((5.0D0/2.0D0)*(raz**(3.0D0))-(3.0D0/2.0D0)*raz)
          T1(0)     = 1.0_wp*R6 * ((-35.0D0/4.0D0)*(raz**(3.0D0))+(15.0D0/4.0D0)*raz)
          T1(3)     = 1.0_wp*R5 * ((15.0D0/2.0D0)*(raz**(2.0D0))-(3.0D0/2.0D0))
!
!                      31c, 00
!
        Else If ((ia .EQ. 3) .AND. (ak .EQ. 1) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R4 * (rt6*((5.0D0/4.0D0)*rax*(raz**(2.0D0))-(1.0D0/4.0D0)*rax))
          T1(0)     = 1.0_wp*R6 * (rt6*((-35.0D0/8.0D0)*rax*(raz**(2.0D0))+(5.0D0/8.0D0)*rax))
          T1(1)     = 1.0_wp*R5 * (rt6*((5.0D0/4.0D0)*(raz**(2.0D0))-(1.0D0/4.0D0)))
          T1(3)     = 1.0_wp*R5 * (rt6*((5.0D0/2.0D0)*rax*raz))
! 
!                      31s, 00
!
        Else If ((ia .EQ. 3) .AND. (ak .EQ. 2) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R4 * (rt6*((5.0D0/4.0D0)*ray*(raz**(2.0D0))-(1.0D0/4.0D0)*ray))
          T1(0)     = 1.0_wp*R6 * (rt6*((-35.0D0/8.0D0)*ray*(raz**(2.0D0))+(5.0D0/8.0D0)*ray))
          T1(2)     = 1.0_wp*R5 * (rt6*((5.0D0/4.0D0)*(raz**(2.0D0))-(1.0D0/4.0D0)))
          T1(3)     = 1.0_wp*R5 * (rt6*((5.0D0/2.0D0)*ray*raz))
!    
!                      32c, 00
!
        Else If ((ia .EQ. 3) .AND. (ak .EQ. 3) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R4 * (rt15*((1.0D0/2.0D0)*(rax**(2.0D0))*raz-(1.0D0/2.0D0)*(ray**(2.0D0))*raz))
          T1(0)     = 1.0_wp*R6 * (rt15*((-7.0D0/4.0D0)*(rax**(2.0D0))*raz+(7.0D0/4.0D0)*(ray**(2.0D0))*raz))
          T1(1)     = 1.0_wp*R5 * (rt15*rax*raz)
          T1(2)     = 1.0_wp*R5 * (rt15*(-1.0D0)*ray*raz)
          T1(3)     = 1.0_wp*R5 * (rt15*((1.0D0/2.0D0)*(rax**(2.0D0))-(1.0D0/2.0D0)*(ray**(2.0D0))))
!
!                      32s, 00
!
        Else If ((ia .EQ. 3) .AND. (ak .EQ. 4) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R4 * (rt15*rax*ray*raz)
          T1(0)     = 1.0_wp*R6 * (rt15*(-7.0D0/2.0D0)*rax*ray*raz)
          T1(1)     = 1.0_wp*R5 * (rt15*ray*raz)
          T1(2)     = 1.0_wp*R5 * (rt15*rax*raz)
          T1(3)     = 1.0_wp*R5 * (rt15*rax*ray)
!
!                      33c, 00
!
        Else If ((ia .EQ. 3) .AND. (ak .EQ. 5) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R4 * (rt10*((1.0D0/4.0D0)*(rax**(3.0D0))-(3.0D0/4.0D0)*rax*(ray**(2.0D0))))
          T1(0)     = 1.0_wp*R6 * (rt10*((-7.0D0/8.0D0)*(rax**(3.0D0))+(21.0D0/8.0D0)*rax*(ray**(2.0D0))))
          T1(1)     = 1.0_wp*R5 * (rt10*((3.0D0/4.0D0)*(rax**(2.0D0))-(3.0D0/4.0D0)*(ray**(2.0D0))))
          T1(2)     = 1.0_wp*R5 * (rt10*(-3.0D0/2.0D0)*rax*ray)
!
!                      33s, 00
!
        Else If ((ia .EQ. 3) .AND. (ak .EQ. 6) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = 1.0_wp*R4 * (rt10*((3.0D0/4.0D0)*ray*(rax**(2.0D0))-(1.0D0/4.0D0)*(ray**(3.0D0))))
          T1(0)     = 1.0_wp*R6 * (rt10*((-21.0D0/8.0D0)*ray*(rax**(2.0D0))+(7.0D0/8.0D0)*(ray**(3.0D0))))
          T1(1)     = 1.0_wp*R5 * (rt10*(3.0D0/2.0D0)*rax*ray)
          T1(2)     = 1.0_wp*R5 * (rt10*((3.0D0/4.0D0)*(rax**(2.0D0))-(3.0D0/4.0D0)*(ray**(2.0D0))))
!
!  Charge-octopole     00, 30
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 3) .AND. (bk .EQ. 0)) Then
          !T0        = R4 * ((5.0D0/2.0D0)*(rbz**(3.0D0))-(3.0D0/2.0D0)*rbz)
          !T1(0)     = R6 * ((-35.0D0/4.0D0)*(rbz**(3.0D0))+(15.0D0/4.0D0)*rbz)
          !T1(6)     = R5 * ((15.0D0/2.0D0)*(rbz**(2.0D0))-(3.0D0/2.0D0))
          T0        = 1.0_wp*R4 * ((-5.0D0/2.0D0)*(rbz**(3.0D0))+(3.0D0/2.0D0)*rbz)
          T1(0)     = 1.0_wp*R6 * ((35.0D0/4.0D0)*(rbz**(3.0D0))-(15.0D0/4.0D0)*rbz)
          T1(6)     = 1.0_wp*R5 * ((-15.0D0/2.0D0)*(rbz**(2.0D0))+(3.0D0/2.0D0))
!
!                      00, 31c
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 3) .AND. (bk .EQ. 1)) Then
          !T0        = 1.0_wp*R4 * (rt6*((5.0D0/4.0D0)*rbx*(rbz**(2.0D0))-(1.0D0/4.0D0)*rbx))
          !T1(0)     = 1.0_wp*R6 * (rt6*((-35.0D0/8.0D0)*rbx*(rbz**(2.0D0))+(5.0D0/8.0D0)*rbx))
          !T1(4)     = 1.0_wp*R5 * (rt6*((5.0D0/4.0D0)*(rbz**(2.0D0))-(1.0D0/4.0D0)))
          !T1(6)     = 1.0_wp*R5 * (rt6*((5.0D0/2.0D0)*rbx*rbz))
          T0        = 1.0_wp*R4 * (rt6*((-5.0D0/4.0D0)*rbx*(rbz**(2.0D0))+(1.0D0/4.0D0)*rbx))
          T1(0)     = 1.0_wp*R6 * (rt6*((35.0D0/8.0D0)*rbx*(rbz**(2.0D0))-(5.0D0/8.0D0)*rbx))
          T1(4)     = 1.0_wp*R5 * (rt6*((-5.0D0/4.0D0)*(rbz**(2.0D0))+(1.0D0/4.0D0)))
          T1(6)     = 1.0_wp*R5 * (rt6*((-5.0D0/2.0D0)*rbx*rbz))
!        
!                      00, 31s
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 3) .AND. (bk .EQ. 2)) Then
          !T0        = 1.0_wp*R4 * (rt6*((5.0D0/4.0D0)*rby*(rbz**(2.0D0))-(1.0D0/4.0D0)*rby))
          !T1(0)     = 1.0_wp*R6 * (rt6*((-35.0D0/8.0D0)*rby*(rbz**(2.0D0))+(5.0D0/8.0D0)*rby))
          !T1(5)     = 1.0_wp*R5 * (rt6*((5.0D0/4.0D0)*(rbz**(2.0D0))-(1.0D0/4.0D0)))
          !T1(6)     = 1.0_wp*R5 * (rt6*((5.0D0/2.0D0)*rby*rbz))
          T0        = 1.0_wp*R4 * (rt6*((-5.0D0/4.0D0)*rby*(rbz**(2.0D0))+(1.0D0/4.0D0)*rby))
          T1(0)     = 1.0_wp*R6 * (rt6*((35.0D0/8.0D0)*rby*(rbz**(2.0D0))-(5.0D0/8.0D0)*rby))
          T1(5)     = 1.0_wp*R5 * (rt6*((-5.0D0/4.0D0)*(rbz**(2.0D0))+(1.0D0/4.0D0)))
          T1(6)     = 1.0_wp*R5 * (rt6*((-5.0D0/2.0D0)*rby*rbz))
!
!                      00, 32c
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 3) .AND. (bk .EQ. 3)) Then
          !T0        = 1.0_wp*R4 * (rt15*((1.0D0/2.0D0)*(rbx**(2.0D0))*rbz-(1.0D0/2.0D0)*(rby**(2.0D0))*rbz))
          !T1(0)     = 1.0_wp*R6 * (rt15*((-7.0D0/4.0D0)*(rbx**(2.0D0))*rbz+(7.0D0/4.0D0)*(rby**(2.0D0))*rbz))
          !T1(4)     = 1.0_wp*R5 * (rt15*rbx*rbz)
          !T1(5)     = 1.0_wp*R5 * (rt15*(-1.0D0)*rby*rbz)
          !T1(6)     = 1.0_wp*R5 * (rt15*((1.0D0/2.0D0)*(rbx**(2.0D0))-(1.0D0/2.0D0)*(rby**(2.0D0))))
          T0        = 1.0_wp*R4 * (rt15*((-1.0D0/2.0D0)*(rbx**(2.0D0))*rbz+(1.0D0/2.0D0)*(rby**(2.0D0))*rbz))
          T1(0)     = 1.0_wp*R6 * (rt15*((7.0D0/4.0D0)*(rbx**(2.0D0))*rbz-(7.0D0/4.0D0)*(rby**(2.0D0))*rbz))
          T1(4)     = 1.0_wp*R5 * (rt15*rbx*rbz)*(-1.0_wp)
          T1(5)     = 1.0_wp*R5 * (rt15*(1.0D0)*rby*rbz)
          T1(6)     = 1.0_wp*R5 * (rt15*((-1.0D0/2.0D0)*(rbx**(2.0D0))+(1.0D0/2.0D0)*(rby**(2.0D0))))
!      
!                      00, 32s
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 3) .AND. (bk .EQ. 4)) Then
          !T0        = 1.0_wp*R4 * (rt15*rbx*rby*rbz)
          !T1(0)     = 1.0_wp*R6 * (rt15*(-7.0D0/2.0D0)*rbx*rby*rbz)
          !T1(4)     = 1.0_wp*R5 * (rt15*rby*rbz)
          !T1(5)     = 1.0_wp*R5 * (rt15*rbx*rbz)
          !T1(6)     = 1.0_wp*R5 * (rt15*rbx*rby)
          T0        = 1.0_wp*R4 * (rt15*rbx*rby*rbz)*(-1.0_wp)
          T1(0)     = 1.0_wp*R6 * (rt15*(7.0D0/2.0D0)*rbx*rby*rbz)
          T1(4)     = 1.0_wp*R5 * (rt15*rby*rbz)*(-1.0_wp)
          T1(5)     = 1.0_wp*R5 * (rt15*rbx*rbz)*(-1.0_wp)
          T1(6)     = 1.0_wp*R5 * (rt15*rbx*rby)*(-1.0_wp)
!
!                      00, 33c
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 3) .AND. (bk .EQ. 5)) Then
          !T0        = 1.0_wp*R4 * (rt10*((1.0D0/4.0D0)*(rbx**(3.0D0))-(3.0D0/4.0D0)*rbx*(rby**(2.0D0))))
          !T1(0)     = 1.0_wp*R6 * (rt10*((-7.0D0/8.0D0)*(rbx**(3.0D0))+(21.0D0/8.0D0)*rbx*(rby**(2.0D0))))
          !T1(4)     = 1.0_wp*R5 * (rt10*((3.0D0/4.0D0)*(rbx**(2.0D0))-(3.0D0/4.0D0)*(rby**(2.0D0))))
          !T1(5)     = 1.0_wp*R5 * (rt10*(-3.0D0/2.0D0)*rbx*rby)
          T0        = 1.0_wp*R4 * (rt10*((-1.0D0/4.0D0)*(rbx**(3.0D0))+(3.0D0/4.0D0)*rbx*(rby**(2.0D0))))
          T1(0)     = 1.0_wp*R6 * (rt10*((7.0D0/8.0D0)*(rbx**(3.0D0))-(21.0D0/8.0D0)*rbx*(rby**(2.0D0))))
          T1(4)     = 1.0_wp*R5 * (rt10*((-3.0D0/4.0D0)*(rbx**(2.0D0))+(3.0D0/4.0D0)*(rby**(2.0D0))))
          T1(5)     = 1.0_wp*R5 * (rt10*(3.0D0/2.0D0)*rbx*rby)
!
!                      00, 33s
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 3) .AND. (bk .EQ. 6)) Then
          !T0        = 1.0_wp*R4 * (rt10*((3.0D0/4.0D0)*rby*(rbx**(2.0D0))-(1.0D0/4.0D0)*(rby**(3.0D0))))
          !T1(0)     = 1.0_wp*R6 * (rt10*((-21.0D0/8.0D0)*rby*(rbx**(2.0D0))+(7.0D0/8.0D0)*(rby**(3.0D0))))
          !T1(4)     = 1.0_wp*R5 * (rt10*(3.0D0/2.0D0)*rbx*rby)
          !T1(5)     = 1.0_wp*R5 * (rt10*((3.0D0/4.0D0)*(rbx**(2.0D0))-(3.0D0/4.0D0)*(rby**(2.0D0))))
          T0        = 1.0_wp*R4 * (rt10*((-3.0D0/4.0D0)*rby*(rbx**(2.0D0))+(1.0D0/4.0D0)*(rby**(3.0D0))))
          T1(0)     = 1.0_wp*R6 * (rt10*((21.0D0/8.0D0)*rby*(rbx**(2.0D0))-(7.0D0/8.0D0)*(rby**(3.0D0))))
          T1(4)     = 1.0_wp*R5 * (rt10*(-3.0D0/2.0D0)*rbx*rby)
          T1(5)     = 1.0_wp*R5 * (rt10*((-3.0D0/4.0D0)*(rbx**(2.0D0))+(3.0D0/4.0D0)*(rby**(2.0D0))))
!
!  Quadrupole-dipole   20, 10
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 0) .AND. (ib .EQ. 1) .AND. (bk .EQ. 0)) Then
          !T0           = 1.0_wp*R4 * (7.5D0*raz*raz*rBeta + 3.0D0*raz*czB - 1.5D0*rBeta)
          !T1(0)        = 1.0_wp*R6 * (-26.25D0*raz*raz*rBeta - 7.5D0*raz*czB + 3.75*rBeta)
          !T1(3)        = 1.0_wp*R5 * (15.0D0*raz*rBeta + 3.0D0*czB)
          !T1(4+ob)     = 1.0_wp*R5 * (7.5D0*(raz*raz) - 1.5D0)
          !T1(9+(3*ob)) = 1.0_wp*R4 * (3.0D0*raz)
          T0           = 1.0_wp*R4 * ((-7.5D0)*raz*raz*rBeta + 3.0D0*raz*czB + 1.5D0*rBeta)
          T1(0)        = 1.0_wp*R6 * (26.25D0*raz*raz*rBeta - 7.5D0*raz*czB - 3.75_wp*rBeta)
          T1(3)        = 1.0_wp*R5 * (-15.0D0*raz*rBeta + 3.0D0*czB)
          T1(4+ob)     = 1.0_wp*R5 * (-7.5D0*(raz*raz) + 1.5D0)
          T1(9+(3*ob)) = 1.0_wp*R4 * (3.0D0*raz)
!
!                      21c, 10
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 1) .AND. (ib .EQ. 1) .AND. (bk .EQ. 0)) Then
          !T0           = 1.0_wp*R4 * rt3 * (rax*czB + raz*cxB + 5.0D0*rax*raz*rBeta)
          !T1(0)        = 0.0_wp*R6 * rt3 * (-2.5D0*rax*czB - 2.5D0*raz*cxB - 17.5D0*rax*raz*rBeta)
          !T1(1)        = 0.0_wp*R5 * rt3 * (czB + 5.0D0*raz*rBeta)
          !T1(3)        = 0.0_wp*R5 * rt3 * (cxB + 5.0D0*rax*rBeta)
          !T1(4+ob)     = 0.0_wp*R5 * rt3 * (5.0D0*rax*raz)
          !T1(7+(3*ob)) = 0.0_wp*R4 * rt3 * raz
          !T1(9+(3*ob)) = 0.0_wp*R4 * rt3 * rax
          T0           = 1.0_wp*R4 * rt3 * (rax*czB + raz*cxB - 5.0D0*rax*raz*rBeta)
          T1(0)        = 1.0_wp*R6 * rt3 * (-2.5D0*rax*czB - 2.5D0*raz*cxB + 17.5D0*rax*raz*rBeta)
          T1(1)        = 1.0_wp*R5 * rt3 * (czB - 5.0D0*raz*rBeta)
          T1(3)        = 1.0_wp*R5 * rt3 * (cxB - 5.0D0*rax*rBeta)
          T1(4+ob)     = 1.0_wp*R5 * rt3 * (-5.0D0*rax*raz)
          T1(7+(3*ob)) = 1.0_wp*R4 * rt3 * raz
          T1(9+(3*ob)) = 1.0_wp*R4 * rt3 * rax
!
!                      21s, 10
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 2) .AND. (ib .EQ. 1) .AND. (bk .EQ. 0)) Then
          !T0           = 1.0_wp*R4 * rt3 * (ray*czB + raz*cyB + 5.0D0*ray*raz*rBeta)
          !T1(0)        = 0.0_wp*R6 * rt3 * (-2.5D0*ray*czB - 2.5D0*raz*cyB - 17.5D0*ray*raz*rBeta)
          !T1(2)        = 0.0_wp*R5 * rt3 * (czB + 5.0D0*raz*rBeta)
          !T1(3)        = 0.0_wp*R5 * rt3 * (cyB + 5.0D0*ray*rBeta)
          !T1(4+ob)     = 0.0_wp*R5 * rt3 * (5.0D0*ray*raz)
          !T1(8+(3*ob)) = 0.0_wp*R4 * rt3 * raz
          !T1(9+(3*ob)) = 0.0_wp*R4 * rt3 * ray
          T0           = 1.0_wp*R4 * rt3 * (ray*czB + raz*cyB - 5.0D0*ray*raz*rBeta)
          T1(0)        = 1.0_wp*R6 * rt3 * (-2.5D0*ray*czB - 2.5D0*raz*cyB + 17.5D0*ray*raz*rBeta)
          T1(2)        = 1.0_wp*R5 * rt3 * (czB - 5.0D0*raz*rBeta)
          T1(3)        = 1.0_wp*R5 * rt3 * (cyB - 5.0D0*ray*rBeta)
          T1(4+ob)     = 1.0_wp*R5 * rt3 * (-5.0D0*ray*raz)
          T1(8+(3*ob)) = 1.0_wp*R4 * rt3 * raz
          T1(9+(3*ob)) = 1.0_wp*R4 * rt3 * ray
!
!                      22c, 10
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 3) .AND. (ib .EQ. 1) .AND. (bk .EQ. 0)) Then
          !T0           = 1.0_wp*R4 * rt3 * (2.5D0*rBeta*rax*rax - 2.5D0*rBeta*ray*ray + rax*cxB - ray*cyB)
          !T1(0)        = 0.0_wp*R6 * rt3 * (8.75D0*rBeta*ray*ray - 8.75D0*rBeta*rax*rax - 2.5D0*rax*cxB + 2.5D0*ray*cyB)
          !T1(1)        = 0.0_wp*R5 * rt3 * ( 5.0D0*rBeta*rax + cxB)
          !T1(2)        = 0.0_wp*R5 * rt3 * (-5.0D0*rBeta*ray - cyB)
          !T1(4+ob)     = 0.0_wp*R5 * rt3 * 2.5D0 * (rax*rax - ray*ray)
          !T1(7+(3*ob)) = 0.0_wp*R4 * rt3 * ( rax)
          !T1(8+(3*ob)) = 0.0_wp*R4 * rt3 * (-ray)
          T0           = 1.0_wp*R4 * rt3 * ((-2.5D0)*rBeta*rax*rax + 2.5D0*rBeta*ray*ray + rax*cxB - ray*cyB)
          T1(0)        = 1.0_wp*R6 * rt3 * ((-8.75D0)*rBeta*ray*ray + 8.75D0*rBeta*rax*rax - 2.5D0*rax*cxB + 2.5D0*ray*cyB)
          T1(1)        = 1.0_wp*R5 * rt3 * ((-5.0D0)*rBeta*rax + cxB)
          T1(2)        = 1.0_wp*R5 * rt3 * ( 5.0D0*rBeta*ray - cyB)
          T1(4+ob)     = 1.0_wp*R5 * rt3 * 2.5D0 * (ray*ray - rax*rax)
          T1(7+(3*ob)) = 1.0_wp*R4 * rt3 * (rax)
          T1(8+(3*ob)) = 1.0_wp*R4 * rt3 * (-ray)
!
!                      22s, 10
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 4) .AND. (ib .EQ. 1) .AND. (bk .EQ. 0)) Then
          !T0           = 1.0_wp*R4 * rt3 * (ray*cxB + rax*cyB + 5.0D0*ray*rax*rBeta)
          !T1(0)        = 0.0_wp*R6 * rt3 * (-2.5D0*ray*cxB - 2.5D0*rax*cyB - 17.5D0*ray*rax*rBeta)
          !T1(1)        = 0.0_wp*R5 * rt3 * (cyB + 5.0D0*ray*rBeta)
          !T1(2)        = 0.0_wp*R5 * rt3 * (cxB + 5.0D0*rax*rBeta)
          !T1(4+ob)     = 0.0_wp*R5 * rt3 * (5.0D0*ray*rax)
          !T1(7+(3*ob)) = 0.0_wp*R4 * rt3 * ray
          !T1(8+(3*ob)) = 0.0_wp*R4 * rt3 * rax
          T0           = 1.0_wp*R4 * rt3 * (ray*cxB + rax*cyB - 5.0D0*ray*rax*rBeta)
          T1(0)        = 1.0_wp*R6 * rt3 * (-2.5D0*ray*cxB - 2.5D0*rax*cyB + 17.5D0*ray*rax*rBeta)
          T1(1)        = 1.0_wp*R5 * rt3 * (cyB - 5.0D0*ray*rBeta)
          T1(2)        = 1.0_wp*R5 * rt3 * (cxB - 5.0D0*rax*rBeta)
          T1(4+ob)     = 1.0_wp*R5 * rt3 * (-5.0D0*ray*rax)
          T1(7+(3*ob)) = 1.0_wp*R4 * rt3 * ray
          T1(8+(3*ob)) = 1.0_wp*R4 * rt3 * rax
!
!  Dipole-quadrupole   10, 20
!
        Else If ((ia .EQ. 1) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 0)) Then
          !T0           = 1.0_wp*R4 * (7.5D0*rbz*rbz*rAlpha + 3.0D0*rbz*czA - 1.5D0*rAlpha)
          !T1(0)        = 1.0_wp*R6 * (-26.25D0*rbz*rbz*rAlpha - 7.5D0*rbz*czA + 3.75*rAlpha)
          !T1(6)        = 1.0_wp*R5 * (15.0D0*rbz*rAlpha + 3.0D0*czA)
          !T1(1+oa)     = 1.0_wp*R5 * (7.5D0*(rbz*rbz) - 1.5D0)
          !T1(9+(3*oa)) = 1.0_wp*R4 * (3.0D0*rbz)
          T0           = 1.0_wp*R4 * (7.5D0*rbz*rbz*rAlpha - 3.0D0*rbz*czA - 1.5D0*rAlpha)
          T1(0)        = 1.0_wp*R6 * (-26.25D0*rbz*rbz*rAlpha + 7.5D0*rbz*czA + 3.75_wp*rAlpha)
          T1(6)        = 1.0_wp*R5 * (15.0D0*rbz*rAlpha - 3.0D0*czA)
          T1(1+oa)     = 1.0_wp*R5 * (7.5D0*(rbz*rbz) - 1.5D0)
          T1(13+oa) = -1.0_wp*R4*(3.0D0*rbz)
!
!                      10, 21c
!
        Else If ((ia .EQ. 1) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 1)) Then
          !T0           = 1.0_wp*R4 * rt3 * (rbx*czA + rbz*cxA + 5.0D0*rbx*rbz*rAlpha)
          !T1(0)        = 0.0_wp*R6 * rt3 * (-2.5D0*rbx*czA - 2.5D0*rbz*cxA - 17.5D0*rbx*rbz*rAlpha)
          !T1(4)        = 0.0_wp*R5 * rt3 * (czA + 5.0D0*rbz*rAlpha)
          !T1(6)        = 0.0_wp*R5 * rt3 * (cxA + 5.0D0*rbx*rAlpha)
          !T1(1+oa)     = 0.0_wp*R5 * rt3 * (5.0D0*rbx*rbz)
          !T1(7+(3*oa)) = 0.0_wp*R4 * rt3 * rbz
          !T1(9+(3*oa)) = 0.0_wp*R4 * rt3 * rbx
          T0           = 1.0_wp*R4 * rt3 * ((-1.0_wp)*rbx*czA - rbz*cxA + 5.0D0*rbx*rbz*rAlpha)
          T1(0)        = 1.0_wp*R6 * rt3 * (2.5D0*rbx*czA + 2.5D0*rbz*cxA - 17.5D0*rbx*rbz*rAlpha)
          T1(4)        = 1.0_wp*R5 * rt3 * ((-1.0_wp)*czA + 5.0D0*rbz*rAlpha)
          T1(6)        = 1.0_wp*R5 * rt3 * ((-1.0_wp)*cxA + 5.0D0*rbx*rAlpha)
          T1(1+oa)     = 1.0_wp*R5 * rt3 * (5.0D0*rbx*rbz)
          T1(7+oa)  = -1.0_wp*R4 * rt3 * rbz
          T1(13+oa) = -1.0_wp*R4 * rt3 * rbx
!
!                      10,21s
!
        Else If ((ia .EQ. 1) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 2)) Then
          !T0           = 1.0_wp*R4 * rt3 * (rby*czA + rbz*cyA + 5.0D0*rby*rbz*rAlpha)
          !T1(0)        = 0.0_wp*R6 * rt3 * (-2.5D0*rby*czA - 2.5D0*rbz*cyA - 17.5D0*rby*rbz*rAlpha)
          !T1(5)        = 0.0_wp*R5 * rt3 * (czA + 5.0D0*rbz*rAlpha)
          !T1(6)        = 0.0_wp*R5 * rt3 * (cyA + 5.0D0*rby*rAlpha)
          !T1(1+oa)     = 0.0_wp*R5 * rt3 * (5.0D0*rby*rbz)
          !T1(8+(3*oa)) = 0.0_wp*R4 * rt3 * rbz
          !T1(9+(3*oa)) = 0.0_wp*R4 * rt3 * rby
          T0           = 1.0_wp*R4 * rt3 * ((-1.0_wp)*rby*czA - rbz*cyA + 5.0D0*rby*rbz*rAlpha)
          T1(0)        = 1.0_wp*R6 * rt3 * (+2.5D0*rby*czA + 2.5D0*rbz*cyA - 17.5D0*rby*rbz*rAlpha)
          T1(5)        = 1.0_wp*R5 * rt3 * ((-1.0_wp)*czA + 5.0D0*rbz*rAlpha)
          T1(6)        = 1.0_wp*R5 * rt3 * ((-1.0_wp)*cyA + 5.0D0*rby*rAlpha)
          T1(1+oa)     = 1.0_wp*R5 * rt3 * (5.0D0*rby*rbz)
          T1(10+oa) = -1.0_wp*R4 * rt3 * rbz
          T1(13+oa) = -1.0_wp*R4 * rt3 * rby

!                      10, 22c
!
        Else If ((ia .EQ. 1) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 3)) Then
          !T0           = 1.0_wp*R4 * rt3 * (2.5D0*rAlpha*rbx*rbx - 2.5D0*rAlpha*rby*rby + rbx*cxA - rby*cyA)
          !T1(0)        = 0.0_wp*R6 * rt3 * (8.75D0*rAlpha*rby*rby - 8.75D0*rAlpha*rbx*rbx - 2.5D0*rbx*cxA + 2.5D0*rby*cyA)
          !T1(4)        = 0.0_wp*R5 * rt3 * ( 5.0D0*rAlpha*rbx + cxA)
          !T1(5)        = 0.0_wp*R5 * rt3 * (-5.0D0*rAlpha*rby - cyA)
          !T1(1+oa)     = 0.0_wp*R5 * rt3 * 2.5D0 * (rbx*rbx - rby*rby)
          !T1(7+(3*oa)) = 0.0_wp*R4 * rt3 * ( rbx)
          !T1(8+(3*oa)) = 0.0_wp*R4 * rt3 * (-rby)
          T0           = 1.0_wp*R4 * rt3 * (2.5D0*rAlpha*rbx*rbx - 2.5D0*rAlpha*rby*rby - rbx*cxA + rby*cyA)
          T1(0)        = 1.0_wp*R6 * rt3 * (8.75D0*rAlpha*rby*rby - 8.75D0*rAlpha*rbx*rbx + 2.5D0*rbx*cxA - 2.5D0*rby*cyA)
          T1(4)        = 1.0_wp*R5 * rt3 * ( 5.0D0*rAlpha*rbx - cxA)
          T1(5)        = 1.0_wp*R5 * rt3 * (-5.0D0*rAlpha*rby + cyA)
          T1(1+oa)     = 1.0_wp*R5 * rt3 * 2.5D0 * (rbx*rbx - rby*rby)
          T1(7+oa)     = 1.0_wp*R4 * rt3 * (-rbx)
          T1(10+oa)    = 1.0_wp*R4 * rt3 * (rby)
!
!                      10, 22s
!
        Else If ((ia .EQ. 1) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 4)) Then
          !T0           = 1.0_wp*R4 * rt3 * (rby*cxA + rbx*cyA + 5.0D0*rby*rbx*rAlpha)
          !T1(0)        = 0.0_wp*R6 * rt3 * (-2.5D0*rby*cxA - 2.5D0*rbx*cyA - 17.5D0*rby*rbx*rAlpha)
          !T1(4)        = 0.0_wp*R5 * rt3 * (cyA + 5.0D0*rby*rAlpha)
          !T1(5)        = 0.0_wp*R5 * rt3 * (cxA + 5.0D0*rbx*rAlpha)
          !T1(1+oa)     = 0.0_wp*R5 * rt3 * (5.0D0*rby*rbx)
          !T1(7+(3*oa)) = 0.0_wp*R4 * rt3 * rby
          !T1(8+(3*oa)) = 0.0_wp* R4 * rt3 * rbx
          T0           = 1.0_wp*R4 * rt3 * ((-1.0_wp)*rby*cxA - rbx*cyA + 5.0D0*rby*rbx*rAlpha)
          T1(0)        = 1.0_wp*R6 * rt3 * (2.5D0*rby*cxA + 2.5D0*rbx*cyA - 17.5D0*rby*rbx*rAlpha)
          T1(4)        = 1.0_wp*R5 * rt3 * ((-1.0_wp)*cyA + 5.0D0*rby*rAlpha)
          T1(5)        = 1.0_wp*R5 * rt3 * ((-1.0_wp)*cxA + 5.0D0*rbx*rAlpha)
          T1(1+oa)     = 1.0_wp*R5 * rt3 * (5.0D0*rby*rbx)
          T1(7+oa)     = -1.0_wp*R4 * rt3 * rby
          T1(10+oa)    = -1.0_wp*R4 * rt3 * rbx
!
!  L = 5
!  Hexadecapole-charge 40, 00
!
        Else If ((ia .EQ. 4) .AND. (ak .EQ. 0) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = R5 * ((35.0D0/8.0D0)*(raz**(4.0D0))-(15.0D0/4.0D0)*(raz**(2.0D0))+(3.0D0/8.0D0))
          T1(0)     = R7 * ((-315.0D0/16.0D0)*(raz**(4.0D0))+(105.0D0/8.0D0)*(raz**(2.0D0))-(15.0D0/16.0D0))
          T1(3)     = R6 * ((35.0D0/2.0D0)*(raz**(3.0D0))-(15.0D0/2.0D0)*raz)
!
!                      41c, 00
!
        Else If ((ia .EQ. 4) .AND. (ak .EQ. 1) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = R5 * (rt10*((7.0D0/4.0D0)*rax*(raz**(3.0D0))-(3.0D0/4.0D0)*rax*raz))
          T1(0)     = R7 * (rt10*((-63.0D0/8.0D0)*rax*(raz**(3.0D0))+(21.0D0/8.0D0)*rax*raz))
          T1(1)     = R6 * (rt10*((7.0D0/4.0D0)*(raz**(3.0D0))-(3.0D0/4.0D0)*raz))
          T1(3)     = R6 * (rt10*((21.0D0/4.0D0)*rax*(raz**(2.0D0))-(3.0D0/4.0D0)*rax))
!
!                      41s, 00
!
        Else If ((ia .EQ. 4) .AND. (ak .EQ. 2) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = R5 * (rt10*((7.0D0/4.0D0)*ray*(raz**(3.0D0))-(3.0D0/4.0D0)*ray*raz))
          T1(0)     = R7 * (rt10*((-63.0D0/8.0D0)*ray*(raz**(3.0D0))+(21.0D0/8.0D0)*ray*raz))
          T1(2)     = R6 * (rt10*((7.0D0/4.0D0)*(raz**(3.0D0))-(3.0D0/4.0D0)*raz))
          T1(3)     = R6 * (rt10*((21.0D0/4.0D0)*ray*(raz**(2.0D0))-(3.0D0/4.0D0)*ray))
!
!                      42c, 00
!
        Else If ((ia .EQ. 4) .AND. (ak .EQ. 3) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = R5 * (rt5*((7.0D0/4.0D0)*(rax**(2.0D0))*(raz**(2.0D0))-(7.0D0/4.0D0)*(ray**(2.0D0))*(raz**(2.0D0)) &
                                             & -(1.0D0/4.0D0)*(rax**(2.0D0))+(1.0D0/4.0D0)*(ray**(2.0D0))))
          T1(0)     = R7 * (rt5*((-63.0D0/8.0D0)*(rax**(2.0D0))*(raz**(2.0D0))+(63.0D0/8.0D0)*(ray**(2.0D0))*(raz**(2.0D0)) &
                                             & +(7.0D0/8.0D0)*(rax**(2.0D0))-(7.0D0/8.0D0)*(ray**(2.0D0))))
          T1(1)     = R6 * (rt5*((7.0D0/2.0D0)*rax*(raz**(2.0D0))-(1.0D0/2.0D0)*rax))
          T1(2)     = R6 * (rt5*((-7.0D0/2.0D0)*ray*(raz**(2.0D0))+(1.0D0/2.0D0)*ray))
          T1(3)     = R6 * (rt5*((7.0D0/2.0D0)*(rax**(2.0D0))*raz-(7.0D0/2.0D0)*(ray**(2.0D0))*raz))
!
!                      42s, 00
!
        Else If ((ia .EQ. 4) .AND. (ak .EQ. 4) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = R5 * (rt5*((7.0D0/2.0D0)*rax*ray*(raz**(2.0D0))-(1.0D0/2.0D0)*rax*ray))
          T1(0)     = R7 * (rt5*((-63.0D0/4.0D0)*rax*ray*(raz**(2.0D0))+(7.0D0/4.0D0)*rax*ray))
          T1(1)     = R6 * (rt5*((7.0D0/2.0D0)*ray*(raz**(2.0D0))-(1.0D0/2.0D0)*ray))
          T1(2)     = R6 * (rt5*((7.0D0/2.0D0)*rax*(raz**(2.0D0))-(1.0D0/2.0D0)*rax))
          T1(3)     = R6 * (rt5*(7.0D0*rax*ray*raz))
!
!                      43c, 00
!
        Else If ((ia .EQ. 4) .AND. (ak .EQ. 5) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = R5 * (rt70*((1.0D0/4.0D0)*(rax**(3.0D0))*raz-(3.0D0/4.0D0)*rax*(ray**(2.0D0))*raz))
          T1(0)     = R7 * (rt70*((-9.0D0/8.0D0)*(rax**(3.0D0))*raz+(27.0D0/8.0D0)*rax*(ray**(2.0D0))*raz))
          T1(1)     = R6 * (rt70*((3.0D0/4.0D0)*(rax**(2.0D0))*raz-(3.0D0/4.0D0)*(ray**(2.0D0))*raz))
          T1(2)     = R6 * (rt70*((-3.0D0/2.0D0)*rax*ray*raz))
          T1(3)     = R6 * (rt70*((1.0D0/4.0D0)*(rax**(3.0D0))-(3.0D0/4.0D0)*rax*(ray**(2.0D0))))
!
!                      43s, 00
!
        Else If ((ia .EQ. 4) .AND. (ak .EQ. 6) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = R5 * (rt70*((3.0D0/4.0D0)*(rax**(2.0D0))*ray*raz-(1.0D0/4.0D0)*(ray**(3.0D0))*raz))
          T1(0)     = R7 * (rt70*((-27.0D0/8.0D0)*(rax**(2.0D0))*ray*raz+(9.0D0/8.0D0)*(ray**(3.0D0))*raz))
          T1(1)     = R6 * (rt70*((3.0D0/2.0D0)*rax*ray*raz))
          T1(2)     = R6 * (rt70*((3.0D0/4.0D0)*(rax**(2.0D0))*raz-(3.0D0/4.0D0)*(ray**(2.0D0))*raz))
          T1(3)     = R6 * (rt70*((3.0D0/4.0D0)*(rax**(2.0D0))*ray-(1.0D0/4.0D0)*(ray**(3.0D0))))
!
!                      44c, 00
!
        Else If ((ia .EQ. 4) .AND. (ak .EQ. 7) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = R5 * (rt35*((1.0D0/8.0D0)*(rax**(4.0D0))-(3.0D0/4.0D0)*(rax**(2.0D0))*(ray**(2.0D0))+(1.0D0/8.0D0)*(ray**(4.0D0))))
          T1(0)     = R7 * (rt35*((-9.0D0/16.0D0)*(rax**(4.0D0))+(27.0D0/8.0D0)*(rax**(2.0D0))*(ray**(2.0D0))-(9.0D0/16.0D0)*(ray**(4.0D0))))
          T1(1)     = R6 * (rt35*((1.0D0/2.0D0)*(rax**(3.0D0))-(3.0D0/2.0D0)*rax*(ray**(2.0D0))))
          T1(2)     = R6 * (rt35*((-3.0D0/2.0D0)*(rax**(2.0D0))*ray+(1.0D0/2.0D0)*(ray**(3.0D0))))
!
!                      44s, 00
!
        Else If ((ia .EQ. 4) .AND. (ak .EQ. 8) .AND. (ib .EQ. 0) .AND. (bk .EQ. 0)) Then
          T0        = R5 * (rt35*((1.0D0/2.0D0)*(rax**(3.0D0))*ray-(1.0D0/2.0D0)*rax*(ray**(3.0D0))))
          T1(0)     = R7 * (rt35*((-9.0D0/4.0D0)*(rax**(3.0D0))*ray+(9.0D0/4.0D0)*rax*(ray**(3.0D0))))
          T1(1)     = R6 * (rt35*((3.0D0/2.0D0)*(rax**(2.0D0))*ray-(1.0D0/2.0D0)*(ray**(3.0D0))))
          T1(2)     = R6 * (rt35*((1.0D0/2.0D0)*(rax**(3.0D0))-(3.0D0/2.0D0)*rax*(ray**(2.0D0))))
!
!  Charge-hexadecapole 00, 40
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 4) .AND. (bk .EQ. 0)) Then
          T0        = R5 * ((35.0D0/8.0D0)*(rbz**(4.0D0))-(15.0D0/4.0D0)*(rbz**(2.0D0))+(3.0D0/8.0D0))
          T1(0)     = R7 * ((-315.0D0/16.0D0)*(rbz**(4.0D0))+(105.0D0/8.0D0)*(rbz**(2.0D0))-(15.0D0/16.0D0))
          T1(6)     = R6 * ((35.0D0/2.0D0)*(rbz**(3.0D0))-(15.0D0/2.0D0)*rbz)
!
!                      00, 41c
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 4) .AND. (bk .EQ. 1)) Then
          T0        = R5 * (rt10*((7.0D0/4.0D0)*rbx*(rbz**(3.0D0))-(3.0D0/4.0D0)*rbx*rbz))
          T1(0)     = R7 * (rt10*((-63.0D0/8.0D0)*rbx*(rbz**(3.0D0))+(21.0D0/8.0D0)*rbx*rbz))
          T1(4)     = R6 * (rt10*((7.0D0/4.0D0)*(rbz**(3.0D0))-(3.0D0/4.0D0)*rbz))
          T1(6)     = R6 * (rt10*((21.0D0/4.0D0)*rbx*(rbz**(2.0D0))-(3.0D0/4.0D0)*rbx))
!
!                      00, 41s
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 4) .AND. (bk .EQ. 2)) Then
          T0        = R5 * (rt10*((7.0D0/4.0D0)*rby*(rbz**(3.0D0))-(3.0D0/4.0D0)*rby*rbz))
          T1(0)     = R7 * (rt10*((-63.0D0/8.0D0)*rby*(rbz**(3.0D0))+(21.0D0/8.0D0)*rby*rbz))
          T1(5)     = R6 * (rt10*((7.0D0/4.0D0)*(rbz**(3.0D0))-(3.0D0/4.0D0)*rbz))
          T1(6)     = R6 * (rt10*((21.0D0/4.0D0)*rby*(rbz**(2.0D0))-(3.0D0/4.0D0)*rby))
!
!                      00, 42c
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 4) .AND. (bk .EQ. 3)) Then
          T0        = R5 * (rt5*((7.0D0/4.0D0)*(rbx**(2.0D0))*(rbz**(2.0D0))-(7.0D0/4.0D0)*(rby**(2.0D0))*(rbz**(2.0D0)) &
                                      & -(1.0D0/4.0D0)*(rbx**(2.0D0))+(1.0D0/4.0D0)*(rby**(2.0D0))))
          T1(0)     = R7 * (rt5*((-63.0D0/8.0D0)*(rbx**(2.0D0))*(rbz**(2.0D0))+(63.0D0/8.0D0)*(rby**(2.0D0))*(rbz**(2.0D0)) &
                                      & +(7.0D0/8.0D0)*(rbx**(2.0D0))-(7.0D0/8.0D0)*(rby**(2.0D0))))
          T1(4)     = R6 * (rt5*((7.0D0/2.0D0)*rbx*(rbz**(2.0D0))-(1.0D0/2.0D0)*rbx))
          T1(5)     = R6 * (rt5*((-7.0D0/2.0D0)*rby*(rbz**(2.0D0))+(1.0D0/2.0D0)*rby))
          T1(6)     = R6 * (rt5*((7.0D0/2.0D0)*(rbx**(2.0D0))*rbz-(7.0D0/2.0D0)*(rby**(2.0D0))*rbz))
!
!                      00, 42s
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 4) .AND. (bk .EQ. 4)) Then
          T0        = R5 * (rt5*((7.0D0/2.0D0)*rbx*rby*(rbz**(2.0D0))-(1.0D0/2.0D0)*rbx*rby))
          T1(0)     = R7 * (rt5*((-63.0D0/4.0D0)*rbx*rby*(rbz**(2.0D0))+(7.0D0/4.0D0)*rbx*rby))
          T1(4)     = R6 * (rt5*((7.0D0/2.0D0)*rby*(rbz**(2.0D0))-(1.0D0/2.0D0)*rby))
          T1(5)     = R6 * (rt5*((7.0D0/2.0D0)*rbx*(rbz**(2.0D0))-(1.0D0/2.0D0)*rbx))
          T1(6)     = R6 * (rt5*(7.0D0*rbx*rby*rbz))
!
!                      00, 43c
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 4) .AND. (bk .EQ. 5)) Then
          T0        = R5 * (rt70*((1.0D0/4.0D0)*(rbx**(3.0D0))*rbz-(3.0D0/4.0D0)*rbx*(rby**(2.0D0))*rbz))
          T1(0)     = R7 * (rt70*((-9.0D0/8.0D0)*(rbx**(3.0D0))*rbz+(27.0D0/8.0D0)*rbx*(rby**(2.0D0))*rbz))
          T1(4)     = R6 * (rt70*((3.0D0/4.0D0)*(rbx**(2.0D0))*rbz-(3.0D0/4.0D0)*(rby**(2.0D0))*rbz))
          T1(5)     = R6 * (rt70*((-3.0D0/2.0D0)*rbx*rby*rbz))
          T1(6)     = R6 * (rt70*((1.0D0/4.0D0)*(rbx**(3.0D0))-(3.0D0/4.0D0)*rbx*(rby**(2.0D0))))
!
!                      00, 43s
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 4) .AND. (bk .EQ. 6)) Then
          T0        = R5 * (rt70*((3.0D0/4.0D0)*(rbx**(2.0D0))*rby*rbz-(1.0D0/4.0D0)*(rby**(3.0D0))*rbz))
          T1(0)     = R7 * (rt70*((-27.0D0/8.0D0)*(rbx**(2.0D0))*rby*rbz+(9.0D0/8.0D0)*(rby**(3.0D0))*rbz))
          T1(4)     = R6 * (rt70*((3.0D0/2.0D0)*rbx*rby*rbz))
          T1(5)     = R6 * (rt70*((3.0D0/4.0D0)*(rbx**(2.0D0))*rbz-(3.0D0/4.0D0)*(rby**(2.0D0))*rbz))
          T1(6)     = R6 * (rt70*((3.0D0/4.0D0)*(rbx**(2.0D0))*rby-(1.0D0/4.0D0)*(rby**(3.0D0))))
!
!                      00, 44c
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 4) .AND. (bk .EQ. 7)) Then
          T0        = R5 * (rt35*((1.0D0/8.0D0)*(rbx**(4.0D0))-(3.0D0/4.0D0)*(rbx**(2.0D0))*(rby**(2.0D0))+(1.0D0/8.0D0)*(rby**(4.0D0))))
          T1(0)     = R7 * (rt35*((-9.0D0/16.0D0)*(rbx**(4.0D0))+(27.0D0/8.0D0)*(rbx**(2.0D0))*(rby**(2.0D0))-(9.0D0/16.0D0)*(rby**(4.0D0))))
          T1(4)     = R6 * (rt35*((1.0D0/2.0D0)*(rbx**(3.0D0))-(3.0D0/2.0D0)*rbx*(rby**(2.0D0))))
          T1(5)     = R6 * (rt35*((-3.0D0/2.0D0)*(rbx**(2.0D0))*rby+(1.0D0/2.0D0)*(rby**(3.0D0))))
!
!                      00, 44s
!
        Else If ((ia .EQ. 0) .AND. (ak .EQ. 0) .AND. (ib .EQ. 4) .AND. (bk .EQ. 8)) Then
          T0        = R5 * (rt35*((1.0D0/2.0D0)*(rbx**(3.0D0))*rby-(1.0D0/2.0D0)*rbx*(rby**(3.0D0))))
          T1(0)     = R7 * (rt35*((-9.0D0/4.0D0)*(rbx**(3.0D0))*rby+(9.0D0/4.0D0)*rbx*(rby**(3.0D0))))
          T1(4)     = R6 * (rt35*((3.0D0/2.0D0)*(rbx**(2.0D0))*rby-(1.0D0/2.0D0)*(rby**(3.0D0))))
          T1(5)     = R6 * (rt35*((1.0D0/2.0D0)*(rbx**(3.0D0))-(3.0D0/2.0D0)*rbx*(rby**(2.0D0))))
!
!  Octopole-dipole     30, 10
!
        Else If ((ia .EQ. 3) .AND. (ak .EQ. 0) .AND. (ib .EQ. 1) .AND. (bk .EQ. 0)) Then
          !T0           = 0.0_wp*R5 * (17.5D0*raz*raz*raz*rBeta + 7.5D0*raz*raz*czB - 7.5D0*raz*rBeta - 1.5D0*czB)
          !T1(0)        = 0.0_wp*R7 * (-78.75D0*raz*raz*raz*rBeta - 26.25D0*raz*raz*czB + 26.25D0*raz*rBeta + 3.75D0*czB)
          !T1(3)        = 0.0_wp*R6 * (52.5D0*raz*raz*rBeta + 15.0D0*raz*czB - 7.5D0*rBeta)
          !T1(4+ob)     = 0.0_wp*R6 * (17.5D0*raz*raz*raz - 7.5D0*raz)
          !T1(9+(3*ob)) = 0.0_wp*R5 * (7.5D0*raz*raz - 1.5D0)
          T0           = 1.0_wp*R5 * (-17.5D0*raz*raz*raz*rBeta + 7.5D0*raz*raz*czB + 7.5D0*raz*rBeta - 1.5D0*czB)
          T1(0)        = 1.0_wp*R7 * (78.75D0*raz*raz*raz*rBeta - 26.25D0*raz*raz*czB - 26.25D0*raz*rBeta + 3.75D0*czB)
          T1(3)        = 1.0_wp*R6 * (-52.5D0*raz*raz*rBeta + 15.0D0*raz*czB + 7.5D0*rBeta)
          T1(4+ob)     = 1.0_wp*R6 * (-17.5D0*raz*raz*raz + 7.5D0*raz)
          T1(9+(3*ob)) = 1.0_wp*R5 * (7.5D0*raz*raz - 1.5D0)

!
!                      31c, 10
!
        Else If ((ia .EQ. 3) .AND. (ak .EQ. 1) .AND. (ib .EQ. 1) .AND. (bk .EQ. 0)) Then
          !T0           = 0.0_wp*R5 * rt6 * (8.75D0*rax*raz*raz*rBeta + 1.25D0*raz*raz*cxB &
          !                      & + 2.5D0*rax*raz*czB - 1.25D0*rax*rBeta - 0.25D0*cxB)
          !T1(0)        = 0.0_wp*R7 * rt6 * (-39.375D0*rax*raz*raz*rBeta - 4.375D0*raz*raz*cxB &
          !                         & - 8.75D0*rax*raz*czB + 4.375D0*rax*rBeta + 0.625D0*cxB)
          !T1(1)        = 0.0_wp*R6 * rt6 * (8.75D0*raz*raz*rBeta + 2.5D0*raz*czB - 1.25D0*rBeta)
          !T1(3)        = 0.0_wp*R6 * rt6 * (17.5D0*rax*raz*rBeta + 2.5D0*raz*cxB + 2.5D0*rax*czB)
          !T1(4+ob)     = 0.0_wp*R6 * rt6 * (8.75D0*rax*raz*raz - 1.25D0*rax)
          !T1(7+(3*ob)) = 0.0_wp*R5 * rt6 * (1.25D0*raz*raz - 0.25D0)
          !T1(9+(3*ob)) = 0.0_wp*R5 * rt6 * (2.5D0*rax*raz)
          T0           = 1.0_wp*R5 * rt6 * (-8.75D0*rax*raz*raz*rBeta + 1.25D0*raz*raz*cxB &
                                & + 2.5D0*rax*raz*czB + 1.25D0*rax*rBeta - 0.25D0*cxB)
          T1(0)        = 1.0_wp*R7 * rt6 * (39.375D0*rax*raz*raz*rBeta - 4.375D0*raz*raz*cxB &
                                   & - 8.75D0*rax*raz*czB - 4.375D0*rax*rBeta + 0.625D0*cxB)
          T1(1)        = 1.0_wp*R6 * rt6 * (-8.75D0*raz*raz*rBeta + 2.5D0*raz*czB + 1.25D0*rBeta)
          T1(3)        = 1.0_wp*R6 * rt6 * (-17.5D0*rax*raz*rBeta + 2.5D0*raz*cxB + 2.5D0*rax*czB)
          T1(4+ob)     = 1.0_wp*R6 * rt6 * (-8.75D0*rax*raz*raz + 1.25D0*rax)
          T1(7+(3*ob)) = 1.0_wp*R5 * rt6 * (1.25D0*raz*raz - 0.25D0)
          T1(9+(3*ob)) = 1.0_wp*R5 * rt6 * (2.5D0*rax*raz)

!
!                      31s, 10
! 
        Else If ((ia .EQ. 3) .AND. (ak .EQ. 2) .AND. (ib .EQ. 1) .AND. (bk .EQ. 0)) Then
          !T0           = 0.0_wp*R5 * rt6 * (8.75D0*ray*raz*raz*rBeta + 1.25D0*raz*raz*cyB &
          !                         & + 2.5D0*ray*raz*czB - 1.25D0*ray*rBeta - 0.25D0*cyB)
          !T1(0)        = 0.0_wp*R7 * rt6 * (-39.375D0*ray*raz*raz*rBeta - 4.375D0*raz*raz*cyB - 8.75D0*ray*raz*czB &
          !                         & + 4.375D0*ray*rBeta + 0.625D0*cyB)
          !T1(2)        = 0.0_wp*R6 * rt6 * (8.75D0*raz*raz*rBeta + 2.5D0*raz*czB - 1.25D0*rBeta)
          !T1(3)        = 0.0_wp*R6 * rt6 * (17.5D0*ray*raz*rBeta + 2.5D0*raz*cyB + 2.5D0*ray*czB)
          !T1(4+ob)     = 0.0_wp*R6 * rt6 * (8.75D0*ray*raz*raz - 1.25D0*ray)
          !T1(8+(3*ob)) = 0.0_wp*R5 * rt6 * (1.25D0*raz*raz - 0.25D0)
          !T1(9+(3*ob)) = 0.0_wp*R5 * rt6 * (2.5D0*ray*raz)
          T0           = 1.0_wp*R5 * rt6 * (-8.75D0*ray*raz*raz*rBeta + 1.25D0*raz*raz*cyB &
                                   & + 2.5D0*ray*raz*czB + 1.25D0*ray*rBeta - 0.25D0*cyB)
          T1(0)        = 1.0_wp*R7 * rt6 * (39.375D0*ray*raz*raz*rBeta - 4.375D0*raz*raz*cyB - 8.75D0*ray*raz*czB &
                                   & - 4.375D0*ray*rBeta + 0.625D0*cyB)
          T1(2)        = 1.0_wp*R6 * rt6 * (-8.75D0*raz*raz*rBeta + 2.5D0*raz*czB + 1.25D0*rBeta)
          T1(3)        = 1.0_wp*R6 * rt6 * (-17.5D0*ray*raz*rBeta + 2.5D0*raz*cyB + 2.5D0*ray*czB)
          T1(4+ob)     = 1.0_wp*R6 * rt6 * (-8.75D0*ray*raz*raz + 1.25D0*ray)
          T1(8+(3*ob)) = 1.0_wp*R5 * rt6 * (1.25D0*raz*raz - 0.25D0)
          T1(9+(3*ob)) = 1.0_wp*R5 * rt6 * (2.5D0*ray*raz)

!
!                      32c, 10
! 
        Else If ((ia .EQ. 3) .AND. (ak .EQ. 3) .AND. (ib .EQ. 1) .AND. (bk .EQ. 0)) Then
          !T0           = 0.0_wp*R5 * rt15 * (3.5d0*raz*rax*rax*rBeta - 3.5d0*raz*ray*ray*rBeta &
          !                          & + 0.5d0*rax*rax*czB - 0.5d0*ray*ray*czB &
          !                          & + raz*rax*cxB - raz*ray*cyB)
          !T1(0)        = 0.0_wp*R7 * rt15 * (-15.75d0*raz*rax*rax*rBeta + 15.75d0*raz*ray*ray*rBeta &
          !                          & -1.75d0*rax*rax*czB + 1.75d0*ray*ray*czB &
          !                          & -3.5d0*raz*rax*cxB + 3.5d0*raz*ray*cyB)
          !T1(3)        = 0.0_wp*R6 * rt15 * (3.5d0*rax*rax*rBeta - 3.5d0*ray*ray*rBeta &
          !                          & + rax*cxB - ray*cyB)
          !T1(1)        = 0.0_wp*R6 * rt15 * (7.0d0*raz*rax*rBeta + rax*czB + raz*cxB)
          !T1(2)        = 0.0_wp*R6 * rt15 * (-7.0d0*raz*ray*rBeta - ray*czB - raz*cyB)
          !T1(4+ob)     = 0.0_wp*R6 * rt15 * (3.5d0*raz*rax*rax - 3.5d0*raz*ray*ray) 
          !T1(7+(3*ob)) = 0.0_wp*R5 * rt15 * ( raz*rax)
          !T1(8+(3*ob)) = 0.0_wp*R5 * rt15 * (-raz*ray)
          !T1(9+(3*ob)) = 0.0_wp*R5 * rt15 * (0.5d0*rax*rax - 0.5d0*ray*ray)
          T0           = 1.0_wp*R5 * rt15 * (-3.5d0*raz*rax*rax*rBeta + 3.5d0*raz*ray*ray*rBeta &
                                    & + 0.5d0*rax*rax*czB - 0.5d0*ray*ray*czB &
                                    & + raz*rax*cxB - raz*ray*cyB)
          T1(0)        = 1.0_wp*R7 * rt15 * (15.75d0*raz*rax*rax*rBeta - 15.75d0*raz*ray*ray*rBeta &
                                    & -1.75d0*rax*rax*czB + 1.75d0*ray*ray*czB &
                                    & -3.5d0*raz*rax*cxB + 3.5d0*raz*ray*cyB)
          T1(3)        = 1.0_wp*R6 * rt15 * (-3.5d0*rax*rax*rBeta + 3.5d0*ray*ray*rBeta &
                                    & + rax*cxB - ray*cyB)
          T1(1)        = 1.0_wp*R6 * rt15 * (-7.0d0*raz*rax*rBeta + rax*czB + raz*cxB)
          T1(2)        = 1.0_wp*R6 * rt15 * (7.0d0*raz*ray*rBeta - ray*czB - raz*cyB)
          T1(4+ob)     = 1.0_wp*R6 * rt15 * (-3.5d0*raz*rax*rax + 3.5d0*raz*ray*ray) 
          T1(7+(3*ob)) = 1.0_wp*R5 * rt15 * ( raz*rax)
          T1(8+(3*ob)) = 1.0_wp*R5 * rt15 * (-raz*ray)
          T1(9+(3*ob)) = 1.0_wp*R5 * rt15 * (0.5d0*rax*rax - 0.5d0*ray*ray)

!
!                      32s, 10
! 
        Else If ((ia .EQ. 3) .AND. (ak .EQ. 4) .AND. (ib .EQ. 1) .AND. (bk .EQ. 0)) Then
          !T0           = 0.0_wp*R5 * rt15 * (7.0d0*rax*ray*raz*rBeta + rax*ray*czB + rax*raz*cyB + ray*raz*cxB)
          !T1(0)        = 0.0_wp*R7 * rt15 * (-31.5d0*rax*ray*raz*rBeta - 3.5d0*rax*ray*czB & 
          !                          & - 3.5d0*rax*raz*cyB - 3.5d0*ray*raz*cxB)
          !T1(1)        = 0.0_wp*R6 * rt15 * (7.0d0*ray*raz*rBeta + ray*czB + raz*cyB)
          !T1(2)        = 0.0_wp*R6 * rt15 * (7.0d0*rax*raz*rBeta + rax*czB + raz*cxB)
          !T1(3)        = 0.0_wp*R6 * rt15 * (7.0d0*rax*ray*rBeta + rax*cyB + ray*cxB)
          !T1(4+ob)     = 0.0_wp*R6 * rt15 * (7.0d0*rax*ray*raz)
          !T1(7+(3*ob)) = 0.0_wp*R5 * rt15 * (ray*raz)
          !T1(8+(3*ob)) = 0.0_wp*R5 * rt15 * (rax*raz)
          !T1(9+(3*ob)) = 0.0_wp*R5 * rt15 * (rax*ray)
          T0           = 1.0_wp*R5 * rt15 * (-7.0d0*rax*ray*raz*rBeta + rax*ray*czB + rax*raz*cyB + ray*raz*cxB)
          T1(0)        = 1.0_wp*R7 * rt15 * (31.5d0*rax*ray*raz*rBeta - 3.5d0*rax*ray*czB & 
                                    & - 3.5d0*rax*raz*cyB - 3.5d0*ray*raz*cxB)
          T1(1)        = 1.0_wp*R6 * rt15 * (-7.0d0*ray*raz*rBeta + ray*czB + raz*cyB)
          T1(2)        = 1.0_wp*R6 * rt15 * (-7.0d0*rax*raz*rBeta + rax*czB + raz*cxB)
          T1(3)        = 1.0_wp*R6 * rt15 * (-7.0d0*rax*ray*rBeta + rax*cyB + ray*cxB)
          T1(4+ob)     = 1.0_wp*R6 * rt15 * (-7.0d0*rax*ray*raz)
          T1(7+(3*ob)) = 1.0_wp*R5 * rt15 * (ray*raz)
          T1(8+(3*ob)) = 1.0_wp*R5 * rt15 * (rax*raz)
          T1(9+(3*ob)) = 1.0_wp*R5 * rt15 * (rax*ray)

!
!                      33c, 10
! 
        Else If ((ia .EQ. 3) .AND. (ak .EQ. 5) .AND. (ib .EQ. 1) .AND. (bk .EQ. 0)) Then
          !T0           = 0.0_wp*R5 * rt10 * (1.75d0*rax*rax*rax*rBeta + 0.75d0*rax*rax*cxB - 0.75d0*ray*ray*cxB &
          !                       & - 5.25d0*rax*ray*ray*rBeta - 1.5d0*rax*ray*cyB)
          !T1(0)        = 0.0_wp*R7 * rt10 * (-7.875d0*rax*rax*rax*rBeta - 2.625d0*rax*rax*cxB + 2.625d0*ray*ray*cxB &
          !                          & + 23.625*rax*ray*ray*rBeta + 5.25*rax*ray*cyB)
          !T1(1)        = 0.0_wp*R6 * rt10 * (5.25d0*rax*rax*rBeta + 1.5d0*rax*cxB - 5.25d0*ray*ray*rBeta - 1.5d0*ray*cyB)
          !T1(2)        = 0.0_wp*R6 * rt10 * ((-1.5d0)*ray*cxB - 10.5d0*rax*ray*rBeta - 1.5d0*rax*cyB)
          !T1(4+ob)     = 0.0_wp*R6 * rt10 * (1.75d0*rax*rax*rax - 5.25d0*rax*ray*ray)
          !T1(7+(3*ob)) = 0.0_wp*R5 * rt10 * (0.75d0*rax*rax - 0.75d0*ray*ray)
          !T1(8+(3*ob)) = 0.0_wp*R5 * rt10 * ((-1.5d0)*rax*ray)
          T0           = 1.0_wp*R5 * rt10 * (-1.75d0*rax*rax*rax*rBeta + 0.75d0*rax*rax*cxB - 0.75d0*ray*ray*cxB &
                                 & + 5.25d0*rax*ray*ray*rBeta - 1.5d0*rax*ray*cyB)
          T1(0)        = 1.0_wp*R7 * rt10 * (7.875d0*rax*rax*rax*rBeta - 2.625d0*rax*rax*cxB + 2.625d0*ray*ray*cxB &
                                    & - 23.625*rax*ray*ray*rBeta + 5.25*rax*ray*cyB)
          T1(1)        = 1.0_wp*R6 * rt10 * (-5.25d0*rax*rax*rBeta + 1.5d0*rax*cxB + 5.25d0*ray*ray*rBeta - 1.5d0*ray*cyB)
          T1(2)        = 1.0_wp*R6 * rt10 * ((-1.5d0)*ray*cxB + 10.5d0*rax*ray*rBeta - 1.5d0*rax*cyB)
          T1(4+ob)     = 1.0_wp*R6 * rt10 * (-1.75d0*rax*rax*rax + 5.25d0*rax*ray*ray)
          T1(7+(3*ob)) = 1.0_wp*R5 * rt10 * (0.75d0*rax*rax - 0.75d0*ray*ray)
          T1(8+(3*ob)) = 1.0_wp*R5 * rt10 * ((-1.5d0)*rax*ray)

!
!                      33s, 10
!
        Else If ((ia .EQ. 3) .AND. (ak .EQ. 6) .AND. (ib .EQ. 1) .AND. (bk .EQ. 0)) Then
          !T0           = 0.0_wp*R5 * rt10 * ((-1.75d0)*ray*ray*ray*rBeta + 0.75d0*rax*rax*cyB - 0.75d0*ray*ray*cyB &
          !                          & + 5.25d0*rax*rax*ray*rBeta + 1.5d0*rax*ray*cxB)
          !T1(0)        = 0.0_wp*R7 * rt10 * (7.875d0*ray*ray*ray*rBeta - 2.625d0*rax*rax*cyB + 2.625d0*ray*ray*cyB &
          !                          & - 23.625d0*rax*rax*ray*rBeta - 5.25d0*rax*ray*cxB)
          !T1(1)        = 0.0_wp*R6 * rt10 * (1.5d0*rax*cyB + 10.5d0*rax*ray*rBeta + 1.5d0*ray*cxB)
          !T1(2)        = 0.0_wp*R6 * rt10 * ((-5.25d0)*ray*ray*rBeta - 1.5d0*ray*cyB + 5.25d0*rax*rax*rBeta + 1.5d0*rax*cxB)
          !T1(4+ob)     = 0.0_wp*R6 * rt10 * ((-1.75d0)*ray*ray*ray + 5.25d0*rax*rax*ray)
          !T1(7+(3*ob)) = 0.0_wp*R5 * rt10 * (1.5d0*rax*ray)
          !T1(8+(3*ob)) = 0.0_wp*R5 * rt10 * (0.75d0*rax*rax - 0.75d0*ray*ray)
          T0           = 1.0_wp*R5 * rt10 * ((1.75d0)*ray*ray*ray*rBeta + 0.75d0*rax*rax*cyB - 0.75d0*ray*ray*cyB &
                                    & - 5.25d0*rax*rax*ray*rBeta + 1.5d0*rax*ray*cxB)
          T1(0)        = 1.0_wp*R7 * rt10 * (-7.875d0*ray*ray*ray*rBeta - 2.625d0*rax*rax*cyB + 2.625d0*ray*ray*cyB &
                                    & + 23.625d0*rax*rax*ray*rBeta - 5.25d0*rax*ray*cxB)
          T1(1)        = 1.0_wp*R6 * rt10 * (1.5d0*rax*cyB - 10.5d0*rax*ray*rBeta + 1.5d0*ray*cxB)
          T1(2)        = 1.0_wp*R6 * rt10 * ((5.25d0)*ray*ray*rBeta - 1.5d0*ray*cyB - 5.25d0*rax*rax*rBeta + 1.5d0*rax*cxB)
          T1(4+ob)     = 1.0_wp*R6 * rt10 * ((1.75d0)*ray*ray*ray - 5.25d0*rax*rax*ray)
          T1(7+(3*ob)) = 1.0_wp*R5 * rt10 * (1.5d0*rax*ray)
          T1(8+(3*ob)) = 1.0_wp*R5 * rt10 * (0.75d0*rax*rax - 0.75d0*ray*ray)

!
!  Dipole-Octopole     10, 30
!
        Else If ((ia .EQ. 1) .AND. (ak .EQ. 0) .AND. (ib .EQ. 3) .AND. (bk .EQ. 0)) Then
          !T0           = 0.0_wp*R5 * (17.5D0*rbz*rbz*rbz*rAlpha + 7.5D0*rbz*rbz*czA - 7.5D0*rbz*rAlpha - 1.5D0*czA)
          !T1(0)        = 0.0_wp*R7 * (-78.75D0*rbz*rbz*rbz*rAlpha - 26.25D0*rbz*rbz*czA + 26.25D0*rbz*rAlpha + 3.75*czA)
          !T1(6)        = 0.0_wp*R6 * (52.5D0*rbz*rbz*rAlpha + 15.0D0*rbz*czA - 7.5D0*rAlpha)
          !T1(1+oa)     = 0.0_wp*R6 * (17.5D0*rbz*rbz*rbz - 7.5D0*rbz)
          !T1(9+(3*oa)) = 0.0_wp*R5 * (7.5D0*rbz*rbz - 1.5D0)
          T0           = 1.0_wp*R5 * (-17.5D0*rbz*rbz*rbz*rAlpha + 7.5D0*rbz*rbz*czA + 7.5D0*rbz*rAlpha - 1.5D0*czA)
          T1(0)        = 1.0_wp*R7 * (78.75D0*rbz*rbz*rbz*rAlpha - 26.25D0*rbz*rbz*czA - 26.25D0*rbz*rAlpha + 3.75*czA)
          T1(6)        = 1.0_wp*R6 * (-52.5D0*rbz*rbz*rAlpha + 15.0D0*rbz*czA + 7.5D0*rAlpha)
          T1(1+oa)     = 1.0_wp*R6 * (-17.5D0*rbz*rbz*rbz + 7.5D0*rbz)
          T1(9+(3*oa)) = 1.0_wp*R5 * (7.5D0*rbz*rbz - 1.5D0)

!
!                      10, 31c
!
        Else If ((ia .EQ. 1) .AND. (ak .EQ. 0) .AND. (ib .EQ. 3) .AND. (bk .EQ. 1)) Then
          !T0           = 0.0_wp*R5 * rt6 * (8.75D0*rbx*rbz*rbz*rAlpha + 1.25D0*rbz*rbz*cxA &
          !                      & + 2.5D0*rbx*rbz*czA - 1.25D0*rbx*rAlpha - 0.25D0*cxA)
          !T1(0)        = 0.0_wp*R7 * rt6 * (-39.375D0*rbx*rbz*rbz*rAlpha - 4.375D0*rbz*rbz*cxA &
          !                         & - 8.75D0*rbx*rbz*czA + 4.375D0*rbx*rAlpha + 0.625D0*cxA)
          !T1(4)        = 0.0_wp*R6 * rt6 * (8.75D0*rbz*rbz*rAlpha + 2.5D0*rbz*czA - 1.25D0*rAlpha)
          !T1(6)        = 0.0_wp*R6 * rt6 * (17.5D0*rbx*rbz*rAlpha + 2.5D0*rbz*cxA + 2.5D0*rbx*czA)
          !T1(1+oa)     = 0.0_wp*R6 * rt6 * (8.75D0*rbx*rbz*rbz - 1.25D0*rbx)
          !T1(7+(3*oa)) = 0.0_wp*R5 * rt6 * (1.25D0*rbz*rbz - 0.25D0)
          !T1(9+(3*oa)) = 0.0_wp*R5 * rt6 * (2.5D0*rbx*rbz)
          T0           = 1.0_wp*R5 * rt6 * (-8.75D0*rbx*rbz*rbz*rAlpha + 1.25D0*rbz*rbz*cxA &
                                & + 2.5D0*rbx*rbz*czA + 1.25D0*rbx*rAlpha - 0.25D0*cxA)
          T1(0)        = 1.0_wp*R7 * rt6 * (39.375D0*rbx*rbz*rbz*rAlpha - 4.375D0*rbz*rbz*cxA &
                                   & - 8.75D0*rbx*rbz*czA - 4.375D0*rbx*rAlpha + 0.625D0*cxA)
          T1(4)        = 1.0_wp*R6 * rt6 * (-8.75D0*rbz*rbz*rAlpha + 2.5D0*rbz*czA + 1.25D0*rAlpha)
          T1(6)        = 1.0_wp*R6 * rt6 * (-17.5D0*rbx*rbz*rAlpha + 2.5D0*rbz*cxA + 2.5D0*rbx*czA)
          T1(1+oa)     = 1.0_wp*R6 * rt6 * (-8.75D0*rbx*rbz*rbz + 1.25D0*rbx)
          T1(7+(3*oa)) = 1.0_wp*R5 * rt6 * (1.25D0*rbz*rbz - 0.25D0)
          T1(9+(3*oa)) = 1.0_wp*R5 * rt6 * (2.5D0*rbx*rbz)

!
!                      10, 31s
!
        Else If ((ia .EQ. 1) .AND. (ak .EQ. 0) .AND. (ib .EQ. 3) .AND. (bk .EQ. 2)) Then
          !T0           = 0.0_wp*R5 * rt6 * (8.75D0*rby*rbz*rbz*rAlpha + 1.25D0*rbz*rbz*cyA &
          !                         & + 2.5D0*rby*rbz*czA - 1.25D0*rby*rAlpha - 0.25D0*cyA)
          !T1(0)        = 0.0_wp*R7 * rt6 * (-39.375D0*rby*rbz*rbz*rAlpha - 4.375D0*rbz*rbz*cyA - 8.75D0*rby*rbz*czA &
          !                         & + 4.375D0*rby*rAlpha +0.625D0*cyA)
          !T1(5)        = 0.0_wp*R6 * rt6 * (8.75D0*rbz*rbz*rAlpha + 2.5D0*rbz*czA - 1.25D0*rAlpha)
          !T1(6)        = 0.0_wp*R6 * rt6 * (17.5D0*rby*rbz*rAlpha + 2.5D0*rbz*cyA + 2.5D0*rby*czA)
          !T1(1+oa)     = 0.0_wp*R6 * rt6 * (8.75D0*rby*rbz*rbz - 1.25D0*rby)
          !T1(8+(3*oa))    = 0.0_wp*R5 * rt6 * (1.25D0*rbz*rbz - 0.25D0)
          !T1(9+(3*oa))    = 0.0_wp*R5 * rt6 * (2.5D0*rby*rbz)
          T0           = 1.0_wp*R5 * rt6 * (-8.75D0*rby*rbz*rbz*rAlpha + 1.25D0*rbz*rbz*cyA &
                                   & + 2.5D0*rby*rbz*czA + 1.25D0*rby*rAlpha - 0.25D0*cyA)
          T1(0)        = 1.0_wp*R7 * rt6 * (39.375D0*rby*rbz*rbz*rAlpha - 4.375D0*rbz*rbz*cyA - 8.75D0*rby*rbz*czA &
                                   & - 4.375D0*rby*rAlpha +0.625D0*cyA)
          T1(5)        = 1.0_wp*R6 * rt6 * (-8.75D0*rbz*rbz*rAlpha + 2.5D0*rbz*czA + 1.25D0*rAlpha)
          T1(6)        = 1.0_wp*R6 * rt6 * (-17.5D0*rby*rbz*rAlpha + 2.5D0*rbz*cyA + 2.5D0*rby*czA)
          T1(1+oa)     = 1.0_wp*R6 * rt6 * (-8.75D0*rby*rbz*rbz + 1.25D0*rby)
          T1(8+(3*oa))    = 1.0_wp*R5 * rt6 * (1.25D0*rbz*rbz - 0.25D0)
          T1(9+(3*oa))    = 1.0_wp*R5 * rt6 * (2.5D0*rby*rbz)
!
!                      10, 32c
! JCRT
        Else If ((ia .EQ. 1) .AND. (ak .EQ. 0) .AND. (ib .EQ. 3) .AND. (bk .EQ. 3)) Then
          !T0           = 0.0_wp*R5 * rt15 * (3.5d0*rbz*rbx*rbx*rAlpha - 3.5d0*rbz*rby*rby*rAlpha &
          !                          & + 0.5d0*rbx*rbx*czA - 0.5d0*rby*rby*czA &
          !                          & + rbz*rbx*cxA - rbz*rby*cyA)
          !T1(0)        = 0.0_wp*R7 * rt15 * (-15.75d0*rbz*rbx*rbx*rAlpha + 15.75d0*rbz*rby*rby*rAlpha &
          !                          & -1.75d0*rbx*rbx*czA + 1.75d0*rby*rby*czA &
          !                          & -3.5d0*rbz*rbx*cxA + 3.5D0*rbz*rby*cyA)
          !T1(6)        = 0.0_wp*R6 * rt15 * (3.5d0*rbx*rbx*rAlpha - 3.5d0*rby*rby*rAlpha &
          !                          & + rbx*cxA - rby*cyA)
          !T1(4)        = 0.0_wp*R6 * rt15 * (7.0d0*rbz*rbx*rAlpha + rbx*czA + rbz*cxA)
          !T1(5)        = 0.0_wp*R6 * rt15 * (-7.0d0*rbz*rby*rAlpha - rby*czA - rbz*cyA)
          !T1(1+oa)     = 0.0_wp*R6 * rt15 * (3.5d0*rbz*rbx*rbx - 3.5d0*rbz*rby*rby)
          !T1(7+(3*oa))     = 0.0_wp*R5 * rt15 * ( rbz*rbx)
          !T1(8+(3*oa))    = 0.0_wp*R5 * rt15 * (-rbz*rby)
          !T1(9+(3*oa))    = 0.0_wp*R5 * rt15 * (0.5d0*rbx*rbx - 0.5d0*rby*rby)
          T0           = 1.0_wp*R5 * rt15 * (-3.5d0*rbz*rbx*rbx*rAlpha + 3.5d0*rbz*rby*rby*rAlpha &
                                    & + 0.5d0*rbx*rbx*czA - 0.5d0*rby*rby*czA &
                                    & + rbz*rbx*cxA - rbz*rby*cyA)
          T1(0)        = 1.0_wp*R7 * rt15 * (15.75d0*rbz*rbx*rbx*rAlpha - 15.75d0*rbz*rby*rby*rAlpha &
                                    & -1.75d0*rbx*rbx*czA + 1.75d0*rby*rby*czA &
                                    & -3.5d0*rbz*rbx*cxA + 3.5D0*rbz*rby*cyA)
          T1(6)        = 1.0_wp*R6 * rt15 * (-3.5d0*rbx*rbx*rAlpha + 3.5d0*rby*rby*rAlpha &
                                    & + rbx*cxA - rby*cyA)
          T1(4)        = 1.0_wp*R6 * rt15 * (-7.0d0*rbz*rbx*rAlpha + rbx*czA + rbz*cxA)
          T1(5)        = 1.0_wp*R6 * rt15 * (7.0d0*rbz*rby*rAlpha - rby*czA - rbz*cyA)
          T1(1+oa)     = 1.0_wp*R6 * rt15 * (-3.5d0*rbz*rbx*rbx + 3.5d0*rbz*rby*rby)
          T1(7+(3*oa))     = 1.0_wp*R5 * rt15 * ( rbz*rbx)
          T1(8+(3*oa))    = 1.0_wp*R5 * rt15 * (-rbz*rby)
          T1(9+(3*oa))    = 1.0_wp*R5 * rt15 * (0.5d0*rbx*rbx - 0.5d0*rby*rby)

!
!                      10, 32s
!
        Else If ((ia .EQ. 1) .AND. (ak .EQ. 0) .AND. (ib .EQ. 3) .AND. (bk .EQ. 4)) Then
          !T0           = 0.0_wp*R5 * rt15 * (7.0d0*rbx*rby*rbz*rAlpha + rbx*rby*czA + rbx*rbz*cyA + rby*rbz*cxA)
          !T1(0)        = 0.0_wp*R7 * rt15 * (-31.5d0*rbx*rby*rbz*rAlpha - 3.5d0*rbx*rby*czA &
          !                          & - 3.5d0*rbx*rbz*cyA - 3.5d0*rby*rbz*cxA)
          !T1(4)        = 0.0_wp*R6 * rt15 * (7.0d0*rby*rbz*rAlpha + rby*czA + rbz*cyA)
          !T1(5)        = 0.0_wp*R6 * rt15 * (7.0d0*rbx*rbz*rAlpha + rbx*czA + rbz*cxA)
          !T1(6)        = 0.0_wp*R6 * rt15 * (7.0d0*rbx*rby*rAlpha + rbx*cyA + rby*cxA)
          !T1(1+oa)     = 0.0_wp*R6 * rt15 * (7.0d0*rbx*rby*rbz)
          !T1(7+(3*oa))     = 0.0_wp*R5 * rt15 * (rby*rbz)
          !T1(8+(3*oa))    = 0.0_wp*R5 * rt15 * (rbx*rbz)
          !T1(9+(3*oa))    = 0.0_wp*R5 * rt15 * (rbx*rby)
          T0           = 1.0_wp*R5 * rt15 * (-7.0d0*rbx*rby*rbz*rAlpha + rbx*rby*czA + rbx*rbz*cyA + rby*rbz*cxA)
          T1(0)        = 1.0_wp*R7 * rt15 * (31.5d0*rbx*rby*rbz*rAlpha - 3.5d0*rbx*rby*czA &
                                    & - 3.5d0*rbx*rbz*cyA - 3.5d0*rby*rbz*cxA)
          T1(4)        = 1.0_wp*R6 * rt15 * (-7.0d0*rby*rbz*rAlpha + rby*czA + rbz*cyA)
          T1(5)        = 1.0_wp*R6 * rt15 * (-7.0d0*rbx*rbz*rAlpha + rbx*czA + rbz*cxA)
          T1(6)        = 1.0_wp*R6 * rt15 * (-7.0d0*rbx*rby*rAlpha + rbx*cyA + rby*cxA)
          T1(1+oa)     = 1.0_wp*R6 * rt15 * (-7.0d0*rbx*rby*rbz)
          T1(7+(3*oa))     = 1.0_wp*R5 * rt15 * (rby*rbz)
          T1(8+(3*oa))    = 1.0_wp*R5 * rt15 * (rbx*rbz)
          T1(9+(3*oa))    = 1.0_wp*R5 * rt15 * (rbx*rby)

!
!                      10, 33c
!
        Else If ((ia .EQ. 1) .AND. (ak .EQ. 0) .AND. (ib .EQ. 3) .AND. (bk .EQ. 5)) Then
          !T0           = 0.0_wp*R5 * rt10 * (1.75d0*rbx*rbx*rbx*rAlpha + 0.75d0*rbx*rbx*cxA - 0.75d0*rby*rby*cxA &
          !                       & - 5.25d0*rbx*rby*rby*rAlpha - 1.5d0*rbx*rby*cyA)
          !T1(0)        = 0.0_wp*R7 * rt10 * (-7.875d0*rbx*rbx*rbx*rAlpha - 2.625d0*rbx*rbx*cxA + 2.625d0*rby*rby*cxA &
          !                          & + 23.625*rbx*rby*rby*rAlpha + 5.25*rbx*rby*cyA)
          !T1(4)        = 0.0_wp*R6 * rt10 * (5.25d0*rbx*rbx*rAlpha + 1.5d0*rbx*cxA - 5.25d0*rby*rby*rAlpha - 1.5d0*rby*cyA)
          !T1(5)        = 0.0_wp*R6 * rt10 * ((-1.5d0)*rby*cxA - 10.5d0*rbx*rby*rAlpha - 1.5d0*rbx*cyA)
          !T1(1+oa)     = 0.0_wp*R6 * rt10 * (1.75d0*rbx*rbx*rbx - 5.25d0*rbx*rby*rby)
          !T1(7+(3*oa))     = 0.0_wp*R5 * rt10 * (0.75d0*rbx*rbx - 0.75d0*rby*rby)
          !T1(8+(3*oa))    = 0.0_wp*R5 * rt10 * ((-1.5d0)*rbx*rby)
          T0           = 1.0_wp*R5 * rt10 * ((-1.75_wp)*rbx*rbx*rbx*rAlpha + 0.75d0*rbx*rbx*cxA - 0.75d0*rby*rby*cxA &
                                 & + 5.25d0*rbx*rby*rby*rAlpha - 1.5d0*rbx*rby*cyA)
          T1(0)        = 1.0_wp*R7 * rt10 * (7.875d0*rbx*rbx*rbx*rAlpha - 2.625d0*rbx*rbx*cxA + 2.625d0*rby*rby*cxA &
                                    & - 23.625*rbx*rby*rby*rAlpha + 5.25*rbx*rby*cyA)
          T1(4)        = 1.0_wp*R6 * rt10 * (-5.25d0*rbx*rbx*rAlpha + 1.5d0*rbx*cxA + 5.25d0*rby*rby*rAlpha - 1.5d0*rby*cyA)
          T1(5)        = 1.0_wp*R6 * rt10 * ((-1.5d0)*rby*cxA + 10.5d0*rbx*rby*rAlpha - 1.5d0*rbx*cyA)
          T1(1+oa)     = 1.0_wp*R6 * rt10 * (-1.75d0*rbx*rbx*rbx + 5.25d0*rbx*rby*rby)
          T1(7+(3*oa))     = 1.0_wp*R5 * rt10 * (0.75d0*rbx*rbx - 0.75d0*rby*rby)
          T1(8+(3*oa))    = 1.0_wp*R5 * rt10 * ((-1.5d0)*rbx*rby)

!
!                      10, 33s
!
        Else If ((ia .EQ. 1) .AND. (ak .EQ. 0) .AND. (ib .EQ. 3) .AND. (bk .EQ. 6)) Then
          !T0        = 0.0_wp*R5 * rt10 * ((-1.75d0)*rby*rby*rby*rAlpha + 0.75d0*rbx*rbx*cyA - 0.75d0*rby*rby*cyA &
          !                       & + 5.25d0*rbx*rbx*rby*rAlpha + 1.5d0*rbx*rby*cxA)
          !T1(0)     = 0.0_wp*R7 * rt10 * (7.875d0*rby*rby*rby*rAlpha - 2.625d0*rbx*rbx*cyA + 2.625d0*rby*rby*cyA &
          !                       & - 23.625d0*rbx*rbx*rby*rAlpha - 5.25d0*rbx*rby*cxA)
          !T1(4)     = 0.0_wp*R6 * rt10 * (1.5d0*rbx*cyA + 10.5d0*rbx*rby*rAlpha + 1.5d0*rby*cxA)
          !T1(5)     = 0.0_wp*R6 * rt10 * ((-5.25d0)*rby*rby*rAlpha - 1.5d0*rby*cyA &
          !                       & + 5.25d0*rbx*rbx*rAlpha + 1.5d0*rbx*cxA)
          !T1(1+oa)  = 0.0_wp*R6 * rt10 * ((-1.75d0)*rby*rby*rby + 5.25d0*rbx*rbx*rby)
          !T1(7+(3*oa))  = 0.0_wp*R5 * rt10 * (1.5d0*rbx*rby)
          !T1(8+(3*oa)) = 0.0_wp*R5 * rt10 * (0.75d0*rbx*rbx - 0.75d0*rby*rby)
          T0        = 1.0_wp*R5 * rt10 * ((1.75d0)*rby*rby*rby*rAlpha + 0.75d0*rbx*rbx*cyA - 0.75d0*rby*rby*cyA &
                                 & - 5.25d0*rbx*rbx*rby*rAlpha + 1.5d0*rbx*rby*cxA)
          T1(0)     = 1.0_wp*R7 * rt10 * (-7.875d0*rby*rby*rby*rAlpha - 2.625d0*rbx*rbx*cyA + 2.625d0*rby*rby*cyA &
                                 & + 23.625d0*rbx*rbx*rby*rAlpha - 5.25d0*rbx*rby*cxA)
          T1(4)     = 1.0_wp*R6 * rt10 * (1.5d0*rbx*cyA - 10.5d0*rbx*rby*rAlpha + 1.5d0*rby*cxA)
          T1(5)     = 1.0_wp*R6 * rt10 * ((5.25d0)*rby*rby*rAlpha - 1.5d0*rby*cyA &
                                 & - 5.25d0*rbx*rbx*rAlpha + 1.5d0*rbx*cxA)
          T1(1+oa)  = 1.0_wp*R6 * rt10 * ((1.75d0)*rby*rby*rby - 5.25d0*rbx*rbx*rby)
          T1(7+(3*oa))  = 1.0_wp*R5 * rt10 * (1.5d0*rbx*rby)
          T1(8+(3*oa)) = 1.0_wp*R5 * rt10 * (0.75d0*rbx*rbx - 0.75d0*rby*rby)

!
!  Quadrupole-quadrupole 20, 20
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 0)) Then
          !T0        = R5 * ((105.0D0/4.0D0)*(raz**(2.0D0))*(rbz**(2.0D0))-(15.0D0/4.0D0)*(raz**(2.0D0))-(15.0D0/4.0D0)*(rbz**(2.0D0))+15.0D0*raz*rbz*Czz  &
          !                                    & +(3.0D0/2.0D0)*(Czz**(2.0D0))+(3.0D0/4.0D0))
          !T1(0)     = R7 * ((-945.0D0/8.0D0)*(raz**(2.0D0))*(rbz**(2.0D0))+(105.0D0/8.0D0)*(raz**(2.0D0))+(105.0D0/8.0D0)*(rbz**(2.0D0))-(105.0D0/2.0D0)*raz*rbz*Czz  &
          !                                    & -(15.0D0/4.0D0)*(Czz**(2.0D0))-(15.0D0/8.0D0))
          !T1(3)     = R6 * ((105.0D0/2.0D0)*raz*(rbz**(2.0D0))-(15.0D0/2.0D0)*raz+15.0D0*rbz*Czz)
          !T1(6)     = R6 * ((105.0D0/2.0D0)*(raz**(2.0D0))*rbz-(15.0D0/2.0D0)*rbz+15.0D0*raz*Czz)
          !T1(15)    = R5 * (15.0D0*raz*rbz+3.0D0*Czz)
          T0        = R5 * ((105.0D0/4.0D0)*(raz**(2.0D0))*(rbz**(2.0D0))-(15.0D0/4.0D0)*(raz**(2.0D0))-(15.0D0/4.0D0)*(rbz**(2.0D0))-15.0D0*raz*rbz*Czz  &
                                              & +(3.0D0/2.0D0)*(Czz**(2.0D0))+(3.0D0/4.0D0))
          T1(0)     = R7 * ((-945.0D0/8.0D0)*(raz**(2.0D0))*(rbz**(2.0D0))+(105.0D0/8.0D0)*(raz**(2.0D0))+(105.0D0/8.0D0)*(rbz**(2.0D0))+(105.0D0/2.0D0)*raz*rbz*Czz  &
                                              & -(15.0D0/4.0D0)*(Czz**(2.0D0))-(15.0D0/8.0D0))
          T1(3)     = R6 * ((105.0D0/2.0D0)*raz*(rbz**(2.0D0))-(15.0D0/2.0D0)*raz-15.0D0*rbz*Czz)
          T1(6)     = R6 * ((105.0D0/2.0D0)*(raz**(2.0D0))*rbz-(15.0D0/2.0D0)*rbz-15.0D0*raz*Czz)
          T1(15)    = R5 * ((-15.0D0)*raz*rbz+3.0D0*Czz)
!
!                      20, 21c
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 1)) Then
          !T0        = R5 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rbx*rbz-(5.0D0/2.0D0)*rbx*rbz+5.0D0*raz*rbx*Czz+5.0D0*raz*rbz*Czx+Czx*Czz))
          !T1(0)     = R7 * (rt3*((-315.0D0/4.0D0)*(raz**(2.0D0))*rbx*rbz+(35.0D0/4.0D0)*rbx*rbz-(35.0D0/2.0D0)*raz*rbx*Czz-(35.0D0/2.0D0)*raz*rbz*Czx &
          !                                    & -(5.0D0/2.0D0)*Czx*Czz))
          !T1(3)     = R6 * (rt3*(35.0D0*raz*rbx*rbz+5.0D0*rbx*Czz+5.0D0*rbz*Czx))
          !T1(4)     = R6 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rbz-(5.0D0/2.0D0)*rbz+5.0D0*raz*Czz))
          !T1(6)     = R6 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rbx-(5.0D0/2.0D0)*rbx+5.0D0*raz*Czx))
          !T1(9)     = R5 * (rt3*(5.0D0*raz*rbz+Czz))
          !T1(15)    = R5 * (rt3*(5.0D0*raz*rbx+Czx))
          T0        = R5 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rbx*rbz-(5.0D0/2.0D0)*rbx*rbz-5.0D0*raz*rbx*Czz-5.0D0*raz*rbz*Czx+Czx*Czz))
          T1(0)     = R7 * (rt3*((-315.0D0/4.0D0)*(raz**(2.0D0))*rbx*rbz+(35.0D0/4.0D0)*rbx*rbz+(35.0D0/2.0D0)*raz*rbx*Czz+(35.0D0/2.0D0)*raz*rbz*Czx &
                                              & -(5.0D0/2.0D0)*Czx*Czz))
          T1(3)     = R6 * (rt3*(35.0D0*raz*rbx*rbz-5.0D0*rbx*Czz-5.0D0*rbz*Czx))
          T1(4)     = R6 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rbz-(5.0D0/2.0D0)*rbz-5.0D0*raz*Czz))
          T1(6)     = R6 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rbx-(5.0D0/2.0D0)*rbx-5.0D0*raz*Czx))
          T1(9)     = R5 * (rt3*((-5.0D0)*raz*rbz+Czz))
          T1(15)    = R5 * (rt3*((-5.0D0)*raz*rbx+Czx))
!
!                      20, 21s
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 2)) Then
          !T0        = R5 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rby*rbz-(5.0D0/2.0D0)*rby*rbz+5.0D0*raz*rby*Czz+5.0D0*raz*rbz*Czy+Czy*Czz))
          !T1(0)     = R7 * (rt3*((-315.0D0/4.0D0)*(raz**(2.0D0))*rby*rbz+(35.0D0/4.0D0)*rby*rbz-(35.0D0/2.0D0)*raz*rby*Czz-(35.0D0/2.0D0)*raz*rbz*Czy &
          !                                    & -(5.0D0/2.0D0)*Czy*Czz))
          !T1(3)     = R6 * (rt3*(35.0D0*raz*rby*rbz+5.0D0*rby*Czz+5.0D0*rbz*Czy))
          !T1(5)     = R6 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rbz-(5.0D0/2.0D0)*rbz+5.0D0*raz*Czz))
          !T1(6)     = R6 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rby-(5.0D0/2.0D0)*rby+5.0D0*raz*Czy))
          !T1(12)    = R5 * (rt3*(5.0D0*raz*rbz+Czz))
          !T1(15)    = R5 * (rt3*(5.0D0*raz*rby+Czy))
          T0        = R5 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rby*rbz-(5.0D0/2.0D0)*rby*rbz-5.0D0*raz*rby*Czz-5.0D0*raz*rbz*Czy+Czy*Czz))
          T1(0)     = R7 * (rt3*((-315.0D0/4.0D0)*(raz**(2.0D0))*rby*rbz+(35.0D0/4.0D0)*rby*rbz+(35.0D0/2.0D0)*raz*rby*Czz+(35.0D0/2.0D0)*raz*rbz*Czy &
                                              & -(5.0D0/2.0D0)*Czy*Czz))
          T1(3)     = R6 * (rt3*(35.0D0*raz*rby*rbz-5.0D0*rby*Czz-5.0D0*rbz*Czy))
          T1(5)     = R6 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rbz-(5.0D0/2.0D0)*rbz-5.0D0*raz*Czz))
          T1(6)     = R6 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rby-(5.0D0/2.0D0)*rby-5.0D0*raz*Czy))
          T1(12)    = R5 * (rt3*((-5.0D0)*raz*rbz+Czz))
          T1(15)    = R5 * (rt3*((-5.0D0)*raz*rby+Czy))
!
!                      20, 22c
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 3)) Then
          !T0        = R5 * (rt3*((35.0D0/4.0D0)*(raz**(2.0D0))*(rbx**(2.0D0))-(35.0D0/4.0D0)*(raz**(2.0D0))*(rby**(2.0D0))-(5.0D0/4.0D0)*(rbx**(2.0D0))+(5.0D0/4.0D0)*(rby**(2.0D0)) &
          !                                    & +5.0D0*raz*rbx*Czx-5.0D0*raz*rby*Czy+(1.0D0/2.0D0)*(Czx**(2.0D0))-(1.0D0/2.0D0)*(Czy**(2.0D0))))
          !T1(0)     = R7 * (rt3*((-315.0D0/8.0D0)*(raz**(2.0D0))*(rbx**(2.0D0))+(315.0D0/8.0D0)*(raz**(2.0D0))*(rby**(2.0D0))+(35.0D0/8.0D0)*(rbx**(2.0D0))-(35.0D0/8.0D0)*(rby**(2.0D0)) &
          !                                    & -(35.0D0/2.0D0)*raz*rbx*Czx+(35.0D0/2.0D0)*raz*rby*Czy-(5.0D0/4.0D0)*(Czx**(2.0D0))+(5.0D0/4.0D0)*(Czy**(2.0D0))))
          !T1(3)     = R6 * (rt3*((35.0D0/2.0D0)*raz*(rbx**(2.0D0))-(35.0D0/2.0D0)*raz*(rby**(2.0D0))+5.0D0*rbx*Czx-5.0D0*rby*Czy))
          !T1(4)     = R6 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rbx-(5.0D0/2.0D0)*rbx+5.0D0*raz*Czx))
          !T1(5)     = R6 * (rt3*((-35.0D0/2.0D0)*(raz**(2.0D0))*rby+(5.0D0/2.0D0)*rby-5.0D0*raz*Czy))
          !T1(9)     = R5 * (rt3*(5.0D0*raz*rbx+Czx))
          !T1(12)    = R5 * (rt3*(-5.0D0*raz*rby-Czy))
          T0        = R5 * (rt3*((35.0D0/4.0D0)*(raz**(2.0D0))*(rbx**(2.0D0))-(35.0D0/4.0D0)*(raz**(2.0D0))*(rby**(2.0D0))-(5.0D0/4.0D0)*(rbx**(2.0D0))+(5.0D0/4.0D0)*(rby**(2.0D0)) &
                                              & -5.0D0*raz*rbx*Czx+5.0D0*raz*rby*Czy+(1.0D0/2.0D0)*(Czx**(2.0D0))-(1.0D0/2.0D0)*(Czy**(2.0D0))))
          T1(0)     = R7 * (rt3*((-315.0D0/8.0D0)*(raz**(2.0D0))*(rbx**(2.0D0))+(315.0D0/8.0D0)*(raz**(2.0D0))*(rby**(2.0D0))+(35.0D0/8.0D0)*(rbx**(2.0D0))-(35.0D0/8.0D0)*(rby**(2.0D0)) &
                                              & +(35.0D0/2.0D0)*raz*rbx*Czx-(35.0D0/2.0D0)*raz*rby*Czy-(5.0D0/4.0D0)*(Czx**(2.0D0))+(5.0D0/4.0D0)*(Czy**(2.0D0))))
          T1(3)     = R6 * (rt3*((35.0D0/2.0D0)*raz*(rbx**(2.0D0))-(35.0D0/2.0D0)*raz*(rby**(2.0D0))-5.0D0*rbx*Czx+5.0D0*rby*Czy))
          T1(4)     = R6 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rbx-(5.0D0/2.0D0)*rbx-5.0D0*raz*Czx))
          T1(5)     = R6 * (rt3*((-35.0D0/2.0D0)*(raz**(2.0D0))*rby+(5.0D0/2.0D0)*rby+5.0D0*raz*Czy))
          T1(9)     = R5 * (rt3*((-5.0D0)*raz*rbx+Czx))
          T1(12)    = R5 * (rt3*(5.0D0*raz*rby-Czy))
!
!                      20, 22s
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 0) .AND. (ib .EQ. 2) .AND. (bk .EQ. 4)) Then
          !T0        = R5 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rbx*rby-(5.0D0/2.0D0)*rbx*rby+5.0D0*raz*rbx*Czy+5.0D0*raz*rby*Czx+Czx*Czy))
          !T1(0)     = R7 * (rt3*((-315.0D0/4.0D0)*(raz**(2.0D0))*rbx*rby+(35.0D0/4.0D0)*rbx*rby-(35.0D0/2.0D0)*raz*rbx*Czy-(35.0D0/2.0D0)*raz*rby*Czx &
          !                                    & -(5.0D0/2.0D0)*Czx*Czy))
          !T1(3)     = R6 * (rt3*(35.0D0*raz*rbx*rby+5.0D0*rbx*Czy+5.0D0*rby*Czx))
          !T1(4)     = R6 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rby-(5.0D0/2.0D0)*rby+5.0D0*raz*Czy))
          !T1(5)     = R6 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rbx-(5.0D0/2.0D0)*rbx+5.0D0*raz*Czx))
          !T1(9)     = R5 * (rt3*(5.0D0*raz*rby+Czy))
          !T1(12)    = R5 * (rt3*(5.0D0*raz*rbx+Czx))
          T0        = R5 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rbx*rby-(5.0D0/2.0D0)*rbx*rby-5.0D0*raz*rbx*Czy-5.0D0*raz*rby*Czx+Czx*Czy))
          T1(0)     = R7 * (rt3*((-315.0D0/4.0D0)*(raz**(2.0D0))*rbx*rby+(35.0D0/4.0D0)*rbx*rby+(35.0D0/2.0D0)*raz*rbx*Czy+(35.0D0/2.0D0)*raz*rby*Czx &
                                              & -(5.0D0/2.0D0)*Czx*Czy))
          T1(3)     = R6 * (rt3*(35.0D0*raz*rbx*rby-5.0D0*rbx*Czy-5.0D0*rby*Czx))
          T1(4)     = R6 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rby-(5.0D0/2.0D0)*rby-5.0D0*raz*Czy))
          T1(5)     = R6 * (rt3*((35.0D0/2.0D0)*(raz**(2.0D0))*rbx-(5.0D0/2.0D0)*rbx-5.0D0*raz*Czx))
          T1(9)     = R5 * (rt3*((-5.0D0)*raz*rby+Czy))
          T1(12)    = R5 * (rt3*((-5.0D0)*raz*rbx+Czx))
!
!                      21c, 20
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 1) .AND. (ib .EQ. 2) .AND. (bk .EQ. 0)) Then
          !T0        = R5 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*rax*raz-(5.0D0/2.0D0)*rax*raz+5.0D0*rbz*rax*Czz+5.0D0*rbz*raz*Cxz+Cxz*Czz))
          !T1(0)     = R7 * (rt3*((-315.0D0/4.0D0)*(rbz**(2.0D0))*rax*raz+(35.0D0/4.0D0)*rax*raz-(35.0D0/2.0D0)*rbz*rax*Czz-(35.0D0/2.0D0)*rbz*raz*Cxz &
          !                                     & -(5.0D0/2.0D0)*Cxz*Czz))
          !T1(6)     = R6 * (rt3*(35.0D0*rbz*rax*raz+5.0D0*rax*Czz+5.0D0*raz*Cxz))
          !T1(1)     = R6 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*raz-(5.0D0/2.0D0)*raz+5.0D0*rbz*Czz))
          !T1(3)     = R6 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*rax-(5.0D0/2.0D0)*rax+5.0D0*rbz*Cxz))
          !T1(13)    = R5 * (rt3*(5.0D0*rbz*raz+Czz))
          !T1(15)    = R5 * (rt3*(5.0D0*rbz*rax+Cxz))
          T0        = R5 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*rax*raz-(5.0D0/2.0D0)*rax*raz-5.0D0*rbz*rax*Czz-5.0D0*rbz*raz*Cxz+Cxz*Czz))
          T1(0)     = R7 * (rt3*((-315.0D0/4.0D0)*(rbz**(2.0D0))*rax*raz+(35.0D0/4.0D0)*rax*raz+(35.0D0/2.0D0)*rbz*rax*Czz+(35.0D0/2.0D0)*rbz*raz*Cxz &
                                               & -(5.0D0/2.0D0)*Cxz*Czz))
          T1(6)     = R6 * (rt3*(35.0D0*rbz*rax*raz-5.0D0*rax*Czz-5.0D0*raz*Cxz))
          T1(1)     = R6 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*raz-(5.0D0/2.0D0)*raz-5.0D0*rbz*Czz))
          T1(3)     = R6 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*rax-(5.0D0/2.0D0)*rax-5.0D0*rbz*Cxz))
          T1(13)    = R5 * (rt3*((-5.0D0)*rbz*raz+Czz))
          T1(15)    = R5 * (rt3*((-5.0D0)*rbz*rax+Cxz))
!
!                      21c, 21c
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 1) .AND. (ib .EQ. 2) .AND. (bk .EQ. 1)) Then
          !T0        = R5 * (35.0D0*rax*raz*rbx*rbz+5.0D0*rax*rbx*Czz+5.0D0*rax*rbz*Czx+5.0D0*raz*rbx*Cxz+5.0D0*raz*rbz*Cxx+Cxx*Czz+Cxz*Czx)
          !T1(0)     = R7 * ((-315.0D0/2.0D0)*rax*raz*rbx*rbz-(35.0D0/2.0D0)*rax*rbx*Czz-(35.0D0/2.0D0)*rax*rbz*Czx-(35.0D0/2.0D0)*raz*rbx*Cxz-(35.0D0/2.0D0)*raz*rbz*Cxx &
          !                                     & -(5.0D0/2.0D0)*Cxx*Czz-(5.0D0/2.0D0)*Cxz*Czx)
          !T1(1)     = R6 * (35.0D0*raz*rbx*rbz+5.0D0*rbx*Czz+5.0D0*rbz*Czx)
          !T1(3)     = R6 * (35.0D0*rax*rbx*rbz+5.0D0*rbx*Cxz+5.0D0*rbz*Cxx)
          !T1(4)     = R6 * (35.0D0*rax*raz*rbz+5.0D0*rax*Czz+5.0D0*raz*Cxz)
          !T1(6)     = R6 * (35.0D0*rax*raz*rbx+5.0D0*rax*Czx+5.0D0*raz*Cxx)
          !T1(7)     = R5 * (5.0D0*raz*rbz+Czz)
          !T1(9)     = R5 * (5.0D0*rax*rbz+Cxz)
          !T1(13)    = R5 * (5.0D0*raz*rbx+Czx)
          !T1(15)    = R5 * (5.0D0*rax*rbx+Cxx)
          T0        = R5 * (35.0D0*rax*raz*rbx*rbz-5.0D0*rax*rbx*Czz-5.0D0*rax*rbz*Czx-5.0D0*raz*rbx*Cxz-5.0D0*raz*rbz*Cxx+Cxx*Czz+Cxz*Czx)
          T1(0)     = R7 * ((-315.0D0/2.0D0)*rax*raz*rbx*rbz+(35.0D0/2.0D0)*rax*rbx*Czz+(35.0D0/2.0D0)*rax*rbz*Czx+(35.0D0/2.0D0)*raz*rbx*Cxz+(35.0D0/2.0D0)*raz*rbz*Cxx &
                                               & -(5.0D0/2.0D0)*Cxx*Czz-(5.0D0/2.0D0)*Cxz*Czx)
          T1(1)     = R6 * (35.0D0*raz*rbx*rbz-5.0D0*rbx*Czz-5.0D0*rbz*Czx)
          T1(3)     = R6 * (35.0D0*rax*rbx*rbz-5.0D0*rbx*Cxz-5.0D0*rbz*Cxx)
          T1(4)     = R6 * (35.0D0*rax*raz*rbz-5.0D0*rax*Czz-5.0D0*raz*Cxz)
          T1(6)     = R6 * (35.0D0*rax*raz*rbx-5.0D0*rax*Czx-5.0D0*raz*Cxx)
          T1(7)     = R5 * ((-5.0D0)*raz*rbz+Czz)
          T1(9)     = R5 * ((-5.0D0)*rax*rbz+Cxz)
          T1(13)    = R5 * ((-5.0D0)*raz*rbx+Czx)
          T1(15)    = R5 * ((-5.0D0)*rax*rbx+Cxx)
!
!                      21c, 21s
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 1) .AND. (ib .EQ. 2) .AND. (bk .EQ. 2)) Then
          !T0        = R5 * (35.0D0*rax*raz*rby*rbz+5.0D0*rax*rby*Czz+5.0D0*rax*rbz*Czy+5.0D0*raz*rby*Cxz+5.0D0*raz*rbz*Cxy+Cxy*Czz+Cxz*Czy)
          !T1(0)     = R7 * ((-315.0D0/2.0D0)*rax*raz*rby*rbz-(35.0D0/2.0D0)*rax*rby*Czz-(35.0D0/2.0D0)*rax*rbz*Czy-(35.0D0/2.0D0)*raz*rby*Cxz-(35.0D0/2.0D0)*raz*rbz*Cxy &
          !                                                     & -(5.0D0/2.0D0)*Cxy*Czz-(5.0D0/2.0D0)*Cxz*Czy)
          !T1(1)     = R6 * (35.0D0*raz*rby*rbz+5.0D0*rby*Czz+5.0D0*rbz*Czy)
          !T1(3)     = R6 * (35.0D0*rax*rby*rbz+5.0D0*rby*Cxz+5.0D0*rbz*Cxy)
          !T1(5)     = R6 * (35.0D0*rax*raz*rbz+5.0D0*rax*Czz+5.0D0*raz*Cxz)
          !T1(6)     = R6 * (35.0D0*rax*raz*rby+5.0D0*rax*Czy+5.0D0*raz*Cxy)
          !T1(10)    = R5 * (5.0D0*raz*rbz+Czz)
          !T1(12)    = R5 * (5.0D0*rax*rbz+Cxz)
          !T1(13)    = R5 * (5.0D0*raz*rby+Czy)
          !T1(15)    = R5 * (5.0D0*rax*rby+Cxy)
          T0        = R5 * (35.0D0*rax*raz*rby*rbz-5.0D0*rax*rby*Czz-5.0D0*rax*rbz*Czy-5.0D0*raz*rby*Cxz-5.0D0*raz*rbz*Cxy+Cxy*Czz+Cxz*Czy)
          T1(0)     = R7 * ((-315.0D0/2.0D0)*rax*raz*rby*rbz+(35.0D0/2.0D0)*rax*rby*Czz+(35.0D0/2.0D0)*rax*rbz*Czy+(35.0D0/2.0D0)*raz*rby*Cxz+(35.0D0/2.0D0)*raz*rbz*Cxy &
                                                               & -(5.0D0/2.0D0)*Cxy*Czz-(5.0D0/2.0D0)*Cxz*Czy)
          T1(1)     = R6 * (35.0D0*raz*rby*rbz-5.0D0*rby*Czz-5.0D0*rbz*Czy)
          T1(3)     = R6 * (35.0D0*rax*rby*rbz-5.0D0*rby*Cxz-5.0D0*rbz*Cxy)
          T1(5)     = R6 * (35.0D0*rax*raz*rbz-5.0D0*rax*Czz-5.0D0*raz*Cxz)
          T1(6)     = R6 * (35.0D0*rax*raz*rby-5.0D0*rax*Czy-5.0D0*raz*Cxy)
          T1(10)    = R5 * ((-5.0D0)*raz*rbz+Czz)
          T1(12)    = R5 * ((-5.0D0)*rax*rbz+Cxz)
          T1(13)    = R5 * ((-5.0D0)*raz*rby+Czy)
          T1(15)    = R5 * ((-5.0D0)*rax*rby+Cxy)
!
!                      21c, 22c
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 1) .AND. (ib .EQ. 2) .AND. (bk .EQ. 3)) Then
          !T0        = R5 * ((35.0D0/2.0D0)*rax*raz*(rbx**(2.0D0))-(35.0D0/2.0D0)*rax*raz*(rby**(2.0D0))+5.0D0*rax*rbx*Czx-5.0D0*rax*rby*Czy+5.0D0*raz*rbx*Cxx-5.0D0*raz*rby*Cxy+Cxx*Czx-Cxy*Czy)
          !T1(0)     = R7 * ((-315.0D0/4.0D0)*rax*raz*(rbx**(2.0D0))+(315.0D0/4.0D0)*rax*raz*(rby**(2.0D0))-(35.0D0/2.0D0)*rax*rbx*Czx+(35.0D0/2.0D0)*rax*rby*Czy &
          !                                                  & -(35.0D0/2.0D0)*raz*rbx*Cxx+(35.0D0/2.0D0)*raz*rby*Cxy-(5.0D0/2.0D0)*Cxx*Czx+(5.0D0/2.0D0)*Cxy*Czy)
          !T1(1)     = R6 * ((35.0D0/2.0D0)*raz*(rbx**(2.0D0))-(35.0D0/2.0D0)*raz*(rby**(2.0D0))+5.0D0*rbx*Czx-5.0D0*rby*Czy)
          !T1(3)     = R6 * ((35.0D0/2.0D0)*rax*(rbx**(2.0D0))-(35.0D0/2.0D0)*rax*(rby**(2.0D0))+5.0D0*rbx*Cxx-5.0D0*rby*Cxy)
          !T1(4)     = R6 * (35.0D0*rax*raz*rbx+5.0D0*rax*Czx+5.0D0*raz*Cxx)
          !T1(5)     = R6 * (-35.0D0*rax*raz*rby-5.0D0*rax*Czy-5.0D0*raz*Cxy)
          !T1(7)     = R5 * (5.0D0*raz*rbx+Czx)
          !T1(9)     = R5 * (5.0D0*rax*rbx+Cxx)
          !T1(10)    = R5 * (-5.0D0*raz*rby-Czy)
          !T1(12)    = R5 * (-5.0D0*rax*rby-Cxy)
          T0        = R5 * ((35.0D0/2.0D0)*rax*raz*(rbx**(2.0D0))-(35.0D0/2.0D0)*rax*raz*(rby**(2.0D0))-5.0D0*rax*rbx*Czx+5.0D0*rax*rby*Czy-5.0D0*raz*rbx*Cxx+5.0D0*raz*rby*Cxy+Cxx*Czx-Cxy*Czy)
          T1(0)     = R7 * ((-315.0D0/4.0D0)*rax*raz*(rbx**(2.0D0))+(315.0D0/4.0D0)*rax*raz*(rby**(2.0D0))+(35.0D0/2.0D0)*rax*rbx*Czx-(35.0D0/2.0D0)*rax*rby*Czy &
                                                            & +(35.0D0/2.0D0)*raz*rbx*Cxx-(35.0D0/2.0D0)*raz*rby*Cxy-(5.0D0/2.0D0)*Cxx*Czx+(5.0D0/2.0D0)*Cxy*Czy)
          T1(1)     = R6 * ((35.0D0/2.0D0)*raz*(rbx**(2.0D0))-(35.0D0/2.0D0)*raz*(rby**(2.0D0))-5.0D0*rbx*Czx+5.0D0*rby*Czy)
          T1(3)     = R6 * ((35.0D0/2.0D0)*rax*(rbx**(2.0D0))-(35.0D0/2.0D0)*rax*(rby**(2.0D0))-5.0D0*rbx*Cxx+5.0D0*rby*Cxy)
          T1(4)     = R6 * (35.0D0*rax*raz*rbx-5.0D0*rax*Czx-5.0D0*raz*Cxx)
          T1(5)     = R6 * (-35.0D0*rax*raz*rby+5.0D0*rax*Czy+5.0D0*raz*Cxy)
          T1(7)     = R5 * ((-5.0D0)*raz*rbx+Czx)
          T1(9)     = R5 * ((-5.0D0)*rax*rbx+Cxx)
          T1(10)    = R5 * (5.0D0*raz*rby-Czy)
          T1(12)    = R5 * (5.0D0*rax*rby-Cxy)
!
!                      21c, 22s
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 1) .AND. (ib .EQ. 2) .AND. (bk .EQ. 4)) Then
          !T0        = R5 * (35.0D0*rax*raz*rbx*rby+5.0D0*rax*rbx*Czy+5.0D0*rax*rby*Czx+5.0D0*raz*rbx*Cxy+5.0D0*raz*rby*Cxx+Cxx*Czy+Cxy*Czx)
          !T1(0)     = R7 * ((-315.0D0/2.0D0)*rax*raz*rbx*rby-(35.0D0/2.0D0)*rax*rbx*Czy-(35.0D0/2.0D0)*rax*rby*Czx-(35.0D0/2.0D0)*raz*rbx*Cxy-(35.0D0/2.0D0)*raz*rby*Cxx &
          !                                               & -(5.0D0/2.0D0)*Cxx*Czy-(5.0D0/2.0D0)*Cxy*Czx)
          !T1(1)     = R6 * (35.0D0*raz*rbx*rby+5.0D0*rbx*Czy+5.0D0*rby*Czx)
          !T1(3)     = R6 * (35.0D0*rax*rbx*rby+5.0D0*rbx*Cxy+5.0D0*rby*Cxx)
          !T1(4)     = R6 * (35.0D0*rax*raz*rby+5.0D0*rax*Czy+5.0D0*raz*Cxy)
          !T1(5)     = R6 * (35.0D0*rax*raz*rbx+5.0D0*rax*Czx+5.0D0*raz*Cxx)
          !T1(7)     = R5 * (5.0D0*raz*rby+Czy)
          !T1(9)     = R5 * (5.0D0*rax*rby+Cxy)
          !T1(10)    = R5 * (5.0D0*raz*rbx+Czx)
          !T1(12)    = R5 * (5.0D0*rax*rbx+Cxx)
          T0        = R5 * (35.0D0*rax*raz*rbx*rby-5.0D0*rax*rbx*Czy-5.0D0*rax*rby*Czx-5.0D0*raz*rbx*Cxy-5.0D0*raz*rby*Cxx+Cxx*Czy+Cxy*Czx)
          T1(0)     = R7 * ((-315.0D0/2.0D0)*rax*raz*rbx*rby+(35.0D0/2.0D0)*rax*rbx*Czy+(35.0D0/2.0D0)*rax*rby*Czx+(35.0D0/2.0D0)*raz*rbx*Cxy+(35.0D0/2.0D0)*raz*rby*Cxx &
                                                         & -(5.0D0/2.0D0)*Cxx*Czy-(5.0D0/2.0D0)*Cxy*Czx)
          T1(1)     = R6 * (35.0D0*raz*rbx*rby-5.0D0*rbx*Czy-5.0D0*rby*Czx)
          T1(3)     = R6 * (35.0D0*rax*rbx*rby-5.0D0*rbx*Cxy-5.0D0*rby*Cxx)
          T1(4)     = R6 * (35.0D0*rax*raz*rby-5.0D0*rax*Czy-5.0D0*raz*Cxy)
          T1(5)     = R6 * (35.0D0*rax*raz*rbx-5.0D0*rax*Czx-5.0D0*raz*Cxx)
          T1(7)     = R5 * ((-5.0D0)*raz*rby+Czy)
          T1(9)     = R5 * ((-5.0D0)*rax*rby+Cxy)
          T1(10)    = R5 * ((-5.0D0)*raz*rbx+Czx)
          T1(12)    = R5 * ((-5.0D0)*rax*rbx+Cxx)
!
!                      21s, 20
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 2) .AND. (ib .EQ. 2) .AND. (bk .EQ. 0)) Then
          !T0        = R5 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*ray*raz-(5.0D0/2.0D0)*ray*raz+5.0D0*rbz*ray*Czz+5.0D0*rbz*raz*Cyz+Cyz*Czz))
          !T1(0)     = R7 * (rt3*((-315.0D0/4.0D0)*(rbz**(2.0D0))*ray*raz+(35.0D0/4.0D0)*ray*raz-(35.0D0/2.0D0)*rbz*ray*Czz-(35.0D0/2.0D0)*rbz*raz*Cyz &
          !                                     & -(5.0D0/2.0D0)*Cyz*Czz))
          !T1(6)     = R6 * (rt3*(35.0D0*rbz*ray*raz+5.0D0*ray*Czz+5.0D0*raz*Cyz))
          !T1(2)     = R6 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*raz-(5.0D0/2.0D0)*raz+5.0D0*rbz*Czz))
          !T1(3)     = R6 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*ray-(5.0D0/2.0D0)*ray+5.0D0*rbz*Cyz))
          !T1(14)    = R5 * (rt3*(5.0D0*rbz*raz+Czz))
          !T1(15)    = R5 * (rt3*(5.0D0*rbz*ray+Cyz))
          T0        = R5 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*ray*raz-(5.0D0/2.0D0)*ray*raz-5.0D0*rbz*ray*Czz-5.0D0*rbz*raz*Cyz+Cyz*Czz))
          T1(0)     = R7 * (rt3*((-315.0D0/4.0D0)*(rbz**(2.0D0))*ray*raz+(35.0D0/4.0D0)*ray*raz+(35.0D0/2.0D0)*rbz*ray*Czz+(35.0D0/2.0D0)*rbz*raz*Cyz &
                                               & -(5.0D0/2.0D0)*Cyz*Czz))
          T1(6)     = R6 * (rt3*(35.0D0*rbz*ray*raz-5.0D0*ray*Czz-5.0D0*raz*Cyz))
          T1(2)     = R6 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*raz-(5.0D0/2.0D0)*raz-5.0D0*rbz*Czz))
          T1(3)     = R6 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*ray-(5.0D0/2.0D0)*ray-5.0D0*rbz*Cyz))
          T1(14)    = R5 * (rt3*((-5.0D0)*rbz*raz+Czz))
          T1(15)    = R5 * (rt3*((-5.0D0)*rbz*ray+Cyz))
!
!                      21s, 21c
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 2) .AND. (ib .EQ. 2) .AND. (bk .EQ. 1)) Then
          !T0        = R5 * (35.0D0*rbx*rbz*ray*raz+5.0D0*rbx*ray*Czz+5.0D0*rbx*raz*Cyz+5.0D0*rbz*ray*Czx+5.0D0*rbz*raz*Cyx+Cyx*Czz+Czx*Cyz)
          !T1(0)     = R7 * ((-315.0D0/2.0D0)*rbx*rbz*ray*raz-(35.0D0/2.0D0)*rbx*ray*Czz-(35.0D0/2.0D0)*rbx*raz*Cyz-(35.0D0/2.0D0)*rbz*ray*Czx-(35.0D0/2.0D0)*rbz*raz*Cyx &
          !                                        & -(5.0D0/2.0D0)*Cyx*Czz-(5.0D0/2.0D0)*Czx*Cyz)
          !T1(4)     = R6 * (35.0D0*rbz*ray*raz+5.0D0*ray*Czz+5.0D0*raz*Cyz)
          !T1(6)     = R6 * (35.0D0*rbx*ray*raz+5.0D0*ray*Czx+5.0D0*raz*Cyx)
          !T1(2)     = R6 * (35.0D0*rbx*rbz*raz+5.0D0*rbx*Czz+5.0D0*rbz*Czx)
          !T1(3)     = R6 * (35.0D0*rbx*rbz*ray+5.0D0*rbx*Cyz+5.0D0*rbz*Cyx)
          !T1(8)     = R5 * (5.0D0*rbz*raz+Czz)
          !T1(14)    = R5 * (5.0D0*rbx*raz+Czx)
          !T1(9)     = R5 * (5.0D0*rbz*ray+Cyz)
          !T1(15)    = R5 * (5.0D0*rbx*ray+Cyx)
          T0        = R5 * (35.0D0*rbx*rbz*ray*raz-5.0D0*rbx*ray*Czz-5.0D0*rbx*raz*Cyz-5.0D0*rbz*ray*Czx-5.0D0*rbz*raz*Cyx+Cyx*Czz+Czx*Cyz)
          T1(0)     = R7 * ((-315.0D0/2.0D0)*rbx*rbz*ray*raz+(35.0D0/2.0D0)*rbx*ray*Czz+(35.0D0/2.0D0)*rbx*raz*Cyz+(35.0D0/2.0D0)*rbz*ray*Czx+(35.0D0/2.0D0)*rbz*raz*Cyx &
                                                  & -(5.0D0/2.0D0)*Cyx*Czz-(5.0D0/2.0D0)*Czx*Cyz)
          T1(4)     = R6 * (35.0D0*rbz*ray*raz-5.0D0*ray*Czz-5.0D0*raz*Cyz)
          T1(6)     = R6 * (35.0D0*rbx*ray*raz-5.0D0*ray*Czx-5.0D0*raz*Cyx)
          T1(2)     = R6 * (35.0D0*rbx*rbz*raz-5.0D0*rbx*Czz-5.0D0*rbz*Czx)
          T1(3)     = R6 * (35.0D0*rbx*rbz*ray-5.0D0*rbx*Cyz-5.0D0*rbz*Cyx)
          T1(8)     = R5 * ((-5.0D0)*rbz*raz+Czz)
          T1(14)    = R5 * ((-5.0D0)*rbx*raz+Czx)
          T1(9)     = R5 * ((-5.0D0)*rbz*ray+Cyz)
          T1(15)    = R5 * ((-5.0D0)*rbx*ray+Cyx)
!
!                      21s, 21s
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 2) .AND. (ib .EQ. 2) .AND. (bk .EQ. 2)) Then
          !T0        = R5 * (35.0D0*ray*raz*rby*rbz+5.0D0*ray*rby*Czz+5.0D0*ray*rbz*Czy+5.0D0*raz*rby*Cyz+5.0D0*raz*rbz*Cyy+Cyy*Czz+Cyz*Czy)
          !T1(0)     = R7 * ((-315.0D0/2.0D0)*ray*raz*rby*rbz-(35.0D0/2.0D0)*ray*rby*Czz-(35.0D0/2.0D0)*ray*rbz*Czy-(35.0D0/2.0D0)*raz*rby*Cyz-(35.0D0/2.0D0)*raz*rbz*Cyy &
          !                               & -(5.0D0/2.0D0)*Cyy*Czz-(5.0D0/2.0D0)*Cyz*Czy)
          !T1(2)     = R6 * (35.0D0*raz*rby*rbz+5.0D0*rby*Czz+5.0D0*rbz*Czy)
          !T1(3)     = R6 * (35.0D0*ray*rby*rbz+5.0D0*rby*Cyz+5.0D0*rbz*Cyy)
          !T1(5)     = R6 * (35.0D0*ray*raz*rbz+5.0D0*ray*Czz+5.0D0*raz*Cyz)
          !T1(6)     = R6 * (35.0D0*ray*raz*rby+5.0D0*ray*Czy+5.0D0*raz*Cyy)
          !T1(11)    = R5 * (5.0D0*raz*rbz+Czz)
          !T1(12)    = R5 * (5.0D0*ray*rbz+Cyz)
          !T1(14)    = R5 * (5.0D0*raz*rby+Czy)
          !T1(15)    = R5 * (5.0D0*ray*rby+Cyy)
          T0        = R5 * (35.0D0*ray*raz*rby*rbz-5.0D0*ray*rby*Czz-5.0D0*ray*rbz*Czy-5.0D0*raz*rby*Cyz-5.0D0*raz*rbz*Cyy+Cyy*Czz+Cyz*Czy)
          T1(0)     = R7 * ((-315.0D0/2.0D0)*ray*raz*rby*rbz+(35.0D0/2.0D0)*ray*rby*Czz+(35.0D0/2.0D0)*ray*rbz*Czy+(35.0D0/2.0D0)*raz*rby*Cyz+(35.0D0/2.0D0)*raz*rbz*Cyy &
                                         & -(5.0D0/2.0D0)*Cyy*Czz-(5.0D0/2.0D0)*Cyz*Czy)
          T1(2)     = R6 * (35.0D0*raz*rby*rbz-5.0D0*rby*Czz-5.0D0*rbz*Czy)
          T1(3)     = R6 * (35.0D0*ray*rby*rbz-5.0D0*rby*Cyz-5.0D0*rbz*Cyy)
          T1(5)     = R6 * (35.0D0*ray*raz*rbz-5.0D0*ray*Czz-5.0D0*raz*Cyz)
          T1(6)     = R6 * (35.0D0*ray*raz*rby-5.0D0*ray*Czy-5.0D0*raz*Cyy)
          T1(11)    = R5 * ((-5.0D0)*raz*rbz+Czz)
          T1(12)    = R5 * ((-5.0D0)*ray*rbz+Cyz)
          T1(14)    = R5 * ((-5.0D0)*raz*rby+Czy)
          T1(15)    = R5 * ((-5.0D0)*ray*rby+Cyy)
!
!                      21s, 22c
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 2) .AND. (ib .EQ. 2) .AND. (bk .EQ. 3)) Then
          !T0        = R5 * ((35.0D0/2.0D0)*ray*raz*(rbx**(2.0D0))-(35.0D0/2.0D0)*ray*raz*(rby**(2.0D0))+5.0D0*ray*rbx*Czx-5.0D0*ray*rby*Czy+5.0D0*raz*rbx*Cyx-5.0D0*raz*rby*Cyy+Cyx*Czx-Cyy*Czy)
          !T1(0)     = R7 * ((-315.0D0/4.0D0)*ray*raz*(rbx**(2.0D0))+(315.0D0/4.0D0)*ray*raz*(rby**(2.0D0))-(35.0D0/2.0D0)*ray*rbx*Czx+(35.0D0/2.0D0)*ray*rby*Czy &
          !                               & -(35.0D0/2.0D0)*raz*rbx*Cyx+(35.0D0/2.0D0)*raz*rby*Cyy-(5.0D0/2.0D0)*Cyx*Czx+(5.0D0/2.0D0)*Cyy*Czy)
          !T1(2)     = R6 * ((35.0D0/2.0D0)*raz*(rbx**(2.0D0))-(35.0D0/2.0D0)*raz*(rby**(2.0D0))+5.0D0*rbx*Czx-5.0D0*rby*Czy)
          !T1(3)     = R6 * ((35.0D0/2.0D0)*ray*(rbx**(2.0D0))-(35.0D0/2.0D0)*ray*(rby**(2.0D0))+5.0D0*rbx*Cyx-5.0D0*rby*Cyy)
          !T1(4)     = R6 * (35.0D0*ray*raz*rbx+5.0D0*ray*Czx+5.0D0*raz*Cyx)
          !T1(5)     = R6 * (-35.0D0*ray*raz*rby-5.0D0*ray*Czy-5.0D0*raz*Cyy)
          !T1(8)     = R5 * (5.0D0*raz*rbx+Czx)
          !T1(9)     = R5 * (5.0D0*ray*rbx+Cyx)
          !T1(11)    = R5 * (-5.0D0*raz*rby-Czy)
          !T1(12)    = R5 * (-5.0D0*ray*rby-Cyy)
          T0        = R5 * ((35.0D0/2.0D0)*ray*raz*(rbx**(2.0D0))-(35.0D0/2.0D0)*ray*raz*(rby**(2.0D0))-5.0D0*ray*rbx*Czx+5.0D0*ray*rby*Czy-5.0D0*raz*rbx*Cyx+5.0D0*raz*rby*Cyy+Cyx*Czx-Cyy*Czy)
          T1(0)     = R7 * ((-315.0D0/4.0D0)*ray*raz*(rbx**(2.0D0))+(315.0D0/4.0D0)*ray*raz*(rby**(2.0D0))+(35.0D0/2.0D0)*ray*rbx*Czx-(35.0D0/2.0D0)*ray*rby*Czy &
                                         & +(35.0D0/2.0D0)*raz*rbx*Cyx-(35.0D0/2.0D0)*raz*rby*Cyy-(5.0D0/2.0D0)*Cyx*Czx+(5.0D0/2.0D0)*Cyy*Czy)
          T1(2)     = R6 * ((35.0D0/2.0D0)*raz*(rbx**(2.0D0))-(35.0D0/2.0D0)*raz*(rby**(2.0D0))-5.0D0*rbx*Czx+5.0D0*rby*Czy)
          T1(3)     = R6 * ((35.0D0/2.0D0)*ray*(rbx**(2.0D0))-(35.0D0/2.0D0)*ray*(rby**(2.0D0))-5.0D0*rbx*Cyx+5.0D0*rby*Cyy)
          T1(4)     = R6 * (35.0D0*ray*raz*rbx-5.0D0*ray*Czx-5.0D0*raz*Cyx)
          T1(5)     = R6 * (-35.0D0*ray*raz*rby+5.0D0*ray*Czy+5.0D0*raz*Cyy)
          T1(8)     = R5 * ((-5.0D0)*raz*rbx+Czx)
          T1(9)     = R5 * ((-5.0D0)*ray*rbx+Cyx)
          T1(11)    = R5 * (5.0D0*raz*rby-Czy)
          T1(12)    = R5 * (5.0D0*ray*rby-Cyy)
!
!                      21s, 22s
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 2) .AND. (ib .EQ. 2) .AND. (bk .EQ. 4)) Then
          !T0        = R5 * (35.0D0*ray*raz*rbx*rby+5.0D0*ray*rbx*Czy+5.0D0*ray*rby*Czx+5.0D0*raz*rbx*Cyy+5.0D0*raz*rby*Cyx+Cyx*Czy+Cyy*Czx)
          !T1(0)     = R7 * ((-315.0D0/2.0D0)*ray*raz*rbx*rby-(35.0D0/2.0D0)*ray*rbx*Czy-(35.0D0/2.0D0)*ray*rby*Czx-(35.0D0/2.0D0)*raz*rbx*Cyy-(35.0D0/2.0D0)*raz*rby*Cyx &
          !                                & -(5.0D0/2.0D0)*Cyx*Czy-(5.0D0/2.0D0)*Cyy*Czx)
          !T1(2)     = R6 * (35.0D0*raz*rbx*rby+5.0D0*rbx*Czy+5.0D0*rby*Czx)
          !T1(3)     = R6 * (35.0D0*ray*rbx*rby+5.0D0*rbx*Cyy+5.0D0*rby*Cyx)
          !T1(4)     = R6 * (35.0D0*ray*raz*rby+5.0D0*ray*Czy+5.0D0*raz*Cyy)
          !T1(5)     = R6 * (35.0D0*ray*raz*rbx+5.0D0*ray*Czx+5.0D0*raz*Cyx)
          !T1(8)     = R5 * (5.0D0*raz*rby+Czy)
          !T1(9)     = R5 * (5.0D0*ray*rby+Cyy)
          !T1(11)    = R5 * (5.0D0*raz*rbx+Czx)
          !T1(12)    = R5 * (5.0D0*ray*rbx+Cyx)
          T0        = R5 * (35.0D0*ray*raz*rbx*rby-5.0D0*ray*rbx*Czy-5.0D0*ray*rby*Czx-5.0D0*raz*rbx*Cyy-5.0D0*raz*rby*Cyx+Cyx*Czy+Cyy*Czx)
          T1(0)     = R7 * ((-315.0D0/2.0D0)*ray*raz*rbx*rby+(35.0D0/2.0D0)*ray*rbx*Czy+(35.0D0/2.0D0)*ray*rby*Czx+(35.0D0/2.0D0)*raz*rbx*Cyy+(35.0D0/2.0D0)*raz*rby*Cyx &
                                          & -(5.0D0/2.0D0)*Cyx*Czy-(5.0D0/2.0D0)*Cyy*Czx)
          T1(2)     = R6 * (35.0D0*raz*rbx*rby-5.0D0*rbx*Czy-5.0D0*rby*Czx)
          T1(3)     = R6 * (35.0D0*ray*rbx*rby-5.0D0*rbx*Cyy-5.0D0*rby*Cyx)
          T1(4)     = R6 * (35.0D0*ray*raz*rby-5.0D0*ray*Czy-5.0D0*raz*Cyy)
          T1(5)     = R6 * (35.0D0*ray*raz*rbx-5.0D0*ray*Czx-5.0D0*raz*Cyx)
          T1(8)     = R5 * ((-5.0D0)*raz*rby+Czy)
          T1(9)     = R5 * ((-5.0D0)*ray*rby+Cyy)
          T1(11)    = R5 * ((-5.0D0)*raz*rbx+Czx)
          T1(12)    = R5 * ((-5.0D0)*ray*rbx+Cyx)
!
!                      22c, 20
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 3) .AND. (ib .EQ. 2) .AND. (bk .EQ. 0)) Then
          !T0        = R5 * (rt3*((35.0D0/4.0D0)*(rbz**(2.0D0))*(rax**(2.0D0))-(35.0D0/4.0D0)*(rbz**(2.0D0))*(ray**(2.0D0))-(5.0D0/4.0D0)*(rax**(2.0D0))+(5.0D0/4.0D0)*(ray**(2.0D0)) &
          !                                    & +5.0D0*rbz*rax*Cxz-5.0D0*rbz*ray*Cyz+(1.0D0/2.0D0)*(Cxz**(2.0D0))-(1.0D0/2.0D0)*(Cyz**(2.0D0))))
          !T1(0)     = R7 * (rt3*((-315.0D0/8.0D0)*(rbz**(2.0D0))*(rax**(2.0D0))+(315.0D0/8.0D0)*(rbz**(2.0D0))*(ray**(2.0D0))+(35.0D0/8.0D0)*(rax**(2.0D0))-(35.0D0/8.0D0)*(ray**(2.0D0)) &
          !                                    & -(35.0D0/2.0D0)*rbz*rax*Cxz+(35.0D0/2.0D0)*rbz*ray*Cyz-(5.0D0/4.0D0)*(Cxz**(2.0D0))+(5.0D0/4.0D0)*(Cyz**(2.0D0))))
          !T1(6)     = R6 * (rt3*((35.0D0/2.0D0)*rbz*(rax**(2.0D0))-(35.0D0/2.0D0)*rbz*(ray**(2.0D0))+5.0D0*rax*Cxz-5.0D0*ray*Cyz))
          !T1(1)     = R6 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*rax-(5.0D0/2.0D0)*rax+5.0D0*rbz*Cxz))
          !T1(2)     = R6 * (rt3*((-35.0D0/2.0D0)*(rbz**(2.0D0))*ray+(5.0D0/2.0D0)*ray-5.0D0*rbz*Cyz))
          !T1(13)    = R5 * (rt3*(5.0D0*rbz*rax+Cxz))
          !T1(14)    = R5 * (rt3*(-5.0D0*rbz*ray-Cyz))
          T0        = R5 * (rt3*((35.0D0/4.0D0)*(rbz**(2.0D0))*(rax**(2.0D0))-(35.0D0/4.0D0)*(rbz**(2.0D0))*(ray**(2.0D0))-(5.0D0/4.0D0)*(rax**(2.0D0))+(5.0D0/4.0D0)*(ray**(2.0D0)) &
                                              & -5.0D0*rbz*rax*Cxz+5.0D0*rbz*ray*Cyz+(1.0D0/2.0D0)*(Cxz**(2.0D0))-(1.0D0/2.0D0)*(Cyz**(2.0D0))))
          T1(0)     = R7 * (rt3*((-315.0D0/8.0D0)*(rbz**(2.0D0))*(rax**(2.0D0))+(315.0D0/8.0D0)*(rbz**(2.0D0))*(ray**(2.0D0))+(35.0D0/8.0D0)*(rax**(2.0D0))-(35.0D0/8.0D0)*(ray**(2.0D0)) &
                                              & +(35.0D0/2.0D0)*rbz*rax*Cxz-(35.0D0/2.0D0)*rbz*ray*Cyz-(5.0D0/4.0D0)*(Cxz**(2.0D0))+(5.0D0/4.0D0)*(Cyz**(2.0D0))))
          T1(6)     = R6 * (rt3*((35.0D0/2.0D0)*rbz*(rax**(2.0D0))-(35.0D0/2.0D0)*rbz*(ray**(2.0D0))-5.0D0*rax*Cxz+5.0D0*ray*Cyz))
          T1(1)     = R6 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*rax-(5.0D0/2.0D0)*rax-5.0D0*rbz*Cxz))
          T1(2)     = R6 * (rt3*((-35.0D0/2.0D0)*(rbz**(2.0D0))*ray+(5.0D0/2.0D0)*ray+5.0D0*rbz*Cyz))
          T1(13)    = R5 * (rt3*((-5.0D0)*rbz*rax+Cxz))
          T1(14)    = R5 * (rt3*(5.0D0*rbz*ray-Cyz))
!
!                      22c, 21c
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 3) .AND. (ib .EQ. 2) .AND. (bk .EQ. 1)) Then
          !T0        = R5 * ((35.0D0/2.0D0)*rbx*rbz*(rax**(2.0D0))-(35.0D0/2.0D0)*rbx*rbz*(ray**(2.0D0))+5.0D0*rbx*rax*Cxz-5.0D0*rbx*ray*Cyz+5.0D0*rbz*rax*Cxx-5.0D0*rbz*ray*Cyx+Cxx*Cxz-Cyx*Cyz)
          !T1(0)     = R7 * ((-315.0D0/4.0D0)*rbx*rbz*(rax**(2.0D0))+(315.0D0/4.0D0)*rbx*rbz*(ray**(2.0D0))-(35.0D0/2.0D0)*rbx*rax*Cxz+(35.0D0/2.0D0)*rbx*ray*Cyz &
          !                                    & -(35.0D0/2.0D0)*rbz*rax*Cxx+(35.0D0/2.0D0)*rbz*ray*Cyx-(5.0D0/2.0D0)*Cxx*Cxz+(5.0D0/2.0D0)*Cyx*Cyz)
          !T1(4)     = R6 * ((35.0D0/2.0D0)*rbz*(rax**(2.0D0))-(35.0D0/2.0D0)*rbz*(ray**(2.0D0))+5.0D0*rax*Cxz-5.0D0*ray*Cyz)
          !T1(6)     = R6 * ((35.0D0/2.0D0)*rbx*(rax**(2.0D0))-(35.0D0/2.0D0)*rbx*(ray**(2.0D0))+5.0D0*rax*Cxx-5.0D0*ray*Cyx)
          !T1(1)     = R6 * (35.0D0*rbx*rbz*rax+5.0D0*rbx*Cxz+5.0D0*rbz*Cxx)
          !T1(2)     = R6 * (-35.0D0*rbx*rbz*ray-5.0D0*rbx*Cyz-5.0D0*rbz*Cyx)
          !T1(7)     = R5 * (5.0D0*rbz*rax+Cxz)
          !T1(13)    = R5 * (5.0D0*rbx*rax+Cxx)
          !T1(8)     = R5 * (-5.0D0*rbz*ray-Cyz)
          !T1(14)    = R5 * (-5.0D0*rbx*ray-Cyx)
          T0        = R5 * ((35.0D0/2.0D0)*rbx*rbz*(rax**(2.0D0))-(35.0D0/2.0D0)*rbx*rbz*(ray**(2.0D0))-5.0D0*rbx*rax*Cxz+5.0D0*rbx*ray*Cyz-5.0D0*rbz*rax*Cxx+5.0D0*rbz*ray*Cyx+Cxx*Cxz-Cyx*Cyz)
          T1(0)     = R7 * ((-315.0D0/4.0D0)*rbx*rbz*(rax**(2.0D0))+(315.0D0/4.0D0)*rbx*rbz*(ray**(2.0D0))+(35.0D0/2.0D0)*rbx*rax*Cxz-(35.0D0/2.0D0)*rbx*ray*Cyz &
                                              & +(35.0D0/2.0D0)*rbz*rax*Cxx-(35.0D0/2.0D0)*rbz*ray*Cyx-(5.0D0/2.0D0)*Cxx*Cxz+(5.0D0/2.0D0)*Cyx*Cyz)
          T1(4)     = R6 * ((35.0D0/2.0D0)*rbz*(rax**(2.0D0))-(35.0D0/2.0D0)*rbz*(ray**(2.0D0))-5.0D0*rax*Cxz+5.0D0*ray*Cyz)
          T1(6)     = R6 * ((35.0D0/2.0D0)*rbx*(rax**(2.0D0))-(35.0D0/2.0D0)*rbx*(ray**(2.0D0))-5.0D0*rax*Cxx+5.0D0*ray*Cyx)
          T1(1)     = R6 * (35.0D0*rbx*rbz*rax-5.0D0*rbx*Cxz-5.0D0*rbz*Cxx)
          T1(2)     = R6 * (-35.0D0*rbx*rbz*ray+5.0D0*rbx*Cyz+5.0D0*rbz*Cyx)
          T1(7)     = R5 * ((-5.0D0)*rbz*rax+Cxz)
          T1(13)    = R5 * ((-5.0D0)*rbx*rax+Cxx)
          T1(8)     = R5 * (5.0D0*rbz*ray-Cyz)
          T1(14)    = R5 * (5.0D0*rbx*ray-Cyx)
!
!                      22c, 21s
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 3) .AND. (ib .EQ. 2) .AND. (bk .EQ. 2)) Then
          !T0        = R5 * ((35.0D0/2.0D0)*rby*rbz*(rax**(2.0D0))-(35.0D0/2.0D0)*rby*rbz*(ray**(2.0D0))+5.0D0*rby*rax*Cxz-5.0D0*rby*ray*Cyz+5.0D0*rbz*rax*Cxy-5.0D0*rbz*ray*Cyy+Cxy*Cxz-Cyy*Cyz)
          !T1(0)     = R7 * ((-315.0D0/4.0D0)*rby*rbz*(rax**(2.0D0))+(315.0D0/4.0D0)*rby*rbz*(ray**(2.0D0))-(35.0D0/2.0D0)*rby*rax*Cxz+(35.0D0/2.0D0)*rby*ray*Cyz &
          !                                    & -(35.0D0/2.0D0)*rbz*rax*Cxy+(35.0D0/2.0D0)*rbz*ray*Cyy-(5.0D0/2.0D0)*Cxy*Cxz+(5.0D0/2.0D0)*Cyy*Cyz)
          !T1(5)     = R6 * ((35.0D0/2.0D0)*rbz*(rax**(2.0D0))-(35.0D0/2.0D0)*rbz*(ray**(2.0D0))+5.0D0*rax*Cxz-5.0D0*ray*Cyz)
          !T1(6)     = R6 * ((35.0D0/2.0D0)*rby*(rax**(2.0D0))-(35.0D0/2.0D0)*rby*(ray**(2.0D0))+5.0D0*rax*Cxy-5.0D0*ray*Cyy)
          !T1(1)     = R6 * (35.0D0*rby*rbz*rax+5.0D0*rby*Cxz+5.0D0*rbz*Cxy)
          !T1(2)     = R6 * (-35.0D0*rby*rbz*ray-5.0D0*rby*Cyz-5.0D0*rbz*Cyy)
          !T1(10)    = R5 * (5.0D0*rbz*rax+Cxz)
          !T1(13)    = R5 * (5.0D0*rby*rax+Cxy)
          !T1(11)    = R5 * (-5.0D0*rbz*ray-Cyz)
          !T1(14)    = R5 * (-5.0D0*rby*ray-Cyy)
          T0        = R5 * ((35.0D0/2.0D0)*rby*rbz*(rax**(2.0D0))-(35.0D0/2.0D0)*rby*rbz*(ray**(2.0D0))-5.0D0*rby*rax*Cxz+5.0D0*rby*ray*Cyz-5.0D0*rbz*rax*Cxy+5.0D0*rbz*ray*Cyy+Cxy*Cxz-Cyy*Cyz)
          T1(0)     = R7 * ((-315.0D0/4.0D0)*rby*rbz*(rax**(2.0D0))+(315.0D0/4.0D0)*rby*rbz*(ray**(2.0D0))+(35.0D0/2.0D0)*rby*rax*Cxz-(35.0D0/2.0D0)*rby*ray*Cyz &
                                              & +(35.0D0/2.0D0)*rbz*rax*Cxy-(35.0D0/2.0D0)*rbz*ray*Cyy-(5.0D0/2.0D0)*Cxy*Cxz+(5.0D0/2.0D0)*Cyy*Cyz)
          T1(5)     = R6 * ((35.0D0/2.0D0)*rbz*(rax**(2.0D0))-(35.0D0/2.0D0)*rbz*(ray**(2.0D0))-5.0D0*rax*Cxz+5.0D0*ray*Cyz)
          T1(6)     = R6 * ((35.0D0/2.0D0)*rby*(rax**(2.0D0))-(35.0D0/2.0D0)*rby*(ray**(2.0D0))-5.0D0*rax*Cxy+5.0D0*ray*Cyy)
          T1(1)     = R6 * (35.0D0*rby*rbz*rax-5.0D0*rby*Cxz-5.0D0*rbz*Cxy)
          T1(2)     = R6 * (-35.0D0*rby*rbz*ray+5.0D0*rby*Cyz+5.0D0*rbz*Cyy)
          T1(10)    = R5 * ((-5.0D0)*rbz*rax+Cxz)
          T1(13)    = R5 * ((-5.0D0)*rby*rax+Cxy)
          T1(11)    = R5 * (5.0D0*rbz*ray-Cyz)
          T1(14)    = R5 * (5.0D0*rby*ray-Cyy)
!
!                      22c, 22c
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 3) .AND. (ib .EQ. 2) .AND. (bk .EQ. 3)) Then
          !T0        = R5 * ((35.0D0/4.0D0)*(rax**(2.0D0))*(rbx**(2.0D0))-(35.0D0/4.0D0)*(rax**(2.0D0))*(rby**(2.0D0))-(35.0D0/4.0D0)*(ray**(2.0D0))*(rbx**(2.0D0))+(35.0D0/4.0D0)*(ray**(2.0D0))*(rby**(2.0D0)) &
          !                      & +5.0D0*rax*rbx*Cxx-5.0D0*rax*rby*Cxy-5.0D0*ray*rbx*Cyx+5.0D0*ray*rby*Cyy &
          !                      & +(1.0D0/2.0D0)*(Cxx**(2.0D0))-(1.0D0/2.0D0)*(Cxy**(2.0D0))-(1.0D0/2.0D0)*(Cyx**(2.0D0))+(1.0D0/2.0D0)*(Cyy**(2.0D0)))
          !T1(0)     = R7 * ((-315.0D0/8.0D0)*(rax**(2.0D0))*(rbx**(2.0D0))+(315.0D0/8.0D0)*(rax**(2.0D0))*(rby**(2.0D0))+(315.0D0/8.0D0)*(ray**(2.0D0))*(rbx**(2.0D0))-(315.0D0/8.0D0)*(ray**(2.0D0))*(rby**(2.0D0)) &
          !                      & -(35.0D0/2.0D0)*rax*rbx*Cxx+(35.0D0/2.0D0)*rax*rby*Cxy+(35.0D0/2.0D0)*ray*rbx*Cyx-(35.0D0/2.0D0)*ray*rby*Cyy &
          !                      & -(5.0D0/4.0D0)*(Cxx**(2.0D0))+(5.0D0/4.0D0)*(Cxy**(2.0D0))+(5.0D0/4.0D0)*(Cyx**(2.0D0))-(5.0D0/4.0D0)*(Cyy**(2.0D0)))
          !T1(1)     = R6 * ((35.0D0/2.0D0)*rax*(rbx**(2.0D0))-(35.0D0/2.0D0)*rax*(rby**(2.0D0))+5.0D0*rbx*Cxx-5.0D0*rby*Cxy)
          !T1(2)     = R6 * ((-35.0D0/2.0D0)*ray*(rbx**(2.0D0))+(35.0D0/2.0D0)*ray*(rby**(2.0D0))-5.0D0*rbx*Cyx+5.0D0*rby*Cyy)
          !T1(4)     = R6 * ((35.0D0/2.0D0)*(rax**(2.0D0))*rbx-(35.0D0/2.0D0)*(ray**(2.0D0))*rbx+5.0D0*rax*Cxx-5.0D0*ray*Cyx)
          !T1(5)     = R6 * ((-35.0D0/2.0D0)*(rax**(2.0D0))*rby+(35.0D0/2.0D0)*(ray**(2.0D0))*rby-5.0D0*rax*Cxy+5.0D0*ray*Cyy)
          !T1(7)     = R5 * (5.0D0*rax*rbx+Cxx)
          !T1(8)     = R5 * (-5.0D0*ray*rbx-Cyx)
          !T1(10)    = R5 * (-5.0D0*rax*rby-Cxy)
          !T1(11)    = R5 * (5.0D0*ray*rby+Cyy)
          T0        = R5 * ((35.0D0/4.0D0)*(rax**(2.0D0))*(rbx**(2.0D0))-(35.0D0/4.0D0)*(rax**(2.0D0))*(rby**(2.0D0))-(35.0D0/4.0D0)*(ray**(2.0D0))*(rbx**(2.0D0))+(35.0D0/4.0D0)*(ray**(2.0D0))*(rby**(2.0D0)) &
                                & -5.0D0*rax*rbx*Cxx+5.0D0*rax*rby*Cxy+5.0D0*ray*rbx*Cyx-5.0D0*ray*rby*Cyy &
                                & +(1.0D0/2.0D0)*(Cxx**(2.0D0))-(1.0D0/2.0D0)*(Cxy**(2.0D0))-(1.0D0/2.0D0)*(Cyx**(2.0D0))+(1.0D0/2.0D0)*(Cyy**(2.0D0)))
          T1(0)     = R7 * ((-315.0D0/8.0D0)*(rax**(2.0D0))*(rbx**(2.0D0))+(315.0D0/8.0D0)*(rax**(2.0D0))*(rby**(2.0D0))+(315.0D0/8.0D0)*(ray**(2.0D0))*(rbx**(2.0D0))-(315.0D0/8.0D0)*(ray**(2.0D0))*(rby**(2.0D0)) &
                                & +(35.0D0/2.0D0)*rax*rbx*Cxx-(35.0D0/2.0D0)*rax*rby*Cxy-(35.0D0/2.0D0)*ray*rbx*Cyx+(35.0D0/2.0D0)*ray*rby*Cyy &
                                & -(5.0D0/4.0D0)*(Cxx**(2.0D0))+(5.0D0/4.0D0)*(Cxy**(2.0D0))+(5.0D0/4.0D0)*(Cyx**(2.0D0))-(5.0D0/4.0D0)*(Cyy**(2.0D0)))
          T1(1)     = R6 * ((35.0D0/2.0D0)*rax*(rbx**(2.0D0))-(35.0D0/2.0D0)*rax*(rby**(2.0D0))-5.0D0*rbx*Cxx+5.0D0*rby*Cxy)
          T1(2)     = R6 * ((-35.0D0/2.0D0)*ray*(rbx**(2.0D0))+(35.0D0/2.0D0)*ray*(rby**(2.0D0))+5.0D0*rbx*Cyx-5.0D0*rby*Cyy)
          T1(4)     = R6 * ((35.0D0/2.0D0)*(rax**(2.0D0))*rbx-(35.0D0/2.0D0)*(ray**(2.0D0))*rbx-5.0D0*rax*Cxx+5.0D0*ray*Cyx)
          T1(5)     = R6 * ((-35.0D0/2.0D0)*(rax**(2.0D0))*rby+(35.0D0/2.0D0)*(ray**(2.0D0))*rby+5.0D0*rax*Cxy-5.0D0*ray*Cyy)
          T1(7)     = R5 * ((-5.0D0)*rax*rbx+Cxx)
          T1(8)     = R5 * (5.0D0*ray*rbx-Cyx)
          T1(10)    = R5 * (5.0D0*rax*rby-Cxy)
          T1(11)    = R5 * ((-5.0D0)*ray*rby+Cyy)
!
!                      22c, 22s
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 3) .AND. (ib .EQ. 2) .AND. (bk .EQ. 4)) Then
          !T0        = R5 * ((35.0D0/2.0D0)*rbx*rby*(rax**(2.0D0))-(35.0D0/2.0D0)*rbx*rby*(ray**(2.0D0))+5.0D0*rbx*rax*Cxy-5.0D0*rbx*ray*Cyy+5.0D0*rby*rax*Cxx-5.0D0*rby*ray*Cyx+Cxx*Cxy-Cyx*Cyy)
          !T1(0)     = R7 * ((-315.0D0/4.0D0)*rbx*rby*(rax**(2.0D0))+(315.0D0/4.0D0)*rbx*rby*(ray**(2.0D0))-(35.0D0/2.0D0)*rbx*rax*Cxy+(35.0D0/2.0D0)*rbx*ray*Cyy &
          !                      & -(35.0D0/2.0D0)*rby*rax*Cxx+(35.0D0/2.0D0)*rby*ray*Cyx-(5.0D0/2.0D0)*Cxx*Cxy+(5.0D0/2.0D0)*Cyx*Cyy)
          !T1(4)     = R6 * ((35.0D0/2.0D0)*rby*(rax**(2.0D0))-(35.0D0/2.0D0)*rby*(ray**(2.0D0))+5.0D0*rax*Cxy-5.0D0*ray*Cyy)
          !T1(5)     = R6 * ((35.0D0/2.0D0)*rbx*(rax**(2.0D0))-(35.0D0/2.0D0)*rbx*(ray**(2.0D0))+5.0D0*rax*Cxx-5.0D0*ray*Cyx)
          !T1(1)     = R6 * (35.0D0*rbx*rby*rax+5.0D0*rbx*Cxy+5.0D0*rby*Cxx)
          !T1(2)     = R6 * (-35.0D0*rbx*rby*ray-5.0D0*rbx*Cyy-5.0D0*rby*Cyx)
          !T1(7)     = R5 * (5.0D0*rby*rax+Cxy)
          !T1(10)    = R5 * (5.0D0*rbx*rax+Cxx)
          !T1(8)     = R5 * (-5.0D0*rby*ray-Cyy)
          !T1(11)    = R5 * (-5.0D0*rbx*ray-Cyx)
          T0        = R5 * ((35.0D0/2.0D0)*rbx*rby*(rax**(2.0D0))-(35.0D0/2.0D0)*rbx*rby*(ray**(2.0D0))-5.0D0*rbx*rax*Cxy+5.0D0*rbx*ray*Cyy-5.0D0*rby*rax*Cxx+5.0D0*rby*ray*Cyx+Cxx*Cxy-Cyx*Cyy)
          T1(0)     = R7 * ((-315.0D0/4.0D0)*rbx*rby*(rax**(2.0D0))+(315.0D0/4.0D0)*rbx*rby*(ray**(2.0D0))+(35.0D0/2.0D0)*rbx*rax*Cxy-(35.0D0/2.0D0)*rbx*ray*Cyy &
                                & +(35.0D0/2.0D0)*rby*rax*Cxx-(35.0D0/2.0D0)*rby*ray*Cyx-(5.0D0/2.0D0)*Cxx*Cxy+(5.0D0/2.0D0)*Cyx*Cyy)
          T1(4)     = R6 * ((35.0D0/2.0D0)*rby*(rax**(2.0D0))-(35.0D0/2.0D0)*rby*(ray**(2.0D0))-5.0D0*rax*Cxy+5.0D0*ray*Cyy)
          T1(5)     = R6 * ((35.0D0/2.0D0)*rbx*(rax**(2.0D0))-(35.0D0/2.0D0)*rbx*(ray**(2.0D0))-5.0D0*rax*Cxx+5.0D0*ray*Cyx)
          T1(1)     = R6 * (35.0D0*rbx*rby*rax-5.0D0*rbx*Cxy-5.0D0*rby*Cxx)
          T1(2)     = R6 * (-35.0D0*rbx*rby*ray+5.0D0*rbx*Cyy+5.0D0*rby*Cyx)
          T1(7)     = R5 * ((-5.0D0)*rby*rax+Cxy)
          T1(10)    = R5 * ((-5.0D0)*rbx*rax+Cxx)
          T1(8)     = R5 * (5.0D0*rby*ray-Cyy)
          T1(11)    = R5 * (5.0D0*rbx*ray-Cyx)!HERE
!
!                      22s, 20
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 4) .AND. (ib .EQ. 2) .AND. (bk .EQ. 0)) Then
          !T0        = R5 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*rax*ray-(5.0D0/2.0D0)*rax*ray+5.0D0*rbz*rax*Cyz+5.0D0*rbz*ray*Cxz+Cxz*Cyz))
          !T1(0)     = R7 * (rt3*((-315.0D0/4.0D0)*(rbz**(2.0D0))*rax*ray+(35.0D0/4.0D0)*rax*ray-(35.0D0/2.0D0)*rbz*rax*Cyz-(35.0D0/2.0D0)*rbz*ray*Cxz &
          !                      & -(5.0D0/2.0D0)*Cxz*Cyz))
          !T1(6)     = R6 * (rt3*(35.0D0*rbz*rax*ray+5.0D0*rax*Cyz+5.0D0*ray*Cxz))
          !T1(1)     = R6 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*ray-(5.0D0/2.0D0)*ray+5.0D0*rbz*Cyz))
          !T1(2)     = R6 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*rax-(5.0D0/2.0D0)*rax+5.0D0*rbz*Cxz))
          !T1(13)    = R5 * (rt3*(5.0D0*rbz*ray+Cyz))
          !T1(14)    = R5 * (rt3*(5.0D0*rbz*rax+Cxz))
          T0        = R5 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*rax*ray-(5.0D0/2.0D0)*rax*ray-5.0D0*rbz*rax*Cyz-5.0D0*rbz*ray*Cxz+Cxz*Cyz))
          T1(0)     = R7 * (rt3*((-315.0D0/4.0D0)*(rbz**(2.0D0))*rax*ray+(35.0D0/4.0D0)*rax*ray+(35.0D0/2.0D0)*rbz*rax*Cyz+(35.0D0/2.0D0)*rbz*ray*Cxz &
                                & -(5.0D0/2.0D0)*Cxz*Cyz))
          T1(6)     = R6 * (rt3*(35.0D0*rbz*rax*ray-5.0D0*rax*Cyz-5.0D0*ray*Cxz))
          T1(1)     = R6 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*ray-(5.0D0/2.0D0)*ray-5.0D0*rbz*Cyz))
          T1(2)     = R6 * (rt3*((35.0D0/2.0D0)*(rbz**(2.0D0))*rax-(5.0D0/2.0D0)*rax-5.0D0*rbz*Cxz))
          T1(13)    = R5 * (rt3*((-5.0D0)*rbz*ray+Cyz))
          T1(14)    = R5 * (rt3*((-5.0D0)*rbz*rax+Cxz))
!
!                      22s, 21c
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 4) .AND. (ib .EQ. 2) .AND. (bk .EQ. 1)) Then
          !T0        = R5 * (35.0D0*rbx*rbz*rax*ray+5.0D0*rbx*rax*Cyz+5.0D0*rbx*ray*Cxz+5.0D0*rbz*rax*Cyx+5.0D0*rbz*ray*Cxx+Cxx*Cyz+Cyx*Cxz)
          !T1(0)     = R7 * ((-315.0D0/2.0D0)*rbx*rbz*rax*ray-(35.0D0/2.0D0)*rbx*rax*Cyz-(35.0D0/2.0D0)*rbx*ray*Cxz-(35.0D0/2.0D0)*rbz*rax*Cyx-(35.0D0/2.0D0)*rbz*ray*Cxx &
          !                       & -(5.0D0/2.0D0)*Cxx*Cyz-(5.0D0/2.0D0)*Cyx*Cxz)
          !T1(4)     = R6 * (35.0D0*rbz*rax*ray+5.0D0*rax*Cyz+5.0D0*ray*Cxz)
          !T1(6)     = R6 * (35.0D0*rbx*rax*ray+5.0D0*rax*Cyx+5.0D0*ray*Cxx)
          !T1(1)     = R6 * (35.0D0*rbx*rbz*ray+5.0D0*rbx*Cyz+5.0D0*rbz*Cyx)
          !T1(2)     = R6 * (35.0D0*rbx*rbz*rax+5.0D0*rbx*Cxz+5.0D0*rbz*Cxx)
          !T1(7)     = R5 * (5.0D0*rbz*ray+Cyz)
          !T1(13)    = R5 * (5.0D0*rbx*ray+Cyx)
          !T1(8)     = R5 * (5.0D0*rbz*rax+Cxz)
          !T1(14)    = R5 * (5.0D0*rbx*rax+Cxx)
          T0        = R5 * (35.0D0*rbx*rbz*rax*ray-5.0D0*rbx*rax*Cyz-5.0D0*rbx*ray*Cxz-5.0D0*rbz*rax*Cyx-5.0D0*rbz*ray*Cxx+Cxx*Cyz+Cyx*Cxz)
          T1(0)     = R7 * ((-315.0D0/2.0D0)*rbx*rbz*rax*ray+(35.0D0/2.0D0)*rbx*rax*Cyz+(35.0D0/2.0D0)*rbx*ray*Cxz+(35.0D0/2.0D0)*rbz*rax*Cyx+(35.0D0/2.0D0)*rbz*ray*Cxx &
                                 & -(5.0D0/2.0D0)*Cxx*Cyz-(5.0D0/2.0D0)*Cyx*Cxz)
          T1(4)     = R6 * (35.0D0*rbz*rax*ray-5.0D0*rax*Cyz-5.0D0*ray*Cxz)
          T1(6)     = R6 * (35.0D0*rbx*rax*ray-5.0D0*rax*Cyx-5.0D0*ray*Cxx)
          T1(1)     = R6 * (35.0D0*rbx*rbz*ray-5.0D0*rbx*Cyz-5.0D0*rbz*Cyx)
          T1(2)     = R6 * (35.0D0*rbx*rbz*rax-5.0D0*rbx*Cxz-5.0D0*rbz*Cxx)
          T1(7)     = R5 * ((-5.0D0)*rbz*ray+Cyz)
          T1(13)    = R5 * ((-5.0D0)*rbx*ray+Cyx)
          T1(8)     = R5 * ((-5.0D0)*rbz*rax+Cxz)
          T1(14)    = R5 * ((-5.0D0)*rbx*rax+Cxx)
!
!                      22s, 21s
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 4) .AND. (ib .EQ. 2) .AND. (bk .EQ. 2)) Then
          !T0        = R5 * (35.0D0*rby*rbz*rax*ray+5.0D0*rby*rax*Cyz+5.0D0*rby*ray*Cxz+5.0D0*rbz*rax*Cyy+5.0D0*rbz*ray*Cxy+Cxy*Cyz+Cyy*Cxz)
          !T1(0)     = R7 * ((-315.0D0/2.0D0)*rby*rbz*rax*ray-(35.0D0/2.0D0)*rby*rax*Cyz-(35.0D0/2.0D0)*rby*ray*Cxz-(35.0D0/2.0D0)*rbz*rax*Cyy-(35.0D0/2.0D0)*rbz*ray*Cxy &
          !                       & -(5.0D0/2.0D0)*Cxy*Cyz-(5.0D0/2.0D0)*Cyy*Cxz)
          !T1(5)     = R6 * (35.0D0*rbz*rax*ray+5.0D0*rax*Cyz+5.0D0*ray*Cxz)
          !T1(6)     = R6 * (35.0D0*rby*rax*ray+5.0D0*rax*Cyy+5.0D0*ray*Cxy)
          !T1(1)     = R6 * (35.0D0*rby*rbz*ray+5.0D0*rby*Cyz+5.0D0*rbz*Cyy)
          !T1(2)     = R6 * (35.0D0*rby*rbz*rax+5.0D0*rby*Cxz+5.0D0*rbz*Cxy)
          !T1(10)    = R5 * (5.0D0*rbz*ray+Cyz)
          !T1(13)    = R5 * (5.0D0*rby*ray+Cyy)
          !T1(11)    = R5 * (5.0D0*rbz*rax+Cxz)
          !T1(14)    = R5 * (5.0D0*rby*rax+Cxy)
          T0        = R5 * (35.0D0*rby*rbz*rax*ray-5.0D0*rby*rax*Cyz-5.0D0*rby*ray*Cxz-5.0D0*rbz*rax*Cyy-5.0D0*rbz*ray*Cxy+Cxy*Cyz+Cyy*Cxz)
          T1(0)     = R7 * ((-315.0D0/2.0D0)*rby*rbz*rax*ray+(35.0D0/2.0D0)*rby*rax*Cyz+(35.0D0/2.0D0)*rby*ray*Cxz+(35.0D0/2.0D0)*rbz*rax*Cyy+(35.0D0/2.0D0)*rbz*ray*Cxy &
                                 & -(5.0D0/2.0D0)*Cxy*Cyz-(5.0D0/2.0D0)*Cyy*Cxz)
          T1(5)     = R6 * (35.0D0*rbz*rax*ray-5.0D0*rax*Cyz-5.0D0*ray*Cxz)
          T1(6)     = R6 * (35.0D0*rby*rax*ray-5.0D0*rax*Cyy-5.0D0*ray*Cxy)
          T1(1)     = R6 * (35.0D0*rby*rbz*ray-5.0D0*rby*Cyz-5.0D0*rbz*Cyy)
          T1(2)     = R6 * (35.0D0*rby*rbz*rax-5.0D0*rby*Cxz-5.0D0*rbz*Cxy)
          T1(10)    = R5 * ((-5.0D0)*rbz*ray+Cyz)
          T1(13)    = R5 * ((-5.0D0)*rby*ray+Cyy)
          T1(11)    = R5 * ((-5.0D0)*rbz*rax+Cxz)
          T1(14)    = R5 * ((-5.0D0)*rby*rax+Cxy)
!
!                      22s, 22c
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 4) .AND. (ib .EQ. 2) .AND. (bk .EQ. 3)) Then
          !T0        = R5 * ((35.0D0/2.0D0)*rax*ray*(rbx**(2.0D0))-(35.0D0/2.0D0)*rax*ray*(rby**(2.0D0))+5.0D0*rax*rbx*Cyx-5.0D0*rax*rby*Cyy+5.0D0*ray*rbx*Cxx-5.0D0*ray*rby*Cxy+Cxx*Cyx-Cxy*Cyy)
          !T1(0)     = R7 * ((-315.0D0/4.0D0)*rax*ray*(rbx**(2.0D0))+(315.0D0/4.0D0)*rax*ray*(rby**(2.0D0))-(35.0D0/2.0D0)*rax*rbx*Cyx+(35.0D0/2.0D0)*rax*rby*Cyy &
          !                      & -(35.0D0/2.0D0)*ray*rbx*Cxx+(35.0D0/2.0D0)*ray*rby*Cxy-(5.0D0/2.0D0)*Cxx*Cyx+(5.0D0/2.0D0)*Cxy*Cyy)
          !T1(1)     = R6 * ((35.0D0/2.0D0)*ray*(rbx**(2.0D0))-(35.0D0/2.0D0)*ray*(rby**(2.0D0))+5.0D0*rbx*Cyx-5.0D0*rby*Cyy)
          !T1(2)     = R6 * ((35.0D0/2.0D0)*rax*(rbx**(2.0D0))-(35.0D0/2.0D0)*rax*(rby**(2.0D0))+5.0D0*rbx*Cxx-5.0D0*rby*Cxy)
          !T1(4)     = R6 * (35.0D0*rax*ray*rbx+5.0D0*rax*Cyx+5.0D0*ray*Cxx)
          !T1(5)     = R6 * (-35.0D0*rax*ray*rby-5.0D0*rax*Cyy-5.0D0*ray*Cxy)
          !T1(7)     = R5 * (5.0D0*ray*rbx+Cyx)
          !T1(8)     = R5 * (5.0D0*rax*rbx+Cxx)
          !T1(10)    = R5 * (-5.0D0*ray*rby-Cyy)
          !T1(11)    = R5 * (-5.0D0*rax*rby-Cxy)
          T0        = R5 * ((35.0D0/2.0D0)*rax*ray*(rbx**(2.0D0))-(35.0D0/2.0D0)*rax*ray*(rby**(2.0D0))-5.0D0*rax*rbx*Cyx+5.0D0*rax*rby*Cyy-5.0D0*ray*rbx*Cxx+5.0D0*ray*rby*Cxy+Cxx*Cyx-Cxy*Cyy)
          T1(0)     = R7 * ((-315.0D0/4.0D0)*rax*ray*(rbx**(2.0D0))+(315.0D0/4.0D0)*rax*ray*(rby**(2.0D0))+(35.0D0/2.0D0)*rax*rbx*Cyx-(35.0D0/2.0D0)*rax*rby*Cyy &
                                & +(35.0D0/2.0D0)*ray*rbx*Cxx-(35.0D0/2.0D0)*ray*rby*Cxy-(5.0D0/2.0D0)*Cxx*Cyx+(5.0D0/2.0D0)*Cxy*Cyy)
          T1(1)     = R6 * ((35.0D0/2.0D0)*ray*(rbx**(2.0D0))-(35.0D0/2.0D0)*ray*(rby**(2.0D0))-5.0D0*rbx*Cyx+5.0D0*rby*Cyy)
          T1(2)     = R6 * ((35.0D0/2.0D0)*rax*(rbx**(2.0D0))-(35.0D0/2.0D0)*rax*(rby**(2.0D0))-5.0D0*rbx*Cxx+5.0D0*rby*Cxy)
          T1(4)     = R6 * (35.0D0*rax*ray*rbx-5.0D0*rax*Cyx-5.0D0*ray*Cxx)
          T1(5)     = R6 * (-35.0D0*rax*ray*rby+5.0D0*rax*Cyy+5.0D0*ray*Cxy)
          T1(7)     = R5 * ((-5.0D0)*ray*rbx+Cyx)
          T1(8)     = R5 * ((-5.0D0)*rax*rbx+Cxx)
          T1(10)    = R5 * (5.0D0*ray*rby-Cyy)
          T1(11)    = R5 * (5.0D0*rax*rby-Cxy)
!
!                      22s, 22s
!
        Else If ((ia .EQ. 2) .AND. (ak .EQ. 4) .AND. (ib .EQ. 2) .AND. (bk .EQ. 4)) Then
          !T0        = R5 * (35.0D0*rax*ray*rbx*rby+5.0D0*rax*rbx*Cyy+5.0D0*rax*rby*Cyx+5.0D0*ray*rbx*Cxy+5.0D0*ray*rby*Cxx+Cxx*Cyy+Cxy*Cyx)
          !T1(0)     = R7 * ((-315.0D0/2.0D0)*rax*ray*rbx*rby-(35.0D0/2.0D0)*rax*rbx*Cyy-(35.0D0/2.0D0)*rax*rby*Cyx-(35.0D0/2.0D0)*ray*rbx*Cxy-(35.0D0/2.0D0)*ray*rby*Cxx &
          !                       & -(5.0D0/2.0D0)*Cxx*Cyy-(5.0D0/2.0D0)*Cxy*Cyx)
          !T1(1)     = R6 * (35.0D0*ray*rbx*rby+5.0D0*rbx*Cyy+5.0D0*rby*Cyx)
          !T1(2)     = R6 * (35.0D0*rax*rbx*rby+5.0D0*rbx*Cxy+5.0D0*rby*Cxx)
          !T1(4)     = R6 * (35.0D0*rax*ray*rby+5.0D0*rax*Cyy+5.0D0*ray*Cxy)
          !T1(5)     = R6 * (35.0D0*rax*ray*rbx+5.0D0*rax*Cyx+5.0D0*ray*Cxx)
          !T1(7)     = R5 * (5.0D0*ray*rby+Cyy)
          !T1(8)     = R5 * (5.0D0*rax*rby+Cxy)
          !T1(10)    = R5 * (5.0D0*ray*rbx+Cyx)
          !T1(11)    = R5 * (5.0D0*rax*rbx+Cxx)
          T0        = R5 * (35.0D0*rax*ray*rbx*rby-5.0D0*rax*rbx*Cyy-5.0D0*rax*rby*Cyx-5.0D0*ray*rbx*Cxy-5.0D0*ray*rby*Cxx+Cxx*Cyy+Cxy*Cyx)
          T1(0)     = R7 * ((-315.0D0/2.0D0)*rax*ray*rbx*rby+(35.0D0/2.0D0)*rax*rbx*Cyy+(35.0D0/2.0D0)*rax*rby*Cyx+(35.0D0/2.0D0)*ray*rbx*Cxy+(35.0D0/2.0D0)*ray*rby*Cxx &
                                 & -(5.0D0/2.0D0)*Cxx*Cyy-(5.0D0/2.0D0)*Cxy*Cyx)
          T1(1)     = R6 * (35.0D0*ray*rbx*rby-5.0D0*rbx*Cyy-5.0D0*rby*Cyx)
          T1(2)     = R6 * (35.0D0*rax*rbx*rby-5.0D0*rbx*Cxy-5.0D0*rby*Cxx)
          T1(4)     = R6 * (35.0D0*rax*ray*rby-5.0D0*rax*Cyy-5.0D0*ray*Cxy)
          T1(5)     = R6 * (35.0D0*rax*ray*rbx-5.0D0*rax*Cyx-5.0D0*ray*Cxx)
          T1(7)     = R5 * ((-5.0D0)*ray*rby+Cyy)
          T1(8)     = R5 * ((-5.0D0)*rax*rby+Cxy)
          T1(10)    = R5 * ((-5.0D0)*ray*rbx+Cyx)
          T1(11)    = R5 * ((-5.0D0)*rax*rbx+Cxx)
         
        End If

  end subroutine
