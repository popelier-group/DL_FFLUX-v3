Subroutine Bohr2Ang(rcut,xxt,yyt,zzt,rrt,ener,vir)

  Use kinds_f90
  Use config_module, Only : xxx,yyy,zzz, &
                          & fxx,fyy,fzz
  Use fflux_module,  Only : conv_B2A,conv_Ha_2_10Jmol,conv_Fau

  Implicit None

  Real( Kind = wp ) :: rcut,xxt,yyt,zzt,rrt,ener,vir

  xxx(:) = xxx(:) * conv_B2A
  yyy(:) = yyy(:) * conv_B2A
  zzz(:) = zzz(:) * conv_B2A

  fxx(:) = fxx(:) * conv_Fau
  fyy(:) = fyy(:) * conv_Fau
  fzz(:) = fzz(:) * conv_Fau

  ener = ener * conv_Ha_2_10Jmol

  vir = vir * conv_Ha_2_10Jmol

  rcut = rcut * conv_B2A

End Subroutine


Subroutine Ang2Bohr(rcut,xxt,yyt,zzt,rrt,ener,vir)

  Use kinds_f90
  Use config_module, Only : xxx,yyy,zzz, &
                          & fxx,fyy,fzz
  Use fflux_module,  Only : conv_B2A,conv_Ha_2_10Jmol,conv_Fau

  Implicit None
  
  Real( Kind = wp ) :: rcut,xxt,yyt,zzt,rrt,ener,vir,InvB2A,InvFau,InvHa_2_10Jmol

  !InvB2A = 1.0_wp / conv_B2A
  !InvFau = 1.0_wp / conv_Fau
  !InvHa_2_10Jmol = 1.0_wp / conv_Ha_2_10Jmol

  xxx(:) = xxx(:) / conv_B2A
  yyy(:) = yyy(:) / conv_B2A
  zzz(:) = zzz(:) / conv_B2A

  fxx(:) = fxx(:) / conv_Fau
  fyy(:) = fyy(:) / conv_Fau
  fzz(:) = fzz(:) / conv_Fau

  ener = ener / conv_Ha_2_10Jmol

  vir = vir / conv_Ha_2_10Jmol

  rcut = rcut / conv_B2A

End Subroutine

!=======================================================

Subroutine Atm2Dlp(rcut,enr1,vir1,enr2,vir2)

! Subroutine to convert from atomic to dlpoly units 

  Use kinds_f90
  Use config_module, Only : xxx,yyy,zzz, &
                          & fxx,fyy,fzz,cell
  Use fflux_module,  Only : conv_B2A,conv_Ha_2_10Jmol,conv_Fau

  Implicit None

  Real( Kind = wp ) :: rcut,enr1,vir1,enr2,vir2

  xxx(:) = xxx(:) * conv_B2A
  yyy(:) = yyy(:) * conv_B2A
  zzz(:) = zzz(:) * conv_B2A

  fxx(:) = fxx(:) * conv_Fau
  fyy(:) = fyy(:) * conv_Fau
  fzz(:) = fzz(:) * conv_Fau

  enr1 = enr1 * conv_Ha_2_10Jmol
  enr2 = enr2 * conv_Ha_2_10Jmol

  vir1 = vir1 * conv_Ha_2_10Jmol
  vir2 = vir2 * conv_Ha_2_10Jmol

  rcut = rcut * conv_B2A

  cell(:) = cell(:) * conv_B2A

End Subroutine

!=======================================================

Subroutine Dlp2Atm(rcut,enr1,vir1,enr2,vir2)

! Subroutine to convert from dlpoly to atomic units 

  Use kinds_f90
  Use config_module, Only : xxx,yyy,zzz, &
                          & fxx,fyy,fzz,cell
  Use fflux_module,  Only : conv_B2A,conv_Ha_2_10Jmol,conv_Fau

  Implicit None
  
  Real( Kind = wp ) :: rcut,enr1,vir1,enr2,vir2

  xxx(:) = xxx(:) / conv_B2A
  yyy(:) = yyy(:) / conv_B2A
  zzz(:) = zzz(:) / conv_B2A

  fxx(:) = fxx(:) / conv_Fau
  fyy(:) = fyy(:) / conv_Fau
  fzz(:) = fzz(:) / conv_Fau

  enr1 = enr1 / conv_Ha_2_10Jmol
  enr2 = enr2 / conv_Ha_2_10Jmol

  vir1 = vir1 / conv_Ha_2_10Jmol
  vir2 = vir2 / conv_Ha_2_10Jmol

  rcut = rcut / conv_B2A

  cell(:) = cell(:) / conv_B2A

End Subroutine
!=======================================================

Subroutine Atm2Dlp_pos()

! Subroutine to convert from atomic to dlpoly units 
! Take no vars - just convert necessary arrays - BCBS 

  Use kinds_f90
  Use config_module, Only : xxx,yyy,zzz,cell
  Use fflux_module,  Only : conv_B2A,conv_Ha_2_10Jmol,conv_Fau

  Implicit None

  xxx(:) = xxx(:) * conv_B2A
  yyy(:) = yyy(:) * conv_B2A
  zzz(:) = zzz(:) * conv_B2A

  cell(:) = cell(:) * conv_B2A

End Subroutine

!=======================================================

Subroutine Dlp2Atm_pos()

! Subroutine to convert from dlpoly to atomic units
! Take no vars - just convert necessary arrays - BCBS 

  Use kinds_f90
  Use config_module, Only : xxx,yyy,zzz,cell
  Use fflux_module,  Only : conv_B2A,conv_Ha_2_10Jmol,conv_Fau

  Implicit None

  xxx(:) = xxx(:) / conv_B2A
  yyy(:) = yyy(:) / conv_B2A
  zzz(:) = zzz(:) / conv_B2A

  cell(:) = cell(:) / conv_B2A

End Subroutine
