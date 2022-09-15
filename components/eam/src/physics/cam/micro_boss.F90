module micro_boss

!--------------------------------------------------------------------------
!
! This module contains process rates calculated by the BOSS microphysics
! framework (two-category version), originally developed by Hugh Morrison, Sean
! Patrick Santos, and Marcus van Lier-Walqui, and re-implemented here for use in
! E3SM by Sean Patrick Santos.
!
!--------------------------------------------------------------------------

use micro_mg_utils, only: r8, qsmall, pi, rhow, var_coef

implicit none
private

public :: &
     micro_boss_init, &
     boss_thermo, &
     boss_m, &
     boss_rain_evap, &
     boss_cloud_selfcol, &
     boss_auto, &
     boss_accr, &
     boss_rain_selfcol, &
     boss_cloud_fall, &
     boss_rain_fall

! Conversions between liquid mass mixing ratio (ql) and third moment (M3).
real(r8), parameter, public :: m3_to_ql = rhow * pi / 6.
real(r8), parameter, public :: ql_to_m3 = 1. / m3_to_ql


! BOSS tunable parameters
real(r8) :: rain_evap_a
real(r8) :: rain_evap_b
real(r8) :: cloud_selfcol_a
real(r8) :: cloud_selfcol_b
real(r8) :: auto_cloud_size
real(r8) :: auto_rain_size
real(r8) :: auto_non_rain_term_a
real(r8) :: auto_non_rain_term_b
real(r8) :: auto_rain_term_a
real(r8) :: auto_rain_term_b_cloud
real(r8) :: auto_rain_term_b_rain
real(r8) :: auto_rain_term_b_num
real(r8) :: accr_a
real(r8) :: accr_b_cloud
real(r8) :: accr_b_rain
real(r8) :: rain_selfcol_a
real(r8) :: rain_selfcol_b
real(r8) :: cloud_fall0_a
real(r8) :: cloud_fall0_b
real(r8) :: cloud_fall3_a
real(r8) :: cloud_fall3_b
real(r8) :: rain_fall0_a
real(r8) :: rain_fall0_b
real(r8) :: rain_fall3_a
real(r8) :: rain_fall3_b

contains

! Initialize BOSS parameters.
subroutine micro_boss_init(rain_evap_a_in, rain_evap_b_in, &
     cloud_selfcol_a_in, cloud_selfcol_b_in, &
     auto_cloud_size_in, auto_rain_size_in, &
     auto_non_rain_term_a_in, auto_non_rain_term_b_in, &
     auto_rain_term_a_in, auto_rain_term_b_cloud_in, &
     auto_rain_term_b_rain_in, auto_rain_term_b_num_in, &
     accr_a_in, accr_b_cloud_in, accr_b_rain_in, &
     rain_selfcol_a_in, rain_selfcol_b_in, cloud_fall0_a_in, &
     cloud_fall0_b_in, cloud_fall3_a_in, cloud_fall3_b_in, &
     rain_fall0_a_in, rain_fall0_b_in, rain_fall3_a_in, &
     rain_fall3_b_in)

  real(r8), intent(in) :: rain_evap_a_in
  real(r8), intent(in) :: rain_evap_b_in
  real(r8), intent(in) :: cloud_selfcol_a_in
  real(r8), intent(in) :: cloud_selfcol_b_in
  real(r8), intent(in) :: auto_cloud_size_in
  real(r8), intent(in) :: auto_rain_size_in
  real(r8), intent(in) :: auto_non_rain_term_a_in
  real(r8), intent(in) :: auto_non_rain_term_b_in
  real(r8), intent(in) :: auto_rain_term_a_in
  real(r8), intent(in) :: auto_rain_term_b_cloud_in
  real(r8), intent(in) :: auto_rain_term_b_rain_in
  real(r8), intent(in) :: auto_rain_term_b_num_in
  real(r8), intent(in) :: accr_a_in
  real(r8), intent(in) :: accr_b_cloud_in
  real(r8), intent(in) :: accr_b_rain_in
  real(r8), intent(in) :: rain_selfcol_a_in
  real(r8), intent(in) :: rain_selfcol_b_in
  real(r8), intent(in) :: cloud_fall0_a_in
  real(r8), intent(in) :: cloud_fall0_b_in
  real(r8), intent(in) :: cloud_fall3_a_in
  real(r8), intent(in) :: cloud_fall3_b_in
  real(r8), intent(in) :: rain_fall0_a_in
  real(r8), intent(in) :: rain_fall0_b_in
  real(r8), intent(in) :: rain_fall3_a_in
  real(r8), intent(in) :: rain_fall3_b_in

  rain_evap_a = rain_evap_a_in
  rain_evap_b = rain_evap_b_in
  cloud_selfcol_a = cloud_selfcol_a_in
  cloud_selfcol_b = cloud_selfcol_b_in
  auto_cloud_size = auto_cloud_size_in
  auto_rain_size = auto_rain_size_in
  auto_non_rain_term_a = auto_non_rain_term_a_in
  auto_non_rain_term_b = auto_non_rain_term_b_in
  auto_rain_term_a = auto_rain_term_a_in
  auto_rain_term_b_cloud = auto_rain_term_b_cloud_in
  auto_rain_term_b_rain = auto_rain_term_b_rain_in
  auto_rain_term_b_num = auto_rain_term_b_num_in
  accr_a = accr_a_in
  accr_b_cloud = accr_b_cloud_in
  accr_b_rain = accr_b_rain_in
  rain_selfcol_a = rain_selfcol_a_in
  rain_selfcol_b = rain_selfcol_b_in
  cloud_fall0_a = cloud_fall0_a_in
  cloud_fall0_b = cloud_fall0_b_in
  cloud_fall3_a = cloud_fall3_a_in
  cloud_fall3_b = cloud_fall3_b_in
  rain_fall0_a = rain_fall0_a_in
  rain_fall0_b = rain_fall0_b_in
  rain_fall3_a = rain_fall3_a_in
  rain_fall3_b = rain_fall3_b_in

end subroutine micro_boss_init

! Calculate thermodynamic (super/sub-saturation) factor used by BOSS.
elemental function boss_thermo(rho, t, qv, qvs, lv, dv)
  use micro_mg_utils, only: calc_ab
  ! Density of air (kg/m^3)
  real(r8), intent(in) :: rho
  ! Temperature (K)
  real(r8), intent(in) :: t
  ! Specific humidity (kg/kg)
  real(r8), intent(in) :: qv
  ! Saturation specific humidity (kg/kg)
  real(r8), intent(in) :: qvs
  ! Latent heat of vaporization (J/kg)
  real(r8), intent(in) :: lv
  ! Water vapor diffusivity (m^2/s)
  real(r8), intent(in) :: dv
  ! Thermodynamic factor (kg/m/s)
  real(r8) :: boss_thermo

  ! correction factor due to latent heat
  real(r8) :: ab

  ab = calc_ab(t, qvs, lv)

  boss_thermo = 2._r8 * pi * rho * dv * (qv/qvs - 1._r8) * qvs / ab

end function boss_thermo

! Calculate m = M3/M0 (proportional to mean particle volume or mass).
elemental function boss_m(ql, nl)
  ! Drop mass (kg/kg)
  real(r8), intent(in) :: ql
  ! Drop number (#/kg)
  real(r8), intent(in) :: nl
  ! Mean size "m" (m^3)
  real(r8) :: boss_m

  if (ql >= qsmall) then
     boss_m = ql_to_m3 * ql / nl
  else
     boss_m = 0._r8
  end if

end function boss_m

! Calculate rain evaporation rate.
subroutine boss_rain_evap(thermo, qr, nr, mr, dqr_dt)
  ! Thermodynamic factor (kg/m/s)
  real(r8), intent(in) :: thermo(:)
  ! Rain mass (kg/kg)
  real(r8), intent(in) :: qr(:)
  ! Rain number (#/kg)
  real(r8), intent(in) :: nr(:)
  ! Mean rain size (m^3)
  real(r8), intent(in) :: mr(:)
  ! Rain mass tendency (kg/kg/s)
  real(r8), intent(out) :: dqr_dt(:)

  ! Rain M3 tendency (m^3/kg/s)
  real(r8) :: dm3_dt(size(dqr_dt))

  where (qr >= qsmall)
     dm3_dt = rain_evap_a * thermo * nr * mr**rain_evap_b
  elsewhere
     dm3_dt = 0._r8
  end where

  dqr_dt = m3_to_ql * dm3_dt

end subroutine boss_rain_evap

! Calculate cloud self-collection rate.
subroutine boss_cloud_selfcol(qc, nc, mc, relvar, dnc_dt)
  ! Cloud mass (kg/kg)
  real(r8), intent(in) :: qc(:)
  ! Cloud number (#/kg)
  real(r8), intent(in) :: nc(:)
  ! Mean cloud size (m^3)
  real(r8), intent(in) :: mc(:)
  ! Relative variance in cloud water (unitless)
  real(r8), intent(in) :: relvar(:)
  ! Cloud number tendency (#/kg/s)
  real(r8), intent(out) :: dnc_dt(:)

  ! Factor accounting for effects of relvar (unitless)
  real(r8) :: relvar_coef(size(dnc_dt))

  where (qc >= qsmall)
     relvar_coef = var_coef(relvar, cloud_selfcol_b)
     dnc_dt = relvar_coef * cloud_selfcol_a * nc**2 * mc**cloud_selfcol_b
  elsewhere
     dnc_dt = 0._r8
  end where

end subroutine boss_cloud_selfcol

! Calculate cloud autoconversion rate.
subroutine boss_auto(qc, nc, mc, qr, nr, mr, relvar, dql_dt, dnc_dt, dnr_dt)
  ! Cloud mass (kg/kg)
  real(r8), intent(in) :: qc(:)
  ! Cloud number (#/kg)
  real(r8), intent(in) :: nc(:)
  ! Mean cloud size (m^3)
  real(r8), intent(in) :: mc(:)
  ! Rain mass (kg/kg)
  real(r8), intent(in) :: qr(:)
  ! Rain number (#/kg)
  real(r8), intent(in) :: nr(:)
  ! Mean rain size (m^3)
  real(r8), intent(in) :: mr(:)
  ! Relative variance in cloud water (unitless)
  real(r8), intent(in) :: relvar(:)
  ! Mass tendency (#/kg/s)
  real(r8), intent(out) :: dql_dt(:)
  ! Cloud number tendency (#/kg/s)
  real(r8), intent(out) :: dnc_dt(:)
  ! Rain number tendency (#/kg/s)
  real(r8), intent(out) :: dnr_dt(:)

  ! Factor accounting for effects of relvar (unitless)
  real(r8) :: relvar_coef(size(dql_dt))
  ! Non-rain-dependent M3 tendency term (m^3/kg/s)
  real(r8) :: dm3_dt_non_rain_term(size(dql_dt))
  ! Rain-dependent M3 tendency term (m^3/kg/s)
  real(r8) :: dm3_dt_rain_term(size(dql_dt))
  ! Total M3 tendency (m^3/kg/s)
  real(r8) :: dm3_dt(size(dql_dt))

  where (qc >= qsmall)
     relvar_coef = var_coef(relvar, auto_non_rain_term_b)
     dm3_dt_non_rain_term = relvar_coef * auto_non_rain_term_a * nc**2 &
          * mc**auto_non_rain_term_b
     where (qr >= qsmall)
        relvar_coef = var_coef(relvar, auto_rain_term_b_cloud)
        dm3_dt_rain_term = relvar_coef * auto_rain_term_a * nc**2 &
             * mc**auto_rain_term_b_cloud * mr**auto_rain_term_b_rain &
             * (nr/nc)**auto_rain_term_b_num
     elsewhere
        dm3_dt_rain_term = 0._r8
     end where
     dm3_dt = dm3_dt_non_rain_term + dm3_dt_rain_term
  elsewhere
     dm3_dt = 0._r8
  end where

  dql_dt = m3_to_ql * dm3_dt
  dnc_dt = auto_cloud_size * dm3_dt
  dnr_dt = auto_rain_size * dm3_dt

end subroutine boss_auto

! Calculate cloud accretion rate.
subroutine boss_accr(qc, nc, mc, qr, nr, mr, relvar, dql_dt, dnc_dt)
  ! Cloud mass (kg/kg)
  real(r8), intent(in) :: qc(:)
  ! Cloud number (#/kg)
  real(r8), intent(in) :: nc(:)
  ! Mean cloud size (m^3)
  real(r8), intent(in) :: mc(:)
  ! Rain mass (kg/kg)
  real(r8), intent(in) :: qr(:)
  ! Rain number (#/kg)
  real(r8), intent(in) :: nr(:)
  ! Mean rain size (m^3)
  real(r8), intent(in) :: mr(:)
  ! Relative variance in cloud water (unitless)
  real(r8), intent(in) :: relvar(:)
  ! Mass tendency (#/kg/s)
  real(r8), intent(out) :: dql_dt(:)
  ! Cloud number tendency (#/kg/s)
  real(r8), intent(out) :: dnc_dt(:)

  ! Factor accounting for effects of relvar (unitless)
  real(r8) :: relvar_coef(size(dql_dt))

  ! M3 tendency (m^3/kg/s)
  real(r8) :: dm3_dt(size(dql_dt))

  where (qc >= qsmall .and. qr >= qsmall)
     relvar_coef = var_coef(relvar, accr_b_cloud)
     dm3_dt = relvar_coef * accr_a * nc * nr * mc**accr_b_cloud &
          * mr**accr_b_rain
     ! Set this inside the where statement to avoid floating point exceptions
     ! when mc is zero.
     dnc_dt = dm3_dt / mc
  elsewhere
     dm3_dt = 0._r8
     dnc_dt = 0._r8
  end where

  dql_dt = m3_to_ql * dm3_dt

end subroutine boss_accr

! Calculate rain self-collection rate.
subroutine boss_rain_selfcol(qr, nr, mr, dnr_dt)
  ! Rain mass (kg/kg)
  real(r8), intent(in) :: qr(:)
  ! Rain number (#/kg)
  real(r8), intent(in) :: nr(:)
  ! Mean rain size (m^3)
  real(r8), intent(in) :: mr(:)
  ! Rain number tendency (#/kg/s)
  real(r8), intent(out) :: dnr_dt(:)

  where (qr >= qsmall)
     dnr_dt = rain_selfcol_a * nr**2 * mr**rain_selfcol_b
  elsewhere
     dnr_dt = 0._r8
  end where

end subroutine boss_rain_selfcol

! Calculate cloud fall speed.
subroutine boss_cloud_fall(stokes_fac, qc, mc, relvar, vq, vn)
  ! Scaling factor calculated using Stokes' law. (1/(m*s)
  real(r8), intent(in) :: stokes_fac
  ! Cloud mass (kg/kg)
  real(r8), intent(in) :: qc
  ! Mean cloud size (m^3)
  real(r8), intent(in) :: mc
  ! Relative variance in cloud water (unitless)
  real(r8), intent(in) :: relvar
  ! Cloud mass fall speed (m/s)
  real(r8), intent(out) :: vq
  ! Cloud number fall speed (m/s)
  real(r8), intent(out) :: vn

  ! Factor accounting for effects of relvar (unitless)
  real(r8) :: relvar_coef

  if (qc >= qsmall) then
     ! Adding 1 here is due to the fact that the actual flux out the bottom of
     ! the grid cell is proportional to vq * qc (or equivalently, higher fall
     ! speeds are correlated with higher mass).
     relvar_coef = var_coef(relvar, cloud_fall3_b + 1.)
     vq = relvar_coef * stokes_fac * cloud_fall3_a * mc**cloud_fall3_b
     relvar_coef = var_coef(relvar, cloud_fall0_b)
     vn = relvar_coef * stokes_fac * cloud_fall0_a * mc**cloud_fall0_b
  else
     vq = 0._r8
     vn = 0._r8
  end if

end subroutine boss_cloud_fall

! Calculate rain fall speed.
subroutine boss_rain_fall(rhofac, qr, mr, vq, vn)
  ! Density correction factor (unitless)
  real(r8), intent(in) :: rhofac
  ! Rain mass (kg/kg)
  real(r8), intent(in) :: qr
  ! Mean rain size (m^3)
  real(r8), intent(in) :: mr
  ! Rain mass fall speed (m/s)
  real(r8), intent(out) :: vq
  ! Rain number fall speed (m/s)
  real(r8), intent(out) :: vn

  if (qr >= qsmall) then
     vq = rhofac * rain_fall3_a * mr**rain_fall3_b
     vn = rhofac * rain_fall0_a * mr**rain_fall0_b
  else
     vq = 0._r8
     vn = 0._r8
  end if

end subroutine boss_rain_fall

end module micro_boss
