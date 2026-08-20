module THC_EOS_Utils_module

! ----------------------------------------------------------------------------
! THC EOS utility routines.
!
! Thin wrappers around the EOS_Water_module extended density and viscosity
! interfaces (procedure pointers set in THCSetup), plus the Somerton
! effective-thermal-conductivity mixing model, providing the property +
! derivative set THC needs for its 3x3 (P,T,C) Jacobian.
!
!   * Density  -> EOSWaterDensityExt
!   * Viscosity-> EOSWaterViscosityExt
!   * The EOS returns MOLAR-density derivative (dwp/dwt are d(rho_kmol)/dP,/dT);
!     they are converted here to MASS-density derivatives via FMWH2O so the
!     outputs match the documented auxvar units kg/(m^3 . X)
!   * d(rho)/dC and d(mu)/dC are NOT exposed by the EOS.
!     They are formed locally by the chain rule  d/dC = d/ds . ds/dC, where the
!     salinity derivative d/ds is obtained by a one-sided finite difference.
!   * No porous-medium thermal-conductivity mixing exists in the EOS
!     Somerton (1974) mixing is implemented here.
!
! Author: Piyoosh Jaysaval
! Date:   06/09/26
! ----------------------------------------------------------------------------

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module, only : FMWH2O
  use EOS_Water_module, only : EOSWaterDensityExt, &
                               EOSWaterViscosityExt, &
                               EOSWaterDensity, &
                               EOSWaterViscosity, &
                               EOSWaterDensityIsBrine, &
                               EOSWaterSaturationPressure

  implicit none

  private

  ! --- module parameters --------------------------------------------------
  ! Effective molar mass of the dissolved solids [kg/mol] (NaCl-equivalent;
  ! FMWNACL = 58.44277 kg/kmol).
  PetscReal, parameter, public :: thc_molar_mass_solute = 0.05844d0
  ! Reference liquid density [kg/m^3] used to convert C [mol/L] to the salinity
  ! mass fraction the Ext EOS routines require. The dependence is weak for the
  ! dilute-to-moderate TDS range; callers needing tight consistency may iterate
  ! using the current-iteration density.
  PetscReal, parameter, public :: thc_density_reference = 1000.d0
  ! Constituent thermal conductivities for the Somerton mixing model [W/(m.K)].
  PetscReal, parameter, public :: thc_kappa_liquid = 0.6d0
  PetscReal, parameter, public :: thc_kappa_air    = 0.025d0

  ! Finite-difference step in salinity mass fraction for the d/dC chain rule.
  PetscReal, parameter :: thc_salinity_perturbation = 1.d-6
  ! Floor on saturation to keep the Somerton derivative finite at S_l -> 0.
  PetscReal, parameter, public :: thc_sat_floor = 1.d-6

  public :: THCConcToMassFraction, &
            THCDensityAndDerivs, &
            THCViscosityAndDerivs, &
            THCThermalConductivityEff

contains

! ****************************************************************************
subroutine THCConcToMassFraction(C, den_kg, M_s, s, ds_dC)
  !
  ! Convert solute concentration C [mol/L] to a salinity mass fraction
  ! s [kg/kg] (the aux(1) input required by the extended EOS routines) and
  ! return ds/dC for the chain rule.
  !
  !   mass_solute_per_m3 = C * M_s * 1000   ! mol/L -> mol/m^3 -> kg/m^3
  !   s     = mass_solute_per_m3 / den_kg   ! [kg/kg]
  !   ds_dC = M_s * 1000 / den_kg           ! [kg/kg / (mol/L)]
  !
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !
  implicit none

  PetscReal, intent(in)  :: C       ! concentration [mol/L]
  PetscReal, intent(in)  :: den_kg  ! liquid density [kg/m^3]
  PetscReal, intent(in)  :: M_s     ! solute molar mass [kg/mol]
  PetscReal, intent(out) :: s       ! salinity mass fraction [kg/kg]
  PetscReal, intent(out) :: ds_dC   ! d(s)/d(C) [kg/kg / (mol/L)]

  ds_dC = M_s * 1.d3 / den_kg
  s     = C * ds_dC

end subroutine THCConcToMassFraction

! ****************************************************************************
subroutine THCDensityAndDerivs(T, P, C, den_kg, den_kmol, &
                                   dden_dp, dden_dT, dden_dC, ierr)
  !
  ! Liquid density rho_l(P,T,C) and all derivatives needed for the THC
  ! Jacobian, via EOSWaterDensityExt.
  !
  ! The EOS returns dwp/dwt = d(rho_kmol)/dP,/dT; these are
  ! converted to MASS-density derivatives via rho_kg = rho_kmol * FMWH2O.
  ! d(rho)/dC is formed locally (the EOS exposes no salinity derivative,
  ! d(rho)/dC = d(rho)/ds . ds/dC, with d(rho)/ds from a
  ! one-sided finite difference on the salinity mass fraction.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !
  implicit none

  PetscReal,      intent(in)    :: T         ! temperature [C]
  PetscReal,      intent(in)    :: P         ! pressure [Pa]
  PetscReal,      intent(in)    :: C         ! concentration [mol/L]
  PetscReal,      intent(out)   :: den_kg    ! rho_l [kg/m^3]
  PetscReal,      intent(out)   :: den_kmol  ! rho_l [kmol/m^3]
  PetscReal,      intent(out)   :: dden_dp   ! d(rho_kg)/dP [kg/(m^3.Pa)]
  PetscReal,      intent(out)   :: dden_dT   ! d(rho_kg)/dT [kg/(m^3.C)]
  PetscReal,      intent(out)   :: dden_dC   ! d(rho_kg)/dC [kg/(m^3.(mol/L))]
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: aux1(1)
  PetscReal :: dw, dwmol, dwp, dwt
  PetscReal :: dw_pert, dwmol_pert, dwp_pert, dwt_pert
  PetscReal :: s, ds_dC, dden_ds

  if (.not.EOSWaterDensityIsBrine()) then
    ! Pure-water correlation selected (IF97, IFC67, ...)
    call EOSWaterDensity(T, P, dw, dwmol, dwp, dwt, ierr)
    den_kg   = dw
    den_kmol = dwmol
    dden_dp  = dwp * FMWH2O
    dden_dT  = dwt * FMWH2O
    dden_dC  = 0.d0
    return
  endif

  ! C [mol/L] -> salinity mass fraction s [kg/kg]; seed with the reference
  ! density, then one Picard pass so s is consistent with the EOS density
  call THCConcToMassFraction(C, thc_density_reference, &
                                 thc_molar_mass_solute, s, ds_dC)
  aux1(1) = s
  call EOSWaterDensityExt(T, P, aux1, dw, dwmol, dwp, dwt, ierr)
  call THCConcToMassFraction(C, dw, thc_molar_mass_solute, s, ds_dC)
  aux1(1) = s

  ! Density with analytic P,T derivatives (molar-density basis).
  ! ds_dC neglects the implicit d(rho)/dC term (second order in s).
  call EOSWaterDensityExt(T, P, aux1, dw, dwmol, dwp, dwt, ierr)
  den_kg   = dw
  den_kmol = dwmol

  ! Convert molar-density derivatives to mass-density derivatives.
  dden_dp = dwp * FMWH2O
  dden_dT = dwt * FMWH2O

  ! d(rho)/dC by chain rule, d(rho)/ds via one-sided finite difference.
  aux1(1) = s + thc_salinity_perturbation
  call EOSWaterDensityExt(T, P, aux1, dw_pert, dwmol_pert, dwp_pert, &
                          dwt_pert, ierr)
  dden_ds = (dw_pert - dw) / thc_salinity_perturbation
  dden_dC = dden_ds * ds_dC

end subroutine THCDensityAndDerivs

! ****************************************************************************
subroutine THCViscosityAndDerivs(T, P, C, den_kg, vis, dvis_dT, dvis_dC, &
                                     ierr)
  !
  ! Liquid dynamic viscosity mu_l(T,C) and the derivatives needed for the
  ! THC Jacobian, via EOSWaterViscosityExt.
  !
  ! The viscosity interface requires saturation pressure inputs (PS, dPS_dT).
  ! Batzle & Wang Eq 32 does not actually use them, but they
  ! are computed here for correctness / pointer-agnosticism.
  !
  ! Pressure dependence is omitted by design (dVW_dP == 0), so no dvis_dP is
  ! returned. d(mu)/dC is formed locally:
  !   d(mu)/dC = d(mu)/ds . ds/dC, d(mu)/ds via one-sided finite difference.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !
  implicit none

  PetscReal,      intent(in)    :: T        ! temperature [C]
  PetscReal,      intent(in)    :: P        ! pressure [Pa]
  PetscReal,      intent(in)    :: C        ! concentration [mol/L]
  PetscReal,      intent(in)    :: den_kg   ! rho_l [kg/m^3]
  PetscReal,      intent(out)   :: vis      ! mu_l [Pa.s]
  PetscReal,      intent(out)   :: dvis_dT  ! d(mu_l)/dT [Pa.s/C]
  PetscReal,      intent(out)   :: dvis_dC  ! d(mu_l)/dC [Pa.s/(mol/L)]
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: aux1(1)
  PetscReal :: PS, dPS_dT
  PetscReal :: VW, dVW_dT, dVW_dP
  PetscReal :: VW_pert, dVW_dT_pert, dVW_dP_pert
  PetscReal :: s, ds_dC, dvis_ds

  ! Saturation pressure (interface input; unused by Batzle & Wang Eq 32).
  call EOSWaterSaturationPressure(T, PS, dPS_dT, ierr)

  if (.not.EOSWaterDensityIsBrine()) then
    ! Pure-water correlation
    call EOSWaterViscosity(T, P, PS, dPS_dT, VW, dVW_dT, dVW_dP, ierr)
    vis     = VW
    dvis_dT = dVW_dT
    dvis_dC = 0.d0
    return
  endif

  ! C [mol/L] -> salinity mass fraction s [kg/kg] using the current density
  call THCConcToMassFraction(C, den_kg, thc_molar_mass_solute, s, ds_dC)
  aux1(1) = s

  ! Viscosity with analytic T derivative (dVW_dP is hard-set to 0 by the EOS).
  call EOSWaterViscosityExt(T, P, PS, dPS_dT, aux1, VW, dVW_dT, dVW_dP, &
                            ierr)
  vis     = VW
  dvis_dT = dVW_dT

  ! d(mu)/dC by chain rule, d(mu)/ds via one-sided finite difference.
  aux1(1) = s + thc_salinity_perturbation
  call EOSWaterViscosityExt(T, P, PS, dPS_dT, aux1, VW_pert, dVW_dT_pert, &
                            dVW_dP_pert, ierr)
  dvis_ds = (VW_pert - VW) / thc_salinity_perturbation
  dvis_dC = dvis_ds * ds_dC

end subroutine THCViscosityAndDerivs

! ****************************************************************************
subroutine THCThermalConductivityEff(S_l, phi, kappa_s, &
                                         kappa_eff, dkappa_dsat)
  !
  ! Effective porous-medium thermal conductivity via the Somerton et al. (1974)
  ! mixing model. No equivalent exists in the EOS
  ! so it is implemented here entirely from local inputs.
  !
  !   kappa_sat = kappa_s^(1-phi) . kappa_l^phi      [geometric mean, saturated]
  !   kappa_dry = kappa_s^(1-phi) . kappa_air^phi    [geometric mean, dry]
  !   kappa_eff = kappa_dry + sqrt(S_l) . (kappa_sat - kappa_dry)
  !   dkappa/dS_l = (kappa_sat - kappa_dry) / (2 . sqrt(S_l))
  !
  ! kappa_l = 0.6 W/(m.K) (constant), kappa_air = 0.025 W/(m.K).
  ! The saturation is floored at thc_sat_floor so the derivative stays
  ! finite as S_l -> 0 (physically S_l >= S_lr > 0).
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !
  implicit none

  PetscReal, intent(in)  :: S_l          ! liquid saturation [-]
  PetscReal, intent(in)  :: phi          ! porosity [-]
  PetscReal, intent(in)  :: kappa_s      ! solid thermal conductivity [W/(m.K)]
  PetscReal, intent(out) :: kappa_eff    ! effective conductivity [W/(m.K)]
  PetscReal, intent(out) :: dkappa_dsat  ! d(kappa_eff)/d(S_l) [W/(m.K)]

  PetscReal :: kappa_sat, kappa_dry, sqrt_sat, sat

  sat = max(S_l, thc_sat_floor)

  kappa_sat = kappa_s**(1.d0 - phi) * thc_kappa_liquid**phi
  kappa_dry = kappa_s**(1.d0 - phi) * thc_kappa_air**phi

  sqrt_sat    = sqrt(sat)
  kappa_eff   = kappa_dry + sqrt_sat * (kappa_sat - kappa_dry)
  dkappa_dsat = (kappa_sat - kappa_dry) / (2.d0 * sqrt_sat)

end subroutine THCThermalConductivityEff

end module THC_EOS_Utils_module
