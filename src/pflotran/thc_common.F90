module THC_Common_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Global_Aux_module
  use THC_Aux_module

  implicit none

  private

  ! Cutoff parameters
  PetscReal, parameter :: eps     = 1.d-8
  PetscReal, parameter :: floweps = 0.d0

  ! The branch below is implemented now so the energy BC Jacobian is complete.
  PetscInt, parameter :: THC_CONVECTIVE_BC = -2

  ! --- public interface ---------------------------------------------------
  public :: THCAccumulation, &
            THCFlux, &
            THCBCFlux, &
            THCSrcSink, &
            THCAccumDerivative, &
            THCFluxDerivative, &
            THCBCFluxDerivative, &
            THCSrcSinkDerivative

contains

! ************************************************************************** !

subroutine THCAccumulation(thc_auxvar,global_auxvar,material_auxvar, &
                               option,Res,Jac,dResdparam,calculate_derivatives)
  !
  ! Computes the non-fixed portion of the accumulation term for the THC
  ! residual (flow, energy, solute) and its 3x3 accumulation Jacobian block.
  !
  ! Argument list mirrors ZFlowAccumulation. dResdparam is retained for
  ! interface compatibility but is a NO-OP here (zeroed, never filled) --
  ! THC does not yet support the adjoint-parameter sensitivities.
  !
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(thc_auxvar_type) :: thc_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(THC_NDOF)
  PetscReal :: Jac(THC_NDOF,THC_NDOF)
  PetscReal :: dResdparam(THC_NDOF,THC_NDOF)
  PetscBool :: calculate_derivatives

  PetscReal :: porosity
  PetscReal :: saturation
  PetscReal :: por_sat
  PetscReal :: dpor_sat_dp
  PetscReal :: volume_over_dt
  PetscReal :: temperature
  PetscReal :: den_kg
  PetscReal :: den_kmol
  PetscReal :: c_p
  PetscReal :: rho_s_c_s
  PetscReal :: tempreal
  PetscReal, parameter :: L_per_m3 = 1.d3

  Res = 0.d0
  Jac = 0.d0
  dResdparam = 0.d0

  ! v_over_t[m^3 bulk/sec] = vol[m^3 bulk] / dt[sec]
  volume_over_dt = material_auxvar%volume / option%flow_dt
  ! must use thc_auxvar%effective_porosity here as it enables numerical
  ! derivatives to be employed
  porosity    = thc_auxvar%effective_porosity
  saturation  = thc_auxvar%sat
  por_sat     = saturation * porosity
  ! d(por_sat)/dP = dS/dP * por + S * dpor/dP
  dpor_sat_dp = thc_auxvar%dsat_dp * porosity + &
                saturation * thc_auxvar%dpor_dp
  temperature = thc_auxvar%temp
  den_kg      = thc_auxvar%den_kg
  den_kmol    = thc_auxvar%den_kmol
  c_p         = thc_specific_heat_liquid       ! c_p,l [J/(kg.K)]
  rho_s_c_s   = thc_auxvar%heat_cap_solid     ! rho_s * c_s [J/(m^3.K)]

  ! --- Flow (liquid mass) equation : units kmol/s -------------------------
  ! Res[kmol/sec] = sat[m^3 liq/m^3 void] * por[m^3 void/m^3 bulk] *
  !                 rho_kmol[kmol/m^3] * vol[m^3 bulk] / dt[sec]
  tempreal = volume_over_dt * den_kmol
  Res(thc_pressure_dof) = por_sat * tempreal

  ! --- Solute equation : units mol/s --------------------------------------
  ! Res[mole/sec] = c[mol/L] * 1000[L/m^3] * sat * por * vol / dt
  tempreal = L_per_m3 * volume_over_dt
  Res(thc_concentration_dof) = thc_auxvar%conc * tempreal * por_sat

  ! --- Energy equation : units J/s = W ------------------------------------
  ! RHO_CP_T     : Res = [phi*S_l*rho_kg*c_p*T + (1-phi)*rho_s*c_s*T]*vol/dt
  ! FULL_EOS     : Res = [phi*S_l*rho_kg*u(P,T) +(1-phi)*rho_s*c_s*T]*vol/dt
  !                      \_liquid thermal energy_/ \__solid thermal energy__/
  if (thc_energy_mode == THC_ENERGY_FULL_EOS) then
    Res(thc_temperature_dof) = volume_over_dt * &
      (por_sat * den_kg * thc_auxvar%u + &
       (1.d0 - porosity) * rho_s_c_s * temperature)
  else
    Res(thc_temperature_dof) = volume_over_dt * temperature * &
      (por_sat * den_kg * c_p + (1.d0 - porosity) * rho_s_c_s)
  endif

  if (calculate_derivatives) then

    ! ===================================================================== !
    ! Flow equation accumulation Jacobian
    !   A_flow = phi * S_l * rho_kmol * V/dt
    ! Note: rho_kmol derivatives are obtained from the stored mass-density
    ! derivatives via d(rho_kmol)/dX = d(rho_kg)/dX / FMWH2O.
    ! ===================================================================== !
    ! dA_flow/dP = V/dt * [rho_kmol * d(por_sat)/dP + por_sat * d(rho_kmol)/dP]
    Jac(thc_pressure_dof,thc_pressure_dof) = volume_over_dt * &
      (den_kmol * dpor_sat_dp + &
       por_sat * thc_auxvar%dden_dp / FMWH2O)
    ! dA_flow/dT = V/dt * por_sat * d(rho_kmol)/dT
    Jac(thc_pressure_dof,thc_temperature_dof) = volume_over_dt * &
      por_sat * thc_auxvar%dden_dT / FMWH2O
    ! dA_flow/dC = V/dt * por_sat * d(rho_kmol)/dC
    Jac(thc_pressure_dof,thc_concentration_dof) = volume_over_dt * &
      por_sat * thc_auxvar%dden_dC / FMWH2O

    ! ===================================================================== !
    ! Solute equation accumulation Jacobian
    !   A_sol = C * 1000 * phi * S_l * V/dt   (density-independent)
    ! ===================================================================== !
    ! dA_sol/dP = C * 1000 * V/dt * d(por_sat)/dP
    Jac(thc_concentration_dof,thc_pressure_dof) = &
      thc_auxvar%conc * L_per_m3 * volume_over_dt * dpor_sat_dp
    ! dA_sol/dT = 0    (no T dependence in solute accumulation)
    Jac(thc_concentration_dof,thc_temperature_dof) = 0.d0
    ! dA_sol/dC = 1000 * phi * S_l * V/dt
    Jac(thc_concentration_dof,thc_concentration_dof) = &
      L_per_m3 * volume_over_dt * por_sat

    ! ===================================================================== !
    ! Energy equation accumulation Jacobian
    !   A_energy = [phi*S_l*rho_kg*c_p*T + (1-phi)*rho_s*c_s*T] * V/dt
    ! ===================================================================== !
    if (thc_energy_mode == THC_ENERGY_FULL_EOS) then
      ! FULL_EOS: A_e = [phi*S_l*rho_kg*u + (1-phi)*rho_s*c_s*T] * V/dt
      ! dA_e/dP = V/dt * [ u*rho_kg*d(por_sat)/dP + u*por_sat*d(rho_kg)/dP
      !                    + por_sat*rho_kg*du/dP - rho_s*c_s*T*dphi/dP ]
      Jac(thc_temperature_dof,thc_pressure_dof) = volume_over_dt * &
        (thc_auxvar%u * den_kg * dpor_sat_dp + &
         thc_auxvar%u * por_sat * thc_auxvar%dden_dp + &
         por_sat * den_kg * thc_auxvar%du_dP - &
         rho_s_c_s * temperature * thc_auxvar%dpor_dp)
      ! dA_e/dT = V/dt * [ u*por_sat*d(rho_kg)/dT + por_sat*rho_kg*du/dT
      !                    + (1-phi)*rho_s*c_s ]
      Jac(thc_temperature_dof,thc_temperature_dof) = volume_over_dt * &
        (thc_auxvar%u * por_sat * thc_auxvar%dden_dT + &
         por_sat * den_kg * thc_auxvar%du_dT + &
         (1.d0 - porosity) * rho_s_c_s)
      ! dA_e/dC = V/dt * [ u*por_sat*d(rho_kg)/dC + por_sat*rho_kg*du/dC ]
      Jac(thc_temperature_dof,thc_concentration_dof) = volume_over_dt *&
        (thc_auxvar%u * por_sat * thc_auxvar%dden_dC + &
         por_sat * den_kg * thc_auxvar%du_dC)
    else
      ! RHO_CP_T formulation
      ! dA_e/dP = V/dt * T * [ rho_kg*c_p*d(por_sat)/dP
      !                        + phi*S_l*c_p*d(rho_kg)/dP
      !                        - rho_s*c_s*dphi/dP ]
      Jac(thc_temperature_dof,thc_pressure_dof) = volume_over_dt * &
        temperature * (den_kg * c_p * dpor_sat_dp + &
                       por_sat * c_p * thc_auxvar%dden_dp - &
                       rho_s_c_s * thc_auxvar%dpor_dp)
      ! dA_e/dT = V/dt * [ phi*S_l*rho_kg*c_p + (1-phi)*rho_s*c_s
      !                    + phi*S_l*c_p*T*d(rho_kg)/dT ]
      !   first two terms = (rho c)_eff (effective volumetric heat capacity);
      !   third term = temperature-dependent-density correction.
      Jac(thc_temperature_dof,thc_temperature_dof) = volume_over_dt * &
        (por_sat * den_kg * c_p + (1.d0 - porosity) * rho_s_c_s + &
         por_sat * c_p * temperature * thc_auxvar%dden_dT)
      ! dA_e/dC = V/dt * phi*S_l*c_p*T*d(rho_kg)/dC
      Jac(thc_temperature_dof,thc_concentration_dof) = volume_over_dt *&
        por_sat * c_p * temperature * thc_auxvar%dden_dC
    endif

  endif

end subroutine THCAccumulation

! ************************************************************************** !

subroutine THCFlux(thc_auxvar_up,global_auxvar_up, &
                       material_auxvar_up, &
                       thc_auxvar_dn,global_auxvar_dn, &
                       material_auxvar_dn, &
                       area, dist, &
                       thc_parameter, &
                       option,v_darcy,Res, &
                       Jup,Jdn, &
                       dResdparamup,dResdparamdn, &
                       calculate_derivatives, &
                       debug_connection)
  !
  ! Computes the internal flux terms (flow / energy / solute) for the THC
  ! residual and the full 3x3 (P,T,C) Jacobian blocks for a single internal
  ! connection.  Mirrors ZFlowFluxHarmonicPermOnly (zflow_common.F90), extended
  ! to variable density / viscosity and the mandatory 3-DOF system.
  !
  ! Harmonic averaging is used for BOTH the intrinsic permeability AND the
  ! effective thermal conductivity.  A single
  ! upwind direction -- the sign of the total potential difference -- governs
  ! the mobility, the transported density/enthalpy, and the advected
  ! concentration for all three equations.  The sign
  ! convention is identical to ZFLOW
  !
  ! dResdparamup / dResdparamdn are retained for interface compatibility but are
  ! a NO-OP here (zeroed, never filled) -- THC does not yet support the
  ! ZFLOW adjoint-parameter sensitivities.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !
  use Option_module
  use Material_Aux_module
  use Connection_module

  implicit none

  type(thc_auxvar_type) :: thc_auxvar_up, thc_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(thc_parameter_type) :: thc_parameter
  PetscReal :: Res(THC_NDOF)
  PetscReal :: Jup(THC_NDOF,THC_NDOF)
  PetscReal :: Jdn(THC_NDOF,THC_NDOF)
  PetscReal :: dResdparamup(THC_NDOF,THC_NDOF)
  PetscReal :: dResdparamdn(THC_NDOF,THC_NDOF)
  PetscBool :: calculate_derivatives
  PetscBool :: debug_connection

  PetscReal, parameter :: L_per_m3 = 1.d3

  PetscReal :: dist_gravity      ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight
  PetscReal :: geometric_weight_up, geometric_weight_dn

  ! permeability / mobility / volumetric flux
  PetscReal :: perm_up, perm_dn
  PetscReal :: perm_ave_over_dist
  PetscReal :: numerator, denominator
  PetscReal :: denom_D, denom_kappa
  PetscReal :: rho_avg, gravity_term, delta_pressure
  PetscReal :: kr, dkr_dpup, dkr_dpdn
  PetscReal :: vis_upw, mobility, kr_over_vis2
  PetscReal :: Gamma             ! geometric transmissibility = k_harm/d * area
  PetscReal :: q
  PetscReal :: wup, wdn          ! upwind weights (1/0) for transported scalars
  PetscReal :: advective_weight_up, advective_weight_dn

  ! potential derivatives
  PetscReal :: dPi_dpup, dPi_dpdn, dPi_dTup, dPi_dTdn, dPi_dCup, dPi_dCdn
  ! mobility derivatives
  PetscReal :: dmob_dpup, dmob_dpdn, dmob_dTup, dmob_dTdn, dmob_dCup, dmob_dCdn
  ! volumetric flux derivatives
  PetscReal :: dq_dpup, dq_dpdn, dq_dTup, dq_dTdn, dq_dCup, dq_dCdn

  ! upwind-selected transported quantities
  PetscReal :: den_kmol_adv, den_kg_adv, temp_upw, conc_upw, h_upw
  PetscReal :: h_spec_upw        ! upwinded specific enthalpy [J/kg] (FULL_EOS)
  ! molar-density derivative routing (kmol = mass / FMWH2O)
  PetscReal :: ddkmol_dp_up, ddkmol_dp_dn
  PetscReal :: ddkmol_dT_up, ddkmol_dT_dn
  PetscReal :: ddkmol_dC_up, ddkmol_dC_dn
  ! mass-density derivative routing (for enthalpy)
  PetscReal :: ddkg_dp_up, ddkg_dp_dn
  PetscReal :: ddkg_dT_up, ddkg_dT_dn
  PetscReal :: ddkg_dC_up, ddkg_dC_dn
  ! enthalpy derivatives
  PetscReal :: dh_dp_up, dh_dp_dn, dh_dT_up, dh_dT_dn, dh_dC_up, dh_dC_dn

  ! solute diffusion (hydrodynamic dispersion)
  PetscReal :: c_p
  PetscReal :: D_hyd_up, D_hyd_dn
  PetscReal :: Deff_over_dist, Gamma_D, delta_conc
  PetscReal :: dD_hyd_up_dp, dD_hyd_up_dT, dD_hyd_up_dC
  PetscReal :: dD_hyd_dn_dp, dD_hyd_dn_dT, dD_hyd_dn_dC
  PetscReal :: facD_up, facD_dn
  PetscReal :: dGamma_D_dpup, dGamma_D_dpdn
  PetscReal :: dGamma_D_dTup, dGamma_D_dTdn
  PetscReal :: dGamma_D_dCup, dGamma_D_dCdn
  PetscBool :: diffusion_on

  ! energy conduction
  PetscReal :: kappa_up, kappa_dn
  PetscReal :: kappa_over_dist, Gamma_kappa, delta_temp
  PetscReal :: dkappa_up_dp, dkappa_dn_dp
  PetscReal :: dkappa_up_dsat, dkappa_dn_dsat
  PetscReal :: dkappa_up_dT, dkappa_dn_dT
  PetscReal :: dGamma_kappa_dpup, dGamma_kappa_dpdn
  PetscReal :: dGamma_kappa_dTup, dGamma_kappa_dTdn
  PetscBool :: conduction_on

  PetscReal :: tempreal

  Res = 0.d0
  Jup = 0.d0
  Jdn = 0.d0
  dResdparamup = 0.d0
  dResdparamdn = 0.d0
  v_darcy = 0.d0

  q = 0.d0
  dq_dpup = 0.d0; dq_dpdn = 0.d0
  dq_dTup = 0.d0; dq_dTdn = 0.d0
  dq_dCup = 0.d0; dq_dCdn = 0.d0

  c_p = thc_specific_heat_liquid       ! c_p,l [J/(kg.K)]

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)
  geometric_weight_up = upweight
  geometric_weight_dn = 1.d0 - upweight

  ! --------------------------------------------------------------------------
  ! Harmonic-mean intrinsic permeability
  !   k_harm / d_total = (k_up*k_dn) / (d_up*k_dn + d_dn*k_up)
  !   Gamma = (k_harm/d_total) * area
  ! --------------------------------------------------------------------------
  call PermeabilityTensorToScalar(material_auxvar_up,dist,perm_up)
  call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

  numerator = perm_up * perm_dn
  denominator = dist_up*perm_dn + dist_dn*perm_up
  perm_ave_over_dist = numerator / denominator
  Gamma = perm_ave_over_dist * area

  rho_avg = geometric_weight_up * thc_auxvar_up%den_kg + &
            geometric_weight_dn * thc_auxvar_dn%den_kg
  gravity_term = rho_avg * dist_gravity
  delta_pressure = thc_auxvar_up%pres - &
                   thc_auxvar_dn%pres + &
                   gravity_term

  ! single upwind direction for all three equations:
  ! sign(q) = sign(delta_pressure) since Gamma, mobility >= 0.
  if (delta_pressure >= 0.d0) then
    wup = 1.d0
    wdn = 0.d0
  else
    wup = 0.d0
    wdn = 1.d0
  endif

  ! --------------------------------------------------------------------------
  ! Upwind mobility lambda = kr_up / vis_up  (variable viscosity)
  ! --------------------------------------------------------------------------
  kr = 0.d0
  dkr_dpup = 0.d0
  dkr_dpdn = 0.d0
  if (thc_auxvar_up%kr + thc_auxvar_dn%kr > floweps) then
    if (thc_tensorial_rel_perm) then
      if (delta_pressure >= 0.d0) then
        call THCAuxTensorialRelPerm(thc_auxvar_up, &
                thc_parameter% &
                  tensorial_rel_perm_exponent(:,material_auxvar_up%id), &
                dist,kr,dkr_dpup,option)
      else
        call THCAuxTensorialRelPerm(thc_auxvar_dn, &
                thc_parameter% &
                  tensorial_rel_perm_exponent(:,material_auxvar_dn%id), &
                dist,kr,dkr_dpdn,option)
      endif
    else
      if (delta_pressure >= 0.d0) then
        kr = thc_auxvar_up%kr
        dkr_dpup = thc_auxvar_up%dkr_dp
      else
        kr = thc_auxvar_dn%kr
        dkr_dpdn = thc_auxvar_dn%dkr_dp
      endif
    endif
  endif

  ! upwind viscosity (mobility uses the SAME upwind cell as kr)
  vis_upw = wup*thc_auxvar_up%vis + wdn*thc_auxvar_dn%vis
  mobility = kr / vis_upw

  ! q[m^3/sec] = Gamma * mobility * delta_pressure
  q = Gamma * mobility * delta_pressure
  v_darcy = perm_ave_over_dist * mobility * delta_pressure

  ! Density may be fully upwinded or geometrically interpolated to
  ! reproduce TH. Mobility, enthalpy, temperature, and concentration remain
  ! flow-upwinded in both formulations.
  advective_weight_up = wup
  advective_weight_dn = wdn
  if (thc_advective_density_mode == &
      THC_ADVECTIVE_DENSITY_TH_COMPATIBLE) then
    advective_weight_up = geometric_weight_up
    advective_weight_dn = geometric_weight_dn
  endif
  den_kmol_adv = advective_weight_up*thc_auxvar_up%den_kmol + &
                 advective_weight_dn*thc_auxvar_dn%den_kmol
  den_kg_adv   = advective_weight_up*thc_auxvar_up%den_kg + &
                 advective_weight_dn*thc_auxvar_dn%den_kg
  temp_upw     = wup*thc_auxvar_up%temp     + wdn*thc_auxvar_dn%temp
  conc_upw     = wup*thc_auxvar_up%conc     + wdn*thc_auxvar_dn%conc
  ! Volumetric enthalpy [J/m^3] carried by the advective flux:
  !   RHO_CP_T: rho_kg * c_p * T   ;   FULL_EOS: rho_kg * h_spec(P,T)
  if (thc_energy_mode == THC_ENERGY_FULL_EOS) then
    h_spec_upw = wup*thc_auxvar_up%h + wdn*thc_auxvar_dn%h
    h_upw      = den_kg_adv * h_spec_upw
  else
    h_spec_upw = c_p * temp_upw
    h_upw      = den_kg_adv * h_spec_upw
  endif

  ! --------------------------------------------------------------------------
  ! Solute hydrodynamic dispersion with harmonic-mean D_eff
  !   D_hyd = phi * S_l * tau * D_mol      (D_mol = D_mol(T,C), per-cell)
  ! --------------------------------------------------------------------------
  D_hyd_up = thc_auxvar_up%effective_porosity * thc_auxvar_up%sat * &
             material_auxvar_up%tortuosity * thc_auxvar_up%diff_mol
  D_hyd_dn = thc_auxvar_dn%effective_porosity * thc_auxvar_dn%sat * &
             material_auxvar_dn%tortuosity * thc_auxvar_dn%diff_mol
  numerator = D_hyd_up * D_hyd_dn
  denom_D = dist_up*D_hyd_dn + dist_dn*D_hyd_up
  diffusion_on = PETSC_TRUE
  if (denom_D == 0.d0) then
    ! turn off diffusion
    numerator = 0.d0
    denom_D = 1.d0
    diffusion_on = PETSC_FALSE
  endif
  Deff_over_dist = numerator / denom_D
  Gamma_D = area * Deff_over_dist
  delta_conc = thc_auxvar_up%conc - thc_auxvar_dn%conc

  ! --------------------------------------------------------------------------
  ! Energy conduction with harmonic-mean effective thermal conductivity
  ! kappa_face/d = (k_up*k_dn)/(d_up*k_dn+d_dn*k_up)
  ! --------------------------------------------------------------------------
  ! evaluated per connection so an anisotropic tensor can be projected onto
  ! this face; the auxvar copy is an unprojected per-cell value for output only
  call THCThermalCondEval(thc_auxvar_up,material_auxvar_up%id,thc_parameter, &
                          dist,PETSC_TRUE,kappa_up,dkappa_up_dsat, &
                          dkappa_up_dT,option)
  call THCThermalCondEval(thc_auxvar_dn,material_auxvar_dn%id,thc_parameter, &
                          dist,PETSC_TRUE,kappa_dn,dkappa_dn_dsat, &
                          dkappa_dn_dT,option)
  numerator = kappa_up * kappa_dn
  denom_kappa = dist_up*kappa_dn + dist_dn*kappa_up
  conduction_on = PETSC_TRUE
  if (denom_kappa == 0.d0) then
    ! turn off conduction
    numerator = 0.d0
    denom_kappa = 1.d0
    conduction_on = PETSC_FALSE
  endif
  kappa_over_dist = numerator / denom_kappa
  Gamma_kappa = area * kappa_over_dist
  delta_temp = thc_auxvar_up%temp - thc_auxvar_dn%temp

  ! --------------------------------------------------------------------------
  ! Residuals (same sign convention as ZFLOW)
  ! --------------------------------------------------------------------------
  ! Flow [kmol/s]
  Res(thc_pressure_dof) = den_kmol_adv * q
  ! Solute [mol/s] : advection + hydrodynamic dispersion
  Res(thc_concentration_dof) = (conc_upw * q + Gamma_D * delta_conc) * &
                                   L_per_m3
  ! Energy [W] : advected enthalpy + Fourier conduction
  Res(thc_temperature_dof) = h_upw * q + Gamma_kappa * delta_temp

  if (calculate_derivatives) then

    ! ===================================================================== !
    ! Building blocks
    ! ===================================================================== !
    ! Potential derivatives  dPi = d(P_up - P_dn + rho_avg*g*dz)
    dPi_dpup =  1.d0 + geometric_weight_up * thc_auxvar_up%dden_dp * &
                        dist_gravity
    dPi_dpdn = -1.d0 + geometric_weight_dn * thc_auxvar_dn%dden_dp * &
                        dist_gravity
    dPi_dTup = geometric_weight_up * thc_auxvar_up%dden_dT * dist_gravity
    dPi_dTdn = geometric_weight_dn * thc_auxvar_dn%dden_dT * dist_gravity
    dPi_dCup = geometric_weight_up * thc_auxvar_up%dden_dC * dist_gravity
    dPi_dCdn = geometric_weight_dn * thc_auxvar_dn%dden_dC * dist_gravity

    ! Mobility derivatives  lambda = kr_upw / vis_upw  (upwind only)
    kr_over_vis2 = kr / (vis_upw*vis_upw)
    dmob_dpup = dkr_dpup / vis_upw
    dmob_dpdn = dkr_dpdn / vis_upw
    dmob_dTup = -kr_over_vis2 * (wup*thc_auxvar_up%dvis_dT)
    dmob_dTdn = -kr_over_vis2 * (wdn*thc_auxvar_dn%dvis_dT)
    dmob_dCup = -kr_over_vis2 * (wup*thc_auxvar_up%dvis_dC)
    dmob_dCdn = -kr_over_vis2 * (wdn*thc_auxvar_dn%dvis_dC)

    ! Volumetric flux derivatives  q = Gamma * lambda * dPi (product rule)
    dq_dpup = Gamma * (dmob_dpup*delta_pressure + mobility*dPi_dpup)
    dq_dpdn = Gamma * (dmob_dpdn*delta_pressure + mobility*dPi_dpdn)
    dq_dTup = Gamma * (dmob_dTup*delta_pressure + mobility*dPi_dTup)
    dq_dTdn = Gamma * (dmob_dTdn*delta_pressure + mobility*dPi_dTdn)
    dq_dCup = Gamma * (dmob_dCup*delta_pressure + mobility*dPi_dCup)
    dq_dCdn = Gamma * (dmob_dCdn*delta_pressure + mobility*dPi_dCdn)

    ! Advective-density derivative routing (kmol = mass/FMWH2O)
    ddkmol_dp_up = advective_weight_up*thc_auxvar_up%dden_dp / FMWH2O
    ddkmol_dp_dn = advective_weight_dn*thc_auxvar_dn%dden_dp / FMWH2O
    ddkmol_dT_up = advective_weight_up*thc_auxvar_up%dden_dT / FMWH2O
    ddkmol_dT_dn = advective_weight_dn*thc_auxvar_dn%dden_dT / FMWH2O
    ddkmol_dC_up = advective_weight_up*thc_auxvar_up%dden_dC / FMWH2O
    ddkmol_dC_dn = advective_weight_dn*thc_auxvar_dn%dden_dC / FMWH2O
    ddkg_dp_up   = advective_weight_up*thc_auxvar_up%dden_dp
    ddkg_dp_dn   = advective_weight_dn*thc_auxvar_dn%dden_dp
    ddkg_dT_up   = advective_weight_up*thc_auxvar_up%dden_dT
    ddkg_dT_dn   = advective_weight_dn*thc_auxvar_dn%dden_dT
    ddkg_dC_up   = advective_weight_up*thc_auxvar_up%dden_dC
    ddkg_dC_dn   = advective_weight_dn*thc_auxvar_dn%dden_dC

    ! Diffusive transmissibility derivatives (harmonic D_eff)
    dGamma_D_dpup = 0.d0; dGamma_D_dpdn = 0.d0
    dGamma_D_dTup = 0.d0; dGamma_D_dTdn = 0.d0
    dGamma_D_dCup = 0.d0; dGamma_D_dCdn = 0.d0
    if (diffusion_on) then
      ! cell-level D_eff derivatives
      dD_hyd_up_dp = (thc_auxvar_up%dsat_dp * &
                      thc_auxvar_up%effective_porosity + &
                      thc_auxvar_up%sat * thc_auxvar_up%dpor_dp) * &
                     material_auxvar_up%tortuosity * thc_auxvar_up%diff_mol
      dD_hyd_up_dT = thc_auxvar_up%effective_porosity * &
                     thc_auxvar_up%sat * material_auxvar_up%tortuosity * &
                     thc_auxvar_up%ddiff_dT
      dD_hyd_up_dC = thc_auxvar_up%effective_porosity * &
                     thc_auxvar_up%sat * material_auxvar_up%tortuosity * &
                     thc_auxvar_up%ddiff_dC
      dD_hyd_dn_dp = (thc_auxvar_dn%dsat_dp * &
                      thc_auxvar_dn%effective_porosity + &
                      thc_auxvar_dn%sat * thc_auxvar_dn%dpor_dp) * &
                     material_auxvar_dn%tortuosity * thc_auxvar_dn%diff_mol
      dD_hyd_dn_dT = thc_auxvar_dn%effective_porosity * &
                     thc_auxvar_dn%sat * material_auxvar_dn%tortuosity * &
                     thc_auxvar_dn%ddiff_dT
      dD_hyd_dn_dC = thc_auxvar_dn%effective_porosity * &
                     thc_auxvar_dn%sat * material_auxvar_dn%tortuosity * &
                     thc_auxvar_dn%ddiff_dC
      ! d(Deff_over_dist)/dD_hyd_up = D_hyd_dn^2 * dist_up / denom^2  (ZFLOW form)
      tempreal = denom_D * denom_D
      facD_up = area * D_hyd_dn * D_hyd_dn * dist_up / tempreal
      facD_dn = area * D_hyd_up * D_hyd_up * dist_dn / tempreal
      dGamma_D_dpup = facD_up * dD_hyd_up_dp
      dGamma_D_dpdn = facD_dn * dD_hyd_dn_dp
      dGamma_D_dTup = facD_up * dD_hyd_up_dT
      dGamma_D_dTdn = facD_dn * dD_hyd_dn_dT
      dGamma_D_dCup = facD_up * dD_hyd_up_dC
      dGamma_D_dCdn = facD_dn * dD_hyd_dn_dC
    endif

    ! Thermal transmissibility derivatives (harmonic kappa; kappa = kappa(S_l))
    dGamma_kappa_dpup = 0.d0
    dGamma_kappa_dpdn = 0.d0
    dGamma_kappa_dTup = 0.d0
    dGamma_kappa_dTdn = 0.d0
    if (conduction_on) then
      dkappa_up_dp = dkappa_up_dsat * thc_auxvar_up%dsat_dp
      dkappa_dn_dp = dkappa_dn_dsat * thc_auxvar_dn%dsat_dp
      tempreal = denom_kappa * denom_kappa
      dGamma_kappa_dpup = area * kappa_dn * kappa_dn * dist_up / tempreal * &
                          dkappa_up_dp
      dGamma_kappa_dpdn = area * kappa_up * kappa_up * dist_dn / tempreal * &
                          dkappa_dn_dp
      ! kappa(T) only via THERMAL_CHARACTERISTIC_CURVES; zero otherwise
      dGamma_kappa_dTup = area * kappa_dn * kappa_dn * dist_up / tempreal * &
                          dkappa_up_dT
      dGamma_kappa_dTdn = area * kappa_up * kappa_up * dist_dn / tempreal * &
                          dkappa_dn_dT
    endif

    ! Enthalpy derivatives  (volumetric enthalpy carried by flux)
    !   RHO_CP_T: h_upw = rho_kg_adv * c_p * T_upw
    !   FULL_EOS: h_upw = rho_kg_adv * h_spec_upw
    if (thc_energy_mode == THC_ENERGY_FULL_EOS) then
      dh_dp_up = ddkg_dp_up*h_spec_upw + den_kg_adv*wup*thc_auxvar_up%dh_dP
      dh_dp_dn = ddkg_dp_dn*h_spec_upw + den_kg_adv*wdn*thc_auxvar_dn%dh_dP
      dh_dT_up = ddkg_dT_up*h_spec_upw + den_kg_adv*wup*thc_auxvar_up%dh_dT
      dh_dT_dn = ddkg_dT_dn*h_spec_upw + den_kg_adv*wdn*thc_auxvar_dn%dh_dT
      dh_dC_up = ddkg_dC_up*h_spec_upw + den_kg_adv*wup*thc_auxvar_up%dh_dC
      dh_dC_dn = ddkg_dC_dn*h_spec_upw + den_kg_adv*wdn*thc_auxvar_dn%dh_dC
    else
      dh_dp_up = c_p * temp_upw * ddkg_dp_up
      dh_dp_dn = c_p * temp_upw * ddkg_dp_dn
      dh_dT_up = c_p * (ddkg_dT_up*temp_upw + den_kg_adv*wup)
      dh_dT_dn = c_p * (ddkg_dT_dn*temp_upw + den_kg_adv*wdn)
      dh_dC_up = c_p * temp_upw * ddkg_dC_up
      dh_dC_dn = c_p * temp_upw * ddkg_dC_dn
    endif

    ! ===================================================================== !
    ! Flow equation flux Jacobian  F_flow = rho_kmol_adv * q
    ! ===================================================================== !
    Jup(thc_pressure_dof,thc_pressure_dof) = &
      ddkmol_dp_up*q + den_kmol_adv*dq_dpup
    Jup(thc_pressure_dof,thc_temperature_dof) = &
      ddkmol_dT_up*q + den_kmol_adv*dq_dTup
    Jup(thc_pressure_dof,thc_concentration_dof) = &
      ddkmol_dC_up*q + den_kmol_adv*dq_dCup
    Jdn(thc_pressure_dof,thc_pressure_dof) = &
      ddkmol_dp_dn*q + den_kmol_adv*dq_dpdn
    Jdn(thc_pressure_dof,thc_temperature_dof) = &
      ddkmol_dT_dn*q + den_kmol_adv*dq_dTdn
    Jdn(thc_pressure_dof,thc_concentration_dof) = &
      ddkmol_dC_dn*q + den_kmol_adv*dq_dCdn

    ! ===================================================================== !
    ! Solute equation flux Jacobian
    !   F_sol = (C_upw*q + Gamma_D*dC) * 1000
    !   d(C_upw)/dC_up = wup, d(C_upw)/dC_dn = wdn ; d(dC)/dC_up = +1, dn = -1
    ! ===================================================================== !
    Jup(thc_concentration_dof,thc_pressure_dof) = &
      (conc_upw*dq_dpup + dGamma_D_dpup*delta_conc) * L_per_m3
    Jup(thc_concentration_dof,thc_temperature_dof) = &
      (conc_upw*dq_dTup + dGamma_D_dTup*delta_conc) * L_per_m3
    Jup(thc_concentration_dof,thc_concentration_dof) = &
      (wup*q + conc_upw*dq_dCup + dGamma_D_dCup*delta_conc + Gamma_D) * L_per_m3
    Jdn(thc_concentration_dof,thc_pressure_dof) = &
      (conc_upw*dq_dpdn + dGamma_D_dpdn*delta_conc) * L_per_m3
    Jdn(thc_concentration_dof,thc_temperature_dof) = &
      (conc_upw*dq_dTdn + dGamma_D_dTdn*delta_conc) * L_per_m3
    Jdn(thc_concentration_dof,thc_concentration_dof) = &
      (wdn*q + conc_upw*dq_dCdn + dGamma_D_dCdn*delta_conc - Gamma_D) * L_per_m3

    ! ===================================================================== !
    ! Energy equation flux Jacobian
    !   F_energy = h_upw*q + Gamma_kappa*dT
    !   d(dT)/dT_up = +1, d(dT)/dT_dn = -1 ; Gamma_kappa indep of C
    ! ===================================================================== !
    Jup(thc_temperature_dof,thc_pressure_dof) = &
      dh_dp_up*q + h_upw*dq_dpup + dGamma_kappa_dpup*delta_temp
    Jup(thc_temperature_dof,thc_temperature_dof) = &
      dh_dT_up*q + h_upw*dq_dTup + Gamma_kappa + &
      dGamma_kappa_dTup*delta_temp
    Jup(thc_temperature_dof,thc_concentration_dof) = &
      dh_dC_up*q + h_upw*dq_dCup
    Jdn(thc_temperature_dof,thc_pressure_dof) = &
      dh_dp_dn*q + h_upw*dq_dpdn + dGamma_kappa_dpdn*delta_temp
    Jdn(thc_temperature_dof,thc_temperature_dof) = &
      dh_dT_dn*q + h_upw*dq_dTdn - Gamma_kappa + &
      dGamma_kappa_dTdn*delta_temp
    Jdn(thc_temperature_dof,thc_concentration_dof) = &
      dh_dC_dn*q + h_upw*dq_dCdn

  endif

  if (debug_connection) then
    write(*,'(9x,"IFF(res_f,kr,dp): ",8es12.4)') Res(thc_pressure_dof), &
                                                 kr, delta_pressure
    write(*,'(9x,"IFE(res_e,Gk,dT): ",8es12.4)') Res(thc_temperature_dof), &
                                                 Gamma_kappa, delta_temp
    write(*,'(9x,"IFS(res_s,q,dC): ",8es12.4)') Res(thc_concentration_dof),&
                                                q, delta_conc
  endif

end subroutine THCFlux

! ************************************************************************** !

subroutine THCBCFlux(ibndtype,auxvar_mapping,auxvars, &
                         thc_auxvar_up,global_auxvar_up, &
                         thc_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         area,dist, &
                         thc_parameter, &
                         option,v_darcy,Res,Jdn, &
                         dResdparamdn, &
                         calculate_derivatives, &
                         debug_connection)
  !
  ! Computes the boundary flux terms (flow / energy / solute) for the THC
  ! residual and the 3x3 boundary Jacobian block (interior cell only).
  ! Mirrors ZFlowBCFluxHarmonicPermOnly (zflow_common.F90), extended to the
  ! mandatory 3-DOF system with variable density / viscosity.
  !
  ! The "up" cell is the boundary ghost (its state is fixed by the BC), the
  ! "dn" cell is the interior cell.  Only Jdn (the interior cell) is filled.
  !
  ! Flow / solute BCs:
  !   DIRICHLET / SEEPAGE / CONDUCTANCE / HYDROSTATIC(+seepage/conductance) /
  !   PONDED_WATER  -- pressure-driven flux with variable-density hydrostatic
  !                    gravity. UPWIND uses arithmetic boundary density;
  !                    TH_COMPATIBLE uses the prescribed boundary density ;
  !   NEUMANN       -- prescribed Darcy velocity.
  !
  ! Thermal BCs, keyed on the energy condition slot:
  !   DIRICHLET             -- prescribed T at ghost: advection + conduction ;
  !   NEUMANN               -- prescribed heat flux q_T [W/m^2] (no Jacobian) ;
  !   THC_CONVECTIVE_BC -- Robin / third-kind: h_conv*(T_ext - T_cell) .
  !
  ! Injection vs extraction for T and C is by the sign of the volumetric flux q
  ! inflow carries the boundary (ghost) state,
  ! outflow carries the interior cell state -- a single upwind direction shared
  ! by all three equations (consistent with THCFlux and ZFLOW).
  !
  ! dResdparamdn is retained for interface compatibility but is a NO-OP here.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !
  use Option_module
  use Material_Aux_module
  use String_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(:)
  PetscInt :: auxvar_mapping(THC_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(thc_auxvar_type) :: thc_auxvar_up, thc_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(thc_parameter_type) :: thc_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(THC_NDOF)
  PetscReal :: Jdn(THC_NDOF,THC_NDOF)
  PetscReal :: dResdparamdn(THC_NDOF,THC_NDOF)
  PetscBool :: calculate_derivatives
  PetscBool :: debug_connection

  PetscReal, parameter :: L_per_m3 = 1.d3

  PetscInt :: bc_type, energy_bc_type
  PetscInt :: idof
  PetscReal :: perm_dn, perm_ave_over_dist
  PetscReal :: dist_gravity
  PetscReal :: gravity_weight_up, gravity_weight_dn
  PetscReal :: rho_avg, gravity_term, delta_pressure, ddelta_pressure_dpdn
  PetscReal :: boundary_pressure
  PetscReal :: kr, dkr_dpdn
  PetscReal :: vis_upw, mobility, kr_over_vis2
  PetscReal :: Gamma_bc, q
  PetscReal :: wup, wdn
  PetscReal :: advective_weight_up, advective_weight_dn
  PetscBool :: derivative_toggle

  ! potential / mobility / flux derivatives (interior-cell only)
  PetscReal :: dPi_dpdn, dPi_dTdn, dPi_dCdn
  PetscReal :: dmob_dpdn, dmob_dTdn, dmob_dCdn
  PetscReal :: dq_dpdn, dq_dTdn, dq_dCdn

  ! transported quantities (density follows ADVECTIVE_DENSITY)
  PetscReal :: c_p
  PetscReal :: den_kmol_adv, den_kg_adv, temp_upw, conc_upw, h_upw
  PetscReal :: h_spec_upw        ! upwinded specific enthalpy [J/kg] (FULL_EOS)
  PetscReal :: ddkmol_dp_dn, ddkmol_dT_dn, ddkmol_dC_dn
  PetscReal :: ddkg_dp_dn, ddkg_dT_dn, ddkg_dC_dn
  PetscReal :: dh_dp_dn, dh_dT_dn, dh_dC_dn

  ! solute diffusion
  PetscReal :: delta_conc, D_hyd_dn, Deff_over_dist, Gamma_D, dispersion_scale
  PetscReal :: dD_hyd_dn_dp, dD_hyd_dn_dT, dD_hyd_dn_dC
  PetscReal :: dGamma_D_dpdn, dGamma_D_dTdn, dGamma_D_dCdn

  ! energy conduction / prescribed-flux / convective
  PetscReal :: Gamma_kappa, delta_temp, dGamma_kappa_dpdn, dGamma_kappa_dTdn
  PetscReal :: kappa_bc, dkappa_bc_dsat, dkappa_bc_dT
  PetscReal :: heat_flux_bc, h_conv, temp_ext

  Res = 0.d0
  Jdn = 0.d0
  dResdparamdn = 0.d0
  v_darcy = 0.d0
  q = 0.d0
  delta_pressure = 0.d0
  delta_conc = 0.d0
  ddelta_pressure_dpdn = 0.d0
  perm_ave_over_dist = 0.d0
  kr = 0.d0
  dkr_dpdn = 0.d0
  dq_dpdn = 0.d0; dq_dTdn = 0.d0; dq_dCdn = 0.d0
  wup = 0.d0; wdn = 0.d0
  Gamma_bc = 0.d0
  mobility = 0.d0
  vis_upw = thc_auxvar_dn%vis

  c_p = thc_specific_heat_liquid       ! c_p,l [J/(kg.K)]

  ! dist(0) = scalar distance to boundary face; dist(1:3) = unit vector
  dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

  gravity_weight_up = 1.d0
  gravity_weight_dn = 0.d0

  ! ========================================================================= !
  ! Flow equation -- determine the volumetric flux q
  ! ========================================================================= !
  call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

  bc_type = ibndtype(auxvar_mapping(THC_COND_WATER_INDEX))
  derivative_toggle = PETSC_TRUE
  select case(bc_type)
    case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
         HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
         PONDED_WATER_BC)
      if (thc_auxvar_up%kr + thc_auxvar_dn%kr > floweps) then

        if (bc_type == DIRICHLET_CONDUCTANCE_BC .or. &
            bc_type == HYDROSTATIC_CONDUCTANCE_BC) then
          ! conductance already includes 1/dist; viscosity applied via mobility
          idof = auxvar_mapping(THC_COND_WATER_AUX_INDEX)
          perm_ave_over_dist = auxvars(idof)
        else
          perm_ave_over_dist = perm_dn / dist(0)
        endif

        ! variable-density hydrostatic gravity term
        rho_avg = gravity_weight_up*thc_auxvar_up%den_kg + &
                  gravity_weight_dn*thc_auxvar_dn%den_kg
        gravity_term = rho_avg * dist_gravity
        boundary_pressure = thc_auxvar_up%pres
        delta_pressure = boundary_pressure - thc_auxvar_dn%pres + &
                         gravity_term
        ddelta_pressure_dpdn = -1.d0

        select case(bc_type)
          case(DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
               HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
               PONDED_WATER_BC)
                ! flow in         ! boundary cell is <= pref
            if (delta_pressure > 0.d0 .and. &
                thc_auxvar_up%pres - &
                  option%flow%reference_pressure < eps) then
              delta_pressure = 0.d0
              ddelta_pressure_dpdn = 0.d0
            endif
        end select

        ! upwind mobility direction (shared by all three equations)
        if (delta_pressure >= 0.d0) then
          wup = 1.d0
          wdn = 0.d0
        else
          wup = 0.d0
          wdn = 1.d0
        endif

        if (thc_tensorial_rel_perm) then
          if (delta_pressure >= 0.d0) then
            call THCAuxTensorialRelPerm(thc_auxvar_up, &
                    thc_parameter% &
                      tensorial_rel_perm_exponent(:,material_auxvar_dn%id), &
                    dist,kr,dkr_dpdn,option)
            ! upwind (boundary) pressure is fixed
            dkr_dpdn = 0.d0
          else
            call THCAuxTensorialRelPerm(thc_auxvar_dn, &
                    thc_parameter% &
                      tensorial_rel_perm_exponent(:,material_auxvar_dn%id), &
                    dist,kr,dkr_dpdn,option)
          endif
        else
          if (delta_pressure >= 0.d0) then
            kr = thc_auxvar_up%kr
            dkr_dpdn = 0.d0
          else
            kr = thc_auxvar_dn%kr
            dkr_dpdn = thc_auxvar_dn%dkr_dp
          endif
        endif

        ! upwind viscosity (boundary state fixed; interior varies)
        vis_upw = wup*thc_auxvar_up%vis + wdn*thc_auxvar_dn%vis
        mobility = kr / vis_upw
        ! v_darcy[m/sec] = perm/dist * kr/mu * dP
        v_darcy = perm_ave_over_dist * mobility * delta_pressure

        if (bc_type == PONDED_WATER_BC) then
          idof = auxvar_mapping(THC_COND_WATER_INDEX)
          if (v_darcy > auxvars(idof) / option%flow_dt) then
            v_darcy = auxvars(idof) / option%flow_dt
            ddelta_pressure_dpdn = 0.d0
            dkr_dpdn = 0.d0
            derivative_toggle = PETSC_FALSE
          endif
        endif
      endif

    case(NEUMANN_BC)
      idof = auxvar_mapping(THC_COND_WATER_INDEX)
      if (dabs(auxvars(idof)) > floweps) then
        v_darcy = auxvars(idof)
      endif
      derivative_toggle = PETSC_FALSE

    case default
      option%io_buffer = &
        'Boundary condition type (' // trim(StringWrite(bc_type)) // &
        ') not recognized in THCBCFlux phase loop.'
      call PrintErrMsg(option)
  end select

  ! q[m^3 liquid/sec] = v_darcy[m/sec] * area[m^2]
  q = v_darcy * area
  Gamma_bc = perm_ave_over_dist * area

  ! Upwind direction for transported scalars. For Dirichlet the wup/wdn set
  ! above (sign of delta_pressure) already matches sign(q); for Neumann (and
  ! the ponded cap) reset it from the prescribed velocity.
  if (.not. derivative_toggle .or. bc_type == NEUMANN_BC) then
    if (q >= 0.d0) then
      wup = 1.d0; wdn = 0.d0
    else
      wup = 0.d0; wdn = 1.d0
    endif
  endif

  ! TH pressure-driven boundaries carry boundary-state density independently
  ! of the enthalpy/concentration upwind direction. Prescribed Neumann flow
  ! already uses the same inflow/outflow density selection in TH and THC.
  advective_weight_up = wup
  advective_weight_dn = wdn
  if (thc_advective_density_mode == &
      THC_ADVECTIVE_DENSITY_TH_COMPATIBLE .and. bc_type /= NEUMANN_BC) then
    advective_weight_up = 1.d0
    advective_weight_dn = 0.d0
  endif
  den_kmol_adv = advective_weight_up*thc_auxvar_up%den_kmol + &
                 advective_weight_dn*thc_auxvar_dn%den_kmol
  den_kg_adv   = advective_weight_up*thc_auxvar_up%den_kg + &
                 advective_weight_dn*thc_auxvar_dn%den_kg
  temp_upw     = wup*thc_auxvar_up%temp     + wdn*thc_auxvar_dn%temp
  conc_upw     = wup*thc_auxvar_up%conc     + wdn*thc_auxvar_dn%conc
  ! Volumetric enthalpy [J/m^3]: RHO_CP_T -> rho*c_p*T ; FULL_EOS -> rho*h_spec
  if (thc_energy_mode == THC_ENERGY_FULL_EOS) then
    h_spec_upw = wup*thc_auxvar_up%h + wdn*thc_auxvar_dn%h
    h_upw      = den_kg_adv * h_spec_upw
  else
    h_spec_upw = c_p * temp_upw
    h_upw      = den_kg_adv * h_spec_upw
  endif

  ! --- advective residual contributions (all three equations) -------------
  ! Flow [kmol/s]
  Res(thc_pressure_dof) = den_kmol_adv * q
  ! Solute advection [mol/s] (diffusion added below)
  Res(thc_concentration_dof) = conc_upw * q * L_per_m3
  ! Energy advection [W] (conduction / prescribed flux added below)
  Res(thc_temperature_dof) = h_upw * q

  ! ========================================================================= !
  ! Solute diffusion at the boundary (zero-gradient on outflow)
  ! ========================================================================= !
  delta_conc = thc_auxvar_up%conc - thc_auxvar_dn%conc
  if (q >= 0.d0) then
    dispersion_scale = 1.d0   ! inflow: full diffusive exchange with ghost
  else
    dispersion_scale = 0.d0   ! outflow: zero-gradient outlet BC
  endif
  D_hyd_dn = thc_auxvar_dn%effective_porosity * thc_auxvar_dn%sat * &
             material_auxvar_dn%tortuosity * thc_auxvar_dn%diff_mol
  Deff_over_dist = dispersion_scale * D_hyd_dn / dist(0)
  Gamma_D = area * Deff_over_dist
  Res(thc_concentration_dof) = Res(thc_concentration_dof) + &
                                   Gamma_D * delta_conc * L_per_m3

  ! ========================================================================= !
  ! Energy conduction / prescribed flux / convective BC
  ! ========================================================================= !
  energy_bc_type = ibndtype(auxvar_mapping(THC_COND_ENERGY_INDEX))
  Gamma_kappa = 0.d0
  dGamma_kappa_dpdn = 0.d0
  dGamma_kappa_dTdn = 0.d0
  delta_temp = thc_auxvar_up%temp - thc_auxvar_dn%temp
  select case(energy_bc_type)
    case(DIRICHLET_BC,HYDROSTATIC_BC,DIRICHLET_SEEPAGE_BC, &
         DIRICHLET_CONDUCTANCE_BC,HYDROSTATIC_SEEPAGE_BC, &
         HYDROSTATIC_CONDUCTANCE_BC)
      ! prescribed T at ghost: Fourier conduction over the half-cell distance,
      ! using the interior effective conductivity.
      call THCThermalCondEval(thc_auxvar_dn,material_auxvar_dn%id, &
                              thc_parameter,dist,PETSC_TRUE,kappa_bc, &
                              dkappa_bc_dsat,dkappa_bc_dT,option)
      Gamma_kappa = kappa_bc / dist(0) * area
      Res(thc_temperature_dof) = Res(thc_temperature_dof) + &
                                     Gamma_kappa * delta_temp
      dGamma_kappa_dpdn = dkappa_bc_dsat * &
                          thc_auxvar_dn%dsat_dp / dist(0) * area
      dGamma_kappa_dTdn = dkappa_bc_dT / dist(0) * area
    case(NEUMANN_BC)
      ! prescribed heat flux q_T [W/m^2]; +q_T means heat enters the domain
      ! No Jacobian (value independent of solution).
      heat_flux_bc = auxvars(auxvar_mapping(THC_COND_ENERGY_INDEX))
      Res(thc_temperature_dof) = Res(thc_temperature_dof) + &
                                     heat_flux_bc * area
    case(THC_CONVECTIVE_BC)
      ! Robin / third-kind: heat flux INTO domain = h_conv*(T_ext - T_cell)
      ! h_conv stored in the water-aux slot, T_ext in the
      ! energy slot (placeholder plumbing -- see module header).
      temp_ext = auxvars(auxvar_mapping(THC_COND_ENERGY_INDEX))
      h_conv = auxvars(auxvar_mapping(THC_COND_WATER_AUX_INDEX))
      Res(thc_temperature_dof) = Res(thc_temperature_dof) + &
        h_conv * area * (temp_ext - thc_auxvar_dn%temp)
    case(ZERO_GRADIENT_BC,DIRICHLET_ZERO_GRADIENT_BC,NULL_CONDITION)
      ! no conductive boundary term -- advection only
    case default
      ! treat any other energy BC type as zero-gradient (advection only)
  end select

  ! ========================================================================= !
  ! Boundary Jacobian (interior cell only)
  ! ========================================================================= !
  if (calculate_derivatives) then

    ! --- potential / mobility / flux derivatives (Dirichlet family only) ---
    if (derivative_toggle) then
      dPi_dpdn = ddelta_pressure_dpdn + gravity_weight_dn * &
                 thc_auxvar_dn%dden_dp * dist_gravity
      dPi_dTdn = gravity_weight_dn*thc_auxvar_dn%dden_dT*dist_gravity
      dPi_dCdn = gravity_weight_dn*thc_auxvar_dn%dden_dC*dist_gravity

      kr_over_vis2 = kr / (vis_upw*vis_upw)
      dmob_dpdn = dkr_dpdn / vis_upw
      dmob_dTdn = -kr_over_vis2 * (wdn*thc_auxvar_dn%dvis_dT)
      dmob_dCdn = -kr_over_vis2 * (wdn*thc_auxvar_dn%dvis_dC)

      dq_dpdn = Gamma_bc * (dmob_dpdn*delta_pressure + mobility*dPi_dpdn)
      dq_dTdn = Gamma_bc * (dmob_dTdn*delta_pressure + mobility*dPi_dTdn)
      dq_dCdn = Gamma_bc * (dmob_dCdn*delta_pressure + mobility*dPi_dCdn)
    else
      dq_dpdn = 0.d0; dq_dTdn = 0.d0; dq_dCdn = 0.d0
    endif

    ! Interior-cell advective-density derivatives. These are zero for
    ! TH-compatible pressure-driven boundaries and upwinded for Neumann flow.
    ddkmol_dp_dn = advective_weight_dn*thc_auxvar_dn%dden_dp / FMWH2O
    ddkmol_dT_dn = advective_weight_dn*thc_auxvar_dn%dden_dT / FMWH2O
    ddkmol_dC_dn = advective_weight_dn*thc_auxvar_dn%dden_dC / FMWH2O
    ddkg_dp_dn   = advective_weight_dn*thc_auxvar_dn%dden_dp
    ddkg_dT_dn   = advective_weight_dn*thc_auxvar_dn%dden_dT
    ddkg_dC_dn   = advective_weight_dn*thc_auxvar_dn%dden_dC

    ! enthalpy derivatives (volumetric); interior (dn) cell only at boundary
    !   RHO_CP_T: h_upw = rho_kg_adv * c_p * T_upw
    !   FULL_EOS: h_upw = rho_kg_adv * h_spec_upw
    if (thc_energy_mode == THC_ENERGY_FULL_EOS) then
      dh_dp_dn = ddkg_dp_dn*h_spec_upw + den_kg_adv*wdn*thc_auxvar_dn%dh_dP
      dh_dT_dn = ddkg_dT_dn*h_spec_upw + den_kg_adv*wdn*thc_auxvar_dn%dh_dT
      dh_dC_dn = ddkg_dC_dn*h_spec_upw + den_kg_adv*wdn*thc_auxvar_dn%dh_dC
    else
      dh_dp_dn = c_p * temp_upw * ddkg_dp_dn
      dh_dT_dn = c_p * (ddkg_dT_dn*temp_upw + den_kg_adv*wdn)
      dh_dC_dn = c_p * temp_upw * ddkg_dC_dn
    endif

    ! diffusive transmissibility derivatives
    dGamma_D_dpdn = 0.d0; dGamma_D_dTdn = 0.d0; dGamma_D_dCdn = 0.d0
    if (dispersion_scale > 0.d0) then
      dD_hyd_dn_dp = (thc_auxvar_dn%dsat_dp * &
                      thc_auxvar_dn%effective_porosity + &
                      thc_auxvar_dn%sat * thc_auxvar_dn%dpor_dp) * &
                     material_auxvar_dn%tortuosity * thc_auxvar_dn%diff_mol
      dD_hyd_dn_dT = thc_auxvar_dn%effective_porosity * &
                     thc_auxvar_dn%sat * material_auxvar_dn%tortuosity * &
                     thc_auxvar_dn%ddiff_dT
      dD_hyd_dn_dC = thc_auxvar_dn%effective_porosity * &
                     thc_auxvar_dn%sat * material_auxvar_dn%tortuosity * &
                     thc_auxvar_dn%ddiff_dC
      dGamma_D_dpdn = area * dispersion_scale * dD_hyd_dn_dp / dist(0)
      dGamma_D_dTdn = area * dispersion_scale * dD_hyd_dn_dT / dist(0)
      dGamma_D_dCdn = area * dispersion_scale * dD_hyd_dn_dC / dist(0)
    endif

    ! ----------------------------------------------------------------------- !
    ! Flow equation  F_flow = rho_kmol_adv * q
    ! ----------------------------------------------------------------------- !
    Jdn(thc_pressure_dof,thc_pressure_dof) = &
      ddkmol_dp_dn*q + den_kmol_adv*dq_dpdn
    Jdn(thc_pressure_dof,thc_temperature_dof) = &
      ddkmol_dT_dn*q + den_kmol_adv*dq_dTdn
    Jdn(thc_pressure_dof,thc_concentration_dof) = &
      ddkmol_dC_dn*q + den_kmol_adv*dq_dCdn

    ! ----------------------------------------------------------------------- !
    ! Solute equation  F_sol = (C_upw*q + Gamma_D*dC)*1000
    !   d(C_upw)/dC_dn = wdn ; d(dC)/dC_dn = -1
    ! ----------------------------------------------------------------------- !
    Jdn(thc_concentration_dof,thc_pressure_dof) = &
      (conc_upw*dq_dpdn + dGamma_D_dpdn*delta_conc) * L_per_m3
    Jdn(thc_concentration_dof,thc_temperature_dof) = &
      (conc_upw*dq_dTdn + dGamma_D_dTdn*delta_conc) * L_per_m3
    Jdn(thc_concentration_dof,thc_concentration_dof) = &
      (wdn*q + conc_upw*dq_dCdn + dGamma_D_dCdn*delta_conc - Gamma_D) * L_per_m3

    ! ----------------------------------------------------------------------- !
    ! Energy equation  F_energy = h_upw*q + Gamma_kappa*dT  (+ Robin)
    !   d(dT)/dT_dn = -1 ; Gamma_kappa indep of C
    ! ----------------------------------------------------------------------- !
    Jdn(thc_temperature_dof,thc_pressure_dof) = &
      dh_dp_dn*q + h_upw*dq_dpdn + dGamma_kappa_dpdn*delta_temp
    Jdn(thc_temperature_dof,thc_temperature_dof) = &
      dh_dT_dn*q + h_upw*dq_dTdn - Gamma_kappa + &
      dGamma_kappa_dTdn*delta_temp
    Jdn(thc_temperature_dof,thc_concentration_dof) = &
      dh_dC_dn*q + h_upw*dq_dCdn
    ! convective (Robin) thermal BC: d/dT_cell of h_conv*(T_ext - T_cell)
    if (energy_bc_type == THC_CONVECTIVE_BC) then
      h_conv = auxvars(auxvar_mapping(THC_COND_WATER_AUX_INDEX))
      Jdn(thc_temperature_dof,thc_temperature_dof) = &
        Jdn(thc_temperature_dof,thc_temperature_dof) - h_conv * area
    endif

  endif

  if (debug_connection) then
    write(*,'(9x,"BFF(res_f,kr,dp): ",8es12.4)') Res(thc_pressure_dof), &
                                                 kr, delta_pressure
    write(*,'(9x,"BFE(res_e,Gk,dT): ",8es12.4)') Res(thc_temperature_dof), &
                                                 Gamma_kappa, delta_temp
    write(*,'(9x,"BFS(res_s,q,dC): ",8es12.4)') Res(thc_concentration_dof),&
                                                q, delta_conc
  endif

end subroutine THCBCFlux

! ************************************************************************** !

subroutine THCSrcSink(option,flow_aux_real_var,flow_src_sink_mapping, &
                          flow_src_sink_type, rate_aux_real, &
                          thc_auxvar,global_auxvar,material_auxvar, &
                          ss_flow_vol_flux,Res,Jdn,dResdparamdn, &
                          calculate_derivatives)
  !
  ! Computes the source/sink terms (flow / energy / solute) for the THC
  ! residual and its 3x3 Jacobian.  Mirrors ZFlowSrcSink (zflow_common.F90),
  ! extended to the energy equation and to a solution-dependent density.
  !
  ! Injection vs extraction is by the sign of the volumetric rate q_vol
  !   q_vol >= 0 (injection):  solute at C_inj, energy at T_inj (both fixed) ;
  !   q_vol <  0 (extraction): solute at C_cell, energy at T_cell (interior) .
  !
  ! dResdparamdn is retained for interface compatibility but is a NO-OP here.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/09/26
  !
  use Option_module
  use Material_Aux_module
  use EOS_Water_module
  use Utility_module, only : Smoothstep

  implicit none

  type(option_type) :: option
  PetscReal :: flow_aux_real_var(:)
  PetscInt :: flow_src_sink_mapping(:)
  PetscInt :: flow_src_sink_type(:)
  PetscReal :: rate_aux_real(:)
  type(thc_auxvar_type) :: thc_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: ss_flow_vol_flux
  PetscReal :: Res(THC_NDOF)
  PetscReal :: Jdn(THC_NDOF,THC_NDOF)
  PetscReal :: dResdparamdn(THC_NDOF,THC_NDOF)
  PetscBool :: calculate_derivatives

  PetscReal, parameter :: L_per_m3 = 1.d3
  PetscReal :: qsrc_m3, qsrc_L
  PetscReal :: qsrc0, qsrc_kmol, qsrc_mass
  PetscReal :: dqsrc_kmol_dp, dqsrc_kmol_dT, dqsrc_kmol_dC
  PetscReal :: dqsrc_mass_dp, dqsrc_mass_dT, dqsrc_mass_dC
  PetscReal :: threshold_pressure, pressure_span
  PetscReal :: min_pressure, max_pressure, volume_scale
  PetscReal :: smoothstep_scale, dsmoothstep_dp, dqsrc_kmol_dsmoothstep
  PetscBool :: inhibit_flow_above_pressure
  PetscReal :: c_p, conc_used, temp_used
  PetscReal :: h_used, dh_used_dP, dh_used_dT, dh_used_dC
  PetscReal :: hw, hw_dp, hw_dT
  PetscErrorCode :: ierr
  PetscInt :: dof_index, index_

  Res = 0.d0
  Jdn = 0.d0
  dResdparamdn = 0.d0

  c_p = thc_specific_heat_liquid       ! c_p,l [J/(kg.K)]

  ! ========================================================================= !
  ! Flow (liquid mass) source/sink : units kmol/s
  !   Volumetric-input cases : qsrc_kmol = q_vol * rho_kmol(P,T,C)
  !                            -> density-dependent (dden_* derivatives)
  !   Mass-input cases       : qsrc_kmol = q_mass / M_w   (fixed molar rate)
  !                            -> NO density derivative (avoids double-density)
  ! ========================================================================= !
  dof_index = flow_src_sink_mapping(THC_COND_WATER_INDEX)
  qsrc0 = flow_aux_real_var(dof_index)
  qsrc_m3 = 0.d0
  qsrc_kmol = 0.d0
  dqsrc_kmol_dp = 0.d0
  dqsrc_kmol_dT = 0.d0
  dqsrc_kmol_dC = 0.d0
  select case(flow_src_sink_type(dof_index))
    case(VOLUMETRIC_RATE_SS)
      ! qsrc0 is the volumetric rate [m^3/s]; convert to kmol/s
      qsrc_m3 = qsrc0
      qsrc_kmol = qsrc_m3 * thc_auxvar%den_kmol
      dqsrc_kmol_dp = qsrc_m3 * thc_auxvar%dden_dp / FMWH2O
      dqsrc_kmol_dT = qsrc_m3 * thc_auxvar%dden_dT / FMWH2O
      dqsrc_kmol_dC = qsrc_m3 * thc_auxvar%dden_dC / FMWH2O
    case(SCALED_VOLUMETRIC_RATE_SS)
      index_ = flow_src_sink_mapping(THC_COND_WATER_AUX_INDEX)
      ! qsrc0 is the volumetric rate [m^3/s]; convert to kmol/s
      qsrc_m3 = qsrc0 * flow_aux_real_var(index_)
      qsrc_kmol = qsrc_m3 * thc_auxvar%den_kmol
      dqsrc_kmol_dp = qsrc_m3 * thc_auxvar%dden_dp / FMWH2O
      dqsrc_kmol_dT = qsrc_m3 * thc_auxvar%dden_dT / FMWH2O
      dqsrc_kmol_dC = qsrc_m3 * thc_auxvar%dden_dC / FMWH2O
    case(MASS_RATE_SS)
      ! qsrc0 is the mass rate [kg/s]; convert to kmol/s
      ! [kg/s -> kmol/s; FMWH2O -> g/mol = kg/kmol]
      qsrc_kmol = qsrc0 / FMWH2O
      qsrc_m3 = qsrc_kmol / thc_auxvar%den_kmol
    case(SCALED_MASS_RATE_SS)
      index_ = flow_src_sink_mapping(THC_COND_WATER_AUX_INDEX)
      ! qsrc0 is the mass rate [kg/s]; convert to kmol/s
      qsrc_kmol = qsrc0 / FMWH2O * flow_aux_real_var(index_)
      qsrc_m3 = qsrc_kmol / thc_auxvar%den_kmol
    case(PRES_REG_MASS_RATE_SS)
      ! pressure-regulated mass rate: a smoothstep ramp on cell pressure scales
      ! the injected/extracted mass rate (mirrors TH's PRES_REG_MASS_RATE_SS).
      threshold_pressure = rate_aux_real(1)
      inhibit_flow_above_pressure = (threshold_pressure > 0.d0)
      threshold_pressure = dabs(threshold_pressure)
      pressure_span = rate_aux_real(2)
      if (inhibit_flow_above_pressure) then
        min_pressure = threshold_pressure - pressure_span
        max_pressure = threshold_pressure
      else
        min_pressure = threshold_pressure
        max_pressure = threshold_pressure + pressure_span
      endif
      index_ = flow_src_sink_mapping(THC_COND_WATER_AUX_INDEX)
      volume_scale = flow_aux_real_var(index_)
      call Smoothstep(thc_auxvar%pres,min_pressure,max_pressure, &
                      smoothstep_scale,dsmoothstep_dp)
      if (inhibit_flow_above_pressure) then
        smoothstep_scale = 1.d0 - smoothstep_scale
        dsmoothstep_dp = -1.d0 * dsmoothstep_dp
      endif
      dqsrc_kmol_dsmoothstep = qsrc0 / FMWH2O * volume_scale
      qsrc_kmol = dqsrc_kmol_dsmoothstep * smoothstep_scale
      dqsrc_kmol_dp = dqsrc_kmol_dsmoothstep * dsmoothstep_dp
      qsrc_m3 = qsrc_kmol / thc_auxvar%den_kmol
    case default
      option%io_buffer = 'src_sink_type not supported in THCSrcSink'
      call PrintErrMsg(option)
  end select
  ss_flow_vol_flux = qsrc_m3
  qsrc_L = qsrc_m3 * L_per_m3

  ! water mass source rate [kg/s] used by the energy equation
  qsrc_mass = qsrc_kmol * FMWH2O
  dqsrc_mass_dp = dqsrc_kmol_dp * FMWH2O
  dqsrc_mass_dT = dqsrc_kmol_dT * FMWH2O
  dqsrc_mass_dC = dqsrc_kmol_dC * FMWH2O

  ! Res[kmol/sec]
  Res(thc_pressure_dof) = qsrc_kmol

  if (calculate_derivatives) then
    ! dQ_flow/dX : density-dependent for volumetric input, smoothstep-dependent
    ! for the pressure-regulated case, zero for the fixed mass-rate cases.
    Jdn(thc_pressure_dof,thc_pressure_dof) = dqsrc_kmol_dp
    Jdn(thc_pressure_dof,thc_temperature_dof) = dqsrc_kmol_dT
    Jdn(thc_pressure_dof,thc_concentration_dof) = dqsrc_kmol_dC
  endif

  ! ========================================================================= !
  ! Solute source/sink : units mol/s
  ! ========================================================================= !
  if (qsrc_m3 >= 0.d0) then
    ! injection at prescribed C_inj
    dof_index = flow_src_sink_mapping(THC_COND_SOLUTE_INDEX)
    conc_used = flow_aux_real_var(dof_index)
  else
    ! extraction carries interior cell concentration
    conc_used = thc_auxvar%conc
  endif
  Res(thc_concentration_dof) = qsrc_L * conc_used
  if (calculate_derivatives) then
    if (qsrc_m3 < 0.d0) then
      ! dQ_sol/dC = q_vol*1000   (C_cell depends on solution); 0 for injection
      Jdn(thc_concentration_dof,thc_concentration_dof) = qsrc_L
    endif
    ! Note: for the mass-input cases qsrc_m3 = qsrc_kmol/rho_kmol carries a weak
    ! density dependence on (P,T,C); those second-order solute cross-derivatives
    ! are intentionally omitted (they only affect Newton convergence, not the
    ! converged solution, and TH carries no solute term at all).
  endif

  ! ========================================================================= !
  ! Energy source/sink : units W
  !   RHO_CP_T: Q_energy = q_mass(P,T,C) * c_p * T_used
  !   FULL_EOS: Q_energy = q_mass(P,T,C) * h_spec_used
  !     injection  -> h_spec at prescribed T_inj and cell pressure
  !     extraction -> h_spec at cell conditions (thc_auxvar%h)
  ! ========================================================================= !
  if (qsrc_m3 >= 0.d0) then
    ! injection at prescribed T_inj (fixed)
    dof_index = flow_src_sink_mapping(THC_COND_ENERGY_INDEX)
    temp_used = flow_aux_real_var(dof_index)
  else
    ! extraction carries interior cell temperature
    temp_used = thc_auxvar%temp
  endif

  if (thc_energy_mode == THC_ENERGY_FULL_EOS) then
    ! specific enthalpy h_used [J/kg] and its derivatives w.r.t. cell DOFs
    if (qsrc_m3 >= 0.d0) then
      ! injection: h(T_inj, P_cell); independent of cell T and C, varies with P
      ierr = 0
      call EOSWaterEnthalpy(temp_used,thc_auxvar%pres,hw,hw_dp,hw_dT,ierr)
      if (ierr /= 0) then
        option%io_buffer = 'Error computing THC injection enthalpy in &
          &THCSrcSink (FULL_EOS).'
        call PrintErrMsg(option)
      endif
      h_used     = hw / FMWH2O
      dh_used_dP = hw_dp / FMWH2O
      dh_used_dT = 0.d0
      dh_used_dC = 0.d0
    else
      ! extraction: cell specific enthalpy and its stored derivatives
      h_used     = thc_auxvar%h
      dh_used_dP = thc_auxvar%dh_dP
      dh_used_dT = thc_auxvar%dh_dT
      dh_used_dC = thc_auxvar%dh_dC
    endif
    Res(thc_temperature_dof) = qsrc_mass * h_used
    if (calculate_derivatives) then
      ! product rule on q_mass(P,T,C) * h_used(P[,T,C])
      Jdn(thc_temperature_dof,thc_pressure_dof) = &
        dqsrc_mass_dp * h_used + qsrc_mass * dh_used_dP
      Jdn(thc_temperature_dof,thc_temperature_dof) = &
        dqsrc_mass_dT * h_used + qsrc_mass * dh_used_dT
      Jdn(thc_temperature_dof,thc_concentration_dof) = &
        dqsrc_mass_dC * h_used + qsrc_mass * dh_used_dC
    endif
  else
    Res(thc_temperature_dof) = qsrc_mass * c_p * temp_used
    if (calculate_derivatives) then
      ! q_mass-coupling terms (present for both injection and extraction)
      Jdn(thc_temperature_dof,thc_pressure_dof) = &
        dqsrc_mass_dp * c_p * temp_used
      Jdn(thc_temperature_dof,thc_concentration_dof) = &
        dqsrc_mass_dC * c_p * temp_used
      if (qsrc_m3 >= 0.d0) then
        ! injection: T_inj fixed -> only q_mass's T-dependence
        Jdn(thc_temperature_dof,thc_temperature_dof) = &
          dqsrc_mass_dT * c_p * temp_used
      else
        ! extraction: T_used = T_cell -> product rule on q_mass*T_cell
        Jdn(thc_temperature_dof,thc_temperature_dof) = &
          c_p * (dqsrc_mass_dT * temp_used + qsrc_mass)
      endif
    endif
  endif

end subroutine THCSrcSink

! ************************************************************************** !

subroutine THCAccumDerivative(thc_auxvar,global_auxvar, &
                                  material_auxvar, &
                                  option,Res,Jac,dResdparam)
  !
  ! Computes derivatives of the accumulation term for the Jacobian.
  !
  ! When thc_numerical_derivatives is .false. (default), the analytical 3x3
  ! accumulation block computed by THCAccumulation is returned unchanged.
  ! When .true., the block is formed by first-order forward finite differences
  ! of the residual w.r.t. each perturbed primary variable (P,T,C):
  !   Jac(ieq,idof) = [Res(x+pert*e_idof) - Res(x)] / pert .
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/23/26
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(thc_auxvar_type) :: thc_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(THC_NDOF)
  PetscReal :: Jac(THC_NDOF,THC_NDOF)
  PetscReal :: dResdparam(THC_NDOF,THC_NDOF)

  PetscReal :: res_pert(THC_NDOF)
  PetscReal :: Jdum(THC_NDOF,THC_NDOF)
  PetscReal :: dJdum(THC_NDOF,THC_NDOF)
  PetscInt :: idof, ieq

  call THCAccumulation(thc_auxvar(ZERO_INTEGER), &
                           global_auxvar, &
                           material_auxvar, &
                           option,Res,Jac,dResdparam, &
                           .not.thc_numerical_derivatives)

  if (thc_numerical_derivatives) then
    do idof = 1, option%nflowdof
      call THCAccumulation(thc_auxvar(idof), &
                               global_auxvar, &
                               material_auxvar, &
                               option,res_pert,Jdum,dJdum, &
                               PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jac(ieq,idof) = (res_pert(ieq)-Res(ieq))/thc_auxvar(idof)%pert
      enddo
    enddo
  endif

end subroutine THCAccumDerivative

! ************************************************************************** !

subroutine THCFluxDerivative(thc_auxvar_up,global_auxvar_up, &
                                 material_auxvar_up, &
                                 thc_auxvar_dn,global_auxvar_dn, &
                                 material_auxvar_dn, &
                                 area,dist,thc_parameter,option,v_darcy, &
                                 Res,Jup,Jdn,dResdparamup,dResdparamdn, &
                                 debug_connection)
  !
  ! Computes the derivatives of the internal flux terms for the Jacobian.
  !
  ! Analytical path (thc_numerical_derivatives = .false.): returns the
  ! analytical Jup/Jdn from THCFlux.  Numerical path: perturbs upgradient
  ! DOFs (-> Jup) and then the downgradient DOFs (-> Jdn) one at a time & forms
  ! forward differences of the residual.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/23/26
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(thc_auxvar_type) :: thc_auxvar_up(0:), thc_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(thc_parameter_type) :: thc_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(THC_NDOF)
  PetscReal :: Jup(THC_NDOF,THC_NDOF)
  PetscReal :: Jdn(THC_NDOF,THC_NDOF)
  PetscReal :: dResdparamup(THC_NDOF,THC_NDOF)
  PetscReal :: dResdparamdn(THC_NDOF,THC_NDOF)
  PetscBool :: debug_connection

  PetscReal :: res_pert(THC_NDOF)
  PetscReal :: v_darcy_dum
  PetscReal :: Jdum(THC_NDOF,THC_NDOF)
  PetscReal :: dJdum(THC_NDOF,THC_NDOF)
  PetscInt :: idof, ieq

  Jup = 0.d0
  Jdn = 0.d0

  call THCFlux(thc_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                   material_auxvar_up, &
                   thc_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                   material_auxvar_dn, &
                   area,dist,thc_parameter,option,v_darcy, &
                   Res,Jup,Jdn,dResdparamup,dResdparamdn, &
                   .not.thc_numerical_derivatives, &
                   debug_connection)

  if (thc_numerical_derivatives) then
    ! upgradient derivatives
    do idof = 1, option%nflowdof
      call THCFlux(thc_auxvar_up(idof),global_auxvar_up, &
                       material_auxvar_up, &
                       thc_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                       material_auxvar_dn, &
                       area,dist,thc_parameter,option,v_darcy_dum, &
                       res_pert,Jdum,Jdum,dJdum,dJdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jup(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        thc_auxvar_up(idof)%pert
      enddo
    enddo
    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call THCFlux(thc_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                       material_auxvar_up, &
                       thc_auxvar_dn(idof),global_auxvar_dn, &
                       material_auxvar_dn, &
                       area,dist,thc_parameter,option,v_darcy_dum, &
                       res_pert,Jdum,Jdum,dJdum,dJdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jdn(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        thc_auxvar_dn(idof)%pert
      enddo
    enddo
  endif

end subroutine THCFluxDerivative

! ************************************************************************** !

subroutine THCBCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                                   thc_auxvar_up, &
                                   global_auxvar_up, &
                                   thc_auxvar_dn,global_auxvar_dn, &
                                   material_auxvar_dn, &
                                   area,dist,thc_parameter,option,v_darcy, &
                                   Res,Jdn,dResdparamdn,debug_connection)
  !
  ! Computes the derivatives of the boundary flux terms for the Jacobian.
  ! Only the interior ("dn") cell DOFs are perturbed; the "up" cell is the fixed
  ! boundary ghost.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/23/26
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(:)
  PetscInt :: auxvar_mapping(THC_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(thc_auxvar_type) :: thc_auxvar_up, thc_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(thc_parameter_type) :: thc_parameter
  PetscReal :: v_darcy
  PetscReal :: Res(THC_NDOF)
  PetscReal :: Jdn(THC_NDOF,THC_NDOF)
  PetscReal :: dResdparamdn(THC_NDOF,THC_NDOF)
  PetscBool :: debug_connection

  PetscReal :: res_pert(THC_NDOF)
  PetscReal :: v_darcy_dum
  PetscReal :: Jdum(THC_NDOF,THC_NDOF)
  PetscReal :: dJdum(THC_NDOF,THC_NDOF)
  PetscInt :: idof, ieq

  Jdn = 0.d0

  call THCBCFlux(ibndtype,auxvar_mapping,auxvars, &
                     thc_auxvar_up,global_auxvar_up, &
                     thc_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                     material_auxvar_dn, &
                     area,dist,thc_parameter,option,v_darcy, &
                     Res,Jdn,dResdparamdn,.not.thc_numerical_derivatives, &
                     debug_connection)

  if (thc_numerical_derivatives) then
    ! downgradient (interior cell) derivatives
    do idof = 1, option%nflowdof
      call THCBCFlux(ibndtype,auxvar_mapping,auxvars, &
                         thc_auxvar_up,global_auxvar_up, &
                         thc_auxvar_dn(idof),global_auxvar_dn, &
                         material_auxvar_dn, &
                         area,dist,thc_parameter,option,v_darcy_dum, &
                         res_pert,Jdum,dJdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jdn(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        thc_auxvar_dn(idof)%pert
      enddo
    enddo
  endif

end subroutine THCBCFluxDerivative

! ************************************************************************** !

subroutine THCSrcSinkDerivative(option,flow_aux_real_var, &
                                    flow_src_sink_mapping, &
                                    flow_src_sink_type, rate_aux_real, &
                                    thc_auxvar,global_auxvar, &
                                    material_auxvar, &
                                    ss_flow_vol_flux, &
                                    Res,Jac,dResdparam)
  !
  ! Computes the derivatives of the source/sink terms for the Jacobian.
  !
  ! Author: Piyoosh Jaysaval
  ! Date:   06/23/26
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  PetscReal :: flow_aux_real_var(:)
  PetscInt :: flow_src_sink_mapping(:)
  PetscInt :: flow_src_sink_type(:)
  PetscReal :: rate_aux_real(:)
  type(thc_auxvar_type) :: thc_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: ss_flow_vol_flux
  PetscReal :: Res(THC_NDOF)
  PetscReal :: Jac(THC_NDOF,THC_NDOF)
  PetscReal :: dResdparam(THC_NDOF,THC_NDOF)

  PetscReal :: res_pert(THC_NDOF)
  PetscReal :: dummy_real
  PetscReal :: Jdum(THC_NDOF,THC_NDOF)
  PetscReal :: dJdum(THC_NDOF,THC_NDOF)
  PetscInt :: idof, ieq

  Jac = 0.d0
  ! unperturbed thc_auxvar's value
  call THCSrcSink(option,flow_aux_real_var,flow_src_sink_mapping, &
                      flow_src_sink_type, rate_aux_real, &
                      thc_auxvar(ZERO_INTEGER),global_auxvar, &
                      material_auxvar, &
                      ss_flow_vol_flux, &
                      Res,Jac,dResdparam,.not.thc_numerical_derivatives)

  if (thc_numerical_derivatives) then
    ! perturbed thc_auxvar's value
    do idof = 1, option%nflowdof
      call THCSrcSink(option,flow_aux_real_var,flow_src_sink_mapping, &
                          flow_src_sink_type, rate_aux_real, &
                          thc_auxvar(idof),global_auxvar, &
                          material_auxvar, &
                          dummy_real, &
                          res_pert,Jdum,dJdum,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jac(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        thc_auxvar(idof)%pert
      enddo
    enddo
  endif

end subroutine THCSrcSinkDerivative

! ************************************************************************** !

end module THC_Common_module
