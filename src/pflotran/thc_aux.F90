module THC_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module
  use Material_Aux_module
  use EOS_Water_module
  use THC_EOS_Utils_module
  use Characteristic_Curves_Thermal_module
  use Utility_module, only : Arrhenius

  implicit none

  private

  ! Solid / fluid thermal properties
  PetscReal, public :: thc_specific_heat_liquid = 4184.d0  ! c_p,l [J/(kg.K)]
  PetscReal, public :: thc_specific_heat_solid = 800.d0   ! c_s   [J/(kg.K)]
  PetscReal, public :: thc_density_solid       = 2650.d0  ! rho_s [kg/m^3]
  ! grain (solid) thermal conductivity used by the porosity-dependent
  ! grain-Somerton (THCThermalConductivityEff)
  PetscReal, public :: thc_kappa_solid         = 3.d0     ! kappa_s [W/(m.K)]
  ! optional OPTIONS-block global BULK dry/wet conductivities; when set they
  ! feed the TH-style interpolation (kappa_eff = ckdry + sqrt(S)*(ckwet-ckdry))
  ! and take precedence over thc_kappa_solid.  Default UNINITIALIZED so the
  ! grain-Somerton path remains the fallback when they are not supplied.
  PetscReal, public :: thc_kappa_dry = UNINITIALIZED_DOUBLE  ! [W/(m.K)]
  PetscReal, public :: thc_kappa_wet = UNINITIALIZED_DOUBLE  ! [W/(m.K)]
  ! reference temperature for the Arrhenius diffusion scaling (matches RT)
  PetscReal, parameter, public :: thc_diffusion_ref_temp = 25.d0  ! [C]

  ! Energy formulation switch (selectable via ENERGY_FORMULATION input card).
  !   THC_ENERGY_RHO_CP_T (default): liquid energy = rho*c_p*T
  !   THC_ENERGY_FULL_EOS (TH-equivalent)  : liquid energy uses EOS u(P,T)
  !                                              in accumulation and h(P,T) in
  !                                              advection (single-phase liquid)
  ! Mass-specific (J/kg, J/m^3) units throughout -- no option%scale applied,
  ! unlike TH which carries molar MJ. The rock term stays rho_s*c_s*T in both.
  PetscInt, parameter, public :: THC_ENERGY_RHO_CP_T = 1
  PetscInt, parameter, public :: THC_ENERGY_FULL_EOS = 2
  PetscInt, public :: thc_energy_mode = THC_ENERGY_RHO_CP_T

  PetscInt, parameter, public :: THC_ADVECTIVE_DENSITY_UPWIND = 1
  PetscInt, parameter, public :: THC_ADVECTIVE_DENSITY_TH_COMPATIBLE = 2
  PetscInt, public :: thc_advective_density_mode = &
                        THC_ADVECTIVE_DENSITY_UPWIND

  ! perturbation controls (mirror zflow_aux.F90)
  PetscReal, public :: thc_rel_pert = 1.d-8
  PetscReal, public :: thc_pres_min_pert = 1.d-2
  PetscReal, public :: thc_temp_min_pert = 1.d-6 ! not based on anything
  PetscReal, public :: thc_conc_min_pert = 1.d-6 ! not based on anything

  PetscReal, pointer, public :: thc_min_pert(:)

  PetscBool, public :: thc_calc_accum = PETSC_TRUE
  PetscBool, public :: thc_calc_flux = PETSC_TRUE
  PetscBool, public :: thc_calc_bcflux = PETSC_TRUE

  PetscBool, public :: thc_numerical_derivatives = PETSC_FALSE
  PetscBool, public :: thc_simult_function_evals = PETSC_TRUE
  PetscBool, public :: thc_tensorial_rel_perm = PETSC_FALSE
  PetscBool, public :: thc_acknowledge_no_compress = PETSC_FALSE

  ! debugging
  PetscInt, public :: thc_ni_count
  PetscInt, public :: thc_ts_cut_count
  PetscInt, public :: thc_ts_count
  PetscInt, public :: thc_debug_cell_id

  ! Fixed DOF indices -- all always active
  PetscInt, parameter, public :: thc_pressure_dof      = 1
  PetscInt, parameter, public :: thc_temperature_dof   = 2
  PetscInt, parameter, public :: thc_concentration_dof = 3
  PetscInt, parameter, public :: THC_NDOF = 3

  PetscInt, parameter, public :: THC_COND_WATER_INDEX = 1
  PetscInt, parameter, public :: THC_COND_ENERGY_INDEX = 2
  PetscInt, parameter, public :: THC_COND_SOLUTE_INDEX = 3
  PetscInt, parameter, public :: THC_COND_WATER_AUX_INDEX = 4
  PetscInt, parameter, public :: THC_MAX_INDEX = 4

  PetscInt, parameter, public :: THC_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: THC_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: THC_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: THC_UPDATE_FOR_BOUNDARY = 2

  PetscInt, parameter, public :: THC_LIQ_SAT_WRT_LIQ_PRES = 1
  PetscInt, parameter, public :: THC_LIQ_PRES_WRT_POROS = 2

  type, public :: thc_auxvar_type
    PetscReal :: pres ! liquid pressure P_l [Pa]
    PetscReal :: sat  ! liquid saturation S_l [-]
    PetscReal :: pc   ! capillary pressure [Pa]
    PetscReal :: kr   ! relative permeability [-]
    PetscReal :: effective_porosity   ! phi_eff [-]
    PetscReal :: dpor_dp              ! dphi/dP_l [1/Pa]
    PetscReal :: dsat_dp ! derivative of saturation wrt pressure [1/Pa]
    PetscReal :: dkr_dp  ! derivative of rel. perm. wrt pressure [1/Pa]
    PetscReal :: effective_saturation
    PetscReal :: deffsat_dp
    PetscReal :: conc ! concentration C [mol/L]
    PetscReal :: pert
    PetscReal :: mat_pert(1)

    PetscReal :: temp             ! T [C]
    PetscReal :: den_kg           ! rho_l [kg/m^3]
    PetscReal :: den_kmol         ! rho_l / M_w [kmol/m^3]
    PetscReal :: dden_dp          ! d(rho_l)/dP_l [kg/(m^3.Pa)]
    PetscReal :: dden_dT          ! d(rho_l)/dT [kg/(m^3.K)]
    PetscReal :: dden_dC          ! d(rho_l)/dC [kg/(m^3.(mol/L))]
    PetscReal :: vis              ! mu_l [Pa.s]
    PetscReal :: dvis_dT          ! d(mu_l)/dT [Pa.s/K]
    PetscReal :: dvis_dC          ! d(mu_l)/dC [Pa.s.L/mol]
    PetscReal :: diff_mol         ! D_mol [m^2/s]
    PetscReal :: ddiff_dT         ! d(D_mol)/dT [m^2/(s.K)]
    PetscReal :: ddiff_dC         ! d(D_mol)/dC [m^2.L/(s.mol)]
    PetscReal :: therm_cond_eff   ! kappa_eff [W/(m.K)]
    PetscReal :: dtherm_cond_dsat ! d(kappa_eff)/dS_l [W/(m.K)]
    PetscReal :: dtherm_cond_dT   ! d(kappa_eff)/dT [W/(m.K.C)]
    PetscReal :: heat_cap_liquid  ! rho_l * c_p,l [J/(m^3.K)]
    PetscReal :: heat_cap_solid   ! rho_s * c_s   [J/(m^3.K)]
    ! Full-EOS energy fields (used only when thc_energy_mode==FULL_EOS).
    ! Mass-specific (per kg) -- converted from molar EOS values via FMWH2O.
    PetscReal :: u                ! specific internal energy u(P,T) [J/kg]
    PetscReal :: h                ! specific enthalpy h(P,T)        [J/kg]
    PetscReal :: du_dP            ! d(u)/dP at const T,C [J/(kg.Pa)]
    PetscReal :: du_dT            ! d(u)/dT at const P,C [J/(kg.K)]
    PetscReal :: du_dC            ! d(u)/dC at const P,T [J.L/(kg.mol)]
    PetscReal :: dh_dP            ! d(h)/dP at const T,C [J/(kg.Pa)]
    PetscReal :: dh_dT            ! d(h)/dT at const P,C [J/(kg.K)]
    PetscReal :: dh_dC            ! d(h)/dC at const P,T [J.L/(kg.mol)]
  end type thc_auxvar_type

  type, public :: thc_parameter_type
    PetscBool :: check_post_converged
    PetscReal, pointer :: tensorial_rel_perm_exponent(:,:)
    PetscReal :: diffusion_coef
    ! FLUID_PROPERTY DIFFUSION_ACTIVATION_ENERGY [J/mol]; when Initialized the
    ! molecular diffusion is scaled by Arrhenius(AE,T,25C), as in RT/NWT
    PetscReal :: diffusion_activation_energy
    ! Per-material thermal properties (indexed by material id).  Populated in
    ! THCSetup from the MATERIAL_PROPERTY card (TH-compatible), falling back
    ! to the MODE THC OPTIONS-block globals when a material omits them.
    PetscReal, pointer :: dencpr(:)  ! rho_s * c_s [J/(m^3.K)]
    PetscReal, pointer :: ckdry(:)   ! bulk dry  thermal conductivity [W/(m.K)]
    PetscReal, pointer :: ckwet(:)   ! bulk wet  thermal conductivity [W/(m.K)]
    ! per-material THERMAL_CHARACTERISTIC_CURVES; takes precedence over
    ! ckdry/ckwet and is the only path supplying d(kappa_eff)/dT
    type(cc_thermal_ptr_type), pointer :: thermal_cc(:)
    ! When ckdry(imat)/ckwet(imat) are UNINITIALIZED_DOUBLE the per-cell
    ! conductivity falls back to the grain-Somerton model using the global
    ! thc_kappa_solid; otherwise the TH-style bulk interpolation
    ! kappa_eff = ckdry + sqrt(S_l)*(ckwet-ckdry) is used.
  end type thc_parameter_type

  type, public :: thc_type
    PetscBool :: auxvars_up_to_date
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(thc_parameter_type), pointer :: thc_parameter
    type(thc_auxvar_type), pointer :: auxvars(:,:)
    type(thc_auxvar_type), pointer :: auxvars_bc(:)
    type(thc_auxvar_type), pointer :: auxvars_ss(:)
    type(material_auxvar_type), pointer :: material_auxvars_pert(:,:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type thc_type

  interface THCAuxVarDestroy
    module procedure THCAuxVarSingleDestroy
    module procedure THCAuxVarArray1Destroy
    module procedure THCAuxVarArray2Destroy
  end interface THCAuxVarDestroy

  interface THCOutputAuxVars
    module procedure THCOutputAuxVars1
  end interface THCOutputAuxVars

  public :: THCAuxCreate, &
            THCAuxDestroy, &
            THCAuxVarCompute, &
            THCAuxVarInit, &
            THCAuxVarCopy, &
            THCAuxVarDestroy, &
            THCAuxVarStrip, &
            THCAuxVarPerturb, &
            THCAuxMapConditionIndices, &
            THCPrintAuxVars, &
            THCOutputAuxVars, &
            THCAuxTensorialRelPerm, &
            THCThermalCondEval

contains

! ************************************************************************** !

function THCAuxCreate(option)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Option_module

  implicit none

  type(option_type) :: option

  type(thc_type), pointer :: THCAuxCreate

  type(thc_type), pointer :: aux

  nullify(thc_min_pert)

  allocate(aux)
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  nullify(aux%material_auxvars_pert)
  nullify(aux%matrix_zeroing)

  allocate(aux%thc_parameter)
  aux%thc_parameter%check_post_converged = PETSC_FALSE
  nullify(aux%thc_parameter%tensorial_rel_perm_exponent)
  aux%thc_parameter%diffusion_coef = 0.d0
  aux%thc_parameter%diffusion_activation_energy = UNINITIALIZED_DOUBLE
  nullify(aux%thc_parameter%dencpr)
  nullify(aux%thc_parameter%ckdry)
  nullify(aux%thc_parameter%ckwet)
  nullify(aux%thc_parameter%thermal_cc)

  THCAuxCreate => aux

end function THCAuxCreate

! ************************************************************************** !

subroutine THCAuxVarInit(auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Option_module

  implicit none

  type(thc_auxvar_type) :: auxvar
  type(option_type) :: option

  auxvar%pres = 0.d0
  auxvar%sat = 0.d0
  auxvar%pc = 0.d0
  auxvar%kr = 0.d0
  auxvar%effective_porosity = 0.d0
  auxvar%dpor_dp = 0.d0
  auxvar%dsat_dp = 0.d0
  auxvar%dkr_dp = 0.d0
  auxvar%effective_saturation = UNINITIALIZED_DOUBLE
  auxvar%deffsat_dp = UNINITIALIZED_DOUBLE
  auxvar%conc = 0.d0
  auxvar%pert = 0.d0
  auxvar%mat_pert = 0.d0

  auxvar%temp = 0.d0
  auxvar%den_kg = 0.d0
  auxvar%den_kmol = 0.d0
  auxvar%dden_dp = 0.d0
  auxvar%dden_dT = 0.d0
  auxvar%dden_dC = 0.d0
  auxvar%vis = 0.d0
  auxvar%dvis_dT = 0.d0
  auxvar%dvis_dC = 0.d0
  auxvar%diff_mol = 0.d0
  auxvar%ddiff_dT = 0.d0
  auxvar%ddiff_dC = 0.d0
  auxvar%therm_cond_eff = 0.d0
  auxvar%dtherm_cond_dsat = 0.d0
  auxvar%dtherm_cond_dT = 0.d0
  auxvar%heat_cap_liquid = 0.d0
  auxvar%heat_cap_solid = 0.d0
  auxvar%u = 0.d0
  auxvar%h = 0.d0
  auxvar%du_dP = 0.d0
  auxvar%du_dT = 0.d0
  auxvar%du_dC = 0.d0
  auxvar%dh_dP = 0.d0
  auxvar%dh_dT = 0.d0
  auxvar%dh_dC = 0.d0

end subroutine THCAuxVarInit

! ************************************************************************** !

subroutine THCAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Option_module

  implicit none

  type(thc_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pres = auxvar%pres
  auxvar2%sat = auxvar%sat
  auxvar2%pc = auxvar%pc
  auxvar2%kr = auxvar%kr
  auxvar2%effective_porosity = auxvar%effective_porosity
  auxvar2%dpor_dp = auxvar%dpor_dp
  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dkr_dp = auxvar%dkr_dp
  auxvar2%effective_saturation = auxvar%effective_saturation
  auxvar2%deffsat_dp = auxvar%deffsat_dp
  auxvar2%conc = auxvar%conc
  auxvar2%pert = auxvar%pert

  auxvar2%temp = auxvar%temp
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%den_kmol = auxvar%den_kmol
  auxvar2%dden_dp = auxvar%dden_dp
  auxvar2%dden_dT = auxvar%dden_dT
  auxvar2%dden_dC = auxvar%dden_dC
  auxvar2%vis = auxvar%vis
  auxvar2%dvis_dT = auxvar%dvis_dT
  auxvar2%dvis_dC = auxvar%dvis_dC
  auxvar2%diff_mol = auxvar%diff_mol
  auxvar2%ddiff_dT = auxvar%ddiff_dT
  auxvar2%ddiff_dC = auxvar%ddiff_dC
  auxvar2%therm_cond_eff = auxvar%therm_cond_eff
  auxvar2%dtherm_cond_dsat = auxvar%dtherm_cond_dsat
  auxvar2%dtherm_cond_dT = auxvar%dtherm_cond_dT
  auxvar2%heat_cap_liquid = auxvar%heat_cap_liquid
  auxvar2%heat_cap_solid = auxvar%heat_cap_solid
  auxvar2%u = auxvar%u
  auxvar2%h = auxvar%h
  auxvar2%du_dP = auxvar%du_dP
  auxvar2%du_dT = auxvar%du_dT
  auxvar2%du_dC = auxvar%du_dC
  auxvar2%dh_dP = auxvar%dh_dP
  auxvar2%dh_dT = auxvar%dh_dT
  auxvar2%dh_dC = auxvar%dh_dC

end subroutine THCAuxVarCopy

! ************************************************************************** !

subroutine THCAuxVarCompute(x,thc_auxvar,global_auxvar, &
                                material_auxvar,characteristic_curves, &
                                thc_parameter, &
                                natural_id,update_porosity,option)
  !
  ! Computes auxiliary variables for each grid cell.
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Option_module
  use Global_Aux_module
  use Characteristic_Curves_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(:)
  type(thc_auxvar_type) :: thc_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(thc_parameter_type) :: thc_parameter
  PetscBool :: update_porosity
  PetscInt :: natural_id

  PetscInt :: imat
  PetscReal :: arrhenius_scale
  PetscReal :: dummy_dist(-1:3)
  PetscBool :: saturated
  PetscReal :: dkr_dsat
  PetscReal :: deffsat_dsat
  PetscReal :: mobility, dmobility_dp, dmobility_dT, dmobility_dC
  PetscErrorCode :: ierr

  ! --------------------------------------------------------------------------
  ! Step 1: read primary variables (P,T,C) from solution vector
  ! --------------------------------------------------------------------------
  thc_auxvar%pres = x(thc_pressure_dof)
  thc_auxvar%temp = x(thc_temperature_dof)
  thc_auxvar%conc = x(thc_concentration_dof)

  ! --------------------------------------------------------------------------
  ! Step 2: effective porosity (soil compressibility)
  ! --------------------------------------------------------------------------
  ! %porosity should never be used. set to bogus value to catch misuse
  material_auxvar%porosity = -888.d0

  if (update_porosity .and. soil_compressibility_index > 0 .and. &
      thc_auxvar%pres > 0.d0) then
    call MaterialCompressSoil(material_auxvar,thc_auxvar%pres, &
                              thc_auxvar%effective_porosity, &
                              thc_auxvar%dpor_dp)
  else
    thc_auxvar%effective_porosity = material_auxvar%porosity_base
    thc_auxvar%dpor_dp = 0.d0
  endif
  material_auxvar%porosity = thc_auxvar%effective_porosity

  ! --------------------------------------------------------------------------
  ! Step 3: capillary pressure and saturation
  ! --------------------------------------------------------------------------
  ! For a very large negative liquid pressure (e.g. -1.d18), the capillary
  ! pressure can go near infinite, resulting in ds_dp being < 1.d-40 below
  ! and flipping the cell to saturated, when it is really far from saturated.
  ! The large negative liquid pressure is then passed to the EOS causing it
  ! to blow up.  Therefore, we truncate to the max capillary pressure here.
  thc_auxvar%pc = min(option%flow%reference_pressure - thc_auxvar%pres, &
                          characteristic_curves%saturation_function%pcmax)

  if (thc_auxvar%pc > 0.d0) then
    saturated = PETSC_FALSE
    call characteristic_curves%saturation_function% &
                               Saturation(thc_auxvar%pc, &
                                          thc_auxvar%sat, &
                                          thc_auxvar%dsat_dp,option)
    ! if ds_dp is 0, we consider the cell saturated.
    if (thc_auxvar%dsat_dp < 1.d-40) then
      saturated = PETSC_TRUE
    else
      ! ----------------------------------------------------------------------
      ! Step 4: relative permeability
      ! ----------------------------------------------------------------------
      dkr_dsat = 0.d0
      call characteristic_curves%liq_rel_perm_function% &
                       RelativePermeability(thc_auxvar%sat, &
                                            thc_auxvar%kr, &
                                            dkr_dsat,option)
      thc_auxvar%dkr_dp = thc_auxvar%dsat_dp * dkr_dsat
      if (thc_tensorial_rel_perm) then
        call characteristic_curves%liq_rel_perm_function% &
          EffectiveSaturation(thc_auxvar%sat, &
                              thc_auxvar%effective_saturation, &
                              deffsat_dsat,option)
        thc_auxvar%deffsat_dp = deffsat_dsat * thc_auxvar%dsat_dp
      endif
    endif
  else
    saturated = PETSC_TRUE
  endif

  ! the purpose for splitting this condition from the 'else' statement
  ! above is due to SaturationFunctionCompute switching a cell to
  ! saturated to prevent unstable (potentially infinite) derivatives when
  ! capillary pressure is very small
  if (saturated) then
    thc_auxvar%pc = 0.d0
    thc_auxvar%sat = 1.d0
    thc_auxvar%kr = 1.d0
    thc_auxvar%dsat_dp = 0.d0
    thc_auxvar%dkr_dp = 0.d0
    thc_auxvar%effective_saturation = 1.d0
    thc_auxvar%deffsat_dp = 0.d0
  endif

  if (option%iflag /= THC_UPDATE_FOR_DERIVATIVE) then
    global_auxvar%sat(1) = thc_auxvar%sat
    if (size(global_auxvar%sat) > 1) then
      global_auxvar%sat(2) = 1.d0 - global_auxvar%sat(1)
    endif
  endif

  ! --------------------------------------------------------------------------
  ! Step 5: liquid density rho_l(P,T,C) and its derivatives
  ! --------------------------------------------------------------------------
  ierr = 0
  call THCDensityAndDerivs(thc_auxvar%temp,thc_auxvar%pres, &
                               thc_auxvar%conc,thc_auxvar%den_kg, &
                               thc_auxvar%den_kmol,thc_auxvar%dden_dp, &
                               thc_auxvar%dden_dT,thc_auxvar%dden_dC, &
                               ierr)
  if (ierr /= 0) then
    option%io_buffer = 'Error computing THC liquid density in &
      &THCAuxVarCompute.'
    call PrintErrMsg(option)
  endif

  ! --------------------------------------------------------------------------
  ! Step 6: liquid viscosity mu_l(T,C) and its derivatives
  ! --------------------------------------------------------------------------
  call THCViscosityAndDerivs(thc_auxvar%temp,thc_auxvar%pres, &
                                 thc_auxvar%conc,thc_auxvar%den_kg, &
                                 thc_auxvar%vis, &
                                 thc_auxvar%dvis_dT,thc_auxvar%dvis_dC, &
                                 ierr)
  if (ierr /= 0) then
    option%io_buffer = 'Error computing THC liquid viscosity in &
      &THCAuxVarCompute.'
    call PrintErrMsg(option)
  endif

  ! --------------------------------------------------------------------------
  ! Step 7: effective porous-medium thermal conductivity
  ! --------------------------------------------------------------------------
  ! Per-material thermal properties are indexed by the cell's material id.
  ! Precedence: THERMAL_CHARACTERISTIC_CURVES (the only path giving kappa(T)),
  ! then the TH-style bulk interpolation kappa_eff = ckdry+sqrt(S)*(ckwet-ckdry)
  ! from THERMAL_CONDUCTIVITY_DRY/WET, then the grain-Somerton model from the
  ! OPTIONS global thc_kappa_solid.
  ! Output only -- the flux kernels re-evaluate per connection via
  ! THCThermalCondEval with the connection direction.  For an anisotropic curve
  ! no single scalar is meaningful, and because TCondTensorToScalar overwrites
  ! kT_dry/wet on the shared curve object this value reflects whichever
  ! connection was projected last.
  imat = material_auxvar%id
  dummy_dist = 0.d0
  call THCThermalCondEval(thc_auxvar,imat,thc_parameter,dummy_dist, &
                          PETSC_FALSE,thc_auxvar%therm_cond_eff, &
                          thc_auxvar%dtherm_cond_dsat, &
                          thc_auxvar%dtherm_cond_dT,option)

  ! --------------------------------------------------------------------------
  ! Step 8: derived quantities
  ! --------------------------------------------------------------------------
  ! volumetric heat capacities
  thc_auxvar%heat_cap_liquid = thc_auxvar%den_kg * &
                                  thc_specific_heat_liquid
  ! solid volumetric heat capacity rho_s*c_s [J/(m^3.K)], per-material when set
  if (associated(thc_parameter%dencpr)) then
    thc_auxvar%heat_cap_solid = thc_parameter%dencpr(imat)
  else
    thc_auxvar%heat_cap_solid = thc_density_solid * &
                                    thc_specific_heat_solid
  endif

  ! molecular diffusion [m^2/s] from FLUID_PROPERTY DIFFUSION_COEFFICIENT,
  ! optionally scaled by Arrhenius(AE,T,25C) when DIFFUSION_ACTIVATION_ENERGY
  ! is supplied (opt-in, as in RT/NWT).  No concentration dependence.
  thc_auxvar%diff_mol = thc_parameter%diffusion_coef
  thc_auxvar%ddiff_dT = 0.d0
  thc_auxvar%ddiff_dC = 0.d0
  if (Initialized(thc_parameter%diffusion_activation_energy)) then
    arrhenius_scale = Arrhenius(thc_parameter%diffusion_activation_energy, &
                                thc_auxvar%temp,thc_diffusion_ref_temp)
    thc_auxvar%diff_mol = thc_auxvar%diff_mol * arrhenius_scale
    thc_auxvar%ddiff_dT = thc_auxvar%diff_mol * &
      thc_parameter%diffusion_activation_energy / IDEAL_GAS_CONSTANT / &
      (thc_auxvar%temp + T273K)**2
  endif

  ! --------------------------------------------------------------------------
  ! Step 8b: full-EOS specific enthalpy / internal energy (FULL_EOS mode only)
  ! --------------------------------------------------------------------------
  ! Mirrors TH (th_aux.F90) but stays in mass-specific J/kg (no option%scale).
  ! Liquid enthalpy h(P,T) from EOS water; internal energy u = h - P/rho.
  if (thc_energy_mode == THC_ENERGY_FULL_EOS) then
    call THCEnergyEOS(thc_auxvar,option)
  else
    thc_auxvar%u = 0.d0
    thc_auxvar%h = 0.d0
    thc_auxvar%du_dP = 0.d0
    thc_auxvar%du_dT = 0.d0
    thc_auxvar%du_dC = 0.d0
    thc_auxvar%dh_dP = 0.d0
    thc_auxvar%dh_dT = 0.d0
    thc_auxvar%dh_dC = 0.d0
  endif

  ! derived mobility derivatives. These are recomputed
  ! by the flux kernel from the stored primitives (kr, dkr_dp, vis, dvis_dT,
  ! dvis_dC); evaluated here to document the chain rule.
  ! Assumes: no pressure dependence of viscosity mu_l, and
  !          no temperature or concentration dependence of kr. Then:
  !   mobility       = kr / mu_l
  !   d(mob)/dP_l    = dkr_dp / mu_l                      (through saturation)
  !   d(mob)/dT      = -kr / mu_l^2 * dmu_l/dT
  !   d(mob)/dC      = -kr / mu_l^2 * dmu_l/dC
  mobility     = thc_auxvar%kr / thc_auxvar%vis
  dmobility_dp = thc_auxvar%dkr_dp / thc_auxvar%vis
  dmobility_dT = -thc_auxvar%kr / thc_auxvar%vis**2 * &
                 thc_auxvar%dvis_dT
  dmobility_dC = -thc_auxvar%kr / thc_auxvar%vis**2 * &
                 thc_auxvar%dvis_dC

end subroutine THCAuxVarCompute

! ************************************************************************** !

subroutine THCThermalCondEval(thc_auxvar,imat,thc_parameter,dist,project, &
                              kappa_eff,dkappa_dsat,dkappa_dT,option)
  !
  ! Effective porous-medium thermal conductivity and its derivatives.
  !
  ! Precedence: THERMAL_CHARACTERISTIC_CURVES (the only path giving kappa(T) or
  ! anisotropy), then the TH-style bulk interpolation
  ! kappa_eff = ckdry+sqrt(S)*(ckwet-ckdry) from THERMAL_CONDUCTIVITY_DRY/WET,
  ! then the grain-Somerton model from the OPTIONS global thc_kappa_solid.
  !
  ! project = .true. applies TCondTensorToScalar(dist) first, which projects an
  ! anisotropic tensor onto the connection direction.  It overwrites kT_dry/wet
  ! on the shared curve object, so it must immediately precede CalculateTCond --
  ! the same ordering TH, GENERAL and SCO2 use.  Callers without a connection
  ! direction (per-cell output) pass .false.
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 08/14/26
  !
  use Option_module

  implicit none

  type(thc_auxvar_type) :: thc_auxvar
  PetscInt :: imat
  type(thc_parameter_type) :: thc_parameter
  PetscReal :: dist(-1:3)
  PetscBool :: project
  PetscReal :: kappa_eff, dkappa_dsat, dkappa_dT
  type(option_type) :: option

  PetscBool :: thermal_cc_set
  PetscReal :: kappa_dry_mat, kappa_wet_mat, sqrt_sat_mat

  dkappa_dT = 0.d0

  thermal_cc_set = PETSC_FALSE
  if (associated(thc_parameter%thermal_cc)) then
    thermal_cc_set = associated(thc_parameter%thermal_cc(imat)%ptr)
  endif
  kappa_dry_mat = UNINITIALIZED_DOUBLE
  kappa_wet_mat = UNINITIALIZED_DOUBLE
  if (associated(thc_parameter%ckdry)) then
    kappa_dry_mat = thc_parameter%ckdry(imat)
    kappa_wet_mat = thc_parameter%ckwet(imat)
  endif

  if (thermal_cc_set) then
    if (project) then
      call thc_parameter%thermal_cc(imat)%ptr% &
             thermal_conductivity_function%TCondTensorToScalar(dist,option)
    endif
    call thc_parameter%thermal_cc(imat)%ptr% &
           thermal_conductivity_function% &
           CalculateTCond(max(thc_auxvar%sat,thc_sat_floor), &
                          thc_auxvar%temp, &
                          thc_auxvar%effective_porosity, &
                          kappa_eff,dkappa_dsat,dkappa_dT,option)
  else if (Initialized(kappa_dry_mat) .and. Initialized(kappa_wet_mat)) then
    ! TH-style bulk dry/wet interpolation (saturation floored for finite deriv)
    sqrt_sat_mat = sqrt(max(thc_auxvar%sat,thc_sat_floor))
    kappa_eff = kappa_dry_mat + sqrt_sat_mat * (kappa_wet_mat - kappa_dry_mat)
    dkappa_dsat = (kappa_wet_mat - kappa_dry_mat) / (2.d0 * sqrt_sat_mat)
  else
    ! grain-Somerton fallback using the OPTIONS-block global solid conductivity
    call THCThermalConductivityEff(thc_auxvar%sat, &
                                       thc_auxvar%effective_porosity, &
                                       thc_kappa_solid, &
                                       kappa_eff,dkappa_dsat)
  endif

end subroutine THCThermalCondEval

! ************************************************************************** !

subroutine THCEnergyEOS(thc_auxvar,option)
  !
  ! Computes mass-specific liquid enthalpy h(P,T) and internal energy u(P,T)
  ! and their derivatives from the water EOS, for the FULL_EOS energy
  ! formulation (single-phase liquid, no ice/vapor).
  !
  ! TH (th_aux.F90) carries molar enthalpy scaled by option%scale (MJ/kmol).
  ! Here we keep mass-specific SI units: h[J/kg] = h_molar[J/kmol]/FMWH2O,
  ! and u = h - P/rho_kg (the P*v work term), consistent with TH's
  !   u_molar = h_molar - P/rho_molar  divided through by FMWH2O.
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/25/26

  use Option_module

  implicit none

  type(thc_auxvar_type) :: thc_auxvar
  type(option_type) :: option

  PetscReal :: hw, hw_dp, hw_dT       ! molar enthalpy [J/kmol] and derivs
  PetscReal :: den_kg, P, inv_den, inv_den2
  PetscErrorCode :: ierr

  ierr = 0
  P = thc_auxvar%pres
  den_kg = thc_auxvar%den_kg

  ! EOS water enthalpy: hw [J/kmol], hw_dp [J/(kmol.Pa)], hw_dT [J/(kmol.K)]
  call EOSWaterEnthalpy(thc_auxvar%temp,P,hw,hw_dp,hw_dT,ierr)
  if (ierr /= 0) then
    option%io_buffer = 'Error computing THC liquid enthalpy in &
      &THCEnergyEOS (FULL_EOS energy formulation).'
    call PrintErrMsg(option)
  endif

  ! Convert molar -> mass-specific [J/kg]
  thc_auxvar%h    = hw    / FMWH2O
  thc_auxvar%dh_dP = hw_dp / FMWH2O
  thc_auxvar%dh_dT = hw_dT / FMWH2O
  thc_auxvar%dh_dC = 0.d0   ! no solute enthalpy correction yet

  ! Internal energy u = h - P/rho_kg   [J/kg]
  !   du/dP = dh/dP - (1/rho - P/rho^2 * drho/dP)
  !   du/dT = dh/dT + P/rho^2 * drho/dT
  !   du/dC = dh/dC + P/rho^2 * drho/dC
  inv_den  = 1.d0 / den_kg
  inv_den2 = inv_den * inv_den
  thc_auxvar%u = thc_auxvar%h - P * inv_den
  thc_auxvar%du_dP = thc_auxvar%dh_dP - inv_den + &
                         P * inv_den2 * thc_auxvar%dden_dp
  thc_auxvar%du_dT = thc_auxvar%dh_dT + &
                         P * inv_den2 * thc_auxvar%dden_dT
  thc_auxvar%du_dC = thc_auxvar%dh_dC + &
                         P * inv_den2 * thc_auxvar%dden_dC

end subroutine THCEnergyEOS

! ************************************************************************** !

subroutine THCAuxVarPerturb(x,thc_auxvar,global_auxvar, &
                                material_auxvar, &
                                material_auxvar_pert, &
                                characteristic_curves, &
                                thc_parameter, &
                                natural_id, &
                                option)
  ! Calculates auxiliary variables for perturbed system
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  PetscReal :: x(:)
  type(option_type) :: option
  PetscInt :: natural_id
  type(thc_auxvar_type) :: thc_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(material_auxvar_type) :: material_auxvar_pert(:)
  class(characteristic_curves_type) :: characteristic_curves
  type(thc_parameter_type) :: thc_parameter

  PetscInt :: idof
  PetscReal :: x_pert(THC_NDOF), pert

  ! THC_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = THC_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    pert = x(idof)*thc_rel_pert+thc_min_pert(idof)
    thc_auxvar(idof)%pert = pert
    x_pert(1:option%nflowdof) = x
    x_pert(idof) = x(idof) + pert
    call THCAuxVarCompute(x_pert,thc_auxvar(idof),global_auxvar, &
                              material_auxvar, &
                              characteristic_curves, &
                              thc_parameter,natural_id, &
                              PETSC_TRUE,option)
  enddo

end subroutine THCAuxVarPerturb

! ************************************************************************** !

subroutine THCAuxTensorialRelPerm(auxvar,tensorial_rel_perm_exponent, &
                                      dist,rel_perm,drel_perm_dp,option)
  !
  ! Computes the tensorial (direction-dependent) relative permeability
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Option_module
  use Utility_module

  implicit none

  type(thc_auxvar_type) :: auxvar
  PetscReal :: tensorial_rel_perm_exponent(3)
  PetscReal :: dist(-1:3)
  PetscReal :: rel_perm
  PetscReal :: drel_perm_dp
  type(option_type) :: option

  PetscReal :: exponent_
  PetscReal :: tensorial_scale

  exponent_ = UtilityTensorToScalar(dist,tensorial_rel_perm_exponent)

  ! remember that the default 0.5 was subtracted from the tensorial value
  ! in THCSetup. If 0.5 is specified for the tensorial exponent in the
  ! input file, this value will be 0.
  tensorial_scale = auxvar%effective_saturation**exponent_
  rel_perm = auxvar%kr * tensorial_scale
  drel_perm_dp = auxvar%dkr_dp * tensorial_scale + &
                 exponent_ * rel_perm / &
                 auxvar%effective_saturation * &
                 auxvar%deffsat_dp

end subroutine THCAuxTensorialRelPerm

! ************************************************************************** !

subroutine THCPrintAuxVars(thc_auxvar,global_auxvar,material_auxvar, &
                               natural_id,string,option)
  !
  ! Prints out the contents of an auxvar
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(thc_auxvar_type) :: thc_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  type(option_type) :: option

  print *, '--------------------------------------------------------'
  print *, trim(string)
  print *, '                 cell id: ', natural_id
  print *, '         liquid pressure: ', thc_auxvar%pres
  print *, '             temperature: ', thc_auxvar%temp
  print *, '           concentration: ', thc_auxvar%conc
  print *, '      capillary pressure: ', thc_auxvar%pc
  print *, '       liquid saturation: ', thc_auxvar%sat
  print *, '      liquid sat (deriv): ', thc_auxvar%dsat_dp
  print *, '         liquid rel perm: ', thc_auxvar%kr
  print *, ' liquid rel perm (deriv): ', thc_auxvar%dkr_dp
  print *, '      effective porosity: ', thc_auxvar%effective_porosity
  print *, '   eff. porosity (deriv): ', thc_auxvar%dpor_dp
  print *, '          liquid density: ', thc_auxvar%den_kg
  print *, '        liquid viscosity: ', thc_auxvar%vis
  print *, '  effective therm. cond.: ', thc_auxvar%therm_cond_eff
  print *, '--------------------------------------------------------'

end subroutine THCPrintAuxVars

! ************************************************************************** !

subroutine THCOutputAuxVars1(thc_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,string,append,option)
  !
  ! Prints out the contents of an auxvar to a file
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(thc_auxvar_type) :: thc_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  PetscBool :: append
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string2

  write(string2,*) natural_id
  string2 = trim(adjustl(string)) // '_' // trim(adjustl(string2)) // '.txt'
  if (append) then
    open(unit=IUNIT_TEMP,file=string2,position='append')
  else
    open(unit=IUNIT_TEMP,file=string2)
  endif

  write(IUNIT_TEMP,*) '--------------------------------------------------------'
  write(IUNIT_TEMP,*) trim(string)
  write(IUNIT_TEMP,*) '                 cell id: ', natural_id
  write(IUNIT_TEMP,*) '         liquid pressure: ', thc_auxvar%pres
  write(IUNIT_TEMP,*) '             temperature: ', thc_auxvar%temp
  write(IUNIT_TEMP,*) '           concentration: ', thc_auxvar%conc
  write(IUNIT_TEMP,*) '      capillary pressure: ', thc_auxvar%pc
  write(IUNIT_TEMP,*) '       liquid saturation: ', thc_auxvar%sat
  write(IUNIT_TEMP,*) '      liquid sat (deriv): ', thc_auxvar%dsat_dp
  write(IUNIT_TEMP,*) '         liquid rel perm: ', thc_auxvar%kr
  write(IUNIT_TEMP,*) ' liquid rel perm (deriv): ', thc_auxvar%dkr_dp
  write(IUNIT_TEMP,*) '      effective porosity: ', thc_auxvar &
                                                      %effective_porosity
  write(IUNIT_TEMP,*) '   eff. porosity (deriv): ', thc_auxvar%dpor_dp
  write(IUNIT_TEMP,*) '          liquid density: ', thc_auxvar%den_kg
  write(IUNIT_TEMP,*) '        liquid viscosity: ', thc_auxvar%vis
  write(IUNIT_TEMP,*) '  effective therm. cond.: ', thc_auxvar &
                                                      %therm_cond_eff
  write(IUNIT_TEMP,*) '--------------------------------------------------------'

  close(IUNIT_TEMP)

end subroutine THCOutputAuxVars1

! ************************************************************************** !

function THCAuxMapConditionIndices(include_water_aux)
  !
  ! Maps indexing of conditions
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Option_module

  implicit none

  PetscBool :: include_water_aux

  PetscInt, pointer :: THCAuxMapConditionIndices(:)

  PetscInt, pointer :: mapping(:)

  PetscInt :: temp_int

  allocate(mapping(THC_MAX_INDEX))
  mapping = UNINITIALIZED_INTEGER

  ! THC always has all three equations active (fixed DOF), so every branch
  ! is taken unconditionally.
  temp_int = 0
  temp_int = temp_int + 1
  mapping(THC_COND_WATER_INDEX) = temp_int
  if (include_water_aux) then
    temp_int = temp_int + 1
    mapping(THC_COND_WATER_AUX_INDEX) = temp_int
  endif
  temp_int = temp_int + 1
  mapping(THC_COND_ENERGY_INDEX) = temp_int
  temp_int = temp_int + 1
  mapping(THC_COND_SOLUTE_INDEX) = temp_int

  THCAuxMapConditionIndices => mapping

end function THCAuxMapConditionIndices

! ************************************************************************** !

subroutine THCAuxVarSingleDestroy(auxvar)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  implicit none

  type(thc_auxvar_type), pointer :: auxvar

  if (associated(auxvar)) then
    call THCAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine THCAuxVarSingleDestroy

! ************************************************************************** !

subroutine THCAuxVarArray1Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  implicit none

  type(thc_auxvar_type), pointer :: auxvars(:)

  PetscInt :: iaux

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call THCAuxVarStrip(auxvars(iaux))
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine THCAuxVarArray1Destroy

! ************************************************************************** !

subroutine THCAuxVarArray2Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  implicit none

  type(thc_auxvar_type), pointer :: auxvars(:,:)

  PetscInt :: iaux, idof

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call THCAuxVarStrip(auxvars(idof-1,iaux))
      enddo
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine THCAuxVarArray2Destroy

! ************************************************************************** !

subroutine THCMaterialAuxVarDestroy(auxvars)
  !
  ! Deallocates material auxiliary object for thc perturbation
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  implicit none

  type(material_auxvar_type), pointer :: auxvars(:,:)

  PetscInt :: iaux, idof

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call MaterialAuxVarStrip(auxvars(idof,iaux))
      enddo
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine THCMaterialAuxVarDestroy

! ************************************************************************** !

subroutine THCAuxVarStrip(auxvar)
  !
  ! THCAuxVarDestroy: Deallocates a thc auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(thc_auxvar_type) :: auxvar

end subroutine THCAuxVarStrip

! ************************************************************************** !

subroutine THCAuxDestroy(aux)
  !
  ! Deallocates a general auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(thc_type), pointer :: aux

  call DeallocateArray(thc_min_pert)

  if (.not.associated(aux)) return

  call THCAuxVarDestroy(aux%auxvars)
  call THCAuxVarDestroy(aux%auxvars_bc)
  call THCAuxVarDestroy(aux%auxvars_ss)
  call THCMaterialAuxVarDestroy(aux%material_auxvars_pert)

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%thc_parameter)) then
    call DeallocateArray(aux%thc_parameter%tensorial_rel_perm_exponent)
    call DeallocateArray(aux%thc_parameter%dencpr)
    call DeallocateArray(aux%thc_parameter%ckdry)
    call DeallocateArray(aux%thc_parameter%ckwet)
    if (associated(aux%thc_parameter%thermal_cc)) then
      deallocate(aux%thc_parameter%thermal_cc)
      nullify(aux%thc_parameter%thermal_cc)
    endif
    deallocate(aux%thc_parameter)
  endif
  nullify(aux%thc_parameter)

  deallocate(aux)
  nullify(aux)

end subroutine THCAuxDestroy

end module THC_Aux_module
