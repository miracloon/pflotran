module TH_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module
  use Geomechanics_Linear_Aux_module

  implicit none

  private

  PetscInt, public :: th_ni_count
  PetscInt, public :: th_ts_cut_count
  PetscInt, public :: th_ts_count

  PetscBool, public :: th_scale_by_volume = PETSC_TRUE
  PetscBool, public :: th_numerical_derivatives = PETSC_FALSE
  PetscReal, public :: th_rel_pert = 1.d-8
  PetscReal, public :: th_pres_min_pert = 1.d-2
  PetscReal, public :: th_temp_min_pert = 1.d-6
  PetscReal, pointer, public :: th_min_pert(:)

  PetscInt, public :: th_ice_model

  PetscInt, parameter, public :: TH_LIQUID_EQUATION_INDEX = 1
  PetscInt, parameter, public :: TH_ENERGY_EQUATION_INDEX = 2
  PetscInt, public :: th_well_eq = UNINITIALIZED_INTEGER

  PetscInt, parameter, public :: TH_PRESSURE_DOF = 1
  PetscInt, parameter, public :: TH_TEMPERATURE_DOF = 2
  PetscInt, public :: th_well_dof = UNINITIALIZED_INTEGER
  PetscInt, public :: th_conductance_dof = UNINITIALIZED_INTEGER

  PetscReal, public :: th_well_init_well_temperature = UNINITIALIZED_DOUBLE

  PetscInt, parameter, public :: TH_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: TH_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: TH_UPDATE_FOR_BOUNDARY = 2
  PetscInt, parameter, public :: TH_UPDATE_FOR_PERTURBATION = 3

  type, public :: th_auxvar_type
    PetscReal :: pres
    PetscReal :: temp
    PetscReal :: sat
    PetscReal :: den
    PetscReal :: den_kg
    PetscReal :: avgmw
    PetscReal :: h
    PetscReal :: u
    PetscReal :: pc
    PetscReal :: vis
    PetscReal :: kr
    PetscReal :: mobility
    PetscReal :: effective_porosity
    PetscReal :: dsat_dp
    PetscReal :: dsat_dT
    PetscReal :: dden_dp
    PetscReal :: dden_dT
    PetscReal :: dmobility_dp
    PetscReal :: dmobility_dT
    PetscReal :: dh_dp
    PetscReal :: dh_dT
    PetscReal :: du_dp
    PetscReal :: du_dT
    PetscReal :: Dk_eff
    PetscReal :: Ke
    PetscReal :: dKe_dp
    PetscReal :: dKe_dT
    PetscReal :: dpor_dp
    PetscReal :: dpor_dT
    PetscReal :: dpres_dtime ! for TS
    PetscReal :: dtemp_dtime  ! for TS
    PetscReal :: d2sat_dp2
    PetscReal :: d2den_dp2
    PetscReal :: d2u_dp2
    PetscReal :: d2sat_dT2
    PetscReal :: d2den_dT2
    PetscReal :: d2u_dT2
    PetscReal :: d2sat_dTdp
    PetscReal :: d2den_dTdp
    PetscReal :: d2u_dTdp

    ! for ice
    type(th_ice_type), pointer :: ice

    ! for well model
    PetscInt :: iwellaux

    ! for numerical derivatives
    PetscReal :: pert
    type(th_auxvar_type), pointer :: auxvar_pert(:)

#if defined(CLM_PFLOTRAN) || defined(CLM_OFFLINE)
    PetscReal :: bc_alpha  ! Brooks Corey parameterization: alpha
    PetscReal :: bc_lambda ! Brooks Corey parameterization: lambda
#endif
  end type th_auxvar_type

  type, public :: th_ice_type
    PetscReal :: Ke_fr
    PetscReal :: dKe_fr_dp
    PetscReal :: dKe_fr_dT
    ! ice
    PetscReal :: sat_ice
    PetscReal :: sat_gas
    PetscReal :: dsat_ice_dp
    PetscReal :: dsat_gas_dp
    PetscReal :: dsat_ice_dT
    PetscReal :: dsat_gas_dT
    PetscReal :: den_ice
    PetscReal :: dden_ice_dp
    PetscReal :: dden_ice_dT
    PetscReal :: u_ice
    PetscReal :: du_ice_dT
    PetscReal :: du_ice_dP
    PetscReal :: den_gas
    PetscReal :: dden_gas_dT
    PetscReal :: u_gas
    PetscReal :: du_gas_dT
    PetscReal :: mol_gas
    PetscReal :: dmol_gas_dT
    ! For DallAmico model
    PetscReal :: pres_fh2o
    PetscReal :: dpres_fh2o_dp
    PetscReal :: dpres_fh2o_dT
  end type th_ice_type

  type, public :: th_well_auxvar_type
    PetscInt :: ghosted_id
    PetscReal :: z
    PetscReal :: temp
    PetscReal :: u
    PetscReal :: h
    PetscReal :: du_dT
    PetscReal :: dh_dT
    PetscReal :: volume
    PetscReal :: segment_length
    PetscReal :: heat_transfer_coef
    PetscReal :: well_index
    PetscReal :: therm_cond_borehole_to_cell
    ! Fluid properties for convective heat transfer
    PetscReal :: den                 ! density [kmol/m^3]
    PetscReal :: den_kg              ! density [kg/m^3]
    PetscReal :: vis                 ! viscosity [Pa-s]
    PetscReal :: therm_cond_fluid    ! thermal conductivity of fluid [W/m/K]
    PetscReal :: spec_heat_fluid     ! specific heat capacity of fluid [J/kg/K]
    ! derivatives
    PetscReal :: dden_dT
    PetscReal :: dtherm_cond_fluid_dT
    ! Dimensionless numbers
    PetscReal :: reynolds            ! Reynolds number [-]
    PetscReal :: prandtl             ! Prandtl number [-]
    PetscReal :: nusselt             ! Nusselt number [-]
    ! for numerical derivatives
    PetscReal :: pert
    type(th_well_auxvar_type), pointer :: auxvar_pert(:)
  end type th_well_auxvar_type

  type, public :: th_parameter_type
    PetscReal, pointer :: dencpr(:)
    PetscReal, pointer :: ckdry(:) ! Thermal conductivity (dry)
    PetscReal, pointer :: ckwet(:) ! Thermal conductivity (wet)
    PetscReal, pointer :: alpha(:)
    PetscReal, pointer :: ckfrozen(:) ! Thermal conductivity (frozen soil)
    PetscReal, pointer :: alpha_fr(:) ! exponent frozen
    PetscReal, pointer :: sir(:,:)
    type(geomech_linear_parameter_type), pointer :: geomech_parameter
  end type th_parameter_type

  type, public :: TH_type
    PetscBool :: auxvars_up_to_date
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(th_parameter_type), pointer :: th_parameter
    type(th_auxvar_type), pointer :: auxvars(:)
    type(th_auxvar_type), pointer :: auxvars_bc(:)
    type(th_auxvar_type), pointer :: auxvars_ss(:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
    PetscInt :: num_well_aux, num_well_aux_bc
    type(th_well_auxvar_type), pointer :: auxvars_well(:)
    type(th_well_auxvar_type), pointer :: auxvars_well_bc(:)
  end type TH_type

  PetscReal, parameter :: epsilon = 1.d-6
  PetscReal, parameter :: perturbation_tolerance = 1.d-8

  public :: THAuxCreate, &
            THAuxDestroy, &
            THAuxVarComputeNoFreezing, &
            THAuxVarInit, &
            THAuxVarCopy, &
            THAuxVarCompute2ndOrderDeriv, &
            THAuxVarPerturb, &
            THWellAuxVarInit, &
            THWellAuxVarCompute, &
            THWellAuxVarPerturb, &
            THWellAuxVarCopyParamsToPert, &
            THAuxVarDestroy

  public :: THAuxVarComputeFreezing

contains

! ************************************************************************** !

function THAuxCreate(option)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: ???
  ! Date: 02/14/08
  !

  use Option_module

  implicit none

  type(option_type) :: option
  type(TH_type), pointer :: THAuxCreate

  type(TH_type), pointer :: aux

  nullify(th_min_pert)

  ! initialize module variables
  th_ice_model = 0

  allocate(aux)
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  aux%num_well_aux = 0
  aux%num_well_aux_bc = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  nullify(aux%matrix_zeroing)
  nullify(aux%auxvars_well)
  nullify(aux%auxvars_well_bc)

  allocate(aux%th_parameter)
  nullify(aux%th_parameter%dencpr)
  nullify(aux%th_parameter%ckdry)
  nullify(aux%th_parameter%ckwet)
  nullify(aux%th_parameter%alpha)
  nullify(aux%th_parameter%ckfrozen)
  nullify(aux%th_parameter%alpha_fr)
  nullify(aux%th_parameter%sir)

  nullify(aux%th_parameter%geomech_parameter)

  THAuxCreate => aux

end function THAuxCreate

! ************************************************************************** !

recursive subroutine THAuxVarInit(auxvar,option,allocate_perturbation)
  !
  ! Initialize auxiliary object
  !
  ! Author: ???
  ! Date: 02/14/08
  !

  use Option_module
  use PFLOTRAN_Constants_module, only : UNINITIALIZED_DOUBLE

  implicit none

  type(th_auxvar_type) :: auxvar
  type(option_type) :: option
  PetscBool :: allocate_perturbation

  PetscReal :: uninit_value
  PetscInt :: i
!  uninit_value     = UNINITIALIZED_DOUBLE
  uninit_value     = 0.d0

  auxvar%pres      = uninit_value
  auxvar%temp      = uninit_value
  auxvar%sat       = uninit_value
  auxvar%den       = uninit_value
  auxvar%den_kg    = uninit_value
  auxvar%avgmw     = uninit_value
  auxvar%h         = uninit_value
  auxvar%u         = uninit_value
  auxvar%pc        = uninit_value
  auxvar%kr        = uninit_value
  !auxvar%dkr_dp   = uninit_value
  auxvar%vis       = uninit_value
  !auxvar%dvis_dp  = uninit_value
  auxvar%mobility  = uninit_value
  auxvar%effective_porosity  = uninit_value
  auxvar%dsat_dp   = uninit_value
  auxvar%dsat_dT   = uninit_value
  auxvar%dden_dp   = uninit_value
  auxvar%dden_dT   = uninit_value
  auxvar%dmobility_dp   = uninit_value
  auxvar%dmobility_dT   = uninit_value
  auxvar%dh_dp     = uninit_value
  auxvar%dh_dT     = uninit_value
  auxvar%du_dp     = uninit_value
  auxvar%du_dT     = uninit_value
  auxvar%Dk_eff    = uninit_value
  auxvar%Ke        = uninit_value
  auxvar%dKe_dp    = uninit_value
  auxvar%dKe_dT    = uninit_value
  auxvar%dpor_dp   = uninit_value
  auxvar%dpor_dT   = uninit_value
 if (option%flow%th_freezing) then
    allocate(auxvar%ice)
    auxvar%ice%Ke_fr     = uninit_value
    auxvar%ice%dKe_fr_dp = uninit_value
    auxvar%ice%dKe_fr_dT = uninit_value
    ! NOTE(bja, 2013-12) always initialize ice variables to zero, even if
    !                    not used!
    auxvar%ice%sat_ice       = uninit_value
    auxvar%ice%sat_gas       = uninit_value
    auxvar%ice%dsat_ice_dp   = uninit_value
    auxvar%ice%dsat_gas_dp   = uninit_value
    auxvar%ice%dsat_ice_dT   = uninit_value
    auxvar%ice%dsat_gas_dT   = uninit_value
    auxvar%ice%den_ice       = uninit_value
    auxvar%ice%dden_ice_dp   = uninit_value
    auxvar%ice%dden_ice_dT   = uninit_value
    auxvar%ice%u_ice         = uninit_value
    auxvar%ice%du_ice_dT     = uninit_value
    auxvar%ice%du_ice_dP     = uninit_value
    auxvar%ice%den_gas       = uninit_value
    auxvar%ice%dden_gas_dT   = uninit_value
    auxvar%ice%u_gas         = uninit_value
    auxvar%ice%du_gas_dT     = uninit_value
    auxvar%ice%mol_gas       = uninit_value
    auxvar%ice%dmol_gas_dT   = uninit_value
    auxvar%ice%pres_fh2o     = uninit_value
    auxvar%ice%dpres_fh2o_dp = uninit_value
    auxvar%ice%dpres_fh2o_dT = uninit_value
  else
    nullify(auxvar%ice)
  endif

  auxvar%iwellaux = UNINITIALIZED_INTEGER

#if defined(CLM_PFLOTRAN) || defined(CLM_OFFLINE)
  auxvar%bc_alpha      = uninit_value
  auxvar%bc_lambda     = uninit_value
#endif

  auxvar%pert = uninit_value
  if (allocate_perturbation) then
    allocate(auxvar%auxvar_pert(option%nflowdof))
    do i = 1, option%nflowdof
      call THAuxVarInit(auxvar%auxvar_pert(i),option,PETSC_FALSE)
    enddo
  else
    nullify(auxvar%auxvar_pert)
  endif

end subroutine THAuxVarInit

! ************************************************************************** !

subroutine THAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: ???
  ! Date: 12/13/07
  !

  use Option_module

  implicit none

  type(th_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pres = auxvar%pres
  auxvar2%temp = auxvar%temp
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%den_kg = auxvar%den_kg

  auxvar2%avgmw = auxvar%avgmw
  auxvar2%h = auxvar%h
  auxvar2%u = auxvar%u
  auxvar2%pc = auxvar%pc
  auxvar2%kr = auxvar%kr
!  auxvar2%dkr_dp = auxvar%dkr_dp
  auxvar2%vis = auxvar%vis
!  auxvar2%dvis_dp = auxvar%dvis_dp
  auxvar2%mobility = auxvar%mobility
  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dsat_dT = auxvar%dsat_dT
  auxvar2%dden_dp = auxvar%dden_dp
  auxvar2%dden_dT = auxvar%dden_dT
  auxvar2%dmobility_dp = auxvar%dmobility_dp
  auxvar2%dmobility_dT = auxvar%dmobility_dT
  auxvar2%dh_dp = auxvar%dh_dp
  auxvar2%dh_dT = auxvar%dh_dT
  auxvar2%du_dp = auxvar%du_dp
  auxvar2%du_dT = auxvar%du_dT
  auxvar2%Dk_eff = auxvar%Dk_eff
  auxvar2%Ke = auxvar%Ke
  auxvar2%dKe_dp = auxvar%dKe_dp
  auxvar2%dKe_dT = auxvar%dKe_dT
  if (associated(auxvar%ice)) then
    auxvar2%ice%Ke_fr = auxvar%ice%Ke_fr
    auxvar2%ice%dKe_fr_dp = auxvar%ice%dKe_fr_dp
    auxvar2%ice%dKe_fr_dT = auxvar%ice%dKe_fr_dT
    auxvar2%ice%sat_ice = auxvar%ice%sat_ice
    auxvar2%ice%sat_gas = auxvar%ice%sat_gas
    auxvar2%ice%dsat_ice_dp = auxvar%ice%dsat_ice_dp
    auxvar2%ice%dsat_gas_dp = auxvar%ice%dsat_gas_dp
    auxvar2%ice%dsat_ice_dT = auxvar%ice%dsat_ice_dT
    auxvar2%ice%dsat_gas_dT = auxvar%ice%dsat_gas_dT
    auxvar2%ice%den_ice = auxvar%ice%den_ice
    auxvar2%ice%dden_ice_dp = auxvar%ice%dden_ice_dp
    auxvar2%ice%dden_ice_dT = auxvar%ice%dden_ice_dT
    auxvar2%ice%u_ice = auxvar%ice%u_ice
    auxvar2%ice%du_ice_dT = auxvar%ice%du_ice_dT
    auxvar2%ice%du_ice_dP = auxvar%ice%du_ice_dP
    auxvar2%ice%pres_fh2o = auxvar%ice%pres_fh2o
    auxvar2%ice%dpres_fh2o_dp = auxvar%ice%dpres_fh2o_dp
    auxvar2%ice%dpres_fh2o_dT = auxvar%ice%dpres_fh2o_dT
    auxvar2%ice%den_gas = auxvar%ice%den_gas
    auxvar2%ice%dden_gas_dT = auxvar%ice%dden_gas_dT
    auxvar2%ice%u_gas = auxvar%ice%u_gas
    auxvar2%ice%du_gas_dT = auxvar%ice%du_gas_dT
    auxvar2%ice%mol_gas = auxvar%ice%mol_gas
    auxvar2%ice%dmol_gas_dT = auxvar%ice%dmol_gas_dT
  endif

#if defined(CLM_PFLOTRAN) || defined(CLM_OFFLINE)
  auxvar2%bc_alpha  = auxvar%bc_alpha
  auxvar2%bc_lambda = auxvar%bc_lambda
#endif

end subroutine THAuxVarCopy

! ************************************************************************** !

subroutine THAuxPorosity(th_auxvar,material_auxvar,global_auxvar, &
                         update_porosity,th_parameter,option)
  !
  ! Calculates the update to porosity for TH
  !
  ! Author: Glenn Hammond
  ! Date: 03/04/26
  !
  use Material_Aux_module
  use Global_Aux_module
  use Option_module
  use Geomechanics_Linear_Aux_module

  implicit none

  type(th_auxvar_type) :: th_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscBool :: update_porosity
  type(th_parameter_type) :: th_parameter
  type(option_type) :: option

  type(geomech_linear_parameter_type), pointer :: geomech_param

  PetscReal :: cell_pressure

  PetscReal :: biot_coeff
  PetscReal :: youngs_mod, poissons_ratio
  PetscReal :: dr_bulk_modulus
  PetscReal :: por

  PetscInt :: id_press_0, id_temp_0
  PetscInt :: id_vstrain_0, id_vstrain, id_press_mech, id_temp_mech
  PetscInt :: id_porosity_mech, id_porosity_flow
  PetscInt :: mat_id
  PetscReal :: porosity_0, C1, C2, press_0, temp_0
  PetscReal :: vstrain_0, vstrain, del_vstrain, alpha
  PetscReal :: press_mech, temp_mech

  geomech_param => th_parameter%geomech_parameter

  if (option%iflag /= TH_UPDATE_FOR_BOUNDARY) then
    if (update_porosity) then
      if (soil_compressibility_index > 0) then
        cell_pressure = th_auxvar%pres
        call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                  th_auxvar%effective_porosity, &
                                  th_auxvar%dpor_dp)
      else
        th_auxvar%effective_porosity = material_auxvar%porosity_base
        th_auxvar%dpor_dp =  0.d0
        th_auxvar%dpor_dT = 0.d0
      endif
    else
      th_auxvar%effective_porosity = material_auxvar%porosity_base
      th_auxvar%dpor_dp = 0.d0
      th_auxvar%dpor_dT = 0.d0
    endif

    if (associated(geomech_param)) then
      select case(option%geomechanics%flow_coupling)
        case(GEOMECH_TWO_WAY_COUPLED)
          select case(option%geomechanics%split_scheme)
            case(GEOMECH_FIXED_STRESS_SPLIT)
              ! geomechanics mat id corresponding to cell
              mat_id = material_auxvar%secondary_material_id
              ! get geomech material properties
              youngs_mod = geomech_param%youngs_modulus(mat_id)
              poissons_ratio = geomech_param%poissons_ratio(mat_id)
              biot_coeff = geomech_param%biot_coeff(mat_id)
              alpha = geomech_param%thermal_exp_coeff(mat_id)
              porosity_0 = material_auxvar%porosity_0
              ! 3D bulk modulus
              dr_bulk_modulus = youngs_mod / &
                                (3.d0 * (1.d0 - (2.d0 * poissons_ratio)))
              ! C1 constant
              C1 = (biot_coeff-porosity_0)*(1.d0 - biot_coeff)/dr_bulk_modulus
              ! C2 constant
              C2 = biot_coeff**2/dr_bulk_modulus
              ! get stored values
              id_press_0 = geomech_param%press_0_id
              id_temp_0 = geomech_param%temp_0_id
              id_vstrain_0 = geomech_param%vol_strain_0_id
              id_vstrain = geomech_param%vol_strain_id
              id_press_mech = geomech_param%stored_pressure_id
              id_temp_mech = geomech_param%stored_temperature_id
              id_porosity_mech = geomech_param%stored_porosity_id
              id_porosity_flow = geomech_param%flow_porosity_id
              alpha = geomech_param%thermal_exp_coeff(mat_id)

              press_0 = global_auxvar%parameters(id_press_0)
              press_mech = global_auxvar%parameters(id_press_mech)
              temp_0 = global_auxvar%parameters(id_temp_0)
              temp_mech = global_auxvar%parameters(id_temp_mech)
              vstrain_0 = global_auxvar%parameters(id_vstrain_0)
              vstrain = global_auxvar%parameters(id_vstrain)
              if (option%iflag == TH_UPDATE_FOR_FIXED_ACCUM) then
                por = global_auxvar%parameters(id_porosity_mech)
                th_auxvar%effective_porosity = por
              elseif (option%iflag == TH_UPDATE_FOR_ACCUM) then
                ! delta vstrain
                del_vstrain = vstrain - vstrain_0
                ! Burghardt 2017 paper
                por = porosity_0 + (biot_coeff * del_vstrain) + &
                      ( C1 * ( press_mech - press_0 )) + &
                      ( (C2 + C1) * ( th_auxvar%pres - press_mech )) + &
                      ! thermal term from continuous strain-based law
                      ( -3.d0*alpha*(biot_coeff - porosity_0) * &
                        ( temp_mech - temp_0 )) + &
                      ! thermal term for fixed-stress increment
                      ( 3.d0*alpha*porosity_0 * &
                        ( th_auxvar%temp - temp_mech ))
                ! store new (mass conserved) flow porosity for threshold check
                global_auxvar%parameters(id_porosity_flow) = por

                th_auxvar%effective_porosity = por
                ! Burghardt 2017 paper
                th_auxvar%dpor_dp = (biot_coeff**2)/dr_bulk_modulus + &
                    ((biot_coeff-porosity_0)*(1.d0-biot_coeff)) / &
                    dr_bulk_modulus
                th_auxvar%dpor_dT = 3.d0*alpha*porosity_0
              endif
          end select
      end select
    else
      material_auxvar%dporosity_dp = UNINITIALIZED_DOUBLE
    endif

    if (option%iflag /= TH_UPDATE_FOR_PERTURBATION) then
      material_auxvar%porosity = th_auxvar%effective_porosity
    endif
  endif

end subroutine THAuxPorosity

! ************************************************************************** !

subroutine THAuxVarComputeNoFreezing(x,auxvar,global_auxvar, &
                                     material_auxvar, &
                                     iphase,characteristic_curves, &
                                     thermal_cc, &
                                     th_parameter, icct, natural_id, &
                                     update_porosity,option)
  !
  ! Computes auxiliary variables for each grid cell
  !
  ! Author: ???
  ! Date: 02/22/08
  !

  use Option_module
  use Global_Aux_module

  use EOS_Water_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module
  use Characteristic_Curves_Thermal_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  class(cc_thermal_type) :: thermal_cc
  PetscReal :: x(option%nflowdof)
  type(th_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscInt :: iphase
  type(th_parameter_type) :: th_parameter
  PetscInt :: icct
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  PetscBool :: update_porosity

  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,hw,sat_pressure,visl
  PetscReal :: ds_dp, dkr_dp
  PetscReal :: dvis_dT, dvis_dp
  PetscReal :: dw_dp, dw_dT, hw_dp, hw_dT
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dT
  PetscReal :: Ke
  PetscReal :: alpha
  PetscReal :: Dk
  PetscReal :: Dk_dry
  PetscReal :: aux(1)
  PetscReal :: dkr_dsat1
  PetscReal :: dk_ds, dk_dT

  ierr = 0
  auxvar%den = 0.d0
  auxvar%den_kg = 0.d0
  auxvar%sat = 0.d0
  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%avgmw = 0.d0
  auxvar%mobility = 0.d0
  auxvar%kr = 0.d0

  auxvar%pres = x(1)
  auxvar%temp = x(2)

  global_auxvar%pres(1) = UNINITIALIZED_DOUBLE
  global_auxvar%temp = UNINITIALIZED_DOUBLE
  global_auxvar%den = UNINITIALIZED_DOUBLE
  global_auxvar%den_kg = UNINITIALIZED_DOUBLE
  global_auxvar%sat = UNINITIALIZED_DOUBLE

  call THAuxPorosity(auxvar,material_auxvar,global_auxvar,update_porosity, &
                     th_parameter,option)

  auxvar%pc = min(option%flow%reference_pressure - auxvar%pres, &
                  characteristic_curves%saturation_function%pcmax)

  if (option%flow%disable_capillarity) then
    auxvar%pc = 0.d0
  endif

!***************  Liquid phase properties **************************
  auxvar%avgmw = FMWH2O

  pw = option%flow%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0

!  if (auxvar%pc > 0.d0) then
  if (auxvar%pc > 1.d0) then
     iphase = 3

#if defined(CLM_PFLOTRAN) || defined(CLM_OFFLINE)
    if (auxvar%bc_alpha > 0.d0) then
      select type(sf => characteristic_curves%saturation_function)
        class is(sat_func_vg_type)
          sf%m     = auxvar%bc_lambda
          sf%alpha = auxvar%bc_alpha
        class is(sat_func_bc_type)
            sf%lambda = auxvar%bc_lambda
            sf%alpha  = auxvar%bc_alpha
        class default
          option%io_buffer = 'CLM-PFLOTRAN only supports ' // &
            'sat_func_vg_type and sat_func_bc_type'
          call printErrMsg(option)
      end select

      select type(rpf => characteristic_curves%liq_rel_perm_function)
        class is(rpf_mualem_vg_liq_type)
          rpf%m = auxvar%bc_lambda
        class is(rpf_burdine_bc_liq_type)
          rpf%lambda = auxvar%bc_lambda
        class is(rpf_mualem_bc_liq_type)
          rpf%lambda = auxvar%bc_lambda
        class is(rpf_burdine_vg_liq_type)
          rpf%m = auxvar%bc_lambda
        class default
          option%io_buffer = 'Unsupported LIQUID-REL-PERM-FUNCTION'
          call printErrMsg(option)
      end select
    endif
#endif

    call characteristic_curves%saturation_function% &
        Saturation(auxvar%pc,auxvar%sat, &
                   ds_dp, option)

    if (ds_dp < 1.d-40) then
      iphase = 1
      auxvar%pc = 0.d0
      auxvar%sat = 1.d0
      auxvar%kr = 1.d0
      pw = max(auxvar%pres,pw)
      dpw_dp = 1.d0
    else
      call characteristic_curves%liq_rel_perm_function% &
             RelativePermeability(auxvar%sat,auxvar%kr, &
                                  dkr_dsat1,option)

      dkr_dp = ds_dp * dkr_dsat1
      dpw_dp = 0.d0
    endif

  else
    iphase = 1
    auxvar%pc = 0.d0
    auxvar%sat = 1.d0
    auxvar%kr = 1.d0
    pw = max(auxvar%pres,pw)
    dpw_dp = 1.d0
  endif

  ! may need to compute dpsat_dT to pass to VISW
  call EOSWaterSaturationPressure(auxvar%temp,sat_pressure,dpsat_dT,ierr)
  call EOSWaterEnthalpy(auxvar%temp,pw,hw,hw_dp,hw_dT,ierr)
  if (.not.option%flow%density_depends_on_salinity) then
    call EOSWaterDensity(auxvar%temp,pw,dw_kg,dw_mol,dw_dp,dw_dT,ierr)
    if (ierr /= 0) then
      call PrintMsgByCell(option,natural_id, &
                       'Error in THAuxVarComputeNoFreezing->EOSWaterDensity')
    endif
    call EOSWaterViscosity(auxvar%temp,pw,sat_pressure,dpsat_dT,visl, &
                           dvis_dT,dvis_dp,ierr)
  else
    if (option%iflag == TH_UPDATE_FOR_FIXED_ACCUM) then
      ! For the computation of fixed accumulation term use NaCl
      ! value, m_nacl(2), from the previous time step.
      aux(1) = global_auxvar%m_nacl(2)
    else
      ! Use NaCl value for the current time step, m_nacl(1), for computing
      ! the accumulation term
      aux(1) = global_auxvar%m_nacl(1)
    endif
    call EOSWaterDensityExt(auxvar%temp,pw,aux, &
                            dw_kg,dw_mol,dw_dp,dw_dT,ierr)
    if (ierr /= 0) then
      call PrintMsgByCell(option,natural_id, &
                     'Error in THAuxVarComputeNoFreezing->EOSWaterDensityExt')
    endif
    call EOSWaterViscosityExt(auxvar%temp,pw,sat_pressure,dpsat_dT,aux, &
                              visl,dvis_dT,dvis_dp,ierr)
  endif
  ! J/kmol -> whatever units (default of option%scale is 1.e-6)
  hw = hw * option%scale
  hw_dp = hw_dp * option%scale
  hw_dT = hw_dT * option%scale

!  call VISW_noderiv(option%temp,pw,sat_pressure,visl,ierr)
  if (iphase == 3) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif

  auxvar%den = dw_mol
  auxvar%den_kg = dw_kg

  auxvar%h = hw
  auxvar%u = auxvar%h - pw / dw_mol * option%scale
  auxvar%mobility = auxvar%kr/visl

  auxvar%vis = visl
!  auxvar%dvis_dp = dvis_dp
!  auxvar%kr = kr
!  auxvar%dkr_dp = dkr_dp
  auxvar%dsat_dp = ds_dp
  auxvar%dden_dT = dw_dT
  auxvar%dsat_dT = 0.d0
  auxvar%dden_dp = dw_dp

!geh: contribution of dvis_dpsat is now added in EOSWaterViscosity
!  auxvar%dmobility_dT = -kr/(visl*visl)*(dvis_dT+dvis_dpsat*dpsat_dT)
  auxvar%dmobility_dT = -auxvar%kr/(visl*visl)*dvis_dT
  auxvar%dmobility_dp = dkr_dp/visl - auxvar%kr/(visl*visl)*dvis_dp
  if (iphase < 3) then !kludge since pw is constant in the unsat zone
    auxvar%dh_dp = hw_dp
    auxvar%du_dp = hw_dp - (dpw_dp/dw_mol-pw/(dw_mol*dw_mol)*dw_dp)*option%scale
  else
    auxvar%dh_dp = 0.d0
    auxvar%du_dp = 0.d0
  endif

  auxvar%dh_dT = hw_dT
  auxvar%du_dT = hw_dT + pw/(dw_mol*dw_mol)*option%scale*dw_dT

  ! Parameters for computation of effective thermal conductivity
  alpha = th_parameter%alpha(icct)
  Dk = th_parameter%ckwet(icct)
  Dk_dry = th_parameter%ckdry(icct)

  !unfrozen soil Kersten number
  Ke = (auxvar%sat + epsilon)**(alpha)
  auxvar%Ke = Ke

  ! Effective thermal conductivity
  call thermal_cc%thermal_conductivity_function%CalculateTCond( &
       auxvar%sat,auxvar%temp, &
       material_auxvar%porosity,auxvar%Dk_eff,dk_ds,dk_dT,option)

  ! Derivative of soil Kersten number
  auxvar%dKe_dp = alpha*(auxvar%sat + epsilon)**(alpha - 1.d0)* &
                  auxvar%dsat_dp
  auxvar%dKe_dT = 0.d0

  if (option%iflag /= TH_UPDATE_FOR_PERTURBATION) then
    global_auxvar%sat(1) = auxvar%sat
  endif

end subroutine THAuxVarComputeNoFreezing

! ************************************************************************** !

subroutine THAuxVarComputeFreezing(x, auxvar, global_auxvar, &
                                   material_auxvar, &
                                   iphase, &
                                   saturation_function, &
                                   thermal_cc, &
                                   th_parameter, icct, natural_id, &
                                   update_porosity,option)
  !
  ! Computes auxillary variables for each grid cell when
  ! ice and vapor phases are present
  !
  ! Author: Satish Karra, LANL
  ! Date: 11/16/11
  !

!sk: Not sure if we need por, perm

  use Option_module
  use Global_Aux_module

  use EOS_Water_module
  use Saturation_Function_module
  use Characteristic_Curves_Thermal_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  class(cc_thermal_type) :: thermal_cc
  PetscReal :: x(option%nflowdof)
  type(th_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(th_parameter_type) :: th_parameter
  PetscInt :: icct
  PetscInt :: iphase
  PetscInt :: natural_id
  PetscBool :: update_porosity

  PetscErrorCode :: ierr
  PetscReal :: pw, dw_kg, dw_mol, hw, sat_pressure, visl
  PetscReal :: ds_dp, dkr_dp, dkr_dT
  PetscReal :: dvis_dT, dvis_dp
  PetscReal :: dw_dp, dw_dT, hw_dp, hw_dT
  PetscReal :: dpw_dp
  PetscReal :: dpsat_dT
  PetscReal :: ice_saturation, gas_saturation
  PetscReal :: dsl_temp
  PetscReal :: dsg_pl, dsg_temp
  PetscReal :: dsi_pl, dsi_temp
  PetscReal :: den_ice, dden_ice_dT, dden_ice_dP
  PetscReal :: u_ice, du_ice_dT, du_ice_dP
  PetscBool :: out_of_table_flag
  PetscReal :: p_th

  PetscReal :: p_g
  PetscReal :: p_sat
  PetscReal :: mol_g
  PetscReal :: C_g
  PetscReal :: dmolg_dT
  PetscReal, parameter :: C_a = 1.86d-3 ! in MJ/kg/K at 300K
  PetscReal, parameter :: C_wv = 1.005d-3 ! in MJ/kg/K

  PetscReal :: Ke
  PetscReal :: Ke_fr
  PetscReal :: alpha
  PetscReal :: alpha_fr
  PetscReal :: Dk
  PetscReal :: Dk_dry
  PetscReal :: Dk_ice
  PetscReal :: dk_ds, dK_di, dk_dT

  ierr = 0
  out_of_table_flag = PETSC_FALSE

  global_auxvar%sat = UNINITIALIZED_DOUBLE
  global_auxvar%den = UNINITIALIZED_DOUBLE
  global_auxvar%den_kg = UNINITIALIZED_DOUBLE

  auxvar%den = 0.d0
  auxvar%den_kg = 0.d0
  auxvar%sat = 0.d0
  auxvar%h = 0.d0
  auxvar%u = 0.d0
  auxvar%avgmw = 0.d0
  auxvar%mobility = 0.d0

  auxvar%pres = x(1)
  auxvar%temp = x(2)

  ! Check if the capillary pressure is less than -100MPa

  if (auxvar%pres - &
      option%flow%reference_pressure < -1.d8 + 1.d0) then
    auxvar%pres = -1.d8 + option%flow%reference_pressure + 1.d0
  endif

  call THAuxPorosity(auxvar,material_auxvar,global_auxvar, &
                     update_porosity,th_parameter,option)

  auxvar%pc = option%flow%reference_pressure - auxvar%pres

!***************  Liquid phase properties **************************
  auxvar%avgmw = FMWH2O

  pw = option%flow%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0
  if (auxvar%pc > 1.d0) then
    iphase = 3
    dpw_dp = 0.d0
  else
    iphase = 1
    auxvar%pc = 0.d0
    pw = auxvar%pres
    dpw_dp = 1.d0
  endif

  call CapillaryPressureThreshold(saturation_function,p_th,option)

#if defined(CLM_PFLOTRAN) || defined(CLM_OFFLINE)
    if(auxvar%bc_alpha > 0.d0) then
       saturation_function%alpha  = auxvar%bc_alpha
       saturation_function%lambda = auxvar%bc_lambda
    endif
#endif

  ! Check if user specified ice model via thermal characteristic curves
  select type(tcf => thermal_cc%thermal_conductivity_function)
  class is (kt_frozen_type)
    if (Initialized(tcf%ice_model)) then
      th_ice_model = tcf%ice_model
    endif
  class default
    option%io_buffer = 'Cannot use thermal characteristic curve "' &
                       // trim(thermal_cc%name) // &
                       '" when FREEZING is active in TH mode.'
    call PrintErrMsg(option)
  end select

  select case (th_ice_model)
    case (PAINTER_EXPLICIT)
      ! Model from Painter, Comp. Geosci. (2011)
      call SatFuncComputeIcePExplicit(auxvar%pres, &
                                      auxvar%temp, ice_saturation, &
                                      auxvar%sat, gas_saturation, &
                                      auxvar%kr, &
                                      ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                      dsi_pl, dsi_temp, dkr_dp, dkr_dT, &
                                      saturation_function, p_th, option)
    case (PAINTER_KARRA_IMPLICIT)
      ! Implicit model from Painter & Karra, VJZ (2013)
      call SatFuncComputeIcePKImplicit(auxvar%pres, &
                                       auxvar%temp, ice_saturation, &
                                       auxvar%sat, gas_saturation, &
                                       auxvar%kr, &
                                       ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                       dsi_pl, dsi_temp, dkr_dp, dkr_dT, &
                                       saturation_function, p_th, option)
    case (PAINTER_KARRA_EXPLICIT)
      ! Explicit model from Painter & Karra, VJZ (2013)
      call SatFuncComputeIcePKExplicit(auxvar%pres, &
                                       auxvar%temp, ice_saturation, &
                                       auxvar%sat, gas_saturation, &
                                       auxvar%kr, &
                                       ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                       dsi_pl, dsi_temp, dkr_dp, dkr_dT, &
                                       saturation_function, p_th, option)
    case (DALL_AMICO)
      ! Model from Dall'Amico (2010) and Dall' Amico et al. (2011)
      call SatFuncComputeIceDallAmico(auxvar%pres, &
                                      auxvar%temp, &
                                      auxvar%ice%pres_fh2o, &
                                      auxvar%ice%dpres_fh2o_dp, &
                                      auxvar%ice%dpres_fh2o_dT, &
                                      ice_saturation, &
                                      auxvar%sat, gas_saturation, &
                                      auxvar%kr, &
                                      ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                      dsi_pl, dsi_temp, dkr_dp, dkr_dT, &
                                      saturation_function, option)
    case (PAINTER_KARRA_EXPLICIT_NOCRYO)
      ! Explicit model from Painter & Karra, VJZ (2013) and removed cryosuction
      call SatFuncComputeIcePKExplicitNoCryo(auxvar%pres, &
                                       auxvar%temp, ice_saturation, &
                                       auxvar%sat, gas_saturation, &
                                       auxvar%kr, &
                                       ds_dp, dsl_temp, dsg_pl, dsg_temp, &
                                       dsi_pl, dsi_temp, dkr_dp, dkr_dT, &
                                       saturation_function, p_th, option)
    case default
      option%io_buffer = 'THCAuxVarComputeIce: Ice model not recognized.'
      call PrintErrMsg(option)
  end select

  call EOSWaterDensity(auxvar%temp,pw,dw_kg,dw_mol,dw_dp,dw_dT,ierr)
  if (ierr /= 0) then
    call PrintMsgByCell(option,natural_id, &
                        'Error in THAuxVarComputeFreezing->EOSWaterDensity')
  endif
  call EOSWaterEnthalpy(auxvar%temp,pw,hw,hw_dp,hw_dT,ierr)
  ! J/kmol -> MJ/kmol
  hw = hw * option%scale
  hw_dp = hw_dp * option%scale
  hw_dT = hw_dT * option%scale

  call EOSWaterSaturationPressure(auxvar%temp, sat_pressure, &
                                  dpsat_dT, ierr)
  call EOSWaterViscosity(auxvar%temp, pw, sat_pressure, dpsat_dT, &
                         visl, dvis_dT,dvis_dp, ierr)

  if (iphase == 3) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif

  auxvar%den = dw_mol
  auxvar%den_kg = dw_kg

  auxvar%h = hw
  auxvar%u = auxvar%h - pw / dw_mol * option%scale
  auxvar%mobility = auxvar%kr/visl
  auxvar%vis = visl
  auxvar%dsat_dp = ds_dp
  auxvar%dden_dT = dw_dT
  auxvar%dden_dp = dw_dp
!geh: contribution of dvis_dpsat is now added in EOSWaterViscosity
!  auxvar%dmobility_dT = -auxvar%kr/(visl*visl)*(dvis_dT + dvis_dpsat*dpsat_dT) + &
!    dkr_dT/visl
  auxvar%dmobility_dT = -auxvar%kr/(visl*visl)*dvis_dT + dkr_dT/visl
  auxvar%dmobility_dp = dkr_dp/visl - auxvar%kr/(visl*visl)*dvis_dp
  auxvar%dh_dp = hw_dp
  auxvar%du_dp = hw_dp - (dpw_dp/dw_mol - pw/(dw_mol*dw_mol)*dw_dp)* &
                  option%scale
  auxvar%dh_dT = hw_dT
  auxvar%du_dT = hw_dT + pw/(dw_mol*dw_mol)*option%scale*dw_dT

  auxvar%ice%sat_ice = ice_saturation
  auxvar%ice%sat_gas = gas_saturation
  auxvar%dsat_dT = dsl_temp
  auxvar%ice%dsat_ice_dp = dsi_pl
  auxvar%ice%dsat_gas_dp = dsg_pl
  auxvar%ice%dsat_ice_dT = dsi_temp
  auxvar%ice%dsat_gas_dT = dsg_temp

  ! Calculate the density, internal energy and derivatives for ice
  call EOSWaterDensityIce(auxvar%temp, auxvar%pres, &
                          den_ice, dden_ice_dT, dden_ice_dP, ierr)
  if (ierr /= 0) then
    call PrintMsgByCell(option,natural_id, &
                      'Error in THAuxVarComputeFreezing->EOSWaterDensityIce')
  endif
  call EOSWaterInternalEnergyIce(auxvar%temp, u_ice, du_ice_dT, du_ice_dP, ierr)

  auxvar%ice%den_ice = den_ice
  auxvar%ice%dden_ice_dT = dden_ice_dT
  auxvar%ice%dden_ice_dp = dden_ice_dP
  auxvar%ice%u_ice = u_ice*1.d-3                  !kJ/kmol --> MJ/kmol
  auxvar%ice%du_ice_dT = du_ice_dT*1.d-3          !kJ/kmol/K --> MJ/kmol/K

  ! Calculate the values and derivatives for density and internal energy
  call EOSWaterSaturationPressure(auxvar%temp, p_sat, ierr)

  p_g            = option%flow%reference_pressure
  auxvar%ice%den_gas = p_g/(IDEAL_GAS_CONSTANT* &
                         (auxvar%temp + T273K))*1.d-3 !in kmol/m3
  mol_g          = p_sat/p_g
  C_g            = C_wv*mol_g*FMWH2O + C_a*(1.d0 - mol_g)*FMWAIR ! in MJ/kmol/K
  auxvar%ice%u_gas   = C_g*(auxvar%temp + T273K)       ! in MJ/kmol
  auxvar%ice%mol_gas = mol_g

  auxvar%ice%dden_gas_dT = -p_g/(IDEAL_GAS_CONSTANT* &
                            (auxvar%temp + T273K)**2)*1.d-3
  dmolg_dT           = dpsat_dT/p_g
  auxvar%ice%du_gas_dT = C_g + (C_wv*dmolg_dT*FMWH2O - C_a*dmolg_dT*FMWAIR)* &
                       (auxvar%temp + T273K)
  auxvar%ice%dmol_gas_dT = dmolg_dT

  ! Parameters for computation of effective thermal conductivity
  alpha = th_parameter%alpha(icct)
  alpha_fr = th_parameter%alpha_fr(icct)
  Dk = th_parameter%ckwet(icct)
  Dk_dry = th_parameter%ckdry(icct)
  Dk_ice = th_parameter%ckfrozen(icct)

  !Soil Kersten number
  Ke = (auxvar%sat + epsilon)**(alpha)
  Ke_fr = (auxvar%ice%sat_ice + epsilon)**(alpha_fr)
  auxvar%Ke = Ke
  auxvar%ice%Ke_fr = Ke_fr

  ! Effective thermal conductivity
  call thermal_cc%thermal_conductivity_function%CalculateFTCond( &
       auxvar%sat,auxvar%ice%sat_ice,auxvar%temp, &
       auxvar%effective_porosity,auxvar%Dk_eff,dk_ds,dK_di,dk_dT,option)

  ! Derivative of Kersten number
  auxvar%dKe_dp = alpha*(auxvar%sat + epsilon)**(alpha - 1.d0)* &
                  auxvar%dsat_dp
  auxvar%dKe_dT = alpha*(auxvar%sat + epsilon)**(alpha - 1.d0)* &
                  auxvar%dsat_dT
  auxvar%ice%dKe_fr_dT = alpha_fr* &
                         (auxvar%ice%sat_ice + epsilon)**(alpha_fr - 1.d0)* &
                         auxvar%ice%dsat_ice_dT
  auxvar%ice%dKe_fr_dp = alpha_fr* &
                         (auxvar%ice%sat_ice + epsilon)**(alpha_fr - 1.d0)* &
                         auxvar%ice%dsat_ice_dp

  if (th_ice_model == DALL_AMICO) then
    auxvar%ice%den_ice = dw_mol
    auxvar%ice%dden_ice_dT = auxvar%dden_dT
    auxvar%ice%dden_ice_dp = auxvar%dden_dp
!    auxvar%ice%u_ice = auxvar%u  ! commented out by S.Karra 06/02/14. setting
!    internal energy of ice and water might not be correct.
!    auxvar%ice%du_ice_dT = auxvar%du_dT

    auxvar%ice%sat_gas       = 0.d0
    auxvar%ice%dsat_gas_dp   = 0.d0
    auxvar%ice%dsat_gas_dT   = 0.d0
    auxvar%ice%den_gas       = 0.d0
    auxvar%ice%dden_gas_dT   = 0.d0
    auxvar%ice%u_gas         = 0.d0
    auxvar%ice%du_gas_dT     = 0.d0
    auxvar%ice%mol_gas       = 0.d0
    auxvar%ice%dmol_gas_dT   = 0.d0
  endif

end subroutine THAuxVarComputeFreezing

! ************************************************************************** !
subroutine THAuxVarCompute2ndOrderDeriv(th_auxvar,global_auxvar, &
                                        material_auxvar,th_parameter, &
                                        icct,characteristic_curves,&
                                        thermal_cc,&
                                        option)

  ! Computes 2nd order derivatives auxiliary variables for each grid cell
  !
  ! Author: Satish Karra
  ! Date: 06/06/2019
  !

  use Option_module
  use Global_Aux_module

  use EOS_Water_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Thermal_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  class(cc_thermal_type) :: thermal_cc
  type(th_auxvar_type) :: th_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: icct

  type(th_parameter_type) :: th_parameter
  PetscInt :: iphase, ideriv
  type(th_auxvar_type) :: th_auxvar_pert
  type(global_auxvar_type) :: global_auxvar_pert
  ! leave as type
  type(material_auxvar_type) :: material_auxvar_pert
  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), pert

  call GlobalAuxVarInit(global_auxvar_pert,option)
  call MaterialAuxVarInit(material_auxvar_pert,option)

  call THAuxVarCopy(th_auxvar,th_auxvar_pert,option)
  call GlobalAuxVarCopy(global_auxvar,global_auxvar_pert,option)
  call MaterialAuxVarCopy(material_auxvar,material_auxvar_pert,option)

  x(1) = th_auxvar%pres
  x(2) = th_auxvar%temp

!  if (option%use_th_freezing) then
!    option%io_buffer = 'ERROR: TH_TS MODE not implemented with freezing'
!    call PrintErrMsg(option)
!  else
!    call THAuxVarComputeNoFreezing(x,th_auxvar,&
!                      global_auxvar,material_auxvar,&
!                      iphase,sat_func, &
!                      th_parameter,icct, &
!                      -999,option)
!  endif

  th_auxvar%d2sat_dp2 = 0.d0
  th_auxvar%d2den_dp2 = 0.d0
  th_auxvar%d2u_dp2 = 0.d0
  th_auxvar%d2sat_dT2 = 0.d0
  th_auxvar%d2den_dT2 = 0.d0
  th_auxvar%d2u_dT2 = 0.d0


  do ideriv = 1,option%nflowdof
    pert = x(ideriv)*perturbation_tolerance
    x_pert = x
    if (option%flow%th_freezing) then
       if (ideriv == 1) then
          if (x_pert(ideriv) < option%flow%reference_pressure) then
             pert = - pert
          endif
          x_pert(ideriv) = x_pert(ideriv) + pert
       endif

       if (ideriv == 2) then
          if (x_pert(ideriv) < 0.d0) then
             pert = - 1.d-8
          else
             pert =  1.d-8
          endif
          x_pert(ideriv) = x_pert(ideriv) + pert
       endif
    else
       x_pert(ideriv) = x_pert(ideriv) + pert
    endif

    if (option%flow%th_freezing) then
      option%io_buffer = 'ERROR: TH_TS MODE not implemented with freezing'
      call PrintErrMsg(option)
    else
      call THAuxVarComputeNoFreezing(x_pert,th_auxvar_pert,&
                            global_auxvar_pert,material_auxvar_pert,&
                            iphase,characteristic_curves, &
                            thermal_cc, &
                            th_parameter,icct, &
                            -999,PETSC_TRUE,option)
    endif


    if (ideriv == 1) then
      th_auxvar%d2sat_dp2 = (th_auxvar_pert%dsat_dp - th_auxvar%dsat_dp)/pert
      th_auxvar%d2sat_dTdp = (th_auxvar_pert%dsat_dT - th_auxvar%dsat_dT)/pert
      th_auxvar%d2den_dp2 = (th_auxvar_pert%dden_dp - th_auxvar%dden_dp)/pert
      th_auxvar%d2den_dTdp = (th_auxvar_pert%dden_dT - th_auxvar%dden_dT)/pert
      th_auxvar%d2u_dp2 = (th_auxvar_pert%du_dp - th_auxvar%du_dp)/pert
      th_auxvar%d2u_dTdp = (th_auxvar_pert%du_dT - th_auxvar%du_dT)/pert
    endif

    if (ideriv == 2) then
      th_auxvar%d2sat_dT2 = (th_auxvar_pert%dsat_dT - th_auxvar%dsat_dT)/pert
      th_auxvar%d2den_dT2 = (th_auxvar_pert%dden_dT - th_auxvar%dden_dT)/pert
      th_auxvar%d2u_dT2 = (th_auxvar_pert%du_dT - th_auxvar%du_dT)/pert
    endif

  enddo

end subroutine THAuxVarCompute2ndOrderDeriv

! ************************************************************************** !

subroutine THAuxVarPerturb(th_auxvar,global_auxvar, &
                           material_auxvar,iphase, &
                           characteristic_curves, &
                           thermal_cc,th_parameter,icct, &
                           natural_id,update_porosity,option)
  ! Calculates auxiliary variables for perturbed TH system
  !
  ! Author: Glenn Hammond
  ! Date: 08/04/25


  use Option_module
  use Global_Aux_module

  use EOS_Water_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module
  use Characteristic_Curves_Thermal_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  class(cc_thermal_type) :: thermal_cc
  type(th_auxvar_type) :: th_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscInt :: iphase
  type(th_parameter_type) :: th_parameter
  PetscInt :: icct
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  PetscBool :: update_porosity

  PetscInt :: idof
  PetscReal :: x(option%nflowdof)
  PetscReal :: x_pert(option%nflowdof)
  PetscReal :: pert

  option%iflag = TH_UPDATE_FOR_PERTURBATION
  x(TH_PRESSURE_DOF) = th_auxvar%pres
  x(TH_TEMPERATURE_DOF) = th_auxvar%temp
  do idof = 1, option%nflowdof
    pert = x(idof)*th_rel_pert+th_min_pert(idof)
    th_auxvar%auxvar_pert(idof)%pert = pert
    x_pert = x
    x_pert(idof) = x(idof) + pert
    call THAuxVarComputeNoFreezing(x_pert, &
                                   th_auxvar%auxvar_pert(idof), &
                                   global_auxvar,material_auxvar,iphase, &
                                   characteristic_curves, &
                                   thermal_cc,th_parameter,icct, &
                                   natural_id,update_porosity,option)
  enddo
  option%iflag = UNINITIALIZED_INTEGER

end subroutine THAuxVarPerturb

! ************************************************************************** !

recursive subroutine THWellAuxVarInit(auxvar,option,allocate_perturbation)
  !
  ! Initialize well auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25
  !

  use Option_module

  implicit none

  type(th_well_auxvar_type) :: auxvar
  type(option_type) :: option
  PetscBool :: allocate_perturbation

  PetscInt :: i

  auxvar%ghosted_id = 0
  auxvar%z = 0.d0
  auxvar%temp = 0.d0
  auxvar%u = 0.d0
  auxvar%h = 0.d0
  auxvar%du_dT = 0.d0
  auxvar%dh_dT = 0.d0
  auxvar%volume = 0.d0
  auxvar%segment_length = 0.d0
  auxvar%heat_transfer_coef = 0.d0
  auxvar%well_index = 0.d0
  auxvar%therm_cond_borehole_to_cell = 0.d0
  auxvar%den = 0.d0
  auxvar%den_kg = 0.d0
  auxvar%vis = 0.d0
  auxvar%therm_cond_fluid = 0.d0
  auxvar%spec_heat_fluid = 0.d0
  auxvar%dden_dT = 0.d0
  auxvar%dtherm_cond_fluid_dT = 0.d0
  auxvar%reynolds = 0.d0
  auxvar%prandtl = 0.d0
  auxvar%nusselt = 0.d0

  auxvar%pert = 0.d0
  if (allocate_perturbation) then
    allocate(auxvar%auxvar_pert(option%nflowdof))
    do i = 1, option%nflowdof
      call THWellAuxVarInit(auxvar%auxvar_pert(i),option,PETSC_FALSE)
    enddo
  else
    nullify(auxvar%auxvar_pert)
  endif

end subroutine THWellAuxVarInit

! ************************************************************************** !

recursive subroutine THWellAuxVarCopyParamsToPert(auxvar)
  !
  ! Copy static parameters to perturbed well auxvar
  !
  ! Author: Glenn Hammond
  ! Date: 09/25/25
  !
  implicit none

  type(th_well_auxvar_type) :: auxvar

  PetscInt :: i

  do i = 1, size(auxvar%auxvar_pert)
    auxvar%auxvar_pert(i)%ghosted_id = auxvar%ghosted_id
    auxvar%auxvar_pert(i)%z = auxvar%z
    auxvar%auxvar_pert(i)%volume = auxvar%volume
    auxvar%auxvar_pert(i)%segment_length = auxvar%segment_length
    auxvar%auxvar_pert(i)%heat_transfer_coef = auxvar%heat_transfer_coef
    auxvar%auxvar_pert(i)%well_index = auxvar%well_index
    auxvar%auxvar_pert(i)%therm_cond_borehole_to_cell = &
      auxvar%therm_cond_borehole_to_cell
    auxvar%auxvar_pert(i)%den = auxvar%den
    auxvar%auxvar_pert(i)%den_kg = auxvar%den_kg
    auxvar%auxvar_pert(i)%vis = auxvar%vis
    auxvar%auxvar_pert(i)%therm_cond_fluid = auxvar%therm_cond_fluid
    auxvar%auxvar_pert(i)%spec_heat_fluid = auxvar%spec_heat_fluid
    auxvar%auxvar_pert(i)%reynolds = auxvar%reynolds
    auxvar%auxvar_pert(i)%prandtl = auxvar%prandtl
    auxvar%auxvar_pert(i)%nusselt = auxvar%nusselt
    auxvar%auxvar_pert(i)%dden_dT = auxvar%dden_dT
    auxvar%auxvar_pert(i)%dtherm_cond_fluid_dT = auxvar%dtherm_cond_fluid_dT
  enddo

end subroutine THWellAuxVarCopyParamsToPert

! ************************************************************************** !

subroutine THWellAuxVarCompute(liquid_pressure,x,auxvar,natural_id,option)
  !
  ! Computes secondary state variables for well auxvar
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  use Option_module
  use EOS_Water_module

  implicit none

  PetscReal :: liquid_pressure
  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  type(th_well_auxvar_type) :: auxvar
  PetscInt :: natural_id

  PetscReal :: pw, dw_kmol, hw
  PetscReal :: dw_dT, dw_dp, dw_kg
  PetscReal :: hw_dp, hw_dT
  PetscReal :: vis, dvis_dp, dvis_dT
  PetscReal :: sat_pressure, dpsat_dT
  PetscReal :: therm_cond, dtherm_cond_dp, dtherm_cond_dT
  PetscErrorCode :: ierr

  pw = liquid_pressure
  auxvar%temp = x(th_well_dof)

  ! Compute density [kg/m^3] and [kmol/m^3]
  call EOSWaterDensity(auxvar%temp,pw,dw_kg,dw_kmol,dw_dp,dw_dT,ierr)
  if (ierr /= 0) then
    call PrintMsgByCell(option,natural_id, &
      'Error in THWellAuxVarCompute->EOSWaterDensity')
  endif
  auxvar%den = dw_kmol
  auxvar%den_kg = dw_kg
  auxvar%dden_dT = dw_dT

  call EOSWaterEnthalpy(auxvar%temp,pw,hw,hw_dp,hw_dT,ierr)
  if (ierr /= 0) then
    call PrintMsgByCell(option,natural_id, &
      'Error in THWellAuxVarCompute->EOSWaterEnthalpy')
  endif

  ! Compute specific heat capacity [J/kg-K] from enthalpy derivative
  ! c_p = dh/dT
  ! We have hw_dT in [J/kmol-K], convert to [J/kg-K]
  auxvar%spec_heat_fluid = hw_dT / FMWH2O

  ! J/kmol -> whatever units (default of option%scale is 1.e-6)
  hw = hw * option%scale
  hw_dp = hw_dp * option%scale
  hw_dT = hw_dT * option%scale
  auxvar%h = hw
  auxvar%u = auxvar%h - pw / dw_kmol * option%scale
  auxvar%dh_dT = hw_dT
  auxvar%du_dT = hw_dT + pw/(dw_kmol*dw_kmol)*dw_dT*option%scale

  ! Compute saturation pressure (needed for viscosity calculation)
  call EOSWaterSaturationPressure(auxvar%temp,sat_pressure,dpsat_dT,ierr)
  if (ierr /= 0) then
    call PrintMsgByCell(option,natural_id, &
      'Error in THWellAuxVarCompute->EOSWaterSaturationPressure')
  endif

  ! Compute viscosity [Pa-s]
  call EOSWaterViscosity(auxvar%temp,pw,sat_pressure,dpsat_dT,vis, &
                         dvis_dT,dvis_dp,ierr)
  if (ierr /= 0) then
    call PrintMsgByCell(option,natural_id, &
      'Error in THWellAuxVarCompute->EOSWaterViscosity')
  endif
  auxvar%vis = vis

  ! Compute thermal conductivity of  water [MW/m-K]
  call EOSWaterThermalConductivityIF97(dw_kg,auxvar%temp,therm_cond, &
                                       PETSC_TRUE,dw_dp,dw_dT, &
                                       dtherm_cond_dp,dtherm_cond_dT)
  auxvar%therm_cond_fluid = therm_cond
  auxvar%dtherm_cond_fluid_dT = dtherm_cond_dT

end subroutine THWellAuxVarCompute

! ************************************************************************** !

subroutine THWellAuxVarPerturb(liquid_pressure,auxvar_well, &
                               natural_id,option)
  !
  ! Calculates well auxiliary variables for perturbed TH system
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25

  use Option_module

  implicit none

  PetscReal :: liquid_pressure
  type(th_well_auxvar_type) :: auxvar_well
  PetscInt :: natural_id
  type(option_type) :: option

  PetscInt :: idof
  PetscReal :: x(option%nflowdof)
  PetscReal :: x_pert(option%nflowdof)
  PetscReal :: pert

  x = UNINITIALIZED_DOUBLE
  x(th_well_dof) = auxvar_well%temp
  do idof = 1, option%nflowdof
    pert = x(idof)*th_rel_pert+th_min_pert(idof)
    auxvar_well%auxvar_pert(idof)%pert = pert
    x_pert(1:option%nflowdof) = x
    x_pert(idof) = x(idof) + pert
    call THWellAuxVarCompute(liquid_pressure,x_pert, &
                             auxvar_well%auxvar_pert(idof), &
                             natural_id,option)
  enddo

end subroutine THWellAuxVarPerturb

! ************************************************************************** !

subroutine THPrintAuxVars(file_unit,th_auxvar,global_auxvar, &
                          material_auxvar,natural_id,string,option)

  ! Prints the content of an TH auxvar to the designated file unit

  ! Author: Glenn Hammond
  ! Date: 07/16/20

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  PetscInt :: file_unit
  type(th_auxvar_type) :: th_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  type(option_type) :: option

  PetscReal ::  liquid_mass

  liquid_mass = material_auxvar%volume*th_auxvar%effective_porosity* &
                th_auxvar%sat*th_auxvar%den

  write(file_unit,*) '--------------------------------------------------------'
  if (len_trim(string) > 0) write(file_unit,*) trim(string)
  write(file_unit,*) '                  cell id: ', natural_id
  write(file_unit,*) '       liquid mass [kmol]: ', liquid_mass
  write(file_unit,*) '          liquid pressure: ', th_auxvar%pres
  write(file_unit,*) '       capillary pressure: ', th_auxvar%pc
  write(file_unit,*) '          temperature [C]: ', th_auxvar%temp
  write(file_unit,*) '          liquid enthalpy: ', th_auxvar%h
  write(file_unit,*) '   liquid internal energy: ', th_auxvar%u
  write(file_unit,*) '         liquid viscosity: ', th_auxvar%vis
  write(file_unit,*) '          liquid mobility: ', th_auxvar%mobility
  write(file_unit,*) '     liquid relative perm: ', th_auxvar%kr
  write(file_unit,*) '        liquid saturation: ', th_auxvar%sat
  write(file_unit,*) '    liquid_density [kmol]: ', th_auxvar%den
  write(file_unit,*) '      liquid_density [kg]: ', th_auxvar%den_kg
  write(file_unit,*) 'eff. thermal conductivity: ', th_auxvar%Dk_eff
  write(file_unit,*) '      effective porosity : ', th_auxvar%effective_porosity
  write(file_unit,*) '             volume [m^3]: ', material_auxvar%volume
  write(file_unit,*) '--------------------------------------------------------'

end subroutine THPrintAuxVars

! ************************************************************************** !

recursive subroutine THWellAuxVarStrip(auxvar_well)
  !
  ! Deallocates a TH auxiliary object
  !
  ! Author: ???
  ! Date: 02/14/08
  !

  implicit none

  type(th_well_auxvar_type) :: auxvar_well

  if (associated(auxvar_well%auxvar_pert)) then
    deallocate(auxvar_well%auxvar_pert)
    nullify(auxvar_well%auxvar_pert)
  endif

end subroutine THWellAuxVarStrip

! ************************************************************************** !

recursive subroutine THAuxVarDestroy(auxvar)
  !
  ! Deallocates a TH auxiliary object
  !
  ! Author: ???
  ! Date: 02/14/08
  !

  implicit none

  type(th_auxvar_type) :: auxvar

  if (associated(auxvar%ice)) deallocate(auxvar%ice)
  nullify(auxvar%ice)
!  call THWellAuxVarDestroy(auxvar%auxvar_well)

end subroutine THAuxVarDestroy

! ************************************************************************** !

subroutine THAuxDestroy(aux)
  !
  ! Deallocates a TH auxiliary object
  !
  ! Author: ???
  ! Date: 02/14/08
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(TH_type), pointer :: aux
  PetscInt :: iaux

  call DeallocateArray(th_min_pert)

  if (.not.associated(aux)) return

  do iaux = 1, aux%num_aux
    call THAuxVarDestroy(aux%auxvars(iaux))
  enddo
  do iaux = 1, aux%num_aux_bc
    call THAuxVarDestroy(aux%auxvars_bc(iaux))
  enddo
  do iaux = 1, aux%num_aux_ss
    call THAuxVarDestroy(aux%auxvars_ss(iaux))
  enddo
  do iaux = 1, aux%num_well_aux
    call THWellAuxVarStrip(aux%auxvars_well(iaux))
  enddo
  do iaux = 1, aux%num_well_aux_bc
    call THWellAuxVarStrip(aux%auxvars_well_bc(iaux))
  enddo

  if (associated(aux%auxvars)) deallocate(aux%auxvars)
  nullify(aux%auxvars)
  if (associated(aux%auxvars_bc)) deallocate(aux%auxvars_bc)
  nullify(aux%auxvars_bc)
  if (associated(aux%auxvars_ss)) deallocate(aux%auxvars_ss)
  nullify(aux%auxvars_ss)
  if (associated(aux%auxvars_well)) deallocate(aux%auxvars_well)
  nullify(aux%auxvars_well)
  if (associated(aux%auxvars_well_bc)) deallocate(aux%auxvars_well_bc)
  nullify(aux%auxvars_well_bc)

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%th_parameter)) then
    if (associated(aux%th_parameter%dencpr)) deallocate(aux%th_parameter%dencpr)
    nullify(aux%th_parameter%dencpr)
    if (associated(aux%th_parameter%ckwet)) deallocate(aux%th_parameter%ckwet)
    nullify(aux%th_parameter%ckwet)
    if (associated(aux%th_parameter%ckdry)) deallocate(aux%th_parameter%ckdry)
    nullify(aux%th_parameter%ckdry)
    if (associated(aux%th_parameter%alpha)) deallocate(aux%th_parameter%alpha)
    nullify(aux%th_parameter%alpha)
    ! ice
    if (associated(aux%th_parameter%ckfrozen)) &
      deallocate(aux%th_parameter%ckfrozen)
    nullify(aux%th_parameter%ckfrozen)
    if (associated(aux%th_parameter%alpha_fr)) &
      deallocate(aux%th_parameter%alpha_fr)
    nullify(aux%th_parameter%alpha_fr)

    if (associated(aux%th_parameter%sir)) deallocate(aux%th_parameter%sir)
    nullify(aux%th_parameter%sir)

    ! solely nullify geomech_parameter as it is destroyed elsewhere
    nullify(aux%th_parameter%geomech_parameter)
  endif
  nullify(aux%th_parameter)

  deallocate(aux)
  nullify(aux)

  end subroutine THAuxDestroy

end module TH_Aux_module

