module THC_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use THC_Aux_module
  use THC_Common_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: THCSetup, &
            THCInitializeTimestep, &
            THCUpdateSolution, &
            THCTimeCut, &
            THCUpdateAuxVars, &
            THCUpdateFixedAccum, &
            THCComputeMassBalance, &
            THCZeroMassBalanceDelta, &
            THCResidual, &
            THCSetPlotVariables, &
            THCMapBCAuxVarsToGlobal, &
            THCDestroy

contains

! ************************************************************************** !

subroutine THCSetup(realization)
  !
  ! Creates arrays for auxiliary variables
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Fluid_module
  use Grid_module
  use Material_Aux_module
  use Output_Aux_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Thermal_module
  use Matrix_Zeroing_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(output_variable_list_type), pointer :: list
  type(material_parameter_type), pointer :: material_parameter
  type(thc_parameter_type), pointer :: thc_parameter
  type(fluid_property_type), pointer :: cur_fluid_property

  PetscInt :: ghosted_id, iconn, sum_connection, local_id
  PetscBool :: error_found
  PetscInt :: flag(10)
  PetscInt :: temp_int, idof, imat, icct
  PetscBool, allocatable :: dof_is_active(:)
  PetscErrorCode :: ierr
                                                ! extra index for derivatives
  type(thc_auxvar_type), pointer :: thc_auxvars(:,:)
  type(thc_auxvar_type), pointer :: thc_auxvars_bc(:)
  type(thc_auxvar_type), pointer :: thc_auxvars_ss(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%THC => THCAuxCreate(option)
  thc_parameter => patch%aux%THC%thc_parameter

  temp_int = size(patch%material_property_array)
  if (thc_tensorial_rel_perm) then
    allocate(thc_parameter%tensorial_rel_perm_exponent(3,temp_int))
    thc_parameter%tensorial_rel_perm_exponent = UNINITIALIZED_DOUBLE
    do imat = 1, temp_int
      if (Initialized(minval(patch%material_property_array(imat)%ptr% &
                                      tensorial_rel_perm_exponent))) then
        thc_parameter%tensorial_rel_perm_exponent(:,imat) = &
          patch%material_property_array(imat)%ptr% &
            ! the tortuosity parameter in hardwired to 0.5 in characteristic
            ! curves. we subtract the default to allow the tensorial value
            ! to override the hardwired default
            tensorial_rel_perm_exponent - 0.5d0
      else
        option%io_buffer = 'A tensorial relative permeability exponent &
          &is not define for material "' // &
          trim(patch%material_property_array(imat)%ptr%name) // '".'
        call PrintErrMsg(option)
      endif
    enddo
  else
    ! check to ensure that user has not parameterized tensorial perm without
    ! adding TENSORIAL_RELATIVE_PERMEABILITY to the simulation OPTIONS block
    do imat = 1, temp_int
      if (Initialized(maxval(patch%material_property_array(imat)%ptr% &
                                      tensorial_rel_perm_exponent))) then
        option%io_buffer = 'A tensorial relative permeability exponent &
          &is define for material "' // &
          trim(patch%material_property_array(imat)%ptr%name) // '" without &
          &TENSORIAL_RELATIVE_PERMEABILITY being defined in the THC &
          &simulation OPTIONS block.'
        call PrintErrMsg(option)
      endif
    enddo
  endif

  ! --------------------------------------------------------------------------
  ! Per-material thermal properties (TH-compatible MATERIAL_PROPERTY card).
  !   dencpr = ROCK_DENSITY * SPECIFIC_HEAT          [J/(m^3.K)]
  !   ckdry  = THERMAL_CONDUCTIVITY_DRY              [W/(m.K)]
  !   ckwet  = THERMAL_CONDUCTIVITY_WET              [W/(m.K)]
  ! Each falls back to the MODE THC OPTIONS-block globals when the material
  ! does not supply it, so existing inputs (OPTIONS-only) are unchanged.  The
  ! conductivity path used per-cell is selected in THCAuxVarCompute: bulk
  ! dry/wet interpolation when both ckdry/ckwet are set, otherwise the
  ! grain-Somerton model from the global thc_kappa_solid.
  ! --------------------------------------------------------------------------
  ! longitudinal dispersivity; transverse dispersion is not implemented
  allocate(thc_parameter%alpha_l(temp_int))
  do imat = 1, temp_int
    if (maxval(patch%material_property_array(imat)%ptr% &
                 dispersivity(2:3)) > 0.d0) then
      option%io_buffer = 'Transverse dispersion is not implemented in &
        &MODE THC; only LONGITUDINAL_DISPERSIVITY is supported (material "' // &
        trim(patch%material_property_array(imat)%ptr%name) // '").'
      call PrintErrMsg(option)
    endif
    thc_parameter%alpha_l(imat) = &
      patch%material_property_array(imat)%ptr%dispersivity(1)
  enddo
  allocate(thc_parameter%dencpr(temp_int))
  allocate(thc_parameter%ckdry(temp_int))
  allocate(thc_parameter%ckwet(temp_int))
  ! per-material THERMAL_CHARACTERISTIC_CURVES, when supplied.  This is the
  ! only conductivity path that returns d(kappa_eff)/dT (POWER,
  ! CUBIC_POLYNOMIAL, LINEAR_RESISTIVITY); DEFAULT reproduces ckdry/ckwet.
  if (associated(patch%char_curves_thermal_array)) then
    allocate(thc_parameter%thermal_cc(temp_int))
    do imat = 1, temp_int
      nullify(thc_parameter%thermal_cc(imat)%ptr)
      icct = patch%material_property_array(imat)%ptr% &
               thermal_conductivity_function_id
      if (Initialized(icct)) then
        thc_parameter%thermal_cc(imat)%ptr => &
          patch%char_curves_thermal_array(icct)%ptr
      endif
    enddo
  endif
  do imat = 1, temp_int
    ! solid volumetric heat capacity rho_s * c_s
    if (Initialized(patch%material_property_array(imat)%ptr%rock_density) .and. &
        Initialized(patch%material_property_array(imat)%ptr%specific_heat)) then
      thc_parameter%dencpr(imat) = &
        patch%material_property_array(imat)%ptr%rock_density * &
        patch%material_property_array(imat)%ptr%specific_heat
    else
      thc_parameter%dencpr(imat) = thc_density_solid * &
                                       thc_specific_heat_solid
    endif
    ! bulk dry/wet thermal conductivity.  Precedence:
    !   1. per-material MATERIAL_PROPERTY THERMAL_CONDUCTIVITY_DRY/WET
    !   2. OPTIONS-block global thc_kappa_dry/wet (bulk)
    !   3. UNINITIALIZED -> grain-Somerton fallback (thc_kappa_solid) in
    !      THCAuxVarCompute
    if (Initialized(patch%material_property_array(imat)%ptr% &
                      thermal_conductivity_dry) .and. &
        Initialized(patch%material_property_array(imat)%ptr% &
                      thermal_conductivity_wet)) then
      thc_parameter%ckdry(imat) = &
        patch%material_property_array(imat)%ptr%thermal_conductivity_dry
      thc_parameter%ckwet(imat) = &
        patch%material_property_array(imat)%ptr%thermal_conductivity_wet
    else if (Initialized(thc_kappa_dry) .and. &
             Initialized(thc_kappa_wet)) then
      thc_parameter%ckdry(imat) = thc_kappa_dry
      thc_parameter%ckwet(imat) = thc_kappa_wet
    else
      thc_parameter%ckdry(imat) = UNINITIALIZED_DOUBLE
      thc_parameter%ckwet(imat) = UNINITIALIZED_DOUBLE
    endif
  enddo

  ! THC always solves the fixed 3-DOF (P,T,C) system; no runtime DOF-index
  ! assignment is needed (contrast ZFLOW, which renumbers its equation indices
  ! here).  thc_pressure_dof / thc_temperature_dof /
  ! thc_concentration_dof are compile-time parameters.

  thc_numerical_derivatives = option%flow%numerical_derivatives
  if (thc_numerical_derivatives) then
    allocate(thc_min_pert(THC_NDOF))
    thc_min_pert(thc_pressure_dof)      = thc_pres_min_pert
    thc_min_pert(thc_temperature_dof)   = thc_temp_min_pert
    thc_min_pert(thc_concentration_dof) = thc_conc_min_pert
    allocate(patch%aux%THC%material_auxvars_pert(ONE_INTEGER,grid%ngmax))
    do ghosted_id = 1, grid%ngmax
      call MaterialAuxVarInit(patch%aux%THC% &
            material_auxvars_pert(ONE_INTEGER,ghosted_id),option)
    enddo
  else
    allocate(patch%aux%THC%material_auxvars_pert(ZERO_INTEGER,grid%ngmax))
  endif

  ! ensure that material properties specific to this module are properly
  ! initialized
  material_parameter => patch%aux%Material%material_parameter
  error_found = PETSC_FALSE

  material_auxvars => patch%aux%Material%auxvars
  flag = 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (material_auxvars(ghosted_id)%volume < 0.d0 .and. flag(1) == 0) then
      flag(1) = 1
      option%io_buffer = 'ERROR: Non-initialized cell volume.'
      call PrintMsgByRank(option)
    endif
    if (material_auxvars(ghosted_id)%porosity_base < 0.d0 .and. &
        flag(2) == 0) then
      flag(2) = 1
      option%io_buffer = 'ERROR: Non-initialized porosity.'
      call PrintMsgByRank(option)
    endif
    if (minval(material_auxvars(ghosted_id)%permeability) < 0.d0 .and. &
        flag(5) == 0) then
      option%io_buffer = 'ERROR: Non-initialized permeability.'
      call PrintMsgByRank(option)
      flag(5) = 1
    endif
  enddo

  error_found = error_found .or. (maxval(flag) > 0)
  call MPI_Allreduce(MPI_IN_PLACE,error_found,ONE_INTEGER_MPI,MPI_C_BOOL, &
                     MPI_LOR,option%mycomm,ierr);CHKERRQ(ierr)
  if (error_found) then
    option%io_buffer = 'Material property errors found in THCSetup.'
    call PrintErrMsg(option)
  endif

  ! initialize parameters
  cur_fluid_property => realization%fluid_properties
  do
    if (.not.associated(cur_fluid_property)) exit
    if (cur_fluid_property%phase_id == LIQUID_PHASE) then
      patch%aux%THC%thc_parameter%diffusion_coef = &
        cur_fluid_property%diffusion_coefficient
      patch%aux%THC%thc_parameter%diffusion_activation_energy = &
        cur_fluid_property%diffusion_activation_energy
      exit
    endif
    cur_fluid_property => cur_fluid_property%next
  enddo

  temp_int = 0
  if (thc_numerical_derivatives) then
    temp_int = option%nflowdof
  endif
  allocate(thc_auxvars(0:temp_int,grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    do idof = 0, temp_int
      call THCAuxVarInit(thc_auxvars(idof,ghosted_id),option)
    enddo
  enddo
  patch%aux%THC%auxvars => thc_auxvars
  patch%aux%THC%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    allocate(thc_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call THCAuxVarInit(thc_auxvars_bc(iconn),option)
    enddo
    patch%aux%THC%auxvars_bc => thc_auxvars_bc
  endif
  patch%aux%THC%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(thc_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call THCAuxVarInit(thc_auxvars_ss(iconn),option)
    enddo
    patch%aux%THC%auxvars_ss => thc_auxvars_ss
  endif
  patch%aux%THC%num_aux_ss = sum_connection

  list => realization%output_option%output_snap_variable_list
  call THCSetPlotVariables(realization,list)
  list => realization%output_option%output_obs_variable_list
  call THCSetPlotVariables(realization,list)

  if (Initialized(thc_debug_cell_id) .and. &
      option%comm%size > 1) then
    option%io_buffer = 'Cannot debug cells in parallel.'
    call PrintErrMsg(option)
  endif

  allocate(dof_is_active(option%nflowdof))
  dof_is_active = PETSC_TRUE
  call PatchCreateZeroArray(patch,dof_is_active, &
                            patch%aux%THC%matrix_zeroing,option)
  deallocate(dof_is_active)

  thc_ts_count = 0
  thc_ts_cut_count = 0
  thc_ni_count = 0

end subroutine THCSetup

! ************************************************************************** !

subroutine THCInitializeTimestep(realization)
  !
  ! Update data in module prior to time step
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class

  implicit none

  class(realization_subsurface_type) :: realization

  call THCUpdateFixedAccum(realization)
  thc_ni_count = 0

end subroutine THCInitializeTimestep

! ************************************************************************** !

subroutine THCUpdateSolution(realization)
  !
  ! Updates data in module after a successful time step
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class

  implicit none

  class(realization_subsurface_type) :: realization

  if (realization%option%compute_mass_balance_new) then
    call THCUpdateMassBalance(realization)
  endif

  thc_ts_count = thc_ts_count + 1
  thc_ts_cut_count = 0
  thc_ni_count = 0

end subroutine THCUpdateSolution

! ************************************************************************** !

subroutine THCTimeCut(realization)
  !
  ! Resets arrays for time step cut
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class

  implicit none

  class(realization_subsurface_type) :: realization

  thc_ts_cut_count = thc_ts_cut_count + 1

  call THCInitializeTimestep(realization)

end subroutine THCTimeCut

! ************************************************************************** !

subroutine THCComputeMassBalance(realization,mass_balance)
  !
  ! Initializes mass balance
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_module

  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(1)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(thc_auxvar_type), pointer :: thc_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt, parameter :: iphase = 1

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  thc_auxvars => patch%aux%THC%auxvars
  material_auxvars => patch%aux%Material%auxvars

  mass_balance = 0.d0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ! volume_phase = saturation*porosity*volume
    mass_balance(iphase) = mass_balance(iphase) + &
        thc_auxvars(ZERO_INTEGER,ghosted_id)%sat* &
        thc_auxvars(ZERO_INTEGER,ghosted_id)%effective_porosity* &
        material_auxvars(ghosted_id)%volume
  enddo

end subroutine THCComputeMassBalance

! ************************************************************************** !

subroutine THCZeroMassBalanceDelta(realization)
  !
  ! Zeros mass balance delta array
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  do iconn = 1, patch%aux%THC%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo
  do iconn = 1, patch%aux%THC%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine THCZeroMassBalanceDelta

! ************************************************************************** !

subroutine THCUpdateMassBalance(realization)
  !
  ! Updates mass balance
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  do iconn = 1, patch%aux%THC%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance(1,1) = &
      global_auxvars_bc(iconn)%mass_balance(1,1) + &
      global_auxvars_bc(iconn)%mass_balance_delta(1,1)*option%flow_dt
  enddo
  do iconn = 1, patch%aux%THC%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance(1,1) = &
      global_auxvars_ss(iconn)%mass_balance(1,1) + &
      global_auxvars_ss(iconn)%mass_balance_delta(1,1)*option%flow_dt
  enddo

end subroutine THCUpdateMassBalance

! ************************************************************************** !

subroutine THCUpdateAuxVars(realization)
  !
  ! Updates the auxiliary variables associated with the THC problem
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Material_Aux_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(thc_auxvar_type), pointer :: thc_auxvars(:,:)
  type(thc_auxvar_type), pointer :: thc_auxvars_bc(:)
  type(thc_auxvar_type), pointer :: thc_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, sum_connection, iconn, natural_id
  PetscInt :: dof_index
  PetscInt :: ghosted_start, ghosted_end, ghosted_offset
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscInt :: water_index, energy_index, solute_index
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  thc_auxvars => patch%aux%THC%auxvars
  thc_auxvars_bc => patch%aux%THC%auxvars_bc
  thc_auxvars_ss => patch%aux%THC%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArray(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ! THC_UPDATE_FOR_ACCUM indicates call from non-perturbation
    option%iflag = THC_UPDATE_FOR_ACCUM
    natural_id = grid%nG2A(ghosted_id)
    ghosted_end = ghosted_id * option%nflowdof
    ghosted_start = ghosted_end - option%nflowdof + 1
    if (grid%nG2L(ghosted_id) == 0) natural_id = -natural_id
    call THCAuxVarCompute(xx_loc_p(ghosted_start:ghosted_end), &
                              thc_auxvars(ZERO_INTEGER,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%cc_id(ghosted_id))%ptr, &
                              patch%aux%THC%thc_parameter, &
                              natural_id, &
                              PETSC_TRUE,option)
  enddo

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    water_index  = boundary_condition%flow_aux_mapping(THC_COND_WATER_INDEX)
    energy_index = boundary_condition%flow_aux_mapping(THC_COND_ENERGY_INDEX)
    solute_index = boundary_condition%flow_aux_mapping(THC_COND_SOLUTE_INDEX)
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      !geh: negate to indicate boundary connection, not actual cell
      natural_id = -grid%nG2A(ghosted_id)
      ghosted_offset = (ghosted_id-1)*option%nflowdof

      ! --- pressure (flow) DOF ----------------------------------------------
      select case(boundary_condition%flow_bc_type(water_index))
        case(DIRICHLET_BC, DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
             HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
          xxbc(thc_pressure_dof) = &
            boundary_condition%flow_aux_real_var(water_index,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC,UNIT_GRADIENT_BC, &
            SURFACE_ZERO_GRADHEIGHT)
          xxbc(thc_pressure_dof) = &
            xx_loc_p(ghosted_offset+thc_pressure_dof)
        case(PONDED_WATER_BC)
          xxbc(thc_pressure_dof) = &
            ! subtract 0.01 to ensure no inflow due to reference pressure
            (option%flow%reference_pressure-0.01d0) + &
            max(boundary_condition%flow_aux_real_var(water_index,iconn)* &
                dot_product(option%gravity, &
                            cur_connection_set%dist(1:3,iconn))* &
                thc_auxvars(ZERO_INTEGER,ghosted_id)%den_kg,0.d0)
        case default
          option%io_buffer = 'flow boundary itype not set up in &
            &THCUpdateAuxVars'
          call PrintErrMsg(option)
      end select

      ! --- temperature (energy) DOF -----------------------------------------
      select case(boundary_condition%flow_bc_type(energy_index))
        case(DIRICHLET_BC, DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
             HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
          xxbc(thc_temperature_dof) = &
            boundary_condition%flow_aux_real_var(energy_index,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC,CONVECTIVE_BC)
          ! prescribed-flux / zero-gradient / Robin: ghost carries interior T
          xxbc(thc_temperature_dof) = &
            xx_loc_p(ghosted_offset+thc_temperature_dof)
        case default
          option%io_buffer = 'energy boundary itype not set up in &
            &THCUpdateAuxVars'
          call PrintErrMsg(option)
      end select

      ! --- concentration (solute) DOF ---------------------------------------
      select case(boundary_condition%flow_bc_type(solute_index))
        case(DIRICHLET_BC, DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
             HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
          xxbc(thc_concentration_dof) = &
            boundary_condition%flow_aux_real_var(solute_index,iconn)
        case(ZERO_GRADIENT_BC)
          xxbc(thc_concentration_dof) = &
            xx_loc_p(ghosted_offset+thc_concentration_dof)
        case default
          option%io_buffer = 'solute boundary itype not set up in &
            &THCUpdateAuxVars'
          call PrintErrMsg(option)
      end select

      ! THC_UPDATE_FOR_BOUNDARY indicates call from non-perturbation
      option%iflag = THC_UPDATE_FOR_BOUNDARY
      call THCAuxVarCompute(xxbc,thc_auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                material_auxvars(ghosted_id), &
                                patch%characteristic_curves_array( &
                                  patch%cc_id(ghosted_id))%ptr, &
                                patch%aux%THC%thc_parameter, &
                                natural_id, &
                                PETSC_FALSE,option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      call THCAuxVarCopy(thc_auxvars(ZERO_INTEGER,ghosted_id), &
                             thc_auxvars_ss(sum_connection),option)
      call GlobalAuxVarCopy(global_auxvars(ghosted_id), &
                            global_auxvars_ss(sum_connection),option)
      ! override concentration from grid cells
      dof_index = source_sink%flow_aux_mapping(THC_COND_SOLUTE_INDEX)
      thc_auxvars_ss(sum_connection)%conc = &
        source_sink%flow_aux_real_var(dof_index,iconn)
    enddo
    source_sink => source_sink%next
  enddo

  call VecRestoreArray(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  patch%aux%THC%auxvars_up_to_date = PETSC_TRUE

end subroutine THCUpdateAuxVars

! ************************************************************************** !

subroutine THCUpdateFixedAccum(realization)
  !
  ! Updates the fixed portion of the accumulation term
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Material_Aux_module
  use Petsc_Utility_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(thc_auxvar_type), pointer :: thc_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, local_id, local_start, local_end, natural_id
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: accum_t_p(:)
  PetscReal :: Res(THC_NDOF)
  PetscReal :: Jdum(THC_NDOF,THC_NDOF)
  PetscReal :: dResdparam(THC_NDOF,THC_NDOF)
  PetscInt :: ndof

  PetscErrorCode :: ierr

  if (.not.thc_calc_accum) return

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  ndof = option%nflowdof

  thc_auxvars => patch%aux%THC%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter

  call VecGetArrayRead(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%flow_accum_t,accum_t_p,ierr);CHKERRQ(ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id * ndof
    local_start = local_end - ndof + 1
    natural_id = grid%nG2A(ghosted_id)
    ! THC_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = THC_UPDATE_FOR_FIXED_ACCUM
    call THCAuxVarCompute(xx_p(local_start:local_end), &
                              thc_auxvars(ZERO_INTEGER,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%cc_id(ghosted_id))%ptr, &
                              patch%aux%THC%thc_parameter, &
                              natural_id, &
                              PETSC_TRUE,option)
    call THCAccumulation(thc_auxvars(ZERO_INTEGER,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             option,Res, &
                             Jdum,dResdparam,PETSC_FALSE)
    call PetUtilVecSVBL(accum_t_p,local_id,Res,ndof,PETSC_TRUE)
  enddo

  call VecRestoreArrayRead(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%flow_accum_t,accum_t_p,ierr);CHKERRQ(ierr)

end subroutine THCUpdateFixedAccum

! ************************************************************************** !

subroutine THCResidual(snes,xx,r,A,realization,debug,ierr)
  !
  ! Computes the residual equation (and 3x3 Jacobian blocks via
  ! the thc_common kernels).
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Debug_module
  use Connection_module
  use Grid_module
  use Coupler_module
  use Material_Aux_module
  use Upwind_Direction_module
  use Petsc_Utility_module
  use Matrix_Zeroing_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  Mat :: A
  class(realization_subsurface_type) :: realization
  type(debug_type) :: debug
  PetscErrorCode :: ierr

  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(material_parameter_type), pointer :: material_parameter
  type(thc_parameter_type), pointer :: thc_parameter
  type(thc_auxvar_type), pointer :: thc_auxvars(:,:)
  type(thc_auxvar_type), pointer :: thc_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars_pert(:,:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: iconn
  PetscReal :: ss_flow_vol_flux
  PetscInt :: sum_connection
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: imat, imat_up, imat_dn

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_t_p(:), accum_t_p2(:)
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: xx_loc_p(:)

  PetscInt :: icc_up, icc_dn
  PetscReal :: Res(THC_NDOF)
  PetscReal :: Jup(THC_NDOF,THC_NDOF)
  PetscReal :: Jdn(THC_NDOF,THC_NDOF)
  PetscReal :: dResdparamup(THC_NDOF,THC_NDOF)
  PetscReal :: dResdparamdn(THC_NDOF,THC_NDOF)
  PetscReal :: dResdparam(THC_NDOF,THC_NDOF)
  PetscReal :: v_darcy

  PetscInt :: ndof
  PetscInt :: istart, iend

  ndof = realization%option%nflowdof

  dResdparamup = 0.d0
  dResdparamdn = 0.d0
  dResdparam = 0.d0

  discretization => realization%discretization
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  thc_auxvars => patch%aux%THC%auxvars
  thc_auxvars_bc => patch%aux%THC%auxvars_bc
  thc_parameter => patch%aux%THC%thc_parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  material_auxvars_pert => patch%aux%THC%material_auxvars_pert

  call MatZeroEntries(A,ierr);CHKERRQ(ierr)

  ! Communication -----------------------------------------
  ! must be called before THCUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call THCUpdateAuxVars(realization)

  ! override flags since they will soon be out of date
  patch%aux%THC%auxvars_up_to_date = PETSC_FALSE

  if (thc_numerical_derivatives) then
    ! Perturb aux vars (fills thc_auxvars(1:ndof,:) consumed by the
    ! *Derivative wrappers below).
    call VecGetArrayRead(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
    do ghosted_id = 1, grid%ngmax  ! For each local node do...
      if (patch%imat(ghosted_id) <= 0) cycle
      iend = ghosted_id*ndof
      istart = iend-ndof+1
      call THCAuxVarPerturb(xx_loc_p(istart:iend), &
                                thc_auxvars(:,ghosted_id), &
                                global_auxvars(ghosted_id), &
                                material_auxvars(ghosted_id), &
                                material_auxvars_pert(:,ghosted_id), &
                                patch%characteristic_curves_array( &
                                  patch%cc_id(ghosted_id))%ptr, &
                                patch%aux%THC%thc_parameter, &
                                grid%nG2A(ghosted_id),option)
    enddo
    call VecRestoreArrayRead(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  endif

  if (option%compute_mass_balance_new) then
    call THCZeroMassBalanceDelta(realization)
  endif

  option%iflag = 1
  ! now assign access pointer to local variables
  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)

  ! Accumulation terms ------------------------------------
  ! accumulation at t(k) (doesn't change during Newton iteration)
  if (thc_calc_accum) then
    call VecGetArrayRead(field%flow_accum_t,accum_t_p,ierr);CHKERRQ(ierr)
    r_p = -accum_t_p
    call VecRestoreArrayRead(field%flow_accum_t,accum_t_p,ierr);CHKERRQ(ierr)

    ! accumulation at t(k+1)
    call VecGetArray(field%flow_accum_tpdt,accum_t_p2,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      imat = patch%imat(ghosted_id)
      if (imat <= 0) cycle
      call THCAccumDerivative(thc_auxvars(:,ghosted_id), &
                                  global_auxvars(ghosted_id), &
                                  material_auxvars(ghosted_id), &
                                  option,Res,Jup,dResdparam)
      call PetUtilVecSVBL(r_p,local_id,Res,ndof,PETSC_FALSE)
      call PetUtilVecSVBL(accum_t_p2,local_id,Res,ndof,PETSC_TRUE)
      call PetUtilMatSVBL(A,ghosted_id,ghosted_id,Jup,ndof)
    enddo
    call VecRestoreArray(field%flow_accum_tpdt,accum_t_p2,ierr);CHKERRQ(ierr)
  else
    r_p = 0.d0
  endif

  if (thc_calc_flux) then
    ! Interior Flux Terms -----------------------------------
    connection_set_list => grid%internal_connection_set_list
    cur_connection_set => connection_set_list%first
    sum_connection = 0
    do
      if (.not.associated(cur_connection_set)) exit
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1

        ghosted_id_up = cur_connection_set%id_up(iconn)
        ghosted_id_dn = cur_connection_set%id_dn(iconn)

        local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
        local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

        imat_up = patch%imat(ghosted_id_up)
        imat_dn = patch%imat(ghosted_id_dn)
        if (imat_up <= 0 .or. imat_dn <= 0) cycle

        icc_up = patch%cc_id(ghosted_id_up)
        icc_dn = patch%cc_id(ghosted_id_dn)

        call THCFluxDerivative(thc_auxvars(:,ghosted_id_up), &
                         global_auxvars(ghosted_id_up), &
                         material_auxvars(ghosted_id_up), &
                         thc_auxvars(:,ghosted_id_dn), &
                         global_auxvars(ghosted_id_dn), &
                         material_auxvars(ghosted_id_dn), &
                         cur_connection_set%area(iconn), &
                         cur_connection_set%dist(:,iconn), &
                         thc_parameter,option,v_darcy, &
                         Res,Jup,Jdn,dResdparamup,dResdparamdn, &
                         PUCast(local_id_up == thc_debug_cell_id .or. &
                                local_id_dn == thc_debug_cell_id))
        ! store velocity for output / post-processing
        patch%internal_velocities(:,sum_connection) = v_darcy
        if (associated(patch%internal_flow_fluxes)) then
          patch%internal_flow_fluxes(1,sum_connection) = &
            Res(thc_pressure_dof)
        endif

        if (local_id_up > 0) then
          call PetUtilVecSVBL(r_p,local_id_up,Res,ndof,PETSC_FALSE)
          call PetUtilMatSVBL(A,ghosted_id_up,ghosted_id_up,Jup,ndof)
          call PetUtilMatSVBL(A,ghosted_id_up,ghosted_id_dn,Jdn,ndof)
        endif

        if (local_id_dn > 0) then
          Res = -Res
          call PetUtilVecSVBL(r_p,local_id_dn,Res,ndof,PETSC_FALSE)
          Jup = -Jup
          Jdn = -Jdn
          call PetUtilMatSVBL(A,ghosted_id_dn,ghosted_id_dn,Jdn,ndof)
          call PetUtilMatSVBL(A,ghosted_id_dn,ghosted_id_up,Jup,ndof)
        endif
      enddo

      cur_connection_set => cur_connection_set%next
    enddo
  endif

  if (thc_calc_bcflux) then
    ! Boundary Flux Terms -----------------------------------
    boundary_condition => patch%boundary_condition_list%first
    sum_connection = 0
    do
      if (.not.associated(boundary_condition)) exit

      cur_connection_set => boundary_condition%connection_set

      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1

        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)

        imat_dn = patch%imat(ghosted_id)
        if (imat_dn <= 0) cycle

        icc_dn = patch%cc_id(ghosted_id)

        call THCBCFluxDerivative(boundary_condition%flow_bc_type, &
                           boundary_condition%flow_aux_mapping, &
                           boundary_condition%flow_aux_real_var(:,iconn), &
                           thc_auxvars_bc(sum_connection), &
                           global_auxvars_bc(sum_connection), &
                           thc_auxvars(:,ghosted_id), &
                           global_auxvars(ghosted_id), &
                           material_auxvars(ghosted_id), &
                           cur_connection_set%area(iconn), &
                           cur_connection_set%dist(:,iconn), &
                           thc_parameter,option, &
                           v_darcy,Res,Jdn,dResdparamdn, &
                           PUCast(local_id == thc_debug_cell_id))
        patch%boundary_velocities(:,sum_connection) = v_darcy
        if (associated(patch%boundary_flow_fluxes)) then
          patch%boundary_flow_fluxes(1,sum_connection) = &
            Res(thc_pressure_dof)
        endif
        if (option%compute_mass_balance_new) then
          ! contribution to boundary
          global_auxvars_bc(sum_connection)%mass_balance_delta(1,1) = &
            global_auxvars_bc(sum_connection)%mass_balance_delta(1,1) - &
            Res(thc_pressure_dof)
        endif
        Res = -Res
        call PetUtilVecSVBL(r_p,local_id,Res,ndof,PETSC_FALSE)
        Jdn = -Jdn
        call PetUtilMatSVBL(A,ghosted_id,ghosted_id,Jdn,ndof)
      enddo
      boundary_condition => boundary_condition%next
    enddo
  endif

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit

    cur_connection_set => source_sink%connection_set

    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      call THCSrcSinkDerivative(option, &
                          source_sink%flow_aux_real_var(:,iconn), &
                          source_sink%flow_aux_mapping, &
                          source_sink%flow_bc_type, &
                          source_sink%flow_condition%rate%aux_real(:), &
                          thc_auxvars(:,ghosted_id), &
                          global_auxvars(ghosted_id), &
                          material_auxvars(ghosted_id), &
                          ss_flow_vol_flux,Res,Jdn, &
                          dResdparamdn)
      if (associated(patch%ss_flow_vol_fluxes)) then
        patch%ss_flow_vol_fluxes(:,sum_connection) = ss_flow_vol_flux
      endif
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(:,sum_connection) = Res(1)
      endif
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_ss(sum_connection)%mass_balance_delta(1,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(1,1) - &
          Res(thc_pressure_dof)
      endif
      Res = -Res
      call PetUtilVecSVBL(r_p,local_id,Res,ndof,PETSC_FALSE)
      Jdn = -Jdn
      call PetUtilMatSVBL(A,ghosted_id,ghosted_id,Jdn,ndof)
    enddo
    source_sink => source_sink%next
  enddo

  call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)

  call MatrixZeroingZeroVecEntries(patch%aux%THC%matrix_zeroing,r)

  if (thc_simult_function_evals) then

    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    ! zero out inactive cells
    call MatrixZeroingZeroMatEntries(patch%aux%THC%matrix_zeroing,A)

    if (debug%matview_Matrix) then
      call DebugMatView(debug,A,'THCFjacobian','', &
                        thc_ts_count,thc_ts_cut_count, &
                        thc_ni_count,option)
    endif
  endif

  ! Mass Transfer
  if (.not.PetscObjectIsNull(field%flow_mass_transfer)) then
    ! scale by -1.d0 for contribution to residual.  A negative contribution
    ! indicates mass being added to system.
    call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
    call VecGetArray(field%flow_mass_transfer,vec_p,ierr);CHKERRQ(ierr)
    ! geh: leave in expanded do loop form instead of VecAXPY for flexibility
    !      in the future
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      imat = patch%imat(ghosted_id)
      if (imat <= 0) cycle
      r_p(local_id) = r_p(local_id) - vec_p(local_id)
    enddo
    call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)
    call VecRestoreArray(field%flow_mass_transfer,vec_p, &
                            ierr);CHKERRQ(ierr)
  endif

  if (debug%vecview_residual) then
    call DebugVecView(debug,r,'THCFresidual','', &
                      thc_ts_count,thc_ts_cut_count, &
                      thc_ni_count,option)
  endif
  if (debug%vecview_solution) then
    call DebugVecView(debug,xx,'THCFxx','', &
                      thc_ts_count,thc_ts_cut_count, &
                      thc_ni_count,option)
  endif

end subroutine THCResidual

! ************************************************************************** !

subroutine THCSetPlotVariables(realization,list)
  !
  ! Adds variables to be printed to list
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class
  use Output_Aux_module
  use Variables_module

  implicit none

  class(realization_subsurface_type) :: realization
  type(output_variable_list_type), pointer :: list

  character(len=MAXWORDLENGTH) :: name, units

  if (associated(list%first)) then
    return
  endif

  if (list%flow_vars) then

    name = 'Liquid Pressure'
    units = 'Pa'
    call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                                 LIQUID_PRESSURE)

    name = 'Temperature'
    units = 'C'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                 TEMPERATURE)

    name = 'Solute Concentration'
    units = 'M'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                 SOLUTE_CONCENTRATION)

    name = 'Liquid Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                 LIQUID_SATURATION)

    name = 'Liquid Density'
    units = 'kg/m^3'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                 LIQUID_DENSITY)

    name = 'Liquid Viscosity'
    units = 'Pa-s'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                 LIQUID_VISCOSITY)

    name = 'Effective Porosity'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                 POROSITY)

    name = 'Thermal Conductivity'
    units = 'W/m-C'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                 THERMAL_CONDUCTIVITY)

    name = 'Capillary Pressure'
    units = 'Pa'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                 CAPILLARY_PRESSURE)

    name = 'Liquid Relative Permeability'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                 LIQUID_RELATIVE_PERMEABILITY)
  endif

end subroutine THCSetPlotVariables

! ************************************************************************** !

subroutine THCMapBCAuxVarsToGlobal(realization)
  !
  ! Maps variables in thc auxvar to global equivalent (for ERT coupling).
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Coupler_module
  use Connection_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(thc_auxvar_type), pointer :: thc_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)

  PetscInt :: sum_connection, iconn

  option => realization%option
  patch => realization%patch

  if (option%ntrandof == 0) return ! no need to update

  thc_auxvars_bc => patch%aux%THC%auxvars_bc
  global_auxvars_bc => patch%aux%Global%auxvars_bc

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      global_auxvars_bc(sum_connection)%sat = &
        thc_auxvars_bc(sum_connection)%sat
    enddo
    boundary_condition => boundary_condition%next
  enddo

end subroutine THCMapBCAuxVarsToGlobal

! ************************************************************************** !

subroutine THCDestroy(realization)
  !
  ! Deallocates variables associated with THC
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class
  use Option_module

  implicit none

  class(realization_subsurface_type) :: realization

  ! place anything that needs to be freed here.
  ! auxvars are deallocated in auxiliary.F90.

end subroutine THCDestroy

end module THC_module
