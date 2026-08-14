module PM_THC_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class

  use PFLOTRAN_Constants_module
  use THC_Aux_module

  implicit none

  private

  PetscInt, parameter :: MAX_CHANGE_LIQ_PRES_NI = 1
  PetscInt, parameter :: MAX_CHANGE_TEMP_NI = 2
  PetscInt, parameter :: MAX_CHANGE_CONC_NI = 3
  PetscInt, parameter :: MAX_RES_LIQ_EQ = 4
  PetscInt, parameter :: MAX_RES_ENERGY_EQ = 5
  PetscInt, parameter :: MAX_RES_SOL_EQ = 6

  type, public, extends(pm_subsurface_flow_type) :: pm_thc_type
    PetscInt, pointer :: max_change_ivar(:)
    PetscReal :: max_allow_liq_pres_change_ni
    PetscReal :: liq_pres_change_ts_governor
    PetscReal :: liq_sat_change_ts_governor
    PetscReal :: temp_change_ts_governor
    PetscReal :: conc_change_ts_governor
    PetscInt :: convergence_flags(MAX_RES_SOL_EQ)
    PetscReal :: convergence_reals(MAX_RES_SOL_EQ)
  contains
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMTHCReadSimOptionsBlock
    procedure, public :: ReadTSBlock => PMTHCReadTSSelectCase
    procedure, public :: ReadNewtonBlock => PMTHCReadNewtonSelectCase
    procedure, public :: Setup => PMTHCSetup
    procedure, public :: InitializeRun => PMTHCInitializeRun
    procedure, public :: InitializeTimestep => PMTHCInitializeTimestep
    procedure, public :: Residual => PMTHCResidual
    procedure, public :: Jacobian => PMTHCJacobian
    procedure, public :: UpdateTimestep => PMTHCUpdateTimestep
    procedure, public :: FinalizeTimestep => PMTHCFinalizeTimestep
    procedure, public :: PreSolve => PMTHCPreSolve
    procedure, public :: PostSolve => PMTHCPostSolve
    procedure, public :: CheckUpdatePost => PMTHCCheckUpdatePost
    procedure, public :: CheckConvergence => PMTHCCheckConvergence
    procedure, public :: TimeCut => PMTHCTimeCut
    procedure, public :: UpdateSolution => PMTHCUpdateSolution
    procedure, public :: UpdateAuxVars => PMTHCUpdateAuxVars
    procedure, public :: MaxChange => PMTHCMaxChange
    procedure, public :: ComputeMassBalance => PMTHCComputeMassBalance
    procedure, public :: InputRecord => PMTHCInputRecord
    procedure, public :: CheckpointBinary => PMTHCCheckpointBinary
    procedure, public :: RestartBinary => PMTHCRestartBinary
    procedure, public :: Destroy => PMTHCDestroy
  end type pm_thc_type

  public :: PMTHCCreate, &
            PMTHCInitObject, &
            PMTHCInitializeRun, &
            PMTHCFinalizeTimestep, &
            PMTHCDestroy

contains

! ************************************************************************** !

function PMTHCCreate()
  !
  ! Creates THC process model shell
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  implicit none

  class(pm_thc_type), pointer :: PMTHCCreate

  class(pm_thc_type), pointer :: thc_pm

  allocate(thc_pm)
  call PMTHCInitObject(thc_pm)

  PMTHCCreate => thc_pm

end function PMTHCCreate

! ************************************************************************** !

subroutine PMTHCInitObject(this)
  !
  ! Initializes THC process model shell
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Option_module
  use String_module
  use Variables_module

  implicit none

  class(pm_thc_type) :: this

  call PMSubsurfaceFlowInit(this)
  this%name = 'THC Flow'
  this%header = 'THC FLOW'

  nullify(this%max_change_ivar)

  ! set to UNINITIALIZED_DOUBLE and report error below if set from input
  this%pressure_change_governor = UNINITIALIZED_DOUBLE
  this%temperature_change_governor = UNINITIALIZED_DOUBLE
  this%saturation_change_governor = UNINITIALIZED_DOUBLE
  this%xmol_change_governor = UNINITIALIZED_DOUBLE

  this%check_post_convergence = PETSC_TRUE

  this%max_allow_liq_pres_change_ni = UNINITIALIZED_DOUBLE
  this%liq_pres_change_ts_governor = 5.d5    ! [Pa]
  this%liq_sat_change_ts_governor = 1.d0
  this%temp_change_ts_governor = 5.d0        ! [C]
  this%conc_change_ts_governor = 1.d0        ! [mol/L]

  this%convergence_flags = 0
  this%convergence_reals = 0.d0

  thc_debug_cell_id = UNINITIALIZED_INTEGER

end subroutine PMTHCInitObject

! ************************************************************************** !

subroutine PMTHCReadSimOptionsBlock(this,input)
  !
  ! Read THC options input block
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use THC_module
  use THC_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module
  use Utility_module

  implicit none

  class(pm_thc_type) :: this
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: word
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option

  error_string = 'THC Options'

  input%ierr = INPUT_ERROR_NONE
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call PMSubsurfFlowReadSimOptionsSC(this,input,keyword,found, &
                                       error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case('NO_ACCUMULATION')
        thc_calc_accum = PETSC_FALSE
      case('NO_FLUX')
        thc_calc_flux = PETSC_FALSE
      case('NO_BCFLUX')
        thc_calc_bcflux = PETSC_FALSE
      case('TENSORIAL_RELATIVE_PERMEABILITY')
        thc_tensorial_rel_perm = PETSC_TRUE
      case('SOLID_DENSITY','SOLID_GRAIN_DENSITY','ROCK_DENSITY')
        call InputReadDouble(input,option,thc_density_solid)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,thc_density_solid,'kg/m^3', &
                                      trim(error_string)//','//keyword,option)
      case('SPECIFIC_HEAT_LIQUID','HEAT_CAPACITY_LIQUID','LIQUID_HEAT_CAPACITY')
        call InputReadDouble(input,option,thc_specific_heat_liquid)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,thc_specific_heat_liquid, &
                                      'J/kg-C',trim(error_string)//','// &
                                      keyword,option)
      case('SPECIFIC_HEAT_SOLID','HEAT_CAPACITY_SOLID','SOLID_HEAT_CAPACITY')
        call InputReadDouble(input,option,thc_specific_heat_solid)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,thc_specific_heat_solid, &
                                      'J/kg-C',trim(error_string)//','// &
                                      keyword,option)
      case('THERMAL_CONDUCTIVITY_SOLID')
        call InputReadDouble(input,option,thc_kappa_solid)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,thc_kappa_solid,'W/m-C', &
                                      trim(error_string)//','//keyword,option)
      case('THERMAL_CONDUCTIVITY_DRY')
        call InputReadDouble(input,option,thc_kappa_dry)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,thc_kappa_dry,'W/m-C', &
                                      trim(error_string)//','//keyword,option)
      case('THERMAL_CONDUCTIVITY_WET')
        call InputReadDouble(input,option,thc_kappa_wet)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,thc_kappa_wet,'W/m-C', &
                                      trim(error_string)//','//keyword,option)
      case('ENERGY_FORMULATION')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,keyword,error_string)
        call StringToUpper(word)
        select case(trim(word))
          case('RHO_CP_T','RHOCPT')
            thc_energy_mode = THC_ENERGY_RHO_CP_T
          case('FULL_EOS','FULLEOS','EOS')
            thc_energy_mode = THC_ENERGY_FULL_EOS
          case default
            call InputKeywordUnrecognized(input,word, &
                   'THC Mode,ENERGY_FORMULATION',option)
        end select
      case('ADVECTIVE_DENSITY')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,keyword,error_string)
        call StringToUpper(word)
        select case(trim(word))
          case('UPWIND')
            thc_advective_density_mode = THC_ADVECTIVE_DENSITY_UPWIND
          case('TH_COMPATIBLE')
            thc_advective_density_mode = &
              THC_ADVECTIVE_DENSITY_TH_COMPATIBLE
          case default
            call InputKeywordUnrecognized(input,word, &
                   'THC Mode,ADVECTIVE_DENSITY',option)
        end select
      case('DEBUG_CELL_ID')
        call InputReadInt(input,option,thc_debug_cell_id)
        call InputErrorMsg(input,option,keyword,error_string)
      case default
        call InputKeywordUnrecognized(input,keyword,'THC Mode',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine PMTHCReadSimOptionsBlock

! ************************************************************************** !

subroutine PMTHCReadTSSelectCase(this,input,keyword,found, &
                                     error_string,option)
  !
  ! Read timestepper settings specific to the THC process model
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  class(pm_thc_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  found = PETSC_TRUE
  call PMSubsurfaceFlowReadTSSelectCase(this,input,keyword,found, &
                                        error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case('LIQ_SAT_CHANGE_TS_GOVERNOR')
      call InputReadDouble(input,option,this%liq_sat_change_ts_governor)
      call InputErrorMsg(input,option,keyword,error_string)
    case('LIQ_PRES_CHANGE_TS_GOVERNOR')
      call InputReadDouble(input,option,this%liq_pres_change_ts_governor)
      call InputErrorMsg(input,option,keyword,error_string)
      ! units conversion since it is absolute
      call InputReadAndConvertUnits(input,this%liq_pres_change_ts_governor, &
                                    'Pa',keyword,option)
    case('TEMP_CHANGE_TS_GOVERNOR')
      call InputReadDouble(input,option,this%temp_change_ts_governor)
      call InputErrorMsg(input,option,keyword,error_string)
    case('CONC_CHANGE_TS_GOVERNOR')
      call InputReadDouble(input,option,this%conc_change_ts_governor)
      call InputErrorMsg(input,option,keyword,error_string)
    case default
      found = PETSC_FALSE
  end select

end subroutine PMTHCReadTSSelectCase

! ************************************************************************** !

subroutine PMTHCReadNewtonSelectCase(this,input,keyword,found, &
                                         error_string,option)
  !
  ! Reads input file parameters associated with the THC process model
  ! Newton solver convergence
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  use THC_Aux_module

  implicit none

  class(pm_thc_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  error_string = 'THC Newton Solver'

  select case(trim(keyword))
    case('ITOL_UPDATE')
      option%io_buffer = 'ITOL_UPDATE not supported with THC. Please &
        &use MAX_ALLOW_LIQ_PRES_CHANGE_NI.'
      call PrintErrMsg(option)
  end select

  found = PETSC_FALSE
  call PMSubsurfaceFlowReadNewtonSelectCase(this,input,keyword,found, &
                                            error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case('REL_PERTURBATION')
      call InputReadDouble(input,option,thc_rel_pert)
      call InputErrorMsg(input,option,keyword,error_string)
      ! no units conversion since it is relative
    case('MIN_LIQ_PRESSURE_PERTURBATION')
      call InputReadDouble(input,option,thc_pres_min_pert)
      call InputErrorMsg(input,option,keyword,error_string)
      call InputReadAndConvertUnits(input,thc_pres_min_pert, &
                                    'Pa',keyword,option)
    case('MIN_TEMPERATURE_PERTURBATION')
      call InputReadDouble(input,option,thc_temp_min_pert)
      call InputErrorMsg(input,option,keyword,error_string)
    case('MIN_CONCENTRATION_PERTURBATION')
      call InputReadDouble(input,option,thc_conc_min_pert)
      call InputErrorMsg(input,option,keyword,error_string)
    case('MAX_ALLOW_LIQ_PRES_CHANGE_NI')
      call InputReadDouble(input,option,this%max_allow_liq_pres_change_ni)
      call InputErrorMsg(input,option,keyword,error_string)
      ! units conversion since it is absolute
      call InputReadAndConvertUnits(input,this%max_allow_liq_pres_change_ni, &
                                    'Pa',keyword,option)
    case default
      found = PETSC_FALSE

  end select

end subroutine PMTHCReadNewtonSelectCase

! ************************************************************************** !

subroutine PMTHCSetup(this)
  !
  ! Sets up auxvars and parameters
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  use THC_module

  implicit none

  class(pm_thc_type) :: this

  call this%SetRealization()
  call THCSetup(this%realization)
  call PMSubsurfaceFlowSetup(this)

end subroutine PMTHCSetup

! ************************************************************************** !

recursive subroutine PMTHCInitializeRun(this)
  !
  ! Initializes the THC mode run.
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  use Realization_Base_class
  use Patch_module
  use Field_module
  use Material_Aux_module
  use Option_module
  use Variables_module

  implicit none

  class(pm_thc_type) :: this

  PetscInt :: i
  PetscErrorCode :: ierr
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option

  patch => this%realization%patch
  field => this%realization%field
  option => this%option

  if (this%steady_state) thc_calc_accum = PETSC_FALSE

  ! THC always solves all three equations; the max change variables
  ! are pressure, saturation, temperature, and concentration.
  allocate(this%max_change_ivar(4))
  this%max_change_ivar(1) = LIQUID_PRESSURE
  this%max_change_ivar(2) = LIQUID_SATURATION
  this%max_change_ivar(3) = TEMPERATURE
  this%max_change_ivar(4) = SOLUTE_CONCENTRATION

  ! need to allocate vectors for max change
  i = size(this%max_change_ivar)
  call VecDuplicateVecs(field%work,i,field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, size(field%max_change_vecs)
    call RealizationGetVariable(this%realization,field%max_change_vecs(i), &
                                this%max_change_ivar(i),ZERO_INTEGER)
  enddo

  ! call parent implementation
  call PMSubsurfaceFlowInitializeRun(this)

  if (soil_compressibility_index == 0 .and. associated(option%inversion)) then
    option%io_buffer = 'Soil compressibility must be employed for THC &
      &when used for inversion.'
    call PrintErrMsg(option)
  endif

end subroutine PMTHCInitializeRun

! ************************************************************************** !

subroutine PMTHCInitializeTimestep(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use THC_module, only : THCInitializeTimestep

  implicit none

  class(pm_thc_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
  call THCInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)

  this%convergence_flags = 0
  this%convergence_reals = 0.d0

end subroutine PMTHCInitializeTimestep

! ************************************************************************** !

subroutine PMTHCFinalizeTimestep(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  implicit none

  class(pm_thc_type) :: this

  call PMSubsurfaceFlowFinalizeTimestep(this)

end subroutine PMTHCFinalizeTimestep

! ************************************************************************** !

subroutine PMTHCPreSolve(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  implicit none

  class(pm_thc_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMTHCPreSolve

! ************************************************************************** !

subroutine PMTHCPostSolve(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  implicit none

  class(pm_thc_type) :: this

end subroutine PMTHCPostSolve

! ************************************************************************** !

subroutine PMTHCUpdateTimestep(this,update_dt, &
                                   dt,dt_min,dt_max,iacceleration, &
                                   num_newton_iterations,tfac, &
                                   time_step_max_growth_factor)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Option_module
  use Utility_module, only : Equal

  implicit none

  class(pm_thc_type) :: this
  PetscBool :: update_dt
  PetscReal :: dt
  PetscReal :: dt_min ! DO NOT USE (see comment below)
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: sat_ratio, pres_ratio, temp_ratio, conc_ratio
  PetscReal :: min_ratio
  PetscReal :: dt_prev

  if (update_dt .and. iacceleration /= 0) then
    dt_prev = dt
    ! calculate the time step ramping factor
    sat_ratio = (2.d0*this%liq_sat_change_ts_governor)/ &
                (this%liq_sat_change_ts_governor+this%max_saturation_change)
    pres_ratio = (2.d0*this%liq_pres_change_ts_governor)/ &
                (this%liq_pres_change_ts_governor+this%max_pressure_change)
    temp_ratio = (2.d0*this%temp_change_ts_governor)/ &
                (this%temp_change_ts_governor+this%max_temperature_change)
    conc_ratio = (2.d0*this%conc_change_ts_governor)/ &
                (this%conc_change_ts_governor+this%max_xmol_change)
    min_ratio = min(sat_ratio,pres_ratio,temp_ratio,conc_ratio)
    ! pick minimum time step from calc'd ramping factor or maximum ramping factor
    dt = min(min_ratio*dt,time_step_max_growth_factor*dt)
    ! make sure time step is within bounds given in the input deck
    dt = min(dt,dt_max)
    if (this%logging_verbosity > 0) then
      if (Equal(dt,dt_max)) then
        string = 'maximum time step size'
      else if (min_ratio > time_step_max_growth_factor) then
        string = 'maximum time step growth factor'
      else if (Equal(min_ratio,sat_ratio)) then
        string = 'liquid saturation governor'
      else if (Equal(min_ratio,pres_ratio)) then
        string = 'liquid pressure governor'
      else if (Equal(min_ratio,temp_ratio)) then
        string = 'temperature governor'
      else
        string = 'concentration governor'
      endif
      string = 'TS update: ' // trim(string)
      call PrintMsg(this%option,string)
    endif
  endif

  if (Initialized(this%cfl_governor)) then
    call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt,dt_max)
  endif

end subroutine PMTHCUpdateTimestep

! ************************************************************************** !

subroutine PMTHCResidual(this,snes,xx,r,ierr)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use THC_module, only : THCResidual
  use Debug_module
  use Grid_module

  implicit none

  class(pm_thc_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  Mat :: M

  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  ! calculate residual
  if (thc_simult_function_evals) then
    call SNESGetJacobian(snes,M,PETSC_NULL_MAT,PETSC_NULL_FUNCTION, &
                         PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
    call THCResidual(snes,xx,r,M,this%realization,this%debug,ierr)
  else
    call THCResidual(snes,xx,r,PETSC_NULL_MAT,this%realization, &
                         this%debug,ierr)
  endif

  if (this%debug%vecview_residual) then
    call DebugVecView(this%debug,r,'THCFresidual.out',this%option)
  endif
  if (this%debug%vecview_solution) then
    call DebugVecView(this%debug,xx,'THCFxx.out',this%option)
  endif

end subroutine PMTHCResidual

! ************************************************************************** !

subroutine PMTHCJacobian(this,snes,xx,A,B,ierr)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Debug_module
  use Option_module

  implicit none

  class(pm_thc_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  PetscReal :: norm

  ! the Jacobian was already calculated in PMTHCResidual

  if (this%debug%matview_Matrix) then
    call DebugMatView(this%debug,A,'THCFjacobian','', &
                      thc_ts_count,thc_ts_cut_count, &
                      thc_ni_count,this%option)
  endif

  if (this%debug%norm_Matrix) then
    call MatNorm(A,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("1 norm: ",es11.4)') norm
    call PrintMsg(this%option)
    call MatNorm(A,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("2 norm: ",es11.4)') norm
    call PrintMsg(this%option)
    call MatNorm(A,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("inf norm: ",es11.4)') norm
    call PrintMsg(this%option)
  endif

  thc_ni_count = thc_ni_count + 1

end subroutine PMTHCJacobian

! ************************************************************************** !

subroutine PMTHCCheckUpdatePost(this,snes,X0,dX,X1,dX_changed, &
                                    X1_changed,ierr)
  !
  ! Tracks the maximum absolute change in each primary variable (P, T, C)
  ! over a Newton iteration so that PMTHCCheckConvergence can enforce a
  ! solution-change convergence criterion.
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Material_Aux_module
  use THC_Aux_module

  implicit none

  class(pm_thc_type) :: this
  SNES :: snes
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch

  PetscReal, pointer :: dX_p(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscReal :: tempreal

  PetscReal :: max_abs_pressure_change_NI
  PetscInt :: max_abs_pressure_change_NI_cell
  PetscReal :: max_abs_temp_change_NI
  PetscInt :: max_abs_temp_change_NI_cell
  PetscReal :: max_abs_conc_change_NI
  PetscInt :: max_abs_conc_change_NI_cell

  grid => this%realization%patch%grid
  option => this%realization%option
  patch => this%realization%patch

  ! If these are changed from true, we must add a global reduction on both
  ! variables to ensure that their values match across all processes. Otherwise
  ! PETSc will throw an error in debug mode or ignore the error in optimized.
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE

  ! reset the solution-change flags each Newton iteration (THC has no
  ! pre-check, which is where ZFLOW performs this reset)
  this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI) = 0
  this%convergence_flags(MAX_CHANGE_TEMP_NI) = 0
  this%convergence_flags(MAX_CHANGE_CONC_NI) = 0

  call VecGetArray(dX,dX_p,ierr);CHKERRQ(ierr)
  max_abs_pressure_change_NI = 0.d0
  max_abs_pressure_change_NI_cell = 0
  max_abs_temp_change_NI = 0.d0
  max_abs_temp_change_NI_cell = 0
  max_abs_conc_change_NI = 0.d0
  max_abs_conc_change_NI_cell = 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle

    offset = (local_id-1)*option%nflowdof
    tempreal = dabs(dX_p(offset+thc_pressure_dof))
    if (tempreal > dabs(max_abs_pressure_change_NI)) then
      max_abs_pressure_change_NI_cell = grid%nG2A(ghosted_id)
      max_abs_pressure_change_NI = tempreal
    endif
    tempreal = dabs(dX_p(offset+thc_temperature_dof))
    if (tempreal > dabs(max_abs_temp_change_NI)) then
      max_abs_temp_change_NI_cell = grid%nG2A(ghosted_id)
      max_abs_temp_change_NI = tempreal
    endif
    tempreal = dabs(dX_p(offset+thc_concentration_dof))
    if (tempreal > dabs(max_abs_conc_change_NI)) then
      max_abs_conc_change_NI_cell = grid%nG2A(ghosted_id)
      max_abs_conc_change_NI = tempreal
    endif
  enddo

  ! flag non-convergence in liquid pressure when a tolerance has been set
  if (Initialized(this%max_allow_liq_pres_change_ni) .and. &
      max_abs_pressure_change_NI > this%max_allow_liq_pres_change_ni) then
    this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI) = &
      max_abs_pressure_change_NI_cell
  endif

  ! the following reals are for REPORTING purposes only
  this%convergence_reals(MAX_CHANGE_LIQ_PRES_NI) = max_abs_pressure_change_NI
  this%convergence_reals(MAX_CHANGE_TEMP_NI) = max_abs_temp_change_NI
  this%convergence_reals(MAX_CHANGE_CONC_NI) = max_abs_conc_change_NI

  call VecRestoreArray(dX,dX_p,ierr);CHKERRQ(ierr)

end subroutine PMTHCCheckUpdatePost

! ************************************************************************** !

subroutine PMTHCCheckConvergence(this,snes,it,xnorm,unorm, &
                                     fnorm,reason,ierr)
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Material_Aux_module
  use THC_Aux_module
  use Convergence_module

  implicit none

  class(pm_thc_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

  Vec :: residual_vec
  PetscReal, pointer :: r_p(:)
  character(len=10) :: reason_string

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(thc_auxvar_type), pointer :: thc_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: converged_flag

  PetscReal :: max_abs_res_liq_
  PetscInt :: max_abs_res_liq_cell
  PetscReal :: max_abs_res_energy_
  PetscInt :: max_abs_res_energy_cell
  PetscReal :: max_abs_res_sol_
  PetscInt :: max_abs_res_sol_cell
  PetscInt :: offset
  PetscMPIInt :: int_mpi

  PetscReal :: residual
  PetscReal :: tempreal

  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch
  thc_auxvars => patch%aux%THC%auxvars
  material_auxvars => patch%aux%Material%auxvars

  max_abs_res_liq_ = 0.d0
  max_abs_res_liq_cell = 0
  max_abs_res_energy_ = 0.d0
  max_abs_res_energy_cell = 0
  max_abs_res_sol_ = 0.d0
  max_abs_res_sol_cell = 0

  residual_vec = field%flow_r
  ! check residual terms
  call VecGetArrayRead(residual_vec,r_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%nflowdof
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    residual = r_p(offset+thc_pressure_dof)
    tempreal = dabs(residual)
    if (tempreal > max_abs_res_liq_) then
      max_abs_res_liq_ = tempreal
      max_abs_res_liq_cell = grid%nG2A(ghosted_id)
    endif
    residual = r_p(offset+thc_temperature_dof)
    tempreal = dabs(residual)
    if (tempreal > max_abs_res_energy_) then
      max_abs_res_energy_ = tempreal
      max_abs_res_energy_cell = grid%nG2A(ghosted_id)
    endif
    residual = r_p(offset+thc_concentration_dof)
    tempreal = dabs(residual)
    if (tempreal > max_abs_res_sol_) then
      max_abs_res_sol_ = tempreal
      max_abs_res_sol_cell = grid%nG2A(ghosted_id)
    endif
  enddo

  ! the following flags are for REPORTING purposes only
  this%convergence_flags(MAX_RES_LIQ_EQ) = max_abs_res_liq_cell
  this%convergence_reals(MAX_RES_LIQ_EQ) = max_abs_res_liq_
  this%convergence_flags(MAX_RES_ENERGY_EQ) = max_abs_res_energy_cell
  this%convergence_reals(MAX_RES_ENERGY_EQ) = max_abs_res_energy_
  this%convergence_flags(MAX_RES_SOL_EQ) = max_abs_res_sol_cell
  this%convergence_reals(MAX_RES_SOL_EQ) = max_abs_res_sol_

  int_mpi = size(this%convergence_flags)
  call MPI_Allreduce(MPI_IN_PLACE,this%convergence_flags,int_mpi,MPIU_INTEGER, &
                     MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
  int_mpi = size(this%convergence_reals)
  call MPI_Allreduce(MPI_IN_PLACE,this%convergence_reals,int_mpi, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                     ierr);CHKERRQ(ierr)

  ! these conditionals cannot change order
  reason_string = '----| '
  converged_flag = CONVERGENCE_CONVERGED
  if (this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI) > 0) then
    reason_string(1:1) = 'P'
    converged_flag = CONVERGENCE_KEEP_ITERATING
  endif
  if (this%convergence_flags(MAX_CHANGE_TEMP_NI) > 0) then
    reason_string(2:2) = 'T'
    converged_flag = CONVERGENCE_KEEP_ITERATING
  endif
  if (this%convergence_flags(MAX_CHANGE_CONC_NI) > 0) then
    reason_string(3:3) = 'C'
    converged_flag = CONVERGENCE_KEEP_ITERATING
  endif

  if (Initialized(this%max_allow_liq_pres_change_ni)) then
    option%convergence = converged_flag
  else
    ! forced standard 2 norms
    option%convergence = CONVERGENCE_OFF
  endif

  call VecRestoreArrayRead(residual_vec,r_p,ierr);CHKERRQ(ierr)

  call PMSubsurfaceFlowCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                        reason,ierr)

end subroutine PMTHCCheckConvergence

! ************************************************************************** !

subroutine PMTHCTimeCut(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use THC_module, only : THCTimeCut

  implicit none

  class(pm_thc_type) :: this

  call PMSubsurfaceFlowTimeCut(this)
  call THCTimeCut(this%realization)

  this%convergence_flags = 0
  this%convergence_reals = 0.d0

end subroutine PMTHCTimeCut

! ************************************************************************** !

subroutine PMTHCUpdateSolution(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use THC_module, only : THCUpdateSolution, &
                             THCMapBCAuxVarsToGlobal

  implicit none

  class(pm_thc_type) :: this

  call PMSubsurfaceFlowUpdateSolution(this)
  call THCUpdateSolution(this%realization)
  call THCMapBCAuxVarsToGlobal(this%realization)

end subroutine PMTHCUpdateSolution

! ************************************************************************** !

subroutine PMTHCUpdateAuxVars(this)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  use THC_module, only : THCUpdateAuxVars

  implicit none

  class(pm_thc_type) :: this

  call THCUpdateAuxVars(this%realization)

end subroutine PMTHCUpdateAuxVars

! ************************************************************************** !

subroutine PMTHCMaxChange(this)
  !
  ! Computes the maximum change in the solution for the timestep governors.
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use String_module
  use THC_Aux_module

  implicit none

  class(pm_thc_type) :: this

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_old_ptr(:), vec_new_ptr(:)
  PetscReal, allocatable :: max_change_global(:)
  PetscReal :: max_change, change
  PetscInt :: i, j
  PetscErrorCode :: ierr

  realization => this%realization
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  allocate(max_change_global(size(this%max_change_ivar)))
  max_change_global = 0.d0

  do i = 1, size(this%max_change_ivar)
    call RealizationGetVariable(realization,field%work, &
                                this%max_change_ivar(i),ZERO_INTEGER)
    ! yes, we could use VecWAXPY and a norm here, but we need the ability
    ! to customize
    call VecGetArray(field%work,vec_new_ptr,ierr);CHKERRQ(ierr)
    call VecGetArray(field%max_change_vecs(i),vec_old_ptr, &
                        ierr);CHKERRQ(ierr)
    max_change = 0.d0
    do j = 1, grid%nlmax
      change = dabs(vec_new_ptr(j)-vec_old_ptr(j))
      max_change = max(max_change,change)
    enddo
    max_change_global(i) = max_change
    call VecRestoreArray(field%work,vec_new_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArray(field%max_change_vecs(i),vec_old_ptr, &
                            ierr);CHKERRQ(ierr)
    call VecCopy(field%work,field%max_change_vecs(i),ierr);CHKERRQ(ierr)
  enddo
  i = size(max_change_global)
  call MPI_Allreduce(MPI_IN_PLACE,max_change_global,i,MPI_DOUBLE_PRECISION, &
                     MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)

  write(option%io_buffer,'("  --> max change: dpl= ",1pe12.4," dsl= ",&
                         &1pe12.4," dt= ",1pe12.4," dc= ",1pe12.4)') &
    max_change_global(1:4)
  call PrintMsg(option)
  this%max_pressure_change = max_change_global(1)
  this%max_saturation_change = max_change_global(2)
  this%max_temperature_change = max_change_global(3)
  ! hijacking xmol_change for concentration
  this%max_xmol_change = max_change_global(4)

  deallocate(max_change_global)

end subroutine PMTHCMaxChange

! ************************************************************************** !

subroutine PMTHCComputeMassBalance(this,mass_balance_array)
  !
  ! Computes the current solute/water/energy mass in the domain
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use THC_module, only : THCComputeMassBalance

  implicit none

  class(pm_thc_type) :: this
  PetscReal :: mass_balance_array(:)

  call THCComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMTHCComputeMassBalance

! ************************************************************************** !

subroutine PMTHCInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  implicit none

  class(pm_thc_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'thc'

end subroutine PMTHCInputRecord

! ************************************************************************** !

subroutine PMTHCCheckpointBinary(this,viewer)
  !
  ! Checkpoints data associated with THC PM
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  use Checkpoint_module
  use Global_module

  implicit none

  class(pm_thc_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowCheckpointBinary(this,viewer)

end subroutine PMTHCCheckpointBinary

! ************************************************************************** !

subroutine PMTHCRestartBinary(this,viewer)
  !
  ! Restarts data associated with THC PM
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26

  use Checkpoint_module
  use Global_module

  implicit none

  class(pm_thc_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowRestartBinary(this,viewer)

end subroutine PMTHCRestartBinary

! ************************************************************************** !

subroutine PMTHCDestroy(this)
  !
  ! Destroys THC process model
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/09/26
  !

  use THC_module, only : THCDestroy
  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_thc_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  call DeallocateArray(this%max_change_ivar)
  call THCDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMTHCDestroy

end module PM_THC_class
