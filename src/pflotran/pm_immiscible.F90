module PM_Immiscible_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class

  use PFLOTRAN_Constants_module
  use Immiscible_Aux_module

  implicit none

  private

  PetscInt, parameter :: MAX_CHANGE_PRES_NI = 1
  PetscInt, parameter :: MAX_CHANGE_SAT_NI = 2
  PetscInt, parameter :: MAX_RES_WET_EQ = 3
  PetscInt, parameter :: MAX_RES_NONWET_EQ = 4

  type, public, extends(pm_subsurface_flow_type) :: pm_immiscible_type
    PetscInt, pointer :: max_change_ivar(:)
    PetscReal :: max_allow_pres_change_ni
    PetscReal :: pres_change_ts_governor
    PetscReal :: sat_change_ts_governor
    PetscInt :: convergence_flags(MAX_RES_NONWET_EQ)
    PetscReal :: convergence_reals(MAX_RES_NONWET_EQ)
    PetscReal :: sat_update_trunc_ni
    PetscInt :: convergence_verbosity
  contains
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMImmiscibleReadSimOptionsBlock
    procedure, public :: ReadTSBlock => PMImmiscibleReadTSSelectCase
    procedure, public :: ReadNewtonBlock => PMImmiscibleReadNewtonSelectCase
    procedure, public :: Setup => PMImmiscibleSetup
    procedure, public :: InitializeRun => PMImmiscibleInitializeRun
    procedure, public :: InitializeTimestep => PMImmiscibleInitializeTimestep
    procedure, public :: Residual => PMImmiscibleResidual
    procedure, public :: Jacobian => PMImmiscibleJacobian
    procedure, public :: UpdateTimestep => PMImmiscibleUpdateTimestep
    procedure, public :: FinalizeTimestep => PMImmiscibleFinalizeTimestep
    procedure, public :: PreSolve => PMImmisciblePreSolve
    procedure, public :: PostSolve => PMImmisciblePostSolve
    procedure, public :: CheckUpdatePre => PMImmiscibleCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMImmiscibleCheckUpdatePost
    procedure, public :: CheckConvergence => PMImmiscibleCheckConvergence
    procedure, public :: TimeCut => PMImmiscibleTimeCut
    procedure, public :: UpdateSolution => PMImmiscibleUpdateSolution
    procedure, public :: UpdateAuxVars => PMImmiscibleUpdateAuxVars
    procedure, public :: MaxChange => PMImmiscibleMaxChange
    procedure, public :: ComputeMassBalance => PMImmiscibleComputeMassBalance
    procedure, public :: InputRecord => PMImmiscibleInputRecord
    procedure, public :: CheckpointBinary => PMImmiscibleCheckpointBinary
    procedure, public :: RestartBinary => PMImmiscibleRestartBinary
    procedure, public :: Destroy => PMImmiscibleDestroy
  end type pm_immiscible_type

  public :: PMImmiscibleCreate, &
            PMImmiscibleInitObject, &
            PMImmiscibleInitializeRun, &
            PMImmiscibleFinalizeTimestep, &
            PMImmiscibleCheckUpdatePre, &
            PMImmiscibleDestroy

contains

! ************************************************************************** !

function PMImmiscibleCreate()
  !
  ! Creates Immiscible process models shell
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !
  implicit none

  class(pm_immiscible_type), pointer :: PMImmiscibleCreate

  class(pm_immiscible_type), pointer :: immis_pm

  allocate(immis_pm)
  call PMImmiscibleInitObject(immis_pm)

  PMImmiscibleCreate => immis_pm

end function PMImmiscibleCreate

! ************************************************************************** !

subroutine PMImmiscibleInitObject(this)
  !
  ! Creates Immiscible process models shell
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !
  use Option_module
  use String_module
  use Variables_module

  implicit none

  class(pm_immiscible_type) :: this

  call PMSubsurfaceFlowInit(this)
  this%name = 'Immiscible Flow'
  this%header = 'IMMISIBLE FLOW'

  nullify(this%max_change_ivar)

  ! set to UNINITIALIZED_DOUBLE and report error below is set from input
  this%pressure_change_governor = UNINITIALIZED_DOUBLE
  this%temperature_change_governor = UNINITIALIZED_DOUBLE
  this%saturation_change_governor = UNINITIALIZED_DOUBLE
  this%xmol_change_governor = UNINITIALIZED_DOUBLE

  this%check_post_convergence = PETSC_TRUE
  this%convergence_verbosity = 0

  this%max_allow_pres_change_ni = UNINITIALIZED_DOUBLE
  this%pres_change_ts_governor = 5.d5    ! [Pa]
  this%sat_change_ts_governor = 1.d0
  this%sat_update_trunc_ni = UNINITIALIZED_DOUBLE

  this%convergence_flags = 0
  this%convergence_reals = 0.d0

  immis_debug_cell_id = UNINITIALIZED_INTEGER

end subroutine PMImmiscibleInitObject

! ************************************************************************** !

subroutine PMImmiscibleReadSimOptionsBlock(this,input)
  !
  ! Read IMMISCIBLE options input block
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !
  use Immiscible_module
  use Immiscible_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module
  use Utility_module
  use EOS_Water_module

  implicit none

  class(pm_immiscible_type) :: this
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: error_string
!  character(len=MAXSTRINGLENGTH) :: local_error_string
  PetscBool :: found
!  PetscReal :: array(1)
  PetscInt :: temp_int

  option => this%option

  error_string = 'Immiscible Options'

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
      case('VERBOSE_CONVERGENCE')
        this%convergence_verbosity = 1
        call InputReadInt(input,option,temp_int)
        if (.not.InputError(input)) then
          this%convergence_verbosity = temp_int
        else
          call InputDefaultMsg(input,option,keyword)
        endif
      case('NO_ACCUMULATION')
        immis_calc_accum = PETSC_FALSE
      case('NO_FLUX')
        immis_calc_flux = PETSC_FALSE
      case('NO_BCFLUX')
        immis_calc_bcflux = PETSC_FALSE
      case('DEBUG_CELL_ID')
        call InputReadInt(input,option,immis_debug_cell_id)
        call InputErrorMsg(input,option,keyword,error_string)
      case default
        call InputKeywordUnrecognized(input,keyword,'Immiscible Mode',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine PMImmiscibleReadSimOptionsBlock

! ************************************************************************** !

subroutine PMImmiscibleReadTSSelectCase(this,input,keyword,found, &
                                        error_string,option)
  !
  ! Read timestepper settings specific to the IMMISCIBLE process model
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  class(pm_immiscible_type) :: this
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
    case('SAT_CHANGE_TS_GOVERNOR')
      call InputReadDouble(input,option,this%sat_change_ts_governor)
      call InputErrorMsg(input,option,keyword,error_string)
    case('PRES_CHANGE_TS_GOVERNOR')
      call InputReadDouble(input,option,this%pres_change_ts_governor)
      call InputErrorMsg(input,option,keyword,error_string)
      ! units conversion since it is absolute
      call InputReadAndConvertUnits(input,this%pres_change_ts_governor, &
                                    'Pa',keyword,option)
    case default
      found = PETSC_FALSE
  end select

end subroutine PMImmiscibleReadTSSelectCase

! ************************************************************************** !

subroutine PMImmiscibleReadNewtonSelectCase(this,input,keyword,found, &
                                            error_string,option)
  !
  ! Reads input file parameters associated with the IMMISCIBLE process model
  ! Newton solver convergence
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  use Immiscible_Aux_module

  implicit none

  class(pm_immiscible_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  error_string = 'IMMISCIBLE Newton Solver'

  select case(trim(keyword))
    case('ITOL_UPDATE')
      option%io_buffer = 'ITOL_UPDATE not supported with IMMISCIBLE. Please &
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
      call InputReadDouble(input,option,immis_rel_pert)
      call InputErrorMsg(input,option,keyword,error_string)
      ! no units conversion since it is relative
    case('MIN_PRESSURE_PERTURBATION')
      call InputReadDouble(input,option,immis_pres_min_pert)
      call InputErrorMsg(input,option,keyword,error_string)
      call InputReadAndConvertUnits(input,immis_pres_min_pert, &
                                    'Pa',keyword,option)
    case('MAX_ALLOW_PRES_CHANGE_NI')
      call InputReadDouble(input,option,this%max_allow_pres_change_ni)
      call InputErrorMsg(input,option,keyword,error_string)
      ! units conversion since it is absolute
      call InputReadAndConvertUnits(input,this%max_allow_pres_change_ni, &
                                    'Pa',keyword,option)
    case('SATURATION_UPDATE_TRUNCATION_NI')
      call InputReadDouble(input,option,this%sat_update_trunc_ni)
      call InputErrorMsg(input,option,keyword,error_string)
    case default
      found = PETSC_FALSE

  end select

end subroutine PMImmiscibleReadNewtonSelectCase

! ************************************************************************** !

subroutine PMImmiscibleSetup(this)
  !
  ! Sets up auxvars and parameters
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24

  use Immiscible_module
  use EOS_Gas_module
  use EOS_module

  implicit none

  class(pm_immiscible_type) :: this

  call EOSGasSetFMW(immis_fmw(2))
  call this%SetRealization()
  call EOSGasTableProcess(this%realization%option)
  call EOSReferenceDensity(this%realization%option)
  call ImmiscibleSetup(this%realization)
  call PMSubsurfaceFlowSetup(this)

end subroutine PMImmiscibleSetup

! ************************************************************************** !

recursive subroutine PMImmiscibleInitializeRun(this)
  !
  ! Initializes the IMMISCIBLE mode run.
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24

  use Realization_Base_class
  use Patch_module
  use Field_module
  use Material_Aux_module
  use Option_module
  use Variables_module

  implicit none

  class(pm_immiscible_type) :: this

  PetscInt :: i
  PetscErrorCode :: ierr
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
!  PetscInt :: ivar

  patch => this%realization%patch
  field => this%realization%field
  option => this%option

  if (this%steady_state) immis_calc_accum = PETSC_FALSE

  allocate(this%max_change_ivar(4))
  this%max_change_ivar(1) = LIQUID_PRESSURE
  this%max_change_ivar(2) = GAS_PRESSURE
  this%max_change_ivar(3) = LIQUID_SATURATION
  this%max_change_ivar(4) = GAS_SATURATION

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

  if (Initialized(this%temperature_change_governor)) then
    option%io_buffer = 'TEMPERATURE_CHANGE_GOVERNOR &
      &may not be used with IMMISCIBLE.'
    call PrintErrMsg(option)
  endif

end subroutine PMImmiscibleInitializeRun

! ************************************************************************** !

subroutine PMImmiscibleInitializeTimestep(this)
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !
  use Immiscible_module, only : ImmiscibleInitializeTimestep

  implicit none

  class(pm_immiscible_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
  call ImmiscibleInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)

  this%convergence_flags = 0
  this%convergence_reals = 0.d0

end subroutine PMImmiscibleInitializeTimestep

! ************************************************************************** !

subroutine PMImmiscibleFinalizeTimestep(this)
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !
  implicit none

  class(pm_immiscible_type) :: this

  call PMSubsurfaceFlowFinalizeTimestep(this)

end subroutine PMImmiscibleFinalizeTimestep

! ************************************************************************** !

subroutine PMImmisciblePreSolve(this)
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24

  implicit none

  class(pm_immiscible_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMImmisciblePreSolve

! ************************************************************************** !

subroutine PMImmisciblePostSolve(this)
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24

  use Upwind_Direction_module
  use Option_module

  implicit none

  class(pm_immiscible_type) :: this

end subroutine PMImmisciblePostSolve

! ************************************************************************** !

subroutine PMImmiscibleUpdateTimestep(this,update_dt, &
                                 dt,dt_min,dt_max,iacceleration, &
                                 num_newton_iterations,tfac, &
                                 time_step_max_growth_factor)
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Option_module
  use Utility_module, only : Equal

  implicit none

  class(pm_immiscible_type) :: this
  PetscBool :: update_dt
  PetscReal :: dt
  PetscReal :: dt_min ! DO NOT USE (see comment below)
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: sat_ratio, pres_ratio
  PetscReal :: dt_prev

  if (update_dt .and. iacceleration /= 0) then
    dt_prev = dt
    ! calculate the time step ramping factor
    sat_ratio = (2.d0*this%sat_change_ts_governor)/ &
                (this%sat_change_ts_governor+this%max_saturation_change)
    pres_ratio = (2.d0*this%pres_change_ts_governor)/ &
                (this%pres_change_ts_governor+this%max_pressure_change)
    ! pick minimum time step from calc'd ramping factor or maximum ramping factor
    dt = min(min(sat_ratio,pres_ratio)*dt,time_step_max_growth_factor*dt)
    ! make sure time step is within bounds given in the input deck
    dt = min(dt,dt_max)
    if (this%logging_verbosity > 0) then
      if (Equal(dt,dt_max)) then
        string = 'maximum time step size'
      else if (min(sat_ratio,pres_ratio) > time_step_max_growth_factor) then
        string = 'maximum time step growth factor'
      else if (sat_ratio < pres_ratio) then
        string = 'liquid saturation governor'
      else
        string = 'liquid pressure governor'
      endif
        string = 'TS update: ' // trim(string)
      call PrintMsg(this%option,string)
    endif
  endif

  if (Initialized(this%cfl_governor)) then
    call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt,dt_max)
  endif

end subroutine PMImmiscibleUpdateTimestep

! ************************************************************************** !

subroutine PMImmiscibleResidual(this,snes,xx,r,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !
  use Immiscible_module, only : ImmiscibleResidual
  use Debug_module
  use Grid_module

  implicit none

  class(pm_immiscible_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  Mat :: M

  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  ! calculate residual
  if (immis_simultaneous_res_jac_calc) then
    call SNESGetJacobian(snes,M,PETSC_NULL_MAT,PETSC_NULL_FUNCTION, &
                         PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
    call ImmiscibleResidual(snes,xx,r,M,this%realization,this%debug,ierr)
  else
    call ImmiscibleResidual(snes,xx,r,PETSC_NULL_MAT,this%realization, &
                            this%debug,ierr)
  endif

end subroutine PMImmiscibleResidual

! ************************************************************************** !

subroutine PMImmiscibleJacobian(this,snes,xx,A,B,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !
  use Debug_module
  use Option_module

  implicit none

  class(pm_immiscible_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  PetscReal :: norm

  ! the Jacobian was already calculated in PMImmiscibleResidual

  if (this%debug%matview_Matrix) then
    call DebugMatView(this%debug,A,'Immis_jacobian','', &
                      immis_ts_count,immis_ts_cut_count, &
                      immis_ni_count,this%option)
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

  immis_ni_count = immis_ni_count + 1

end subroutine PMImmiscibleJacobian

! ************************************************************************** !

subroutine PMImmiscibleCheckUpdatePre(this,snes,X,dX,changed,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !
  use Grid_module
  use Material_Aux_module
  use Option_module
  use Patch_module

  implicit none

  class(pm_immiscible_type) :: this
  SNES :: snes
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr

#if 0
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscInt :: p_index
  PetscReal :: pc, p_target
  PetscReal :: p0, p1, dp
  PetscReal :: sl
  PetscReal :: tempreal
  PetscReal, pointer :: X_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscBool :: sat_update_trunc_flag
  type(immiscible_auxvar_type), pointer :: immiscible_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  patch => this%realization%patch
  grid => patch%grid
  option => this%realization%option

  immiscible_auxvars => patch%aux%Immiscible%auxvars
  material_auxvars => patch%aux%Material%auxvars

  this%convergence_flags = 0
  this%convergence_reals = 0.d0
  changed = PETSC_FALSE

  sat_update_trunc_flag = Initialized(this%sat_update_trunc_ni)

  call VecGetArray(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(X,X_p,ierr);CHKERRQ(ierr)
#if 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle

    offset = (local_id-1)*option%nflowdof
    if (immis_liq_flow_eq > 0) then
      p_index = offset+immis_liq_flow_eq
      dp = -dX_p(p_index)
      p0 = X_p(p_index)
      p1 = p0+dp
      if (unsat_to_sat_damping_flag) then
        sl = immiscible_auxvars(ZERO_INTEGER,ghosted_id)%sat + &
            sign(this%sat_update_trunc_ni,dp)
        call patch%characteristic_curves_array( &
              patch%cc_id(ghosted_id))%ptr%saturation_function% &
                CapillaryPressure(sl,pc,tempreal,option)
        if (pc > 0.d0) then
          p_target = p_ref-pc
          if ((dp >= 0.d0 .and. p1 > p_target) .or. &
              (dp < 0.d0 .and. p1 < p_target)) then
            dX_p(p_index) = p0-p_target ! p1 = p0 - dX_p()
            changed = PETSC_TRUE
          endif
        endif
        ! update these incase used below
        dp = -dX_p(p_index)
        p1 = p0+dp
      endif
      if (unsat_to_sat_damping_flag) then
        ! the following initiate damping when transitioning from
        ! unsaturated to saturated state
        if (p0 < p_ref .and. p1 > p_ref) then
          dX_p(p_index) = this%unsat_to_sat_pres_damping_ni*dX_p(p_index)
          changed = PETSC_TRUE
        endif
      endif
    endif
  enddo
#endif
  call VecRestoreArray(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(X,X_p,ierr);CHKERRQ(ierr)
#endif

end subroutine PMImmiscibleCheckUpdatePre

! ************************************************************************** !

subroutine PMImmiscibleCheckUpdatePost(this,snes,X0,dX,X1,dX_changed, &
                                       X1_changed,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Material_Aux_module
  use Immiscible_Aux_module

  implicit none

  class(pm_immiscible_type) :: this
  SNES :: snes
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr

#if 0
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch

  PetscReal, pointer :: X0_p(:)
  PetscReal, pointer :: X1_p(:)
  PetscReal, pointer :: dX_p(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscReal :: tempreal

  PetscBool :: converged_wetting_pressure
  PetscBool :: converged_nonwetting_pressure
  PetscReal :: max_abs_pressure_change_NI
  PetscInt :: max_abs_pressure_change_NI_cell

  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch

  ! If these are changed from true, we must add a global reduction on both
  ! variables to ensure that their values match across all processes. Otherwise
  ! PETSc will throw an error in debug mode or ignore the error in optimized.
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE

  call VecGetArray(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(X0,X0_p,ierr);CHKERRQ(ierr)
  call VecGetArray(X1,X1_p,ierr);CHKERRQ(ierr)
  converged_wetting_pressure = PETSC_TRUE
  converged_nonwetting_pressure = PETSC_TRUE
  max_abs_pressure_change_NI = 0.d0
  max_abs_pressure_change_NI_cell = 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle

#if 0
    offset = (local_id-1)*option%nflowdof
    if (immis_liq_flow_eq > 0) then
      ! maximum absolute change in liquid pressure over Newton iteration
      tempreal = dabs(dX_p(offset+immis_liq_flow_eq))
      if (tempreal > dabs(max_abs_pressure_change_NI)) then
        max_abs_pressure_change_NI_cell = grid%nG2A(ghosted_id)
        max_abs_pressure_change_NI = tempreal
      endif
    endif
    if (immis_sol_tran_eq > 0) then
      ! maximum absolute change in liquid pressure over Newton iteration
      tempreal = dabs(dX_p(offset+immis_sol_tran_eq))
      if (tempreal > dabs(max_abs_conc_change_NI)) then
        max_abs_conc_change_NI_cell = grid%nG2A(ghosted_id)
        max_abs_conc_change_NI = tempreal
      endif
    endif
#endif
  enddo

  if (Initialized(this%max_allow_pres_change_ni) .and. &
      max_abs_pressure_change_NI > this%max_allow_pres_change_ni) then
    converged_wetting_pressure = PETSC_FALSE
  endif

  ! the following flags are used in detemining convergence
  if (.not.converged_wetting_pressure) then
    this%convergence_flags(MAX_CHANGE_PRES_NI) = &
      max_abs_pressure_change_NI_cell
  endif

  ! the following flags are for REPORTING purposes only
  this%convergence_reals(MAX_CHANGE_PRES_NI) = max_abs_pressure_change_NI

  call VecRestoreArray(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(X0,X0_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(X1,X1_p,ierr);CHKERRQ(ierr)
#endif

end subroutine PMImmiscibleCheckUpdatePost

! ************************************************************************** !

subroutine PMImmiscibleCheckConvergence(this,snes,it,xnorm,unorm, &
                                        fnorm,reason,ierr)
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Material_Aux_module
  use Immiscible_Aux_module
  use Convergence_module

  implicit none

  class(pm_immiscible_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

  Vec :: residual_vec
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_tpdt_p(:)
  PetscReal, pointer :: X1_p(:)
  character(len=10) :: reason_string

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(immiscible_auxvar_type), pointer :: immiscible_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

!  PetscInt :: local_id, ghosted_id
  PetscInt :: converged_flag

  PetscReal :: max_abs_res_wet_
  PetscInt :: max_abs_res_wet_cell
  PetscReal :: max_abs_res_nonwet_
  PetscInt :: max_abs_res_nonwet_cell
!  PetscInt :: offset
  PetscMPIInt :: int_mpi

!  PetscReal :: accumulation
!  PetscReal :: residual
!  PetscReal :: tempreal

  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch
  immiscible_auxvars => patch%aux%Immiscible%auxvars
  material_auxvars => patch%aux%Material%auxvars

  residual_vec = field%flow_r
  ! check residual terms

  call VecGetArrayRead(residual_vec,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%flow_accum_tpdt,accum_tpdt_p,ierr);CHKERRQ(ierr)
  call VecGetArrayRead(field%flow_xx,X1_p,ierr);CHKERRQ(ierr)

  max_abs_res_wet_ = UNINITIALIZED_DOUBLE
  max_abs_res_wet_cell = UNINITIALIZED_INTEGER
  max_abs_res_nonwet_ = UNINITIALIZED_DOUBLE
  max_abs_res_nonwet_cell = UNINITIALIZED_INTEGER

#if 0
  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%nflowdof
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (immis_liq_flow_eq > 0) then
      residual = r_p(offset+immis_liq_flow_eq)
      accumulation = accum2_p(offset+immis_liq_flow_eq)
      ! residual
      tempreal = dabs(residual)
      if (tempreal > max_abs_res_wet_) then
        max_abs_res_wet_ = tempreal
        max_abs_res_wet_cell = grid%nG2A(ghosted_id)
      endif
    endif
    if (immis_sol_tran_eq > 0) then
      residual = r_p(offset+immis_sol_tran_eq)
      accumulation = accum2_p(offset+immis_sol_tran_eq)
      ! residual
      tempreal = dabs(residual)
      if (tempreal > max_abs_res_nonwet_) then
        max_abs_res_nonwet_ = tempreal
        max_abs_res_nonwet_cell = grid%nG2A(ghosted_id)
      endif
    endif
  enddo
#endif

  ! the following flags are used in detemining convergence
  ! currently none

  ! the following flags are for REPORTING purposes only
  this%convergence_flags(MAX_RES_WET_EQ) = max_abs_res_wet_cell
  this%convergence_reals(MAX_RES_WET_EQ) = max_abs_res_wet_
  this%convergence_flags(MAX_RES_NONWET_EQ) = max_abs_res_nonwet_cell
  this%convergence_reals(MAX_RES_NONWET_EQ) = max_abs_res_nonwet_

  if (this%convergence_verbosity >= 10) then
    print *, option%myrank, &
      this%convergence_flags(MAX_RES_WET_EQ), &
      this%convergence_reals(MAX_RES_WET_EQ), &
      this%convergence_flags(MAX_RES_NONWET_EQ), &
      this%convergence_reals(MAX_RES_NONWET_EQ), &
      this%convergence_flags(MAX_CHANGE_PRES_NI), &
      this%convergence_reals(MAX_CHANGE_PRES_NI)
  endif

  int_mpi = size(this%convergence_flags)
  call MPI_Allreduce(MPI_IN_PLACE,this%convergence_flags,int_mpi,MPIU_INTEGER, &
                     MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
  int_mpi = size(this%convergence_reals)
  call MPI_Allreduce(MPI_IN_PLACE,this%convergence_reals,int_mpi, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                     ierr);CHKERRQ(ierr)

  ! these conditionals cannot change order
  reason_string = '---| '
  converged_flag = CONVERGENCE_CONVERGED
  if (this%convergence_flags(MAX_CHANGE_PRES_NI) > 0) then
    reason_string(1:1) = 'P'
    converged_flag = CONVERGENCE_KEEP_ITERATING
  endif
  if (this%convergence_flags(MAX_CHANGE_SAT_NI) > 0) then
    reason_string(2:2) = 'S'
    converged_flag = CONVERGENCE_KEEP_ITERATING
  endif

  if (this%convergence_verbosity > 0 .and. &
      OptionPrintToScreen(option)) then
    if (option%comm%size > 1) then
      write(*,'(4x,"Rsn: ",a10,2es10.2)') reason_string, &
        this%convergence_reals(MAX_RES_WET_EQ), &
        this%convergence_reals(MAX_CHANGE_PRES_NI)
    else if (grid%nmax > 9999) then
      write(*,'(4x,"Rsn: ",a10,2(i8,es10.2))') reason_string, &
        this%convergence_flags(MAX_RES_WET_EQ), &
        this%convergence_reals(MAX_RES_WET_EQ), &
        this%convergence_flags(MAX_CHANGE_PRES_NI), &
        this%convergence_reals(MAX_CHANGE_PRES_NI)
    else
      write(*,'(4x,"Rsn: ",a10,2(i5,es10.2))') reason_string, &
        this%convergence_flags(MAX_RES_WET_EQ), &
        this%convergence_reals(MAX_RES_WET_EQ), &
        this%convergence_flags(MAX_CHANGE_PRES_NI), &
        this%convergence_reals(MAX_CHANGE_PRES_NI)
    endif
  endif

  if (Initialized(this%max_allow_pres_change_ni)) then
    option%convergence = converged_flag
  else
    ! forced standard 2 norms
    option%convergence = CONVERGENCE_OFF
  endif

  call VecRestoreArrayRead(residual_vec,r_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(field%flow_accum_tpdt,accum_tpdt_p, &
                           ierr);CHKERRQ(ierr)
  call VecRestoreArrayRead(field%flow_xx,X1_p,ierr);CHKERRQ(ierr)

  call PMSubsurfaceFlowCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                        reason,ierr)

end subroutine PMImmiscibleCheckConvergence

! ************************************************************************** !

subroutine PMImmiscibleTimeCut(this)
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !

  use Immiscible_module, only : ImmiscibleTimeCut

  implicit none

  class(pm_immiscible_type) :: this

  call PMSubsurfaceFlowTimeCut(this)
  call ImmiscibleTimeCut(this%realization)

  this%convergence_flags = 0
  this%convergence_reals = 0.d0

end subroutine PMImmiscibleTimeCut

! ************************************************************************** !

subroutine PMImmiscibleUpdateSolution(this)
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !

  use Immiscible_module, only : ImmiscibleUpdateSolution, &
                           ImmiscibleMapBCAuxVarsToGlobal

  implicit none

  class(pm_immiscible_type) :: this

  call PMSubsurfaceFlowUpdateSolution(this)
  call ImmiscibleUpdateSolution(this%realization)
  call ImmiscibleMapBCAuxVarsToGlobal(this%realization)

end subroutine PMImmiscibleUpdateSolution

! ************************************************************************** !

subroutine PMImmiscibleUpdateAuxVars(this)
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  use Immiscible_module, only : ImmiscibleUpdateAuxVars

  implicit none

  class(pm_immiscible_type) :: this

  call ImmiscibleUpdateAuxVars(this%realization)

end subroutine PMImmiscibleUpdateAuxVars

! ************************************************************************** !

subroutine PMImmiscibleMaxChange(this)
  !
  ! Not needed given ImmiscibleMaxChange is called in PostSolve
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use String_module
  use Immiscible_Aux_module

  implicit none

  class(pm_immiscible_type) :: this

#if 0
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_old_ptr(:), vec_new_ptr(:)
  PetscReal, allocatable :: max_change_global(:)
  PetscReal :: max_change, change
  PetscInt :: i, j
  PetscInt :: ivar
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
#if 0
  ivar = 1
  if (immis_liq_flow_eq > 0) then
    write(option%io_buffer,'("  --> max change: dpl= ",1pe12.4, " dsl= ",&
                           &1pe12.4)') &
      max_change_global(ivar:ivar+1)
    this%max_pressure_change = max_change_global(ivar)
    this%max_saturation_change = max_change_global(ivar+1)
    ivar = ivar+2
    call PrintMsg(option)
  endif
  if (immis_sol_tran_eq > 0) then
    write(option%io_buffer,'(19x,"dc= ",1pe12.4)') max_change_global(ivar)
    ! hijacking xmol_change
    this%max_xmol_change = max_change_global(ivar)
    ivar = ivar+1
    call PrintMsg(option)
  endif
#endif
  deallocate(max_change_global)
#endif

end subroutine PMImmiscibleMaxChange

! ************************************************************************** !

subroutine PMImmiscibleComputeMassBalance(this,mass_balance_array)
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !

  use Immiscible_module, only : ImmiscibleComputeMassBalance

  implicit none

  class(pm_immiscible_type) :: this
  PetscReal :: mass_balance_array(:)

  call ImmiscibleComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMImmiscibleComputeMassBalance

! ************************************************************************** !

subroutine PMImmiscibleInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 05/17/24
  !

  implicit none

  class(pm_immiscible_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'immiscible'

end subroutine PMImmiscibleInputRecord

! ************************************************************************** !

subroutine PMImmiscibleCheckpointBinary(this,viewer)
  !
  ! Checkpoints data associated with Immiscible PM
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24

  use Checkpoint_module
  use Global_module

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_immiscible_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowCheckpointBinary(this,viewer)

end subroutine PMImmiscibleCheckpointBinary

! ************************************************************************** !

subroutine PMImmiscibleRestartBinary(this,viewer)
  !
  ! Restarts data associated with Immiscible PM
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24

  use Checkpoint_module
  use Global_module

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_immiscible_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowRestartBinary(this,viewer)

end subroutine PMImmiscibleRestartBinary

! ************************************************************************** !

subroutine PMImmiscibleDestroy(this)
  !
  ! Destroys Immiscible process model
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/24
  !

  use Immiscible_module, only : ImmiscibleDestroy
  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_immiscible_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  call DeallocateArray(this%max_change_ivar)
  call ImmiscibleDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMImmiscibleDestroy

end module PM_Immiscible_class
