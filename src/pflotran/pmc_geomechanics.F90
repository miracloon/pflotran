module PMC_Geomechanics_class

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PMC_Base_class
  use Realization_Subsurface_class
  use Geomechanics_Realization_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pmc_base_type) :: pmc_geomechanics_type
    class(realization_subsurface_type), pointer :: subsurf_realization
    class(realization_geomech_type), pointer :: geomech_realization
  contains
    procedure, public :: Init => PMCGeomechanicsInit
    procedure, public :: SetupSolvers => PMCGeomechanicsSetupSolvers
    procedure, public :: InitializeRun => PMCGeomechanicsInitializeRun
    procedure, public :: RunToTime => PMCGeomechanicsRunToTime
    procedure, public :: GetAuxData => PMCGeomechanicsGetAuxData
    procedure, public :: SetAuxData => PMCGeomechanicsSetAuxData
    procedure, public :: CheckpointBinary => PMCGeomechanicsCheckpointBinary
    procedure, public :: RestartBinary => PMCGeomechanicsRestartBinary
    procedure, public :: CheckpointHDF5 => PMCGeomechanicsCheckpointHDF5
    procedure, public :: RestartHDF5 => PMCGeomechanicsRestartHDF5
    procedure, public :: Destroy => PMCGeomechanicsDestroy
  end type pmc_geomechanics_type

  public :: PMCGeomechanicsCreate

contains

! ************************************************************************** !

function PMCGeomechanicsCreate()
  !
  ! This routine allocates and initializes a new object.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !

  implicit none

  class(pmc_geomechanics_type), pointer :: PMCGeomechanicsCreate

  class(pmc_geomechanics_type), pointer :: pmc

#ifdef DEBUG
  print *, 'PMCGeomechanicsCreate%Create()'
#endif

  allocate(pmc)
  call pmc%Init()

  PMCGeomechanicsCreate => pmc

end function PMCGeomechanicsCreate

! ************************************************************************** !

subroutine PMCGeomechanicsInit(this)
  !
  ! This routine initializes a new process model coupler object.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !

  implicit none

  class(pmc_geomechanics_type) :: this

#ifdef DEBUG
  print *, 'PMCGeomechanics%Init()'
#endif

  call PMCBaseInit(this)
  nullify(this%subsurf_realization)
  nullify(this%geomech_realization)

end subroutine PMCGeomechanicsInit

! ************************************************************************** !

subroutine PMCGeomechanicsSetupSolvers(this)
  !
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  !
  use Geomechanics_Discretization_module
  use Timestepper_KSP_class
  use PM_Base_class
  use PM_Base_Pointer_module
  use Option_module
  use Solver_module

  implicit none

  class(pmc_geomechanics_type) :: this

  class(realization_geomech_type), pointer :: geomech_realization
  class(geomech_discretization_type), pointer :: geomech_discretization
  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: dm_mat_type_found
  PetscErrorCode :: ierr

#ifdef DEBUG
  call PrintMsg(this%option,'PMCGeomechanicsSetupSolvers')
#endif

  option => this%option
  solver => this%timestepper%solver
  geomech_realization => this%geomech_realization
  geomech_discretization => geomech_realization%geomech_discretization

  select type(ts=>this%timestepper)
    class is(timestepper_KSP_type)
    class default
      option%io_buffer = 'A KSP timestepper must be used for geomechanics.'
      call PrintErrMsg(option)
  end select

  call SolverCreateKSP(solver,option%mycomm)

  call PrintMsg(option,"  Beginning setup of GEOMECH KSP")
  call KSPSetOptionsPrefix(solver%ksp,"geomech_",ierr);CHKERRQ(ierr)
  call SolverCheckCommandLine(solver)

  if (Uninitialized(solver%Mpre_mat_type) .and. &
      Uninitialized(solver%M_mat_type)) then
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,'geomech_', &
                               '-dm_mat_type',string, &
                               dm_mat_type_found,ierr);CHKERRQ(ierr)
    if (dm_mat_type_found) then
      solver%Mpre_mat_type = trim(string)
    else
      solver%Mpre_mat_type = MATBAIJ
    endif
    solver%M_mat_type = solver%Mpre_mat_type
  else if (Uninitialized(solver%Mpre_mat_type)) then
    if (solver%M_mat_type == MATMFFD) then
      solver%Mpre_mat_type = MATBAIJ
    else
      solver%Mpre_mat_type = solver%M_mat_type
    endif
  else if (Uninitialized(solver%M_mat_type)) then
    solver%M_mat_type = solver%Mpre_mat_type
  endif

  call GeomechDiscretizationCreateMatrix(geomech_realization% &
                                         geomech_discretization,NGEODOF, &
                                         solver%Mpre_mat_type, &
                                         solver%Mpre,option)

  call MatSetOptionsPrefix(solver%Mpre,"geomech_",ierr);CHKERRQ(ierr)
  solver%M = solver%Mpre
  geomech_realization%geomech_field%A = solver%M


  ! Have PETSc do a SNES_View() at the end of each solve if verbosity > 0.
  if (option%verbosity >= 2) then
    string = '-geomech_ksp_view'
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS,string, &
                                  ierr);CHKERRQ(ierr)
    string = '-geomech_ksp_monitor'
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS,string, &
                                  ierr);CHKERRQ(ierr)
  endif


  ! call KSPSetOperators(solver%ksp,solver%M,solver%Mpre,ierr);CHKERRQ(ierr)

  call SolverSetKSPOptions(solver,option)

  call PrintMsg(option,"  Finished setting up GEOMECH KSP ")

end subroutine PMCGeomechanicsSetupSolvers

! ************************************************************************** !

recursive subroutine PMCGeomechanicsInitializeRun(this)
  !
  ! Initializes the geomechanics process model coupler
  !
  ! Author: Glenn Hammond
  ! Date: 03/31/25

  use Timestepper_Base_class

  implicit none

  class(pmc_geomechanics_type) :: this

  PetscReal :: target_time
  PetscInt :: local_stop_flag

  call PMCBaseInitializeRun(this)

  target_time = 0.d0
  local_stop_flag = TS_CONTINUE

  if (this%option%geomechanics%split_scheme == GEOMECH_DRAINED_SPLIT) &
    return

  ! On restart, the geomech state is already restored from checkpoint.
  ! The RunToTime(0.0) initialization solve must be skipped because the
  ! KSP timestepper's target_time is already past 0 (= the restart time),
  ! and calling SetTargetTime(0.0) would produce a negative dt.
  ! Also seed option%flow_dt from the restored flow timestepper so that
  ! the first post-restart coupling step has a valid value for the
  ! overwrite in PMCGeomechanicsRunToTime.
  if (this%option%restart_flag) then
    if (associated(this%peer)) then
      this%option%flow_dt = this%peer%timestepper%dt
    endif
    return
  endif

  ! continue to initialize geomechanics
  ! used for GEOMECH_FIXED_STRESS_SPLIT and
  ! GEOMECH_FIXED_STRAIN_SPLIT
  this%timestepper%steps = -1

  call this%RunToTime(target_time,local_stop_flag)

end subroutine PMCGeomechanicsInitializeRun

! ************************************************************************** !

recursive subroutine PMCGeomechanicsRunToTime(this,sync_time,stop_flag)
  !
  ! This routine runs the geomechanics simulation.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !

  use Timestepper_Base_class
  use Option_module
  use PM_Base_class
  use Output_Geomechanics_module

  implicit none

  class(pmc_geomechanics_type), target :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  PetscInt :: local_stop_flag
  PetscBool :: sync_flag
  PetscBool :: snapshot_plot_flag
  PetscBool :: observation_plot_flag
  PetscBool :: massbal_plot_flag
  PetscBool :: conserv_plot_flag
  PetscBool :: checkpoint_flag
  PetscBool :: peer_already_run_to_time

  class(pm_base_type), pointer :: cur_pm
  PetscReal :: geomech_target_time
  PetscReal :: geomech_dt
  PetscErrorCode :: ierr

  if (stop_flag == TS_STOP_FAILURE) return

  if (this%stage /= 0) then
    call PetscLogStagePush(this%stage,ierr);CHKERRQ(ierr)
  endif

  this%option%io_buffer = trim(this%name) // ':' // trim(this%pm_list%name)
  call PrintVerboseMsg(this%option)

  ! Get data of other process-model
  call this%GetAuxData()

  local_stop_flag = 0

  call SetOutputFlags(this)
  sync_flag = PETSC_FALSE
  snapshot_plot_flag = PETSC_FALSE
  observation_plot_flag = PETSC_FALSE
  massbal_plot_flag = PETSC_FALSE
  conserv_plot_flag = PETSC_FALSE

  ! jaa: assignment of target time and dt here before calling
  ! SetTargetTime ensures snapshot_plot_flag is set correctly
  ! especially for cases when xmf output is periodic
  select case(this%option%geomechanics%split_scheme)
    case(GEOMECH_FIXED_STRAIN_SPLIT, GEOMECH_FIXED_STRESS_SPLIT)
      this%timestepper%target_time = this%option%time
      this%timestepper%dt = this%option%flow_dt
  end select

  call this%timestepper%SetTargetTime(sync_time,this%option,local_stop_flag, &
                                      sync_flag, &
                                      snapshot_plot_flag, &
                                      observation_plot_flag, &
                                      massbal_plot_flag, &
                                      conserv_plot_flag,checkpoint_flag)

  ! Save the correct geomech target_time and dt before the overwrite below.
  ! For fixed-stress/fixed-strain splits, lines 271-272 overwrite the KSP
  ! timestepper's dt and target_time with flow values needed for the solve.
  ! We restore the correct geomech values at the end of this routine so that
  ! checkpoint writes consistent state.
  geomech_target_time = this%timestepper%target_time
  geomech_dt = this%timestepper%dt

  ! overwrites target time and dt for geomech when flow is the master pm
  select case(this%option%geomechanics%split_scheme)
    case(GEOMECH_FIXED_STRAIN_SPLIT, GEOMECH_FIXED_STRESS_SPLIT)
      this%timestepper%dt = this%option%flow_dt
      this%timestepper%target_time = this%option%time
    case default
      this%option%dt = this%timestepper%dt
      this%option%time = this%timestepper%target_time-this%timestepper%dt
  end select

  call this%StepDT(local_stop_flag)

  ! Have to loop over all process models coupled in this object and update
  ! the time step size.  Still need code to force all process models to
  ! use the same time step size if tightly or iteratively coupled.
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    ! have to update option%time for conditions
    this%option%time = this%timestepper%target_time
    call cur_pm%UpdateSolution()
    ! Geomechanics PM does not have an associate time
    !call this%timestepper%UpdateDT(cur_pm)
    cur_pm => cur_pm%next
  enddo

  ! Run underlying process model couplers
  if (associated(this%child)) then
    ! Set data needed by process-model
    call this%SetAuxData()
    call this%child%RunToTime(this%timestepper%target_time,local_stop_flag)
    call this%GetAuxData()
  endif

  peer_already_run_to_time = PETSC_FALSE
  if (sync_flag .and. associated(this%peer)) then
    ! synchronize peers
    call this%SetAuxData()
    ! Run neighboring process model couplers
    call this%peer%RunToTime(this%timestepper%target_time, &
                             local_stop_flag)
    peer_already_run_to_time = PETSC_TRUE
    call this%GetAuxData()
  endif

  if (this%timestepper%time_step_cut_flag) then
    snapshot_plot_flag = PETSC_FALSE
  endif
  ! however, if we are using the modulus of the output_option%imod, we may
  ! still print
  if (mod(this%timestepper%steps,this%pm_list% &
          output_option%periodic_snap_output_ts_imod) == 0) then
    snapshot_plot_flag = PETSC_TRUE
  endif
  if (mod(this%timestepper%steps,this%pm_list%output_option% &
          periodic_obs_output_ts_imod) == 0) then
    observation_plot_flag = PETSC_TRUE
  endif
  if (mod(this%timestepper%steps,this%pm_list%output_option% &
          periodic_msbl_output_ts_imod) == 0) then
    massbal_plot_flag = PETSC_TRUE
  endif
  if (mod(this%timestepper%steps,this%pm_list%output_option% &
          periodic_cons_output_ts_imod) == 0) then
    conserv_plot_flag = PETSC_TRUE
  endif

  call OutputGeomechanics(this%geomech_realization,snapshot_plot_flag, &
                          observation_plot_flag,massbal_plot_flag, &
                          conserv_plot_flag)
  ! Set data needed by process-model
  call this%SetAuxData()

  ! Run neighboring process model couplers
  if (associated(this%peer) .and. .not.peer_already_run_to_time) then
    call this%SetAuxData()
    call this%peer%RunToTime(sync_time,local_stop_flag)
    call this%GetAuxData()
  endif

  ! Restore the correct geomech target_time and dt for the KSP timestepper.
  ! The overwrite at lines 271-272 set the KSP timestepper to flow values
  ! needed during StepDT, but the timestepper must reflect the actual geomech
  ! coupling step for checkpoint correctness.
  select case(this%option%geomechanics%split_scheme)
    case(GEOMECH_FIXED_STRAIN_SPLIT, GEOMECH_FIXED_STRESS_SPLIT)
      this%timestepper%target_time = geomech_target_time
      this%timestepper%dt = geomech_dt
  end select

  stop_flag = max(stop_flag,local_stop_flag)

  if (this%stage /= 0) then
    call PetscLogStagePop(ierr);CHKERRQ(ierr)
  endif

end subroutine PMCGeomechanicsRunToTime

! ************************************************************************** !

subroutine PMCGeomechanicsSetAuxData(this)
  !
  ! This routine updates data in simulation_aux that is required by other
  ! process models.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !
  use Option_module
  use Grid_module
  use Discretization_module
  use Geomechanics_Subsurface_Properties_module
  use Parameter_module
  use Global_Aux_module

  implicit none

  class(pmc_geomechanics_type) :: this

  type(grid_type), pointer :: grid
  PetscInt :: local_id, ghosted_id
  PetscScalar, pointer :: por0_p(:)
  PetscScalar, pointer :: por_p(:)
  PetscScalar, pointer :: perm0_p(:)
  PetscScalar, pointer :: perm_p(:)
  PetscScalar, pointer :: strain_p(:)
  PetscScalar, pointer :: stress_p(:)
  PetscScalar, pointer :: press_p(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscReal :: local_stress(6), local_strain(6)
  PetscReal :: local_pressure, local_temp
  PetscErrorCode :: ierr
  PetscReal :: por_new
  PetscReal :: perm_new
  PetscInt :: i
  PetscInt :: parameter_index

  PetscInt :: strainv_0_id, press_0_id
  PetscInt :: strainv_id
  PetscInt :: flow_porosity_id
  PetscInt :: temp_0_id

  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvar
  type(option_type), pointer :: option

  Vec :: geomech_vec
  Vec :: subsurf_vec

  PetscReal :: strain_vol
  PetscInt :: imat
  PetscInt :: id_porosity_mech
  PetscInt :: id_pressure_mech

#if GEOMECH_DEBUG
  PetscViewer :: viewer
  print *, 'PMCGeomechSetAuxData'
#endif

  option => this%option

  ! If at initialization stage, do nothing
!  if (this%timestepper%steps == 0) return

  select type(pmc => this)
    class is(pmc_geomechanics_type)
      if (option%geomechanics%flow_coupling == &
            GEOMECH_TWO_WAY_COUPLED .or. &
          option%geomechanics%geophysics_coupling == &
            GEOMECH_ERT_COUPLING) then

        grid => pmc%subsurf_realization%patch%grid
        global_auxvars => pmc%subsurf_realization% &
                          patch%aux%Global%auxvars

        ! Find the number of geomech grid nodes for each flow cell
        call VecDuplicate(pmc%geomech_realization%geomech_field%strain, &
                          geomech_vec,ierr);CHKERRQ(ierr)
        call VecSet(geomech_vec,1.d0,ierr);CHKERRQ(ierr)

        call VecDuplicate(pmc%sim_aux%subsurf_strain,subsurf_vec, &
                          ierr);CHKERRQ(ierr)
        call VecSet(subsurf_vec,0.d0,ierr);CHKERRQ(ierr)

        call VecScatterBegin(pmc%sim_aux%geomechanics_to_subsurf,geomech_vec, &
                             subsurf_vec,ADD_VALUES,SCATTER_FORWARD, &
                             ierr);CHKERRQ(ierr)
        call VecScatterEnd(pmc%sim_aux%geomechanics_to_subsurf,geomech_vec, &
                           subsurf_vec,ADD_VALUES,SCATTER_FORWARD, &
                           ierr);CHKERRQ(ierr)


#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'subsurf_vec_adjacency_count.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call VecView(subsurf_vec,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

       ! Save strain dataset in sim_aux%subsurf_strain
        call VecSet(pmc%sim_aux%subsurf_strain,0.d0,ierr);CHKERRQ(ierr)
        call VecScatterBegin(pmc%sim_aux%geomechanics_to_subsurf, &
                             pmc%geomech_realization%geomech_field%strain, &
                             pmc%sim_aux%subsurf_strain,ADD_VALUES, &
                             SCATTER_FORWARD,ierr);CHKERRQ(ierr)
        call VecScatterEnd(pmc%sim_aux%geomechanics_to_subsurf, &
                           pmc%geomech_realization%geomech_field%strain, &
                           pmc%sim_aux%subsurf_strain,ADD_VALUES, &
                           SCATTER_FORWARD,ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'subsurf_strain_vector_before_averaging.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call VecView(pmc%sim_aux%subsurf_strain,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

        ! Save stress dataset in sim_aux%subsurf_stress
        call VecSet(pmc%sim_aux%subsurf_stress,0.d0,ierr);CHKERRQ(ierr)
        call VecScatterBegin(pmc%sim_aux%geomechanics_to_subsurf, &
                             pmc%geomech_realization%geomech_field%stress, &
                             pmc%sim_aux%subsurf_stress,ADD_VALUES, &
                             SCATTER_FORWARD,ierr);CHKERRQ(ierr)
        call VecScatterEnd(pmc%sim_aux%geomechanics_to_subsurf, &
                           pmc%geomech_realization%geomech_field%stress, &
                           pmc%sim_aux%subsurf_stress,ADD_VALUES, &
                           SCATTER_FORWARD,ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'subsurf_stress_vector_before_averaging.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call VecView(pmc%sim_aux%subsurf_stress,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

        ! Calculate the average stress and strain
        call VecPointwiseDivide(pmc%sim_aux%subsurf_strain, &
                                pmc%sim_aux%subsurf_strain,subsurf_vec, &
                                ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'subsurf_strain_vector_after_averaging.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call VecView(pmc%sim_aux%subsurf_strain,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

        call VecPointwiseDivide(pmc%sim_aux%subsurf_stress, &
                                pmc%sim_aux%subsurf_stress,subsurf_vec, &
                                ierr);CHKERRQ(ierr)

#if GEOMECH_DEBUG
  call PetscViewerASCIIOpen(option%mycomm, &
                            'subsurf_stress_vector_after_averaging.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call VecView(pmc%sim_aux%subsurf_stress,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

        ! Strain
        call VecGetArray(pmc%sim_aux%subsurf_strain,strain_p, &
                            ierr);CHKERRQ(ierr)
        ! Stress
        call VecGetArray(pmc%sim_aux%subsurf_stress,stress_p, &
                            ierr);CHKERRQ(ierr)

        if (option%geomechanics%geophysics_coupling == &
                                GEOMECH_ERT_COUPLING) then
          ! stress
          call VecGetArray(pmc%subsurf_realization%field%work, &
                              vec_ptr,ierr);CHKERRQ(ierr)
          do local_id = 1, grid%nlmax
            vec_ptr(local_id) = 0.d0
            do i = 1, 3
              vec_ptr(local_id) = vec_ptr(local_id) + &
                              stress_p((local_id - 1)*SIX_INTEGER + i)
            enddo
            vec_ptr(local_id) = vec_ptr(local_id) / 3.d0
          enddo
          call VecRestoreArray(pmc%subsurf_realization%field%work, &
                                  vec_ptr,ierr);CHKERRQ(ierr)
          call DiscretizationGlobalToLocal( &
                 pmc%subsurf_realization%discretization, &
                 pmc%subsurf_realization%field%work, &
                 pmc%subsurf_realization%field%work_loc,ONEDOF)
          call VecGetArray(pmc%subsurf_realization%field%work_loc, &
                              vec_ptr,ierr);CHKERRQ(ierr)
          parameter_index = ParameterGetIDFromName('geomechanics_stress', &
                                                   option)
          do ghosted_id = 1, grid%ngmax
            pmc%subsurf_realization%patch%aux% &
              Global%auxvars(ghosted_id)%parameters(parameter_index) = &
                vec_ptr(ghosted_id)
          enddo
          call VecRestoreArray(pmc%subsurf_realization%field%work_loc, &
                                  vec_ptr,ierr);CHKERRQ(ierr)
          ! strain
          call VecGetArray(pmc%subsurf_realization%field%work, &
                              vec_ptr,ierr);CHKERRQ(ierr)
          do local_id = 1, grid%nlmax
            vec_ptr(local_id) = 0.d0
            do i = 1, 3
              vec_ptr(local_id) = vec_ptr(local_id) + &
                              strain_p((local_id - 1)*SIX_INTEGER + i)
            enddo
            !vec_ptr(local_id) = vec_ptr(local_id) / 3.d0
          enddo
          call VecRestoreArray(pmc%subsurf_realization%field%work, &
                                  vec_ptr,ierr);CHKERRQ(ierr)
          call DiscretizationGlobalToLocal( &
                 pmc%subsurf_realization%discretization, &
                 pmc%subsurf_realization%field%work, &
                 pmc%subsurf_realization%field%work_loc,ONEDOF)
          call VecGetArray(pmc%subsurf_realization%field%work_loc, &
                              vec_ptr,ierr);CHKERRQ(ierr)
          parameter_index = ParameterGetIDFromName('geomechanics_strain', &
                                                   option)
          do ghosted_id = 1, grid%ngmax
            pmc%subsurf_realization%patch%aux% &
              Global%auxvars(ghosted_id)%parameters(parameter_index) = &
                vec_ptr(ghosted_id)
          enddo
          call VecRestoreArray(pmc%subsurf_realization%field%work_loc, &
                                  vec_ptr,ierr);CHKERRQ(ierr)
        endif

        if (option%geomechanics%flow_coupling == &
            GEOMECH_TWO_WAY_COUPLED) then

          ! Update porosity dataset in sim_aux%subsurf_por
          call VecGetArray(pmc%sim_aux%subsurf_por0,por0_p, &
                              ierr);CHKERRQ(ierr)
          call VecGetArray(pmc%sim_aux%subsurf_por,por_p,ierr);CHKERRQ(ierr)
          ! Perm
          call VecGetArray(pmc%sim_aux%subsurf_perm0,perm0_p, &
                              ierr);CHKERRQ(ierr)
          call VecGetArray(pmc%sim_aux%subsurf_perm,perm_p, &
                              ierr);CHKERRQ(ierr)
          ! Flow
          call VecGetArray(pmc%subsurf_realization%field%flow_xx,press_p, &
                              ierr);CHKERRQ(ierr)

          ! Hoist constant parameter ID lookups outside the element loop
          if (option%geomechanics%split_scheme == &
              GEOMECH_FIXED_STRESS_SPLIT) then
            strainv_0_id = ParameterGetIDFromName('vol_strain_0',option)
            strainv_id = ParameterGetIDFromName('vol_strain', option)
            press_0_id = ParameterGetIDFromName('press_0', option)
            flow_porosity_id = ParameterGetIDFromName('flow_porosity',option)
            temp_0_id = ParameterGetIDFromName('temp_0',option)
            id_porosity_mech=ParameterGetIDFromName('stored_porosity',option)
            id_pressure_mech=ParameterGetIDFromName('stored_pressure', option)
          endif

          do local_id = 1, grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            do i = 1, SIX_INTEGER
              local_stress(i) = stress_p((local_id - 1)*SIX_INTEGER + i)
              local_strain(i) = strain_p((local_id - 1)*SIX_INTEGER + i)
            enddo
            if (option%iflowmode == TH_MODE) then
              local_pressure = press_p(local_id*option%nflowdof-1)
              local_temp = press_p(local_id*option%nflowdof)
            else ! richards mode
              local_pressure = press_p(local_id)
            endif
            ! Update porosity based on stress/strain
            call GeomechanicsSubsurfacePropsPoroEvaluate( &
                  grid, &
                  pmc%subsurf_realization%patch%aux%Material%auxvars(ghosted_id), &
                  por0_p(local_id),local_stress,local_strain,local_pressure, &
                  por_new,option)
            por_p(local_id) = por_new
            ! Update permeability based on stress/strain
            call GeomechanicsSubsurfacePropsPermEvaluate( &
                  grid, &
                  pmc%subsurf_realization%patch%aux%Material%auxvars(ghosted_id), &
                  perm0_p(local_id),local_stress,local_strain,local_pressure, &
                  perm_new,option)
            perm_p(local_id) = perm_new

            select case(option%geomechanics%split_scheme)
            case(GEOMECH_FIXED_STRESS_SPLIT)
              ! ID's are registered in factory_subsurface.F90
              ! (lookups hoisted before loop)
              global_auxvar => global_auxvars(ghosted_id)
              imat = pmc%subsurf_realization%patch%aux%Material% &
                     auxvars(ghosted_id)%id
              ! calc volumetric strain
              strain_vol = local_strain(1) + local_strain(2) + local_strain(3)

              if (option%geomechanics%initial_flag) then ! part of pre-processing
                global_auxvar%parameters(press_0_id) = local_pressure
                global_auxvar%parameters(flow_porosity_id) = &
                 pmc%subsurf_realization%patch%aux%Material% &
                  auxvars(ghosted_id)%porosity_0
                global_auxvar%parameters(strainv_0_id) = 0.d0
                global_auxvar%parameters(strainv_id) = 0.d0
                global_auxvar%parameters(id_porosity_mech) = &
                  pmc%subsurf_realization%patch%aux%Material% &
                  auxvars(ghosted_id)%porosity_0

              else ! after initial solve
                if (this%timestepper%steps == 0) then ! at init
                  global_auxvar%parameters(strainv_0_id) = strain_vol
                  global_auxvar%parameters(strainv_id) = strain_vol
                  global_auxvar%parameters(flow_porosity_id) = &
                    pmc%subsurf_realization%patch%aux%Material% &
                    auxvars(ghosted_id)%porosity_0
                  global_auxvar%parameters(id_porosity_mech) = &
                    pmc%subsurf_realization%patch%aux%Material% &
                    auxvars(ghosted_id)%porosity_0
                  if (this%option%iflowmode == TH_MODE) then
                    global_auxvar%parameters(temp_0_id) = local_temp
                  endif

                else ! for steps >= 1
                  global_auxvar%parameters(strainv_id) = strain_vol
                  global_auxvar%parameters(id_porosity_mech) = &
                  global_auxvar%parameters(flow_porosity_id)
                endif
              endif

              global_auxvar%parameters(id_pressure_mech) = &
                local_pressure
            end select
          enddo

          call VecRestoreArray(pmc%sim_aux%subsurf_por0,por0_p, &
                                  ierr);CHKERRQ(ierr)
          call VecRestoreArray(pmc%sim_aux%subsurf_por,por_p, &
                                  ierr);CHKERRQ(ierr)
          call VecRestoreArray(pmc%subsurf_realization%field%flow_xx,press_p, &
                                  ierr);CHKERRQ(ierr)

          call VecRestoreArray(pmc%sim_aux%subsurf_perm0,perm0_p, &
                                  ierr);CHKERRQ(ierr)
          call VecRestoreArray(pmc%sim_aux%subsurf_perm,perm_p, &
                                  ierr);CHKERRQ(ierr)
        endif

        call VecRestoreArray(pmc%sim_aux%subsurf_stress,stress_p, &
                                ierr);CHKERRQ(ierr)
        call VecRestoreArray(pmc%sim_aux%subsurf_strain,strain_p, &
                                ierr);CHKERRQ(ierr)

        call VecDestroy(geomech_vec,ierr);CHKERRQ(ierr)
        call VecDestroy(subsurf_vec,ierr);CHKERRQ(ierr)

      endif

  end select

end subroutine PMCGeomechanicsSetAuxData

! ************************************************************************** !

subroutine PMCGeomechanicsGetAuxData(this)
  !
  ! This routine updates data for geomechanics simulation from other process
  ! models.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !

  use Option_module
  use Geomechanics_Discretization_module
  use Geomechanics_Force_module

  implicit none

  class(pmc_geomechanics_type) :: this

  PetscErrorCode :: ierr

#if GEOMECH_DEBUG
print *, 'PMCGeomechanicsGetAuxData'
#endif

  select type(pmc => this)
    class is(pmc_geomechanics_type)

      call VecScatterBegin(pmc%sim_aux%subsurf_to_geomechanics, &
                           pmc%sim_aux%subsurf_pres, &
                           pmc%geomech_realization%geomech_field%press, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      call VecScatterEnd(pmc%sim_aux%subsurf_to_geomechanics, &
                         pmc%sim_aux%subsurf_pres, &
                         pmc%geomech_realization%geomech_field%press, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

      call VecScatterBegin(pmc%sim_aux%subsurf_to_geomechanics, &
                           pmc%sim_aux%subsurf_temp, &
                           pmc%geomech_realization%geomech_field%temp, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      call VecScatterEnd(pmc%sim_aux%subsurf_to_geomechanics, &
                         pmc%sim_aux%subsurf_temp, &
                         pmc%geomech_realization%geomech_field%temp, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

      call VecScatterBegin(pmc%sim_aux%subsurf_to_geomechanics, &
                      pmc%sim_aux%subsurf_fluid_den, &
                      pmc%geomech_realization%geomech_field%fluid_density, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      call VecScatterEnd(pmc%sim_aux%subsurf_to_geomechanics, &
                      pmc%sim_aux%subsurf_fluid_den, &
                      pmc%geomech_realization%geomech_field%fluid_density, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

      call VecScatterBegin(pmc%sim_aux%subsurf_to_geomechanics, &
                      pmc%sim_aux%subsurf_por, &
                      pmc%geomech_realization%geomech_field%porosity, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      call VecScatterEnd(pmc%sim_aux%subsurf_to_geomechanics, &
                      pmc%sim_aux%subsurf_por, &
                      pmc%geomech_realization%geomech_field%porosity, &
                      INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

      call GeomechDiscretizationGlobalToLocal( &
                            pmc%geomech_realization%geomech_discretization, &
                            pmc%geomech_realization%geomech_field%press, &
                            pmc%geomech_realization%geomech_field%press_loc, &
                            ONEDOF)

      call GeomechDiscretizationGlobalToLocal( &
                            pmc%geomech_realization%geomech_discretization, &
                            pmc%geomech_realization%geomech_field%temp, &
                            pmc%geomech_realization%geomech_field%temp_loc, &
                            ONEDOF)

      call GeomechDiscretizationGlobalToLocal( &
                    pmc%geomech_realization%geomech_discretization, &
                    pmc%geomech_realization%geomech_field%fluid_density, &
                    pmc%geomech_realization%geomech_field%fluid_density_loc, &
                    ONEDOF)

      call GeomechDiscretizationGlobalToLocal( &
                    pmc%geomech_realization%geomech_discretization, &
                    pmc%geomech_realization%geomech_field%porosity, &
                    pmc%geomech_realization%geomech_field%porosity_loc, &
                    ONEDOF)

  end select

end subroutine PMCGeomechanicsGetAuxData

! ************************************************************************** !

recursive subroutine PMCGeomechanicsCheckpointBinary(this,viewer,append_name)
  !
  ! Checkpoints both flow and geomechanics PMC state to a binary file.
  ! Overrides PMCBaseCheckpointBinary to handle the ordering problem:
  ! geomech is the list head but is_master=FALSE, so the base class would
  ! try to write geomech data before the file is opened by the master (flow).
  ! This override opens the file, writes flow data first (preserving
  ! flow-as-master ordering), then writes geomech data, then closes the file.
  !
  ! Author: Satish Karra, PNNL
  ! Date: 07/09/2025
  !
#include <petsc/finclude/petscbag.h>
  use petscbag

  use Logging_module
  use Checkpoint_module, only : CheckpointOpenFileForWriteBinary, &
                                CheckPointWriteCompatibilityBinary
  use PM_Base_class
  use Option_module, only : PrintMsg

  implicit none

  class(pmc_geomechanics_type) :: this
  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: append_name

  class(pm_base_type), pointer :: cur_pm
  class(pmc_base_header_type), pointer :: header
  type(pmc_base_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize
  PetscLogDouble :: tstart, tend
  PetscErrorCode :: ierr

  bagsize = size(transfer(dummy_header,dummy_char))

  ! --- Open file (taking over the role of is_master) ---
  call PetscLogStagePush(logging%stage(OUTPUT_STAGE),ierr);CHKERRQ(ierr)
  call PetscLogEventBegin(logging%event_checkpoint,ierr);CHKERRQ(ierr)
  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  call CheckpointOpenFileForWriteBinary(viewer,append_name,this%option)
  call CheckPointWriteCompatibilityBinary(viewer,this%option)

  ! --- Write PMC header using flow (peer) PMC's output info ---
  ! Use flow PMC's header info since flow is the conceptual master.
  call PetscBagCreate(this%option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call PMCBaseRegisterHeader(this%peer,bag,header)
  call PMCBaseSetHeader(this%peer,bag,header)
  call PetscBagView(bag,viewer,ierr);CHKERRQ(ierr)
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

  ! --- Write flow (peer) PMC data first ---
  ! Flow timestepper
  if (associated(this%peer%timestepper)) then
    call this%peer%timestepper%CheckpointBinary(viewer,this%option)
  endif
  ! Flow PM(s)
  cur_pm => this%peer%pm_list
  do
    if (.not.associated(cur_pm)) exit
    call cur_pm%CheckpointBinary(viewer)
    cur_pm => cur_pm%next
  enddo
  ! Flow children (if any)
  if (associated(this%peer%child)) then
    call this%peer%child%CheckpointBinary(viewer,append_name)
  endif

  ! --- Write geomechanics (this) PMC data second ---
  ! Geomech timestepper
  if (associated(this%timestepper)) then
    call this%timestepper%CheckpointBinary(viewer,this%option)
  endif
  ! Geomech PM(s)
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    call cur_pm%CheckpointBinary(viewer)
    cur_pm => cur_pm%next
  enddo
  ! Geomech children (if any)
  if (associated(this%child)) then
    call this%child%CheckpointBinary(viewer,append_name)
  endif

  ! --- Close file ---
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  PetscObjectNullify(viewer)
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(this%option%io_buffer, &
        '(6x,"Seconds to write to checkpoint file: ", f10.2)') &
    tend-tstart
  call PrintMsg(this%option)
  call PetscLogEventEnd(logging%event_checkpoint,ierr);CHKERRQ(ierr)
  call PetscLogStagePop(ierr);CHKERRQ(ierr)

end subroutine PMCGeomechanicsCheckpointBinary

! ************************************************************************** !

recursive subroutine PMCGeomechanicsRestartBinary(this,viewer)
  !
  ! Restarts both flow and geomechanics PMC state from a binary file.
  ! Overrides PMCBaseRestartBinary to handle the ordering problem.
  ! Opens the file, reads flow data first (matching checkpoint order),
  ! then reads geomech data, then closes the file.
  !
  ! Author: Satish Karra, PNNL
  ! Date: 07/09/2025
  !
#include <petsc/finclude/petscbag.h>
  use petscbag

  use Logging_module
  use Checkpoint_module, only : CheckPointReadCompatibilityBinary
  use Petsc_Utility_module, only : PUTestFile
  use PM_Base_class
  use Option_module, only : PrintMsg, PrintErrMsg
  use Waypoint_module, only : WaypointSkipToTime

  implicit none

  class(pmc_geomechanics_type) :: this
  PetscViewer :: viewer

  class(pm_base_type), pointer :: cur_pm
  class(pmc_base_header_type), pointer :: header
  type(pmc_base_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize
  PetscLogDouble :: tstart, tend
  PetscBool :: flag
  PetscErrorCode :: ierr

  bagsize = size(transfer(dummy_header,dummy_char))

  ! --- Open file ---
  call PUTestFile(this%option%restart_filename,'r',flag, &
                  ierr);CHKERRQ(ierr)
  if (.not.flag) then
    this%option%io_buffer = 'Restart file "' // &
      trim(this%option%restart_filename) // '" not found.'
    call PrintErrMsg(this%option)
  endif
  this%option%io_buffer = 'Restarting with checkpoint file "' // &
    trim(this%option%restart_filename) // '".'
  call PrintMsg(this%option)
  call PetscLogEventBegin(logging%event_restart,ierr);CHKERRQ(ierr)
  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  call PetscViewerBinaryOpen(this%option%mycomm, &
                             this%option%restart_filename,FILE_MODE_READ, &
                             viewer,ierr);CHKERRQ(ierr)
  call PetscViewerBinarySetSkipOptions(viewer,PETSC_TRUE, &
                                       ierr);CHKERRQ(ierr)
  call CheckPointReadCompatibilityBinary(viewer,this%option)

  ! --- Read PMC header (flow's header, since flow wrote it) ---
  call PetscBagCreate(this%option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call PMCBaseRegisterHeader(this%peer,bag,header)
  call PetscBagLoad(viewer,bag,ierr);CHKERRQ(ierr)
  call PMCBaseGetHeader(this%peer,header)
  if (Initialized(this%option%restart_time)) then
    this%peer%pm_list%realization_base%output_option%plot_number = 0
  endif
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

  ! --- Read flow (peer) PMC data first ---
  ! Flow timestepper
  if (associated(this%peer%timestepper)) then
    call this%peer%timestepper%RestartBinary(viewer,this%option)
    if (Initialized(this%option%restart_time)) then
      call this%peer%timestepper%Reset()
    endif
    call WaypointSkipToTime(this%peer%timestepper%cur_waypoint, &
                            this%peer%timestepper%target_time)
    if (.not.associated(this%peer%timestepper%cur_waypoint)) then
      write(this%option%io_buffer,*) this%peer%timestepper%target_time* &
        this%peer%pm_list%realization_base%output_option%tconv
      this%option%io_buffer = 'Simulation is being restarted at a time that &
        &is at or beyond the end of checkpointed simulation (' // &
        trim(adjustl(this%option%io_buffer)) // &
        trim(this%peer%pm_list%realization_base%output_option%tunit) // ').'
      call PrintMsg(this%option)
    endif
    ! dt_max is not saved in the checkpoint header, so restore it from the
    ! current waypoint to avoid UNINITIALIZED_DOUBLE (-999) corrupting dt
    ! when revert_dt is TRUE in SetTargetTime.
    if (associated(this%peer%timestepper%cur_waypoint)) then
      this%peer%timestepper%dt_max = &
        this%peer%timestepper%cur_waypoint%dt_max
    endif
    this%option%time = this%peer%timestepper%target_time
  endif
  ! Flow PM(s)
  cur_pm => this%peer%pm_list
  do
    if (.not.associated(cur_pm)) exit
    if (cur_pm%skip_restart) then
      this%option%io_buffer = 'Due to sequential nature of binary files, &
        &skipping restart for binary formatted files is not allowed.'
      call PrintErrMsg(this%option)
    endif
    call cur_pm%RestartBinary(viewer)
    cur_pm => cur_pm%next
  enddo
  ! Flow children (if any)
  if (associated(this%peer%child)) then
    call this%peer%child%RestartBinary(viewer)
  endif

  ! --- Read geomechanics (this) PMC data second ---
  ! Geomech timestepper
  if (associated(this%timestepper)) then
    call this%timestepper%RestartBinary(viewer,this%option)
    if (Initialized(this%option%restart_time)) then
      call this%timestepper%Reset()
    endif
    call WaypointSkipToTime(this%timestepper%cur_waypoint, &
                            this%timestepper%target_time)
    if (.not.associated(this%timestepper%cur_waypoint)) then
      write(this%option%io_buffer,*) this%timestepper%target_time* &
        this%pm_list%realization_base%output_option%tconv
      this%option%io_buffer = 'Simulation is being restarted at a time that &
        &is at or beyond the end of checkpointed simulation (' // &
        trim(adjustl(this%option%io_buffer)) // &
        trim(this%pm_list%realization_base%output_option%tunit) // ').'
      call PrintMsg(this%option)
    endif
    ! dt_max is not saved in the checkpoint header, so restore it from the
    ! current waypoint to avoid UNINITIALIZED_DOUBLE (-999) corrupting dt
    ! when revert_dt is TRUE in SetTargetTime.
    if (associated(this%timestepper%cur_waypoint)) then
      this%timestepper%dt_max = this%timestepper%cur_waypoint%dt_max
    endif
    this%option%time = this%timestepper%target_time
  endif
  ! Geomech PM(s)
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    if (cur_pm%skip_restart) then
      this%option%io_buffer = 'Due to sequential nature of binary files, &
        &skipping restart for binary formatted files is not allowed.'
      call PrintErrMsg(this%option)
    endif
    call cur_pm%RestartBinary(viewer)
    cur_pm => cur_pm%next
  enddo
  ! Geomech children (if any)
  if (associated(this%child)) then
    call this%child%RestartBinary(viewer)
  endif

  ! --- Close file ---
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(this%option%io_buffer, &
        '("      Seconds to read from restart file: ", f10.2)') &
    tend-tstart
  call PrintMsg(this%option)
  call PetscLogEventEnd(logging%event_restart,ierr);CHKERRQ(ierr)

end subroutine PMCGeomechanicsRestartBinary

! ************************************************************************** !

recursive subroutine PMCGeomechanicsCheckpointHDF5(this,h5_chk_grp_id, &
                                                    append_name)
  !
  ! Checkpoints both flow and geomechanics PMC state to an HDF5 file.
  ! Overrides PMCBaseCheckpointHDF5 to handle the ordering problem.
  ! Opens the file, writes flow data first, then geomech data.
  !
  ! Author: Satish Karra, PNNL
  ! Date: 07/09/2025
  !

  use Logging_module
  use Checkpoint_module, only : CheckpointOpenFileForWriteHDF5, &
                                CheckPointWriteCompatibilityHDF5
  use PM_Base_class
  use hdf5
  use HDF5_Aux_module
  use Option_module, only : PrintMsg

  implicit none

  class(pmc_geomechanics_type) :: this
  integer(HID_T) :: h5_chk_grp_id
  character(len=MAXSTRINGLENGTH) :: append_name

  integer(HID_T) :: h5_file_id
  integer(HID_T) :: h5_pmc_grp_id
  integer(HID_T) :: h5_pm_grp_id

  class(pm_base_type), pointer :: cur_pm
  PetscLogDouble :: tstart, tend
  PetscErrorCode :: ierr

  ! --- Open file ---
  call PetscLogStagePush(logging%stage(OUTPUT_STAGE),ierr);CHKERRQ(ierr)
  call PetscLogEventBegin(logging%event_checkpoint,ierr);CHKERRQ(ierr)
  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  call CheckpointOpenFileForWriteHDF5(h5_file_id, &
                                      h5_chk_grp_id, &
                                      append_name,this%option)
  call CheckPointWriteCompatibilityHDF5(h5_chk_grp_id, &
                                        this%option)

  ! --- Write flow (peer) PMC data first ---
  call HDF5GroupCreate(h5_chk_grp_id,trim(this%peer%name),h5_pmc_grp_id, &
                       this%option)
  call PMCBaseSetHeaderHDF5(this%peer,h5_pmc_grp_id,this%option)
  if (associated(this%peer%timestepper)) then
    call this%peer%timestepper%CheckpointHDF5(h5_pmc_grp_id,this%option)
  endif
  cur_pm => this%peer%pm_list
  do
    if (.not.associated(cur_pm)) exit
    call HDF5GroupCreate(h5_pmc_grp_id,trim(cur_pm%name),h5_pm_grp_id, &
                         this%option)
    call cur_pm%CheckpointHDF5(h5_pm_grp_id)
    call HDF5GroupClose(h5_pm_grp_id,this%option)
    cur_pm => cur_pm%next
  enddo
  call HDF5GroupClose(h5_pmc_grp_id,this%option)
  ! Flow children (if any)
  if (associated(this%peer%child)) then
    call this%peer%child%CheckpointHDF5(h5_chk_grp_id,append_name)
  endif

  ! --- Write geomechanics (this) PMC data second ---
  call HDF5GroupCreate(h5_chk_grp_id,trim(this%name),h5_pmc_grp_id, &
                       this%option)
  if (associated(this%timestepper)) then
    call this%timestepper%CheckpointHDF5(h5_pmc_grp_id,this%option)
  endif
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    call HDF5GroupCreate(h5_pmc_grp_id,trim(cur_pm%name),h5_pm_grp_id, &
                         this%option)
    call cur_pm%CheckpointHDF5(h5_pm_grp_id)
    call HDF5GroupClose(h5_pm_grp_id,this%option)
    cur_pm => cur_pm%next
  enddo
  call HDF5GroupClose(h5_pmc_grp_id,this%option)
  ! Geomech children (if any)
  if (associated(this%child)) then
    call this%child%CheckpointHDF5(h5_chk_grp_id,append_name)
  endif

  ! --- Close file ---
  call HDF5GroupClose(h5_chk_grp_id,this%option)
  call HDF5FileClose(h5_file_id,this%option)
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(this%option%io_buffer, &
        '("      Seconds to write to checkpoint file: ", f10.2)') &
    tend-tstart
  call PrintMsg(this%option)
  call PetscLogEventEnd(logging%event_checkpoint,ierr);CHKERRQ(ierr)
  call PetscLogStagePop(ierr);CHKERRQ(ierr)

end subroutine PMCGeomechanicsCheckpointHDF5

! ************************************************************************** !

recursive subroutine PMCGeomechanicsRestartHDF5(this,h5_chk_grp_id)
  !
  ! Restarts both flow and geomechanics PMC state from an HDF5 file.
  ! Overrides PMCBaseRestartHDF5 to handle the ordering problem.
  ! Opens the file, reads flow data first, then geomech data.
  !
  ! Author: Satish Karra, PNNL
  ! Date: 07/09/2025
  !

  use Logging_module
  use PM_Base_class
  use hdf5
  use HDF5_Aux_module
  use Checkpoint_module, only : CheckPointReadCompatibilityHDF5, &
                                CheckpointOpenFileForReadHDF5
  use Option_module, only : PrintMsg, PrintErrMsg
  use Waypoint_module, only : WaypointSkipToTime

  implicit none

  class(pmc_geomechanics_type) :: this
  integer(HID_T) :: h5_chk_grp_id

  class(pm_base_type), pointer :: cur_pm
  PetscLogDouble :: tstart, tend
  PetscErrorCode :: ierr

  integer(HID_T) :: h5_file_id
  integer(HID_T) :: h5_pmc_grp_id
  integer(HID_T) :: h5_pm_grp_id

  ! --- Open file ---
  this%option%io_buffer = 'Restarting with checkpoint file: ' // &
    trim(this%option%restart_filename)
  call PrintMsg(this%option)
  call PetscLogEventBegin(logging%event_restart,ierr);CHKERRQ(ierr)
  call PetscTime(tstart,ierr);CHKERRQ(ierr)

  call CheckpointOpenFileForReadHDF5(this%option%restart_filename, &
                                     h5_file_id, &
                                     h5_chk_grp_id, &
                                     this%option)

  call CheckPointReadCompatibilityHDF5(h5_chk_grp_id, &
                                       this%option)

  ! --- Read flow (peer) PMC data first ---
  call HDF5GroupOpen(h5_chk_grp_id,this%peer%name,h5_pmc_grp_id, &
                     this%option)
  call PMCBaseGetHeaderHDF5(this%peer,h5_pmc_grp_id,this%option)
  if (Initialized(this%option%restart_time)) then
    this%peer%pm_list%realization_base%output_option%plot_number = 0
  endif
  if (associated(this%peer%timestepper)) then
    call this%peer%timestepper%RestartHDF5(h5_pmc_grp_id,this%option)
    if (Initialized(this%option%restart_time)) then
      call this%peer%timestepper%Reset()
    endif
    call WaypointSkipToTime(this%peer%timestepper%cur_waypoint, &
                            this%peer%timestepper%target_time)
    if (.not.associated(this%peer%timestepper%cur_waypoint)) then
      write(this%option%io_buffer,*) this%peer%timestepper%target_time/ &
        this%peer%pm_list%realization_base%output_option%tconv
      this%option%io_buffer = 'Simulation is being restarted at a time that &
        &is at or beyond the end of checkpointed simulation (' // &
        trim(adjustl(this%option%io_buffer)) // &
        trim(this%peer%pm_list%realization_base%output_option%tunit) // ').'
      call PrintErrMsg(this%option)
    endif
    ! dt_max is not saved in the checkpoint header, so restore it from the
    ! current waypoint to avoid UNINITIALIZED_DOUBLE (-999) corrupting dt
    ! when revert_dt is TRUE in SetTargetTime.
    if (associated(this%peer%timestepper%cur_waypoint)) then
      this%peer%timestepper%dt_max = &
        this%peer%timestepper%cur_waypoint%dt_max
    endif
    this%option%time = this%peer%timestepper%target_time
  endif
  cur_pm => this%peer%pm_list
  do
    if (.not.associated(cur_pm)) exit
    call HDF5GroupOpen(h5_pmc_grp_id,cur_pm%name,h5_pm_grp_id, &
                       this%option)
    call cur_pm%RestartHDF5(h5_pm_grp_id)
    call HDF5GroupClose(h5_pm_grp_id,this%option)
    call this%peer%SetAuxData()
    cur_pm => cur_pm%next
  enddo
  call HDF5GroupClose(h5_pmc_grp_id,this%option)
  ! Flow children (if any)
  if (associated(this%peer%child)) then
    call this%peer%child%RestartHDF5(h5_chk_grp_id)
  endif

  ! --- Read geomechanics (this) PMC data second ---
  call HDF5GroupOpen(h5_chk_grp_id,this%name,h5_pmc_grp_id, &
                     this%option)
  if (associated(this%timestepper)) then
    call this%timestepper%RestartHDF5(h5_pmc_grp_id,this%option)
    if (Initialized(this%option%restart_time)) then
      call this%timestepper%Reset()
    endif
    call WaypointSkipToTime(this%timestepper%cur_waypoint, &
                            this%timestepper%target_time)
    if (.not.associated(this%timestepper%cur_waypoint)) then
      write(this%option%io_buffer,*) this%timestepper%target_time/ &
        this%pm_list%realization_base%output_option%tconv
      this%option%io_buffer = 'Simulation is being restarted at a time that &
        &is at or beyond the end of checkpointed simulation (' // &
        trim(adjustl(this%option%io_buffer)) // &
        trim(this%pm_list%realization_base%output_option%tunit) // ').'
      call PrintErrMsg(this%option)
    endif
    ! dt_max is not saved in the checkpoint header, so restore it from the
    ! current waypoint to avoid UNINITIALIZED_DOUBLE (-999) corrupting dt
    ! when revert_dt is TRUE in SetTargetTime.
    if (associated(this%timestepper%cur_waypoint)) then
      this%timestepper%dt_max = this%timestepper%cur_waypoint%dt_max
    endif
    this%option%time = this%timestepper%target_time
  endif
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    call HDF5GroupOpen(h5_pmc_grp_id,cur_pm%name,h5_pm_grp_id, &
                       this%option)
    call cur_pm%RestartHDF5(h5_pm_grp_id)
    call HDF5GroupClose(h5_pm_grp_id,this%option)
    call this%SetAuxData()
    cur_pm => cur_pm%next
  enddo
  call HDF5GroupClose(h5_pmc_grp_id,this%option)
  ! Geomech children (if any)
  if (associated(this%child)) then
    call this%child%RestartHDF5(h5_chk_grp_id)
  endif

  ! --- Close file ---
  call HDF5GroupClose(h5_chk_grp_id,this%option)
  call HDF5FileClose(h5_file_id,this%option)
  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(this%option%io_buffer, &
        '("      Seconds to read from restart file: ", f10.2)') &
    tend-tstart
  call PrintMsg(this%option)
  call PetscLogEventEnd(logging%event_restart,ierr);CHKERRQ(ierr)

end subroutine PMCGeomechanicsRestartHDF5

! ************************************************************************** !

subroutine PMCGeomechanicsStrip(this)
  !
  ! Deallocates members of PMC Geomechanics.
  !
  ! Author: Satish Karra
  ! Date: 06/01/16

  implicit none

  class(pmc_geomechanics_type) :: this

  call PMCBaseStrip(this)
  ! realizations destroyed elsewhere
  nullify(this%subsurf_realization)
  nullify(this%geomech_realization)

end subroutine PMCGeomechanicsStrip

! ************************************************************************** !

recursive subroutine PMCGeomechanicsDestroy(this)
  !
  ! Author: Satish Karra
  ! Date: 06/01/16
  !
  use Option_module

  implicit none

  class(pmc_geomechanics_type) :: this

#ifdef DEBUG
  call PrintMsg(this%option,'PMCGeomechanics%Destroy()')
#endif

  if (associated(this%child)) then
    call this%child%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%child)
    nullify(this%child)
  endif

  if (associated(this%peer)) then
    call this%peer%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%peer)
    nullify(this%peer)
  endif

  call PMCGeomechanicsStrip(this)

end subroutine PMCGeomechanicsDestroy

end module PMC_Geomechanics_class
