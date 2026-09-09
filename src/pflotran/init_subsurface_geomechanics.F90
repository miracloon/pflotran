module Init_Subsurface_Geomech_module

#include "petsc/finclude/petscmat.h"
  use petscmat

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: InitSubsurfGeomechReadRequiredCards, &
            InitSubsurfGeomechReadInput, &
            InitSubsurfGeomechJumpStart, & ! remove later
            InitSubsurfGeomechSetupRealization, &
            InitSubsurfGeomechInitSimulation, &
            InitSubsurfGeomechSetGeomechMode, &
            InitSubsurfGeomechChkInactiveCells, &
            InitSubsurfGeomechSetupPMC, &
            InitSubsurfGeomechReadSimBlock
contains

! ************************************************************************** !

subroutine InitSubsurfGeomechReadRequiredCards(geomech_realization,input)
  !
  ! Reads the required input file cards
  ! related to geomechanics
  !
  ! Author: Satish Karra, LANL
  ! Date: 05/23/13
  !
  ! jaa: moved from factory_geomechanics.F90 on 1/28/25

  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Input_Aux_module
  use String_module
  use Patch_module
  use Option_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  type(input_type), pointer :: input

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option

  option => geomech_realization%option

! Read in select required cards
!.........................................................................

  ! GEOMECHANICS information
  string = "GEOMECHANICS"
  call InputFindStringInFile(input,option,string)
  if (InputError(input)) return

  string = "GEOMECHANICS_GRID"
  call InputFindStringInFile(input,option,string)
  call InitSubsurfGeomechReadGridBlock(geomech_realization,input,option)

end subroutine InitSubsurfGeomechReadRequiredCards

! ************************************************************************** !

subroutine InitSubsurfGeomechReadInput(geomech,geomech_solver, &
                                     input,option,output_option)
  !
  ! Reads the geomechanics input data
  !
  ! Author: Satish Karra, LANL
  ! Date: 05/23/13
  !
  ! jaa: moved from factory_geomechanics.F90 on 1/28/25

  use Option_module
  use Input_Aux_module
  use String_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Material_module
  use Geomechanics_Region_module
  use Geomechanics_Debug_module
  use Geomechanics_Strata_module
  use Geomechanics_Condition_module
  use Geomechanics_Coupler_module
  use Geomechanics_Regression_module
  use Option_Geomechanics_module
  use Output_Aux_module
  use Output_Tecplot_module
  use Realization_Base_class
  use Solver_module
  use Units_module
  use Waypoint_module
  use Dataset_Base_class
  use Dataset_module
  use Dataset_Common_HDF5_class
  use Utility_module, only : DeallocateArray, UtilityReadArray
  use Geomechanics_Attr_module

  ! Still need to add other geomech modules for output, etc once created

  implicit none

  type(solver_type), pointer :: geomech_solver
  type(input_type), pointer :: input
  type(geomechanics_attr_type), pointer:: geomech
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option

  class(realization_geomech_type), pointer :: geomech_realization
  class(dataset_base_type), pointer :: dataset

  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_material_property_type),pointer :: geomech_material_property
  type(geomech_grid_type), pointer :: grid
  type(gm_region_type), pointer :: region
  type(geomech_strata_type), pointer :: strata
  type(geomech_condition_type), pointer :: condition
  type(geomech_coupler_type), pointer :: coupler
  type(geomech_observation_type), pointer :: geomech_observation
  type(waypoint_list_type), pointer :: waypoint_list
  type(waypoint_list_type), pointer :: waypoint_list_dt_coupling
  type(waypoint_type), pointer :: waypoint

  character(len=MAXWORDLENGTH) :: word, internal_units
  character(len=MAXWORDLENGTH) :: card
  character(len=1) :: backslash
  PetscReal :: temp_real

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  input%ierr = 0
! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  waypoint_list => geomech%waypoint_list
  waypoint_list_dt_coupling => geomech%waypoint_list_dt_coupling
  geomech_realization => geomech%realization
  geomech_discretization => geomech_realization%geomech_discretization
  geomech_realization%body_force(:) = option%geomechanics%gravity(:)

  if (associated(geomech_realization%geomech_patch)) grid => &
    geomech_realization%geomech_patch%geomech_grid

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS')
    call StringToUpper(word)
    option%io_buffer = 'word :: ' // trim(word)
    call PrintMsg(option)

    select case(trim(word))

      !.........................................................................
      ! Read geomechanics grid information
      case ('GEOMECHANICS_GRID')
        call InputSkipToEND(input,option,trim(word))

      !.........................................................................
      ! Read geomechanics material information
      case ('GEOMECHANICS_MATERIAL_PROPERTY')
        geomech_material_property => GeomechanicsMaterialPropertyCreate()

        call InputReadWord(input,option,geomech_material_property%name, &
                           PETSC_TRUE)

        call InputErrorMsg(input,option,'name','GEOMECHANICS_MATERIAL_PROPERTY')
        call GeomechanicsMaterialPropertyRead(geomech_material_property,input, &
                                              option)
        call GeomechanicsMaterialPropertyAddToList(geomech_material_property, &
                                geomech_realization%geomech_material_properties)
        nullify(geomech_material_property)

      case ('GEOMECHANICS_SET_REF_P_T_TO_IC')
        option%geomechanics%set_ref_pres_and_temp_to_IC = PETSC_TRUE

      !.........................................................................
      ! Read geomechanics datasets
      case ('GEOMECHANICS_DATASET')
          nullify(dataset)
          call DatasetRead(input,dataset,option)
          call DatasetBaseAddToList(dataset,geomech_realization%geomech_datasets)
          nullify(dataset)

      !.........................................................................
      case ('GEOMECHANICS_BODY_FORCE')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','GEOMECHANICS_BODY_FORCE')
          call StringToUpper(word)
          select case(trim(word))
            case('X','BODY_FORCE_X')
              call DatasetReadDoubleorDataset(input, &
                        geomech_realization%body_force(X_DIRECTION), &
                        geomech_realization%body_force_x_dataset, &
                        'X','GEOMECHANICS_BODY_FORCE',option)
              call InputErrorMsg(input,option,'X','GEOMECHANICS_BODY_FORCE')
            case('Y','BODY_FORCE_Y')
              call DatasetReadDoubleorDataset(input, &
                        geomech_realization%body_force(Y_DIRECTION), &
                        geomech_realization%body_force_y_dataset, &
                        'Y','GEOMECHANICS_BODY_FORCE',option)
              call InputErrorMsg(input,option,'Y','GEOMECHANICS_BODY_FORCE')
            case('Z','BODY_FORCE_Z')
              call DatasetReadDoubleorDataset(input, &
                        geomech_realization%body_force(Z_DIRECTION), &
                        geomech_realization%body_force_z_dataset, &
                        'Z','GEOMECHANICS_BODY_FORCE',option)
              call InputErrorMsg(input,option,'Z','GEOMECHANICS_BODY_FORCE')
            case default
              call InputKeywordUnrecognized(input,word, &
                                            'GEOMECHANICS_BODY_FORCE',option)
          end select
        enddo
        call InputPopBlock(input,option)

      !.........................................................................
      case ('GEOMECHANICS_REGION')
        region => GeomechRegionCreate()
        call InputReadWord(input,option,region%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','GEOMECHANICS_REGION')
        call PrintMsg(option,region%name)
        call GeomechRegionRead(region,input,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call GeomechRegionAddToList(region,geomech_realization%geomech_region_list)
        nullify(region)

      !.........................................................................
      case ('GEOMECHANICS_CONDITION')
        condition => GeomechConditionCreate(option)
        call InputReadWord(input,option,condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'GEOMECHANICS_CONDITION','name')
        call PrintMsg(option,condition%name)
        call GeomechConditionRead(condition,input,option)
        call GeomechConditionAddToList(condition,geomech_realization%geomech_conditions)
        nullify(condition)

     !.........................................................................
      case ('GEOMECHANICS_BOUNDARY_CONDITION')
        coupler =>  GeomechCouplerCreate(GM_BOUNDARY_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Geomech Boundary Condition name')
        call GeomechCouplerRead(coupler,input,option)
        call GeomechRealizAddGeomechCoupler(geomech_realization,coupler)
        nullify(coupler)

      !.........................................................................
      case ('GEOMECHANICS_SRC_SINK')
        coupler => GeomechCouplerCreate(GM_SRC_SINK_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Source Sink name')
        call GeomechCouplerRead(coupler,input,option)
        call GeomechRealizAddGeomechCoupler(geomech_realization,coupler)
        nullify(coupler)

     !....................
      case ('GEOMECHANICS_LINEAR_SOLVER')
        call SolverReadLinear(geomech_solver,input,option)

      !.....................
      case ('GEOMECHANICS_REGRESSION')
        call GeomechanicsRegressionRead(geomech%regression,input,option)

      !.........................................................................
      case ('GEOMECHANICS_TIME')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'word','GEOMECHANICS_TIME')
          select case(trim(word))
            case('COUPLING_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option, &
                                 'Coupling Timestep Size','GEOMECHANICS_TIME')
              internal_units = 'sec'
              call InputReadAndConvertUnits(input, temp_real, &
                                            internal_units,'GEOMECHANICS_TIME,&
                                            &COUPLING_TIMESTEP_SIZE',option)
              geomech_realization%dt_coupling = temp_real
              call InputReadCard(input,option,word)
              if (input%ierr == 0) then
                  call StringToUpper(word)
                  if (StringCompare(word, 'AT', TWO_INTEGER)) then
                      waypoint => WaypointCreate()
                      waypoint%dt_max = temp_real
                      call InputReadDouble(input,option,waypoint%time)
                      call InputErrorMsg(input,option,'COUPLING_TIMESTEP_SIZE &
                                                      &Update Time',card)
                      call InputReadAndConvertUnits(input,waypoint%time, &
                                                    internal_units, &
                                                    'GEOMECHANICS_TIME,COUPLING_TIMESTEP_SIZE,&
                                                    &Update Time',option)
                  else
                      option%io_buffer = 'Keyword under "COUPLING_TIMESTEP_SIZE" &
                                        &after coupling timestep size should be "AT".'
                      call PrintErrMsg(option)
                  endif
                  if (.not.associated(waypoint_list_dt_coupling)) then
                      waypoint_list_dt_coupling => WaypointListCreate()
                  endif
                  call WaypointInsertInList(waypoint, &
                                            waypoint_list_dt_coupling, &
                                            option)
                  if (waypoint_list_dt_coupling%first%time > 0.d0) then
                      option%io_buffer = 'First time after keyword "AT" under &
                                         &"COUPLING_TIMESTEP_SIZE" must be zero (0.d0).'
                      call PrintErrMsg(option)
                  endif
              endif
            case('SYNC_FLOW_TIMESTEP_SIZE') ! for sync drained split
              option%geomechanics%sync_flow_dt = PETSC_TRUE
            case default
              call InputKeywordUnrecognized(input,word, &
                                            'GEOMECHANICS_TIME',option)
            end select
        enddo
        call InputPopBlock(input,option)

      !.........................................................................
      case ('GEOMECHANICS_FLOW_INTERPOLATION')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,word,'GEOMECHANICS_FLOW_INTERPOLATION')
          call StringToUpper(word)
          select case(trim(word))
            case('ORDER')
              call InputReadInt(input,option, &
                                option%geomechanics%flow_interp_order)
              call InputErrorMsg(input,option,'ORDER', &
                                 'GEOMECHANICS_FLOW_INTERPOLATION')
              select case(option%geomechanics%flow_interp_order)
                case(GEOMECH_FLOW_INTERP_ORDER_0TH,GEOMECH_FLOW_INTERP_ORDER_1ST)
                case default
                  option%io_buffer = 'GEOMECHANICS_FLOW_INTERPOLATION ORDER ' // &
                    'must be 0 (0th order) or 1 (1st order).'
                  call PrintErrMsg(option)
              end select
            case default
              call InputKeywordUnrecognized(input,word, &
                                           'GEOMECHANICS_FLOW_INTERPOLATION', &
                                           option)
          end select
        enddo
        call InputPopBlock(input,option)

      !.........................................................................
      case ('GEOMECHANICS_DEBUG')
        call GeomechDebugRead(geomech_realization%geomech_debug,input,option)

      !.........................................................................
      case ('GEOMECHANICS_MAPPING_FILE')
        call InputReadFilename(input,option,grid%mapping_filename)
        call InputErrorMsg(input,option,'keyword','mapping_file')
        call GeomechSubsurfMapFromFilename(grid,grid%mapping_filename,option)

      !.........................................................................
      case ('GEOMECHANICS_OUTPUT')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','GEOMECHANICS_OUTPUT')
          call StringToUpper(word)
          select case(trim(word))
            case('TIMES')
              option%io_buffer = 'Subsurface times are now used for ' // &
              'geomechanics as well. No need for TIMES keyword under ' // &
              'GEOMECHANICS_OUTPUT.'
              call PrintWrnMsg(option)
            case('FORMAT')
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'keyword','GEOMECHANICS_OUTPUT,&
                                                         &FORMAT')
              call StringToUpper(word)
              select case(trim(word))
                case ('HDF5')
                  output_option%print_hdf5 = PETSC_TRUE
                  call InputReadCard(input,option,word)
                  call InputDefaultMsg(input,option, &
                                       'GEOMECHANICS_OUTPUT,FORMAT,HDF5,&
                                        &# FILES')
                  if (len_trim(word) > 1) then
                    call StringToUpper(word)
                    select case(trim(word))
                      case('SINGLE_FILE')
                        output_option%print_single_h5_file = PETSC_TRUE
                      case('MULTIPLE_FILES')
                        output_option%print_single_h5_file = PETSC_FALSE
                      case default
                        option%io_buffer = 'HDF5 keyword (' // trim(word) // &
                          ') not recongnized.  Use "SINGLE_FILE" or ' // &
                          '"MULTIPLE_FILES".'
                        call PrintErrMsg(option)
                    end select
                  endif
                case ('TECPLOT')
                  output_option%print_tecplot = PETSC_TRUE
                  call InputReadCard(input,option,word)
                  call InputErrorMsg(input,option,'TECPLOT','GEOMECHANICS_OUTPUT,FORMAT')
                  call StringToUpper(word)
                  output_option%tecplot_format = TECPLOT_FEQUADRILATERAL_FORMAT ! By default it is unstructured
                case ('VTK')
                  output_option%print_vtk = PETSC_TRUE
                case default
                  call InputKeywordUnrecognized(input,word, &
                                 'GEOMECHANICS_OUTPUT,FORMAT',option)
              end select
            case('PERIODIC_OBSERVATION_TIMESTEP')
              call InputReadInt(input,option, &
                   output_option%periodic_obs_output_ts_imod)
              call InputErrorMsg(input,option,'timestep modulus', &
                   'GEOMECHANICS_OUTPUT,PERIODIC_OBSERVATION_TIMESTEP')
            case default
              call InputKeywordUnrecognized(input,word, &
                             'GEOMECHANICS_OUTPUT',option)
          end select
        enddo
        call InputPopBlock(input,option)

      !.........................................................................
      case ('GEOMECHANICS_STRATIGRAPHY','GEOMECHANICS_STRATA')
        strata => GeomechStrataCreate()
        call GeomechStrataRead(strata,input,option)
        call GeomechRealizAddStrata(geomech_realization,strata)
        nullify(strata)

      !.........................................................................
      case ('GEOMECHANICS_OBSERVATION')
        geomech_observation => GeomechObservationCreate()
        call InputReadWord(input,option,geomech_observation%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','GEOMECHANICS_OBSERVATION')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','GEOMECHANICS_OBSERVATION')
          call StringToUpper(word)
          select case(trim(word))
            case('GEOMECHANICS_REGION','REGION')
              call InputReadWord(input,option, &
                                 geomech_observation%region_name,PETSC_TRUE)
              call InputErrorMsg(input,option,'REGION', &
                                 'GEOMECHANICS_OBSERVATION')
            case default
              call InputKeywordUnrecognized(input,word, &
                                 'GEOMECHANICS_OBSERVATION',option)
          end select
        enddo
        call InputPopBlock(input,option)
        if (len_trim(geomech_observation%region_name) < 1) then
          option%io_buffer = 'A GEOMECHANICS_REGION must be specified ' // &
            'in GEOMECHANICS_OBSERVATION: ' // &
            trim(geomech_observation%name) // '.'
          call PrintErrMsg(option)
        endif
        call GeomechObservationAddToList(geomech_observation, &
          geomech_realization%geomech_patch%geomech_observation_list)
        nullify(geomech_observation)

      !.........................................................................
      case ('GEOMECHANICS_IMPROVE_TET_WEIGHT')
        option%geomechanics%improve_tet_weighting = PETSC_TRUE

      !.........................................................................
      case ('END_GEOMECHANICS')
        exit

      !.........................................................................
      case default
        call InputKeywordUnrecognized(input,word, &
                                 'GeomechanicsInitReadInput',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine InitSubsurfGeomechReadInput

! ************************************************************************** !

subroutine InitSubsurfGeomechJumpStart(geomech)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !
  ! jaa: moved from factory_geomechanics.F90 on 1/28/25

  use Geomechanics_Realization_class
  use Option_module
  use Timestepper_KSP_class
  use Output_Aux_module
  use Output_module, only : Output, OutputPrintCouplers
  use Output_Geomechanics_module
  use Logging_module
  use Condition_Control_module
  use Simulation_Subsurface_class
  use PMC_Geomechanics_class
  use Geomechanics_Attr_module

  implicit none

  type(geomechanics_attr_type) :: geomech

  class(realization_geomech_type), pointer :: geomech_realization
  class(pmc_geomechanics_type), pointer :: geomech_pmc
  class(timestepper_ksp_type), pointer :: geomech_timestepper
  type(option_type), pointer :: option

  PetscBool :: snapshot_plot_flag,observation_plot_flag,massbal_plot_flag
  PetscBool :: geomech_read
  PetscBool :: failure
  PetscErrorCode :: ierr

  geomech_realization => geomech%realization
  geomech_pmc => geomech%process_model_coupler
  geomech_timestepper => TimestepperKSPCast(geomech_pmc%timestepper)
  option => geomech_pmc%option

  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           "-vecload_block_size",failure,ierr);CHKERRQ(ierr)

  geomech_timestepper%name = 'GEOMECHANICS'

  snapshot_plot_flag = PETSC_FALSE
  observation_plot_flag = PETSC_FALSE
  massbal_plot_flag = PETSC_FALSE
  geomech_read = PETSC_FALSE
  failure = PETSC_FALSE

  call OutputGeomechInit(geomech_timestepper%steps)

end subroutine InitSubsurfGeomechJumpStart

! ************************************************************************** !

subroutine InitSubsurfGeomechReadGridBlock(geomech_realization,input,option)
  !
  ! Reads the required geomechanics data from input file
  !
  ! Author: Satish Karra, LANL
  ! Date: 05/23/13
  !
  ! jaa: moved from factory_geomechanics.F90 on 1/28/25

  use Option_module
  use Input_Aux_module
  use String_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  type(input_type), pointer :: input
  type(option_type), pointer :: option

  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_patch_type), pointer :: patch
  character(len=MAXWORDLENGTH) :: word
  type(grid_unstructured_type), pointer :: ugrid
  character(len=MAXWORDLENGTH) :: card
  PetscBool :: read_external_unstructured

  geomech_discretization => geomech_realization%geomech_discretization

  input%ierr = 0
  ! we initialize the word to blanks to avoid error reported by valgrind
  word = ''
  read_external_unstructured = PETSC_FALSE

  geomech_discretization%grid  => GMGridCreate()
  ugrid => UGridCreate()

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,card)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS')
    call StringToUpper(word)

    select case(trim(word))
      case ('TYPE')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword','TYPE')
        call StringToUpper(word)

        select case(trim(word))
          case ('UNSTRUCTURED')
            geomech_discretization%itype = UNSTRUCTURED_GRID
            geomech_discretization%ctype = 'UNSTRUCTURED'
            call InputReadFilename(input,option,geomech_discretization%filename)
            call InputErrorMsg(input,option,'keyword','filename')
            read_external_unstructured = PETSC_TRUE
          case ('STRUCTURED_INTERNAL','STRUCTURED')
            ! Geomechanics currently uses the unstructured DM path.  For
            ! structured flow grids we build an internal hex mesh at setup.
            geomech_discretization%itype = UNSTRUCTURED_GRID
            geomech_discretization%ctype = 'STRUCTURED_INTERNAL'
          case default
            option%io_buffer = 'GEOMECHANICS_GRID TYPE must be UNSTRUCTURED or STRUCTURED_INTERNAL.'
            call PrintErrMsg(option)
        end select
      case ('GRAVITY')
        call InputReadDouble(input,option,option%geomechanics% &
                             gravity(X_DIRECTION))
        call InputErrorMsg(input,option,'x-direction','GEOMECH GRAVITY')
        call InputReadDouble(input,option,option%geomechanics% &
                             gravity(Y_DIRECTION))
        call InputErrorMsg(input,option,'y-direction','GEOMECH GRAVITY')
        call InputReadDouble(input,option,option%geomechanics% &
                             gravity(Z_DIRECTION))
        call InputErrorMsg(input,option,'z-direction','GEOMECH GRAVITY')
        if (OptionIsIORank(option) .and. OptionPrintToScreen(option)) &
            write(option%fid_out,'(/," *GEOMECH_GRAV",/, &
            & "  gravity    = "," [m/s^2]",3x,1p3e12.4 &
            & )') option%geomechanics%gravity(1:3)
      case ('MAX_CELLS_SHARING_A_VERTEX')
        call InputReadInt(input,option,ugrid%max_cells_sharing_a_vertex)
        call InputErrorMsg(input,option,'max_cells_sharing_a_vertex', &
                           'GEOMECHANICS_GRID')
      case default
        call InputKeywordUnrecognized(input,word,'GEOMECHANICS_GRID',option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (read_external_unstructured) then
    call UGridRead(ugrid,geomech_discretization%filename,option)
    call UGridDecompose(ugrid,option)
    call CopySubsurfaceGridtoGeomechGrid(ugrid, &
                                         geomech_discretization%grid, &
                                         option)
  endif
  patch => GeomechanicsPatchCreate()
  patch%geomech_grid => geomech_discretization%grid
  geomech_realization%geomech_patch => patch

end subroutine InitSubsurfGeomechReadGridBlock

! ************************************************************************** !

subroutine InitSubsurfGeomechBuildStructuredInternalGrid(subsurf_realization, &
                                                         geomech_realization)
  !
  ! Builds an internal HEX geomechanics mesh from the structured flow grid
  ! and initializes flow-geomech coupling maps.
  !
  ! Author: Satish Karra
  ! Date: 02/23/26
  !

#include "petsc/finclude/petscao.h"
  use petscao
#include "petsc/finclude/petscis.h"
  use petscis

  use Option_module
  use Realization_Subsurface_class
  use Geomechanics_Realization_class
  use Geomechanics_Discretization_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Option_Geomechanics_module
  use Grid_module
  use Grid_Structured_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Cell_module

  implicit none

  class(realization_subsurface_type), pointer :: subsurf_realization
  class(realization_geomech_type), pointer :: geomech_realization

  type(option_type), pointer :: option
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_grid_type), pointer :: geomech_grid
  type(grid_type), pointer :: flow_grid
  type(grid_structured_type), pointer :: sgrid
  type(grid_unstructured_type), pointer :: ugrid
  PetscInt :: nx, ny, nz, nxy
  PetscInt :: nxv, nyv, nzv, nxyv
  PetscInt :: local_id, ghosted_id
  PetscInt :: ic, jc, kc
  PetscInt :: iv(8), vid
  PetscInt :: count, vertex_count
  PetscInt, allocatable :: all_vertex_ids(:)
  PetscInt, allocatable :: sorted_perm(:)
  PetscInt, allocatable :: unique_vertex_ids(:)
  PetscInt, allocatable :: occ_to_unique(:)
  PetscInt, allocatable :: ao_natural(:)
  PetscInt, allocatable :: ao_petsc(:)
  IS :: is_ao_natural
  IS :: is_ao_petsc
  PetscReal, allocatable :: xv(:), yv(:), zv(:)
  PetscInt :: i, j, k
  PetscInt :: ix, iy, iz
  PetscInt :: owner_i, owner_j, owner_k
  PetscInt :: i_start, i_end, j_start, j_end, k_start, k_end
  PetscInt :: ci, cj, ck
  PetscInt :: map_count
  PetscInt, allocatable :: tmp_map_cells(:)
  PetscInt, allocatable :: tmp_map_vertices(:)
  PetscErrorCode :: ierr

  option => geomech_realization%option
  geomech_discretization => geomech_realization%geomech_discretization

  if (trim(geomech_discretization%ctype) /= 'STRUCTURED_INTERNAL') return

  flow_grid => subsurf_realization%discretization%grid
  if (.not.associated(flow_grid)) then
    option%io_buffer = 'Structured-internal geomechanics mesh requested, but subsurface grid is not available.'
    call PrintErrMsg(option)
  endif
  if (flow_grid%itype /= STRUCTURED_GRID) then
    option%io_buffer = 'GEOMECHANICS_GRID TYPE STRUCTURED_INTERNAL requires a structured flow GRID.'
    call PrintErrMsg(option)
  endif

  sgrid => flow_grid%structured_grid
  if (.not.associated(sgrid)) then
    option%io_buffer = 'Structured flow grid metadata is missing.'
    call PrintErrMsg(option)
  endif
  if (.not.associated(sgrid%dx_global) .or. .not.associated(sgrid%dy_global) .or. &
      .not.associated(sgrid%dz_global)) then
    option%io_buffer = 'Structured flow grid spacing has not been initialized before geomechanics setup.'
    call PrintErrMsg(option)
  endif

  nx = sgrid%nx
  ny = sgrid%ny
  nz = sgrid%nz
  nxy = nx*ny
  nxv = nx + 1
  nyv = ny + 1
  nzv = nz + 1
  nxyv = nxv*nyv

  allocate(xv(0:nx))
  allocate(yv(0:ny))
  allocate(zv(0:nz))
  xv(0) = sgrid%bounds(X_DIRECTION,LOWER)
  yv(0) = sgrid%bounds(Y_DIRECTION,LOWER)
  zv(0) = sgrid%bounds(Z_DIRECTION,LOWER)
  do i = 1, nx
    xv(i) = xv(i-1) + sgrid%dx_global(i)
  enddo
  do j = 1, ny
    yv(j) = yv(j-1) + sgrid%dy_global(j)
  enddo
  do k = 1, nz
    zv(k) = zv(k-1) + sgrid%dz_global(k)
  enddo

  ugrid => UGridCreate()
  ugrid%nmax = flow_grid%nmax
  ugrid%nlmax = flow_grid%nlmax
  ugrid%ngmax = flow_grid%ngmax
  ugrid%global_offset = flow_grid%global_offset
  ugrid%num_vertices_global = nxv*nyv*nzv
  ugrid%max_ndual_per_cell = 8
  ugrid%max_nvert_per_cell = 8
  ugrid%max_cells_sharing_a_vertex = max(8,ugrid%max_cells_sharing_a_vertex)

  allocate(ugrid%cell_ids_natural(ugrid%nlmax))
  allocate(ugrid%cell_ids_petsc(ugrid%nlmax))
  allocate(ugrid%cell_type(ugrid%nlmax))
  allocate(ugrid%cell_vertices(0:ugrid%max_nvert_per_cell,ugrid%nlmax))
  ugrid%cell_type = HEX_TYPE
  ugrid%cell_vertices = 0
  ugrid%cell_vertices(0,:) = 8

  allocate(all_vertex_ids(ugrid%nlmax*8))
  count = 0
  do local_id = 1, ugrid%nlmax
    ghosted_id = flow_grid%nL2G(local_id)
    ugrid%cell_ids_natural(local_id) = flow_grid%nG2A(ghosted_id)
    ugrid%cell_ids_petsc(local_id) = flow_grid%global_offset + local_id

    ic = mod(ugrid%cell_ids_natural(local_id)-1,nx)
    jc = mod((ugrid%cell_ids_natural(local_id)-1)/nx,ny)
    kc = (ugrid%cell_ids_natural(local_id)-1)/nxy

    iv(1) = 1 + (ic)   + (jc)*nxv + (kc)*nxyv
    iv(2) = 1 + (ic+1) + (jc)*nxv + (kc)*nxyv
    iv(3) = 1 + (ic+1) + (jc+1)*nxv + (kc)*nxyv
    iv(4) = 1 + (ic)   + (jc+1)*nxv + (kc)*nxyv
    iv(5) = 1 + (ic)   + (jc)*nxv + (kc+1)*nxyv
    iv(6) = 1 + (ic+1) + (jc)*nxv + (kc+1)*nxyv
    iv(7) = 1 + (ic+1) + (jc+1)*nxv + (kc+1)*nxyv
    iv(8) = 1 + (ic)   + (jc+1)*nxv + (kc+1)*nxyv

    do i = 1, 8
      count = count + 1
      all_vertex_ids(count) = iv(i)
    enddo
  enddo

  allocate(sorted_perm(count))
  do i = 1, count
    sorted_perm(i) = i
  enddo
  sorted_perm = sorted_perm - 1
  call PetscSortIntWithPermutation(count,all_vertex_ids,sorted_perm, &
                                   ierr);CHKERRQ(ierr)
  sorted_perm = sorted_perm + 1

  allocate(unique_vertex_ids(count))
  allocate(occ_to_unique(count))
  unique_vertex_ids = 0
  occ_to_unique = 0
  unique_vertex_ids(1) = all_vertex_ids(sorted_perm(1))
  vertex_count = 1
  occ_to_unique(sorted_perm(1)) = vertex_count
  do i = 2, count
    vid = all_vertex_ids(sorted_perm(i))
    if (vid > unique_vertex_ids(vertex_count)) then
      vertex_count = vertex_count + 1
      unique_vertex_ids(vertex_count) = vid
    endif
    occ_to_unique(sorted_perm(i)) = vertex_count
  enddo

  count = 0
  do local_id = 1, ugrid%nlmax
    do i = 1, 8
      count = count + 1
      ugrid%cell_vertices(i,local_id) = occ_to_unique(count)
    enddo
  enddo

  ugrid%num_vertices_local = vertex_count
  allocate(ugrid%vertex_ids_natural(vertex_count))
  allocate(ugrid%vertices(vertex_count))
  do local_id = 1, vertex_count
    ugrid%vertex_ids_natural(local_id) = unique_vertex_ids(local_id)
    vid = unique_vertex_ids(local_id) - 1
    iz = vid/nxyv
    iy = mod(vid,nxyv)/nxv
    ix = mod(vid,nxv)
    ugrid%vertices(local_id)%id = unique_vertex_ids(local_id)
    ugrid%vertices(local_id)%x = xv(ix)
    ugrid%vertices(local_id)%y = yv(iy)
    ugrid%vertices(local_id)%z = zv(iz)
  enddo

  allocate(ao_natural(ugrid%nlmax))
  allocate(ao_petsc(ugrid%nlmax))
  ao_natural = ugrid%cell_ids_natural - 1
  ao_petsc = ugrid%cell_ids_petsc - 1
  call ISCreateGeneral(option%mycomm,ugrid%nlmax,ao_natural, &
                       PETSC_COPY_VALUES,is_ao_natural,ierr);CHKERRQ(ierr)
  call ISCreateGeneral(option%mycomm,ugrid%nlmax,ao_petsc, &
                       PETSC_COPY_VALUES,is_ao_petsc,ierr);CHKERRQ(ierr)
  call AOCreateMappingIS(is_ao_natural,is_ao_petsc, &
                         ugrid%ao_natural_to_petsc,ierr);CHKERRQ(ierr)
  call ISDestroy(is_ao_natural,ierr);CHKERRQ(ierr)
  call ISDestroy(is_ao_petsc,ierr);CHKERRQ(ierr)

  call CopySubsurfaceGridtoGeomechGrid(ugrid,geomech_discretization%grid,option)
  geomech_grid => geomech_discretization%grid

  ! Internal map for structured HEX coupling.
  ! 0th order: one owner flow cell per geomech vertex (latest-cell-wins).
  ! 1st order: one geomech vertex mapped to all adjacent flow cells.
  if (associated(geomech_grid%mapping_cell_ids_flow)) then
    deallocate(geomech_grid%mapping_cell_ids_flow)
  endif
  if (associated(geomech_grid%mapping_vertex_ids_geomech)) then
    deallocate(geomech_grid%mapping_vertex_ids_geomech)
  endif
  select case(option%geomechanics%flow_interp_order)
    case(GEOMECH_FLOW_INTERP_ORDER_0TH)
      geomech_grid%mapping_num_cells = geomech_grid%nlmax_node
      allocate(geomech_grid%mapping_cell_ids_flow(geomech_grid%mapping_num_cells))
      allocate(geomech_grid%mapping_vertex_ids_geomech(geomech_grid%mapping_num_cells))

      do local_id = 1, geomech_grid%nlmax_node
        vid = geomech_grid%node_ids_local_natural(local_id)
        geomech_grid%mapping_vertex_ids_geomech(local_id) = vid

        ! Decode 1-based structured vertex id -> (ix,iy,iz), zero-based.
        count = vid - 1
        iz = count/nxyv
        iy = mod(count,nxyv)/nxv
        ix = mod(count,nxv)

        ! Latest-cell-wins owner for this vertex in global natural ordering.
        owner_i = min(ix + 1, nx)
        owner_j = min(iy + 1, ny)
        owner_k = min(iz + 1, nz)
        geomech_grid%mapping_cell_ids_flow(local_id) = &
          owner_i + (owner_j-1)*nx + (owner_k-1)*nxy
      enddo

    case(GEOMECH_FLOW_INTERP_ORDER_1ST)
      ! Upper bound: 8 adjacent cells per local geomech node.
      allocate(geomech_grid%mapping_cell_ids_flow(8*geomech_grid%nlmax_node))
      allocate(geomech_grid%mapping_vertex_ids_geomech(8*geomech_grid%nlmax_node))

      map_count = 0
      do local_id = 1, geomech_grid%nlmax_node
        vid = geomech_grid%node_ids_local_natural(local_id)

        count = vid - 1
        iz = count/nxyv
        iy = mod(count,nxyv)/nxv
        ix = mod(count,nxv)

        i_start = max(1,ix)
        i_end   = min(nx,ix+1)
        j_start = max(1,iy)
        j_end   = min(ny,iy+1)
        k_start = max(1,iz)
        k_end   = min(nz,iz+1)

        do ck = k_start, k_end
          do cj = j_start, j_end
            do ci = i_start, i_end
              map_count = map_count + 1
              geomech_grid%mapping_vertex_ids_geomech(map_count) = vid
              geomech_grid%mapping_cell_ids_flow(map_count) = &
                ci + (cj-1)*nx + (ck-1)*nxy
            enddo
          enddo
        enddo
      enddo

      geomech_grid%mapping_num_cells = map_count
      if (map_count < size(geomech_grid%mapping_cell_ids_flow)) then
        allocate(tmp_map_cells(map_count))
        allocate(tmp_map_vertices(map_count))
        tmp_map_cells = geomech_grid%mapping_cell_ids_flow(1:map_count)
        tmp_map_vertices = geomech_grid%mapping_vertex_ids_geomech(1:map_count)
        deallocate(geomech_grid%mapping_cell_ids_flow)
        deallocate(geomech_grid%mapping_vertex_ids_geomech)
        allocate(geomech_grid%mapping_cell_ids_flow(map_count))
        allocate(geomech_grid%mapping_vertex_ids_geomech(map_count))
        geomech_grid%mapping_cell_ids_flow = tmp_map_cells
        geomech_grid%mapping_vertex_ids_geomech = tmp_map_vertices
        deallocate(tmp_map_cells)
        deallocate(tmp_map_vertices)
      endif

    case default
      option%io_buffer = 'Invalid GEOMECHANICS_FLOW_INTERPOLATION ORDER for structured-internal coupling. Use 0 or 1.'
      call PrintErrMsg(option)
  end select

  call UGridDestroy(ugrid)

  deallocate(all_vertex_ids)
  deallocate(sorted_perm)
  deallocate(unique_vertex_ids)
  deallocate(occ_to_unique)
  deallocate(ao_natural)
  deallocate(ao_petsc)
  deallocate(xv)
  deallocate(yv)
  deallocate(zv)

end subroutine InitSubsurfGeomechBuildStructuredInternalGrid

! ************************************************************************** !

subroutine InitSubsurfGeomechSetupRealization(subsurf_realization, &
                                              geomech_realization)
  !
  ! Initializes material property data structres and assign them to the domain.
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  !
  ! jaa: moved from factory_geomechanics.F90 on 1/28/25

  use Geomechanics_Realization_class
  use Geomechanics_Global_module
  use Geomechanics_Force_module
  use Realization_Subsurface_class
  use Simulation_Subsurface_class

  use Option_module

  implicit none

  class(realization_subsurface_type), pointer :: subsurf_realization
  class(realization_geomech_type), pointer :: geomech_realization

  type(option_type), pointer :: option

  option => subsurf_realization%option

  call InitSubsurfGeomechBuildStructuredInternalGrid(subsurf_realization, &
                                                     geomech_realization)

  call GeomechRealizCreateDiscretization(geomech_realization)

  if (option%geomechanics%flow_coupling /= 0 .or. &
      option%geomechanics%geophysics_coupling /= 0) then
    call GeomechCreateGeomechSubsurfVec(subsurf_realization, &
                                        geomech_realization)
    call GeomechCreateSubsurfStressStrainVec(subsurf_realization, &
                                              geomech_realization)

    call GeomechRealizMapSubsurfGeomechGrid(subsurf_realization, &
                                            geomech_realization, &
                                            option)
  endif
  call GeomechRealizLocalizeRegions(geomech_realization)
  call GeomechRealizPassFieldPtrToPatch(geomech_realization)
  call GeomechRealizProcessMatProp(geomech_realization)
  call GeomechRealizProcessGeomechCouplers(geomech_realization)
  call GeomechRealizProcessGeomechConditions(geomech_realization)
  call InitMatPropToGeomechRegions(geomech_realization)
  call GeomechRealizInitAllCouplerAuxVars(geomech_realization)
  call GeomechRealizPrintCouplers(geomech_realization)
  call GeomechGridElemSharedByNodes(geomech_realization,option)
  call GeomechForceSetup(geomech_realization)
  call GeomechGlobalSetup(geomech_realization)

  ! SK: We are solving quasi-steady state solution for geomechanics.
  ! Initial condition is not needed, hence CondControlAssignFlowInitCondGeomech
  ! is not needed, at this point.
  call GeomechForceUpdateAuxVars(geomech_realization)

end subroutine InitSubsurfGeomechSetupRealization

! ************************************************************************** !

subroutine InitMatPropToGeomechRegions(geomech_realization)
  !
  ! This routine assigns geomech material
  ! properties to associated regions
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  !
  ! jaa: moved from factory_geomechanics.F90 on 1/28/25

  use Geomechanics_Realization_class
  use Geomechanics_Discretization_module
  use Geomechanics_Strata_module
  use Geomechanics_Region_module
  use Geomechanics_Material_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Field_module
  use Geomechanics_Patch_module
  use Realization_Subsurface_class, only : MATERIAL_ID_ARRAY
  use Option_module

  implicit none

  class(realization_geomech_type) :: geomech_realization

  PetscInt :: ivertex, local_id, ghosted_id, geomech_material_id
  PetscInt :: istart, iend
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscErrorCode :: ierr

  Vec :: temp_vec_loc
  PetscReal, pointer :: temp_vec_loc_p(:)

  type(option_type), pointer :: option
  type(geomech_grid_type), pointer :: grid
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_field_type), pointer :: field
  type(geomech_strata_type), pointer :: strata
  type(geomech_patch_type), pointer :: patch

  type(geomech_material_property_type), pointer :: geomech_material_property
  type(geomech_material_property_type), pointer :: null_geomech_material_property
  type(gm_region_type), pointer :: region
  PetscBool :: update_ghosted_material_ids
  PetscReal, pointer :: imech_loc_p(:)

  option => geomech_realization%option
  geomech_discretization => geomech_realization%geomech_discretization
  field => geomech_realization%geomech_field
  patch => geomech_realization%geomech_patch

  ! loop over all patches and allocation material id arrays
  if (.not.associated(patch%imat)) then
    allocate(patch%imat(patch%geomech_grid%ngmax_node))
    ! initialize to "unset"
    patch%imat = UNINITIALIZED_INTEGER
  endif

  ! if material ids are set based on region, as opposed to being read in
  ! we must communicate the ghosted ids.  This flag toggles this operation.
  update_ghosted_material_ids = PETSC_FALSE
  grid => patch%geomech_grid
  strata => patch%geomech_strata_list%first
  do
    if (.not.associated(strata)) exit
    ! Read in cell by cell material ids if they exist
    if (.not.associated(strata%region) .and. strata%active) then
      option%io_buffer = 'Reading of material prop from file for' // &
        ' geomech is not implemented.'
      call PrintErrMsgByRank(option)
    ! Otherwise, set based on region
    else if (strata%active) then
      update_ghosted_material_ids = PETSC_TRUE
      region => strata%region
      geomech_material_property => strata%material_property
      if (associated(region)) then
        istart = 1
        iend = region%num_verts
      else
        istart = 1
        iend = grid%nlmax_node
      endif
      do ivertex = istart, iend
        if (associated(region)) then
          local_id = region%vertex_ids(ivertex)
        else
          local_id = ivertex
        endif
        ghosted_id = grid%nL2G(local_id)
        patch%imat(ghosted_id) = geomech_material_property%id
      enddo
    endif
    strata => strata%next
  enddo

  if (update_ghosted_material_ids) then
    ! update ghosted material ids
    call GeomechRealizLocalToLocalWithArray(geomech_realization, &
                                            MATERIAL_ID_ARRAY)
  endif

  ! set cell by cell material properties
  ! create null material property for inactive cells
  null_geomech_material_property => GeomechanicsMaterialPropertyCreate()
  call VecGetArray(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax_node
    ghosted_id = grid%nL2G(local_id)
    geomech_material_id = patch%imat(ghosted_id)
    if (geomech_material_id == 0) then ! accomodate inactive cells
      geomech_material_property = null_geomech_material_property
    else if ( geomech_material_id > 0 .and. &
              geomech_material_id <= &
              size(geomech_realization%geomech_material_property_array)) then
      geomech_material_property => &
         geomech_realization% &
           geomech_material_property_array(geomech_material_id)%ptr
      if (.not.associated(geomech_material_property)) then
        write(dataset_name,*) geomech_material_id
        option%io_buffer = 'No material property for geomech material id ' // &
                            trim(adjustl(dataset_name)) &
                            //  ' defined in input file.'
        call PrintErrMsgByRank(option)
      endif
    else if (Uninitialized(geomech_material_id)) then
      write(dataset_name,*) grid%nG2A(ghosted_id)
      option%io_buffer = 'Uninitialized geomech material id in patch at cell ' // &
                         trim(adjustl(dataset_name))
      call PrintErrMsgByRank(option)
    else if (geomech_material_id > size(geomech_realization% &
      geomech_material_property_array)) then
      write(option%io_buffer,*) geomech_material_id
      option%io_buffer = 'Unmatched geomech material id in patch:' // &
        adjustl(trim(option%io_buffer))
      call PrintErrMsgByRank(option)
    else
      option%io_buffer = 'Something messed up with geomech material ids. ' // &
        ' Possibly material ids not assigned to all grid cells. ' // &
        ' Contact Glenn/Satish!'
      call PrintErrMsgByRank(option)
    endif
    imech_loc_p(ghosted_id) = geomech_material_property%id
  enddo ! local_id - loop
  call VecRestoreArray(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)

  ! read in any user-defined geomech property fields
  call GeomechDiscretizationDuplicateVector(geomech_discretization, &
                                            field%press_loc, &
                                            temp_vec_loc)
  do geomech_material_id = 1, size(patch%geomech_material_property_array)
    geomech_material_property => &
            patch%geomech_material_property_array(geomech_material_id)%ptr
    if (.not.associated(geomech_material_property)) cycle
    ! Young's modulus
    if (associated(geomech_material_property%youngs_modulus_dataset)) then
      ! Set the value of field%youngs_modulus for this material to the dataset values
      call GeomechReadDatasetToVecWithMask(geomech_realization, &
            geomech_material_property%youngs_modulus_dataset, &
            geomech_material_property%id,PETSC_FALSE,field%youngs_modulus_loc,temp_vec_loc)
    else
      ! Set the value of field%youngs_modulus for this material to the constant value
      call VecGetArray(field%youngs_modulus_loc,temp_vec_loc_p,ierr);CHKERRQ(ierr)
      do local_id = 1, grid%nlmax_node
        if (patch%imat(grid%nL2G(local_id)) == geomech_material_property%id) then
          temp_vec_loc_p(local_id) = geomech_material_property%youngs_modulus
        endif
      enddo
      call VecRestoreArray(field%youngs_modulus_loc,temp_vec_loc_p,ierr);CHKERRQ(ierr)
    endif
    ! Poisson's ratio
    if (associated(geomech_material_property%poissons_ratio_dataset)) then
      ! Set the value of field%poissons_ratio for this material to the dataset values
      call GeomechReadDatasetToVecWithMask(geomech_realization, &
            geomech_material_property%poissons_ratio_dataset, &
            geomech_material_property%id,PETSC_FALSE,field%poissons_ratio_loc,temp_vec_loc)
    else
      ! Set the value of field%poissons_ratio for this material to the constant value
      call VecGetArray(field%poissons_ratio_loc,temp_vec_loc_p,ierr);CHKERRQ(ierr)
      do local_id = 1, grid%nlmax_node
        if (patch%imat(grid%nL2G(local_id)) == geomech_material_property%id) then
          temp_vec_loc_p(local_id) = geomech_material_property%poissons_ratio
        endif
      enddo
      call VecRestoreArray(field%poissons_ratio_loc,temp_vec_loc_p,ierr);CHKERRQ(ierr)
    endif
    ! Density
    if (associated(geomech_material_property%density_dataset)) then
      ! Set the value of field%density for this material to the dataset values
      call GeomechReadDatasetToVecWithMask(geomech_realization, &
            geomech_material_property%density_dataset, &
            geomech_material_property%id,PETSC_FALSE,field%density_loc,temp_vec_loc)
    else
      ! Set the value of field%density for this material to the constant value
      call VecGetArray(field%density_loc,temp_vec_loc_p,ierr);CHKERRQ(ierr)
      do local_id = 1, grid%nlmax_node
        if (patch%imat(grid%nL2G(local_id)) == geomech_material_property%id) then
          temp_vec_loc_p(local_id) = geomech_material_property%density
        endif
      enddo
      call VecRestoreArray(field%density_loc,temp_vec_loc_p,ierr);CHKERRQ(ierr)
    endif
    ! Biot's coefficient
    if (associated(geomech_material_property%biot_coeff_dataset)) then
      ! Set the value of field%biot_coeff for this material to the dataset values
      call GeomechReadDatasetToVecWithMask(geomech_realization, &
            geomech_material_property%biot_coeff_dataset, &
            geomech_material_property%id,PETSC_FALSE,field%biot_coeff_loc,temp_vec_loc)
    else
      ! Set the value of field%biot_coeff for this material to the constant value
      call VecGetArray(field%biot_coeff_loc,temp_vec_loc_p,ierr);CHKERRQ(ierr)
      do local_id = 1, grid%nlmax_node
        if (patch%imat(grid%nL2G(local_id)) == geomech_material_property%id) then
          temp_vec_loc_p(local_id) = geomech_material_property%biot_coeff
        endif
      enddo
      call VecRestoreArray(field%biot_coeff_loc,temp_vec_loc_p,ierr);CHKERRQ(ierr)
    endif
    ! Thermal expansion coefficient
    if (associated(geomech_material_property%thermal_exp_coeff_dataset)) then
      ! Set the value of field%thermal_exp_coeff for this material to the dataset values
      call GeomechReadDatasetToVecWithMask(geomech_realization, &
            geomech_material_property%thermal_exp_coeff_dataset, &
            geomech_material_property%id,PETSC_FALSE,field%thermal_exp_coeff_loc,temp_vec_loc)
    else
      ! Set the value of field%thermal_exp_coeff for this material to the constant value
      call VecGetArray(field%thermal_exp_coeff_loc,temp_vec_loc_p,ierr);CHKERRQ(ierr)
      do local_id = 1, grid%nlmax_node
        if (patch%imat(grid%nL2G(local_id)) == geomech_material_property%id) then
          temp_vec_loc_p(local_id) = geomech_material_property%thermal_exp_coeff
        endif
      enddo
      call VecRestoreArray(field%thermal_exp_coeff_loc,temp_vec_loc_p,ierr);CHKERRQ(ierr)
    endif
  enddo

  ! Scatter owned-node property values to ghost nodes
  call GeomechDiscretizationLocalToLocal(geomech_discretization, &
                                         field%youngs_modulus_loc, &
                                         field%youngs_modulus_loc,ONEDOF)
  call GeomechDiscretizationLocalToLocal(geomech_discretization, &
                                         field%poissons_ratio_loc, &
                                         field%poissons_ratio_loc,ONEDOF)
  call GeomechDiscretizationLocalToLocal(geomech_discretization, &
                                         field%density_loc, &
                                         field%density_loc,ONEDOF)
  call GeomechDiscretizationLocalToLocal(geomech_discretization, &
                                         field%biot_coeff_loc, &
                                         field%biot_coeff_loc,ONEDOF)
  call GeomechDiscretizationLocalToLocal(geomech_discretization, &
                                         field%thermal_exp_coeff_loc, &
                                         field%thermal_exp_coeff_loc,ONEDOF)

  call VecSet(field%body_force_x_loc,geomech_realization%body_force(X_DIRECTION),ierr);CHKERRQ(ierr)
  call VecSet(field%body_force_y_loc,geomech_realization%body_force(Y_DIRECTION),ierr);CHKERRQ(ierr)
  call VecSet(field%body_force_z_loc,geomech_realization%body_force(Z_DIRECTION),ierr);CHKERRQ(ierr)

  if (associated(geomech_realization%body_force_x_dataset)) then
    call GeomechReadDatasetToVecWithMask(geomech_realization, &
          geomech_realization%body_force_x_dataset, &
          0,PETSC_TRUE,field%body_force_x_loc,temp_vec_loc)
  endif
  if (associated(geomech_realization%body_force_y_dataset)) then
    call GeomechReadDatasetToVecWithMask(geomech_realization, &
          geomech_realization%body_force_y_dataset, &
          0,PETSC_TRUE,field%body_force_y_loc,temp_vec_loc)
  endif
  if (associated(geomech_realization%body_force_z_dataset)) then
    call GeomechReadDatasetToVecWithMask(geomech_realization, &
          geomech_realization%body_force_z_dataset, &
          0,PETSC_TRUE,field%body_force_z_loc,temp_vec_loc)
  endif

  call VecDestroy(temp_vec_loc,ierr);CHKERRQ(ierr)

  call GeomechanicsMaterialPropertyDestroy(null_geomech_material_property)
  nullify(null_geomech_material_property)

  call GeomechDiscretizationLocalToLocal(geomech_discretization,field%imech_loc, &
                                         field%imech_loc,ONEDOF)

end subroutine InitMatPropToGeomechRegions

! ************************************************************************** !

subroutine InitSubsurfGeomechInitSimulation(simulation, pm_geomech)
  !
  ! This routine initializes geomechanics process
  ! model components for factory subsurface
  ! after linkages have been established
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 2/10/25
  !
  use Simulation_Subsurface_class
  use Init_Common_module
  use Option_module
  use PM_Base_class
  use PM_Base_Pointer_module
  use PM_Geomechanics_Force_class
  use PMC_Base_class
  use PMC_Geomechanics_class
  use PFLOTRAN_Constants_module
  use Geomechanics_Discretization_module
  use Geomechanics_Grid_Aux_module
  use GMDM_Pointer_module
  use Geomechanics_Force_module
  use Geomechanics_Realization_class
  use Geomechanics_Regression_module
  use Simulation_Aux_module
  use Realization_Subsurface_class
  use Realization_Base_class
  use Timestepper_KSP_class
  use Logging_module
  use Output_Aux_module
  use Waypoint_module
  use PM_Richards_class
  use PM_TH_class

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_geomech_force_type), pointer :: pm_geomech

  type(option_type), pointer :: option
  class(realization_subsurface_type), pointer :: subsurf_realization
  class(realization_geomech_type), pointer :: geomech_realization
  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: pmc_dummy
  type(gmdm_ptr_type), pointer :: dm_ptr
  class(pmc_geomechanics_type), pointer :: pmc_geomech
  class(timestepper_ksp_type), pointer :: timestepper
  type(geomechanics_regression_type), pointer :: geomech_regression
  PetscErrorCode :: ierr

  if (.not. associated(pm_geomech)) return

  nullify(pmc_dummy)

  option => simulation%option
  geomech_realization => simulation%geomech%realization
  subsurf_realization => simulation%realization
  pmc_geomech => simulation%geomech%process_model_coupler
  timestepper => TimestepperKSPCast(pmc_geomech%timestepper)

  geomech_regression => simulation%geomech%regression

  if (option%flow%creep_closure_on) then
    option%io_buffer = 'Creep closure is not supported with geomechanics&
      &-subsurface coupling. Disable creep closure or geomechanics &
      &coupling in InitSubsurfGeomechInitSimulation.'
    call PrintErrMsg(option)
  endif

  ! initialize geomech realization
  call InitSubsurfGeomechSetupRealization(simulation%realization,&
                                          simulation%geomech%realization)

  call pm_geomech%PMGeomechForceSetRealization(geomech_realization, &
                                               subsurf_realization)
  call pm_geomech%Setup()

  call pmc_geomech%SetupSolvers()

  ! Here I first calculate the linear part of the jacobian and store it
  ! since the jacobian is always linear with geomech (even when coupled with
  ! flow since we are performing sequential coupling). Although
  ! SNESSetJacobian is called, nothing is done there and PETSc just
  ! re-uses the linear Jacobian at all iterations and times
  call MatSetOption(timestepper%solver%M,MAT_NEW_NONZERO_ALLOCATION_ERR, &
                    PETSC_FALSE,ierr);CHKERRQ(ierr)
  call GeomechForceAssembleCoeffMatrix(timestepper%solver%M, &
                                      geomech_realization)
  call MatSetOption(timestepper%solver%M,MAT_NEW_NONZERO_ALLOCATION_ERR, &
                    PETSC_TRUE,ierr);CHKERRQ(ierr)
  nullify(simulation%process_model_coupler_list)

  ! sim_aux: Create PETSc Vectors and VectorScatters
  call SimAuxCopySubsurfVec(simulation%sim_aux,subsurf_realization%field%work)

  call SimAuxCopySubsurfGeomechVec(simulation%sim_aux, &
        geomech_realization%geomech_field%strain_subsurf)

  call GeomechRealizMapSubsurfGeomechGrid(subsurf_realization, &
                                          geomech_realization, &
                                          option)

  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex( &
              geomech_realization%geomech_discretization, ONEDOF)

  call SimAuxCopyVecScatter(simulation%sim_aux, &
                            dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                            SUBSURF_TO_GEOMECHANICS)
  call SimAuxCopyVecScatter(simulation%sim_aux, &
                            dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                            GEOMECHANICS_TO_SUBSURF)

  call GeomechanicsRegressionCreateMapping(geomech_regression, &
                                           geomech_realization)

  ! sim_aux: Set pointer
  simulation%flow_process_model_coupler%sim_aux => simulation%sim_aux
  if (associated(simulation%tran_process_model_coupler)) &
    simulation%tran_process_model_coupler%sim_aux => simulation%sim_aux
  if (option%ngeomechdof>0 .and. associated(pmc_geomech)) &
    pmc_geomech%sim_aux => simulation%sim_aux

  ! set geomech as not master
  pmc_geomech%is_master = PETSC_FALSE
  ! link geomech and master
  ! jaa: set flow as the master
  simulation%process_model_coupler_list => &
    simulation%flow_process_model_coupler
  ! link subsurface flow as peer
  ! jaa: set geomech as a child
!  simulation%process_model_coupler_list%child => &
!    pmc_geomech
  call PMCBaseSetChildPeerPtr(pmc_geomech%CastToBase(),PM_CHILD, &
                    simulation%flow_process_model_coupler%CastToBase(), &
                    pmc_dummy,PM_APPEND)

  call InitSubsurfGeomechChkInactiveCells(geomech_realization, &
                                          subsurf_realization)

  ! Set data in sim_aux
  cur_process_model_coupler => simulation%process_model_coupler_list
  call simulation%flow_process_model_coupler%SetAuxData()
  call simulation%geomech%process_model_coupler%GetAuxData()
  call simulation%geomech%process_model_coupler%SetAuxData()

  ! Set the flow pointer to the geomech_parameter
  select type(pm => simulation%process_model_list)
    class is(pm_th_type)
      pm_geomech%subsurf_realization%patch%aux%th%th_parameter% &
        geomech_parameter => pm_geomech%geomech_realization% &
        geomech_patch%geomech_aux%Linear%linear_parameter

    class is(pm_richards_type)
      pm_geomech%subsurf_realization%patch%aux%richards% &
        richards_parameter%geomech_parameter => &
        pm_geomech%geomech_realization%geomech_patch% &
        geomech_aux%Linear%linear_parameter

    class default
      option%io_buffer = 'Only RICHARDS and TH modes can be coupled'// &
                      'to geomechanics InitSubsurfGeomechInitSimulation'
      call PrintErrMsg(option)
  end select

  ! Map geomechanics material ids onto flow cells for fixed-stress approach
  call GeomechMapMaterialIdsToFlow(subsurf_realization,geomech_realization)

  ! update auxvars once pointer is set
  call simulation%process_model_list%UpdateAuxVars()

  ! this is solely for casting to pmc geomech
  select type(pmc => simulation%geomech%process_model_coupler)
    class is(pmc_geomechanics_type)
      call GeomechStoreInitialPressTemp(pmc%geomech_realization)
  end select

  call InitSubsurfGeomechJumpStart(simulation%geomech)

end subroutine InitSubsurfGeomechInitSimulation

! ************************************************************************** !

subroutine InitSubsurfGeomechSetGeomechMode(pm_geomech,option)
  !
  ! sets geomech options
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 2/10/25
  !
  use Option_module
  use PM_Geomechanics_Force_class

  implicit none

  type(option_type) :: option
  class(pm_geomech_force_type), pointer :: pm_geomech

  if (.not.associated(pm_geomech)) then
    return
  endif

  select type(pm_geomech)
    class is (pm_geomech_force_type)
      option%igeommode = LINEAR_ELASTICITY_MODE
      option%geommode = "GEOMECHANICS"
      option%ngeomechdof = 3 ! displacements in x, y, z directions
      option%n_stress_strain_dof = 6
    class default
      option%io_buffer = 'Unrecognized geomechanics class in '// &
                          'InitSubsurfGeomechSetGeomechMode'
      call PrintErrMsg(option)
  end select

end subroutine InitSubsurfGeomechSetGeomechMode

! ************************************************************************** !

subroutine InitSubsurfGeomechChkInactiveCells(geomech_realization, &
                                             subsurf_realization)
  !
  ! checks if geomech nodes were mapped to inactive flow cells
  !
  ! Author: Glenn, Jumanah
  ! Date: 2/10/25
  !
  use Realization_Subsurface_class
  use Geomechanics_Realization_class
  use Geomechanics_Discretization_module
  use Geomechanics_Grid_Aux_module
  use GMDM_Pointer_module
  use Option_module

  implicit none

  class(realization_subsurface_type) :: subsurf_realization
  class(realization_geomech_type) :: geomech_realization

  type(option_type), pointer :: option
  type(gmdm_ptr_type), pointer :: dm_ptr

  PetscErrorCode :: ierr
  PetscBool :: error_found
  PetscInt :: geomech_local_id, subsurf_local_id, geomech_ghosted_id
  PetscInt :: subsurf_ghosted_id
  PetscReal, pointer :: subsurf_vec_1dof(:)

  error_found = PETSC_FALSE

  option => geomech_realization%option
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex( &
            geomech_realization%geomech_discretization, ONEDOF)

  call VecSet(geomech_realization%geomech_field%subsurf_vec_1dof,-777.d0, &
              ierr);CHKERRQ(ierr)
  call VecSet(geomech_realization%geomech_field%press,-888.d0, &
              ierr);CHKERRQ(ierr)
  call VecGetArray(geomech_realization%geomech_field%subsurf_vec_1dof, &
              subsurf_vec_1dof,ierr);CHKERRQ(ierr)
  do subsurf_local_id = 1, subsurf_realization%patch%grid%nlmax
    subsurf_ghosted_id = subsurf_realization%patch%grid%nL2G(subsurf_local_id)
    subsurf_vec_1dof(subsurf_local_id) = subsurf_realization%patch%imat( &
                                                  subsurf_ghosted_id)
  enddo
  call VecRestoreArray(geomech_realization%geomech_field%subsurf_vec_1dof, &
                          subsurf_vec_1dof,ierr);CHKERRQ(ierr)
  ! Scatter the data
  call VecScatterBegin(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                       geomech_realization%geomech_field%subsurf_vec_1dof, &
                       geomech_realization%geomech_field%press, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                     geomech_realization%geomech_field%subsurf_vec_1dof, &
                     geomech_realization%geomech_field%press, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecGetArray(geomech_realization%geomech_field%press, &
                      subsurf_vec_1dof, ierr);CHKERRQ(ierr)
  do geomech_local_id = 1, geomech_realization%geomech_patch%geomech_grid% &
                           nlmax_node
    geomech_ghosted_id = geomech_realization%geomech_patch%geomech_grid% &
                         nL2G(geomech_local_id)
    if (nint(subsurf_vec_1dof(geomech_local_id)) <= 0) error_found = PETSC_TRUE
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,error_found,ONE_INTEGER_MPI,MPI_C_BOOL, &
                     MPI_LOR,option%mycomm,ierr);CHKERRQ(ierr)
  if (error_found)then
    option%io_buffer = 'Cannot map inactive flow cell to geomechanics '//&
                       'node in the GEOMECHANICS_MAPPING_FILE! '
    call PrintErrMsg(option)
  endif

end subroutine InitSubsurfGeomechChkInactiveCells

! ************************************************************************** !

subroutine InitSubsurfGeomechSetupPMC(simulation,pm_geomech, &
                                     pmc_name,input)
  !
  ! refactored from factory_geomechanics.F90
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 2/10/25
  !
  use Realization_Subsurface_class
  use Option_module
  use Logging_module
  use Input_Aux_module
  use PM_Geomechanics_Force_class
  use Geomechanics_Realization_class
  use Timestepper_KSP_class
  use PMC_Geomechanics_class
  use Output_Aux_module
  use Waypoint_module
  use Simulation_Subsurface_class

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_geomech_force_type), pointer :: pm_geomech
  character(len=*) :: pmc_name
  type(input_type), pointer :: input

  character(len=MAXSTRINGLENGTH) :: string
  class(realization_subsurface_type), pointer :: subsurf_realization
  type(option_type), pointer :: option

  class(pmc_geomechanics_type), pointer :: pmc_geomech
  class(realization_geomech_type), pointer :: geomech_realization
  class(timestepper_ksp_type), pointer :: timestepper

  subsurf_realization => simulation%realization
  option => subsurf_realization%option
  subsurf_realization%output_option => simulation%output_option

  geomech_realization => GeomechRealizCreate(option)
  simulation%geomech%realization => geomech_realization

  input => InputCreate(IN_UNIT,option%input_filename,option)
  call InitSubsurfGeomechReadRequiredCards(geomech_realization,input)
  pmc_geomech => PMCGeomechanicsCreate()

  call pmc_geomech%SetName(pmc_name)
  call pmc_geomech%SetOption(option)
  simulation%geomech%process_model_coupler => pmc_geomech
  pmc_geomech%waypoint_list => simulation%waypoint_list_subsurface
  pmc_geomech%pm_list => pm_geomech
  pmc_geomech%pm_ptr%pm => pm_geomech
  pmc_geomech%geomech_realization => geomech_realization
  pm_geomech%geomech_realization => geomech_realization
  pmc_geomech%subsurf_realization => simulation%realization
  pm_geomech%subsurf_realization => simulation%realization

  ! add time integrator
  timestepper => TimestepperKSPCreate()
  pmc_geomech%timestepper => timestepper

  ! add solver
  call pm_geomech%InitializeSolver()
  timestepper%solver => pm_geomech%solver

  ! set up logging stage
  string = trim(pm_geomech%name)
  call LoggingCreateStage(string,pmc_geomech%stage)

  string = 'GEOMECHANICS'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)
  geomech_realization%output_option => &
    OutputOptionDuplicate(simulation%output_option)
  call OutputVariableListDestroy(geomech_realization%output_option% &
                                   output_snap_variable_list)
  call OutputVariableListDestroy(geomech_realization%output_option% &
                                   output_obs_variable_list)
  geomech_realization%output_option%output_snap_variable_list => &
    OutputVariableListCreate()
  geomech_realization%output_option%output_obs_variable_list => &
    OutputVariableListCreate()
  call InitSubsurfGeomechReadInput(simulation%geomech, &
                                   timestepper%solver, &
                                   input,option, &
                                   geomech_realization%output_option)
  pm_geomech%output_option => geomech_realization%output_option

  ! Hijack subsurface waypoint to geomechanics waypoint
  ! Subsurface controls the output now
  ! Always have snapshot on at t=0
  pmc_geomech%waypoint_list%first%print_snap_output = PETSC_TRUE

  ! link geomech and flow timestepper waypoints to geomech way point list
  if (associated(pmc_geomech)) then
    call pmc_geomech%SetWaypointPtr(pmc_geomech%waypoint_list)
    if (associated(simulation%flow_process_model_coupler)) then
      call simulation%flow_process_model_coupler% &
             SetWaypointPtr(pmc_geomech%waypoint_list)
    endif
  endif

  ! print the waypoints when debug flag is on
  if (geomech_realization%geomech_debug%print_waypoints) then
    call WaypointListPrint(pmc_geomech%waypoint_list,option, &
                           geomech_realization%output_option)
  endif

end subroutine InitSubsurfGeomechSetupPMC

! ************************************************************************** !

subroutine InitSubsurfGeomechReadSimBlock(input,pm)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/25/21
  !
  ! jaa: moved from factory_geomechanics.F90
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_ERT_class

  implicit none

  type(input_type), pointer :: input
  class(pm_base_type), pointer :: pm

  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  option => pm%option

  error_string = 'SIMULATION,PROCESS_MODELS,SUBSURFACE_GEOMECHANICS'

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case('OPTIONS')
        call pm%ReadSimulationOptionsBlock(input)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine InitSubsurfGeomechReadSimBlock

! ************************************************************************** !

subroutine GeomechReadDatasetToVecWithMask(geomech_realization,dataset, &
                                           geomech_material_id,read_all_values,vec,temp_vec)
  !
  ! Reads a geomechanics dataset into a PETSc Vec
  ! (based on SubsurfReadDatasetToVecWithMask)
  !
  ! Author: Kyle Mosley, WSP
  ! Date: 07/2025
  !
  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Patch_module
  use Option_module
  use Input_Aux_module
  Use HDF5_module
  Use Dataset_Base_class
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class

  implicit none

  ! input/output declarations
  class(realization_geomech_type) :: geomech_realization
  class(dataset_base_type) :: dataset
  PetscInt :: geomech_material_id
  PetscBool :: read_all_values
  Vec :: vec
  Vec :: temp_vec

  type(geomech_field_type), pointer :: field
  type(geomech_patch_type), pointer :: patch
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscInt :: local_id
  PetscErrorCode :: ierr
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: work_p(:)

  field => geomech_realization%geomech_field
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option

  call VecGetArray(vec,vec_p,ierr);CHKERRQ(ierr)
  if (index(dataset%filename,'.h5') > 0) then
    group_name = ''
    dataset_name = dataset%name
    select type(dataset)
    class is(dataset_common_hdf5_type)


      dataset_name = dataset%hdf5_dataset_name
      call HDF5ReadCellIndexedRealArrayGM(geomech_realization,temp_vec, &
                                        dataset%filename, &
                                        group_name,dataset_name, &
                                        dataset%realization_dependent)
      call VecGetArray(temp_vec,work_p,ierr);CHKERRQ(ierr)
      if (read_all_values) then
        do local_id = 1, grid%nlmax_node
          vec_p(local_id) = work_p(local_id)
        enddo
      else
        do local_id = 1, grid%nlmax_node
          if (patch%imat(grid%nL2G(local_id)) == geomech_material_id) then
            vec_p(local_id) = work_p(local_id)
          endif
        enddo
      endif
      call VecRestoreArray(temp_vec,work_p,ierr);CHKERRQ(ierr)
    class default
        option%io_buffer = 'Dataset "' // trim(dataset%name) // '" is of the &
          &wrong type for GeomechReadDatasetToVecWithMask()'
        call PrintErrMsg(option)
    end select
  endif
  call VecRestoreArray(vec,vec_p,ierr);CHKERRQ(ierr)

  end subroutine GeomechReadDatasetToVecWithMask

! ************************************************************************** !

end module Init_Subsurface_Geomech_module
