module TH_Well_module

#include "petsc/finclude/petscmat.h"
  use petscmat

  use PFLOTRAN_Constants_module
  use TH_Aux_module
  use PM_Well_class

  implicit none

  private

  PetscReal :: liquid_pressure = 1.d6
  PetscReal :: thermal_diffusivity = 0.143d-6

  PetscInt, parameter :: INLET_BC = 1
  PetscInt, parameter :: OUTLET_BC = 2

  public :: THWellSetup, &
            THWellUpdateAuxVars, &
            THWellPerturb, &
            THWell, &
            THWellAccumulationTerms, &
            THWellMatrixZeroing, &
            THWellUpdateHeatFlux

contains

! ************************************************************************** !

subroutine THWellSetup(pm_well_base,realization)
  !
  ! Configure the well
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25

  use Option_module
  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Field_module
  use Grid_module
  use Matrix_Zeroing_module
  use Grid_Unstructured_Aux_module, only : UCellGetApproxDXYZ

  implicit none

  class(pm_well_type), pointer :: pm_well_base
  class(realization_subsurface_type), pointer :: realization

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(connection_set_type), pointer :: connection_set
  type(connection_set_type), pointer :: bc_connection_set
  type(th_parameter_type), pointer :: th_parameter
  type(th_auxvar_type), pointer :: auxvars(:)
  type(th_well_auxvar_type), pointer :: auxvars_well(:)
  PetscInt, allocatable :: well_conn_to_grid_conn(:)
  PetscInt, allocatable :: global_well_cell_to_ghosted_id(:)
  PetscInt, allocatable :: local_well_cell_to_ghosted_id(:)
  PetscInt, allocatable :: ghosted_id_to_local_well_cell(:)
  PetscInt, allocatable :: ghosted_id_to_global_well_cell(:)
  PetscInt :: iconn_well, iconn_grid
  PetscReal :: surface_area
  PetscInt :: local_id
  PetscInt :: ghosted_id, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, iconn, idof
  PetscInt :: num_connections
  PetscInt :: num_ghosted_well_cells
  PetscReal :: pipe_cross_sectional_area
  PetscReal :: rock_thermal_conductivity
  PetscInt, allocatable :: inactive_cells(:)
  PetscReal, allocatable :: len_(:,:)
  PetscReal :: segment_length
  PetscReal :: dx, dy, dz
  PetscReal :: tempreal
  PetscErrorCode :: ierr
  PetscScalar, pointer :: vec_ptr(:)
  type(connection_set_type), pointer :: grid_connection_set
  PetscInt :: iwell_cell
  PetscInt :: local_well_cell_up, local_well_cell_dn
  PetscInt :: global_well_cell_up, global_well_cell_dn, global_well_cell
  PetscInt :: local_id_bc_in, local_id_bc_out
  PetscInt :: num_bc_connections

  if (.not.associated(pm_well_base)) return

  th_scale_by_volume = PETSC_FALSE
  pm_well => PMWellCastToClosedLoop(pm_well_base)

  field => pm_well%realization%field
  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  pm_well%realization_base => realization
  call pm_well%Setup()

  th_parameter => patch%aux%TH%th_parameter
  auxvars => patch%aux%TH%auxvars

  if (Uninitialized(pm_well%inlet_temperature)) then
    option%io_buffer = 'Well INLET_TEMPERATURE must be specified.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(pm_well%pipe_inner_diameter)) then
    option%io_buffer = 'Well PIPE_INNER_DIAMETER must be specified.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(pm_well%pipe_wall_thickness) .and. &
      Initialized(pm_well%pipe_wall_thermal_conductivity)) then
    option%io_buffer = 'Well PIPE_WALL_THICKNESS must be specified when &
      &PIPE_WALL_THERMAL_CONDUCTIVITY is specified.'
    call PrintErrMsg(option)
  endif
  if (Initialized(pm_well%pipe_wall_thickness) .and. &
      Uninitialized(pm_well%pipe_wall_thermal_conductivity)) then
    option%io_buffer = 'Well PIPE_WALL_THERMAL_CONDUCTIVITY must be &
      &specified when PIPE_WALL_THICKNESS is specified.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(pm_well%insulator_thickness) .and. &
      Initialized(pm_well%insulator_thermal_conductivity)) then
    option%io_buffer = 'Well INSULATOR_THICKNESS must be specified when &
      &INSULATOR_THERMAL_CONDUCTIVITY is specified.'
    call PrintErrMsg(option)
  endif
  if (Initialized(pm_well%insulator_thickness) .and. &
      Uninitialized(pm_well%insulator_thermal_conductivity)) then
    option%io_buffer = 'Well INSULATOR_THERMAL_CONDUCTIVITY must be &
      &specified when INSULATOR_THICKNESS is specified.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(pm_well%mass_flow_rate)) then
    option%io_buffer = 'Well PIPE_MASS_FLOW_RATE must be specified.'
    call PrintErrMsg(option)
  elseif (pm_well%mass_flow_rate < 0.d0) then
    option%io_buffer = 'Well PIPE_MASS_FLOW_RATE must be non-negative.'
    call PrintErrMsg(option)
  endif

  pipe_cross_sectional_area = PI*(0.5d0*pm_well%pipe_inner_diameter)**2

  allocate(len_(3,grid%ngmax))
  len_ = 0.d0
  local_id_bc_in = UNINITIALIZED_INTEGER
  local_id_bc_out = UNINITIALIZED_INTEGER
  num_bc_connections = 0
  num_ghosted_well_cells = 0
  num_connections = 0

  call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
  call VecGetArray(field%work,vec_ptr,ierr);CHKERRQ(ierr)
  do i = 1, size(pm_well%well_grid%h_local_id)
    local_id = pm_well%well_grid%h_local_id(i)
    if (local_id > 0) then
      vec_ptr(local_id) = dble(grid%nG2A(grid%nL2G(local_id)))
    endif
  enddo

  if (pm_well%well_grid%h_rank_id(1) == option%myrank) then
    local_id_bc_in = pm_well%well_grid%h_local_id(1)
  endif
  if (pm_well%well_grid%h_rank_id(pm_well%well_grid%nsegments) == &
      option%myrank) then
    local_id_bc_out = pm_well%well_grid%h_local_id( &
                        pm_well%well_grid%nsegments)
  endif

  call VecRestoreArray(field%work,vec_ptr,ierr);CHKERRQ(ierr)
  call pm_well%realization%comm1%GlobalToLocal(field%work,field%work_loc)
  call VecGetArray(field%work_loc,vec_ptr,ierr);CHKERRQ(ierr)
  allocate(global_well_cell_to_ghosted_id(pm_well%well_grid%nsegments))
  global_well_cell_to_ghosted_id = 0
  allocate(local_well_cell_to_ghosted_id(pm_well%well_grid%nsegments))
  local_well_cell_to_ghosted_id = 0
  allocate(ghosted_id_to_local_well_cell(grid%ngmax))
  ghosted_id_to_local_well_cell = 0
  allocate(ghosted_id_to_global_well_cell(grid%ngmax))
  ghosted_id_to_global_well_cell = 0
  num_ghosted_well_cells = 0
  do ghosted_id = 1, grid%ngmax
    if (vec_ptr(ghosted_id) > 0.d0) then
      num_ghosted_well_cells = num_ghosted_well_cells + 1
      local_well_cell_to_ghosted_id(num_ghosted_well_cells) = ghosted_id
      ghosted_id_to_local_well_cell(ghosted_id) = num_ghosted_well_cells
    endif
  enddo
  do iwell_cell = 1, pm_well%well_grid%nsegments
    if (pm_well%well_grid%h_ghosted_id(iwell_cell) > 0) then
      ghosted_id = pm_well%well_grid%h_ghosted_id(iwell_cell)
      global_well_cell_to_ghosted_id(iwell_cell) = ghosted_id
      ghosted_id_to_global_well_cell(ghosted_id) = iwell_cell
    endif
  enddo
  num_connections = 0
  grid_connection_set => grid%internal_connection_set_list%first
  do
    if (.not.associated(grid_connection_set)) exit
    do iconn = 1, grid_connection_set%num_connections
      ghosted_id_up = grid_connection_set%id_up(iconn)
      ghosted_id_dn = grid_connection_set%id_dn(iconn)
      if (vec_ptr(ghosted_id_up) > 0.d0 .and. &
          vec_ptr(ghosted_id_dn) > 0.d0) then
        num_connections = num_connections + 1
      endif
    enddo
    grid_connection_set => grid_connection_set%next
  enddo
  allocate(well_conn_to_grid_conn(max(num_connections,1)))
  well_conn_to_grid_conn = 0
  num_connections = 0
  grid_connection_set => grid%internal_connection_set_list%first
  do
    if (.not.associated(grid_connection_set)) exit
    do iconn = 1, grid_connection_set%num_connections
      ghosted_id_up = grid_connection_set%id_up(iconn)
      ghosted_id_dn = grid_connection_set%id_dn(iconn)
      if (vec_ptr(ghosted_id_up) > 0.d0 .and. &
          vec_ptr(ghosted_id_dn) > 0.d0) then
        num_connections = num_connections + 1
        well_conn_to_grid_conn(num_connections) = iconn
      endif
    enddo
    grid_connection_set => grid_connection_set%next
  enddo
  call VecRestoreArray(field%work_loc,vec_ptr,ierr);CHKERRQ(ierr)

  connection_set => &
    ConnectionCreate(num_connections,INTERNAL_FACE_CONNECTION_TYPE, &
                     NULL_GRID)

  grid_connection_set => grid%internal_connection_set_list%first
  do iconn_well = 1, num_connections
    iconn_grid = well_conn_to_grid_conn(iconn_well)
    ghosted_id_up = grid_connection_set%id_up(iconn_grid)
    ghosted_id_dn = grid_connection_set%id_dn(iconn_grid)
    local_well_cell_up = ghosted_id_to_local_well_cell(ghosted_id_up)
    local_well_cell_dn = ghosted_id_to_local_well_cell(ghosted_id_dn)
    global_well_cell_up = ghosted_id_to_global_well_cell(ghosted_id_up)
    global_well_cell_dn = ghosted_id_to_global_well_cell(ghosted_id_dn)
    if (global_well_cell_up < 1 .or. global_well_cell_dn < 1) then
      option%io_buffer = 'THWellSetup: well connection between cells that &
        &are not mapped to well segments on this process.'
      call PrintErrMsgByRank(option)
    endif
    ! Orient along the well path (increasing segment index = flow dir).
    if (global_well_cell_up > global_well_cell_dn) then
      i = local_well_cell_up
      local_well_cell_up = local_well_cell_dn
      local_well_cell_dn = i
      i = global_well_cell_up
      global_well_cell_up = global_well_cell_dn
      global_well_cell_dn = i
    endif
    connection_set%id_up(iconn_well) = local_well_cell_up
    connection_set%id_dn(iconn_well) = local_well_cell_dn
    dx = pm_well%well_grid%h(global_well_cell_dn)%x - &
         pm_well%well_grid%h(global_well_cell_up)%x
    dy = pm_well%well_grid%h(global_well_cell_dn)%y - &
         pm_well%well_grid%h(global_well_cell_up)%y
    dz = pm_well%well_grid%h(global_well_cell_dn)%z - &
         pm_well%well_grid%h(global_well_cell_up)%z
    tempreal = sqrt(dx*dx+dy*dy+dz*dz)
    connection_set%dist(-1,iconn_well) = UNINITIALIZED_DOUBLE
    connection_set%dist(0,iconn_well) = tempreal
    if (tempreal > 0.d0) then
      connection_set%dist(1,iconn_well) = dx/tempreal
      connection_set%dist(2,iconn_well) = dy/tempreal
      connection_set%dist(3,iconn_well) = dz/tempreal
    else
      connection_set%dist(1:3,iconn_well) = 0.d0
    endif
    connection_set%area(iconn_well) = pipe_cross_sectional_area
  enddo

  num_bc_connections = 0
  if (local_id_bc_in > 0) num_bc_connections = num_bc_connections + 1
  if (local_id_bc_out > 0) num_bc_connections = num_bc_connections + 1
  bc_connection_set => &
    ConnectionCreate(num_bc_connections,BOUNDARY_FACE_CONNECTION_TYPE, &
                     NULL_GRID)
  num_bc_connections = 0
  if (local_id_bc_in > 0) then
    num_bc_connections = num_bc_connections + 1
    i = 1
    dx = pm_well%well_grid%h(i)%x - pm_well%well_grid%tophole(1)
    dy = pm_well%well_grid%h(i)%y - pm_well%well_grid%tophole(2)
    dz = pm_well%well_grid%h(i)%z - pm_well%well_grid%tophole(3)
    tempreal = sqrt(dx*dx+dy*dy+dz*dz)
    bc_connection_set%id_dn(num_bc_connections) = &
      ghosted_id_to_local_well_cell(grid%nL2G(local_id_bc_in))
    bc_connection_set%dist(-1,num_bc_connections) = UNINITIALIZED_DOUBLE
    bc_connection_set%dist(0,num_bc_connections) = tempreal
    if (tempreal > 0.d0) then
      bc_connection_set%dist(1:3,num_bc_connections) = [dx,dy,dz]/tempreal
    else
      bc_connection_set%dist(1:3,num_bc_connections) = 0.d0
    endif
    bc_connection_set%area(num_bc_connections) = pipe_cross_sectional_area
  endif
  if (local_id_bc_out > 0) then
    num_bc_connections = num_bc_connections + 1
    i = pm_well%well_grid%nsegments
    dx = pm_well%well_grid%bottomhole(1) - pm_well%well_grid%h(i)%x
    dy = pm_well%well_grid%bottomhole(2) - pm_well%well_grid%h(i)%y
    dz = pm_well%well_grid%bottomhole(3) - pm_well%well_grid%h(i)%z
    tempreal = sqrt(dx*dx+dy*dy+dz*dz)
    bc_connection_set%id_dn(num_bc_connections) = &
      ghosted_id_to_local_well_cell(grid%nL2G(local_id_bc_out))
    bc_connection_set%dist(-1,num_bc_connections) = UNINITIALIZED_DOUBLE
    bc_connection_set%dist(0,num_bc_connections) = tempreal
    if (tempreal > 0.d0) then
      bc_connection_set%dist(1:3,num_bc_connections) = [dx,dy,dz]/tempreal
    else
      bc_connection_set%dist(1:3,num_bc_connections) = 0.d0
    endif
    bc_connection_set%area(num_bc_connections) = pipe_cross_sectional_area
  endif
  do i = 1, num_ghosted_well_cells
    ghosted_id = local_well_cell_to_ghosted_id(i)
    global_well_cell = ghosted_id_to_global_well_cell(ghosted_id)
    if (global_well_cell < 1) then
      option%io_buffer = 'THWellSetup: ghosted well cell is not mapped to &
        &a well segment. Ensure well segment data is available on all ranks.'
      call PrintErrMsgByRank(option)
    endif
    len_(:,ghosted_id) = [pm_well%well_grid%dx(global_well_cell), &
                          pm_well%well_grid%dy(global_well_cell), &
                          pm_well%well_grid%dz(global_well_cell)]
  enddo
  allocate(pm_well%well_cells_ghosted(num_ghosted_well_cells))
  do iwell_cell = 1, num_ghosted_well_cells
    ghosted_id = local_well_cell_to_ghosted_id(iwell_cell)
    pm_well%well_cells_ghosted(iwell_cell) = ghosted_id
  enddo
  allocate(pm_well%well_cells_bc_type(num_bc_connections))
  do i = 1, num_bc_connections
    if (local_id_bc_in > 0 .and. i == 1) then
      pm_well%well_cells_bc_type(i) = INLET_BC
    endif
    if (local_id_bc_out > 0 .and. &
        ((local_id_bc_in < 0 .and. i == 1) .or. &
         (local_id_bc_in > 0 .and. i == 2))) then
      pm_well%well_cells_bc_type(i) = OUTLET_BC
    endif
  enddo
  deallocate(global_well_cell_to_ghosted_id)
  deallocate(local_well_cell_to_ghosted_id)
  deallocate(ghosted_id_to_local_well_cell)
  deallocate(ghosted_id_to_global_well_cell)
  deallocate(well_conn_to_grid_conn)


  pm_well%connection_set => connection_set
  pm_well%bc_connection_set => bc_connection_set

  allocate(auxvars_well(num_ghosted_well_cells))
  do i = 1, num_ghosted_well_cells
    ghosted_id = pm_well%well_cells_ghosted(i)
    auxvars(ghosted_id)%iwellaux = i
    call THWellAuxVarInit(auxvars_well(i),option, &
                          th_numerical_derivatives)
    auxvars_well(i)%z = grid%z(ghosted_id)
    auxvars_well(i)%ghosted_id = ghosted_id
    segment_length = sqrt(len_(1,ghosted_id)**2 + &
                          len_(2,ghosted_id)**2 + &
                          len_(3,ghosted_id)**2)
    auxvars_well(i)%segment_length = segment_length
    auxvars_well(i)%volume = segment_length * &
                             pipe_cross_sectional_area
    surface_area = segment_length * pm_well%pipe_inner_diameter * PI
    rock_thermal_conductivity = th_parameter%ckwet(patch%cct_id(ghosted_id))
    select case(grid%itype)
      case(STRUCTURED_GRID)
        dx = grid%structured_grid%dx(ghosted_id)
        dy = grid%structured_grid%dy(ghosted_id)
        dz = grid%structured_grid%dz(ghosted_id)
      case(IMPLICIT_UNSTRUCTURED_GRID)
        call UCellGetApproxDXYZ(grid%unstructured_grid,ghosted_id, &
                                option,dx,dy,dz)
      case default
        option%io_buffer = 'Unsupported grid type in THWellSetup.'
        call PrintErrMsg(option)
    end select
    call THWellWI(rock_thermal_conductivity, &
                  rock_thermal_conductivity, &
                  rock_thermal_conductivity, &
                  dx,dy,dz, &
                  len_(1,ghosted_id), &
                  len_(2,ghosted_id), &
                  len_(3,ghosted_id), &
                  pm_well%pipe_inner_diameter, &
                  0.d0, & ! s
                  auxvars_well(i)%well_index)
    auxvars_well(i)%therm_cond_borehole_to_cell = auxvars_well(i)%well_index / &
                                                  surface_area
    if (th_numerical_derivatives) then
      call THWellAuxVarCopyParamsToPert(auxvars_well(i))
    endif
  enddo
  deallocate(len_)

  patch%aux%TH%auxvars_well => auxvars_well
  patch%aux%TH%num_well_aux = num_ghosted_well_cells
  allocate(auxvars_well(num_bc_connections))
  do i = 1, num_bc_connections
    call THWellAuxVarInit(auxvars_well(i),option,th_numerical_derivatives)
  enddo
  patch%aux%TH%auxvars_well_bc => auxvars_well
  patch%aux%TH%num_well_aux_bc = num_bc_connections

  pm_well%matrix_zeroing => MatrixZeroingCreate()
  allocate(inactive_cells(grid%nlmax))
  inactive_cells = 1
  do i = 1, size(pm_well%well_cells_ghosted)
    ghosted_id = pm_well%well_cells_ghosted(i)
    local_id = grid%nG2L(ghosted_id)
    if (local_id > 0) inactive_cells(local_id) = 0
  enddo
  i = 0
  do local_id = 1, grid%nlmax
    if (inactive_cells(local_id) == 1) then
      i = i + 1
      inactive_cells(i) = local_id
    endif
  enddo
  if (i < grid%nlmax) then
    inactive_cells(i+1:) = UNINITIALIZED_INTEGER
  endif
  call MatrixZeroingAllocateArray(pm_well%matrix_zeroing,i,option)
  do i = 1, pm_well%matrix_zeroing%n_zero_rows
    local_id = inactive_cells(i)
    ghosted_id = grid%nL2G(local_id)
    idof = (local_id-1)*option%nflowdof + th_well_dof
    pm_well%matrix_zeroing%zero_rows_local(i) = idof ! one-based
    idof = (ghosted_id-1)*option%nflowdof + th_well_dof
    pm_well%matrix_zeroing%zero_rows_local_ghosted(i) = idof-1 ! zero-based
  enddo
  deallocate(inactive_cells)

end subroutine THWellSetup

! ************************************************************************** !

subroutine THWellUpdateAuxVars(xx,pm_well_base)
  !
  ! Configure the well
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25

  use Option_module
  use Connection_module
  use Patch_module
  use Grid_module

  implicit none

  PetscReal :: xx(:) ! ghosted solution vec
  class(pm_well_type), pointer :: pm_well_base

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(th_well_auxvar_type), pointer :: auxvars_well(:)
  type(th_well_auxvar_type), pointer :: auxvars_well_bc(:)
  PetscInt :: ghosted_id, natural_id
  PetscInt :: i, ibc
  PetscInt :: istart, iend
  PetscInt :: num_well_cells
  PetscReal :: xxbc(3)

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  auxvars_well => patch%aux%TH%auxvars_well
  auxvars_well_bc => patch%aux%TH%auxvars_well_bc

  num_well_cells = size(pm_well%well_cells_ghosted)
  do i = 1, num_well_cells
    ghosted_id = pm_well%well_cells_ghosted(i)
    natural_id = grid%nG2A(ghosted_id)
    iend = ghosted_id * option%nflowdof
    istart = iend - option%nflowdof + 1
    call THWellAuxVarCompute(liquid_pressure,xx(istart:iend), &
                             auxvars_well(i), &
                             natural_id,option)
  enddo

  xxbc = UNINITIALIZED_DOUBLE
  do ibc = 1, size(auxvars_well_bc)
    if (pm_well%well_cells_bc_type(ibc) == INLET_BC) then
      xxbc(th_well_dof) = pm_well%inlet_temperature
    elseif (pm_well%well_cells_bc_type(ibc) == OUTLET_BC) then
      ! to guarantee zero conduction at boundary
      xxbc(th_well_dof) = auxvars_well(pm_well%bc_connection_set%id_dn(ibc))%temp
    endif
    call THWellAuxVarCompute(liquid_pressure,xxbc, &
                             auxvars_well_bc(ibc), &
                             UNINITIALIZED_INTEGER,option)
  enddo

end subroutine THWellUpdateAuxVars

! ************************************************************************** !

subroutine THWellAccumulationTerms(r,J,pm_well_base,calculate_derivatives)
  !
  ! Calculates contribution of well model to the TH residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25
  !
  use Option_module
  use Petsc_Utility_module
  use Connection_module
  use Patch_module
  use Grid_module

  implicit none

  Vec :: r
  Mat :: J
  class(pm_well_type), pointer :: pm_well_base
  PetscBool :: calculate_derivatives

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  PetscReal, pointer :: r_p(:)
  type(th_well_auxvar_type), pointer :: auxvars_well(:)
  PetscReal :: Res(3)
  PetscReal :: Jblock(3,3)
  PetscReal :: dt
  PetscInt :: ndof
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: i
  PetscErrorCode :: ierr

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  ndof = option%nflowdof
  auxvars_well => patch%aux%TH%auxvars_well
  dt = option%flow_dt

  if (.not.calculate_derivatives) then
    call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  endif

  do i = 1, size(pm_well%well_cells_ghosted)
    ghosted_id = pm_well%well_cells_ghosted(i)
    local_id = grid%nG2L(ghosted_id)
    ! only owned cells contribute to residual/Jacobian rows
    if (local_id <= 0) cycle
    call THWellAccumulationDerivative(auxvars_well(i),Res,Jblock, &
                                      ndof,calculate_derivatives)
    Res = Res/dt
    Jblock = Jblock/dt
    if (calculate_derivatives) then
      call PetUtilMatSVBL(J,ghosted_id,ghosted_id,Jblock,ndof)
    else
      call PetUtilVecSVBL(r_p,local_id,Res,ndof,PETSC_FALSE)
    endif
  enddo

  ! calculate heat transfer between well cells and subsurface cells
  if (.not.calculate_derivatives) then
    call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine THWellAccumulationTerms

! ************************************************************************** !

subroutine THWell(r,J,pm_well_base,calculate_derivatives)
  !
  ! Calculates contribution of well model to the TH residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25
  !
  use Option_module
  use Petsc_Utility_module
  use Connection_module
  use Patch_module
  use Grid_module

  implicit none

  Vec :: r
  Mat :: J
  class(pm_well_type), pointer :: pm_well_base
  PetscBool :: calculate_derivatives

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  PetscReal, pointer :: r_p(:)
  type(th_auxvar_type), pointer :: auxvars(:)
  type(th_well_auxvar_type), pointer :: auxvars_well(:)
  type(th_well_auxvar_type), pointer :: auxvars_well_bc(:)
  type(connection_set_type), pointer :: connection_set
  type(th_well_auxvar_type), pointer :: auxvar_well_up, auxvar_well_dn
  PetscReal :: Res(3)
  PetscReal :: Jup(3,3), Jdn(3,3)
  PetscInt :: ndof
  PetscInt :: local_id, local_id_up, local_id_dn
  PetscInt :: ghosted_id, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, iconn, ibc
  PetscInt :: ghosted_well_id, ghosted_well_id_up, ghosted_well_id_dn
  PetscReal :: velocity_scale
  PetscErrorCode :: ierr

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  ndof = option%nflowdof
  auxvars => patch%aux%TH%auxvars
  auxvars_well => patch%aux%TH%auxvars_well
  auxvars_well_bc => patch%aux%TH%auxvars_well_bc

  call THWellAccumulationTerms(r,J,pm_well_base,calculate_derivatives)

  if (.not.calculate_derivatives) then
    call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  endif

#if 1
  ! calculate heat transfer between well cells
  connection_set => pm_well%connection_set
  i = 0
  do iconn = 1, connection_set%num_connections
    ! cells are ghosted in connection set
    ghosted_well_id_up = connection_set%id_up(iconn)
    ghosted_well_id_dn = connection_set%id_dn(iconn)
    i = i + 1
    call THWellFluxDerivative(auxvars_well(ghosted_well_id_up), &
                    auxvars_well(ghosted_well_id_dn), &
                    pm_well%mass_flow_rate/FMWH2O, &
                    connection_set%dist(0,iconn), &
                    connection_set%area(iconn), &
                    Res,Jup,Jdn,ndof,calculate_derivatives, &
                    PETSC_FALSE)
    ghosted_id_up = auxvars_well(ghosted_well_id_up)%ghosted_id
    ghosted_id_dn = auxvars_well(ghosted_well_id_dn)%ghosted_id
    local_id_up = grid%nG2L(ghosted_id_up)
    local_id_dn = grid%nG2L(ghosted_id_dn)
    if (calculate_derivatives) then
      if (local_id_up > 0) then
        call PetUtilMatSVBL(J,ghosted_id_up,ghosted_id_up,Jup,ndof)
        call PetUtilMatSVBL(J,ghosted_id_up,ghosted_id_dn,Jdn,ndof)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
        call PetUtilMatSVBL(J,ghosted_id_dn,ghosted_id_dn,Jdn,ndof)
        call PetUtilMatSVBL(J,ghosted_id_dn,ghosted_id_up,Jup,ndof)
      endif
    else
      if (local_id_up > 0) then
        call PetUtilVecSVBL(r_p,local_id_up,Res,ndof,PETSC_FALSE)
      endif
      if (local_id_dn > 0) then
        Res = -Res
        call PetUtilVecSVBL(r_p,local_id_dn,Res,ndof,PETSC_FALSE)
      endif
    endif

  enddo
#endif

#if 1
  connection_set => pm_well%bc_connection_set
  do ibc = 1, connection_set%num_connections
    auxvar_well_up => auxvars_well_bc(ibc)
    ghosted_well_id = connection_set%id_dn(ibc)
    select case(pm_well%well_cells_bc_type(ibc))
      case(INLET_BC)
        auxvar_well_dn => auxvars_well(ghosted_well_id)
        velocity_scale = 1.d0
      case(OUTLET_BC) ! outflow
        ! note that the temperature is currently equal to the last well cell,
        ! so this is effectively a zero conduction bc
        auxvar_well_dn => auxvars_well(ghosted_well_id)
        velocity_scale = -1.d0
      case default
        cycle
    end select
    call THWellFluxDerivative(auxvar_well_up, &
                    auxvar_well_dn, &
                    pm_well%mass_flow_rate/FMWH2O*velocity_scale, &
                    connection_set%dist(0,ibc), &
                    connection_set%area(ibc), &
                    Res,Jup,Jdn,ndof,calculate_derivatives, &
                    PETSC_TRUE)
    ghosted_id_dn = auxvar_well_dn%ghosted_id
    local_id_dn = grid%nG2L(ghosted_id_dn)
    if (local_id_dn <= 0) cycle
    if (calculate_derivatives) then
      Jdn = -Jdn
      call PetUtilMatSVBL(J,ghosted_id_dn,ghosted_id_dn,Jdn,ndof)
    else
      Res = -Res
      pm_well%heat_flux(1,ibc) = -Res(3)
      call PetUtilVecSVBL(r_p,local_id_dn,Res,ndof,PETSC_FALSE)
    endif
  enddo
#endif

#if 1
  do i = 1, size(auxvars_well)
    ghosted_id = auxvars_well(i)%ghosted_id
    local_id = grid%nG2L(ghosted_id)
    if (local_id <= 0) cycle
    call THWellHeatExchangeDerivative(auxvars(ghosted_id), &
                                      auxvars_well(i), &
                                      pm_well, &
                                      Res,Jdn, &
                                      ndof,calculate_derivatives)
    if (calculate_derivatives) then
      call PetUtilMatSVBL(J,ghosted_id,ghosted_id,Jdn,ndof)
    else
      call PetUtilVecSVBL(r_p,local_id,Res,ndof,PETSC_FALSE)
    endif
  enddo
#endif

  ! calculate heat transfer between well cells and subsurface cells
  if (.not.calculate_derivatives) then
    call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine THWell

! ************************************************************************** !

subroutine THWellMatrixZeroing(pm_well_base,r,J,calculate_derivatives)
  !
  ! Zeros dof for non-well grid cells
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/25
  !
  use Matrix_Zeroing_module

  implicit none

  class(pm_well_type), pointer :: pm_well_base
  Vec :: r
  Mat :: J
  PetscBool :: calculate_derivatives

  class(pm_well_closed_loop_type), pointer :: pm_well

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  if (calculate_derivatives) then
    call MatrixZeroingZeroMatEntries(pm_well%matrix_zeroing,J)
  else
    call MatrixZeroingZeroVecEntries(pm_well%matrix_zeroing,r)
  endif

end subroutine THWellMatrixZeroing

! ************************************************************************** !

subroutine THWellPerturb(pm_well_base)
  !
  ! Configure the well
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25

  use Option_module
  use Patch_module
  use Grid_module

  implicit none

  class(pm_well_type), pointer :: pm_well_base

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(th_well_auxvar_type), pointer :: auxvars_well(:)
  PetscInt :: ghosted_id, natural_id
  PetscInt :: i
  PetscInt :: num_well_cells

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  auxvars_well => patch%aux%TH%auxvars_well

  num_well_cells = size(pm_well%well_cells_ghosted)
  do i = 1, num_well_cells
    ghosted_id = pm_well%well_cells_ghosted(i)
    natural_id = grid%nG2A(ghosted_id)
    call THWellAuxVarPerturb(liquid_pressure,auxvars_well(i), &
                             natural_id,option)
  enddo

end subroutine THWellPerturb

! ************************************************************************** !

subroutine THWellAccumulation(auxvar_well,Res,Jac,ndof)
  !
  ! Calculates accumulation term for TH well residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  type(th_well_auxvar_type) :: auxvar_well
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: Jac(ndof,ndof)

  PetscReal :: tempreal

  ! Res [MJ] = liquid_density_kmol [kmol liquid/m^3 liquid] *
  !            internal energy [MJ/kmol liquid] *
  !            volume [m^3 liquid]
  Res = 0.d0
  Jac = 0.d0
  tempreal = auxvar_well%volume * auxvar_well%den
  Res(ndof) = tempreal * auxvar_well%u
  Jac(ndof,ndof) = tempreal * auxvar_well%du_dT

end subroutine THWellAccumulation

! ************************************************************************** !

subroutine THWellAccumulationDerivative(auxvar_well,Res,Jac,ndof, &
                                        calculate_derivatives)
  !
  ! Calculates accumulation term for TH well residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  type(th_well_auxvar_type) :: auxvar_well
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: Jac(ndof,ndof)
  PetscBool :: calculate_derivatives

  PetscReal :: Res_pert(ndof)
  PetscReal :: Jdum(ndof,ndof)
  PetscInt :: idof, ieq

  call THWellAccumulation(auxvar_well,Res,Jac,ndof)
  if (calculate_derivatives .and. th_numerical_derivatives) then
    do idof = 1, ndof
      call THWellAccumulation(auxvar_well%auxvar_pert(idof), &
                              Res_pert,Jdum,ndof)
      do ieq = 1, ndof
        Jac(ieq,idof) = (Res_pert(ieq) - Res(ieq)) / &
                        auxvar_well%auxvar_pert(idof)%pert
      enddo
    enddo
  endif

end subroutine THWellAccumulationDerivative

! ************************************************************************** !

subroutine THWellFlux(auxvar_well_up,auxvar_well_dn, &
                      q_kmol,dist,area,Res,Jup,Jdn,ndof, &
                      calculate_derivatives)
  !
  ! Calculates internal flux term for TH well residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  type(th_well_auxvar_type) :: auxvar_well_up, auxvar_well_dn
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: Jup(ndof,ndof)
  PetscReal :: Jdn(ndof,ndof)
  PetscReal :: q_kmol ! [kmol/s]
  PetscReal :: dist
  PetscReal :: area
  PetscBool :: calculate_derivatives

  PetscReal :: tup, tdn
  PetscReal :: dh_dTup, dh_dTdn
  PetscReal :: enthalpy
  PetscReal :: diffusive_energy_flux, convective_energy_flux
  PetscReal :: tempreal

  Res = 0.d0

  tup = auxvar_well_up%temp
  tdn = auxvar_well_dn%temp

  ! convection
  dh_dTup = 0.d0
  dh_dTdn = 0.d0
  if (q_kmol > 0.d0) then
    enthalpy = auxvar_well_up%h
    dh_dTup = auxvar_well_up%dh_dT
  else
    enthalpy = auxvar_well_dn%h
    dh_dTdn = auxvar_well_dn%dh_dT
  endif

  convective_energy_flux = q_kmol * enthalpy

  ! conduction
  ! [MW] = [MW/m-K] * [m^2 bulk] [C] / [m bulk]
  tempreal = thermal_diffusivity * area / dist
  diffusive_energy_flux = tempreal * (tup-tdn)

  Res(th_well_eq) = convective_energy_flux + diffusive_energy_flux

  if (calculate_derivatives) then
    Jup = 0.d0
    Jup(th_well_eq,th_well_dof) = q_kmol * dh_dTup + tempreal
    Jdn = 0.d0
    Jdn(th_well_eq,th_well_dof) = q_kmol * dh_dTdn - tempreal
  endif

end subroutine THWellFlux

! ************************************************************************** !

subroutine THWellFluxDerivative(auxvar_well_up,auxvar_well_dn, &
                                q_kmol,dist,area,Res,Jup,Jdn,ndof, &
                                calculate_derivatives,is_bc)
  !
  ! Calculates internal flux term for TH well residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  type(th_well_auxvar_type) :: auxvar_well_up, auxvar_well_dn
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: Jup(ndof,ndof)
  PetscReal :: Jdn(ndof,ndof)
  PetscReal :: q_kmol ! [kmol/s]
  PetscReal :: dist
  PetscReal :: area
  PetscBool :: calculate_derivatives
  PetscBool :: is_bc

  PetscReal :: Res_pert(ndof)
  PetscReal :: Jdum(ndof,ndof), Jdum1(ndof,ndof)
  PetscInt :: idof, ieq

  call THWellFlux(auxvar_well_up,auxvar_well_dn, &
                  q_kmol,dist,area,Res,Jup,Jdn,ndof, &
                  .not.th_numerical_derivatives)
  if (calculate_derivatives .and. th_numerical_derivatives) then
    Jup = 0.d0
    if (.not.is_bc) then
      do idof = 1, ndof
        call THWellFlux(auxvar_well_up%auxvar_pert(idof),auxvar_well_dn, &
                        q_kmol,dist,area,Res_pert,Jdum,Jdum1,ndof,PETSC_FALSE)
        do ieq = 1, ndof
          Jup(ieq,idof) = (Res_pert(ieq) - Res(ieq)) / &
                          auxvar_well_up%auxvar_pert(idof)%pert
        enddo
      enddo
    endif
    Jdn = 0.d0
    do idof = 1, ndof
      call THWellFlux(auxvar_well_up,auxvar_well_dn%auxvar_pert(idof), &
                      q_kmol,dist,area,Res_pert,Jdum,Jdum1,ndof,PETSC_FALSE)
      do ieq = 1, ndof
        Jdn(ieq,idof) = (Res_pert(ieq) - Res(ieq)) / &
                        auxvar_well_dn%auxvar_pert(idof)%pert
      enddo
    enddo
  endif

end subroutine THWellFluxDerivative

! ************************************************************************** !

function THWellResistanceWellIndex(auxvar_well,pm_well)
  !
  ! Calculates heat flux between pipe wall and insulator
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/26
  !
  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well

  PetscReal :: THWellResistanceWellIndex

  if (pm_well%use_well_index) then
    ! [K/MW] = 1 / [MW/K]
    THWellResistanceWellIndex = 1.d0 / auxvar_well%well_index
  else
    THWellResistanceWellIndex = 0.d0
  endif

end function THWellResistanceWellIndex

! ************************************************************************** !

function THWellResistanceInsulator(auxvar_well,pm_well)
  !
  ! Calculates heat flux between pipe wall and insulator
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/26
  !
  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well

  PetscReal :: THWellResistanceInsulator

  PetscReal :: r_inner
  PetscReal :: r_outer

  r_inner = 0.5d0 * pm_well%pipe_inner_diameter + pm_well%pipe_wall_thickness
  r_outer = r_inner + pm_well%insulator_thickness

  if (pm_well%insulator_thermal_conductivity > 0.d0 .and. &
       r_inner > 0.d0 .and. auxvar_well%segment_length > 0.d0) then
    ! Cylindrical conduction resistance: R = ln(r_o/r_i)/(2*pi*k*L) [K/MW]
    THWellResistanceInsulator = &
      log(r_outer/r_inner) / (2.d0 * PI * &
      pm_well%insulator_thermal_conductivity * auxvar_well%segment_length)
  elseif (pm_well%insulator_thermal_conductivity < 0.d0) then
    THWellResistanceInsulator = 0.d0
  else
    THWellResistanceInsulator = MAX_DOUBLE
  endif

end function THWellResistanceInsulator

! ************************************************************************** !

function THWellResistancePipeWall(auxvar_well,pm_well)
  !
  ! Calculates heat flux between fluid centerline and pipe wall
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/26
  !
  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well

  PetscReal :: THWellResistancePipeWall

  PetscReal :: r_inner
  PetscReal :: r_outer

  r_inner = 0.5d0 * pm_well%pipe_inner_diameter
  r_outer = r_inner + pm_well%pipe_wall_thickness

  if (pm_well%pipe_wall_thermal_conductivity > 0.d0 .and. &
      r_inner > 0.d0 .and. auxvar_well%segment_length > 0.d0) then
    ! Cylindrical conduction resistance: R = ln(r_o/r_i)/(2*pi*k*L) [K/MW]
    THWellResistancePipeWall = &
      log(r_outer/r_inner) / (2.d0 * PI * &
      pm_well%pipe_wall_thermal_conductivity * auxvar_well%segment_length)
  elseif (pm_well%pipe_wall_thermal_conductivity < 0.d0) then
    THWellResistancePipeWall = 0.d0
  else
    THWellResistancePipeWall = MAX_DOUBLE
  endif

end function THWellResistancePipeWall

! ************************************************************************** !

function THWellResistanceFluid(auxvar_well,pm_well)
  !
  ! Calculates heat flux between fluid centerline and pipe wall
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/26
  !
  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well

  PetscReal :: THWellResistanceFluid
  PetscReal :: thermal_coeff_fluid
  PetscReal :: surface_area

  thermal_coeff_fluid = THWellComputeConvCoefficient(auxvar_well,pm_well)
  ! Inner wall surface area for convection
  surface_area = PI * pm_well%pipe_inner_diameter * auxvar_well%segment_length

  ! Convective resistance: R = 1/(h*A) [K/MW] as h [MW/m^2-K] and A [m^2]
  if (thermal_coeff_fluid > 0.d0 .and. surface_area > 0.d0) then
    THWellResistanceFluid = 1.d0 / (thermal_coeff_fluid * surface_area)
  elseif (thermal_coeff_fluid < 0.d0) then
    THWellResistanceFluid = 0.d0
  else
    THWellResistanceFluid = MAX_DOUBLE
  endif

end function THWellResistanceFluid

! ************************************************************************** !

subroutine THWellHeatExchange(auxvar,auxvar_well,pm_well, &
                              Res,J,ndof, &
                              calculate_derivatives)
  !
  ! Calculates heat flux between grid cell and pipe
  !
  ! Author: Glenn Hammond
  ! Date: 09/18/25
  !
  type(th_auxvar_type) :: auxvar
  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: J(ndof,ndof)
  PetscBool :: calculate_derivatives

  PetscReal :: tcell, tpipe
  PetscReal :: diffusive_energy_flux
  PetscReal :: total_thermal_resistance
  PetscReal :: tempreal

  Res = 0.d0

  tcell = auxvar%temp
  tpipe = auxvar_well%temp

  total_thermal_resistance = &  ! [K/MW]
    THWellResistanceFluid(auxvar_well,pm_well) + &
    THWellResistancePipeWall(auxvar_well,pm_well) + &
    THWellResistanceInsulator(auxvar_well,pm_well) + &
    THWellResistanceWellIndex(auxvar_well,pm_well)

  diffusive_energy_flux = (tcell-tpipe) / total_thermal_resistance

  Res(TH_ENERGY_EQUATION_INDEX) = diffusive_energy_flux
  Res(th_well_eq) = -diffusive_energy_flux

  if (calculate_derivatives) then
    J = 0.d0
    tempreal = 1.d0 / total_thermal_resistance
    J(TH_ENERGY_EQUATION_INDEX,th_well_dof) = -tempreal
    J(TH_ENERGY_EQUATION_INDEX,TH_TEMPERATURE_DOF) = tempreal
    J(th_well_eq,th_well_dof) = tempreal
    J(th_well_eq,TH_TEMPERATURE_DOF) = -tempreal
  endif

end subroutine THWellHeatExchange

! ************************************************************************** !
subroutine THWellHeatExchangeDerivative(auxvar,auxvar_well,pm_well, &
                                        Res,J,ndof, &
                                        calculate_derivatives)
  !
  ! Calculates heat flux between grid cell and pipe
  !
  ! Author: Glenn Hammond
  ! Date: 09/18/25
  !
  type(th_auxvar_type) :: auxvar
  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: J(ndof,ndof)
  PetscBool :: calculate_derivatives

  PetscReal :: Res_pert(ndof)
  PetscReal :: Jdum(ndof,ndof)
  PetscReal :: pert
  PetscInt :: idof, ieq

  call THWellHeatExchange(auxvar,auxvar_well,pm_well, &
                          Res,J,ndof, &
                          .not.th_numerical_derivatives)
  if (calculate_derivatives .and. th_numerical_derivatives) then
    J = 0.d0
    do idof = 1, ndof
      call THWellHeatExchange(auxvar%auxvar_pert(idof), &
                              auxvar_well%auxvar_pert(idof),pm_well, &
                              Res_pert,Jdum,ndof,PETSC_FALSE)
      if (idof == th_well_dof) then
        pert = auxvar_well%auxvar_pert(idof)%pert
      else
        pert = auxvar%auxvar_pert(idof)%pert
      endif
      do ieq = 1, ndof
        J(ieq,idof) = (Res_pert(ieq) - Res(ieq)) / pert
      enddo
    enddo
  endif

end subroutine THWellHeatExchangeDerivative

! ************************************************************************** !

subroutine THWellWI(kx,ky,kz,dx,dy,dz,lenx,leny,lenz,pipe_diameter,s,wi)
  !
  ! Calculates well index for a well segment
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/25
  !
  implicit none

  PetscReal :: kx   ! conductivity in x
  PetscReal :: ky
  PetscReal :: kz
  PetscReal :: dx   ! grid spacing in x
  PetscReal :: dy
  PetscReal :: dz
  PetscReal :: lenx ! length of well segement in x
  PetscReal :: leny
  PetscReal :: lenz
  PetscReal :: pipe_diameter
  PetscReal :: s
  PetscReal :: wi

  PetscReal :: sqrt_kx_over_ky
  PetscReal :: sqrt_ky_over_kx
  PetscReal :: sqrt_ky_over_kz
  PetscReal :: sqrt_kz_over_ky
  PetscReal :: sqrt_kx_over_kz
  PetscReal :: sqrt_kz_over_kx
  PetscReal :: ro_x
  PetscReal :: ro_y
  PetscReal :: ro_z
  PetscReal :: wi_x
  PetscReal :: wi_y
  PetscReal :: wi_z
  PetscReal :: rw

  rw = 0.5d0 * pipe_diameter

  sqrt_kx_over_ky = sqrt(kx/ky)
  sqrt_ky_over_kx = sqrt(ky/kx)
  sqrt_ky_over_kz = sqrt(ky/kz)
  sqrt_kz_over_ky = sqrt(kz/ky)
  sqrt_kx_over_kz = sqrt(kx/kz)
  sqrt_kz_over_kx = sqrt(kz/kx)

  ro_x = 0.28d0 * sqrt(sqrt_ky_over_kz*dz**2 + sqrt_kz_over_ky*dy**2) / &
                  (sqrt(sqrt_ky_over_kz) + sqrt(sqrt_kz_over_ky))
  ro_y = 0.28d0 * sqrt(sqrt_kz_over_kx*dx**2 + sqrt_kx_over_kz*dz**2) / &
                  (sqrt(sqrt_kz_over_kx) + sqrt(sqrt_kx_over_kz))
  ro_z = 0.28d0 * sqrt(sqrt_ky_over_kx*dx**2 + sqrt_kx_over_ky*dy**2) / &
                  (sqrt(sqrt_ky_over_kx) + sqrt(sqrt_kx_over_ky))

  wi_x = 2.d0*PI*sqrt(ky*kz)*lenx/(log(ro_x/rw)+s)
  wi_y = 2.d0*PI*sqrt(kx*kz)*leny/(log(ro_y/rw)+s)
  wi_z = 2.d0*PI*sqrt(kx*ky)*lenz/(log(ro_z/rw)+s)

  wi = sqrt(wi_x**2 + wi_y**2 + wi_z**2)

end subroutine THWellWI

! ************************************************************************** !

subroutine THWellComputePipeCoefficient(auxvar_well,pm_well)
  !
  ! Computes pipe wall conduction coefficient for cylindrical wall
  ! h_pipe = k_wall / (r_outer * ln(r_outer/r_inner))
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/18/26
  !
  implicit none

  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well

  PetscReal :: thermal_cond_pipe_wall ! wall thermal conductivity [MW/m-K]
  PetscReal :: r_inner                ! inner radius [m]
  PetscReal :: r_outer                ! outer radius [m]
  PetscReal :: h_pipe_wall            ! [MW/m^2-K]

  thermal_cond_pipe_wall = pm_well%pipe_wall_thermal_conductivity
  r_inner = 0.5d0 * pm_well%pipe_inner_diameter
  r_outer = r_inner + pm_well%pipe_wall_thickness

  if (r_outer > r_inner .and. r_inner > 0.d0) then
    ! For cylindrical wall: R_wall = ln(r_o/r_i)/(2*pi*L*k)
    ! Per unit area (dividing by 2*pi*r_o*L): h = k/(r_o*ln(r_o/r_i))
    h_pipe_wall = thermal_cond_pipe_wall / (r_outer * log(r_outer/r_inner))
  else if (r_outer <= r_inner) then
    ! No wall thickness - infinite conductivity (no resistance)
    h_pipe_wall = 1.d20
  else
    h_pipe_wall = 0.d0
  endif

end subroutine THWellComputePipeCoefficient

! ************************************************************************** !

function THWellComputeConvCoefficient(auxvar_well,pm_well)
  !
  ! Computes convection heat transfer coefficient from fluid flow in well
  ! h = Nu * k / D
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/18/26
  !
  implicit none

  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well

  PetscReal :: density
  PetscReal :: velocity
  PetscReal :: pipe_diameter
  PetscReal :: viscosity
  PetscReal :: heat_capacity
  PetscReal :: thermal_cond
  PetscReal :: THWellComputeConvCoefficient

  PetscReal :: reynolds
  PetscReal :: prandtl
  PetscReal :: nusselt

  reynolds = 0.d0
  prandtl = 0.d0
  nusselt = 0.d0

  pipe_diameter = pm_well%pipe_inner_diameter

  density = auxvar_well%den_kg
  viscosity = auxvar_well%vis
  heat_capacity = auxvar_well%spec_heat_fluid
  thermal_cond = auxvar_well%therm_cond_fluid

  ! local mean velocity derived from the constant mass flow rate
  ! [m/s] = [kg/s] / ([kg/m^3] * [m^2])
  if (density > 0.d0 .and. pipe_diameter > 0.d0) then
    velocity = pm_well%mass_flow_rate / &
               (density * PI * (0.5d0 * pipe_diameter)**2)
  else
    velocity = 0.d0
  endif

  !Compute Reynolds number
  reynolds = THWellComputeReynoldNumber(density,velocity,pipe_diameter, &
                                        viscosity)
  !Compute Prandtl number
  prandtl = THWellComputePrandtlNumber(heat_capacity,viscosity,thermal_cond)
  !Compute Nusselt number
  nusselt = THWellComputeNusseltNumber(reynolds,prandtl,pm_well)

  if (pipe_diameter > 0.d0) then
    THWellComputeConvCoefficient = nusselt * thermal_cond / pipe_diameter
  else
    THWellComputeConvCoefficient = 0.d0
  endif

end function THWellComputeConvCoefficient

! ************************************************************************** !

function THWellComputeReynoldNumber(density,velocity,pipe_diameter, &
                                    viscosity)
  !
  ! Computes Reynolds number for borehole flow
  ! Re = rho * u * D / mu
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/16/26
  !
  implicit none

  PetscReal :: density
  PetscReal :: velocity
  PetscReal :: pipe_diameter
  PetscReal :: viscosity
  PetscReal :: THWellComputeReynoldNumber

  if (viscosity > 0.d0) then
    THWellComputeReynoldNumber = density * dabs(velocity) * pipe_diameter / &
                                 viscosity
  else
    THWellComputeReynoldNumber = 0.d0
  endif

end function THWellComputeReynoldNumber

! ************************************************************************** !

function THWellComputePrandtlNumber(heat_capacity,viscosity,thermal_cond)
  !
  ! Computes Prandtl number for fluid
  ! Pr = c_p * mu / k
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/17/26
  !
  implicit none

  PetscReal :: heat_capacity
  PetscReal :: viscosity
  PetscReal :: thermal_cond
  PetscReal :: THWellComputePrandtlNumber

  if (thermal_cond > 0.d0) then
    THWellComputePrandtlNumber = heat_capacity * viscosity / thermal_cond
  else
    THWellComputePrandtlNumber = 0.d0
  endif

end function THWellComputePrandtlNumber

! ************************************************************************** !

function THWellComputeNusseltNumber(reynolds,prandtl,pm_well)
  !
  ! Computes Nusselt number using the configured correlation
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/18/26
  !
  implicit none

  PetscReal :: reynolds
  PetscReal :: prandtl
  class(pm_well_closed_loop_type) :: pm_well
  PetscReal :: THWellComputeNusseltNumber

  PetscReal :: xi
  PetscReal :: entrance
  PetscReal :: nu_lam
  PetscReal :: nu_turb
  PetscReal :: nu_at_1e4
  PetscReal :: gamma
  PetscReal :: zref
  PetscReal :: r_pi

  nu_lam = pm_well%nusselt_laminar

  select case(pm_well%nusselt_mode)
    case(TH_WELL_NUSSELT_LAMINAR)
      THWellComputeNusseltNumber = nu_lam
    case(TH_WELL_NUSSELT_CHEN_2021)
      if (reynolds <= 0.d0) then
        THWellComputeNusseltNumber = 0.d0
        return
      endif

      xi = (1.8d0 * log(max(reynolds,1.d0)) - 1.5d0)**(-2.d0)
      if (pm_well%use_nusselt_entrance_factor) then
        zref = max(pm_well%nusselt_z_ref,1.d-10)
        r_pi = 0.5d0 * pm_well%pipe_inner_diameter
        entrance = 1.d0 + (r_pi/zref)**(2.d0/3.d0)
      else
        entrance = 1.d0
      endif

      nu_turb = ((xi/8.d0) * reynolds * prandtl) / &
                (1.d0 + 12.7d0 * sqrt(xi/8.d0) * (prandtl**(2.d0/3.d0)-1.d0))
      nu_turb = nu_turb * entrance

      if (reynolds < 2300.d0) then
        THWellComputeNusseltNumber = nu_lam
      elseif (reynolds > 10000.d0) then
        THWellComputeNusseltNumber = nu_turb
      else
        gamma = (reynolds - 2300.d0) / (10000.d0 - 2300.d0)
        gamma = min(1.d0,max(0.d0,gamma))
        nu_at_1e4 = ((0.0308d0/8.d0) * 1.d4 * prandtl) / &
                    (1.d0 + 12.7d0 * sqrt(0.0308d0/8.d0) * &
                    (prandtl**(2.d0/3.d0)-1.d0))
        nu_at_1e4 = nu_at_1e4 * entrance
        THWellComputeNusseltNumber = (1.d0-gamma) * nu_lam + gamma * nu_at_1e4
      endif
    case(TH_WELL_NUSSELT_ZHANG_2015)
      ! Zhang et al. (2015) correlation remains the default.
      if (reynolds > 2300.d0) then
        ! turbulent flow,
        ! Zhang et al. 2015. "The analytical solution of the water-rock heat
        ! transfer coefficient and sensitivity analyses of parameters."
        ! Proceedings World Geothermal Congress 2015, Melbourne, Australia,
        ! 19-25, April 2015, pg 6
        THWellComputeNusseltNumber = 8.7d-5 * (reynolds**0.92d0) * &
                                     (prandtl**1.89d0)
      elseif (reynolds > 0.d0) then
        ! Laminar flow - use Nu = 4.36 for constant heat flux boundary condition
        ! (Nu = 3.66 for constant wall temperature)
        THWellComputeNusseltNumber = nu_lam
      else
        THWellComputeNusseltNumber = 0.d0
      endif
    case default
      ! Default to laminar value if unknown mode
      THWellComputeNusseltNumber = nu_lam
  end select

end function THWellComputeNusseltNumber

! ************************************************************************** !

subroutine THWellUpdateHeatFlux(pm_well_base)
  !
  ! Update the cumulative heat flux for the well model
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/25
  !
  use Option_module

  implicit none

  class(pm_well_type), pointer :: pm_well_base

  class(pm_well_closed_loop_type), pointer :: pm_well

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  pm_well%heat_flux(2,:) = pm_well%heat_flux(2,:) + &
                           pm_well%heat_flux(1,:) * pm_well%option%flow_dt

end subroutine THWellUpdateHeatFlux

end module TH_Well_module

