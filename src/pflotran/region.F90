module Region_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Geometry_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter, public :: DEFINED_BY_BLOCK = 1
  PetscInt, parameter, public :: DEFINED_BY_COORD = 2
  PetscInt, parameter, public :: DEFINED_BY_CELL_IDS = 3
  PetscInt, parameter, public :: DEFINED_BY_CELL_AND_FACE_IDS = 4
  PetscInt, parameter, public :: DEFINED_BY_VERTEX_IDS = 5
  PetscInt, parameter, public :: DEFINED_BY_SIDESET_UGRID = 6
  PetscInt, parameter, public :: DEFINED_BY_FACE_UGRID_EXP = 7
  PetscInt, parameter, public :: DEFINED_BY_POLY_BOUNDARY_FACE = 8
  PetscInt, parameter, public :: DEFINED_BY_POLY_CELL_CENTER = 9
  PetscInt, parameter, public :: DEFINED_BY_CARTESIAN_BOUNDARY = 10
  PetscInt, parameter, public :: DEFINED_BY_PLANAR_PATCH = 11
  PetscInt, parameter, public :: PLANAR_PATCH_ELLIPSE = 1
  PetscInt, parameter, public :: PLANAR_PATCH_RECTANGLE = 2

  type, public :: planar_patch_type
    PetscReal :: centroid(3)
    PetscReal :: normal(3)
    PetscReal :: axis1(3)
    PetscReal :: axis2(3)
    PetscReal :: radii(2)
    PetscReal :: half_thickness
    PetscInt :: shape
  end type planar_patch_type

  type, public :: block_type
    PetscInt :: i1,i2,j1,j2,k1,k2
    type(block_type), pointer :: next
  end type block_type

  type, public :: region_type
    PetscInt :: id
    PetscInt :: def_type
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXSTRINGLENGTH) :: filename
    PetscInt :: i1,i2,j1,j2,k1,k2
    type(point3d_type), pointer :: coordinates(:)
    PetscInt :: iface
    PetscInt :: num_cells
    PetscInt, pointer :: cell_ids(:)
    PetscInt, pointer :: faces(:)
    !TODO(geh): Tear anything to do with structured/unstructured grids other
    !           than cell id ane face id out of region.
    PetscInt, pointer :: vertex_ids(:,:) ! For Unstructured mesh
    PetscInt :: num_verts              ! For Unstructured mesh
    type(region_sideset_type), pointer :: sideset
    type(region_explicit_face_type), pointer :: explicit_faceset
    type(polygonal_volume_type), pointer :: polygonal_volume
    type(planar_patch_type), pointer :: planar_patch
    type(region_type), pointer :: next
  end type region_type

  type, public :: region_ptr_type
    type(region_type), pointer :: ptr
  end type region_ptr_type

  type, public :: region_list_type
    PetscInt :: num_regions
    type(region_type), pointer :: first
    type(region_type), pointer :: last
    type(region_type), pointer :: array(:)
  end type region_list_type

  type, public :: region_sideset_type
    PetscInt :: nfaces
    PetscInt, pointer :: face_vertices(:,:)
  end type region_sideset_type

  type, public :: region_explicit_face_type
    type(point3d_type), pointer :: face_centroids(:)
    PetscReal, pointer :: face_areas(:)
  end type region_explicit_face_type

  interface RegionCreate
    module procedure RegionCreateWithBlock
    module procedure RegionCreateWithList
    module procedure RegionCreateWithNothing
    module procedure RegionCreateWithRegion
  end interface RegionCreate

  interface RegionReadFromFile
    module procedure RegionReadFromFileId
    module procedure RegionReadFromFilename
    module procedure RegionReadSideSet
    module procedure RegionReadExplicitFaceSet
  end interface RegionReadFromFile

  public :: RegionCreate, &
            RegionRead, &
            RegionReadFromFile, &
            RegionInitList, &
            RegionAddToList, &
            RegionGetPtrFromList, &
            RegionDestroyList, &
            RegionReadSideSet, &
            RegionCreateSideset, &
            RegionCheckCellIndexBounds, &
            RegionInputRecord, &
            RegionPointInPlanarPatch, &
            RegionPatchCellMayHit, &
            RegionPatchHitsPolyhedron, &
            RegionPatchHitsSphere, &
            RegionDestroy

contains

! ************************************************************************** !

function RegionCreateWithNothing()
  !
  ! Creates a region with no arguments
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  !

  implicit none

  type(region_type), pointer :: RegionCreateWithNothing

  type(region_type), pointer :: region

  allocate(region)
  region%id = 0
  region%def_type = 0
  region%name = ""
  region%filename = ""
  region%i1 = 0
  region%i2 = 0
  region%j1 = 0
  region%j2 = 0
  region%k1 = 0
  region%k2 = 0
  region%iface = 0
  region%num_cells = 0
  ! By default it is assumed that the region is applicable to strucutred grid,
  ! unless explicitly stated in pflotran input file
  region%num_verts = 0
  nullify(region%coordinates)
  nullify(region%cell_ids)
  nullify(region%faces)
  nullify(region%vertex_ids)
  nullify(region%sideset)
  nullify(region%explicit_faceset)
  nullify(region%polygonal_volume)
  nullify(region%planar_patch)
  nullify(region%next)

  RegionCreateWithNothing => region

end function RegionCreateWithNothing

! ************************************************************************** !

function RegionCreateSideset()
  !
  ! Creates a sideset
  !
  ! Author: Glenn Hammond
  ! Date: 12/19/11
  !

  implicit none

  type(region_sideset_type), pointer :: RegionCreateSideset

  type(region_sideset_type), pointer :: sideset

  allocate(sideset)
  sideset%nfaces = 0
  nullify(sideset%face_vertices)

  RegionCreateSideset => sideset

end function RegionCreateSideset

! ************************************************************************** !

function RegionCreateExplicitFaceSet()
  !
  ! Creates a sideset
  !
  ! Author: Glenn Hammond
  ! Date: 12/19/11
  !

  implicit none

  type(region_explicit_face_type), pointer :: RegionCreateExplicitFaceSet

  type(region_explicit_face_type), pointer :: explicit_faceset

  allocate(explicit_faceset)
  nullify(explicit_faceset%face_centroids)
  nullify(explicit_faceset%face_areas)

  RegionCreateExplicitFaceSet => explicit_faceset

end function RegionCreateExplicitFaceSet

! ************************************************************************** !

function RegionCreateWithBlock(i1,i2,j1,j2,k1,k2)
  !
  ! Creates a region with i,j,k indices for arguments
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  !

  implicit none

  PetscInt :: i1, i2, j1, j2, k1, k2

  type(region_type), pointer :: RegionCreateWithBlock

  type(region_type), pointer :: region

  region => RegionCreateWithNothing()
  region%i1 = i1
  region%i2 = i2
  region%j1 = j2
  region%j2 = j2
  region%k1 = k1
  region%k2 = k2
  region%num_cells = (abs(i2-i1)+1)*(abs(j2-j1)+1)* &
                                    (abs(k2-k1)+1)

  RegionCreateWithBlock => region

end function RegionCreateWithBlock

! ************************************************************************** !

function RegionCreateWithList(list)
  !
  ! RegionCreate: Creates a region from a list of cells
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  !

  implicit none

  PetscInt :: list(:)

  type(region_type), pointer :: RegionCreateWithList

  type(region_type), pointer :: region

  region => RegionCreateWithNothing()
  region%num_cells = size(list)
  allocate(region%cell_ids(region%num_cells))
  region%cell_ids = list

  RegionCreateWithList => region

end function RegionCreateWithList

! ************************************************************************** !

function RegionCreateWithRegion(region)
  !
  ! Creates a copy of a region
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Grid_Unstructured_Cell_module

  implicit none

  type(region_type), pointer :: RegionCreateWithRegion
  type(region_type), pointer :: region

  type(region_type), pointer :: new_region
  PetscInt :: icount

  new_region => RegionCreateWithNothing()

  new_region%id = region%id
  new_region%def_type = region%def_type
  new_region%name = region%name
  new_region%filename = region%filename
  new_region%i1 = region%i1
  new_region%i2 = region%i2
  new_region%j1 = region%j1
  new_region%j2 = region%j2
  new_region%k1 = region%k1
  new_region%k2 = region%k2
  new_region%iface = region%iface
  new_region%num_cells = region%num_cells
  new_region%num_verts = region%num_verts
  if (associated(region%coordinates)) then
    call GeometryCopyCoordinates(region%coordinates, &
                                 new_region%coordinates)
  endif
  if (associated(region%cell_ids)) then
    allocate(new_region%cell_ids(new_region%num_cells))
    new_region%cell_ids(1:new_region%num_cells) = &
      region%cell_ids(1:region%num_cells)
  endif
  if (associated(region%faces)) then
    allocate(new_region%faces(new_region%num_cells))
    new_region%faces(1:new_region%num_cells) = &
      region%faces(1:region%num_cells)
  endif
  if (associated(region%vertex_ids)) then
    allocate(new_region%vertex_ids(0:MAX_VERT_PER_FACE,1:new_region%num_verts))
    new_region%vertex_ids(0:MAX_VERT_PER_FACE,1:new_region%num_verts) = &
    region%vertex_ids(0:MAX_VERT_PER_FACE,1:new_region%num_verts)
  endif
  if (associated(region%sideset)) then
    new_region%sideset => RegionCreateSideSet()
    new_region%sideset%nfaces = region%sideset%nfaces
    allocate(new_region%sideset%face_vertices( &
               size(region%sideset%face_vertices,1), &
               size(region%sideset%face_vertices,2)))
    new_region%sideset%face_vertices = region%sideset%face_vertices
  endif
  if (associated(region%explicit_faceset)) then
    new_region%explicit_faceset => RegionCreateExplicitFaceSet()
    allocate(new_region%explicit_faceset%face_centroids( &
               size(region%explicit_faceset%face_centroids)))
    new_region%explicit_faceset%face_centroids = &
      region%explicit_faceset%face_centroids
    do icount = 1, size(region%explicit_faceset%face_centroids)
      new_region%explicit_faceset%face_centroids(icount)%x = &
        region%explicit_faceset%face_centroids(icount)%x
      new_region%explicit_faceset%face_centroids(icount)%y = &
        region%explicit_faceset%face_centroids(icount)%y
      new_region%explicit_faceset%face_centroids(icount)%z = &
        region%explicit_faceset%face_centroids(icount)%z
    enddo
    allocate(new_region%explicit_faceset%face_areas( &
               size(region%explicit_faceset%face_areas)))
    new_region%explicit_faceset%face_areas = &
      region%explicit_faceset%face_areas
  endif
  if (associated(region%polygonal_volume)) then
    new_region%polygonal_volume => GeometryCreatePolygonalVolume()
    if (associated(region%polygonal_volume%xy_coordinates)) then
      call GeometryCopyCoordinates(region%polygonal_volume%xy_coordinates, &
                                   new_region%polygonal_volume%xy_coordinates)
    endif
    if (associated(region%polygonal_volume%xz_coordinates)) then
      call GeometryCopyCoordinates(region%polygonal_volume%xz_coordinates, &
                                   new_region%polygonal_volume%xz_coordinates)
    endif
    if (associated(region%polygonal_volume%yz_coordinates)) then
      call GeometryCopyCoordinates(region%polygonal_volume%yz_coordinates, &
                                   new_region%polygonal_volume%yz_coordinates)
    endif
  endif
  if (associated(region%planar_patch)) then
    call RegionCopyPlanarPatch(region%planar_patch, &
                               new_region%planar_patch)
  endif

  RegionCreateWithRegion => new_region

end function RegionCreateWithRegion

! ************************************************************************** !

subroutine RegionInitList(list)
  !
  ! Initializes a region list
  !
  ! Author: Glenn Hammond
  ! Date: 10/29/07
  !

  implicit none

  type(region_list_type) :: list

  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_regions = 0

end subroutine RegionInitList

! ************************************************************************** !

subroutine RegionAddToList(new_region,list)
  !
  ! Adds a new region to a region list
  !
  ! Author: Glenn Hammond
  ! Date: 10/29/07
  !

  implicit none

  type(region_type), pointer :: new_region
  type(region_list_type) :: list

  list%num_regions = list%num_regions + 1
  new_region%id = list%num_regions
  if (.not.associated(list%first)) list%first => new_region
  if (associated(list%last)) list%last%next => new_region
  list%last => new_region

end subroutine RegionAddToList

! ************************************************************************** !

subroutine RegionRead(region,input,option)
  !
  ! Reads a region from the input file
  !
  ! Author: Glenn Hammond
  ! Date: 02/20/08
  !
  use Input_Aux_module
  use String_module
  use Option_module
  use Grid_Structured_module

  implicit none

  type(option_type) :: option
  type(region_type) :: region
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: icount

  input%ierr = INPUT_ERROR_NONE
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword','REGION')
    call StringToUpper(keyword)

    select case(trim(keyword))

      case('BLOCK')
        region%def_type = DEFINED_BY_BLOCK
        call InputReadInt(input,option,region%i1)
        if (InputError(input)) then
          input%ierr = INPUT_ERROR_NONE
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'REGION')
          call InputReadInt(input,option,region%i1)
        endif
        call InputErrorMsg(input,option,'i1','REGION')
        call InputReadInt(input,option,region%i2)
        call InputErrorMsg(input,option,'i2','REGION')
        call InputReadInt(input,option,region%j1)
        call InputErrorMsg(input,option,'j1','REGION')
        call InputReadInt(input,option,region%j2)
        call InputErrorMsg(input,option,'j2','REGION')
        call InputReadInt(input,option,region%k1)
        call InputErrorMsg(input,option,'k1','REGION')
        call InputReadInt(input,option,region%k2)
        call InputErrorMsg(input,option,'k2','REGION')
      case('CARTESIAN_BOUNDARY')
        region%def_type = DEFINED_BY_CARTESIAN_BOUNDARY
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'cartesian boundary face','REGION')
        call StringToUpper(word)
        select case(word)
          case('WEST')
            region%iface = WEST_FACE
          case('EAST')
            region%iface = EAST_FACE
          case('NORTH')
            region%iface = NORTH_FACE
          case('SOUTH')
            region%iface = SOUTH_FACE
          case('BOTTOM')
            region%iface = BOTTOM_FACE
          case('TOP')
            region%iface = TOP_FACE
          case default
            option%io_buffer = 'Cartesian boundary face "' // trim(word) // &
              '" not recognized.'
            call PrintErrMsg(option)
        end select
      case('COORDINATE')
        region%def_type = DEFINED_BY_COORD
        allocate(region%coordinates(1))
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%x)
        if (InputError(input)) then
          input%ierr = INPUT_ERROR_NONE
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'REGION')
          call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%x)
        endif
        call InputErrorMsg(input,option,'x-coordinate','REGION')
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%y)
        call InputErrorMsg(input,option,'y-coordinate','REGION')
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%z)
        call InputErrorMsg(input,option,'z-coordinate','REGION')
      case('COORDINATES')
        region%def_type = DEFINED_BY_COORD
        call GeometryReadCoordinates(input,option,region%name, &
                                     region%coordinates)
      case('INFINITE')
        region%def_type = DEFINED_BY_COORD
        allocate(region%coordinates(2))
        region%coordinates(1)%x = -1.d20
        region%coordinates(1)%y = -1.d20
        region%coordinates(1)%z = -1.d20
        region%coordinates(2)%x = 1.d20
        region%coordinates(2)%y = 1.d20
        region%coordinates(2)%z = 1.d20
      case('POLYGON')
        if (.not.associated(region%polygonal_volume)) then
          region%polygonal_volume => GeometryCreatePolygonalVolume()
          region%def_type = DEFINED_BY_POLY_CELL_CENTER
        endif
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','REGION')
          call StringToUpper(word)
          select case(trim(word))
            case('TYPE')
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'polygon type','REGION')
              call StringToUpper(word)
              select case(word)
                case('BOUNDARY_FACES_IN_VOLUME')
                  region%def_type = DEFINED_BY_POLY_BOUNDARY_FACE
                case('CELL_CENTERS_IN_VOLUME')
                  region%def_type = DEFINED_BY_POLY_CELL_CENTER
                case default
                  option%io_buffer = 'REGION->POLYGON->"' // trim(word) // &
                    '" not recognized.'
                  call PrintErrMsg(option)
              end select
            case('XY')
              call GeometryReadCoordinates(input,option,region%name, &
                                         region%polygonal_volume%xy_coordinates)
            case('XZ')
              call GeometryReadCoordinates(input,option,region%name, &
                                         region%polygonal_volume%xz_coordinates)
            case('YZ')
              call GeometryReadCoordinates(input,option,region%name, &
                                         region%polygonal_volume%yz_coordinates)
            case default
              call InputKeywordUnrecognized(input,word, &
                                            'REGION POLYGON',option)
          end select
        enddo
        call InputPopBlock(input,option)
      case('PLANAR_PATCH')
        call RegionReadPlanarPatch(region,input,option)
      case('FILE')
        call InputReadFilename(input,option,region%filename)
        call InputErrorMsg(input,option,'filename','REGION')
      case('LIST')
        call InputReadPflotranString(input, option)
        string = trim(input%buf)
        icount = InputCountWordsInBuffer(input,option)
        input%buf = trim(string)
        select case(icount)
          case(1)
            call RegionReadCellList(region,input,PETSC_FALSE,PETSC_FALSE, &
                                    option)
          case(2)
            call RegionReadCellList(region,input,PETSC_TRUE,PETSC_FALSE, &
                                    option)
          case default
            option%io_buffer = 'REGION LIST format only supported for &
              &CELL_ID or CELL_ID FACE_ID format. One or two integers &
              &per line.'
            call PrintErrMsg(option)
        end select
      case('FACE')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'face','REGION')
        call StringToUpper(word)
        select case(word)
          case('WEST')
            region%iface = WEST_FACE
          case('EAST')
            region%iface = EAST_FACE
          case('NORTH')
            region%iface = NORTH_FACE
          case('SOUTH')
            region%iface = SOUTH_FACE
          case('BOTTOM')
            region%iface = BOTTOM_FACE
          case('TOP')
            region%iface = TOP_FACE
          case default
            option%io_buffer = 'FACE "' // trim(word) // &
              '" not recognized.'
            call PrintErrMsg(option)
        end select
      case default
        call InputKeywordUnrecognized(input,keyword,'REGION',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine RegionRead

! ************************************************************************** !

subroutine RegionReadPlanarPatch(region,input,option)
  !
  ! Reads a PLANAR_PATCH sub-block and sets up the orthonormal
  ! frame.
  !
  ! Author: Glenn Hammond
  ! Date: 09/11/26
  !

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(region_type) :: region
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  PetscReal :: centroid(3)
  PetscReal :: radii(2)
  PetscReal :: half_thickness
  PetscReal :: angles(2)
  PetscReal :: normal(3)
  PetscReal :: axis(3)
  PetscBool :: found_centroid
  PetscBool :: found_radii
  PetscBool :: found_thickness
  PetscBool :: found_angles
  PetscBool :: found_normal
  PetscBool :: found_axis
  PetscInt :: patch_shape

  region%def_type = DEFINED_BY_PLANAR_PATCH
  region%planar_patch => RegionCreatePlanarPatch()
  found_centroid = PETSC_FALSE
  found_radii = PETSC_FALSE
  found_thickness = PETSC_FALSE
  found_angles = PETSC_FALSE
  found_normal = PETSC_FALSE
  found_axis = PETSC_FALSE
  patch_shape = PLANAR_PATCH_ELLIPSE
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','REGION PLANAR_PATCH')
    call StringToUpper(word)
    select case(trim(word))
      case('CENTROID')
        call InputReadNDoubles(input,option,centroid,THREE_INTEGER)
        call InputErrorMsg(input,option,'CENTROID', &
                           'REGION PLANAR_PATCH')
        found_centroid = PETSC_TRUE
      case('ANGLES')
        call InputReadNDoubles(input,option,angles,TWO_INTEGER)
        call InputErrorMsg(input,option,'ANGLES', &
                           'REGION PLANAR_PATCH')
        found_angles = PETSC_TRUE
      case('NORMAL')
        call InputReadNDoubles(input,option,normal,THREE_INTEGER)
        call InputErrorMsg(input,option,'NORMAL', &
                           'REGION PLANAR_PATCH')
        found_normal = PETSC_TRUE
      case('AXIS')
        call InputReadNDoubles(input,option,axis,THREE_INTEGER)
        call InputErrorMsg(input,option,'AXIS', &
                           'REGION PLANAR_PATCH')
        found_axis = PETSC_TRUE
      case('RADII')
        call InputReadNDoubles(input,option,radii,TWO_INTEGER)
        call InputErrorMsg(input,option,'RADII', &
                           'REGION PLANAR_PATCH')
        found_radii = PETSC_TRUE
      case('HALF_THICKNESS')
        call InputReadDouble(input,option,half_thickness)
        call InputErrorMsg(input,option,'HALF_THICKNESS', &
                           'REGION PLANAR_PATCH')
        found_thickness = PETSC_TRUE
      case('SHAPE')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'SHAPE', &
                           'REGION PLANAR_PATCH')
        call StringToUpper(word)
        select case(trim(word))
          case('ELLIPSE')
            patch_shape = PLANAR_PATCH_ELLIPSE
          case('RECTANGLE')
            patch_shape = PLANAR_PATCH_RECTANGLE
          case default
            option%io_buffer = 'REGION PLANAR_PATCH SHAPE "' // &
              trim(word) // '" not recognized. Use ELLIPSE or &
              &RECTANGLE.'
            call PrintErrMsg(option)
        end select
      case default
        call InputKeywordUnrecognized(input,word, &
                                      'REGION PLANAR_PATCH',option)
    end select
  enddo
  call InputPopBlock(input,option)
  if (.not.found_centroid) then
    option%io_buffer = 'REGION PLANAR_PATCH requires CENTROID.'
    call PrintErrMsg(option)
  endif
  if (.not.found_radii) then
    option%io_buffer = 'REGION PLANAR_PATCH requires RADII.'
    call PrintErrMsg(option)
  endif
  if (.not.found_thickness) then
    option%io_buffer = 'REGION PLANAR_PATCH requires HALF_THICKNESS.'
    call PrintErrMsg(option)
  endif
  if (radii(1) <= 0.d0 .or. radii(2) <= 0.d0) then
    option%io_buffer = 'REGION PLANAR_PATCH RADII must be positive.'
    call PrintErrMsg(option)
  endif
  if (half_thickness <= 0.d0) then
    option%io_buffer = 'REGION PLANAR_PATCH HALF_THICKNESS must &
      &be greater than zero.'
    call PrintErrMsg(option)
  endif
  if (found_angles .and. found_normal) then
    option%io_buffer = 'REGION PLANAR_PATCH accepts ANGLES or &
      &NORMAL, not both.'
    call PrintErrMsg(option)
  endif
  if (.not.found_angles .and. .not.found_normal) then
    option%io_buffer = 'REGION PLANAR_PATCH requires ANGLES or &
      &NORMAL.'
    call PrintErrMsg(option)
  endif
  if (found_normal .and. .not.found_axis) then
    option%io_buffer = 'REGION PLANAR_PATCH NORMAL requires AXIS.'
    call PrintErrMsg(option)
  endif
  if (found_angles) then
    call RegionSetupPatchAngles(region%planar_patch,centroid, &
                                angles(1),angles(2),radii, &
                                half_thickness,patch_shape,option)
  else
    call RegionSetupPatchNormal(region%planar_patch,centroid, &
                                normal,axis,radii,half_thickness, &
                                patch_shape,option)
  endif

end subroutine RegionReadPlanarPatch

! ************************************************************************** !

subroutine RegionReadFromFilename(region,option,filename)
  !
  ! Reads a list of cells from a file named filename
  !
  ! Author: Glenn Hammond
  ! Date: 10/29/07
  !

  use Input_Aux_module
  use Option_module
  use Utility_module

  implicit none

  type(region_type) :: region
  type(option_type) :: option
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: filename

  input => InputCreate(IUNIT_TEMP,filename,option)
  call RegionReadFromFileId(region,input,option)
  call InputDestroy(input)

end subroutine RegionReadFromFilename

! ************************************************************************** !

subroutine RegionReadFromFileId(region,input,option)
  !
  ! Reads a list of cells from an open file
  !
  ! Author: Glenn Hammond
  ! Date: 10/29/07
  !

  use Input_Aux_module
  use Option_module
  use Utility_module
  use Logging_module
  use Grid_Unstructured_Cell_module

  implicit none

  type(region_type) :: region
  type(option_type) :: option
  type(input_type), pointer :: input

  character(len=1) :: backslash
  character(len=MAXSTRINGLENGTH) :: string

  PetscInt, pointer :: temp_int_array(:)
  PetscInt, pointer :: vert_id_0_p(:)
  PetscInt, pointer :: vert_id_1_p(:)
  PetscInt, pointer :: vert_id_2_p(:)
  PetscInt, pointer :: vert_id_3_p(:)
  PetscInt, pointer :: vert_id_4_p(:)
  PetscInt :: max_size
  PetscInt :: count
  PetscInt :: temp_int
  PetscInt :: ii
  PetscInt :: istart
  PetscInt :: iend
  PetscInt :: remainder
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_region_read_ascii,ierr);CHKERRQ(ierr)

  max_size = 1000
  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++

  allocate(temp_int_array(max_size))
  allocate(vert_id_0_p(max_size))
  allocate(vert_id_0_p(max_size))
  allocate(vert_id_1_p(max_size))
  allocate(vert_id_2_p(max_size))
  allocate(vert_id_3_p(max_size))
  allocate(vert_id_4_p(max_size))

  temp_int_array = 0
  vert_id_0_p = 0
  vert_id_1_p = -1
  vert_id_2_p = -1
  vert_id_3_p = -1
  vert_id_4_p = -1

  count = 0

  ! Determine if region definition in the input data is one of the following:
  !  1) Contains cell ids only : Only ONE entry per line
  !  2) Contains cell ids and face ids: TWO entries per line
  !  3) Contains vertex ids that make up the face: MORE than two entries per
  !     line
  call InputReadPflotranString(input, option)
  string = trim(input%buf)
  count = InputCountWordsInBuffer(input,option)
  input%buf = trim(string)

  if (count == 0) then
     option%io_buffer = 'ERROR while reading the region "' // &
          trim(region%name) // '" from file (zero entries in first row)'
     call PrintErrMsg(option)
  else if (count == 1) then
    ! Input data contains cell ids
    call RegionReadCellList(region,input,PETSC_FALSE,PETSC_TRUE,option)
  else if (count == 2) then
    ! Input data contains cell ids + face ids
    call RegionReadCellList(region,input,PETSC_TRUE,PETSC_TRUE,option)
  else
    option%io_buffer = 'The number of integers per line listed in &
      &region "' // trim(region%filename) // '" suggests that &
      &unstructured grid vertices are specified. Please use a .ss &
      &or .ex file for unstructured regions.'
!    call PrintErrMsg(option)
    ! Input data contains vertices
    vert_id_0_p(1) = temp_int_array(1)
    vert_id_1_p(1) = temp_int_array(2)
    vert_id_2_p(1) = temp_int_array(3)
    vert_id_3_p(1) = temp_int_array(4)
    if (vert_id_0_p(1) == 4 ) vert_id_4_p(1) = temp_int_array(5)
    count = 1 ! reset the counter to represent the num of rows read
    region%def_type = DEFINED_BY_VERTEX_IDS

    ! Read the data
    do
      ! InputReadPflotranString is at bottom since string has
      ! been read on first pass
      if (InputError(input)) exit
      call InputReadInt(input,option,temp_int)
      if (InputError(input)) exit
      count = count + 1
      vert_id_0_p(count) = temp_int

      vert_id_4_p(count) = UNINITIALIZED_INTEGER
      do ii = 1, vert_id_0_p(count)
        call InputReadInt(input,option,temp_int)
        if (InputError(input)) then
          option%io_buffer = 'ERROR while reading the region "' // &
            trim(region%name) // '" from file'
          call PrintErrMsg(option)
        endif

        select case(ii)
          case(1)
            vert_id_1_p(count) = temp_int
          case(2)
            vert_id_2_p(count) = temp_int
          case(3)
            vert_id_3_p(count) = temp_int
          case(4)
            vert_id_4_p(count) = temp_int
        end select

        if (count+1 > max_size) then ! resize temporary array
          call ReallocateArray(vert_id_0_p,max_size)
          ! since ReallocateArray doubles max_size, we need to divide by 2
          ! before calling again
          max_size = max_size / 2
          call ReallocateArray(vert_id_1_p,max_size)
          max_size = max_size / 2
          call ReallocateArray(vert_id_2_p,max_size)
          max_size = max_size / 2
          call ReallocateArray(vert_id_3_p,max_size)
          max_size = max_size / 2
          call ReallocateArray(vert_id_4_p,max_size)
        endif
      enddo
      call InputReadPflotranString(input,option)
    enddo

    ! Depending on processor rank, save only a portion of data
    region%num_verts = count/option%comm%size
      remainder = count - region%num_verts*option%comm%size
    if (option%myrank < remainder) region%num_verts = region%num_verts + 1
    istart = 0
    iend   = 0
    call MPI_Exscan(region%num_verts,istart,ONE_INTEGER_MPI,MPIU_INTEGER, &
                    MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
    call MPI_Scan(region%num_verts,iend,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                  option%mycomm,ierr);CHKERRQ(ierr)

    ! Allocate memory and save the data
    region%num_verts = iend - istart
    allocate(region%vertex_ids(0:MAX_VERT_PER_FACE,1:region%num_verts))
    region%vertex_ids(0,1:region%num_verts) = vert_id_0_p(istart + 1: iend)
    region%vertex_ids(1,1:region%num_verts) = vert_id_1_p(istart + 1: iend)
    region%vertex_ids(2,1:region%num_verts) = vert_id_2_p(istart + 1: iend)
    region%vertex_ids(3,1:region%num_verts) = vert_id_3_p(istart + 1: iend)
    region%vertex_ids(4,1:region%num_verts) = vert_id_4_p(istart + 1: iend)
    deallocate(vert_id_0_p)
    deallocate(vert_id_1_p)
    deallocate(vert_id_2_p)
    deallocate(vert_id_3_p)
    deallocate(vert_id_4_p)

  endif

#if 0
  count = 1
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    call InputReadInt(input,option,temp_int)
    if (.not.InputError(input)) then
      count = count + 1
      temp_int_array(count) = temp_int
      write(*,*) count,temp_int
    endif
    if (count+1 > max_size) then ! resize temporary array
      call ReallocateArray(temp_int_array,max_size)
    endif
  enddo

  if (count > 0) then
    region%num_cells = count
    allocate(region%cell_ids(count))
    region%cell_ids(1:count) = temp_int_array(1:count)
  else
    region%num_cells = 0
    nullify(region%cell_ids)
  endif
#endif
  deallocate(temp_int_array)

  call PetscLogEventEnd(logging%event_region_read_ascii,ierr);CHKERRQ(ierr)

end subroutine RegionReadFromFileId

! ************************************************************************** !

subroutine RegionReadSideSet(sideset,filename,option)
  !
  ! Reads an unstructured grid sideset
  !
  ! Author: Glenn Hammond
  ! Date: 12/19/11
  !

  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  type(region_sideset_type) :: sideset
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option

  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string, hint
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: num_faces_local_save
  PetscInt :: num_faces_local
  PetscInt :: num_to_read
  PetscInt, parameter :: max_nvert_per_face = 4
  PetscInt, allocatable :: temp_int_array(:,:)

  PetscInt :: iface, ivertex, irank, num_vertices
  PetscInt :: remainder
  PetscErrorCode :: ierr
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscMPIInt :: int_mpi
  PetscInt :: fileid

  fileid = 86
  input => InputCreate(fileid,filename,option)

! Format of sideset file
! type: T=triangle, Q=quadrilateral
! vertn(Q) = 4
! vertn(T) = 3
! -----------------------------------------------------------------
! num_faces  (integer)
! type vert1 vert2 ... vertn  ! for face 1 (integers)
! type vert1 vert2 ... vertn  ! for face 2
! ...
! ...
! type vert1 vert2 ... vertn  ! for face num_faces
! -----------------------------------------------------------------

  hint = 'Unstructured Sideset'

  call InputReadPflotranString(input,option)
  string = 'unstructured sideset'
  call InputReadStringErrorMsg(input,option,hint)

  ! read num_faces
  call InputReadInt(input,option,sideset%nfaces)
  call InputErrorMsg(input,option,'number of faces',hint)

  ! divide faces across ranks
  num_faces_local = sideset%nfaces/option%comm%size
  num_faces_local_save = num_faces_local
  remainder = sideset%nfaces - num_faces_local*option%comm%size
  if (option%myrank < remainder) num_faces_local = &
                                 num_faces_local + 1

  ! allocate array to store vertices for each faces
  allocate(sideset%face_vertices(max_nvert_per_face, &
                                 num_faces_local))
  sideset%face_vertices = UNINITIALIZED_INTEGER

  ! for now, read all faces from ASCII file through io_rank and communicate
  ! to other ranks
  call OptionSetBlocking(option,PETSC_FALSE)
  if (OptionIsIORank(option)) then
    allocate(temp_int_array(max_nvert_per_face, &
                            num_faces_local_save+1))
    ! read for other processors
    do irank = 0, option%comm%size-1
      temp_int_array = UNINITIALIZED_INTEGER
      num_to_read = num_faces_local_save
      if (irank < remainder) num_to_read = num_to_read + 1

      do iface = 1, num_to_read
        ! read in the vertices defining the cell face
        call InputReadPflotranString(input,option)
        call InputReadStringErrorMsg(input,option,hint)
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'face type',hint)
        call StringToUpper(word)
        select case(word)
          case('Q')
            num_vertices = 4
          case('T')
            num_vertices = 3
          case('L')
            num_vertices = 2
          case default
            option%io_buffer = 'Unknown face type "' // trim(word) // &
              '" in sideset file "' // trim(filename) // '". Please use &
              &"Q" (quadrilateral) or "T" (triangle).'
            call PrintErrMsgByRank(option)
        end select
        do ivertex = 1, num_vertices
          call InputReadInt(input,option,temp_int_array(ivertex,iface))
          call InputErrorMsg(input,option,'vertex id',hint)
        enddo
      enddo
      ! if the faces reside on io_rank
      if (OptionIsIORank(option,irank)) then
#if UGRID_DEBUG
        write(string,*) num_faces_local
        string = trim(adjustl(string)) // ' faces stored on p0'
        print *, trim(string)
#endif
        sideset%nfaces = num_faces_local
        sideset%face_vertices(:,1:num_faces_local) = &
          temp_int_array(:,1:num_faces_local)
      else
        ! otherwise communicate to other ranks
#if UGRID_DEBUG
        write(string,*) num_to_read
        write(word,*) irank
        string = trim(adjustl(string)) // ' faces sent from p0 to p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
        int_mpi = num_to_read*max_nvert_per_face
        call MPI_Send(temp_int_array,int_mpi,MPIU_INTEGER,irank,num_to_read, &
                      option%mycomm,ierr);CHKERRQ(ierr)
      endif
    enddo
    deallocate(temp_int_array)
  else
    ! other ranks post the recv
#if UGRID_DEBUG
        write(string,*) num_faces_local
        write(word,*) option%myrank
        string = trim(adjustl(string)) // ' faces received from p0 at p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
    sideset%nfaces = num_faces_local
    int_mpi = num_faces_local*max_nvert_per_face
    call MPI_Recv(sideset%face_vertices,int_mpi,MPIU_INTEGER, &
                  option%comm%io_rank,MPI_ANY_TAG,option%mycomm,status_mpi, &
                  ierr);CHKERRQ(ierr)
  endif
  call OptionSetBlocking(option,PETSC_TRUE)
  call OptionCheckNonBlockingError(option)

!  unstructured_grid%nlmax = num_faces_local
!  unstructured_grid%num_vertices_local = num_vertices_local

  call InputDestroy(input)

end subroutine RegionReadSideSet

! ************************************************************************** !

subroutine RegionReadExplicitFaceSet(explicit_faceset,cell_ids,filename,option)
  !
  ! Reads an unstructured grid explicit region
  !
  ! Author: Glenn Hammond
  ! Date: 05/18/12
  !
  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  type(region_explicit_face_type), pointer :: explicit_faceset
  PetscInt, pointer :: cell_ids(:)
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option

  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: hint
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: fileid

  PetscInt :: num_connections
  PetscInt :: iconn

  explicit_faceset => RegionCreateExplicitFaceSet()

  fileid = 86
  input => InputCreate(fileid,filename,option)

! Format of explicit unstructured grid file
! id_ = integer
! x_, y_, z_, area_ = real
! definitions
! id_ = id of grid cell
! x_ = x coordinate of cell face
! y_ = y coordinate of cell face
! z_ = z coordinate of cell face
! area_ = area of grid cell face
! -----------------------------------------------------------------
! CONNECTIONS <integer>   integer = # connections (M)
! id_1 x_1 y_1 z_1 area_1
! id_2 x_2 y_2 z_2 area_2
! ...
! ...
! id_M x_M y_M z_M area_M
! -----------------------------------------------------------------

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit

    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    hint = trim(word)

    select case(word)
      case('CONNECTIONS')
        hint = 'Explicit Unstructured Grid CONNECTIONS in file: ' // &
          trim(adjustl(filename))
        call InputReadInt(input,option,num_connections)
        call InputErrorMsg(input,option,'number of connections',hint)

        allocate(cell_ids(num_connections))
        cell_ids = 0
        allocate(explicit_faceset%face_areas(num_connections))
        explicit_faceset%face_areas = 0
        allocate(explicit_faceset%face_centroids(num_connections))
        do iconn = 1, num_connections
          explicit_faceset%face_centroids(iconn)%x = 0.d0
          explicit_faceset%face_centroids(iconn)%y = 0.d0
          explicit_faceset%face_centroids(iconn)%z = 0.d0
        enddo
        do iconn = 1, num_connections
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,hint)
          call InputReadInt(input,option,cell_ids(iconn))
          call InputErrorMsg(input,option,'cell id',hint)
          call InputReadDouble(input,option, &
                               explicit_faceset%face_centroids(iconn)%x)
          call InputErrorMsg(input,option,'face x coordinate',hint)
          call InputReadDouble(input,option, &
                               explicit_faceset%face_centroids(iconn)%y)
          call InputErrorMsg(input,option,'face y coordinate',hint)
          call InputReadDouble(input,option, &
                               explicit_faceset%face_centroids(iconn)%z)
          call InputErrorMsg(input,option,'face z coordinate',hint)
          call InputReadDouble(input,option, &
                               explicit_faceset%face_areas(iconn))
          call InputErrorMsg(input,option,'face area',hint)
        enddo
      case default
        call InputKeywordUnrecognized(input,word, &
               'REGION (explicit unstructured grid)',option)
    end select
  enddo
  call InputPopBlock(input,option)

  call InputDestroy(input)

end subroutine RegionReadExplicitFaceSet

! ************************************************************************** !

subroutine RegionReadCellList(region,input,read_faces,from_file,option)
  !
  ! Reads a list of cells (and optional faces) from an ASCII file
  !
  ! Author: Glenn Hammond
  ! Date: 05/19/23
  !
  use Input_Aux_module
  use Option_module
  use Utility_module

  implicit none

  type(region_type) :: region
  type(input_type), pointer :: input
  PetscBool :: read_faces
  PetscBool :: from_file ! true if read from a file (no block terminator)
  type(option_type) :: option

  PetscInt :: array_size
  PetscInt, pointer :: cell_ids(:)
  PetscInt, pointer :: face_ids(:)
  PetscInt :: temp_int
  PetscInt :: icount
  PetscInt :: istart
  PetscInt :: iend
  PetscInt :: remainder
  PetscErrorCode :: ierr

  if (read_faces) then
    region%def_type = DEFINED_BY_CELL_AND_FACE_IDS
  else
    region%def_type = DEFINED_BY_CELL_IDS
  endif

  array_size = 1000
  allocate(cell_ids(array_size))
  cell_ids(:) = 0
  if (read_faces) then
    allocate(face_ids(array_size))
    face_ids(:) = 0
  else
    nullify(face_ids)
  endif

  ! Read the data
  icount = 0
  input%ierr = INPUT_ERROR_NONE ! first pass must be success
  do
    ! InputReadPflotranString is at bottom since string has
    ! been read on first pass
    if (from_file) then
      if (InputError(input)) exit
    else
      if (InputCheckExit(input,option)) exit
    endif
    call InputReadInt(input, option, temp_int)
    if (InputError(input)) then
      if (from_file) then
        exit
      else
        option%io_buffer = 'ERROR reading cell ID in REGION "' // &
          trim(region%name) // '".'
        call PrintErrMsg(option)
      endif
    endif
    icount = icount + 1
    cell_ids(icount) = temp_int

    if (read_faces) then
      call InputReadInt(input,option,temp_int)
      if (InputError(input)) then
        option%io_buffer = 'ERROR reading face ID in REGION "' // &
          trim(region%name) // '".'
        call PrintErrMsg(option)
      endif
      face_ids(icount) = temp_int
    endif
    if (icount+1 > array_size) then ! resize temporary array
      call ReallocateArray(cell_ids,array_size)
      ! since ReallocateArray doubles max_size, we need to divide by 2
      ! before calling again
      if (read_faces) then
        array_size = array_size / 2
        call ReallocateArray(face_ids,array_size)
      endif
    endif
    call InputReadPflotranString(input, option)
  enddo

  ! Depending on processor rank, save only a portion of data
  region%num_cells = icount/option%comm%size
  remainder = icount - region%num_cells*option%comm%size
  if (option%myrank < remainder) region%num_cells = region%num_cells + 1
  istart = 0
  iend   = 0
  call MPI_Exscan(region%num_cells,istart,ONE_INTEGER_MPI,MPIU_INTEGER, &
                  MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)
  call MPI_Scan(region%num_cells,iend,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                option%mycomm,ierr);CHKERRQ(ierr)

  ! Allocate memory and save the data
  allocate(region%cell_ids(region%num_cells))
  region%cell_ids(1:region%num_cells) = cell_ids(istart+1:iend)
  if (read_faces) then
    allocate(region%faces(region%num_cells))
    region%faces(1:region%num_cells) = face_ids(istart+1:iend)
  endif
  call DeallocateArray(cell_ids)
  call DeallocateArray(face_ids)

end subroutine RegionReadCellList

! ************************************************************************** !

function RegionGetPtrFromList(region_name,region_list)
  !
  ! Returns a pointer to the region matching region_name
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  !

  use String_module

  implicit none

  type(region_type), pointer :: RegionGetPtrFromList
  character(len=MAXWORDLENGTH) :: region_name
  PetscInt :: length
  type(region_list_type) :: region_list

  type(region_type), pointer :: region

  nullify(RegionGetPtrFromList)
  region => region_list%first

  do
    if (.not.associated(region)) exit
    length = len_trim(region_name)
    if (length == len_trim(region%name) .and. &
        StringCompare(region%name,region_name,length)) then
      RegionGetPtrFromList => region
      return
    endif
    region => region%next
  enddo

end function RegionGetPtrFromList

! ************************************************************************** !

subroutine RegionCheckCellIndexBounds(region,num_cells,option)
  !
  ! Checks to ensure that cell ids listed in a region are within the bounds
  ! of 1 and the maximum cell id.
  !
  ! Author: Glenn Hammond
  ! Date: 06/27/19
  !
  use Option_module
  use String_module

  implicit none

  type(region_type) :: region
  PetscInt :: num_cells
  type(option_type) :: option

  PetscInt :: cell_id_extremes(2)
  PetscErrorCode :: ierr

  cell_id_extremes(1) = 999999999
  cell_id_extremes(2) = -cell_id_extremes(1)
  if (region%num_cells > 0) then
    cell_id_extremes(1) = minval(region%cell_ids)
    cell_id_extremes(2) = maxval(region%cell_ids)
  endif

  ! invert for MPI max below
  cell_id_extremes(1) = -cell_id_extremes(1)
  call MPI_Allreduce(MPI_IN_PLACE,cell_id_extremes,TWO_INTEGER_MPI, &
                     MPI_INTEGER,MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
  ! invert back
  cell_id_extremes(1) = -cell_id_extremes(1)

  if (cell_id_extremes(1) < 1 .or. cell_id_extremes(2) > num_cells) then
    option%io_buffer = 'The minimum cell ID (' // &
      trim(StringWrite(cell_id_extremes(1))) // &
      ') and/or maximum cell ID (' // &
      trim(StringWrite(cell_id_extremes(2))) // &
      ') for REGION "' // trim(region%name) // &
      '" is outside the GRID cell ID bounds of (1 - ' // &
      trim(StringWrite(num_cells)) // ').'
    call PrintErrMsg(option)
  endif

end subroutine RegionCheckCellIndexBounds

! **************************************************************************** !

subroutine RegionInputRecord(region_list)
  !
  ! Prints ingested region information to the input record file
  !
  ! Author: Jenn Frederick
  ! Date: 03/30/2016
  !
  use Grid_Structured_module

  implicit none

  type(region_list_type), pointer :: region_list

  type(region_type), pointer :: cur_region
  character(len=MAXWORDLENGTH) :: word1, word2
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: k
  PetscInt :: id = INPUT_RECORD_UNIT
  character(len=10) :: sFormat, iFormat

  sFormat = '(ES14.7)'
  iFormat = '(I10)'

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'REGIONS'

  cur_region => region_list%first
  do
    if (.not.associated(cur_region)) exit
    write(id,'(a29)',advance='no') 'region: '
    write(id,'(a)') adjustl(trim(cur_region%name))
    if (len_trim(cur_region%filename) > 0) then
      write(id,'(a29)',advance='no') 'from file: '
      write(id,'(a)') adjustl(trim(cur_region%filename))
    endif

    select case (cur_region%def_type)
    !--------------------------------
      case (DEFINED_BY_BLOCK)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'BLOCK'
        write(id,'(a29)',advance='no') 'I indices: '
        write(word1,iFormat) cur_region%i1
        write(word2,iFormat) cur_region%i2
        write(id,'(a)') adjustl(trim(word1)) // ' ' // adjustl(trim(word2))
        write(id,'(a29)',advance='no') 'J indices: '
        write(word1,iFormat) cur_region%j1
        write(word2,iFormat) cur_region%j2
        write(id,'(a)') adjustl(trim(word1)) // ' ' // adjustl(trim(word2))
        write(id,'(a29)',advance='no') 'K indices: '
        write(word1,iFormat) cur_region%k1
        write(word2,iFormat) cur_region%k2
        write(id,'(a)') adjustl(trim(word1)) // ' ' // adjustl(trim(word2))
    !--------------------------------
      case (DEFINED_BY_CARTESIAN_BOUNDARY)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'CARTESIAN BOUNDARY'
    !--------------------------------
      case (DEFINED_BY_COORD)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'COORDINATE(S)'
        write(id,'(a29)',advance='no') 'X coordinate(s): '
        string = ''
        do k = 1,size(cur_region%coordinates)
         write(word1,sFormat) cur_region%coordinates(k)%x
         string = adjustl(trim(string)) // ' ' // adjustl(trim(word1))
        enddo
        write(id,'(a)') adjustl(trim(string)) // ' m'
        write(id,'(a29)',advance='no') 'Y coordinate(s): '
        string = ''
        do k = 1,size(cur_region%coordinates)
          write(word1,sFormat) cur_region%coordinates(k)%y
          string = adjustl(trim(string)) // ' ' // adjustl(trim(word1))
        enddo
        write(id,'(a)') adjustl(trim(string)) // ' m'
        write(id,'(a29)',advance='no') 'Z coordinate(s): '
        string = ''
        do k = 1,size(cur_region%coordinates)
          write(word1,sFormat) cur_region%coordinates(k)%z
          string = adjustl(trim(string)) // ' ' // adjustl(trim(word1))
        enddo
        write(id,'(a)') adjustl(trim(string)) // ' m'
    !--------------------------------
      case (DEFINED_BY_CELL_AND_FACE_IDS)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'CELL AND FACE IDS'
    !--------------------------------
      case (DEFINED_BY_CELL_IDS)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'CELL IDS'
    !--------------------------------
      case (DEFINED_BY_VERTEX_IDS)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'VERTEX IDS'
    !--------------------------------
      case (DEFINED_BY_FACE_UGRID_EXP)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'FACE UNSTRUCTURED GRID EXPLICIT'
    !--------------------------------
      case (DEFINED_BY_POLY_BOUNDARY_FACE)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'POLYGON BOUNDARY FACES IN VOLUME'
    !--------------------------------
      case (DEFINED_BY_POLY_CELL_CENTER)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'POLYGON CELL CENTERS IN VOLUME'
    !--------------------------------
      case (DEFINED_BY_PLANAR_PATCH)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'PLANAR PATCH'
    !--------------------------------
    end select

    if (cur_region%iface /= 0) then
      write(id,'(a29)',advance='no') 'face: '
      select case (cur_region%iface)
        case (WEST_FACE)
          write(id,'(a)') 'west'
        case (EAST_FACE)
          write(id,'(a)') 'east'
        case (NORTH_FACE)
          write(id,'(a)') 'north'
        case (SOUTH_FACE)
          write(id,'(a)') 'south'
        case (BOTTOM_FACE)
          write(id,'(a)') 'bottom'
        case (TOP_FACE)
          write(id,'(a)') 'top'
      end select
    endif

    write(id,'(a29)') '---------------------------: '
    cur_region => cur_region%next
  enddo

end subroutine RegionInputRecord

! **************************************************************************** !

subroutine RegionDestroySideset(sideset)
  !
  ! Deallocates a unstructured grid side set
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/09
  !

  implicit none

  type(region_sideset_type), pointer :: sideset

  if (.not.associated(sideset)) return

  if (associated(sideset%face_vertices)) deallocate(sideset%face_vertices)
  nullify(sideset%face_vertices)

  deallocate(sideset)
  nullify(sideset)

end subroutine RegionDestroySideset

! ************************************************************************** !

subroutine RegionDestroyExplicitFaceSet(explicit_faceset)
  !
  ! Deallocates a unstructured grid explicit grid
  !
  ! Author: Glenn Hammond
  ! Date: 05/18/12
  !

  use Utility_module, only : DeallocateArray

  implicit none

  type(region_explicit_face_type), pointer :: explicit_faceset

  if (.not.associated(explicit_faceset)) return

  if (associated(explicit_faceset%face_centroids)) &
    deallocate(explicit_faceset%face_centroids)
  nullify(explicit_faceset%face_centroids)
  call DeallocateArray(explicit_faceset%face_areas)

  deallocate(explicit_faceset)
  nullify(explicit_faceset)

end subroutine RegionDestroyExplicitFaceSet

! ************************************************************************** !

subroutine RegionDestroyList(region_list)
  !
  ! Deallocates a list of regions
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  !

  implicit none

  type(region_list_type), pointer :: region_list

  type(region_type), pointer :: region, prev_region

  if (.not.associated(region_list)) return

  region => region_list%first
  do
    if (.not.associated(region)) exit
    prev_region => region
    region => region%next
    call RegionDestroy(prev_region)
  enddo

  region_list%num_regions = 0
  nullify(region_list%first)
  nullify(region_list%last)
  if (associated(region_list%array)) deallocate(region_list%array)
  nullify(region_list%array)

  deallocate(region_list)
  nullify(region_list)

end subroutine RegionDestroyList

! ************************************************************************** !

function RegionCreatePlanarPatch()
  !
  ! Creates a planar patch (ellipse or rectangle) with finite
  ! thickness along the plane normal.
  !
  ! Author: Glenn Hammond
  ! Date: 09/11/26
  !

  implicit none

  type(planar_patch_type), pointer :: RegionCreatePlanarPatch

  type(planar_patch_type), pointer :: patch

  allocate(patch)
  patch%centroid = UNINITIALIZED_DOUBLE
  patch%normal = UNINITIALIZED_DOUBLE
  patch%axis1 = UNINITIALIZED_DOUBLE
  patch%axis2 = UNINITIALIZED_DOUBLE
  patch%radii = UNINITIALIZED_DOUBLE
  patch%half_thickness = UNINITIALIZED_DOUBLE
  patch%shape = PLANAR_PATCH_ELLIPSE

  RegionCreatePlanarPatch => patch

end function RegionCreatePlanarPatch

! ************************************************************************** !

subroutine RegionSetupPatchAngles(patch,centroid,angle_xy,angle_xz, &
                                  radii,half_thickness,shape,option)
  !
  ! Builds the orthonormal frame from XY- and XZ-trace angles
  ! (degrees from +X).
  !
  ! Author: Glenn Hammond
  ! Date: 09/11/26
  !

  use Option_module

  implicit none

  type(planar_patch_type) :: patch
  PetscReal :: centroid(3)
  PetscReal :: angle_xy
  PetscReal :: angle_xz
  PetscReal :: radii(2)
  PetscReal :: half_thickness
  PetscInt :: shape
  type(option_type) :: option

  PetscReal :: u(3)
  PetscReal :: v(3)
  PetscReal :: n(3)
  PetscReal :: e1(3)
  PetscReal :: e2(3)
  PetscReal :: nrm
  PetscReal :: deg_to_rad

  deg_to_rad = PI / 180.d0
  u(1) = cos(angle_xy*deg_to_rad)
  u(2) = sin(angle_xy*deg_to_rad)
  u(3) = 0.d0
  v(1) = cos(angle_xz*deg_to_rad)
  v(2) = 0.d0
  v(3) = sin(angle_xz*deg_to_rad)
  call PatchCrossProduct(u,v,n)
  nrm = sqrt(n(1)**2+n(2)**2+n(3)**2)
  if (nrm < 1.d-12) then
    option%io_buffer = 'PLANAR_PATCH ANGLES produce parallel traces.'
    call PrintErrMsg(option)
  endif
  n = n / nrm
  nrm = sqrt(u(1)**2+u(2)**2+u(3)**2)
  e1 = u / nrm
  call PatchCrossProduct(n,e1,e2)
  nrm = sqrt(e2(1)**2+e2(2)**2+e2(3)**2)
  e2 = e2 / nrm
  patch%centroid = centroid
  patch%normal = n
  patch%axis1 = e1
  patch%axis2 = e2
  patch%radii = radii
  patch%half_thickness = half_thickness
  patch%shape = shape

end subroutine RegionSetupPatchAngles

! ************************************************************************** !

subroutine RegionSetupPatchNormal(patch,centroid,normal,axis, &
                                  radii,half_thickness,shape,option)
  !
  ! Builds the orthonormal frame from a plane normal and an
  ! in-plane axis (projected onto the plane).
  !
  ! Author: Glenn Hammond
  ! Date: 09/11/26
  !

  use Option_module

  implicit none

  type(planar_patch_type) :: patch
  PetscReal :: centroid(3)
  PetscReal :: normal(3)
  PetscReal :: axis(3)
  PetscReal :: radii(2)
  PetscReal :: half_thickness
  PetscInt :: shape
  type(option_type) :: option

  PetscReal :: n(3)
  PetscReal :: e1(3)
  PetscReal :: e2(3)
  PetscReal :: nrm
  PetscReal :: n_dot_a

  nrm = sqrt(normal(1)**2+normal(2)**2+normal(3)**2)
  if (nrm < 1.d-12) then
    option%io_buffer = 'PLANAR_PATCH NORMAL has zero length.'
    call PrintErrMsg(option)
  endif
  n = normal / nrm
  n_dot_a = axis(1)*n(1) + axis(2)*n(2) + axis(3)*n(3)
  e1(1) = axis(1) - n_dot_a*n(1)
  e1(2) = axis(2) - n_dot_a*n(2)
  e1(3) = axis(3) - n_dot_a*n(3)
  nrm = sqrt(e1(1)**2+e1(2)**2+e1(3)**2)
  if (nrm < 1.d-12) then
    option%io_buffer = 'PLANAR_PATCH AXIS is parallel to NORMAL.'
    call PrintErrMsg(option)
  endif
  e1 = e1 / nrm
  call PatchCrossProduct(n,e1,e2)
  nrm = sqrt(e2(1)**2+e2(2)**2+e2(3)**2)
  e2 = e2 / nrm
  patch%centroid = centroid
  patch%normal = n
  patch%axis1 = e1
  patch%axis2 = e2
  patch%radii = radii
  patch%half_thickness = half_thickness
  patch%shape = shape

end subroutine RegionSetupPatchNormal

! ************************************************************************** !

function RegionPointInPlanarPatch(x,y,z,patch)
  !
  ! True if (x,y,z) lies in the finite-thickness planar patch.
  !
  ! Author: Glenn Hammond
  ! Date: 09/11/26
  !

  implicit none

  PetscReal :: x
  PetscReal :: y
  PetscReal :: z
  type(planar_patch_type) :: patch

  PetscBool :: RegionPointInPlanarPatch

  PetscReal :: dx
  PetscReal :: dy
  PetscReal :: dz
  PetscReal :: dist_n
  PetscReal :: xi
  PetscReal :: eta

  RegionPointInPlanarPatch = PETSC_FALSE
  dx = x - patch%centroid(1)
  dy = y - patch%centroid(2)
  dz = z - patch%centroid(3)
  dist_n = dx*patch%normal(1) + dy*patch%normal(2) + dz*patch%normal(3)
  if (abs(dist_n) > patch%half_thickness) return
  xi = dx*patch%axis1(1) + dy*patch%axis1(2) + dz*patch%axis1(3)
  eta = dx*patch%axis2(1) + dy*patch%axis2(2) + dz*patch%axis2(3)
  select case(patch%shape)
    case(PLANAR_PATCH_ELLIPSE)
      if ((xi/patch%radii(1))**2 + (eta/patch%radii(2))**2 <= 1.d0) then
        RegionPointInPlanarPatch = PETSC_TRUE
      endif
    case(PLANAR_PATCH_RECTANGLE)
      if (abs(xi) <= patch%radii(1) .and. abs(eta) <= patch%radii(2)) then
        RegionPointInPlanarPatch = PETSC_TRUE
      endif
  end select

end function RegionPointInPlanarPatch

! ************************************************************************** !

function RegionPatchCellMayHit(patch,cx,cy,cz,radius)
  !
  ! Conservative circumsphere test: false means the cell cannot
  ! meet the patch. True means it might; run the exact test.
  !
  ! Author: Glenn Hammond
  ! Date: 09/11/26
  !

  implicit none

  type(planar_patch_type) :: patch
  PetscReal :: cx
  PetscReal :: cy
  PetscReal :: cz
  PetscReal :: radius

  PetscBool :: RegionPatchCellMayHit

  PetscReal :: dx
  PetscReal :: dy
  PetscReal :: dz
  PetscReal :: dist_n
  PetscReal :: xi
  PetscReal :: eta

  RegionPatchCellMayHit = PETSC_FALSE
  if (radius < 0.d0) return
  dx = cx - patch%centroid(1)
  dy = cy - patch%centroid(2)
  dz = cz - patch%centroid(3)
  dist_n = dx*patch%normal(1) + dy*patch%normal(2) + dz*patch%normal(3)
  if (abs(dist_n) > radius + patch%half_thickness) return
  xi = dx*patch%axis1(1) + dy*patch%axis1(2) + dz*patch%axis1(3)
  eta = dx*patch%axis2(1) + dy*patch%axis2(2) + dz*patch%axis2(3)
  select case(patch%shape)
    case(PLANAR_PATCH_ELLIPSE)
      RegionPatchCellMayHit = PatchCircleHitsEllipse(xi,eta,radius, &
                                      patch%radii(1),patch%radii(2))
    case(PLANAR_PATCH_RECTANGLE)
      RegionPatchCellMayHit = PatchCircleHitsRectangle(xi,eta,radius, &
                                      patch%radii(1),patch%radii(2))
  end select

end function RegionPatchCellMayHit

! ************************************************************************** !

function RegionPatchHitsPolyhedron(patch,x,y,z,nvert, &
                                   dist,qx,qy,qz,xi,eta,lwork)
  !
  ! True if a convex cell (vertex list) intersects the planar patch.
  ! dist(nvert) and qx,qy,qz,xi,eta(lwork) are caller scratch.
  !
  ! Author: Glenn Hammond
  ! Date: 09/11/26
  !

  implicit none

  type(planar_patch_type) :: patch
  PetscInt :: nvert
  PetscReal :: x(nvert)
  PetscReal :: y(nvert)
  PetscReal :: z(nvert)
  PetscReal :: dist(nvert)
  PetscInt :: lwork
  PetscReal :: qx(lwork)
  PetscReal :: qy(lwork)
  PetscReal :: qz(lwork)
  PetscReal :: xi(lwork)
  PetscReal :: eta(lwork)

  PetscBool :: RegionPatchHitsPolyhedron

  PetscInt :: i
  PetscInt :: j
  PetscInt :: nq
  PetscReal :: di
  PetscReal :: dj
  PetscReal :: t
  PetscReal :: h
  PetscReal :: px
  PetscReal :: py
  PetscReal :: pz

  RegionPatchHitsPolyhedron = PETSC_FALSE
  if (nvert < 1) return

  h = patch%half_thickness
  do i = 1, nvert
    dist(i) = (x(i)-patch%centroid(1))*patch%normal(1) + &
              (y(i)-patch%centroid(2))*patch%normal(2) + &
              (z(i)-patch%centroid(3))*patch%normal(3)
  enddo
  if (maxval(dist(1:nvert)) < -h .or. minval(dist(1:nvert)) > h) return

  nq = 0
  do i = 1, nvert
    if (abs(dist(i)) <= h) then
      nq = nq + 1
      qx(nq) = x(i)
      qy(nq) = y(i)
      qz(nq) = z(i)
    endif
  enddo
  do i = 1, nvert-1
    do j = i+1, nvert
      di = dist(i)
      dj = dist(j)
      if ((di-h)*(dj-h) < 0.d0 .and. nq < lwork) then
        t = (h-di)/(dj-di)
        nq = nq + 1
        qx(nq) = x(i) + t*(x(j)-x(i))
        qy(nq) = y(i) + t*(y(j)-y(i))
        qz(nq) = z(i) + t*(z(j)-z(i))
      endif
      if ((di+h)*(dj+h) < 0.d0 .and. nq < lwork) then
        t = (-h-di)/(dj-di)
        nq = nq + 1
        qx(nq) = x(i) + t*(x(j)-x(i))
        qy(nq) = y(i) + t*(y(j)-y(i))
        qz(nq) = z(i) + t*(z(j)-z(i))
      endif
    enddo
  enddo
  if (nq < 1) return

  do i = 1, nq
    px = qx(i) - patch%centroid(1)
    py = qy(i) - patch%centroid(2)
    pz = qz(i) - patch%centroid(3)
    xi(i) = px*patch%axis1(1) + py*patch%axis1(2) + pz*patch%axis1(3)
    eta(i) = px*patch%axis2(1) + py*patch%axis2(2) + pz*patch%axis2(3)
  enddo
  RegionPatchHitsPolyhedron = PatchPolygonHitsPatch(xi,eta,nq,patch)

end function RegionPatchHitsPolyhedron

! ************************************************************************** !

function RegionPatchHitsSphere(patch,cx,cy,cz,volume)
  !
  ! True if a sphere of equal volume at (cx,cy,cz) meets the patch:
  ! the sphere is projected onto the patch plane and the resulting
  ! circle is tested against the ellipse or rectangle.
  !
  ! Author: Glenn Hammond
  ! Date: 09/11/26
  !

  implicit none

  type(planar_patch_type) :: patch
  PetscReal :: cx
  PetscReal :: cy
  PetscReal :: cz
  PetscReal :: volume

  PetscBool :: RegionPatchHitsSphere

  PetscReal :: radius

  RegionPatchHitsSphere = PETSC_FALSE
  if (volume <= 0.d0) return
  radius = (0.75d0*volume/PI)**(1.d0/3.d0)
  RegionPatchHitsSphere = RegionPatchCellMayHit(patch,cx,cy,cz,radius)

end function RegionPatchHitsSphere

! ************************************************************************** !

function PatchPolygonHitsPatch(xi,eta,n,patch)
  !
  ! True if a convex 2D polygon overlaps the in-plane ellipse or
  ! rectangle.
  !

  implicit none

  PetscInt :: n
  PetscReal :: xi(n)
  PetscReal :: eta(n)
  type(planar_patch_type) :: patch

  PetscBool :: PatchPolygonHitsPatch

  PetscInt :: i
  PetscReal :: a
  PetscReal :: b

  PatchPolygonHitsPatch = PETSC_FALSE
  if (n < 1) return
  a = patch%radii(1)
  b = patch%radii(2)
  do i = 1, n
    select case(patch%shape)
      case(PLANAR_PATCH_ELLIPSE)
        if ((xi(i)/a)**2 + (eta(i)/b)**2 <= 1.d0) then
          PatchPolygonHitsPatch = PETSC_TRUE
          return
        endif
      case(PLANAR_PATCH_RECTANGLE)
        if (abs(xi(i)) <= a .and. abs(eta(i)) <= b) then
          PatchPolygonHitsPatch = PETSC_TRUE
          return
        endif
    end select
  enddo
  if (n >= 3) then
    if (PatchOriginInPolygon(xi,eta,n)) then
      PatchPolygonHitsPatch = PETSC_TRUE
      return
    endif
  endif
  do i = 1, n
    select case(patch%shape)
      case(PLANAR_PATCH_ELLIPSE)
        if (PatchSegmentHitsEllipse(xi(i),eta(i), &
              xi(mod(i,n)+1),eta(mod(i,n)+1),a,b)) then
          PatchPolygonHitsPatch = PETSC_TRUE
          return
        endif
      case(PLANAR_PATCH_RECTANGLE)
        if (PatchSegmentHitsRectangle(xi(i),eta(i), &
              xi(mod(i,n)+1),eta(mod(i,n)+1),a,b)) then
          PatchPolygonHitsPatch = PETSC_TRUE
          return
        endif
    end select
  enddo

end function PatchPolygonHitsPatch

! ************************************************************************** !

function PatchOriginInPolygon(xi,eta,n)
  !
  ! Even-odd test of whether (0,0) is inside a 2D polygon.
  !

  implicit none

  PetscInt :: n
  PetscReal :: xi(n)
  PetscReal :: eta(n)

  PetscBool :: PatchOriginInPolygon

  PetscInt :: i
  PetscInt :: j
  PetscReal :: xhit

  PatchOriginInPolygon = PETSC_FALSE
  j = n
  do i = 1, n
    if ((eta(i) < 0.d0 .and. eta(j) >= 0.d0) .or. &
        (eta(j) < 0.d0 .and. eta(i) >= 0.d0)) then
      xhit = xi(i) + (0.d0-eta(i))/(eta(j)-eta(i))*(xi(j)-xi(i))
      if (xhit < 0.d0) then
        PatchOriginInPolygon = .not.PatchOriginInPolygon
      endif
    endif
    j = i
  enddo

end function PatchOriginInPolygon

! ************************************************************************** !

function PatchSegmentHitsEllipse(x1,y1,x2,y2,a,b)
  !
  ! True if segment (x1,y1)-(x2,y2) intersects the filled ellipse.
  !

  implicit none

  PetscReal :: x1, y1, x2, y2
  PetscReal :: a, b

  PetscBool :: PatchSegmentHitsEllipse

  PetscReal :: dx
  PetscReal :: dy
  PetscReal :: aa
  PetscReal :: bb
  PetscReal :: cc
  PetscReal :: disc
  PetscReal :: t
  PetscReal :: sqrt_disc

  PatchSegmentHitsEllipse = PETSC_FALSE
  if ((x1/a)**2+(y1/b)**2 <= 1.d0) then
    PatchSegmentHitsEllipse = PETSC_TRUE
    return
  endif
  if ((x2/a)**2+(y2/b)**2 <= 1.d0) then
    PatchSegmentHitsEllipse = PETSC_TRUE
    return
  endif
  dx = x2 - x1
  dy = y2 - y1
  aa = dx*dx/(a*a) + dy*dy/(b*b)
  bb = 2.d0*(x1*dx/(a*a) + y1*dy/(b*b))
  cc = x1*x1/(a*a) + y1*y1/(b*b) - 1.d0
  if (aa < 1.d-30) return
  disc = bb*bb - 4.d0*aa*cc
  if (disc < 0.d0) return
  sqrt_disc = sqrt(disc)
  t = (-bb - sqrt_disc)/(2.d0*aa)
  if (t >= 0.d0 .and. t <= 1.d0) then
    PatchSegmentHitsEllipse = PETSC_TRUE
    return
  endif
  t = (-bb + sqrt_disc)/(2.d0*aa)
  if (t >= 0.d0 .and. t <= 1.d0) then
    PatchSegmentHitsEllipse = PETSC_TRUE
  endif

end function PatchSegmentHitsEllipse

! ************************************************************************** !

function PatchSegmentHitsRectangle(x1,y1,x2,y2,a,b)
  !
  ! True if segment intersects the filled rectangle [-a,a] x [-b,b].
  !

  implicit none

  PetscReal :: x1, y1, x2, y2
  PetscReal :: a, b

  PetscBool :: PatchSegmentHitsRectangle

  PetscReal :: dx
  PetscReal :: dy
  PetscReal :: t
  PetscReal :: u
  PetscReal :: qx
  PetscReal :: qy
  PetscReal :: rx
  PetscReal :: ry
  PetscReal :: sx
  PetscReal :: sy
  PetscReal :: rxs
  PetscInt :: k
  PetscReal :: ex1(4)
  PetscReal :: ey1(4)
  PetscReal :: ex2(4)
  PetscReal :: ey2(4)

  PatchSegmentHitsRectangle = PETSC_FALSE
  if (abs(x1) <= a .and. abs(y1) <= b) then
    PatchSegmentHitsRectangle = PETSC_TRUE
    return
  endif
  if (abs(x2) <= a .and. abs(y2) <= b) then
    PatchSegmentHitsRectangle = PETSC_TRUE
    return
  endif
  ex1 = (/-a, a, -a, -a/)
  ey1 = (/-b, -b, -b, b/)
  ex2 = (/a, a, -a, a/)
  ey2 = (/-b, b, b, b/)
  dx = x2 - x1
  dy = y2 - y1
  do k = 1, 4
    rx = dx
    ry = dy
    sx = ex2(k) - ex1(k)
    sy = ey2(k) - ey1(k)
    rxs = rx*sy - ry*sx
    if (abs(rxs) < 1.d-30) cycle
    qx = ex1(k) - x1
    qy = ey1(k) - y1
    t = (qx*sy - qy*sx)/rxs
    u = (qx*ry - qy*rx)/rxs
    if (t >= 0.d0 .and. t <= 1.d0 .and. u >= 0.d0 .and. u <= 1.d0) then
      PatchSegmentHitsRectangle = PETSC_TRUE
      return
    endif
  enddo

end function PatchSegmentHitsRectangle

! ************************************************************************** !

function PatchCircleHitsEllipse(xi,eta,radius,a,b)
  !
  ! True if a circle overlaps a filled axis-aligned ellipse.
  !

  implicit none

  PetscReal :: xi
  PetscReal :: eta
  PetscReal :: radius
  PetscReal :: a
  PetscReal :: b

  PetscBool :: PatchCircleHitsEllipse

  PetscReal :: px
  PetscReal :: py
  PetscReal :: u
  PetscReal :: v

  PatchCircleHitsEllipse = PETSC_FALSE
  if ((xi/a)**2 + (eta/b)**2 <= 1.d0) then
    PatchCircleHitsEllipse = PETSC_TRUE
    return
  endif
  u = abs(xi)
  v = abs(eta)
  call PatchClosestOnEllipse(u,v,a,b,px,py)
  if ((px-u)**2 + (py-v)**2 <= radius*radius) then
    PatchCircleHitsEllipse = PETSC_TRUE
  endif

end function PatchCircleHitsEllipse

! ************************************************************************** !

function PatchCircleHitsRectangle(xi,eta,radius,a,b)
  !
  ! True if a circle overlaps the rectangle [-a,a] x [-b,b].
  !

  implicit none

  PetscReal :: xi
  PetscReal :: eta
  PetscReal :: radius
  PetscReal :: a
  PetscReal :: b

  PetscBool :: PatchCircleHitsRectangle

  PetscReal :: qx
  PetscReal :: qy

  qx = min(a,max(-a,xi))
  qy = min(b,max(-b,eta))
  PatchCircleHitsRectangle = ((xi-qx)**2+(eta-qy)**2 <= radius*radius)

end function PatchCircleHitsRectangle

! ************************************************************************** !

subroutine PatchClosestOnEllipse(u,v,a,b,x,y)
  !
  ! Closest point (x,y) on the ellipse x^2/a^2+y^2/b^2=1 to (u,v)
  ! in the first quadrant (u>=0, v>=0).
  !

  implicit none

  PetscReal :: u
  PetscReal :: v
  PetscReal :: a
  PetscReal :: b
  PetscReal :: x
  PetscReal :: y

  PetscReal :: t
  PetscReal :: f
  PetscReal :: df
  PetscReal :: ta
  PetscReal :: tb
  PetscInt :: iter

  if (v < 1.d-14) then
    x = a
    y = 0.d0
    return
  endif
  if (u < 1.d-14) then
    x = 0.d0
    y = b
    return
  endif
  t = b*v - b*b
  do iter = 1, 25
    ta = t + a*a
    tb = t + b*b
    if (abs(ta) < 1.d-30) ta = 1.d-30
    if (abs(tb) < 1.d-30) tb = 1.d-30
    f = (a*u/ta)**2 + (b*v/tb)**2 - 1.d0
    df = -2.d0*((a*u)**2/ta**3 + (b*v)**2/tb**3)
    if (abs(df) < 1.d-30) exit
    t = t - f/df
    if (t < -b*b + 1.d-8) t = -b*b + 1.d-8
    if (abs(f) < 1.d-12) exit
  enddo
  x = a*a*u/(t+a*a)
  y = b*b*v/(t+b*b)

end subroutine PatchClosestOnEllipse

! ************************************************************************** !

subroutine RegionCopyPlanarPatch(patch_in,patch_out)
  !
  ! Deep-copies a planar patch object.
  !
  ! Author: Glenn Hammond
  ! Date: 09/11/26
  !

  implicit none

  type(planar_patch_type) :: patch_in
  type(planar_patch_type), pointer :: patch_out

  patch_out => RegionCreatePlanarPatch()
  patch_out%centroid = patch_in%centroid
  patch_out%normal = patch_in%normal
  patch_out%axis1 = patch_in%axis1
  patch_out%axis2 = patch_in%axis2
  patch_out%radii = patch_in%radii
  patch_out%half_thickness = patch_in%half_thickness
  patch_out%shape = patch_in%shape

end subroutine RegionCopyPlanarPatch

! ************************************************************************** !

subroutine RegionDestroyPlanarPatch(patch)
  !
  ! Deallocates a planar patch object.
  !
  ! Author: Glenn Hammond
  ! Date: 09/11/26
  !

  implicit none

  type(planar_patch_type), pointer :: patch

  if (.not.associated(patch)) return
  deallocate(patch)
  nullify(patch)

end subroutine RegionDestroyPlanarPatch

! ************************************************************************** !

subroutine PatchCrossProduct(a,b,c)
  !
  ! c = a x b for length-3 vectors.
  !
  implicit none

  PetscReal :: a(3)
  PetscReal :: b(3)
  PetscReal :: c(3)

  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)

end subroutine PatchCrossProduct

! ************************************************************************** !

subroutine RegionDestroy(region)
  !
  ! Deallocates a region
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(region_type), pointer :: region

  if (.not.associated(region)) return

  call DeallocateArray(region%cell_ids)
  call DeallocateArray(region%faces)
  if (associated(region%coordinates)) deallocate(region%coordinates)
  nullify(region%coordinates)
  call RegionDestroySideset(region%sideset)
  call RegionDestroyExplicitFaceSet(region%explicit_faceset)
  call GeometryDestroyPolygonalVolume(region%polygonal_volume)
  call RegionDestroyPlanarPatch(region%planar_patch)

  if (associated(region%vertex_ids)) deallocate(region%vertex_ids)
  nullify(region%vertex_ids)

  nullify(region%next)

  deallocate(region)
  nullify(region)

end subroutine RegionDestroy

end module Region_module
