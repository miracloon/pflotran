module GMDM_Pointer_module

#include "petsc/finclude/petscdm.h"
  use petscdm

  implicit none

  private

  type, public :: gmdm_type
    ! local: included both local (non-ghosted) and ghosted nodes
    ! global: includes only local (non-ghosted) nodes
    PetscInt :: ndof
    ! for the below
    ! ghosted = local (non-ghosted) and ghosted nodes
    ! local = local (non-ghosted) nodes
    IS :: is_ghosted_local ! ghosted nodes, local on-processor numbering
    IS :: is_local_local ! local nodes, local on-processor numbering
    IS :: is_ghosted_petsc ! ghosted nodes, petsc numbering
    IS :: is_local_petsc ! local nodes, petsc numbering
    IS :: is_ghosts_local ! ghost nodes, local on-processor numbering
    IS :: is_ghosts_petsc ! ghost nodes, petsc numbering
    IS :: is_local_natural ! local nodes, natural (global) numbering
    VecScatter :: scatter_ltog ! local to global
    VecScatter :: scatter_gtol ! global to local
    VecScatter :: scatter_ltol ! local to local
    VecScatter :: scatter_gton ! global to natural
    VecScatter :: scatter_gton_elem ! global to natural for elements
    ISLocalToGlobalMapping :: mapping_ltog ! petsc vec local to global mapping
    ISLocalToGlobalMapping :: mapping_ltog_elem ! for elements
    Vec :: global_vec ! global vec (no ghost nodes), petsc-ordering
    Vec :: local_vec ! local vec (local + ghosted nodes), local ordering
    Vec :: global_vec_elem
    VecScatter :: scatter_subsurf_to_geomech_ndof ! subsurface to geomech
    VecScatter :: scatter_geomech_to_subsurf_ndof ! geomech to subsurface
  end type gmdm_type

  ! PETSc DM plus the PFLOTRAN geomech scatter/mapping object (gmdm).
  ! %dm may be PETSC_NULL (PetscObjectIsNull) when only gmdm is needed.
  type, public :: gmdm_ptr_type
    DM :: dm
    type(gmdm_type), pointer :: gmdm
  end type gmdm_ptr_type

  public :: GMDMCreate, &
            GMDMDestroy

contains

! ************************************************************************** !

function GMDMCreate()
  !
  ! Creates a geomech grid distributed mesh object
  !
  ! Author: Satish Karra, LANL
  ! Date: 05/22/13
  !
  implicit none

  type(gmdm_type), pointer :: GMDMCreate
  type(gmdm_type), pointer :: gmdm

  allocate(gmdm)
  PetscObjectNullify(gmdm%is_ghosted_local)
  PetscObjectNullify(gmdm%is_local_local)
  PetscObjectNullify(gmdm%is_ghosted_petsc)
  PetscObjectNullify(gmdm%is_local_petsc)
  PetscObjectNullify(gmdm%is_ghosts_local)
  PetscObjectNullify(gmdm%is_ghosts_petsc)
  PetscObjectNullify(gmdm%is_local_natural)
  PetscObjectNullify(gmdm%scatter_ltog)
  PetscObjectNullify(gmdm%scatter_gtol)
  PetscObjectNullify(gmdm%scatter_ltol)
  PetscObjectNullify(gmdm%scatter_gton)
  PetscObjectNullify(gmdm%scatter_gton_elem)
  PetscObjectNullify(gmdm%mapping_ltog)
  PetscObjectNullify(gmdm%mapping_ltog_elem)
  PetscObjectNullify(gmdm%global_vec)
  PetscObjectNullify(gmdm%local_vec)
  PetscObjectNullify(gmdm%global_vec_elem)
  PetscObjectNullify(gmdm%scatter_subsurf_to_geomech_ndof)
  PetscObjectNullify(gmdm%scatter_geomech_to_subsurf_ndof)

  GMDMCreate => gmdm

end function GMDMCreate

! ************************************************************************** !

subroutine GMDMDestroy(gmdm)
  !
  ! Deallocates a geomechanics grid distributed mesh object
  !
  ! Author: Satish Karra, LANL
  ! Date: 05/22/13
  !
  use Petsc_Utility_module

  implicit none

  type(gmdm_type), pointer :: gmdm

  if (.not.associated(gmdm)) return

  call PUISDestroy(gmdm%is_ghosted_local)
  call PUISDestroy(gmdm%is_local_local)
  call PUISDestroy(gmdm%is_ghosted_petsc)
  call PUISDestroy(gmdm%is_local_petsc)
  call PUISDestroy(gmdm%is_ghosts_local)
  call PUISDestroy(gmdm%is_ghosts_petsc)
  call PUISDestroy(gmdm%is_local_natural)
  call PUVecScatterDestroy(gmdm%scatter_ltog)
  call PUVecScatterDestroy(gmdm%scatter_gtol)
  call PUVecScatterDestroy(gmdm%scatter_ltol)
  call PUVecScatterDestroy(gmdm%scatter_gton)
  call PUVecScatterDestroy(gmdm%scatter_gton_elem)
  call PUISLocalToGlobalMappingDestroy(gmdm%mapping_ltog)
  call PUISLocalToGlobalMappingDestroy(gmdm%mapping_ltog_elem)
  call PUVecDestroy(gmdm%global_vec)
  call PUVecDestroy(gmdm%local_vec)
  call PUVecDestroy(gmdm%global_vec_elem)
  call PUVecScatterDestroy(gmdm%scatter_subsurf_to_geomech_ndof)
  call PUVecScatterDestroy(gmdm%scatter_geomech_to_subsurf_ndof)

  deallocate(gmdm)
  nullify(gmdm)

end subroutine GMDMDestroy

end module GMDM_Pointer_module
