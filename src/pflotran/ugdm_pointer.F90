module UGDM_Pointer_module

#include "petsc/finclude/petscdm.h"
#include "petsc/finclude/petscao.h"
  use petscdm
  use petscao

  implicit none

  private

  type, public :: ugdm_type
    ! local: included both local (non-ghosted) and ghosted cells
    ! global: includes only local (non-ghosted) cells
    PetscInt :: ndof
    ! for the below
    ! ghosted = local (non-ghosted) and ghosted cells
    ! local = local (non-ghosted) cells
    IS :: is_ghosted_local ! ghosted cells, local on-processor numbering
    IS :: is_local_local ! local cells, local on-processor numbering
    IS :: is_ghosted_petsc ! ghosted cells, petsc numbering
    IS :: is_local_petsc ! local cells, petsc numbering
    IS :: is_ghosts_local ! ghost cells, local on-processor numbering
    IS :: is_ghosts_petsc ! ghost cells, petsc numbering
    IS :: is_local_natural ! local cells, natural (global) numbering
    VecScatter :: scatter_ltog ! local to global
    VecScatter :: scatter_gtol ! global to local
    VecScatter :: scatter_ltol ! local to local
    VecScatter :: scatter_gton ! global to natural
    ISLocalToGlobalMapping :: mapping_ltog ! petsc vec local to global mapping
    Vec :: global_vec ! global vec (no ghost cells), petsc-ordering
    Vec :: local_vec ! local vec (local + ghosted cells), local ordering
    VecScatter :: scatter_bet_grids ! scatter between surface and subsurface
    VecScatter :: scatter_bet_grids_1dof ! scatter between surface and
                                         ! subsurface grids for 1-DOF
    VecScatter :: scatter_bet_grids_ndof ! scatter between surface and
                                         ! subsurface grids for N-DOFs
    AO :: ao_natural_to_petsc
  end type ugdm_type

  ! Subsurface PETSc DM plus optional unstructured scatter object (ugdm).
  ! Structured uses %dm only (%ugdm null). Unstructured may use both.
  ! %dm may be PETSC_NULL (PetscObjectIsNull) when only ugdm is needed.
  type, public :: dm_ptr_type
    DM :: dm
    type(ugdm_type), pointer :: ugdm
  end type dm_ptr_type

  public :: UGDMCreate, &
            UGridDMDestroy

contains

! ************************************************************************** !

function UGDMCreate()
  !
  ! Creates an unstructured grid distributed mesh object
  !
  ! Author: Glenn Hammond
  ! Date: 10/21/09
  !
  implicit none

  type(ugdm_type), pointer :: UGDMCreate

  type(ugdm_type), pointer :: ugdm

  allocate(ugdm)
  PetscObjectNullify(ugdm%is_ghosted_local)
  PetscObjectNullify(ugdm%is_local_local)
  PetscObjectNullify(ugdm%is_ghosted_petsc)
  PetscObjectNullify(ugdm%is_local_petsc)
  PetscObjectNullify(ugdm%is_ghosts_local)
  PetscObjectNullify(ugdm%is_ghosts_petsc)
  PetscObjectNullify(ugdm%is_local_natural)
  PetscObjectNullify(ugdm%scatter_ltog)
  PetscObjectNullify(ugdm%scatter_gtol)
  PetscObjectNullify(ugdm%scatter_ltol)
  PetscObjectNullify(ugdm%scatter_gton)
  PetscObjectNullify(ugdm%mapping_ltog)
  PetscObjectNullify(ugdm%global_vec)
  PetscObjectNullify(ugdm%local_vec)
  PetscObjectNullify(ugdm%scatter_bet_grids)
  PetscObjectNullify(ugdm%scatter_bet_grids_1dof)
  PetscObjectNullify(ugdm%scatter_bet_grids_ndof)
  ! this is solely a pointer, do not destroy
  PetscObjectNullify(ugdm%ao_natural_to_petsc)
  UGDMCreate => ugdm

end function UGDMCreate

! ************************************************************************** !

subroutine UGridDMDestroy(ugdm)
  !
  ! Deallocates a unstructured grid distributed mesh
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/09
  !
  use Petsc_Utility_module

  implicit none

  type(ugdm_type), pointer :: ugdm

  if (.not.associated(ugdm)) return

  call PUISDestroy(ugdm%is_ghosted_local)
  call PUISDestroy(ugdm%is_local_local)
  call PUISDestroy(ugdm%is_ghosted_petsc)
  call PUISDestroy(ugdm%is_local_petsc)
  call PUISDestroy(ugdm%is_ghosts_local)
  call PUISDestroy(ugdm%is_ghosts_petsc)
  call PUISDestroy(ugdm%is_local_natural)
  call PUVecScatterDestroy(ugdm%scatter_ltog)
  call PUVecScatterDestroy(ugdm%scatter_gtol)
  call PUVecScatterDestroy(ugdm%scatter_ltol)
  call PUVecScatterDestroy(ugdm%scatter_gton)
  call PUISLocalToGlobalMappingDestroy(ugdm%mapping_ltog)
  call PUVecDestroy(ugdm%global_vec)
  call PUVecDestroy(ugdm%local_vec)
  call PUVecScatterDestroy(ugdm%scatter_bet_grids)
  call PUVecScatterDestroy(ugdm%scatter_bet_grids_1dof)
  call PUVecScatterDestroy(ugdm%scatter_bet_grids_ndof)
  ! ugdm%ao_natural_to_petsc is a pointer to ugrid%ao_natural_to_petsc.  Do
  ! not destroy here.
  PetscObjectNullify(ugdm%ao_natural_to_petsc)
  deallocate(ugdm)
  nullify(ugdm)

end subroutine UGridDMDestroy

end module UGDM_Pointer_module
