module Communicator_Base_class

#include "petsc/finclude/petscvec.h"
   use petscvec
   use petscsys
   use UGDM_Pointer_module, only : dm_ptr_type

   use PFLOTRAN_Constants_module

  implicit none

  private

  type, abstract, public :: communicator_type
    type(dm_ptr_type) :: dm_ptr
  contains
    procedure, public :: SetDM => CommunicatorBaseSetDM
    procedure(VecToVec), public, deferred :: GlobalToLocal
    procedure(VecToVec), public, deferred :: LocalToGlobal
    procedure(VecToVec), public, deferred :: LocalToLocal
    procedure(VecToVec), public, deferred :: GlobalToNatural
    procedure(VecToVec), public, deferred :: NaturalToGlobal
    procedure(MapArray), public, deferred :: AONaturalToPetsc
    procedure(BaseDestroy), public, deferred :: Destroy
  end type communicator_type

  abstract interface

    subroutine VecToVec(this,source,destination)
      use petscvec
      import communicator_type
      implicit none
      class(communicator_type) :: this
      Vec :: source
      Vec :: destination
    end subroutine VecToVec

    subroutine MapArray(this,array)
      use petscvec
      import communicator_type
      implicit none
      class(communicator_type) :: this
      PetscInt :: array(:)
    end subroutine MapArray

    subroutine BaseDestroy(this)
      use petscvec
      import communicator_type
      implicit none
      class(communicator_type) :: this
    end subroutine BaseDestroy

  end interface

contains

! ************************************************************************** !

subroutine CommunicatorBaseSetDM(this,dm_ptr)
  implicit none

  class(communicator_type) :: this
  type(dm_ptr_type) :: dm_ptr

  this%dm_ptr = dm_ptr

end subroutine CommunicatorBaseSetDM

end module Communicator_Base_class
