module Geomechanics_Linear_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: geomech_linear_parameter_type

    PetscBool :: youngs_modulus_spatially_varying
    PetscBool :: poissons_ratio_spatially_varying
    PetscBool :: density_spatially_varying
    PetscBool :: biot_coeff_spatially_varying
    PetscBool :: thermal_exp_coeff_spatially_varying

    PetscInt :: press_0_id
    PetscInt :: temp_0_id
    PetscInt :: vol_strain_0_id
    PetscInt :: vol_strain_id
    PetscInt :: stored_pressure_id
    PetscInt :: stored_temperature_id
    PetscInt :: stored_porosity_id
    PetscInt :: flow_porosity_id

    PetscReal, pointer :: youngs_modulus(:)
    PetscReal, pointer :: poissons_ratio(:)
    PetscReal, pointer :: biot_coeff(:)
    PetscReal, pointer :: thermal_exp_coeff(:)
    PetscReal, pointer :: density(:)

  end type geomech_linear_parameter_type

  type, public :: geomech_linear_type
    type(geomech_linear_parameter_type), pointer :: linear_parameter
  end type geomech_linear_type

  public :: GeomechLinearAuxCreate, &
            GeomechLinearAuxDestroy

contains

! ************************************************************************** !

function GeomechLinearAuxCreate()
  !
  ! Nullifies pointers in geomech linear auxiliary type
  !
  ! Author: Glenn Hammond
  ! Date: 01/22/26

  implicit none

  type(geomech_linear_type), pointer :: GeomechLinearAuxCreate

  type(geomech_linear_type), pointer :: aux
  type(geomech_linear_parameter_type), pointer :: linear_parameter

  allocate(aux)
  allocate(linear_parameter)
  linear_parameter%youngs_modulus_spatially_varying = PETSC_FALSE
  linear_parameter%poissons_ratio_spatially_varying = PETSC_FALSE
  linear_parameter%density_spatially_varying = PETSC_FALSE
  linear_parameter%biot_coeff_spatially_varying = PETSC_FALSE
  linear_parameter%thermal_exp_coeff_spatially_varying = PETSC_FALSE
  linear_parameter%press_0_id = UNINITIALIZED_INTEGER
  linear_parameter%temp_0_id = UNINITIALIZED_INTEGER
  linear_parameter%vol_strain_0_id = UNINITIALIZED_INTEGER
  linear_parameter%vol_strain_id = UNINITIALIZED_INTEGER
  linear_parameter%stored_pressure_id = UNINITIALIZED_INTEGER
  linear_parameter%stored_temperature_id = UNINITIALIZED_INTEGER
  linear_parameter%stored_porosity_id = UNINITIALIZED_INTEGER
  linear_parameter%flow_porosity_id = UNINITIALIZED_INTEGER
  nullify(linear_parameter%youngs_modulus)
  nullify(linear_parameter%poissons_ratio)
  nullify(linear_parameter%biot_coeff)
  nullify(linear_parameter%thermal_exp_coeff)
  nullify(linear_parameter%density)

  aux%linear_parameter => linear_parameter

  GeomechLinearAuxCreate => aux

end function GeomechLinearAuxCreate

! ************************************************************************** !

subroutine GeomechLinearAuxDestroy(aux)
  !
  ! Deallocates a geomech Linear auxiliary type
  !
  ! Author: Glenn Hammond
  ! Date: 01/22/26
  !
  use Utility_module
  use Petsc_Utility_module

  implicit none

  type(geomech_linear_type), pointer :: aux

  type(geomech_linear_parameter_type), pointer :: linear_parameter

  if (.not.associated(aux)) return

  linear_parameter => aux%linear_parameter

  call DeallocateArray(linear_parameter%youngs_modulus)
  call DeallocateArray(linear_parameter%poissons_ratio)
  call DeallocateArray(linear_parameter%biot_coeff)
  call DeallocateArray(linear_parameter%thermal_exp_coeff)
  call DeallocateArray(linear_parameter%density)

  deallocate(aux%linear_parameter)
  nullify(aux%linear_parameter)

end subroutine GeomechLinearAuxDestroy

end module Geomechanics_Linear_Aux_module
