module Immiscible_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module
  use Material_Aux_module

  implicit none

  private

  PetscReal, public :: immis_fmw(2) = [FMWH2O,FMWAIR]

  PetscReal, public :: immis_rel_pert = 1.d-8
  PetscReal, public :: immis_pres_min_pert = 1.d-2
  PetscReal, public :: immis_sat_min_pert = 1.d-10

  PetscReal, pointer, public :: immis_min_pert(:)

  PetscBool, public :: immis_calc_accum = PETSC_TRUE
  PetscBool, public :: immis_calc_flux = PETSC_TRUE
  PetscBool, public :: immis_calc_bcflux = PETSC_TRUE
  PetscBool, public :: immis_simultaneous_res_jac_calc = PETSC_TRUE

  PetscBool, public :: immis_numerical_derivatives = PETSC_FALSE

  ! debugging
  PetscInt, public :: immis_ni_count
  PetscInt, public :: immis_ts_cut_count
  PetscInt, public :: immis_ts_count
  PetscInt, public :: immis_debug_cell_id

  ! process models
  PetscInt, parameter, public :: IMMIS_GAS_PRESSURE_DOF = 1
  PetscInt, parameter, public :: IMMIS_GAS_SATURATION_DOF = 2
  PetscInt, parameter, public :: IMMIS_MAX_DOF = IMMIS_GAS_SATURATION_DOF

  PetscInt, parameter, public :: IMMIS_LIQUID_EQUATION_INDEX = 1
  PetscInt, parameter, public :: IMMIS_GAS_EQUATION_INDEX = 2

  PetscInt, parameter, public :: IMMIS_GAS_PRESSURE_INDEX = 1
  PetscInt, parameter, public :: IMMIS_GAS_SATURATION_INDEX = 2
  PetscInt, parameter, public :: IMMIS_LIQUID_FLUX_INDEX = 3
  PetscInt, parameter, public :: IMMIS_GAS_FLUX_INDEX = 4
  PetscInt, parameter, public :: IMMIS_RATE_SCALE_INDEX = 5
  PetscInt, parameter, public :: IMMIS_MAX_INDEX = IMMIS_RATE_SCALE_INDEX

  PetscInt, parameter, public :: IMMIS_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: IMMIS_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: IMMIS_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: IMMIS_UPDATE_FOR_BOUNDARY = 2

  PetscInt, public :: dof_to_primary_variable(2)

  type, public :: immiscible_auxvar_type
    PetscReal :: pres(3)
    PetscReal :: sat(2)
    PetscReal :: den_kg(2)
    PetscReal :: mu(2)
    PetscReal :: kr(2) ! relative permeability
    PetscReal :: pert
  end type immiscible_auxvar_type

  type, public :: immiscible_parameter_type
    PetscBool :: check_post_converged
  end type immiscible_parameter_type

  type, public :: immiscible_type
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(immiscible_parameter_type), pointer :: immiscible_parameter
    type(immiscible_auxvar_type), pointer :: auxvars(:,:)
    type(immiscible_auxvar_type), pointer :: auxvars_bc(:)
    type(immiscible_auxvar_type), pointer :: auxvars_ss(:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type immiscible_type

  interface ImmiscibleAuxVarDestroy
    module procedure ImmiscibleAuxVarSingleDestroy
    module procedure ImmiscibleAuxVarArray1Destroy
    module procedure ImmiscibleAuxVarArray2Destroy
  end interface ImmiscibleAuxVarDestroy

  interface ImmiscibleOutputAuxVars
    module procedure ImmiscibleOutputAuxVars1
  end interface ImmiscibleOutputAuxVars

  public :: ImmiscibleAuxCreate, &
            ImmiscibleAuxDestroy, &
            ImmiscibleAuxVarCompute, &
            ImmiscibleAuxVarInit, &
            ImmiscibleAuxVarCopy, &
            ImmiscibleAuxVarDestroy, &
            ImmiscibleAuxVarStrip, &
            ImmiscibleAuxVarPerturb, &
            ImmisciblePrintAuxVars, &
            ImmiscibleOutputAuxVars

contains

! ************************************************************************** !

function ImmiscibleAuxCreate(option)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24

  use Option_module

  implicit none

  type(option_type) :: option

  type(immiscible_type), pointer :: ImmiscibleAuxCreate

  type(immiscible_type), pointer :: aux

  nullify(immis_min_pert)

  dof_to_primary_variable(1:2) = &
    reshape([IMMIS_GAS_PRESSURE_DOF,IMMIS_GAS_SATURATION_DOF], &
             shape(dof_to_primary_variable))

  allocate(aux)
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  nullify(aux%matrix_zeroing)

  allocate(aux%immiscible_parameter)
  aux%immiscible_parameter%check_post_converged = PETSC_FALSE

  ImmiscibleAuxCreate => aux

end function ImmiscibleAuxCreate

! ************************************************************************** !

subroutine ImmiscibleAuxVarInit(auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24

  use Option_module

  implicit none

  type(immiscible_auxvar_type) :: auxvar
  type(option_type) :: option

  auxvar%pres = 0.d0
  auxvar%sat = 0.d0
  auxvar%den_kg = 0.d0
  auxvar%mu = 0.d0
  auxvar%kr = 0.d0

  auxvar%pert = 0.d0

end subroutine ImmiscibleAuxVarInit

! ************************************************************************** !

subroutine ImmiscibleAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24

  use Option_module

  implicit none

  type(immiscible_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pres = auxvar%pres
  auxvar2%sat = auxvar%sat
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%mu = auxvar%mu
  auxvar2%kr = auxvar%kr

  auxvar2%pert = auxvar%pert

end subroutine ImmiscibleAuxVarCopy

! ************************************************************************** !

subroutine ImmiscibleAuxVarCompute(x,immis_auxvar,global_auxvar, &
                              material_auxvar,characteristic_curves, &
                              natural_id,update_porosity,option)
  !
  ! Computes auxiliary variables for each grid cell
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24

  use Characteristic_Curves_module
  use EOS_Gas_module
  use EOS_Water_module
  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(:)
  type(immiscible_auxvar_type) :: immis_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscBool :: update_porosity
  PetscInt :: natural_id

  PetscInt, parameter :: lid = 1
  PetscInt, parameter :: gid = 2
  PetscInt, parameter :: cpid = 3
  PetscReal :: dpc_dsatl
  PetscReal :: den_air_mol
  PetscReal :: cell_pressure, saturation_pressure
  PetscReal :: dummy, dkrl_dsatl, dkrg_dsatl
  PetscInt :: ierr

  immis_auxvar%pres = 0.d0
  immis_auxvar%sat = 0.d0
  immis_auxvar%den_kg = 0.d0
  immis_auxvar%mu = 0.d0
  immis_auxvar%kr = 0.d0

  global_auxvar%temp = option%flow%reference_temperature

  immis_auxvar%pres(gid) = x(IMMIS_GAS_PRESSURE_DOF)
  immis_auxvar%sat(gid) = x(IMMIS_GAS_SATURATION_DOF)
  immis_auxvar%sat(lid) = 1.d0 - immis_auxvar%sat(gid)
  call characteristic_curves%saturation_function% &
         CapillaryPressure(immis_auxvar%sat(lid),immis_auxvar%pres(cpid), &
                           dpc_dsatl,option)
  immis_auxvar%pres(lid) = immis_auxvar%pres(gid) - immis_auxvar%pres(cpid)

  cell_pressure = immis_auxvar%pres(gid)
  call EOSWaterDensity(global_auxvar%temp,cell_pressure, &
                       immis_auxvar%den_kg(lid),dummy,ierr)
  call EOSGasDensity(global_auxvar%temp,cell_pressure,den_air_mol,ierr)
  immis_auxvar%den_kg(gid) = den_air_mol*immis_fmw(gid)

  call EOSWaterSaturationPressure(global_auxvar%temp, &
                                  saturation_pressure,dummy,ierr)
  call EOSWaterViscosity(global_auxvar%temp,cell_pressure, &
                          saturation_pressure,immis_auxvar%mu(lid),ierr)
  call EOSGasViscosity(global_auxvar%temp,immis_auxvar%pres(gid), &
                       immis_auxvar%pres(gid),den_air_mol, &
                       immis_auxvar%mu(gid),ierr)

  call characteristic_curves%liq_rel_perm_function% &
           RelativePermeability(immis_auxvar%sat(lid),immis_auxvar%kr(lid), &
                                dkrl_dsatl,option)
  call characteristic_curves%gas_rel_perm_function% &
           RelativePermeability(immis_auxvar%sat(lid),immis_auxvar%kr(gid), &
                                dkrg_dsatl,option)

end subroutine ImmiscibleAuxVarCompute

! ************************************************************************** !

subroutine ImmiscibleAuxVarPerturb(x,immis_auxvar,global_auxvar, &
                              material_auxvar, &
                              characteristic_curves, &
                              natural_id, &
                              option)
  ! Calculates auxiliary variables for perturbed system
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  PetscReal :: x(:)
  type(option_type) :: option
  PetscInt :: natural_id
  type(immiscible_auxvar_type) :: immis_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves

  PetscInt :: idof
  PetscReal :: x_pert(IMMIS_MAX_DOF), pert

  ! IMMIS_UPDATE_FOR_DERIVATIVE indicates call from perturbation
  option%iflag = IMMIS_UPDATE_FOR_DERIVATIVE
  do idof = 1, option%nflowdof
    pert = x(idof)*immis_rel_pert+immis_min_pert(idof)
    if (dabs(pert) < 1.d-40) pert = immis_min_pert(idof)
    immis_auxvar(idof)%pert = pert
    x_pert(1:option%nflowdof) = x
    x_pert(idof) = x(idof) + pert
    call ImmiscibleAuxVarCompute(x_pert,immis_auxvar(idof),global_auxvar, &
                            material_auxvar, &
                            characteristic_curves,natural_id, &
                            PETSC_TRUE,option)
  enddo

end subroutine ImmiscibleAuxVarPerturb

! ************************************************************************** !

subroutine ImmisciblePrintAuxVars(immis_auxvar,global_auxvar,material_auxvar, &
                             natural_id,string,option)
  !
  ! Prints out the contents of an auxvar
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(immiscible_auxvar_type) :: immis_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  type(option_type) :: option

  print *, '--------------------------------------------------------'
  print *, trim(string)
  print *, '                 cell id: ', natural_id
  print *, '         liquid pressure: ', immis_auxvar%pres(1)
  print *, '            gas pressure: ', immis_auxvar%pres(2)
  print *, '      capillary pressure: ', immis_auxvar%pres(3)
  print *, '       liquid saturation: ', immis_auxvar%sat(1)
  print *, '          gas saturation: ', immis_auxvar%sat(2)
  print *, '        liquid viscosity: ', immis_auxvar%mu(1)
  print *, '           gas viscosity: ', immis_auxvar%mu(2)
  print *, ' liquid rel permeability: ', immis_auxvar%kr(1)
  print *, '    gas rel permeability: ', immis_auxvar%kr(2)
  print *, '--------------------------------------------------------'

end subroutine ImmisciblePrintAuxVars

! ************************************************************************** !

subroutine ImmiscibleOutputAuxVars1(immis_auxvar,global_auxvar,material_auxvar,&
                               natural_id,string,append,option)
  !
  ! Prints out the contents of an auxvar to a file
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Global_Aux_module
  use Material_Aux_module
  use Option_module

  implicit none

  type(immiscible_auxvar_type) :: immis_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  PetscBool :: append
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string2

  write(string2,*) natural_id
  string2 = trim(adjustl(string)) // '_' // trim(adjustl(string2)) // '.txt'
  if (append) then
    open(unit=IUNIT_TEMP,file=string2,position='append')
  else
    open(unit=IUNIT_TEMP,file=string2)
  endif

  write(IUNIT_TEMP,*) '-------------------------------------------------------'
  write(IUNIT_TEMP,*) trim(string)
  write(IUNIT_TEMP,*) '                 cell id: ', natural_id
  write(IUNIT_TEMP,*) '         liquid pressure: ', immis_auxvar%pres(1)
  write(IUNIT_TEMP,*) '            gas pressure: ', immis_auxvar%pres(2)
  write(IUNIT_TEMP,*) '      capillary pressure: ', immis_auxvar%pres(3)
  write(IUNIT_TEMP,*) '       liquid saturation: ', immis_auxvar%sat(1)
  write(IUNIT_TEMP,*) '          gas saturation: ', immis_auxvar%sat(2)
  write(IUNIT_TEMP,*) '        liquid viscosity: ', immis_auxvar%mu(1)
  write(IUNIT_TEMP,*) '           gas viscosity: ', immis_auxvar%mu(2)
  write(IUNIT_TEMP,*) ' liquid rel permeability: ', immis_auxvar%kr(1)
  write(IUNIT_TEMP,*) '    gas rel permeability: ', immis_auxvar%kr(2)
  write(IUNIT_TEMP,*) '-------------------------------------------------------'
  close(IUNIT_TEMP)

end subroutine ImmiscibleOutputAuxVars1

! ************************************************************************** !

subroutine ImmiscibleAuxVarSingleDestroy(auxvar)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  implicit none

  type(immiscible_auxvar_type), pointer :: auxvar

  if (associated(auxvar)) then
    call ImmiscibleAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine ImmiscibleAuxVarSingleDestroy

! ************************************************************************** !

subroutine ImmiscibleAuxVarArray1Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  implicit none

  type(immiscible_auxvar_type), pointer :: auxvars(:)

  PetscInt :: iaux

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call ImmiscibleAuxVarStrip(auxvars(iaux))
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine ImmiscibleAuxVarArray1Destroy

! ************************************************************************** !

subroutine ImmiscibleAuxVarArray2Destroy(auxvars)
  !
  ! Deallocates a mode auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  implicit none

  type(immiscible_auxvar_type), pointer :: auxvars(:,:)

  PetscInt :: iaux, idof

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call ImmiscibleAuxVarStrip(auxvars(idof-1,iaux))
      enddo
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine ImmiscibleAuxVarArray2Destroy

! ************************************************************************** !

subroutine ImmiscibleMaterialAuxVarDestroy(auxvars)
  !
  ! Deallocates material auxiliary object for immiscible perturbation
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  implicit none

  type(material_auxvar_type), pointer :: auxvars(:,:)

  PetscInt :: iaux, idof

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call MaterialAuxVarStrip(auxvars(idof,iaux))
      enddo
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine ImmiscibleMaterialAuxVarDestroy

! ************************************************************************** !

subroutine ImmiscibleAuxVarStrip(auxvar)
  !
  ! ImmiscibleAuxVarDestroy: Deallocates a general auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(immiscible_auxvar_type) :: auxvar

end subroutine ImmiscibleAuxVarStrip

! ************************************************************************** !

subroutine ImmiscibleAuxDestroy(aux)
  !
  ! Deallocates a general auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(immiscible_type), pointer :: aux

  call DeallocateArray(immis_min_pert)

  if (.not.associated(aux)) return

  call ImmiscibleAuxVarDestroy(aux%auxvars)
  call ImmiscibleAuxVarDestroy(aux%auxvars_bc)
  call ImmiscibleAuxVarDestroy(aux%auxvars_ss)

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%immiscible_parameter)) then
    deallocate(aux%immiscible_parameter)
  endif
  nullify(aux%immiscible_parameter)

  deallocate(aux)
  nullify(aux)

end subroutine ImmiscibleAuxDestroy

end module Immiscible_Aux_module
