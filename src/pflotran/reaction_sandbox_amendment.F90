module Reaction_Sandbox_Amendment_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  PetscInt, parameter :: MODEL_CONSTANT = 0
  PetscInt, parameter :: MODEL_FILTRATION = 1
  PetscInt, parameter :: MODEL_TEMPERATURE = 2

  private

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_amendment_type

    !ID number of amendment species
    PetscInt :: species_Vaq_id ! Aqueous species
    PetscInt :: species_Vim_id ! Immobile species
    PetscInt :: viscosity_id   ! index in global%parameters for viscosity

    !Name of amendment species
    character(len=MAXWORDLENGTH) :: name_aqueous
    character(len=MAXWORDLENGTH) :: name_immobile

    PetscInt :: aqueous_decay_model
    PetscReal :: aqueous_decay_rate
    PetscReal :: aqueous_logDref
    PetscReal :: aqueous_Tref
    PetscReal :: aqueous_zT
    PetscReal :: aqueous_n

    PetscInt :: adsorbed_decay_model
    PetscReal :: adsorbed_decay_rate
    PetscReal :: adsorbed_logDref
    PetscReal :: adsorbed_Tref
    PetscReal :: adsorbed_zT
    PetscReal :: adsorbed_n

    PetscInt :: attachment_model
    PetscReal :: attachment_rate_constant
    PetscReal :: collector_diameter
    PetscReal :: particle_diameter
    PetscReal :: hamaker_constant
    PetscReal :: BOLTZMANN_CONSTANT = 1.380649d-23 !J/K
                 !(Not found on pflotran_constants.f90)
    PetscReal :: particle_density
    PetscReal :: attachment_efficiency

    PetscInt :: detachment_model
    PetscReal :: detachment_rate_constant

    !Debug?
    PetscBool :: debug_option

  contains
    procedure, public :: ReadInput => AmendmentRead
    procedure, public :: Setup => AmendmentSetup
    procedure, public :: Evaluate => AmendmentReact
    procedure, public :: Destroy => AmendmentDestroy

  end type reaction_sandbox_amendment_type

  public :: AmendmentCreate

contains

! ************************************************************************** !

function AmendmentCreate()
  !
  ! Allocates particle transport variables.
  !
  ! Author: Edwin Saavedra C.
  ! Date: 10/01/2020
  !

  implicit none

  class(reaction_sandbox_amendment_type), pointer :: AmendmentCreate

  allocate(AmendmentCreate)

  !ID number of amendment species
  AmendmentCreate%species_Vaq_id = UNINITIALIZED_INTEGER
  AmendmentCreate%species_Vim_id = UNINITIALIZED_INTEGER
  AmendmentCreate%viscosity_id = UNINITIALIZED_INTEGER

  !Name of amendment species
  AmendmentCreate%name_aqueous = ''
  AmendmentCreate%name_immobile = ''

  AmendmentCreate%aqueous_decay_model = UNINITIALIZED_INTEGER
  AmendmentCreate%aqueous_decay_rate = UNINITIALIZED_DOUBLE
  AmendmentCreate%aqueous_logDref = UNINITIALIZED_DOUBLE
  AmendmentCreate%aqueous_Tref = UNINITIALIZED_DOUBLE
  AmendmentCreate%aqueous_zT = UNINITIALIZED_DOUBLE
  AmendmentCreate%aqueous_n = UNINITIALIZED_DOUBLE

  AmendmentCreate%adsorbed_decay_model = UNINITIALIZED_INTEGER
  AmendmentCreate%adsorbed_decay_rate = UNINITIALIZED_DOUBLE
  AmendmentCreate%adsorbed_logDref = UNINITIALIZED_DOUBLE
  AmendmentCreate%adsorbed_Tref = UNINITIALIZED_DOUBLE
  AmendmentCreate%adsorbed_zT = UNINITIALIZED_DOUBLE
  AmendmentCreate%adsorbed_n = UNINITIALIZED_DOUBLE

  AmendmentCreate%attachment_model = UNINITIALIZED_INTEGER
  AmendmentCreate%attachment_rate_constant = UNINITIALIZED_DOUBLE
  AmendmentCreate%collector_diameter = UNINITIALIZED_DOUBLE
  AmendmentCreate%particle_diameter = UNINITIALIZED_DOUBLE
  AmendmentCreate%hamaker_constant = UNINITIALIZED_DOUBLE
  AmendmentCreate%particle_density = UNINITIALIZED_DOUBLE
  AmendmentCreate%attachment_efficiency = UNINITIALIZED_DOUBLE

  AmendmentCreate%detachment_model = UNINITIALIZED_INTEGER
  AmendmentCreate%detachment_rate_constant = UNINITIALIZED_DOUBLE

  AmendmentCreate%debug_option = PETSC_FALSE

  nullify(AmendmentCreate%next)

end function AmendmentCreate

! ************************************************************************** !
subroutine AmendmentRead(this,input,option)
  !
  ! Reads input deck for reaction sandbox parameters
  !
  ! Author: Edwin Saavedra C.
  ! Date: 10/01/2020 - created
  !       05/18/2020 - fixed formatting
  !

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_amendment_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, internal_units, units
  character(len=MAXSTRINGLENGTH) :: error_string, error_string2, error_string3

  error_string = 'CHEMISTRY,REACTION_SANDBOX,AMENDMENT'
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(keyword)
      case('PARTICLE_NAME_AQ')
        ! Bioparticle name while in suspension
        call InputReadWord(input,option,this%name_aqueous,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)

      case('PARTICLE_NAME_IM')
        ! Bioparticle name while immobilized
        call InputReadWord(input,option,this%name_immobile,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)

      case('AQUEOUS_DECAY_MODEL')
        ! Decay rate while in the aqueous phase
        call InputReadCard(input,option,keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'AQUEOUS_DECAY_MODEL',error_string)
        error_string2 = trim(error_string) // ',AQUEOUS_DECAY_MODEL'
        call StringToUpper(keyword)
        select case(keyword)
          case('CONSTANT')
            this%aqueous_decay_model = MODEL_CONSTANT
            error_string3 = trim(error_string2) //',CONSTANT'
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit

              call InputReadCard(input,option,keyword)
              call InputErrorMsg(input,option,'keyword',error_string3)
              call StringToUpper(keyword)
              select case(keyword)
                case('RATE')
                ! Read the double precision rate constant
                  call InputReadDouble(input,option,this%aqueous_decay_rate)
                  call InputErrorMsg(input,option,keyword,error_string3)
                  ! Read the units
                  call InputReadWord(input,option,units,PETSC_TRUE)
                  if (InputError(input)) then
                  ! If units do not exist, assume default units of 1/s which are the
                  ! standard internal PFLOTRAN units for this rate constant.
                    input%err_buf = trim(error_string3) // ',RATE UNITS'
                    call InputDefaultMsg(input,option)
                  else
                    ! If units exist, convert to internal units of 1/s
                    internal_units = 'unitless/sec'
                    this%aqueous_decay_rate = this%aqueous_decay_rate * &
                      UnitsConvertToInternal(units,internal_units, &
                                  trim(error_string3)//',RATE UNITS', &
                                  option)
                  endif

                case default
                  call InputKeywordUnrecognized(input,keyword,error_string3, &
                                                option)
              end select
            end do
            call InputPopBlock(input,option)

          case('TEMPERATURE')
            this%aqueous_decay_model = MODEL_TEMPERATURE
            error_string3 = trim(error_string2) // ',TEMPERATURE'
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit

              call InputReadCard(input,option,keyword)
              call InputErrorMsg(input,option,'keyword', &
                                 error_string3)
              call StringToUpper(keyword)
              select case(keyword)
                case('TREF')
                ! Reference temperature (Probably 4°C)
                  call InputReadDouble(input,option,this%aqueous_Tref)
                  call InputErrorMsg(input,option,keyword,error_string3)
                case('ZT')
                ! Model parameter zT
                  call InputReadDouble(input,option,this%aqueous_zT)
                  call InputErrorMsg(input,option,keyword,error_string3)
                case('N')
                  ! Model parameter n (Probably 1.0 or 2.0 )
                  call InputReadDouble(input,option,this%aqueous_n)
                  call InputErrorMsg(input,option,keyword,error_string3)
                case('LOGDREF')
                ! D reference value (Probably 2.3)
                  call InputReadDouble(input,option,this%aqueous_logDref)
                  call InputErrorMsg(input,option,keyword,error_string3)
                case default
                  call InputKeywordUnrecognized(input,keyword,error_string3, &
                                                option)
              end select
            end do
            call InputPopBlock(input,option)

          case default
            call InputKeywordUnrecognized(input,keyword,error_string2,option)
        end select

      case('ADSORBED_DECAY_MODEL')
        ! Decay rate while in the immobile phase
        call InputReadCard(input,option,keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'ADSORBED_DECAY_MODEL',error_string)
        error_string2 = trim(error_string) // ',ADSORBED_DECAY_MODEL'
        call StringToUpper(keyword)
        select case(keyword)
          case('CONSTANT')
            this%adsorbed_decay_model = MODEL_CONSTANT
            error_string3 = trim(error_string2) // ',CONSTANT'
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit
              call InputReadCard(input,option,keyword)
              call InputErrorMsg(input,option,'keyword',error_string3)
              call StringToUpper(keyword)

              select case(keyword)
                case('RATE')
                ! Read the double precision rate constant
                  call InputReadDouble(input,option,this%adsorbed_decay_rate)
                  call InputErrorMsg(input,option,keyword,error_string3)
                  ! Read the units
                  call InputReadWord(input,option,units,PETSC_TRUE)
                  if (InputError(input)) then
                  ! If units do not exist, assume default units of 1/s which are the
                  ! standard internal PFLOTRAN units for this rate constant.
                    input%err_buf = trim(error_string3) // ',RATE CONSTANT UNITS'
                    call InputDefaultMsg(input,option)
                  else
                    ! If units exist, convert to internal units of 1/s
                    internal_units = 'unitless/sec'
                    this%adsorbed_decay_rate = this%adsorbed_decay_rate * &
                      UnitsConvertToInternal(units,internal_units, &
                             trim(error_string3) // &
                             ',RATE CONSTANT UNITS',option)
                  endif
                case default
                  call InputKeywordUnrecognized(input,keyword,error_string3, &
                                                option)
              end select
            end do
            call InputPopBlock(input,option)

          case('TEMPERATURE')
            this%adsorbed_decay_model = MODEL_TEMPERATURE
            error_string3 = trim(error_string2) // ',TEMPERATURE'
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit

              call InputReadCard(input,option,keyword)
              call InputErrorMsg(input,option,'keyword',error_string3)
              call StringToUpper(keyword)
              select case(keyword)
                case('TREF')
                  call InputReadDouble(input,option,this%adsorbed_Tref)
                  call InputErrorMsg(input,option,'TREF',error_string3)
                case('ZT')
                  call InputReadDouble(input,option,this%adsorbed_zT)
                  call InputErrorMsg(input,option,'ZT',error_string3)
                case('N')
                  call InputReadDouble(input,option,this%adsorbed_n)
                  call InputErrorMsg(input,option,'N',error_string3)
                case('LOGDREF')
                  call InputReadDouble(input,option,this%adsorbed_logDref)
                  call InputErrorMsg(input,option,'LOGDREF',error_string3)
                case default
                  call InputKeywordUnrecognized(input,keyword,error_string3, &
                                                option)
              end select
            end do
            call InputPopBlock(input,option)
          case default
            call InputKeywordUnrecognized(input,keyword,error_string2,option)
        end select

      case('ATTACHMENT_MODEL')
        ! Decay rate while in the aqueous phase
        call InputReadCard(input,option,keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'ATTACHMENT_MODEL', &
               error_string)
        error_string2 = trim(error_string) // ',ATTACHMENT_MODEL'
        call StringToUpper(keyword)

        select case(keyword)
          case('CONSTANT')
            this%attachment_model = MODEL_CONSTANT
            error_string3 = trim(error_string2) // ',CONSTANT'
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit

              call InputReadCard(input,option,keyword)
              call InputErrorMsg(input,option,'keyword',error_string3)
              call StringToUpper(keyword)

              select case(keyword)
                case('RATE')
                ! Read the double precision rate constant
                  call InputReadDouble(input,option,this%attachment_rate_constant)
                  call InputErrorMsg(input,option,'RATE',error_string3)
                  ! Read the units
                  call InputReadWord(input,option,units,PETSC_TRUE)
                  if (InputError(input)) then
                  ! If units do not exist, assume default units of 1/s which are the
                  ! standard internal PFLOTRAN units for this rate constant.
                    input%err_buf = trim(error_string3) // &
                                    ',RATE UNITS assumed as 1/s'
                    call InputDefaultMsg(input,option)
                  else
                    ! If units exist, convert to internal units of 1/s
                    internal_units = 'unitless/sec'
                    this%attachment_rate_constant = this%attachment_rate_constant * &
                      UnitsConvertToInternal(units,internal_units, &
                                             trim(error_string3) // &
                                               ',RATE UNITS', &
                                             option)
                  endif
                case default
                  call InputKeywordUnrecognized(input,keyword,error_string2, &
                                                option)
              end select
            end do
            call InputPopBlock(input,option)

          case('FILTRATION')
            this%attachment_model = MODEL_FILTRATION
            call InputPushBlock(input,option)
            ! Gotta turn this on so Darcy velocity is stored in global
            ! and can be used in the reaction sandbox
            option%flow%store_darcy_vel = PETSC_TRUE
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit

              call InputReadCard(input,option,keyword)
              call InputErrorMsg(input,option,'FILTRATION', &
                     error_string)
              error_string3 = trim(error_string2) // ',FILTRATION'
              call StringToUpper(keyword)
              select case(keyword)
                case('COLLECTOR_DIAMETER')
                ! Diameter of the collector, i.e., soil grain size [m]
                  call InputReadDouble(input,option,this%collector_diameter)
                  call InputErrorMsg(input,option,keyword,error_string3)

                case('PARTICLE_DIAMETER')
                ! Diameter of the amendment [m]
                  call InputReadDouble(input,option,this%particle_diameter)
                  call InputErrorMsg(input,option,keyword,error_string3)

                case('HAMAKER_CONSTANT')
                ! Hamaker constant particle-soil pair [Joules]
                  call InputReadDouble(input,option,this%hamaker_constant)
                  call InputErrorMsg(input,option,keyword,error_string3)

                case('PARTICLE_DENSITY')
                ! Density of the particulates [kg/m3]
                  call InputReadDouble(input,option,this%particle_density)
                  call InputErrorMsg(input,option,keyword,error_string3)

                case('ATTACHMENT_EFFICIENCY')
                ! Collision/attachment efficiency [-]
                  call InputReadDouble(input,option,this%attachment_efficiency)
                  call InputErrorMsg(input,option,keyword,error_string3)

                case('DEBUG')
                ! Print stuff on screen
                  print *, "Will debug -> print CFT stuff: "
                  this%debug_option = PETSC_TRUE ! Edwin debugging

                ! Something else
                case default
                  call InputKeywordUnrecognized(input,keyword, &
                                                error_string3,option)
              end select
            end do
            call InputPopBlock(input,option)

          case default
            call InputKeywordUnrecognized(input,keyword, &
                                          error_string2,option)
        end select

      ! Detachment rate
      case('DETACHMENT_MODEL')
        ! Decay rate while in the aqueous phase
        call InputReadCard(input,option,keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'DETACHMENT_MODEL', &
               error_string)
        error_string2 = trim(error_string) // ',DETACHMENT_MODEL'
        call StringToUpper(keyword)
        select case(keyword)
          case('CONSTANT')
            this%detachment_model = MODEL_CONSTANT
            error_string3 = trim(error_string2) // ',CONSTANT'
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit

              call InputReadCard(input,option,keyword)
              call InputErrorMsg(input,option,'CONSTANT',error_string2)
              call StringToUpper(keyword)

              select case(keyword)
                case('RATE')
                ! Read the double precision rate constant
                  call InputReadDouble(input,option,this%detachment_rate_constant)
                  call InputErrorMsg(input,option,'RATE',error_string3)
                  ! Read the units
                  call InputReadWord(input,option,units,PETSC_TRUE)
                  if (InputError(input)) then
                  ! If units do not exist, assume default units of 1/s which are the
                  ! standard internal PFLOTRAN units for this rate constant.
                    input%err_buf = trim(error_string3) // &
                                    'RATE UNITS assumed as 1/s'
                    call InputDefaultMsg(input,option)
                  else
                    ! If units exist, convert to internal units of 1/s
                    internal_units = 'unitless/sec'
                    this%detachment_rate_constant = this%detachment_rate_constant * &
                      UnitsConvertToInternal(units,internal_units, &
                                             error_string3,option)
                  endif
                case default
                  call InputKeywordUnrecognized(input,keyword, &
                                                error_string3,option)
              end select
            end do
            call InputPopBlock(input,option)

          case default
            call InputKeywordUnrecognized(input,keyword,error_string2,option)
        end select
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  end do

  call InputPopBlock(input,option)

end subroutine AmendmentRead

! ************************************************************************** !

subroutine AmendmentSetup(this,reaction,option)
  !
  ! Sets up the kinetic attachment/dettachment reactions
  !
  ! Author: Edwin Saavedra C.
  ! Date: 10/01/2020
  !

  use Reaction_Aux_module
  use Reaction_Immobile_Aux_module
  use Option_module
  use Parameter_module
  use Material_Aux_module, only : mean_soil_grain_size_index

  implicit none

  class(reaction_sandbox_amendment_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  this%species_Vaq_id = &
    ReactionAuxGetPriSpecIDFromName(this%name_aqueous,reaction,option)

  this%species_Vim_id = &
    ReactionImGetSpeciesIDFromName(this%name_immobile,reaction%immobile,option)

  this%viscosity_id = ParameterGetIDFromName('Viscosity',option)

  select case(this%attachment_model)
    case(MODEL_CONSTANT)
      if (Uninitialized(this%attachment_rate_constant)) then
        option%io_buffer = 'RATE must be defined for ATTACHMENT_MODEL CONSTANT.'
        call PrintErrMsg(option)
      endif
    case(MODEL_FILTRATION)
      if (Uninitialized(this%particle_diameter)) then
        option%io_buffer = 'ATTACHMENT_MODEL FILTRATION requires a PARTICLE_DIAMETER.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(this%hamaker_constant)) then
        option%io_buffer = 'ATTACHMENT_MODEL FILTRATION requires a HAMAKER_CONSTANT.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(this%particle_density)) then
        option%io_buffer = 'ATTACHMENT_MODEL FILTRATION requires a PARTICLE_DENSITY.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(this%attachment_efficiency)) then
        option%io_buffer = 'ATTACHMENT_MODEL FILTRATION requires an ATTACHMENT_EFFICIENCY.'
        call PrintErrMsg(option)
      endif
      if (mean_soil_grain_size_index == 0 .and. &
          Uninitialized(this%collector_diameter)) then
        option%io_buffer = 'FILTRATION_MODEL requires that COLLECTOR_DIAMETER be defined &
          &or MEAN_SOIL_GRAIN_SIZE be specified under MATERIAL_PROPERTY.'
        call PrintErrMsg(option)
      endif
      if (mean_soil_grain_size_index > 0 .and. &
          Initialized(this%collector_diameter)) then
        option%io_buffer = 'Both MEAN_SOIL_GRAIN_SIZE and COLLECTOR_DIAMETER cannot &
          &be specified for Reaction Sandbox AMENDMENT.'
        call PrintErrMsg(option)
      endif
    case default
      if (Initialized(this%attachment_model)) then
        option%io_buffer = 'A model must be specified for ATTACHMENT_MODEL.'
        call PrintErrMsg(option)
      endif
  end select

  select case(this%detachment_model)
    case(MODEL_CONSTANT)
      if (Uninitialized(this%detachment_rate_constant)) then
        option%io_buffer = 'RATE must be defined for DETACHMENT_MODEL CONSTANT.'
        call PrintErrMsg(option)
      endif
    case default
      if (Initialized(this%detachment_model)) then
        option%io_buffer = 'A model must be specified for DETACHMENT_MODEL.'
        call PrintErrMsg(option)
      endif
  end select

  select case(this%aqueous_decay_model)
    case(MODEL_CONSTANT)
      if (Uninitialized(this%aqueous_decay_rate)) then
        option%io_buffer = 'RATE must be defined for AQUEOUS_DECAY_MODEL CONSTANT.'
        call PrintErrMsg(option)
      endif
    case(MODEL_TEMPERATURE)
      if (Uninitialized(this%aqueous_logDref)) then
        option%io_buffer = 'AQUEOUS_DECAY_MODEL TEMPERATURE requires AQUEOUS_LOGDREF.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(this%aqueous_Tref)) then
        option%io_buffer = 'AQUEOUS_DECAY_MODEL TEMPERATURE requires AQUEOUS_TREF.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(this%aqueous_zT)) then
        option%io_buffer = 'AQUEOUS_DECAY_MODEL TEMPERATURE requires AQUEOUS_ZT.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(this%aqueous_N)) then
        option%io_buffer = 'AQUEOUS_DECAY_MODEL TEMPERATURE requires AQUEOUS_N.'
        call PrintErrMsg(option)
      endif
    case default
      if (Initialized(this%aqueous_decay_model)) then
        option%io_buffer = 'A model must be specified for AQUEOUS_DECAY_MODEL.'
        call PrintErrMsg(option)
      endif
  end select

  select case(this%adsorbed_decay_model)
    case(MODEL_CONSTANT)
      if (Uninitialized(this%adsorbed_decay_rate)) then
        option%io_buffer = 'RATE must be defined for ADSORBED_DECAY_MODEL CONSTANT.'
        call PrintErrMsg(option)
      endif
    case(MODEL_TEMPERATURE)
      if (Uninitialized(this%adsorbed_logDref)) then
        option%io_buffer = 'ADSORBED_DECAY_MODEL TEMPERATURE requires ADSORBED_LOGDREF.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(this%adsorbed_Tref)) then
        option%io_buffer = 'ADSORBED_DECAY_MODEL TEMPERATURE requires ADSORBED_TREF.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(this%adsorbed_zT)) then
        option%io_buffer = 'ADSORBED_DECAY_MODEL TEMPERATURE requires ADSORBED_ZT.'
        call PrintErrMsg(option)
      endif
      if (Uninitialized(this%adsorbed_N)) then
        option%io_buffer = 'ADSORBED_DECAY_MODEL TEMPERATURE requires ADSORBED_N.'
        call PrintErrMsg(option)
      endif
    case default
      if (Initialized(this%adsorbed_decay_model)) then
        option%io_buffer = 'A model must be specified for ADSORBED_DECAY_MODEL.'
        call PrintErrMsg(option)
      endif
  end select

end subroutine AmendmentSetup

! ************************************************************************** !

subroutine AmendmentReact(this,Residual,Jacobian,compute_derivative, &
                        rt_auxvar,global_auxvar,material_auxvar, &
                        reaction, option)
  !
  ! Evaluates reaction
  !
  ! Author: Edwin
  ! Date: 04/09/2020 - created
  !       05/18/2021 - remove truncate concentrations
  !

  use Option_module
  use String_module
  use Reaction_Aux_module, only : reaction_rt_type
  use Reaction_Immobile_Aux_module
  use Material_Aux_module, only : material_auxvar_type, &
                                  mean_soil_grain_size_index

  implicit none

  class(reaction_sandbox_amendment_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative

  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: volume                 ! m^3 bulk
  PetscReal :: porosity               ! m^3 pore space / m^3 bulk
  PetscReal :: liquid_saturation      ! m^3 water / m^3 pore space
  PetscReal :: L_water                ! L water
  PetscReal :: temperature
  PetscReal :: rho_f
  PetscReal :: g
  PetscReal :: viscosity

  PetscReal :: Vaq  ! mol/L
  PetscReal :: Vim  ! mol/m^3
  PetscReal :: Rate
  PetscReal :: RateAtt, RateDet  ! mol/sec
  PetscReal :: RateDecayAq, RateDecayIm !Check units
  PetscReal :: stoichVaq
  PetscReal :: stoichVim
  PetscReal :: decayAq
  PetscReal :: decayIm

  ! Filtration model parameters for attachment
  PetscReal :: katt
  PetscReal :: dc,dp,Hamaker,rho_p, alpha, qMag, diffusionCoeff, kB
  PetscReal :: Gm, Gm5, Happel
  PetscReal :: NR, NPe, NvdW, Ngr
  PetscReal :: Eta_D, Eta_I, Eta_G, Eta_0
  PetscReal :: zero

  ! Detachment rate
  PetscReal :: kdet

  PetscInt :: idof_vaq, idof_vim

  ! Global stuff (Check global_aux.F90)
  volume = material_auxvar%volume
  porosity = material_auxvar%porosity
  liquid_saturation = global_auxvar%sat(iphase)
  L_water = porosity*liquid_saturation*volume*1.d3
  ! 1.d3 converts m^3 water -> L water

  zero = 0.d0
  viscosity = 0.0008891 !Pa-s
  if (this%viscosity_id > 0) then
    viscosity = global_auxvar%parameters(this%viscosity_id)
  endif

  temperature = global_auxvar%temp
  rho_f = global_auxvar%den_kg(iphase)
  g = EARTH_GRAVITY

! Assign concentrations of Vaq and Vim
  Vaq = rt_auxvar%total(this%species_Vaq_id,iphase)
  Vim = rt_auxvar%immobile(this%species_Vim_id)

  ! initialize all rates to zero
  Rate = 0.d0
  RateAtt = 0.d0
  RateDet = 0.d0
  RateDecayAq = 0.d0
  RateDecayIm = 0.d0

  ! stoichiometries
  ! reactants have negative stoichiometry
  ! products have positive stoichiometry
  stoichVaq = -1.d0
  stoichVim = 1.d0

  !!!!!!!!!!!!!!!!!!!
  ! Decay rate - Aqueous phase
  !
  !  Check Guillier et al. for model equation (2020)
  !  decay = ln(10)/D
  !  logD = logDref - [(T-Tref)/zT]^n
  !
  !!!!!!!!!!!!!!!!!!!
  select case(this%aqueous_decay_model)
    case(MODEL_CONSTANT)
      decayAq = this%aqueous_decay_rate
    case(MODEL_TEMPERATURE)
      decayAq = (2.302585d0/(10.d0 ** (this%aqueous_logDref - &
                                    (((temperature - this%aqueous_Tref)/ &
                                      this%aqueous_zT)**this%aqueous_n))))/3600.d0  ! 1/s
    case default
      decayAq = 0.d0
  end select

!!!!!!!!!!!!!!!!!!!
! Decay rate - Immobile phase
!!!!!!!!!!!!!!!!!!!

  select case(this%adsorbed_decay_model)
    case(MODEL_CONSTANT)
      decayIm = this%adsorbed_decay_rate
    case(MODEL_TEMPERATURE)
      decayIm = (2.302585d0/(10.d0 ** (this%adsorbed_logDref - &
                                    (((temperature - this%adsorbed_Tref)/ &
                                      this%adsorbed_zT)**this%adsorbed_n))))/3600.d0  ! 1/s
    case default
      decayIm = 0.d0
  end select

  select case(this%attachment_model)
    case(MODEL_CONSTANT)
      katt = this%attachment_rate_constant
    case(MODEL_FILTRATION)
      ! Calculate attachment rate using filtration theory
      ! See Tufenkji & Elimelech 2014 DOI : 10.1021/es034049r
      ! and Saavedra et al. 2020 DOI : 10.1016/j.jconhyd.2020.103565

      kB = this%BOLTZMANN_CONSTANT
      if (mean_soil_grain_size_index > 0) then
        dc = material_auxvar%soil_properties(mean_soil_grain_size_index)
      else
        dc = this%collector_diameter
      endif
      dp = this%particle_diameter
      rho_p = this%particle_density
      alpha = this%attachment_efficiency
      Hamaker = this%hamaker_constant

      qMag = MAX(global_auxvar%darcy_vel(iphase),1.0d-20)

      diffusionCoeff = this%BOLTZMANN_CONSTANT*(temperature+T273K) / &
                    (3.0 * PI * viscosity * dp)

      ! Non-dimensional parameters
      !! Happel parameter As
      Gm = (1.0 - porosity)**(1./3.)
      Gm5 = Gm*Gm*Gm*Gm*Gm
      Happel = (2.0 * (1.0 - Gm5)) / (2.0 - (3.0*Gm) + (3.0*Gm5) - (2.0*Gm*Gm5))

      !! Aspect ratio
      NR = dp/dc

      !! Péclet number
      NPe = (qMag * dc) / (diffusionCoeff)

      !! van der Waals number
      NvdW = Hamaker/(kB*(temperature + T273K))

      !! Gravitational number
      NGr = PI/12.0 * (dp*dp*dp*dp) * (rho_p - rho_f) * g /&
            (kB*(temperature + T273K))

      ! Collector efficiencies
      ! ( see Tufenkji & Elimelech 2014
      !   DOI : 10.1021/es034049r )
      !! Transport by diffusion
      Eta_D = 2.4 &
              * Happel**(1./3.) &
              * NR**(-0.081) &
              * NPe**(-0.715) &
              * NvdW**(0.052)

      !! Transport by interception
      Eta_I = 0.55 &
              * Happel &
              * NR**(1.55) &
              * NPe**(-0.125) &
              * NvdW**(0.125)

      !! Transport due to gravity
      Eta_G = 0.475 &
              * NR**(-1.35) &
              * NPe**(-1.11) &
              * NvdW**(0.053) &
              * NGr**(1.11)

      !! Single collector efficiency
      Eta_0 = Eta_D + Eta_I + Eta_G

      ! Rate of attachment according to CFT
      katt = 1.5 * (1 - porosity) * qMag * alpha * Eta_0 &
            / (dc * porosity)

      ! Edwin debugging
      if(this%debug_option) then
        print '(3x,"porosity = ", ES12.4)', porosity
        print '(3x,"temp C   = ", ES12.4)', Temperature
        print '(3x,"viscosit = ", ES12.4)', viscosity
        print '(3x,"densityF = ", ES12.4)', rho_f
        print '(3x,"densityP = ", ES12.4)', rho_p

        print '(3x,"DarcyqMa = ", ES12.4)', qMag
        print '(3x,"diffCoef = ", ES12.4)', diffusionCoeff
        print '(3x,"Gm       = ", ES12.4)', Gm
        print '(3x,"Gm5      = ", ES12.4)', Gm5
        print '(3x,"Happel   = ", ES12.4)', Happel
        print '(3x,"NR       = ", ES12.4)', NR
        print '(3x,"NPe      = ", ES12.4)', NPe
        print '(3x,"NvW      = ", ES12.4)', NvdW
        print '(3x,"NGr      = ", ES12.4)', NGr
        print '(3x,"EtaD     = ", ES12.4)', Eta_D
        print '(3x,"EtaI     = ", ES12.4)', Eta_I
        print '(3x,"EtaG     = ", ES12.4)', Eta_G
        print '(3x,"Eta0     = ", ES12.4)', Eta_0
        print '(3x,"katt     = ", ES12.4)', katt
        print *, "--------------------"
      endif
    case default
      katt = 0.d0
  end select

  select case(this%detachment_model)
    case(MODEL_CONSTANT)
      kdet = this%detachment_rate_constant
    case default
      kdet = 0.d0
  end select

  RateAtt = 0.0
  RateDet = 0.0
  RateDecayAq = 0.0
  RateDecayIm = 0.0

  ! Build here for attachment/detachment
  ! first-order forward - reverse (A <-> C)
  Rate = katt * Vaq * L_water - kdet * Vim * volume
  RateAtt = stoichVaq * Rate
  RateDet = stoichVim * Rate

  ! Build here for inactivation reactions
  ! first-order (A -> X)
  Rate = decayAq * Vaq * L_water
  RateDecayAq = - Rate

  Rate = decayIm * Vim * volume
  RateDecayIm = - Rate

! This awful block just tries to
! avoid concentrations below 1E-50
! (Is this avoided with TRUNCATE_CONCENTRATION ?) > sure it does
  ! if ( Vaq > 0.0 ) then
  !   if ( Vim > 0.0 ) then
  !     !Do nothing
  !   else if ( Vim <= 0.0 ) then
  !     Vim = 1.0d-50
  !     RateDet = 0.0
  !     RateDecayIm = 0.0
  !   end if
  ! else if ( Vaq <= 0.0 ) then
  !   if ( Vim > 0.0 ) then
  !     Vaq = 1.0d-50
  !     RateAtt = 0.0
  !     RateDecayAq = 0.0
  !   else if ( Vim <= 0.0 ) then
  !     Vim = 1.0d-50
  !     Vaq = 1.0d-50
  !     RateAtt = 0.0
  !     RateDet = 0.0
  !     RateDecayAq = 0.0
  !     RateDecayIm = 0.0
  !   end if
  ! end if

  ! The actual calculation:

  idof_vaq = this%species_Vaq_id
  idof_vim = this%species_Vim_id + reaction%offset_immobile

  Residual(idof_vaq) = Residual(idof_vaq) - RateAtt - RateDecayAq
  Residual(idof_vim) = Residual(idof_vim) - RateDet - RateDecayIm

  if (compute_derivative) then

    ! Residual(idof_vaq) -= RateAtt + RateDecayAq
    ! RateAtt  = -katt*Vaq*L_water + kdet*Vim*volume
    ! RateDecayAq = -decayAq*Vaq*L_water
    ! units = (mol/sec)*(kg water/mol) = kg water/sec
    Jacobian(idof_vaq,idof_vaq) = Jacobian(idof_vaq,idof_vaq) - &
      (-katt*L_water - decayAq*L_water) * &
      rt_auxvar%aqueous%dtotal(this%species_Vaq_id, &
                               this%species_Vaq_id,iphase)

    ! units = (mol/sec)*(m^3 bulk/mol) = m^3 bulk/sec
    Jacobian(idof_vaq,idof_vim) = Jacobian(idof_vaq,idof_vim) - &
      kdet*volume

    ! Residual(idof_vim) -= RateDet + RateDecayIm
    ! RateDet     =  katt*Vaq*L_water - kdet*Vim*volume
    ! RateDecayIm = -decayIm*Vim*volume
    Jacobian(idof_vim,idof_vaq) = Jacobian(idof_vim,idof_vaq) - &
      katt*L_water * &
      rt_auxvar%aqueous%dtotal(this%species_Vaq_id, &
                               this%species_Vaq_id,iphase)

    Jacobian(idof_vim,idof_vim) = Jacobian(idof_vim,idof_vim) - &
      (-kdet*volume - decayIm*volume)

  endif

  ! NOTES
  ! 1. Always subtract contribution from residual
  ! 2. Units of residual are moles/second

end subroutine AmendmentReact

! ************************************************************************** !
subroutine AmendmentDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this
  ! module
  !
  ! Author: Edwin S
  ! Date: 10/01/2020
  !

  implicit none

  class(reaction_sandbox_amendment_type) :: this

end subroutine AmendmentDestroy

end module Reaction_Sandbox_Amendment_class
