module Output_Conservation_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Output_Common_module

  implicit none

  private

  ! flags signifying the first time a routine is called during a given
  ! simulation
  PetscBool :: conservation_first

  public :: OutputConservationInit, &
            OutputConservationMapCouplers, &
            OutputConservation

contains

! ************************************************************************** !

subroutine OutputConservationInit(num_steps)
  !
  ! Initializes module variables for mass and energy conservation
  !
  ! Author: Glenn Hammond
  ! Date: 12/02/25
  !
  use Option_module

  implicit none

  PetscInt :: num_steps

  if (num_steps == 0) then
    conservation_first = PETSC_TRUE
  else
    conservation_first = PETSC_FALSE
  endif

end subroutine OutputConservationInit

! ************************************************************************** !

subroutine OutputConservationMapCouplers(realization_base)
  !
  ! Creates an integral flux for each coupler and maps the flux indices to
  ! the coupler connections
  !
  ! Author: Glenn Hammond
  ! Date: 12/03/25
  !
  use Coupler_module
  use Integral_Flux_module
  use Option_module
  use Output_Aux_module
  use Patch_module
  use Realization_Base_class

  implicit none

  class(realization_base_type) :: realization_base

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  type(coupler_type), pointer :: cur_coupler
  type(integral_flux_list_type), pointer :: integral_flux_list
  type(integral_flux_type), pointer :: new_integral_flux
  PetscInt, pointer :: connections(:)
  PetscInt :: num_connections
  PetscInt :: connection_itype
  PetscInt :: i

  PetscInt, parameter :: boundary_itype = 1
  PetscInt, parameter :: srcsink_itype = 2
  PetscInt, parameter :: prescribed_boundary_itype = 3

  patch => realization_base%patch
  option => realization_base%option
  output_option => realization_base%output_option

  if (.not.realization_base%output_option%print_conservation) return

  integral_flux_list => realization_base%patch%conservation_integral_flux_list

  do connection_itype = boundary_itype, prescribed_boundary_itype
    nullify(cur_coupler)
    select case(connection_itype)
      case(boundary_itype)
        cur_coupler => patch%boundary_condition_list%first
      case(srcsink_itype)
        cur_coupler => patch%source_sink_list%first
      case(prescribed_boundary_itype)
       ! this loop currently skips prescribed conditions as their mapping
       ! occurs earlier in RealizGetPrescribedConnections
    end select
    do
      if (.not.associated(cur_coupler)) exit
      new_integral_flux => IntegralFluxCreate()
      new_integral_flux%name = cur_coupler%name
      call IntegralFluxSizeStorage(new_integral_flux,option)
      num_connections = cur_coupler%connection_set%num_connections
      nullify(connections)
      if (num_connections > 0) then
        allocate(connections(num_connections))
        connections(:) = [(i,i=1,num_connections)]
        connections = connections + cur_coupler%connection_set%offset
      endif
      select case(connection_itype)
        case(boundary_itype)
          new_integral_flux%boundary_connections => connections
        case(srcsink_itype)
          new_integral_flux%source_sink_connections => connections
        case default
          option%io_buffer = 'Unmapped connections in &
            &OutputConservationMapCouplers.'
          call PrintErrMsg(option)
      end select
      call IntegralFluxAddToList(new_integral_flux,integral_flux_list)
      nullify(new_integral_flux)
      cur_coupler => cur_coupler%next
    enddo
  enddo

end subroutine OutputConservationMapCouplers

! ************************************************************************** !

subroutine OutputConservation(realization_base)
  !
  ! Write mass and energy conservation data in column delineated format
  !
  ! Author: Glenn Hammond
  ! Date: 12/02/25
  !
  use Realization_Subsurface_class, only : realization_subsurface_type
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Utility_module
  use Output_Aux_module
  use String_module
  use Integral_Flux_module

  use Richards_module, only : RichardsComputeMassBalance
  use Mphase_module, only : MphaseComputeMassBalance
  use TH_module, only : THComputeMassBalance
  use THC_module, only : THCComputeConservation
  use Reactive_Transport_module, only : RTComputeMassBalance
  use NW_Transport_module, only : NWTComputeMassBalance
  use General_module, only : GeneralComputeMassBalance
  use Hydrate_module, only : HydrateComputeMassBalance
  use WIPP_Flow_module, only : WIPPFloComputeMassBalance
  use ZFlow_module, only : ZFlowComputeMassBalance
  use SCO2_module, only : SCO2ComputeMassBalance, &
                          SCO2ComputeComponentMassBalance
  use SCO2_Aux_module, only : SCO2_ENERGY_EQUATION_INDEX, &
                              SCO2_SALT_EQUATION_INDEX, &
                              sco2_thermal

  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use Reaction_Aux_module
  use NW_Transport_Aux_module
  use Material_Aux_module

  implicit none

  class(realization_base_type), target :: realization_base

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(output_option_type), pointer :: output_option

  class(reaction_rt_type), pointer :: reaction
  class(reaction_nw_type), pointer :: reaction_nw

  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: fid = 86
  PetscInt :: i,icol
  PetscReal :: sum_eq(20)
  PetscReal :: sum_eq_global(20)
  PetscReal :: sum_kg(realization_base%option%nflowspec, &
                      realization_base%option%nphase)
  PetscReal :: sum_kg_global(realization_base%option%nflowspec, &
                             realization_base%option%nphase)
  PetscReal :: sum_trapped(realization_base%option%nphase)
  PetscReal :: sum_trapped_global(realization_base%option%nphase)
  PetscReal :: energy_balance
  PetscReal :: energy_balance_global

  type(integral_flux_type), pointer :: integral_flux

  PetscReal, allocatable :: array(:,:)
  PetscReal, allocatable :: array_global(:,:)
  PetscReal, allocatable :: instantaneous_array(:)
  PetscReal :: flow_dof_scale(10)
  PetscReal :: tempreal
  PetscInt :: istart, iend
  PetscInt :: icomp, iphase, ispec
  PetscInt :: max_tran_size
  PetscInt :: j
  PetscInt :: nflowdof_print
  PetscInt :: nvalues
  PetscInt :: offset
  PetscInt :: cell_ids(realization_base%patch%grid%nlmax)

  PetscMPIInt :: int_mpi
  PetscErrorCode :: ierr

  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option
  reaction => ReactionAuxCast(realization_base%reaction_base)
  reaction_nw => NWTReactionCast(realization_base%reaction_base)
  output_option => realization_base%output_option

  select case(option%iflowmode)
    case(NULL_MODE)
    case(RICHARDS_MODE,TH_MODE,THC_MODE,SCO2_MODE,G_MODE)
    case default
      option%io_buffer = 'Mass and energy conservation must be set up &
        &for the "' // trim(option%flowmode) // '" flow mode.'
      call PrintErrMsg(option)
  end select
  select case(option%itranmode)
    case(NULL_MODE)
    case(RT_MODE)
    case default
      option%io_buffer = 'Mass conservation must be set up for the "' // &
      trim(option%tranmode) // '" transport mode.'
      call PrintErrMsg(option)
  end select

  call OutputCommonMapFlowFormulaWeight(option,flow_dof_scale)

  if (len_trim(output_option%plot_name) > 2) then
    filename = trim(output_option%plot_name) // '-con.dat'
  else
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-con.dat'
  endif

  ! open file
  if (OptionIsIORank(option)) then

    if (conservation_first .or. .not.FileExists(filename)) then

      if (output_option%print_column_ids) then
        icol = 1
      else
        icol = -1
      endif

      open(unit=fid,file=filename,action="write",status="replace")

      call OutputCommonGlobalMassHeader(realization_base,fid,icol,PETSC_TRUE)

      integral_flux => patch%conservation_integral_flux_list%first
      do
        if (.not.associated(integral_flux)) exit
        call OutputCommonFluxHeader(realization_base, &
                                    integral_flux%name,fid,icol,PETSC_TRUE)
        integral_flux => integral_flux%next
      enddo

      call WriteNewLine(fid)
    else
      open(unit=fid,file=filename,action="write",status="old",position="append")
    endif

  endif

  ! write time
  if (OptionIsIORank(option)) then
    call WriteRealFixedNoAdv(fid,option%time/output_option%tconv)
  endif

  if (option%nflowdof > 0) then
    if (OptionIsIORank(option)) then
      call WriteRealFixedNoAdv(fid,option%flow_dt/output_option%tconv)
    endif
  endif
  if (option%ntrandof > 0) then
    if (OptionIsIORank(option)) then
      call WriteRealFixedNoAdv(fid,option%tran_dt/output_option%tconv)
    endif
  endif

! print out global mass and energy balance

  if (option%nflowdof > 0) then
    sum_eq = 0.d0
    sum_kg = 0.d0
    sum_trapped = 0.d0
    energy_balance = 0.d0
    select type(realization_base)
      class is(realization_subsurface_type)
        select case(option%iflowmode)
          case(RICHARDS_MODE,RICHARDS_TS_MODE,PNF_MODE)
            call RichardsComputeMassBalance(realization_base,sum_eq(1))
          case(ZFLOW_MODE)
!            call ZFlowComputeMassBalance(realization_base,sum_kg(1,:))
          case(THC_MODE)
            call THCComputeConservation(realization_base,sum_eq(1), &
                                        sum_eq(2),sum_eq(3))
          case(TH_MODE,TH_TS_MODE)
            call THComputeMassBalance(realization_base,sum_eq(1),sum_eq(2))
          case(MPH_MODE)
!            call MphaseComputeMassBalance(realization_base,sum_kg(:,:), &
!                                          sum_trapped(:))
         case(G_MODE)
            cell_ids = (/ (i, i=1, realization_base%patch%grid%nlmax) /)
            call GeneralComputeMassBalance(realization_base, &
                                           realization_base%patch%grid%nlmax, &
                                           cell_ids,sum_kg(:,:), &
                                           energy_balance)
          case(H_MODE)
!            call HydrateComputeMassBalance(realization_base,sum_kg(:,:))
          case(WF_MODE)
!            call WIPPFloComputeMassBalance(realization_base,sum_kg(:,1))
          case(SCO2_MODE)
            call SCO2ComputeMassBalance(realization_base,sum_kg(:,:), &
                                        sum_trapped(:),energy_balance)
        end select
      class default
        option%io_buffer = 'Unrecognized realization class in &
                           &OutputConservation().'
        call PrintErrMsg(option)
    end select

    select case(option%iflowmode)
      case(RICHARDS_MODE,RICHARDS_TS_MODE,TH_MODE,TH_TS_MODE,THC_MODE)
        int_mpi = option%nflowdof
        call MPI_Reduce(sum_eq,sum_eq_global,int_mpi,MPI_DOUBLE_PRECISION, &
                        MPI_SUM,option%comm%io_rank,option%mycomm, &
                        ierr);CHKERRQ(ierr)
      case(G_MODE)
        ! pack: sum_kg then energy
        nvalues = option%nflowspec*option%nphase
        sum_eq(1:nvalues) = reshape(sum_kg,(/nvalues/))
        nvalues = nvalues + 1
        sum_eq(nvalues) = energy_balance
        int_mpi = nvalues
        call MPI_Reduce(sum_eq,sum_eq_global,int_mpi,MPI_DOUBLE_PRECISION, &
                        MPI_SUM,option%comm%io_rank,option%mycomm, &
                        ierr);CHKERRQ(ierr)
        offset = option%nflowspec*option%nphase
        sum_kg_global = reshape(sum_eq_global(1:offset), &
                                (/option%nflowspec,option%nphase/))
        energy_balance_global = sum_eq_global(nvalues)
      case(SCO2_MODE)
        ! pack: sum_kg, sum_trapped, energy
        nvalues = option%nflowspec*option%nphase
        sum_eq(1:nvalues) = reshape(sum_kg,(/nvalues/))
        offset = nvalues
        nvalues = nvalues + option%nphase
        sum_eq(offset+1:nvalues) = sum_trapped(1:option%nphase)
        nvalues = nvalues + 1
        sum_eq(nvalues) = energy_balance
        int_mpi = nvalues
        call MPI_Reduce(sum_eq,sum_eq_global,int_mpi,MPI_DOUBLE_PRECISION, &
                        MPI_SUM,option%comm%io_rank,option%mycomm, &
                        ierr);CHKERRQ(ierr)
        offset = option%nflowspec*option%nphase
        sum_kg_global = reshape(sum_eq_global(1:offset), &
                                (/option%nflowspec,option%nphase/))
        sum_trapped_global(1:option%nphase) = &
          sum_eq_global(offset+1:offset+option%nphase)
        energy_balance_global = sum_eq_global(nvalues)
    end select

    if (OptionIsIORank(option)) then
      select case(option%iflowmode)
        case(RICHARDS_MODE,RICHARDS_TS_MODE)
          call WriteRealNoAdv(fid,sum_eq_global(1:option%nflowdof))
        case(TH_MODE,TH_TS_MODE,THC_MODE)
          call WriteRealNoAdv(fid,sum_eq_global(1:option%nflowdof))
        case(G_MODE)
          do iphase = 1, option%nphase
            do ispec = 1, option%nflowspec
              if ((iphase == 1) .or. &
                  (iphase == 2 .and. ispec <= 2) .or. &
                  (iphase == 3 .and. ispec == 3)) then
                call WriteRealNoAdv(fid,sum_kg_global(ispec,iphase))
              endif
            enddo
          enddo
          call WriteRealNoAdv(fid,energy_balance_global)
        case(SCO2_MODE)
          do iphase = 1, option%nphase
            do ispec = 1, option%nflowspec
              if (iphase == 1) then
                call WriteRealNoAdv(fid,sum_kg_global(ispec,iphase))
              elseif (iphase == 2 .and. ispec < 3) then
                call WriteRealNoAdv(fid,sum_kg_global(ispec,iphase))
              endif
            enddo
          enddo
          call WriteRealNoAdv(fid,sum_trapped_global(TWO_INTEGER))
          if (sco2_thermal) then
            call WriteRealNoAdv(fid,energy_balance_global)
          endif
        case default
          option%io_buffer = 'Select case statement in OutputConservation &
            &needs to be extended for current flow mode: ' // &
            option%flowmode
          call PrintErrMsg(option)
      end select
    endif
  endif

  if (option%ntrandof > 0) then
     select case(option%itranmode)
       case(RT_MODE)
         if (option%transport%nphase > 1) then
         !TODO(geh): Within RTComputeMassBalance() all the mass is lumped into the
         !           liquid phase.  This need to be split out.  Also, the below
         !           where mass is summed across minerals needs to be moved
         !           to reactive_transport.F90.
           option%io_buffer = 'OutputConservation() needs to be refactored to &
             &consider species in the gas phase.'
     !      call PrintErrMsg(option)
        endif
        max_tran_size = max(reaction%naqcomp,reaction%mineral%nkinmnrl, &
                          reaction%immobile%nimmobile,reaction%gas%nactive_gas)
        ! see RTComputeMassBalance for indexing used below
        allocate(array(max_tran_size,8))
        allocate(array_global(max_tran_size,8))
        array = 0.d0
        select type(realization_base)
          class is(realization_subsurface_type)
            call RTComputeMassBalance(realization_base, &
                                      realization_base%patch%grid%nlmax, &
                                      max_tran_size,array)
          class default
            option%io_buffer = 'Unrecognized realization class in MassBalance().'
            call PrintErrMsg(option)
        end select
        int_mpi = max_tran_size*8
        call MPI_Reduce(array,array_global,int_mpi,MPI_DOUBLE_PRECISION, &
                        MPI_SUM,option%comm%io_rank,option%mycomm, &
                        ierr);CHKERRQ(ierr)

        if (OptionIsIORank(option)) then
          ! total across all phases
          do icomp = 1, reaction%naqcomp
            if (reaction%primary_species_print(icomp)) then
              call WriteRealNoAdv(fid,array_global(icomp,1))
            endif
          enddo
          ! immobile species
          do i = 1, reaction%immobile%nimmobile
            if (reaction%immobile%print_me(i)) then
              call WriteRealNoAdv(fid,array_global(i,7))
            endif
          enddo
          ! gas species
          do i = 1, reaction%gas%nactive_gas
            if (reaction%gas%active_print_me(i)) then
              call WriteRealNoAdv(fid,array_global(i,8))
            endif
          enddo
        endif
    !   print out mineral contribution to mass balance
        if (option%mass_bal_detailed) then
          if (OptionIsIORank(option)) then
            do i = 1, reaction%mineral%nkinmnrl
              if (reaction%mineral%kinmnrl_print(i)) then
                call WriteRealNoAdv(fid,array_global(i,6))
              endif
            enddo
          endif
        endif
        deallocate(array,array_global)
      case(NWT_MODE)
    end select
  endif

  allocate(array(option%nflowdof + option%ntrandof,2))
  allocate(array_global(option%nflowdof + option%ntrandof,2))
  allocate(instantaneous_array(max(option%nflowdof,option%ntrandof)))
  integral_flux => patch%conservation_integral_flux_list%first
  do
    if (.not.associated(integral_flux)) exit
    array = 0.d0
    array_global = 0.d0
    if (option%nflowdof > 0) then
      istart = 1
      iend = option%nflowdof
      instantaneous_array = 0.d0
      call IntegralFluxGetInstantaneous(integral_flux, &
                                        patch%internal_flow_fluxes, &
                                        patch%boundary_flow_fluxes, &
                                        patch%ss_flow_fluxes, &
                                        option%nflowdof, &
                                        instantaneous_array,option)
      array(istart:iend,1) = &
        integral_flux%unscaled_integral_value(istart:iend)
      array(istart:iend,2) = &
        instantaneous_array(1:option%nflowdof)
    endif
    if (option%ntrandof > 0) then
      istart = option%nflowdof+1
      iend = option%nflowdof+option%ntrandof
      instantaneous_array = 0.d0
      call IntegralFluxGetInstantaneous(integral_flux, &
                                        patch%internal_tran_fluxes, &
                                        patch%boundary_tran_fluxes, &
                                        patch%ss_tran_fluxes, &
                                        option%ntrandof, &
                                        instantaneous_array,option)
      array(istart:iend,1) = &
        integral_flux%unscaled_integral_value(istart:iend)
      array(istart:iend,2) = &
        instantaneous_array(1:option%ntrandof)
    endif
    int_mpi = size(array)
    call MPI_Reduce(array,array_global,int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                    option%comm%io_rank,option%mycomm,ierr);CHKERRQ(ierr)
    ! time units conversion
    array_global(:,2) = array_global(:,2) * output_option%tconv
    if (OptionIsIORank(option)) then
      if (option%nflowdof > 0) then
        ! for SCO2, skip well DOF (beyond energy); residual order is
        ! water, CO2, salt, energy, [well]
        nflowdof_print = option%nflowdof
        select case(option%iflowmode)
          case(SCO2_MODE)
            nflowdof_print = SCO2_SALT_EQUATION_INDEX
            if (sco2_thermal) then
              nflowdof_print = SCO2_ENERGY_EQUATION_INDEX
            endif
        end select
        do i = 1, nflowdof_print
          do j = 1, 2  ! 1 = integral, 2 = instantaneous
            tempreal = array_global(i,j)*flow_dof_scale(i)
            call WriteRealNoAdv(fid,tempreal)
          enddo
        enddo
      endif
      if (option%ntrandof > 0) then
        istart = option%nflowdof
        select case(option%itranmode)
          case(RT_MODE)
            do i=1,reaction%naqcomp
              do j = 1, 2  ! 1 = integral, 2 = instantaneous
                if (reaction%primary_species_print(i)) then
                  tempreal = array_global(istart+i,j)
                  call WriteRealNoAdv(fid,tempreal)
                endif
              enddo
            enddo
          case(NWT_MODE)
            do i=1,reaction_nw%params%nspecies
              do j = 1, 2  ! 1 = integral, 2 = instantaneous
                if (reaction_nw%species_print(i)) then
                  tempreal = array_global(istart+i,j)
                  call WriteRealNoAdv(fid,tempreal)
                endif
              enddo
            enddo
        end select
      endif
    endif
    integral_flux => integral_flux%next
  enddo
  deallocate(array)
  deallocate(array_global)
  deallocate(instantaneous_array)

  if (OptionIsIORank(option)) then
    call WriteNewLine(fid)
    close(fid)
  endif

  conservation_first = PETSC_FALSE

end subroutine OutputConservation

end module Output_Conservation_module
