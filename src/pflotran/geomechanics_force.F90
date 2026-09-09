module Geomechanics_Force_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Geomechanics_Auxiliary_module
  use Geomechanics_Global_Aux_module
  use Geomechanics_Linear_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-12
  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public :: GeomechForceSetup, &
            GeomechForceUpdateAuxVars, &
            GeomechanicsForceInitialGuess, &
            GeomechUpdateFromSubsurf, &
            GeomechUpdateSubsurfFromGeomech, &
            GeomechCreateGeomechSubsurfVec, &
            GeomechCreateSubsurfStressStrainVec, &
            GeomechUpdateSolution, &
            GeomechStoreInitialPressTemp, &
            GeomechStoreInitialDisp, &
            GeomechStoreInitialPorosity, &
            GeomechForceSetupLinearSystem, &
            GeomechForceAssembleCoeffMatrix, &
            GeomechMapMaterialIdsToFlow

contains

! ************************************************************************** !

subroutine GeomechForceSetup(geomech_realization)
  !
  ! Sets up the geomechanics calculations
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  !

  use Geomechanics_Realization_class
  use Output_Aux_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  type(output_variable_list_type), pointer :: list

  call GeomechForceSetupPatch(geomech_realization)

  list => geomech_realization%output_option%output_snap_variable_list
  call GeomechForceSetPlotVariables(list)
  list => geomech_realization%output_option%output_obs_variable_list
  call GeomechForceSetPlotVariables(list)

end subroutine GeomechForceSetup

! ************************************************************************** !

subroutine GeomechForceSetupPatch(geomech_realization)
  !
  ! Sets up the arrays for geomech parameters
  !
  ! Author: Satish Karra, LANL
  ! Date: 09/11/13
  !
  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Option_module
  use Parameter_module

  implicit none

  class(realization_geomech_type) :: geomech_realization

  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch
  type(geomech_linear_parameter_type), pointer :: geomech_parameter

  PetscInt :: i

  option => geomech_realization%option
  patch => geomech_realization%geomech_patch

  patch%geomech_aux%Linear => GeomechLinearAuxCreate()

  geomech_parameter => patch%geomech_aux%Linear%linear_parameter

  geomech_parameter%youngs_modulus_spatially_varying = PETSC_FALSE
  geomech_parameter%poissons_ratio_spatially_varying = PETSC_FALSE
  geomech_parameter%density_spatially_varying = PETSC_FALSE
  geomech_parameter%biot_coeff_spatially_varying = PETSC_FALSE
  geomech_parameter%thermal_exp_coeff_spatially_varying = PETSC_FALSE
  do i = 1, size(geomech_realization%geomech_material_property_array)
    geomech_parameter%youngs_modulus_spatially_varying = &
      geomech_parameter%youngs_modulus_spatially_varying .or. &
      associated(geomech_realization%geomech_material_property_array(i)%ptr%youngs_modulus_dataset)
    geomech_parameter%poissons_ratio_spatially_varying = &
      geomech_parameter%poissons_ratio_spatially_varying .or. &
      associated(geomech_realization%geomech_material_property_array(i)%ptr%poissons_ratio_dataset)
    geomech_parameter%density_spatially_varying = &
      geomech_parameter%density_spatially_varying .or. &
      associated(geomech_realization%geomech_material_property_array(i)%ptr%density_dataset)
    geomech_parameter%biot_coeff_spatially_varying = &
      geomech_parameter%biot_coeff_spatially_varying .or. &
      associated(geomech_realization%geomech_material_property_array(i)%ptr%biot_coeff_dataset)
    geomech_parameter%thermal_exp_coeff_spatially_varying = &
      geomech_parameter%thermal_exp_coeff_spatially_varying .or. &
      associated(geomech_realization%geomech_material_property_array(i)%ptr%thermal_exp_coeff_dataset)
  enddo

  ! Young's modulus
  if (geomech_parameter%youngs_modulus_spatially_varying) then
  else
    allocate(geomech_parameter%youngs_modulus &
      (size(geomech_realization%geomech_material_property_array)))
    do i = 1, size(geomech_realization%geomech_material_property_array)
      geomech_parameter%youngs_modulus(geomech_realization% &
        geomech_material_property_array(i)%ptr%id) = geomech_realization% &
        geomech_material_property_array(i)%ptr%youngs_modulus
    enddo
  endif
  ! Poisson's ratio
  if (geomech_parameter%poissons_ratio_spatially_varying) then
  else
    allocate(geomech_parameter%poissons_ratio &
      (size(geomech_realization%geomech_material_property_array)))
    do i = 1, size(geomech_realization%geomech_material_property_array)
      geomech_parameter%poissons_ratio(geomech_realization% &
        geomech_material_property_array(i)%ptr%id) = geomech_realization% &
        geomech_material_property_array(i)%ptr%poissons_ratio
    enddo
  endif
  ! Density
  if (geomech_parameter%density_spatially_varying) then
  else
    allocate(geomech_parameter%density &
      (size(geomech_realization%geomech_material_property_array)))
    do i = 1, size(geomech_realization%geomech_material_property_array)
      geomech_parameter%density(geomech_realization% &
        geomech_material_property_array(i)%ptr%id) = geomech_realization% &
        geomech_material_property_array(i)%ptr%density
    enddo
  endif
  ! Biot's coefficient
  if (geomech_parameter%biot_coeff_spatially_varying) then
  else
    allocate(geomech_parameter%biot_coeff &
      (size(geomech_realization%geomech_material_property_array)))
    do i = 1, size(geomech_realization%geomech_material_property_array)
      geomech_parameter%biot_coeff(geomech_realization% &
        geomech_material_property_array(i)%ptr%id) = geomech_realization% &
        geomech_material_property_array(i)%ptr%biot_coeff
    enddo
  endif
  ! Thermal expansion coefficient
  if (geomech_parameter%thermal_exp_coeff_spatially_varying) then
  else
    allocate(geomech_parameter%thermal_exp_coeff &
      (size(geomech_realization%geomech_material_property_array)))
    do i = 1, size(geomech_realization%geomech_material_property_array)
      geomech_parameter%thermal_exp_coeff(geomech_realization% &
        geomech_material_property_array(i)%ptr%id) = geomech_realization% &
        geomech_material_property_array(i)%ptr%thermal_exp_coeff
    enddo
  endif

  select case(option%geomechanics%flow_coupling)
    case(GEOMECH_TWO_WAY_COUPLED)
      select case(option%geomechanics%split_scheme)
        case(GEOMECH_FIXED_STRESS_SPLIT)
          geomech_parameter%press_0_id = &
            ParameterGetIDFromName('press_0',option)
          geomech_parameter%temp_0_id = &
            ParameterGetIDFromName('temp_0',option)
          geomech_parameter%vol_strain_0_id = &
            ParameterGetIDFromName('vol_strain_0',option)
          geomech_parameter%vol_strain_id = &
            ParameterGetIDFromName('vol_strain',option)
          geomech_parameter%stored_pressure_id = &
            ParameterGetIDFromName('stored_pressure',option)
          geomech_parameter%stored_porosity_id = &
            ParameterGetIDFromName('stored_porosity',option)
          geomech_parameter%flow_porosity_id = &
            ParameterGetIDFromName('flow_porosity',option)
      end select
  end select

end subroutine GeomechForceSetupPatch

! ************************************************************************** !

subroutine GeomechForceSetPlotVariables(list)
  !
  ! Set up of geomechanics plot variables
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  !

  use Output_Aux_module
  use Variables_module

  implicit none

  type(output_variable_list_type), pointer :: list
  type(output_variable_type), pointer :: output_variable

  character(len=MAXWORDLENGTH) :: name, units

  if (associated(list%first)) then
    return
  endif

  name = 'displacement_x'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_DISP_X)

  name = 'displacement_y'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_DISP_Y)

  name = 'displacement_z'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_DISP_Z)

  units = ''
  name = 'Material ID'
  output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE, &
                                          units,GEOMECH_MATERIAL_ID)
  output_variable%iformat = 1 ! integer
  call OutputVariableAddToList(list,output_variable)

  name = 'volumetric strain'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               GEOMECH_VOLUMETRIC_STRAIN)

  name = 'strain_xx'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_XX)

  name = 'strain_yy'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_YY)

  name = 'strain_zz'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_ZZ)

  name = 'strain_xy'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_XY)

  name = 'strain_yz'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_YZ)

  name = 'strain_zx'
  units = ''
  call OutputVariableAddToList(list,name,OUTPUT_STRAIN,units, &
                               STRAIN_ZX)

  name = 'stress_xx'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_XX)

  name = 'stress_yy'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_YY)

  name = 'stress_zz'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_ZZ)

  name = 'stress_xy'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_XY)

  name = 'stress_yz'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_YZ)

  name = 'stress_zx'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_ZX)

  name = 'stress_total_xx'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_TOTAL_XX)

  name = 'stress_total_yy'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_TOTAL_YY)

  name = 'stress_total_zz'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_TOTAL_ZZ)

  name = 'stress_total_xy'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_TOTAL_XY)

  name = 'stress_total_yz'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_TOTAL_YZ)

  name = 'stress_total_zx'
  units = 'Pa'
  call OutputVariableAddToList(list,name,OUTPUT_STRESS,units, &
                               STRESS_TOTAL_ZX)

  name = 'relative_displacement_x'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_REL_DISP_X)

  name = 'relative_displacement_y'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_REL_DISP_Y)

  name = 'relative_displacement_z'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_DISPLACEMENT,units, &
                               GEOMECH_REL_DISP_Z)


end subroutine GeomechForceSetPlotVariables

! ************************************************************************** !

subroutine GeomechanicsForceInitialGuess(geomech_realization)
  !
  ! Sets up the inital guess for the solution
  ! The boundary conditions are set here
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/19/13
  !

  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Option_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Geomechanics_Patch_module
  use Geomechanics_Coupler_module
  use Geomechanics_Region_module

  implicit none

  class(realization_geomech_type) :: geomech_realization

  type(option_type), pointer :: option
  type(geomech_field_type), pointer :: field
  type(geomech_patch_type), pointer :: patch
  type(geomech_coupler_type), pointer :: boundary_condition
  type(geomech_grid_type), pointer :: grid
  type(gm_region_type), pointer :: region

  PetscInt :: ghosted_id,local_id,total_verts,ivertex
  PetscReal, pointer :: xx_p(:)
  PetscErrorCode :: ierr

  option => geomech_realization%option
  field => geomech_realization%geomech_field
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid

  call VecGetArray(field%disp_xx,xx_p,ierr);CHKERRQ(ierr)

  boundary_condition => patch%geomech_boundary_condition_list%first
  total_verts = 0
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      total_verts = total_verts + 1
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      ! X displacement
      if (associated(boundary_condition%geomech_condition%displacement_x)) then
        select case(boundary_condition%geomech_condition%displacement_x%itype)
          case(DIRICHLET_BC)
            xx_p(THREE_INTEGER*(local_id-1) + GEOMECH_DISP_X_DOF) = &
            boundary_condition%geomech_aux_real_var(GEOMECH_DISP_X_DOF,ivertex)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif

      ! Y displacement
      if (associated(boundary_condition%geomech_condition%displacement_y)) then
        select case(boundary_condition%geomech_condition%displacement_y%itype)
          case(DIRICHLET_BC)
            xx_p(THREE_INTEGER*(local_id-1) + GEOMECH_DISP_Y_DOF) = &
            boundary_condition%geomech_aux_real_var(GEOMECH_DISP_Y_DOF,ivertex)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif

      ! Z displacement
      if (associated(boundary_condition%geomech_condition%displacement_z)) then
        select case(boundary_condition%geomech_condition%displacement_z%itype)
          case(DIRICHLET_BC)
            xx_p(THREE_INTEGER*(local_id-1) + GEOMECH_DISP_Z_DOF) = &
            boundary_condition%geomech_aux_real_var(GEOMECH_DISP_Z_DOF,ivertex)
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArray(field%disp_xx,xx_p,ierr);CHKERRQ(ierr)

end subroutine GeomechanicsForceInitialGuess

! ************************************************************************** !

subroutine GeomechForceUpdateAuxVars(geomech_realization)
  !
  ! Updates the geomechanics variables
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/18/13
  !

  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Option_module
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Coupler_module
  use Geomechanics_Material_module
  use Geomechanics_Global_Aux_module
  use Geomechanics_Region_module

  implicit none

  class(realization_geomech_type) :: geomech_realization

  type(option_type), pointer :: option
  type(geomech_patch_type), pointer :: patch
  type(geomech_grid_type), pointer :: grid
  type(geomech_field_type), pointer :: geomech_field
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)

  PetscInt :: ghosted_id
  PetscReal, pointer :: xx_loc_p(:), xx_init_loc_p(:)
  PetscErrorCode :: ierr

  option => geomech_realization%option
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  geomech_field => geomech_realization%geomech_field
  geomech_global_aux_vars => patch%geomech_aux%Global%aux_vars


  ! Communication -----------------------------------------
  call GeomechDiscretizationGlobalToLocal(geomech_realization% &
                                          geomech_discretization, &
                                          geomech_realization% &
                                          geomech_field%disp_xx, &
                                          geomech_realization% &
                                          geomech_field%disp_xx_loc,NGEODOF)

  call VecGetArray(geomech_field%disp_xx_loc,xx_loc_p,ierr)
  call VecGetArray(geomech_field%disp_xx_init_loc,xx_init_loc_p,ierr)

  ! Internal aux vars
  do ghosted_id = 1, grid%ngmax_node
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_X_DOF) = &
      xx_loc_p(GEOMECH_DISP_X_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_Y_DOF) = &
      xx_loc_p(GEOMECH_DISP_Y_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%disp_vector(GEOMECH_DISP_Z_DOF) = &
      xx_loc_p(GEOMECH_DISP_Z_DOF + (ghosted_id-1)*THREE_INTEGER)

    geomech_global_aux_vars(ghosted_id)%rel_disp_vector(GEOMECH_DISP_X_DOF) = &
      xx_loc_p(GEOMECH_DISP_X_DOF + (ghosted_id-1)*THREE_INTEGER) - &
      xx_init_loc_p(GEOMECH_DISP_X_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%rel_disp_vector(GEOMECH_DISP_Y_DOF) = &
      xx_loc_p(GEOMECH_DISP_Y_DOF + (ghosted_id-1)*THREE_INTEGER) - &
      xx_init_loc_p(GEOMECH_DISP_Y_DOF + (ghosted_id-1)*THREE_INTEGER)
    geomech_global_aux_vars(ghosted_id)%rel_disp_vector(GEOMECH_DISP_Z_DOF) = &
      xx_loc_p(GEOMECH_DISP_Z_DOF + (ghosted_id-1)*THREE_INTEGER) - &
      xx_init_loc_p(GEOMECH_DISP_Z_DOF + (ghosted_id-1)*THREE_INTEGER)
 enddo

  call VecRestoreArray(geomech_field%disp_xx_loc,xx_loc_p,ierr)
  call VecRestoreArray(geomech_field%disp_xx_init_loc,xx_init_loc_p,ierr)


end subroutine GeomechForceUpdateAuxVars

! ************************************************************************** !

subroutine ComputeTetVolAtVertex(vert_0, &
                                 vert_1, &
                                 vert_2, &
                                 vert_3, &
                                 volume)
  !
  ! Returns the volume of the specified tetrahedron after
  ! being clipped on the positive side of the planes computed
  ! in this function.
  ! This works by calling ClippedVolume() recursively on one plane
  ! at a time.
  !
  ! Author: Joe Eyles, WSP
  ! Date: 29/04/2025
  !


  PetscReal :: vert_0(3)
  PetscReal :: vert_1(3)
  PetscReal :: vert_2(3)
  PetscReal :: vert_3(3)

  PetscReal :: midPoints(3, 3), normals(3, 3)
  PetscReal :: volume

  midPoints(1, :) = (vert_0(:) + vert_1(:))*0.5
  midPoints(2, :) = (vert_0(:) + vert_2(:))*0.5
  midPoints(3, :) = (vert_0(:) + vert_3(:))*0.5

  normals(1, :) = vert_0(:) - vert_1(:)
  normals(2, :) = vert_0(:) - vert_2(:)
  normals(3, :) = vert_0(:) - vert_3(:)

  normals(1, :) = normals(1, :) / sqrt(dot_product(normals(1, :), normals(1, :)))
  normals(2, :) = normals(2, :) / sqrt(dot_product(normals(2, :), normals(2, :)))
  normals(3, :) = normals(3, :) / sqrt(dot_product(normals(3, :), normals(3, :)))

  volume = ClippedVolume(vert_0, vert_1, vert_2, vert_3, midPoints, normals, 3)
end subroutine ComputeTetVolAtVertex

! ************************************************************************** !

recursive function ClippedVolume(vert_0, vert_1, vert_2, vert_3, midPoints, normals, nSize_in) result(nRet)
  !
  ! Clips on the given plane, splits the resulting
  ! polyhedron into tetrahedra, and then calls itself recursively,
  ! ultimately returning the volume.
  !
  ! Author: Joe Eyles, WSP
  ! Date: 29/04/2025
  !

  PetscReal :: vert_0(3), vert_1(3), vert_2(3), vert_3(3)
  PetscReal :: midPoints(3, 3), normals(3, 3)
  PetscInt :: nSize_in
  PetscInt :: nSize
  PetscReal :: nRet
  PetscReal :: vP(3), vN(3), vvert_0(3), vvert_1(3), vvert_2(3), vvert_3(3), vvert_4(3)
  PetscReal :: V(4, 3)
  PetscReal :: d(4)

  nRet = 0
  nSize = nSize_in

  if (nSize == 0) then
    nRet = abs(dot_product(cross_product(vert_1 - vert_0, vert_2 - vert_0), vert_3 - vert_0) / 6.0)
  else
    nSize = nSize - 1
    vP = midPoints(nSize + 1, :)
    vN = normals(nSize + 1, :)
    V(1, :) = vert_0
    d(1) = dot_product(vert_0 - vP, vN)
    V(2, :) = vert_1
    d(2) = dot_product(vert_1 - vP, vN)
    V(3, :) = vert_2
    d(3) = dot_product(vert_2 - vP, vN)
    V(4, :) = vert_3
    d(4) = dot_product(vert_3 - vP, vN)

    call sort(V, d)

    if (d(1) <= 0.0) then
      nRet = 0.0
    else if (d(4) >= 0.0) then
      nRet = nRet + ClippedVolume(vert_0, vert_1, vert_2, vert_3, midPoints, normals, nSize)
    else if (d(3) > 0.0) then
      vvert_0 = Intersect(V(1, :), V(4, :), d(1), d(4))
      vvert_1 = Intersect(V(2, :), V(4, :), d(2), d(4))
      vvert_2 = Intersect(V(3, :), V(4, :), d(3), d(4))
      nRet = nRet + ClippedVolume(V(1, :), V(2, :), V(3, :), vvert_2, midPoints, normals, nSize) + &
             ClippedVolume(V(1, :), vvert_0, V(2, :), vvert_2, midPoints, normals, nSize) + &
             ClippedVolume(V(2, :), vvert_0, vvert_1, vvert_2, midPoints, normals, nSize)
    else if (d(3) == 0.0) then
      vvert_0 = Intersect(V(1, :), V(4, :), d(1), d(4))
      vvert_1 = Intersect(V(2, :), V(4, :), d(2), d(4))
      nRet = nRet + ClippedVolume(V(1, :), vvert_0, V(2, :), V(3, :), midPoints, normals, nSize) + &
             ClippedVolume(V(2, :), vvert_0, vvert_1, V(3, :), midPoints, normals, nSize)
    else if (d(2) > 0.0) then
      vvert_1 = Intersect(V(1, :), V(3, :), d(1), d(3))
      vvert_2 = Intersect(V(1, :), V(4, :), d(1), d(4))
      vvert_3 = Intersect(V(2, :), V(3, :), d(2), d(3))
      vvert_4 = Intersect(V(2, :), V(4, :), d(2), d(4))
      nRet = nRet + ClippedVolume(V(1, :), vvert_1, vvert_2, vvert_3, midPoints, normals, nSize) + &
             ClippedVolume(V(1, :), vvert_2, vvert_3, vvert_4, midPoints, normals, nSize) + &
             ClippedVolume(V(1, :), V(2, :), vvert_3, vvert_4, midPoints, normals, nSize)
    else
      vvert_1 = Intersect(V(1, :), V(2, :), d(1), d(2))
      vvert_2 = Intersect(V(1, :), V(3, :), d(1), d(3))
      vvert_3 = Intersect(V(1, :), V(4, :), d(1), d(4))
      nRet = nRet + ClippedVolume(V(1, :), vvert_1, vvert_2, vvert_3, midPoints, normals, nSize)
    end if
  end if
end function ClippedVolume

! ************************************************************************** !

function Intersect(v1, v2, d1, d2) result(v)
  !
  ! Computes the intesection
  !
  ! Author: Joe Eyles, WSP
  ! Date: 29/04/2025
  !

  PetscReal :: v1(3), v2(3)
  PetscReal :: d1, d2
  PetscReal :: v(3)

  v = v1 * (-d2 / (d1 - d2)) + v2 * (d1 / (d1 - d2))
end function Intersect

! ************************************************************************** !

function cross_product(v1, v2) result(v)
  !
  ! Computes the cross product between v1 and v2
  !
  ! Author: Joe Eyles, WSP
  ! Date: 29/04/2025
  !

  PetscReal :: v1(3), v2(3)
  PetscReal :: v(3)

  v(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v(3) = v1(1) * v2(2) - v1(2) * v2(1)
end function cross_product

! ************************************************************************** !

subroutine sort(V, d)
  !
  ! Sorts vectors V and d based on d
  !
  ! Author: Joe Eyles, WSP
  ! Date: 29/04/2025
  !

  PetscReal :: V(4, 3)
  PetscReal :: d(4)
  integer :: i, j
  PetscReal :: temp_d
  PetscReal :: temp_V(3)

  do i = 1, 3
    do j = i + 1, 4
      if (d(i) < d(j)) then
        temp_d = d(i)
        d(i) = d(j)
        d(j) = temp_d
        temp_V = V(i, :)
        V(i, :) = V(j, :)
        V(j, :) = temp_V
      end if
    end do
  end do
end subroutine sort

! ************************************************************************** !

function face_unitnormal(v1, v2, v3, option) result(n)
  !
  ! calculates the outward pointing normal vector
  ! for tri face, pass coordinates in the same order
  ! for quad face, pass (v1, v2, v4)
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 09/04/2025
  !

  use Option_module

  PetscReal :: v1(THREE_INTEGER), v2(THREE_INTEGER), v3(THREE_INTEGER)
  type(option_type), optional, intent(inout) :: option
  PetscReal :: n(THREE_INTEGER)

  PetscReal :: e1(THREE_INTEGER), e2(THREE_INTEGER)
  PetscReal :: nmag2

  e1 = v2 - v1
  e2 = v3 - v1
  n = cross_product(e1,e2)
  nmag2 = dot_product(n,n)
  if (nmag2 <= 0.d0) then
    call GeomechForceError('GEOMECHANICS: face normal undefined for degenerate face.', option)
    n = 0.d0
    return
  endif
  n = n/sqrt(nmag2)

end function face_unitnormal

! ************************************************************************** !

subroutine GeomechForceError(message, option)
  use Option_module
  implicit none

  character(len=*), intent(in) :: message
  type(option_type), optional, intent(inout) :: option

  if (present(option)) then
    option%io_buffer = message
    call PrintErrMsg(option)
  else
    print *, trim(message)
    stop
  endif
end subroutine GeomechForceError

! ************************************************************************** !

subroutine GeomechForceApplyTractionBCtoRHS(local_coordinates, &
                                            facetype, &
                                            stress_bc, &
                                            r,w,rhs_vec,option)
  !
  ! Computes the traction contribution for the linear system
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 09/05/2025
  !


  use Grid_Unstructured_Cell_module
  use Shape_Function_module
  use Option_module

  implicit none

  type(option_type) :: option
  PetscInt, intent(in) :: facetype
  PetscReal, intent(in) :: local_coordinates(:,:)  ! (nen,3)
  PetscReal, intent(in) :: stress_bc(SIX_INTEGER)
  PetscReal, pointer, intent(in) :: r(:,:), w(:)
  PetscReal, intent(inout) :: rhs_vec(:)           ! (3*nen)

  type(shapefunction_type) :: shapefunction
  PetscInt :: igpt, len_w, nen, eletype, a, ia
  PetscReal :: J_map(3,2), xp_J(3), surf_J
  PetscReal :: boundary_stress(3,3), normal_vec(3), traction(3)
  PetscReal :: wsurf

  rhs_vec = 0.0d0
  len_w = size(w)
  nen = size(local_coordinates,1)

  eletype = merge(TRI_TYPE, QUAD_TYPE, facetype == TRI_FACE_TYPE)

  ! Expand Voigt stress_bc = [sxx, syy, szz, sxy, syz, szx] into full tensor.
  boundary_stress = 0.0d0
  boundary_stress(1,1) = stress_bc(1)
  boundary_stress(2,2) = stress_bc(2)
  boundary_stress(3,3) = stress_bc(3)
  boundary_stress(1,2) = stress_bc(4); boundary_stress(2,1) = stress_bc(4) ! Symmetric shear stress
  boundary_stress(2,3) = stress_bc(5); boundary_stress(3,2) = stress_bc(5) ! Symmetric shear stress
  boundary_stress(3,1) = stress_bc(6); boundary_stress(1,3) = stress_bc(6) ! Symmetric shear stress

  if (facetype == TRI_FACE_TYPE) then
    normal_vec = face_unitnormal(local_coordinates(1,:), local_coordinates(2,:), local_coordinates(3,:), option)
  else
    normal_vec = face_unitnormal(local_coordinates(1,:), local_coordinates(2,:), local_coordinates(4,:), option)
  end if
  if (dot_product(normal_vec,normal_vec) <= 0.d0) return

  ! Traction vector t = sigma * n.
  traction = matmul(boundary_stress, normal_vec)

  shapefunction%element_type = eletype
  call ShapeFunctionInitialize(shapefunction,option)

  do igpt = 1, len_w
    shapefunction%zeta = r(igpt,:)
    call ShapeFunctionCalculate(shapefunction,option)

    ! Surface mapping Jacobian and area scaling for this quadrature point.
    J_map = matmul(transpose(local_coordinates), shapefunction%DN)  ! (3x2)
    xp_J  = cross_product(J_map(:,1), J_map(:,2))
    surf_J = sqrt(dot_product(xp_J,xp_J))
    if (surf_J <= 0.0d0) then
      call GeomechForceError('GEOMECHANICS: surface jacobian must be positive!', option)
      call ShapeFunctionDestroy(shapefunction)
      return
    end if

    ! Accumulate nodal equivalent traction load: integral(N_a * t dGamma).
    wsurf = w(igpt) * surf_J
    do a = 1, nen
      ia = 3*(a-1)
      rhs_vec(ia+1) = rhs_vec(ia+1) + wsurf * shapefunction%N(a) * traction(1)
      rhs_vec(ia+2) = rhs_vec(ia+2) + wsurf * shapefunction%N(a) * traction(2)
      rhs_vec(ia+3) = rhs_vec(ia+3) + wsurf * shapefunction%N(a) * traction(3)
    end do
  end do

  call ShapeFunctionDestroy(shapefunction)

end subroutine GeomechForceApplyTractionBCtoRHS

! ************************************************************************** !

subroutine GeomechForceSetupLinearSystem(A,solution,rhs,geomech_realization, &
                                         ierr)
  !
  ! Computes the Coefficient matrix and the right hand side of the system
  !
  ! Author: Satish Karra
  ! Date: 06/04/2025
  !

  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Geomechanics_Patch_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Grid_Unstructured_Cell_module
  use Geomechanics_Region_module
  use Geomechanics_Coupler_module
  use Option_module
  use Geomechanics_Auxiliary_module

  implicit none

  Vec :: solution
  Vec :: rhs
  Mat :: A
  class(realization_geomech_type) :: geomech_realization
  PetscErrorCode :: ierr

  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_patch_type), pointer :: patch
  type(geomech_field_type), pointer :: field
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(gm_region_type), pointer :: region
  type(geomech_coupler_type), pointer :: boundary_condition
  type(geomech_linear_parameter_type), pointer :: geomech_parameter

  PetscInt, allocatable :: elenodes(:)
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: local_press(:), local_temp(:)
  PetscInt, allocatable :: petsc_ids(:)
  PetscInt, allocatable :: ids(:)
  PetscReal, allocatable :: rhs_local_vec(:)
  PetscReal, pointer :: press_loc_p(:), temp_loc_p(:)
  PetscReal, pointer :: fluid_density_loc_p(:), porosity_loc_p(:)
  PetscReal, pointer :: press_init_loc_p(:), temp_init_loc_p(:)
  PetscReal, pointer :: fluid_density_init_loc_p(:)
  PetscReal, allocatable :: beta_vec(:), alpha_vec(:)
  PetscReal, allocatable :: density_rock_vec(:), density_fluid_vec(:)
  PetscReal, allocatable :: density_bulk_vec(:)
  PetscReal, allocatable :: youngs_vec(:), poissons_vec(:)
  PetscReal, allocatable :: porosity_vec(:)
  PetscReal, allocatable :: local_body_force(:,:)
  PetscInt :: ielem, ivertex
  PetscInt :: ghosted_id
  PetscInt :: eletype
  PetscInt :: petsc_id, local_id
  PetscReal, pointer :: imech_loc_p(:)
  PetscInt :: size_elenodes, max_elem_nodes, max_elem_dofs, ndofs, idof

  PetscInt :: facetype, nfaces
  PetscInt :: iface, num_vertices, max_face_vertices, face_ndofs
  PetscReal :: stress_bc(SIX_INTEGER)
  PetscInt, allocatable :: face_vertices(:), face_petsc_ids(:), face_ids(:)
  PetscReal, allocatable :: face_local_coordinates(:,:), face_rhs_local_vec(:)

  PetscReal, pointer :: temp_youngs_modulus_loc_p(:)
  PetscReal, pointer :: temp_poissons_ratio_loc_p(:)
  PetscReal, pointer :: temp_density_loc_p(:)
  PetscReal, pointer :: temp_biot_coeff_loc_p(:)
  PetscReal, pointer :: temp_thermal_exp_coeff_loc_p(:)
  PetscReal, pointer :: temp_body_force_x_loc_p(:)
  PetscReal, pointer :: temp_body_force_y_loc_p(:)
  PetscReal, pointer :: temp_body_force_z_loc_p(:)

  field => geomech_realization%geomech_field
  geomech_discretization => geomech_realization%geomech_discretization
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  geomech_parameter => patch%geomech_aux%Linear%linear_parameter


  solution = field%disp_xx

  ! use the stored matrix
  A = field%A
  rhs = field%rhs

  call VecZeroEntries(rhs,ierr);CHKERRQ(ierr)

  ! Get pressure and temperature from subsurface
  call VecGetArray(field%press_loc,press_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%temp_loc,temp_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%porosity_loc,porosity_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%fluid_density_loc,fluid_density_loc_p,ierr);CHKERRQ(ierr)

  ! Get geomech properties
  if (geomech_parameter%youngs_modulus_spatially_varying) then
    call VecGetArray(field%youngs_modulus_loc,temp_youngs_modulus_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%poissons_ratio_spatially_varying) then
    call VecGetArray(field%poissons_ratio_loc,temp_poissons_ratio_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%density_spatially_varying) then
    call VecGetArray(field%density_loc,temp_density_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%biot_coeff_spatially_varying) then
    call VecGetArray(field%biot_coeff_loc,temp_biot_coeff_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%thermal_exp_coeff_spatially_varying) then
    call VecGetArray(field%thermal_exp_coeff_loc,temp_thermal_exp_coeff_loc_p,ierr);CHKERRQ(ierr)
  endif
  call VecGetArray(field%body_force_x_loc,temp_body_force_x_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%body_force_y_loc,temp_body_force_y_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%body_force_z_loc,temp_body_force_z_loc_p,ierr);CHKERRQ(ierr)

  ! Get initial pressure and temperature
  call VecGetArray(field%press_init_loc,press_init_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%temp_init_loc,temp_init_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%fluid_density_init_loc,fluid_density_init_loc_p, &
                      ierr);CHKERRQ(ierr)

  max_elem_nodes = maxval(grid%elem_nodes(0,1:grid%nlmax_elem))
  max_elem_dofs = max_elem_nodes*option%ngeomechdof
  allocate(elenodes(max_elem_nodes))
  allocate(local_coordinates(max_elem_nodes,THREE_INTEGER))
  allocate(local_press(max_elem_nodes))
  allocate(local_temp(max_elem_nodes))
  allocate(petsc_ids(max_elem_nodes))
  allocate(ids(max_elem_dofs))
  allocate(rhs_local_vec(max_elem_dofs))
  allocate(beta_vec(max_elem_nodes))
  allocate(alpha_vec(max_elem_nodes))
  allocate(density_rock_vec(max_elem_nodes))
  allocate(density_fluid_vec(max_elem_nodes))
  allocate(youngs_vec(max_elem_nodes))
  allocate(poissons_vec(max_elem_nodes))
  allocate(porosity_vec(max_elem_nodes))
  allocate(density_bulk_vec(max_elem_nodes))
  allocate(local_body_force(max_elem_nodes,THREE_INTEGER))

  ! Loop over elements on a processor
  do ielem = 1, grid%nlmax_elem
    size_elenodes = grid%elem_nodes(0,ielem)
    ndofs = size_elenodes*option%ngeomechdof
    elenodes(1:size_elenodes) = grid%elem_nodes(1:size_elenodes,ielem)
    eletype = grid%gauss_node(ielem)%entity_type
    do ivertex = 1, size_elenodes
      ghosted_id = elenodes(ivertex)
      local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = grid%nodes(ghosted_id)%x
      local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = grid%nodes(ghosted_id)%y
      local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = grid%nodes(ghosted_id)%z
      petsc_ids(ivertex) = grid%node_ids_ghosted_petsc(ghosted_id)
    enddo
    do ivertex = 1, size_elenodes
      ghosted_id = elenodes(ivertex)
      do idof = 1, option%ngeomechdof
        ids(idof + (ivertex-1)*option%ngeomechdof) = &
          (petsc_ids(ivertex)-1)*option%ngeomechdof + (idof-1)
      enddo
      local_press(ivertex) = press_loc_p(ghosted_id) - press_init_loc_p(ghosted_id)  ! p - p_0
      local_temp(ivertex) = temp_loc_p(ghosted_id) - temp_init_loc_p(ghosted_id)     ! T - T_0
      if (geomech_parameter%thermal_exp_coeff_spatially_varying) then
        alpha_vec(ivertex) = temp_thermal_exp_coeff_loc_p(ghosted_id)
      else
        alpha_vec(ivertex) = &
          geomech_parameter%thermal_exp_coeff(nint(imech_loc_p(ghosted_id)))
      endif
      if (geomech_parameter%biot_coeff_spatially_varying) then
        beta_vec(ivertex) = temp_biot_coeff_loc_p(ghosted_id)
      else
        beta_vec(ivertex) = &
          geomech_parameter%biot_coeff(nint(imech_loc_p(ghosted_id)))
      endif
      if (geomech_parameter%density_spatially_varying) then
        density_rock_vec(ivertex) = temp_density_loc_p(ghosted_id)
      else
        density_rock_vec(ivertex) = &
          geomech_parameter%density(nint(imech_loc_p(ghosted_id)))
      endif
      if (geomech_parameter%youngs_modulus_spatially_varying) then
        youngs_vec(ivertex) = temp_youngs_modulus_loc_p(ghosted_id)
      else
        youngs_vec(ivertex) = &
          geomech_parameter%youngs_modulus(nint(imech_loc_p(ghosted_id)))
      endif
      if (geomech_parameter%poissons_ratio_spatially_varying) then
        poissons_vec(ivertex) = temp_poissons_ratio_loc_p(ghosted_id)
      else
        poissons_vec(ivertex) = &
          geomech_parameter%poissons_ratio(nint(imech_loc_p(ghosted_id)))
      endif
      density_fluid_vec(ivertex) = fluid_density_loc_p(ghosted_id)
      porosity_vec(ivertex) = porosity_loc_p(ghosted_id)
      density_bulk_vec(ivertex) = (porosity_vec(ivertex) * &
                                   density_fluid_vec(ivertex)) + &
                                  ((1.d0 - porosity_vec(ivertex)) * &
                                   density_rock_vec(ivertex))
      local_body_force(ivertex,GEOMECH_DISP_X_DOF) = &
        temp_body_force_x_loc_p(ghosted_id)
      local_body_force(ivertex,GEOMECH_DISP_Y_DOF) = &
        temp_body_force_y_loc_p(ghosted_id)
      local_body_force(ivertex,GEOMECH_DISP_Z_DOF) = &
        temp_body_force_z_loc_p(ghosted_id)
    enddo
    call GeomechForceLocalElemRHS(size_elenodes,local_coordinates(1:size_elenodes,:), &
       local_press(1:size_elenodes),local_temp(1:size_elenodes), &
       youngs_vec(1:size_elenodes),poissons_vec(1:size_elenodes), &
       density_bulk_vec(1:size_elenodes),beta_vec(1:size_elenodes), &
       local_body_force(1:size_elenodes,:), &
       alpha_vec(1:size_elenodes),eletype, &
       grid%gauss_node(ielem)%dim,grid%gauss_node(ielem)%r, &
       grid%gauss_node(ielem)%w,rhs_local_vec(1:ndofs),option)
    call VecSetValues(rhs,ndofs,ids(1:ndofs),rhs_local_vec(1:ndofs),ADD_VALUES, &
                      ierr);CHKERRQ(ierr)
  enddo

  deallocate(elenodes)
  deallocate(local_coordinates)
  deallocate(petsc_ids)
  deallocate(ids)
  deallocate(rhs_local_vec)
  deallocate(local_press)
  deallocate(local_temp)
  deallocate(beta_vec)
  deallocate(alpha_vec)
  deallocate(density_rock_vec)
  deallocate(density_fluid_vec)
  deallocate(youngs_vec)
  deallocate(poissons_vec)
  deallocate(porosity_vec)
  deallocate(density_bulk_vec)
  deallocate(local_body_force)

  call VecRestoreArray(field%press_loc,press_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%temp_loc,temp_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%porosity_loc,porosity_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%fluid_density_loc,fluid_density_loc_p, &
                          ierr);CHKERRQ(ierr)

  if (geomech_parameter%youngs_modulus_spatially_varying) then
    call VecRestoreArray(field%youngs_modulus_loc,temp_youngs_modulus_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%poissons_ratio_spatially_varying) then
    call VecRestoreArray(field%poissons_ratio_loc,temp_poissons_ratio_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%density_spatially_varying) then
    call VecRestoreArray(field%density_loc,temp_density_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%biot_coeff_spatially_varying) then
    call VecRestoreArray(field%biot_coeff_loc,temp_biot_coeff_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%thermal_exp_coeff_spatially_varying) then
    call VecRestoreArray(field%thermal_exp_coeff_loc,temp_thermal_exp_coeff_loc_p,ierr);CHKERRQ(ierr)
  endif
  call VecRestoreArray(field%body_force_x_loc,temp_body_force_x_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%body_force_y_loc,temp_body_force_y_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%body_force_z_loc,temp_body_force_z_loc_p,ierr);CHKERRQ(ierr)

  call VecRestoreArray(field%press_init_loc,press_init_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%temp_init_loc,temp_init_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%fluid_density_init_loc,fluid_density_init_loc_p, &
                          ierr);CHKERRQ(ierr)

  call VecAssemblyBegin(rhs,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(rhs,ierr);CHKERRQ(ierr)

  ! Pre-allocate traction-face work arrays once.
  max_face_vertices = 0
  boundary_condition => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    if (associated(boundary_condition%geomech_condition%traction)) then
      if (boundary_condition%geomech_condition%traction%itype == NEUMANN_BC) then
        nfaces = boundary_condition%region%sideset%nfaces
        do iface = 1, nfaces
          num_vertices = size(boundary_condition%region%sideset%face_vertices(:,iface))
          max_face_vertices = max(max_face_vertices, num_vertices)
        enddo
      endif
    endif
    boundary_condition => boundary_condition%next
  enddo
  max_face_vertices = max(1,max_face_vertices)
  allocate(face_vertices(max_face_vertices))
  allocate(face_local_coordinates(max_face_vertices,THREE_INTEGER))
  allocate(face_petsc_ids(max_face_vertices))
  allocate(face_ids(max_face_vertices*option%ngeomechdof))
  allocate(face_rhs_local_vec(max_face_vertices*option%ngeomechdof))

  ! Traction
  boundary_condition => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    if (associated(boundary_condition%geomech_condition%traction)) then
      select case(boundary_condition%geomech_condition%traction%itype)
        case(DIRICHLET_BC)
          call GeomechForceError('Dirichlet BC for traction not available.', option)
        case(ZERO_GRADIENT_BC)
         ! do nothing
        case(NEUMANN_BC)
          stress_bc = boundary_condition%geomech_condition%traction% &
                        dataset%rarray
          nfaces = boundary_condition%region%sideset%nfaces
          do iface = 1, nfaces
            num_vertices = size(boundary_condition%region%sideset% &
              face_vertices(:,iface))
            if (num_vertices == THREE_INTEGER) then
              facetype = TRI_FACE_TYPE
            else
              facetype = QUAD_FACE_TYPE
            endif
            face_ndofs = num_vertices*option%ngeomechdof
            face_vertices(1:num_vertices) = boundary_condition%region%sideset% &
              face_vertices(1:num_vertices,iface)
            do ivertex = 1, num_vertices
              ghosted_id = face_vertices(ivertex)
              face_local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = &
                grid%nodes(ghosted_id)%x
              face_local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = &
                grid%nodes(ghosted_id)%y
              face_local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = &
                grid%nodes(ghosted_id)%z
              face_petsc_ids(ivertex) = grid%node_ids_ghosted_petsc(ghosted_id)
              do idof = 1, option%ngeomechdof
                face_ids(idof + (ivertex-1)*option%ngeomechdof) = &
                  (face_petsc_ids(ivertex)-1)*option%ngeomechdof + (idof-1)
              enddo
            enddo
            call GeomechForceApplyTractionBCtoRHS( &
                   face_local_coordinates(1:num_vertices,:), &
                   facetype, &
                   stress_bc, &
                   grid%gauss_surf_node(facetype)%r, &
                   grid%gauss_surf_node(facetype)%w, &
                   face_rhs_local_vec(1:face_ndofs),option)
            call VecSetValues(rhs,face_ndofs,face_ids(1:face_ndofs), &
                   face_rhs_local_vec(1:face_ndofs),ADD_VALUES, &
                   ierr);CHKERRQ(ierr)
          enddo
      end select
    endif
    boundary_condition => boundary_condition%next
  enddo
  deallocate(face_vertices)
  deallocate(face_local_coordinates)
  deallocate(face_petsc_ids)
  deallocate(face_ids)
  deallocate(face_rhs_local_vec)
  call VecAssemblyBegin(rhs,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(rhs,ierr);CHKERRQ(ierr)

  ! Force boundary conditions
  boundary_condition => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      petsc_id = grid%node_ids_ghosted_petsc(ghosted_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      ! X force
      if (associated(boundary_condition%geomech_condition%force_x)) then
        select case(boundary_condition%geomech_condition%force_x%itype)
          case(DIRICHLET_BC)
            call VecSetValue(rhs, &
                         (petsc_id-1)*option% &
                           ngeomechdof+GEOMECH_DISP_X_DOF-1, &
                         boundary_condition% &
                           geomech_aux_real_var(GEOMECH_DISP_X_DOF,ivertex), &
                         ADD_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            call GeomechForceError('Neumann BC for force not available.', option)
        end select
      endif

       ! Y force
      if (associated(boundary_condition%geomech_condition%force_y)) then
        select case(boundary_condition%geomech_condition%force_y%itype)
          case(DIRICHLET_BC)
            call VecSetValue(rhs, &
                         (petsc_id-1)*option% &
                           ngeomechdof+GEOMECH_DISP_Y_DOF-1, &
                         boundary_condition% &
                           geomech_aux_real_var(GEOMECH_DISP_Y_DOF,ivertex), &
                         ADD_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            call GeomechForceError('Neumann BC for force not available.', option)

        end select
      endif

       ! Z force
      if (associated(boundary_condition%geomech_condition%force_z)) then
        select case(boundary_condition%geomech_condition%force_z%itype)
          case(DIRICHLET_BC)
            call VecSetValue(rhs, &
                        (petsc_id-1)*option% &
                          ngeomechdof+GEOMECH_DISP_Z_DOF-1, &
                         boundary_condition% &
                           geomech_aux_real_var(GEOMECH_DISP_Z_DOF,ivertex), &
                         ADD_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            call GeomechForceError('Neumann BC for force not available.', option)
        end select
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecAssemblyBegin(rhs,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(rhs,ierr);CHKERRQ(ierr)

  ! Find the boundary nodes with dirichlet and set the residual at those nodes
  ! to zero, later set the Jacobian to 1

  ! displacement boundary conditions
  ! jaa: apply displacement to RHS last
  boundary_condition => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      petsc_id = grid%node_ids_ghosted_petsc(ghosted_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      ! X displacement
      if (associated(boundary_condition%geomech_condition%displacement_x)) then
        select case(boundary_condition%geomech_condition%displacement_x%itype)
          case(DIRICHLET_BC)
            call VecSetValue(rhs, &
                             (petsc_id-1)*option% &
                               ngeomechdof+GEOMECH_DISP_X_DOF-1, &
                             boundary_condition% &
                             geomech_aux_real_var(GEOMECH_DISP_X_DOF, &
                             ivertex),INSERT_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            call GeomechForceError('Neumann BC for displacement not available.', option)
        end select
      endif

      ! Y displacement
      if (associated(boundary_condition%geomech_condition%displacement_y)) then
        select case(boundary_condition%geomech_condition%displacement_y%itype)
          case(DIRICHLET_BC)
            call VecSetValue(rhs, &
                             (petsc_id-1)*option% &
                               ngeomechdof+GEOMECH_DISP_Y_DOF-1, &
                             boundary_condition% &
                             geomech_aux_real_var(GEOMECH_DISP_Y_DOF, &
                             ivertex),INSERT_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            call GeomechForceError('Neumann BC for displacement not available.', option)
        end select
      endif

      ! Z displacement
      if (associated(boundary_condition%geomech_condition%displacement_z)) then
        select case(boundary_condition%geomech_condition%displacement_z%itype)
          case(DIRICHLET_BC)
            call VecSetValue(rhs, &
                             (petsc_id-1)*option% &
                               ngeomechdof+GEOMECH_DISP_Z_DOF-1, &
                             boundary_condition% &
                             geomech_aux_real_var(GEOMECH_DISP_Z_DOF, &
                             ivertex),INSERT_VALUES,ierr);CHKERRQ(ierr)
          case(ZERO_GRADIENT_BC)
           ! do nothing
          case(NEUMANN_BC)
            call GeomechForceError('Neumann BC for displacement not available.', option)
        end select
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Need to assemby here since one cannot mix INSERT_VALUES
  ! and ADD_VALUES
  call VecAssemblyBegin(rhs,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(rhs,ierr);CHKERRQ(ierr)

end subroutine GeomechForceSetupLinearSystem

! ************************************************************************** !

subroutine GeomechForceLocalElemRHS(size_elenodes,local_coordinates, &
                                         local_press,local_temp, &
                                         local_youngs,local_poissons, &
                                         local_density,local_beta, &
                                         local_body_force, &
                                         local_alpha, &
                                         eletype,dim,r,w,rhs_vec,option)
  !
  ! Computes the RHS for a local element
  !
  ! Author: Satish Karra
  ! Date: 06/24/13
  !

  use Grid_Unstructured_Cell_module
  use Shape_Function_module
  use Option_module
  use Utility_module

  implicit none

  type(shapefunction_type) :: shapefunction
  type(option_type) :: option

  PetscInt,  intent(in) :: size_elenodes, eletype, dim
  PetscReal, intent(in) :: local_coordinates(:,:)        ! (nen,3)
  PetscReal, intent(in) :: local_press(:), local_temp(:)
  PetscReal, intent(in) :: local_youngs(:), local_poissons(:)
  PetscReal, intent(in) :: local_density(:), local_beta(:), local_alpha(:)
  PetscReal, intent(in) :: local_body_force(:,:)
  PetscReal, pointer, intent(in) :: r(:,:), w(:)
  PetscReal, intent(inout) :: rhs_vec(:)                 ! (3*nen)

  PetscInt  :: igpt, len_w, a, ia
  PetscReal :: J_map(3,3), inv_J_map(3,3), detJ_map
  PetscReal :: dNdx(size_elenodes,3)
  PetscReal :: lambda, mu, beta, alpha, density, youngs_mod, poissons_ratio
  PetscReal :: bf(3)
  PetscReal :: wdet, dp, dT, coefP, coefT
  PetscReal :: gauss_tet_vol_weight(4)
  PetscReal :: gauss_tot_weight
  integer :: i

  rhs_vec = 0.0d0
  len_w   = size(w)

  if (dim /= 3) then
    call GeomechForceError('GEOMECHANICS: RHS routine expects dim=3.', option)
    return
  end if

  ! Optional improved tet weighting (unchanged logic)
  if(option%geomechanics%improve_tet_weighting .and. eletype == TET_TYPE .and. len_w == 4) then
    call ComputeTetVolAtVertex(local_coordinates(1,:), &
                               local_coordinates(2,:), &
                               local_coordinates(3,:), &
                               local_coordinates(4,:), &
                               gauss_tet_vol_weight(1))
    call ComputeTetVolAtVertex(local_coordinates(2,:), &
                               local_coordinates(3,:), &
                               local_coordinates(4,:), &
                               local_coordinates(1,:), &
                               gauss_tet_vol_weight(2))
    call ComputeTetVolAtVertex(local_coordinates(3,:), &
                               local_coordinates(4,:), &
                               local_coordinates(1,:), &
                               local_coordinates(2,:), &
                               gauss_tet_vol_weight(3))
    call ComputeTetVolAtVertex(local_coordinates(4,:), &
                               local_coordinates(1,:), &
                               local_coordinates(2,:), &
                               local_coordinates(3,:), &
                               gauss_tet_vol_weight(4))
    gauss_tot_weight = sum(gauss_tet_vol_weight)
    do i=1,4
      gauss_tet_vol_weight(i) = gauss_tet_vol_weight(i)/gauss_tot_weight
    end do
  end if

  shapefunction%element_type = eletype
  call ShapeFunctionInitialize(shapefunction,option)

  do igpt = 1, len_w
    shapefunction%zeta = r(igpt,:)
    call ShapeFunctionCalculate(shapefunction,option)

    ! Map reference derivatives to physical space via J = x^T * dN/dzeta.
    ! Geometric mapping and derivative transform at this quadrature point.
    ! Vertex-point evaluation uses same Jacobian/gradient transform as quadrature.
    J_map = matmul(transpose(local_coordinates), shapefunction%DN)

    call MatInv3WithDet(J_map, inv_J_map, detJ_map)
    if (detJ_map <= 0.0d0) then
      call GeomechForceError('GEOMECHANICS: Determinant of J_map has to be positive!', option)
      call ShapeFunctionDestroy(shapefunction)
      return
    end if

    ! Physical shape gradients: dN/dx = dN/dzeta * J^{-1}.
    dNdx = matmul(shapefunction%DN, inv_J_map)

    ! Interpolate material and coupling properties at the quadrature point.
    youngs_mod     = dot_product(shapefunction%N, local_youngs)
    poissons_ratio = dot_product(shapefunction%N, local_poissons)
    alpha          = dot_product(shapefunction%N, local_alpha)
    beta           = dot_product(shapefunction%N, local_beta)
    density        = dot_product(shapefunction%N, local_density)

    call GeomechGetLambdaMu(lambda, mu, youngs_mod, poissons_ratio)
    bf(GEOMECH_DISP_X_DOF) = dot_product(shapefunction%N, &
                       local_body_force(1:size_elenodes,GEOMECH_DISP_X_DOF))
    bf(GEOMECH_DISP_Y_DOF) = dot_product(shapefunction%N, &
                       local_body_force(1:size_elenodes,GEOMECH_DISP_Y_DOF))
    bf(GEOMECH_DISP_Z_DOF) = dot_product(shapefunction%N, &
                       local_body_force(1:size_elenodes,GEOMECH_DISP_Z_DOF))

    ! Effective volume weight for this quadrature point.
    ! Integrate C:(grad N_a, grad N_b) with isotropic Lamé parameters.
    wdet = w(igpt) * detJ_map

    ! Body force term: rhs += wdet*density*(N ⊗ I)*bf
    if(option%geomechanics%improve_tet_weighting .and. eletype == TET_TYPE .and. len_w == 4) then
      wdet = wdet * 4.0d0 * gauss_tet_vol_weight(igpt)
    end if

    ! Body-force contribution: rhs += int( rho * N_a * b dV ).
    do a = 1, size_elenodes
      ia = 3*(a-1)
      rhs_vec(ia+1) = rhs_vec(ia+1) + wdet * density * shapefunction%N(a) * bf(1)
      rhs_vec(ia+2) = rhs_vec(ia+2) + wdet * density * shapefunction%N(a) * bf(2)
      rhs_vec(ia+3) = rhs_vec(ia+3) + wdet * density * shapefunction%N(a) * bf(3)
    end do

    wdet = w(igpt) * detJ_map ! reset the wdet for the coupling terms since it may have been modified by the tet weighting
    ! Pressure and thermal coupling terms: old code uses + beta*dp * vecB + alpha*(3λ+2μ)*dT * vecB
    dp    = dot_product(shapefunction%N, local_press)
    dT    = dot_product(shapefunction%N, local_temp)
    coefP = wdet * beta * dp
    coefT = wdet * alpha * (3.0d0*lambda + 2.0d0*mu) * dT

    ! Combine pressure and thermal isotropic source terms before projection.
    coefP = coefP + coefT
    do a = 1, size_elenodes
      ia = 3*(a-1)
      rhs_vec(ia+1) = rhs_vec(ia+1) + coefP * dNdx(a,1)
      rhs_vec(ia+2) = rhs_vec(ia+2) + coefP * dNdx(a,2)
      rhs_vec(ia+3) = rhs_vec(ia+3) + coefP * dNdx(a,3)
    end do
  end do

  call ShapeFunctionDestroy(shapefunction)

end subroutine GeomechForceLocalElemRHS


! ************************************************************************** !

subroutine GeomechForceAssembleCoeffMatrixLocal(size_elenodes, &
                                                local_coordinates, &
                                                local_youngs, &
                                                local_poissons, &
                                                petsc_ids, &
                                                eletype,dim,r,w,Amat,option)
  !
  ! Forms the local element stiffness matrix and adds it to the global matrix.
  !
  ! Author: Satish Karra
  ! Date: 06/24/13
  ! Updated: 02/24/26

  use Shape_Function_module
  use Option_module
  use Utility_module
  use Petsc_Utility_module

  implicit none

  type(shapefunction_type) :: shapefunction
  type(option_type) :: option
  PetscInt, intent(in) :: size_elenodes, eletype, dim
  PetscReal, intent(in) :: local_coordinates(:,:)     ! (nen,3)
  PetscReal, intent(in) :: local_youngs(:), local_poissons(:)
  PetscInt, intent(in) :: petsc_ids(:)
  PetscReal, pointer, intent(in) :: r(:,:), w(:)
  Mat, intent(inout) :: Amat

  PetscInt :: igpt, len_w, ia, ib
  PetscReal :: J_map(3,3), inv_J_map(3,3), detJ_map
  PetscReal :: dNdx(size_elenodes,3)
  PetscReal :: youngs_mod, poissons_ratio, lambda, mu, wdet
  PetscReal :: gax,gay,gaz, gbx,gby,gbz, gg
  PetscReal :: Kblock(3,3,size_elenodes,size_elenodes)
  PetscReal :: Jac_sub_mat(3,3)
  PetscErrorCode :: ierr

  if (dim /= 3) then
    call GeomechForceError('GEOMECHANICS: this stiffness routine expects dim=3.', option)
    return
  end if

  Kblock = 0.0d0
  len_w = size(w)

  shapefunction%element_type = eletype
  call ShapeFunctionInitialize(shapefunction,option)

  do igpt = 1, len_w
    shapefunction%zeta = r(igpt,:)
    call ShapeFunctionCalculate(shapefunction,option)

    J_map = matmul(transpose(local_coordinates), shapefunction%DN)
    call MatInv3WithDet(J_map, inv_J_map, detJ_map)
    if (detJ_map <= 0.0d0) then
      call GeomechForceError('GEOMECHANICS: det(J) must be positive!', option)
      call ShapeFunctionDestroy(shapefunction)
      return
    end if

    dNdx = matmul(shapefunction%DN, inv_J_map)

    youngs_mod     = dot_product(shapefunction%N, local_youngs)
    poissons_ratio = dot_product(shapefunction%N, local_poissons)
    call GeomechGetLambdaMu(lambda, mu, youngs_mod, poissons_ratio)

    wdet = w(igpt) * detJ_map

    do ia = 1, size_elenodes
      gax = dNdx(ia,1); gay = dNdx(ia,2); gaz = dNdx(ia,3)
      do ib = 1, size_elenodes
        gbx = dNdx(ib,1); gby = dNdx(ib,2); gbz = dNdx(ib,3)

        gg = gax*gbx + gay*gby + gaz*gbz

        Kblock(1,1,ia,ib) = Kblock(1,1,ia,ib) + wdet*( lambda*gax*gbx + mu*gg + mu*gax*gbx )
        Kblock(1,2,ia,ib) = Kblock(1,2,ia,ib) + wdet*( lambda*gax*gby          + mu*gay*gbx )
        Kblock(1,3,ia,ib) = Kblock(1,3,ia,ib) + wdet*( lambda*gax*gbz          + mu*gaz*gbx )

        Kblock(2,1,ia,ib) = Kblock(2,1,ia,ib) + wdet*( lambda*gay*gbx          + mu*gax*gby )
        Kblock(2,2,ia,ib) = Kblock(2,2,ia,ib) + wdet*( lambda*gay*gby + mu*gg + mu*gay*gby )
        Kblock(2,3,ia,ib) = Kblock(2,3,ia,ib) + wdet*( lambda*gay*gbz          + mu*gaz*gby )

        Kblock(3,1,ia,ib) = Kblock(3,1,ia,ib) + wdet*( lambda*gaz*gbx          + mu*gax*gbz )
        Kblock(3,2,ia,ib) = Kblock(3,2,ia,ib) + wdet*( lambda*gaz*gby          + mu*gay*gbz )
        Kblock(3,3,ia,ib) = Kblock(3,3,ia,ib) + wdet*( lambda*gaz*gbz + mu*gg + mu*gaz*gbz )
      end do
    end do
  end do

  call ShapeFunctionDestroy(shapefunction)

  do ia = 1, size_elenodes
    do ib = 1, size_elenodes
      Jac_sub_mat = Kblock(:,:,ia,ib)
      call PUMSetValuesBlocked(Amat,1,petsc_ids(ia)-1,1,petsc_ids(ib)-1, &
                               Jac_sub_mat,ADD_VALUES,ierr);CHKERRQ(ierr)
    end do
  end do
end subroutine GeomechForceAssembleCoeffMatrixLocal

! ************************************************************************** !

subroutine GeomechGetLambdaMu(lambda,mu,E,nu)
  !
  ! Gets the material properties given the position
  ! of the point
  !
  ! Author: Satish Karra
  ! Date: 06/24/13
  !

  PetscReal :: lambda, mu
  PetscReal :: E, nu

  lambda = E*nu/(1.d0+nu)/(1.d0-2.d0*nu)
  mu = E/2.d0/(1.d0+nu)


end subroutine GeomechGetLambdaMu

! ************************************************************************** !

subroutine GeomechGetBodyForce(bf,option)
  !
  ! Gets the body force vector
  !
  ! Author: Satish Karra
  ! Date: 06/24/13
  !

  use Option_module

  implicit none

  type(option_type) :: option

  PetscReal :: bf(THREE_INTEGER)

  bf = 0.d0

  bf(GEOMECH_DISP_X_DOF) = option%geomechanics%gravity(X_DIRECTION)
  bf(GEOMECH_DISP_Y_DOF) = option%geomechanics%gravity(Y_DIRECTION)
  bf(GEOMECH_DISP_Z_DOF) = option%geomechanics%gravity(Z_DIRECTION)

end subroutine GeomechGetBodyForce

! ************************************************************************** !

subroutine GeomechForceAssembleCoeffMatrix(A,geomech_realization)
  !
  ! Computes the Assembled Coefficient Matrix
  !
  ! Author: Satish Karra
  ! Date: 06/21/13
  ! Modified: 07/12/16, 06/12/25
  !

  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Coupler_module
  use Geomechanics_Field_module
  use Geomechanics_Debug_module
  use Geomechanics_Discretization_module
  use Option_module
  use Grid_Unstructured_Cell_module
  use Geomechanics_Region_module
  use Geomechanics_Auxiliary_module
  use Petsc_Utility_module

  implicit none

  Mat :: A

  PetscErrorCode :: ierr

  class(realization_geomech_type) :: geomech_realization
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_patch_type), pointer :: patch
  type(geomech_field_type), pointer :: field
  type(geomech_grid_type), pointer :: grid
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)
  type(option_type), pointer :: option
  type(gm_region_type), pointer :: region
  type(geomech_coupler_type), pointer :: boundary_condition
  type(geomech_linear_parameter_type), pointer :: geomech_parameter

  PetscInt, allocatable :: elenodes(:)
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscInt, allocatable :: ghosted_ids(:), petsc_ids(:)
  PetscInt, allocatable :: rows(:)
  PetscReal, allocatable :: youngs_vec(:), poissons_vec(:)
  PetscInt :: ielem,ivertex
  PetscInt :: ghosted_id
  PetscInt :: eletype
  PetscInt :: local_id, petsc_id
  PetscInt :: vertex_count, count
  PetscReal, pointer :: imech_loc_p(:)
  PetscInt :: size_elenodes, max_elem_nodes

  PetscReal, pointer :: temp_youngs_modulus_loc_p(:)
  PetscReal, pointer :: temp_poissons_ratio_loc_p(:)

  field => geomech_realization%geomech_field
  geomech_discretization => geomech_realization%geomech_discretization
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  geomech_global_aux_vars => patch%geomech_aux%Global%aux_vars
  geomech_parameter => patch%geomech_aux%Linear%linear_parameter

  call MatZeroEntries(A,ierr);CHKERRQ(ierr)
  call VecGetArray(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)

  if (geomech_parameter%youngs_modulus_spatially_varying) then
    call VecGetArray(field%youngs_modulus_loc,temp_youngs_modulus_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%poissons_ratio_spatially_varying) then
    call VecGetArray(field%poissons_ratio_loc,temp_poissons_ratio_loc_p,ierr);CHKERRQ(ierr)
  endif

  max_elem_nodes = maxval(grid%elem_nodes(0,1:grid%nlmax_elem))
  allocate(elenodes(max_elem_nodes))
  allocate(local_coordinates(max_elem_nodes,THREE_INTEGER))
  allocate(ghosted_ids(max_elem_nodes))
  allocate(petsc_ids(max_elem_nodes))
  allocate(youngs_vec(max_elem_nodes))
  allocate(poissons_vec(max_elem_nodes))

  ! Loop over elements on a processor
  do ielem = 1, grid%nlmax_elem
    size_elenodes = grid%elem_nodes(0,ielem)
    elenodes(1:size_elenodes) = grid%elem_nodes(1:size_elenodes,ielem)
    eletype = grid%gauss_node(ielem)%entity_type
    do ivertex = 1, size_elenodes
      ghosted_id = elenodes(ivertex)
      local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = grid%nodes(ghosted_id)%x
      local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = grid%nodes(ghosted_id)%y
      local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = grid%nodes(ghosted_id)%z
      ghosted_ids(ivertex) = ghosted_id
      petsc_ids(ivertex) = grid%node_ids_ghosted_petsc(ghosted_id)
    enddo
    do ivertex = 1, size_elenodes
      ghosted_id = elenodes(ivertex)
      if (geomech_parameter%youngs_modulus_spatially_varying) then
        youngs_vec(ivertex) = temp_youngs_modulus_loc_p(ghosted_id)
      else
        youngs_vec(ivertex) = &
          geomech_parameter%youngs_modulus(nint(imech_loc_p(ghosted_id)))
      endif
      if (geomech_parameter%poissons_ratio_spatially_varying) then
        poissons_vec(ivertex) = temp_poissons_ratio_loc_p(ghosted_id)
      else
        poissons_vec(ivertex) = &
          geomech_parameter%poissons_ratio(nint(imech_loc_p(ghosted_id)))
      endif
    enddo
    call GeomechForceAssembleCoeffMatrixLocal(size_elenodes, &
       local_coordinates(1:size_elenodes,:), &
       youngs_vec(1:size_elenodes),poissons_vec(1:size_elenodes), &
       petsc_ids(1:size_elenodes),eletype,grid%gauss_node(ielem)%dim, &
       grid%gauss_node(ielem)%r,grid%gauss_node(ielem)%w,A,option)
  enddo

  deallocate(elenodes)
  deallocate(local_coordinates)
  deallocate(ghosted_ids)
  deallocate(petsc_ids)
  deallocate(youngs_vec)
  deallocate(poissons_vec)

  call VecRestoreArray(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)

  if (geomech_parameter%youngs_modulus_spatially_varying) then
    call VecRestoreArray(field%youngs_modulus_loc,temp_youngs_modulus_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%poissons_ratio_spatially_varying) then
    call VecRestoreArray(field%poissons_ratio_loc,temp_poissons_ratio_loc_p,ierr);CHKERRQ(ierr)
  endif

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! Find the boundary nodes with dirichlet and set the residual at those nodes
  ! to zero, later set the Jacobian to 1

  ! Find the number of boundary vertices
  vertex_count = 0
  boundary_condition => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    vertex_count = vertex_count + region%num_verts
    boundary_condition => boundary_condition%next
  enddo

  allocate(rows(vertex_count*option%ngeomechdof))
  count = 0

  boundary_condition => patch%geomech_boundary_condition_list%first
  do
    if (.not.associated(boundary_condition)) exit
    region => boundary_condition%region
    do ivertex = 1, region%num_verts
      local_id = region%vertex_ids(ivertex)
      ghosted_id = grid%nL2G(local_id)
      petsc_id = grid%node_ids_ghosted_petsc(ghosted_id)
      if (associated(patch%imat)) then
        if (patch%imat(ghosted_id) <= 0) cycle
      endif

      ! X displacement
      if (associated(boundary_condition%geomech_condition%displacement_x)) then
        select case(boundary_condition%geomech_condition%displacement_x%itype)
          case(DIRICHLET_BC)
            count = count + 1
            rows(count) = (ghosted_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_X_DOF-1
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif

      ! Y displacement
      if (associated(boundary_condition%geomech_condition%displacement_y)) then
        select case(boundary_condition%geomech_condition%displacement_y%itype)
          case(DIRICHLET_BC)
            count = count + 1
            rows(count) = (ghosted_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_Y_DOF-1
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif

      ! Z displacement
      if (associated(boundary_condition%geomech_condition%displacement_z)) then
        select case(boundary_condition%geomech_condition%displacement_z%itype)
          case(DIRICHLET_BC)
            count = count + 1
            rows(count) = (ghosted_id-1)*option%ngeomechdof + &
              GEOMECH_DISP_Z_DOF-1
          case(ZERO_GRADIENT_BC,NEUMANN_BC)
           ! do nothing
        end select
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo

  call MatZeroRowsLocal(A,count,rows,1.d0,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                        ierr);CHKERRQ(ierr)
  call MatSetOption(A,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE, &
                    ierr);CHKERRQ(ierr)
  call MatStoreValues(A,ierr);CHKERRQ(ierr)

  deallocate(rows)

end subroutine GeomechForceAssembleCoeffMatrix

! ************************************************************************** !

subroutine GeomechUpdateFromSubsurf(realization,geomech_realization)
  !
  ! The pressure/temperature from subsurface are
  ! mapped to geomech
  !
  ! Author: Satish Karra, LANL
  ! Date: 09/10/13
  !

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Geomechanics_Realization_class
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use GMDM_Pointer_module
  use Option_module
  use Option_Geomechanics_module
  use Grid_Structured_module

  implicit none

  class(realization_subsurface_type) :: realization
  class(realization_geomech_type) :: geomech_realization
  type(grid_type), pointer :: grid
  type(geomech_grid_type), pointer :: geomech_grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(geomech_field_type), pointer :: geomech_field
  type(gmdm_ptr_type), pointer :: dm_ptr

  PetscErrorCode :: ierr
  PetscReal, pointer :: vec_p(:), xx_loc_p(:), press_p(:), temp_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: nx, ny, nz, nxv, nyv, nxyv
  PetscInt :: vid, ix, iy, iz
  PetscInt :: i_start, i_end, j_start, j_end, k_start, k_end
  PetscInt :: n_adj
  PetscInt :: press_lb, temp_lb
  PetscInt :: idx
  type(grid_structured_type), pointer :: sgrid
  InsertMode :: scatter_insert_mode

  option        => realization%option
  grid          => realization%discretization%grid
  field         => realization%field
  geomech_grid  => geomech_realization%geomech_discretization%grid
  geomech_field => geomech_realization%geomech_field

  ! use the subsurface output option parameters for geomechanics as well
  geomech_realization%output_option%tunit = realization%output_option%tunit
  geomech_realization%output_option%tconv = realization%output_option%tconv

  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_realization% &
                                                   geomech_discretization, &
                                                   ONEDOF)

  scatter_insert_mode = INSERT_VALUES
  if (option%geomechanics%flow_interp_order == GEOMECH_FLOW_INTERP_ORDER_1ST) then
    scatter_insert_mode = ADD_VALUES
  endif


  ! pressure
  call VecGetArray(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(geomech_field%subsurf_vec_1dof,vec_p,ierr)
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    vec_p(local_id) = xx_loc_p(option%nflowdof*(ghosted_id-1)+1)
  enddo
  call VecRestoreArray(geomech_field%subsurf_vec_1dof,vec_p,ierr)
  call VecRestoreArray(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  ! Scatter the data
  if (scatter_insert_mode == ADD_VALUES) then
    call VecSet(geomech_field%press,0.d0,ierr);CHKERRQ(ierr)
  endif
  call VecScatterBegin(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                       geomech_field%subsurf_vec_1dof,geomech_field%press, &
                       scatter_insert_mode,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                     geomech_field%subsurf_vec_1dof,geomech_field%press, &
                     scatter_insert_mode,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  ! temperature
  if (option%nflowdof > 1) then
    call VecGetArray(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
    call VecGetArray(geomech_field%subsurf_vec_1dof,vec_p,ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      vec_p(local_id) = xx_loc_p(option%nflowdof*(ghosted_id-1)+2)
    enddo
    call VecRestoreArray(geomech_field%subsurf_vec_1dof,vec_p,ierr)
    call VecRestoreArray(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

    ! Scatter the data
    if (scatter_insert_mode == ADD_VALUES) then
      call VecSet(geomech_field%temp,0.d0,ierr);CHKERRQ(ierr)
    endif
    call VecScatterBegin(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                         geomech_field%subsurf_vec_1dof,geomech_field%temp, &
                         scatter_insert_mode,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterEnd(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                       geomech_field%subsurf_vec_1dof,geomech_field%temp, &
                       scatter_insert_mode,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  endif

  if (option%geomechanics%flow_interp_order == GEOMECH_FLOW_INTERP_ORDER_1ST .and. &
      trim(geomech_realization%geomech_discretization%ctype) == 'STRUCTURED_INTERNAL') then

    if (grid%itype /= STRUCTURED_GRID) then
      option%io_buffer = 'GEOMECHANICS_FLOW_INTERPOLATION ORDER 1 requires a structured flow grid.'
      call PrintErrMsg(option)
    endif

    sgrid => grid%structured_grid
    nx = sgrid%nx
    ny = sgrid%ny
    nz = sgrid%nz
    nxv = nx + 1
    nyv = ny + 1
    nxyv = nxv*nyv

    call VecGetArray(geomech_field%press,press_p,ierr);CHKERRQ(ierr)
    press_lb = lbound(press_p,ONE_INTEGER)
    if (option%nflowdof > 1) call VecGetArray(geomech_field%temp,temp_p,ierr);CHKERRQ(ierr)
    if (option%nflowdof > 1) temp_lb = lbound(temp_p,ONE_INTEGER)
    do local_id = 1, geomech_grid%nlmax_node
      vid = geomech_grid%node_ids_local_natural(local_id)
      vid = vid - 1
      iz = vid/nxyv
      iy = mod(vid,nxyv)/nxv
      ix = mod(vid,nxv)

      i_start = max(1,ix)
      i_end   = min(nx,ix+1)
      j_start = max(1,iy)
      j_end   = min(ny,iy+1)
      k_start = max(1,iz)
      k_end   = min(nz,iz+1)
      n_adj = (i_end - i_start + 1)*(j_end - j_start + 1)*(k_end - k_start + 1)

      idx = press_lb + local_id - 1
      press_p(idx) = press_p(idx)/dble(n_adj)
      if (option%nflowdof > 1) then
        idx = temp_lb + local_id - 1
        temp_p(idx) = temp_p(idx)/dble(n_adj)
      endif
    enddo
    call VecRestoreArray(geomech_field%press,press_p,ierr);CHKERRQ(ierr)
    if (option%nflowdof > 1) call VecRestoreArray(geomech_field%temp,temp_p,ierr);CHKERRQ(ierr)
  endif

  call GeomechDiscretizationGlobalToLocal(&
                                geomech_realization%geomech_discretization, &
                                geomech_field%press, &
                                geomech_field%press_loc,ONEDOF)

  if (option%nflowdof > 1) &
    call GeomechDiscretizationGlobalToLocal(&
                                geomech_realization%geomech_discretization, &
                                geomech_field%temp, &
                                geomech_field%temp_loc,ONEDOF)

end subroutine GeomechUpdateFromSubsurf

! ************************************************************************** !

subroutine GeomechUpdateSubsurfFromGeomech(realization,geomech_realization)
  !
  ! The stresses and strains from geomech
  ! are mapped to subsurf.
  !
  ! Author: Satish Karra, LANL
  ! Date: 10/10/13
  !

  use Realization_Subsurface_class
  use Discretization_module
  use Grid_module
  use Field_module
  use Geomechanics_Realization_class
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use GMDM_Pointer_module
  use Option_module
  use Option_Geomechanics_module

  implicit none

  class(realization_subsurface_type) :: realization
  class(realization_geomech_type) :: geomech_realization
  type(grid_type), pointer :: grid
  type(geomech_grid_type), pointer :: geomech_grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(geomech_field_type), pointer :: geomech_field
  type(gmdm_ptr_type), pointer :: dm_ptr

  PetscErrorCode :: ierr
  PetscReal, pointer :: strain_subsurf_p(:), stress_subsurf_p(:)
  PetscInt :: local_id, idof
  PetscInt :: strain_lb, stress_lb
  PetscInt :: base, idx
  InsertMode :: scatter_insert_mode

  option        => realization%option
  grid          => realization%discretization%grid
  field         => realization%field
  geomech_grid  => geomech_realization%geomech_discretization%grid
  geomech_field => geomech_realization%geomech_field

  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_realization% &
                                                   geomech_discretization, &
                                                   ONEDOF)

  scatter_insert_mode = INSERT_VALUES
  if (option%geomechanics%flow_interp_order == GEOMECH_FLOW_INTERP_ORDER_1ST) then
    scatter_insert_mode = ADD_VALUES
    call VecSet(geomech_field%strain_subsurf,0.d0,ierr);CHKERRQ(ierr)
    call VecSet(geomech_field%stress_subsurf,0.d0,ierr);CHKERRQ(ierr)
  endif

  ! Scatter the strains
  call VecScatterBegin(dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                       geomech_field%strain,geomech_field%strain_subsurf, &
                       scatter_insert_mode,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                     geomech_field%strain,geomech_field%strain_subsurf, &
                     scatter_insert_mode,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  ! Scatter the stresses
  call VecScatterBegin(dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                       geomech_field%stress,geomech_field%stress_subsurf, &
                       scatter_insert_mode,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                     geomech_field%stress,geomech_field%stress_subsurf, &
                     scatter_insert_mode,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

  if (option%geomechanics%flow_interp_order == GEOMECH_FLOW_INTERP_ORDER_1ST .and. &
      trim(geomech_realization%geomech_discretization%ctype) == 'STRUCTURED_INTERNAL') then
    call VecGetArray(geomech_field%strain_subsurf,strain_subsurf_p,ierr);CHKERRQ(ierr)
    call VecGetArray(geomech_field%stress_subsurf,stress_subsurf_p,ierr);CHKERRQ(ierr)
    strain_lb = lbound(strain_subsurf_p,ONE_INTEGER)
    stress_lb = lbound(stress_subsurf_p,ONE_INTEGER)
    do local_id = 1, grid%nlmax
      base = (local_id-1)*SIX_INTEGER
      do idof = 1, SIX_INTEGER
        idx = strain_lb + base + idof - 1
        strain_subsurf_p(idx) = strain_subsurf_p(idx)/8.d0
        idx = stress_lb + base + idof - 1
        stress_subsurf_p(idx) = stress_subsurf_p(idx)/8.d0
      enddo
    enddo
    call VecRestoreArray(geomech_field%strain_subsurf,strain_subsurf_p,ierr);CHKERRQ(ierr)
    call VecRestoreArray(geomech_field%stress_subsurf,stress_subsurf_p,ierr);CHKERRQ(ierr)
  endif

  ! Scatter from global to local vectors
  call DiscretizationGlobalToLocal(realization%discretization, &
                                   geomech_field%strain_subsurf, &
                                   geomech_field%strain_subsurf_loc, &
                                   NGEODOF)
  call DiscretizationGlobalToLocal(realization%discretization, &
                                   geomech_field%stress_subsurf, &
                                   geomech_field%stress_subsurf_loc, &
                                   NGEODOF)

end subroutine GeomechUpdateSubsurfFromGeomech

! ************************************************************************** !

subroutine GeomechCreateGeomechSubsurfVec(realization,geomech_realization)
  !
  ! Creates the MPI vector that stores the
  ! variables from subsurface
  !
  ! Author: Satish Karra, LANL
  ! Date: 09/10/13
  !

  use Grid_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Geomechanics_Field_module
  use String_module
  use Realization_Subsurface_class
  use Option_module

  implicit none

  class(realization_subsurface_type) :: realization
  class(realization_geomech_type) :: geomech_realization

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(geomech_field_type), pointer :: geomech_field

  PetscErrorCode :: ierr

  option     => realization%option
  grid       => realization%discretization%grid
  geomech_field => geomech_realization%geomech_field

  call VecCreate(option%mycomm,geomech_field%subsurf_vec_1dof, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%subsurf_vec_1dof,grid%nlmax,PETSC_DECIDE, &
                   ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%subsurf_vec_1dof,ierr);CHKERRQ(ierr)

end subroutine GeomechCreateGeomechSubsurfVec

! ************************************************************************** !

subroutine GeomechCreateSubsurfStressStrainVec(realization,geomech_realization)
  !
  ! Creates the subsurface stress and strain
  ! MPI vectors to store information from geomechanics
  !
  ! Author: Satish Karra, LANL
  ! Date: 10/10/13
  !

  use Grid_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Geomechanics_Field_module
  use String_module
  use Realization_Subsurface_class
  use Option_module

  implicit none

  class(realization_subsurface_type) :: realization
  class(realization_geomech_type) :: geomech_realization

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(geomech_field_type), pointer :: geomech_field

  PetscErrorCode :: ierr

  option     => realization%option
  grid       => realization%discretization%grid
  geomech_field => geomech_realization%geomech_field

  ! strain
  call VecCreate(option%mycomm,geomech_field%strain_subsurf, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%strain_subsurf,grid%nlmax*SIX_INTEGER, &
                   PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(geomech_field%strain_subsurf,SIX_INTEGER, &
                       ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%strain_subsurf,ierr);CHKERRQ(ierr)

  ! stress
  call VecCreate(option%mycomm,geomech_field%stress_subsurf, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%stress_subsurf,grid%nlmax*SIX_INTEGER, &
                   PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(geomech_field%stress_subsurf,SIX_INTEGER, &
                       ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%stress_subsurf,ierr);CHKERRQ(ierr)

  ! strain_loc
  call VecCreate(PETSC_COMM_SELF,geomech_field%strain_subsurf_loc, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%strain_subsurf_loc,grid%ngmax*SIX_INTEGER, &
                   PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(geomech_field%strain_subsurf_loc,SIX_INTEGER, &
                       ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%strain_subsurf_loc,ierr);CHKERRQ(ierr)

  ! stress_loc
  call VecCreate(PETSC_COMM_SELF,geomech_field%stress_subsurf_loc, &
                 ierr);CHKERRQ(ierr)
  call VecSetSizes(geomech_field%stress_subsurf_loc,grid%ngmax*SIX_INTEGER, &
                   PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call VecSetBlockSize(geomech_field%stress_subsurf_loc,SIX_INTEGER, &
                       ierr);CHKERRQ(ierr)
  call VecSetFromOptions(geomech_field%stress_subsurf_loc,ierr);CHKERRQ(ierr)

end subroutine GeomechCreateSubsurfStressStrainVec

! ************************************************************************** !

subroutine GeomechForceStressStrain(geomech_realization)
  !
  ! Computes the stress strain on a patch
  !
  ! Author: Satish Karra
  ! Date: 09/17/13
  !

  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Geomechanics_Discretization_module
  use Geomechanics_Patch_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Grid_module
  use Grid_Unstructured_Cell_module
  use Geomechanics_Region_module
  use Geomechanics_Coupler_module
  use Option_module
  use Geomechanics_Auxiliary_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_patch_type), pointer :: patch
  type(geomech_field_type), pointer :: field
  type(geomech_grid_type), pointer :: grid
  type(geomech_global_auxvar_type), pointer :: geomech_global_aux_vars(:)
  type(option_type), pointer :: option
  type(geomech_linear_parameter_type), pointer :: geomech_parameter

  PetscInt, allocatable :: elenodes(:)
  PetscReal, allocatable :: local_coordinates(:,:)
  PetscReal, allocatable :: local_disp(:,:)
  PetscInt, allocatable :: petsc_ids(:)
  PetscInt, allocatable :: ids(:)
  PetscReal, allocatable :: youngs_vec(:), poissons_vec(:)
  PetscReal, allocatable :: alpha_vec(:)  ! thermal expansion coefficient
  PetscReal, allocatable :: beta_vec(:)  ! Biot's coefficient
  PetscReal, allocatable :: local_temp(:)
  PetscReal, allocatable :: local_press(:)
  PetscReal, allocatable :: strain(:,:), stress(:,:), stress_total(:,:)
  PetscInt :: ielem, ivertex
  PetscInt :: ghosted_id
  PetscInt :: eletype, idof
  PetscInt :: local_id
  PetscInt :: size_elenodes, max_elem_nodes
  PetscReal, pointer :: imech_loc_p(:)
  PetscReal, pointer :: strain_loc_p(:)
  PetscReal, pointer :: stress_loc_p(:)
  PetscReal, pointer :: stress_total_loc_p(:)
  PetscReal, pointer :: temp_loc_p(:), temp_init_loc_p(:)
  PetscReal, pointer :: press_loc_p(:), press_init_loc_p(:)
  PetscReal, pointer :: strain_p(:), stress_p(:), stress_total_p(:)
  PetscReal, pointer :: no_elems_p(:)

  PetscReal, pointer :: temp_youngs_modulus_loc_p(:)
  PetscReal, pointer :: temp_poissons_ratio_loc_p(:)
  PetscReal, pointer :: temp_biot_coeff_loc_p(:)
  PetscReal, pointer :: temp_thermal_exp_coeff_loc_p(:)

  PetscErrorCode :: ierr

  field => geomech_realization%geomech_field
  geomech_discretization => geomech_realization%geomech_discretization
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option
  geomech_global_aux_vars => patch%geomech_aux%Global%aux_vars
  geomech_parameter => patch%geomech_aux%Linear%linear_parameter

  call VecSet(field%strain,0.d0,ierr);CHKERRQ(ierr)
  call VecSet(field%stress,0.d0,ierr);CHKERRQ(ierr)
  call VecSet(field%stress_total,0.d0,ierr);CHKERRQ(ierr)

  call VecGetArray(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%strain_loc,strain_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%stress_loc,stress_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%stress_total_loc,stress_total_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%temp_loc,temp_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%temp_init_loc,temp_init_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%press_loc,press_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%press_init_loc,press_init_loc_p,ierr);CHKERRQ(ierr)

  if (geomech_parameter%youngs_modulus_spatially_varying) then
    call VecGetArray(field%youngs_modulus_loc,temp_youngs_modulus_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%poissons_ratio_spatially_varying) then
    call VecGetArray(field%poissons_ratio_loc,temp_poissons_ratio_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%biot_coeff_spatially_varying) then
    call VecGetArray(field%biot_coeff_loc,temp_biot_coeff_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%thermal_exp_coeff_spatially_varying) then
    call VecGetArray(field%thermal_exp_coeff_loc,temp_thermal_exp_coeff_loc_p,ierr);CHKERRQ(ierr)
  endif

  strain_loc_p = 0.d0
  stress_loc_p = 0.d0
  stress_total_loc_p = 0.d0

  max_elem_nodes = maxval(grid%elem_nodes(0,1:grid%nlmax_elem))
  allocate(elenodes(max_elem_nodes))
  allocate(local_coordinates(max_elem_nodes,THREE_INTEGER))
  allocate(local_disp(max_elem_nodes,option%ngeomechdof))
  allocate(petsc_ids(max_elem_nodes))
  allocate(local_temp(max_elem_nodes))
  allocate(local_press(max_elem_nodes))
  allocate(ids(max_elem_nodes*option%ngeomechdof))
  allocate(youngs_vec(max_elem_nodes))
  allocate(poissons_vec(max_elem_nodes))
  allocate(alpha_vec(max_elem_nodes))
  allocate(beta_vec(max_elem_nodes))
  allocate(strain(max_elem_nodes,SIX_INTEGER))
  allocate(stress(max_elem_nodes,SIX_INTEGER))
  allocate(stress_total(max_elem_nodes,SIX_INTEGER))

   ! Loop over elements on a processor
  do ielem = 1, grid%nlmax_elem
    size_elenodes = grid%elem_nodes(0,ielem)
    elenodes(1:size_elenodes) = grid%elem_nodes(1:size_elenodes,ielem)
    eletype = grid%gauss_node(ielem)%entity_type
    do ivertex = 1, size_elenodes
      ghosted_id = elenodes(ivertex)
      local_coordinates(ivertex,GEOMECH_DISP_X_DOF) = grid%nodes(ghosted_id)%x
      local_coordinates(ivertex,GEOMECH_DISP_Y_DOF) = grid%nodes(ghosted_id)%y
      local_coordinates(ivertex,GEOMECH_DISP_Z_DOF) = grid%nodes(ghosted_id)%z
      petsc_ids(ivertex) = grid%node_ids_ghosted_petsc(ghosted_id)
    enddo
    do ivertex = 1, size_elenodes
      ghosted_id = elenodes(ivertex)
      do idof = 1, option%ngeomechdof
        local_disp(ivertex,idof) = &
          geomech_global_aux_vars(ghosted_id)%disp_vector(idof)
        ids(idof + (ivertex-1)*option%ngeomechdof) = &
          (petsc_ids(ivertex)-1)*option%ngeomechdof + (idof-1)
      enddo
      local_temp(ivertex) = &
        temp_loc_p(ghosted_id) - temp_init_loc_p(ghosted_id)
      local_press(ivertex) = &
        press_loc_p(ghosted_id) - press_init_loc_p(ghosted_id)
      if (geomech_parameter%youngs_modulus_spatially_varying) then
        youngs_vec(ivertex) = temp_youngs_modulus_loc_p(ghosted_id)
      else
        youngs_vec(ivertex) = &
          geomech_parameter%youngs_modulus(nint(imech_loc_p(ghosted_id)))
      endif
      if (geomech_parameter%poissons_ratio_spatially_varying) then
        poissons_vec(ivertex) = temp_poissons_ratio_loc_p(ghosted_id)
      else
        poissons_vec(ivertex) = &
          geomech_parameter%poissons_ratio(nint(imech_loc_p(ghosted_id)))
      endif
      if (geomech_parameter%biot_coeff_spatially_varying) then
        beta_vec(ivertex) = temp_biot_coeff_loc_p(ghosted_id)
      else
        beta_vec(ivertex) = &
          geomech_parameter%biot_coeff(nint(imech_loc_p(ghosted_id)))
      endif
      if (geomech_parameter%thermal_exp_coeff_spatially_varying) then
        alpha_vec(ivertex) = temp_thermal_exp_coeff_loc_p(ghosted_id)
      else
        alpha_vec(ivertex) = &
          geomech_parameter%thermal_exp_coeff(nint(imech_loc_p(ghosted_id)))
      endif
    enddo
    call GeomechForceLocalElemStressStrain(size_elenodes, &
       local_coordinates(1:size_elenodes,:), &
       local_disp(1:size_elenodes,:),local_temp(1:size_elenodes), &
       local_press(1:size_elenodes),youngs_vec(1:size_elenodes), &
       poissons_vec(1:size_elenodes),alpha_vec(1:size_elenodes), &
       beta_vec(1:size_elenodes),eletype,grid%gauss_node(ielem)%dim, &
       strain(1:size_elenodes,:),stress(1:size_elenodes,:), &
       stress_total(1:size_elenodes,:),option)

    do ivertex = 1, size_elenodes
      ghosted_id = elenodes(ivertex)
      do idof = 1, SIX_INTEGER
        strain_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) = &
          strain_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) + &
          strain(ivertex,idof)
        stress_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) = &
          stress_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) + &
          stress(ivertex,idof)
        stress_total_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) = &
          stress_total_loc_p(idof + (ghosted_id-1)*SIX_INTEGER) + &
          stress_total(ivertex,idof)
      enddo
    enddo
  enddo

  deallocate(elenodes)
  deallocate(local_coordinates)
  deallocate(local_disp)
  deallocate(local_temp)
  deallocate(local_press)
  deallocate(petsc_ids)
  deallocate(ids)
  deallocate(youngs_vec)
  deallocate(poissons_vec)
  deallocate(alpha_vec)
  deallocate(beta_vec)
  deallocate(strain)
  deallocate(stress)
  deallocate(stress_total)

  call VecRestoreArray(field%temp_loc,temp_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%temp_init_loc,temp_init_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%press_loc,press_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%press_init_loc,press_init_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%strain_loc,strain_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%stress_loc,stress_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%stress_total_loc,stress_total_loc_p,ierr);CHKERRQ(ierr)

  if (geomech_parameter%youngs_modulus_spatially_varying) then
    call VecRestoreArray(field%youngs_modulus_loc,temp_youngs_modulus_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%poissons_ratio_spatially_varying) then
    call VecRestoreArray(field%poissons_ratio_loc,temp_poissons_ratio_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%biot_coeff_spatially_varying) then
    call VecRestoreArray(field%biot_coeff_loc,temp_biot_coeff_loc_p,ierr);CHKERRQ(ierr)
  endif
  if (geomech_parameter%thermal_exp_coeff_spatially_varying) then
    call VecRestoreArray(field%thermal_exp_coeff_loc,temp_thermal_exp_coeff_loc_p,ierr);CHKERRQ(ierr)
  endif

  call GeomechDiscretizationLocalToGlobalAdd(geomech_discretization, &
                                             field%strain_loc,field%strain, &
                                             SIX_INTEGER)
  call GeomechDiscretizationLocalToGlobalAdd(geomech_discretization, &
                                             field%stress_loc,field%stress, &
                                             SIX_INTEGER)
  call GeomechDiscretizationLocalToGlobalAdd(geomech_discretization, &
                                             field%stress_total_loc, &
                                             field%stress_total, &
                                             SIX_INTEGER)

! Now take the average at each node for elements sharing the node
  call VecGetArray(grid%no_elems_sharing_node,no_elems_p, &
                      ierr);CHKERRQ(ierr)
  call VecGetArray(field%strain,strain_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%stress,stress_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%stress_total,stress_total_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax_node
    do idof = 1, SIX_INTEGER
      strain_p(idof + (local_id-1)*SIX_INTEGER) = &
        strain_p(idof + (local_id-1)*SIX_INTEGER)/nint(no_elems_p(local_id))
      stress_p(idof + (local_id-1)*SIX_INTEGER) = &
        stress_p(idof + (local_id-1)*SIX_INTEGER)/nint(no_elems_p(local_id))
      stress_total_p(idof + (local_id-1)*SIX_INTEGER) = &
        stress_total_p(idof + (local_id-1)*SIX_INTEGER)/nint(no_elems_p(local_id))
    enddo
  enddo
  call VecRestoreArray(field%strain,strain_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%stress,stress_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%stress_total,stress_total_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(grid%no_elems_sharing_node,no_elems_p, &
                          ierr);CHKERRQ(ierr)

! Now scatter back to local domains
  call GeomechDiscretizationGlobalToLocal(geomech_discretization, &
                                          field%strain,field%strain_loc, &
                                          SIX_INTEGER)
  call GeomechDiscretizationGlobalToLocal(geomech_discretization, &
                                          field%stress,field%stress_loc, &
                                          SIX_INTEGER)
  call GeomechDiscretizationGlobalToLocal(geomech_discretization, &
                                          field%stress_total, &
                                          field%stress_total_loc, &
                                          SIX_INTEGER)

  call VecGetArray(field%strain_loc,strain_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%stress_loc,stress_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%stress_total_loc,stress_total_loc_p,ierr);CHKERRQ(ierr)
! Copy them to global_aux_vars
  do ghosted_id = 1, grid%ngmax_node
    do idof = 1, SIX_INTEGER
      geomech_global_aux_vars(ghosted_id)%strain(idof) = &
        strain_loc_p(idof + (ghosted_id-1)*SIX_INTEGER)
      geomech_global_aux_vars(ghosted_id)%stress(idof) = &
        stress_loc_p(idof + (ghosted_id-1)*SIX_INTEGER)
      geomech_global_aux_vars(ghosted_id)%stress_total(idof) = &
        stress_total_loc_p(idof + (ghosted_id-1)*SIX_INTEGER)
    enddo
  enddo
  call VecRestoreArray(field%strain_loc,strain_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%stress_loc,stress_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%stress_total_loc,stress_total_loc_p,ierr);CHKERRQ(ierr)

end subroutine GeomechForceStressStrain

! ************************************************************************** !

subroutine GeomechForceLocalElemStressStrain(size_elenodes,local_coordinates, &
                                             local_disp,local_temp, &
                                             local_press, &
                                             local_youngs,local_poissons, &
                                             local_alpha, local_beta, &
                                             eletype,dim,strain,stress, &
                                             stress_total,option)
  !
  ! Computes the stress-strain for a local
  ! element
  !
  ! Author: Satish Karra
  ! Date: 09/17/13
  !

  use Grid_Unstructured_Cell_module
  use Shape_Function_module
  use Option_module
  use Utility_module

  implicit none

  type(shapefunction_type) :: shapefunction
  type(option_type) :: option

  PetscInt,  intent(in) :: size_elenodes, eletype, dim
  PetscReal, intent(in) :: local_coordinates(:,:)      ! (nen,3)
  PetscReal, intent(in) :: local_disp(:,:)             ! (nen,3)
  PetscReal, intent(in) :: local_temp(:), local_press(:)
  PetscReal, intent(in) :: local_youngs(:), local_poissons(:)
  PetscReal, intent(in) :: local_alpha(:), local_beta(:)
  PetscReal, intent(out) :: strain(:,:), stress(:,:), stress_total(:,:) ! (nen,6)

  PetscInt :: a, v
  PetscReal :: J_map(3,3), inv_J_map(3,3), detJ_map
  PetscReal :: dNdx(size_elenodes,3)
  PetscReal :: gradU(3,3), eps(3,3), sig(3,3), sig_total_mat(3,3)
  PetscReal :: lambda, mu, E, nu, alpha, beta, dT, dP
  PetscReal :: trE

  if (dim /= 3) then
    call GeomechForceError('GEOMECHANICS: stress/strain routine expects dim=3.', option)
    return
  end if

  strain = 0.0d0
  stress = 0.0d0
  stress_total = 0.0d0

  shapefunction%element_type = eletype
  call ShapeFunctionInitialize(shapefunction,option)

  do v = 1, size_elenodes
    ! Evaluate at vertex location in reference element
    shapefunction%zeta = shapefunction%coord(v,:)
    call ShapeFunctionCalculate(shapefunction,option)

    J_map = matmul(transpose(local_coordinates), shapefunction%DN)
    call MatInv3WithDet(J_map, inv_J_map, detJ_map)
    if (detJ_map <= 0.0d0) then
      call GeomechForceError('GEOMECHANICS: det(J) must be positive in stress/strain evaluation.', option)
      call ShapeFunctionDestroy(shapefunction)
      return
    end if
    dNdx = matmul(shapefunction%DN, inv_J_map)

    ! Interpolate constitutive coefficients at this evaluation point.
    E     = dot_product(shapefunction%N, local_youngs)
    nu    = dot_product(shapefunction%N, local_poissons)
    alpha = dot_product(shapefunction%N, local_alpha)
    beta  = dot_product(shapefunction%N, local_beta)
    call GeomechGetLambdaMu(lambda, mu, E, nu)

    ! gradU(i,j) = sum_a dNdx(a,i) * u_a(j)
    gradU = 0.0d0
    do a = 1, size_elenodes
      gradU(1,1) = gradU(1,1) + dNdx(a,1)*local_disp(a,1)
      gradU(1,2) = gradU(1,2) + dNdx(a,1)*local_disp(a,2)
      gradU(1,3) = gradU(1,3) + dNdx(a,1)*local_disp(a,3)

      gradU(2,1) = gradU(2,1) + dNdx(a,2)*local_disp(a,1)
      gradU(2,2) = gradU(2,2) + dNdx(a,2)*local_disp(a,2)
      gradU(2,3) = gradU(2,3) + dNdx(a,2)*local_disp(a,3)

      gradU(3,1) = gradU(3,1) + dNdx(a,3)*local_disp(a,1)
      gradU(3,2) = gradU(3,2) + dNdx(a,3)*local_disp(a,2)
      gradU(3,3) = gradU(3,3) + dNdx(a,3)*local_disp(a,3)
    end do

    ! Small-strain tensor and Cauchy stress for isotropic elasticity.
    eps = 0.5d0 * (gradU + transpose(gradU))
    trE = eps(1,1) + eps(2,2) + eps(3,3)

    sig = 2.0d0*mu*eps
    sig(1,1) = sig(1,1) + lambda*trE
    sig(2,2) = sig(2,2) + lambda*trE
    sig(3,3) = sig(3,3) + lambda*trE

    ! Total stress = effective stress - thermal term - Biot pore-pressure term.
    dT = local_temp(v)
    dP = local_press(v)
    sig_total_mat = sig
    sig_total_mat(1,1) = sig_total_mat(1,1) - alpha*dT*(3.0d0*lambda + 2.0d0*mu) - beta*dP
    sig_total_mat(2,2) = sig_total_mat(2,2) - alpha*dT*(3.0d0*lambda + 2.0d0*mu) - beta*dP
    sig_total_mat(3,3) = sig_total_mat(3,3) - alpha*dT*(3.0d0*lambda + 2.0d0*mu) - beta*dP

    ! Store in your 6-component convention (same as your old mapping)
    strain(v,1) = eps(1,1)
    strain(v,2) = eps(2,2)
    strain(v,3) = eps(3,3)
    strain(v,4) = eps(1,2)
    strain(v,5) = eps(2,3)
    strain(v,6) = eps(3,1)

    stress(v,1) = sig(1,1)
    stress(v,2) = sig(2,2)
    stress(v,3) = sig(3,3)
    stress(v,4) = sig(1,2)
    stress(v,5) = sig(2,3)
    stress(v,6) = sig(3,1)

    stress_total(v,1) = sig_total_mat(1,1)
    stress_total(v,2) = sig_total_mat(2,2)
    stress_total(v,3) = sig_total_mat(3,3)
    stress_total(v,4) = sig_total_mat(1,2)
    stress_total(v,5) = sig_total_mat(2,3)
    stress_total(v,6) = sig_total_mat(3,1)
  end do

  call ShapeFunctionDestroy(shapefunction)

end subroutine GeomechForceLocalElemStressStrain

! ************************************************************************** !

subroutine GeomechUpdateSolution(geomech_realization)
  !
  ! Updates data in module after a successful time
  ! step
  !
  ! Author: Satish Karra, LANL
  ! Date: 09/17/13
  !

  use Geomechanics_Realization_class
  use Geomechanics_Field_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  type(geomech_field_type), pointer :: field

  field => geomech_realization%geomech_field

  call GeomechForceUpdateAuxVars(geomech_realization)
  call GeomechUpdateSolutionPatch(geomech_realization)

end subroutine GeomechUpdateSolution

! ************************************************************************** !

subroutine GeomechUpdateSolutionPatch(geomech_realization)
  !
  ! updates data in module after a successful time
  ! step
  !
  ! Author: satish karra, lanl
  ! Date: 09/17/13
  !

  use Geomechanics_Realization_class

  implicit none

  class(realization_geomech_type) :: geomech_realization

  call GeomechForceStressStrain(geomech_realization)

end subroutine GeomechUpdateSolutionPatch

! ************************************************************************** !

subroutine GeomechStoreInitialPressTemp(geomech_realization)
  !
  ! Stores initial pressure and temperature from
  ! subsurface
  !
  ! Author: Satish Karra, LANL
  ! Date: 09/24/13
  !

  use Geomechanics_Realization_class

  implicit none

  class(realization_geomech_type) :: geomech_realization

  PetscErrorCode :: ierr

  call VecCopy(geomech_realization%geomech_field%press_loc, &
               geomech_realization%geomech_field%press_init_loc, &
               ierr);CHKERRQ(ierr)

  call VecCopy(geomech_realization%geomech_field%temp_loc, &
               geomech_realization%geomech_field%temp_init_loc, &
               ierr);CHKERRQ(ierr)

  call VecCopy(geomech_realization%geomech_field%fluid_density_loc, &
               geomech_realization%geomech_field%fluid_density_init_loc, &
               ierr);CHKERRQ(ierr)

end subroutine GeomechStoreInitialPressTemp

! ************************************************************************** !

subroutine GeomechStoreInitialPorosity(realization,geomech_realization)
  !
  ! Stores initial porosity from
  ! subsurface
  !
  ! Author: Satish Karra, LANL
  ! Date: 10/22/13
  !

  use Geomechanics_Realization_class
  use Realization_Subsurface_class
  use Discretization_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  class(realization_subsurface_type) :: realization
  type(discretization_type) :: discretization

  call DiscretizationDuplicateVector(discretization, &
                                     realization%field%work_loc, &
                                     geomech_realization%geomech_field% &
                                     porosity_init_loc)

end subroutine GeomechStoreInitialPorosity

! ************************************************************************** !

subroutine GeomechMapMaterialIdsToFlow(realization,geomech_realization)
  !
  ! Maps geomechanics material ids from geomech nodes onto the
  ! corresponding flow cells via the flow-geomech scatter. Each flow
  ! cell stores the geomech material id on material_auxvar%secondary_material_id
  ! so porosity updates (fixed-stress) use geomech imat
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 07/16/26
  !
  use Realization_Subsurface_class
  use Geomechanics_Realization_class
  use Geomechanics_Discretization_module
  use Geomechanics_Field_module
  use Geomechanics_Grid_Aux_module
  use GMDM_Pointer_module
  use Grid_module
  use Material_Aux_module
  use Option_module
  use Option_Geomechanics_module
  use Patch_module
  use Discretization_module

  implicit none

  class(realization_subsurface_type) :: realization
  class(realization_geomech_type) :: geomech_realization

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(geomech_field_type), pointer :: geomech_field
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscReal, pointer :: subsurf_p(:)
  PetscReal, pointer :: flow_work_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: mat_id
  PetscErrorCode :: ierr

  option => realization%option
  if (option%geomechanics%flow_coupling /= GEOMECH_TWO_WAY_COUPLED) return
  if (option%geomechanics%split_scheme /= GEOMECH_FIXED_STRESS_SPLIT) return

  grid => realization%discretization%grid
  patch => realization%patch
  material_auxvars => patch%aux%Material%auxvars
  geomech_field => geomech_realization%geomech_field
  geomech_discretization => geomech_realization%geomech_discretization

  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   ONEDOF)

  ! Promote local imech ids to global work
  call GeomechDiscretizationLocalToGlobal(geomech_discretization, &
                                          geomech_field%imech_loc, &
                                          geomech_field%work,ONEDOF)

  ! Reverse of subsurf->geomech scatter maps geomech node data onto flow cells
  call VecScatterBegin(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                       geomech_field%work, &
                       geomech_field%subsurf_vec_1dof, &
                       INSERT_VALUES,SCATTER_REVERSE, &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                     geomech_field%work, &
                     geomech_field%subsurf_vec_1dof, &
                     INSERT_VALUES,SCATTER_REVERSE, &
                     ierr);CHKERRQ(ierr)

  ! Promote to ghosted material auxvars via flow work / work_loc
  call VecGetArray(geomech_field%subsurf_vec_1dof,subsurf_p, &
                      ierr);CHKERRQ(ierr)
  call VecGetArray(realization%field%work,flow_work_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    flow_work_p(local_id) = subsurf_p(local_id)
  enddo
  call VecRestoreArray(geomech_field%subsurf_vec_1dof,subsurf_p, &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArray(realization%field%work,flow_work_p,ierr);CHKERRQ(ierr)

  call DiscretizationGlobalToLocal(realization%discretization, &
                                   realization%field%work, &
                                   realization%field%work_loc,ONEDOF)

  call VecGetArray(realization%field%work_loc,flow_work_p,ierr);CHKERRQ(ierr)
  do ghosted_id = 1, grid%ngmax
    if (patch%imat(ghosted_id) <= 0) cycle
    mat_id = nint(flow_work_p(ghosted_id))
    if (mat_id <= 0) then
      write(option%io_buffer,*) grid%nG2A(ghosted_id)
      option%io_buffer = 'No geomechanics material id mapped to flow &
        &cell id ' // trim(adjustl(option%io_buffer)) // &
        ' in GeomechMapMaterialIdsToFlow.'
      call PrintErrMsgByRank(option)
    endif
    material_auxvars(ghosted_id)%secondary_material_id = mat_id
  enddo
  call VecRestoreArray(realization%field%work_loc,flow_work_p, &
                          ierr);CHKERRQ(ierr)

end subroutine GeomechMapMaterialIdsToFlow

! ************************************************************************** !

subroutine GeomechStoreInitialDisp(geomech_realization)
  !
  ! Stores initial displacement for calculating
  ! relative displacements
  !
  ! Author: Satish Karra, LANL
  ! Date: 09/30/13
  !

  use Geomechanics_Realization_class

  implicit none

  class(realization_geomech_type) :: geomech_realization

  PetscErrorCode :: ierr

  call VecCopy(geomech_realization%geomech_field%disp_xx_loc, &
               geomech_realization%geomech_field%disp_xx_init_loc, &
               ierr);CHKERRQ(ierr)

end subroutine GeomechStoreInitialDisp

end module Geomechanics_Force_module
