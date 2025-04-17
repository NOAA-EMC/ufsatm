! ###########################################################################################
!> \file MPAS_init.F90
!>
! ###########################################################################################
module MPAS_init
  use machine,       only: kind_phys
  use MPAS_typedefs, only: MPAS_statein_type, MPAS_stateint_type, MPAS_stateout_type
  use MPAS_typedefs, only: MPAS_init_type
  use GFS_typedefs,  only: GFS_control_type, GFS_diag_type

  implicit none

  public  MPAS_initialize

contains
  ! #########################################################################################
  ! MPAS_initialize
  ! #########################################################################################
  subroutine MPAS_initialize (Model, Diag, Statein, Stateint, Stateout, Init)
#ifdef _OPENMP
    use omp_lib
#endif

    ! Inputs
    type(GFS_control_type),      intent(inout) :: Model
    type(GFS_diag_type),         intent(inout) :: Diag
    type(MPAS_statein_type),     intent(inout) :: Statein
    type(MPAS_stateint_type),    intent(inout) :: Stateint
    type(MPAS_stateout_type),    intent(inout) :: Stateout
    type(MPAS_init_type),        intent(in   ) :: Init

    ! Locals
    integer :: nb
    integer :: nblks
    integer :: nt
    integer :: nthrds
    logical :: non_uniform_blocks
    integer :: ix

    nblks = size(Init%blksz)

#ifdef _OPENMP
    nthrds = omp_get_max_threads()
#else
    nthrds = 1
#endif
    !--- set control properties (including namelist read)
    call Model%init(Init%nlunit, Init%fn_nml, Init%me, Init%master, Init%logunit, Init%isc, &
         Init%jsc, Init%nx, Init%ny, Init%levs, Init%cnx, Init%cny, Init%gnx, Init%gny,     &
         Init%dt_dycore, Init%dt_phys,  Init%iau_offset, Init%bdat, Init%cdat, Init%nwat,   &
         Init%tracer_names, Init%tracer_types, Init%input_nml_file, Init%tile_num,          &
         Init%blksz, Init%ak, Init%bk, Init%restart, Init%hydrostatic, Init%mpi_comm,       &
         Init%fcst_ntasks, nthrds)

    call Statein%create(Model)
    call Stateout%create(Model)
    call Stateint%create(Model)
    call Diag%create(Model)

  end subroutine MPAS_initialize

end module MPAS_init
