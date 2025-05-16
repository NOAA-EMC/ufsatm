! ###########################################################################################
!> \file MPAS_init.F90
!>
! ###########################################################################################
module MPAS_init
  use machine,       only: kind_phys
  use MPAS_typedefs, only: MPAS_statein_type, MPAS_stateint_type, MPAS_stateout_type
  use MPAS_typedefs, only: MPAS_control_type
  use GFS_typedefs,  only: GFS_control_type, GFS_diag_type

  implicit none

  public  MPAS_initialize

contains
  !> #########################################################################################
  !> Procedure to initialize MPAS interface to CCPP Physics.
  !>
  !> #########################################################################################
  subroutine MPAS_initialize (Model, Diag, Statein, Stateint, Stateout, MPAS)
#ifdef _OPENMP
    use omp_lib
#endif

    ! Inputs
    type(GFS_control_type),      intent(inout) :: Model
    type(GFS_diag_type),         intent(inout) :: Diag
    type(MPAS_statein_type),     intent(inout) :: Statein
    type(MPAS_stateint_type),    intent(inout) :: Stateint
    type(MPAS_stateout_type),    intent(inout) :: Stateout
    type(MPAS_control_type),     intent(in   ) :: MPAS

    ! Locals
    integer :: nb
    integer :: nblks
    integer :: nt
    integer :: nthrds
    logical :: non_uniform_blocks
    integer :: ix

    nblks = size(MPAS%blksz)

#ifdef _OPENMP
    nthrds = omp_get_max_threads()
#else
    nthrds = 1
#endif
    !--- set control properties (including namelist read)
    call Model%init(MPAS%nlunit, MPAS%fn_nml, MPAS%me, MPAS%master, MPAS%logunit, &
         MPAS%levs,  MPAS%gnx, MPAS%gny,     &
         MPAS%dt_dycore, MPAS%dt_phys,  MPAS%iau_offset, MPAS%bdat, MPAS%cdat, MPAS%nwat,   &
         MPAS%tracer_names, MPAS%tracer_types, MPAS%input_nml_file, MPAS%tile_num,          &
         MPAS%blksz, MPAS%restart, MPAS%hydrostatic, MPAS%mpi_comm, MPAS%fcst_ntasks, nthrds)

    call Diag%create(Model)

  end subroutine MPAS_initialize

end module MPAS_init
