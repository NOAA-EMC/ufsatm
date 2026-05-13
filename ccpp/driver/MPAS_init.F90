! ###########################################################################################
!> \file MPAS_init.F90
!>
! ###########################################################################################
module MPAS_init
  use machine,            only : kind_phys
  use ufs_mpas_subdriver, only : MPAS_control_type
  use GFS_typedefs,       only : GFS_control_type, GFS_diag_type, GFS_grid_type, GFS_tbd_type
  use GFS_typedefs,       only : GFS_sfcprop_type, GFS_statein_type, GFS_stateout_type, GFS_cldprop_type
  use GFS_typedefs,       only : GFS_radtend_type
  use GFS_typedefs,       only : GFS_coupling_type

  implicit none

  public :: MPAS_initialize

contains
  !> #########################################################################################
  !> Procedure to initialize MPAS interface to CCPP Physics.
  !>
  !> #########################################################################################
  subroutine MPAS_initialize (Model, Diag, Grid, Tbd, SfcProp, Statein, Stateout, CldProp, RadTend,    &
                              Coupling, Init_parm)
#ifdef _OPENMP
    use omp_lib
#endif

    ! Inputs
    type(GFS_control_type),      intent(inout) :: Model
    type(GFS_diag_type),         intent(inout) :: Diag
    type(GFS_grid_type),         intent(inout) :: Grid
    type(GFS_tbd_type),          intent(inout) :: Tbd
    type(GFS_sfcprop_type),      intent(inout) :: SfcProp
    type(GFS_statein_type),      intent(inout) :: Statein
    type(GFS_stateout_type),     intent(inout) :: Stateout
    type(GFS_cldprop_type),      intent(inout) :: Cldprop
    type(GFS_radtend_type),      intent(inout) :: Radtend
    type(GFS_coupling_type),     intent(inout) :: Coupling
    type(MPAS_control_type),     intent(inout) :: Init_parm
    
    ! Locals
    integer :: nb
    integer :: nblks
    integer :: nt
    integer :: nthrds
    logical :: non_uniform_blocks
    integer :: ix

    nblks = size(Init_parm%blksz)

#ifdef _OPENMP
    nthrds = omp_get_max_threads()
#else
    nthrds = 1
#endif

    ! Set control properties (including physics namelist read)
    Model%dycore_active = Model%dycore_mpas
    call Model%init(Init_parm%nlunit, Init_parm%fn_nml, Init_parm%me, Init_parm%master,      &
         Init_parm%logunit, Init_parm%levs, real(Init_parm%dt_dycore, kind_phys),            &
         real(Init_parm%dt_phys, kind_phys), Init_parm%iau_offset, Init_parm%bdat,           &
         Init_parm%cdat, Init_parm%nwat, Init_parm%tracer_names, Init_parm%tracer_types,     &
         Init_parm%input_nml_file, Init_parm%blksz, Init_parm%restart, Init_parm%mpi_comm,   &
         Init_parm%fcst_ntasks, nthrds)

    ! Allocate data containers for physics.
    call Grid%create(Model)
    call Diag%create(Model)
    call Tbd%create(Model)
    call SfcProp%create(Model)
    call Statein%create(Model)
    call Stateout%create(Model)
    call Cldprop%create(Model)
    call Radtend%create(Model)
    call Coupling%create(Model)
    
  end subroutine MPAS_initialize

end module MPAS_init
