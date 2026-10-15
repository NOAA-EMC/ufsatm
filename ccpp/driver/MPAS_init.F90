! ###########################################################################################
!> \file MPAS_init.F90
!>
! ###########################################################################################
module MPAS_init
  use machine,            only : kind_phys
  use GFS_typedefs,       only : GFS_control_type, GFS_diag_type, GFS_grid_type, GFS_tbd_type
  use GFS_typedefs,       only : GFS_sfcprop_type, GFS_statein_type, GFS_stateout_type, GFS_cldprop_type
  use GFS_typedefs,       only : GFS_radtend_type
  use GFS_typedefs,       only : GFS_coupling_type
  use mpi_f08

  implicit none

  public :: MPAS_initialize

contains
  !> #########################################################################################
  !> Procedure to initialize MPAS interface to CCPP Physics.
  !>
  !> #########################################################################################
  subroutine MPAS_initialize (Model, Diag, Grid, Tbd, SfcProp, Statein, Stateout, CldProp,   &
       RadTend, Coupling, me, master, mpicomm, levs, dt_dyn, dt_phys, nml_funit,             &
       nml_filename, bdat, cdat, nwat, fcst_ntasks, blksz, input_nml_file, constituent_name, &
       constituent_type, restart, gnx, gny, p_ref)
#ifdef _OPENMP
    use omp_lib
#endif
    ! Inputs
    integer,                     intent(in   ) :: me
    integer,                     intent(in   ) :: master
    integer,                     intent(in   ) :: levs
    integer,                     intent(in   ) :: dt_dyn
    integer,                     intent(in   ) :: dt_phys
    integer,                     intent(in   ) :: nml_funit
    integer,                     intent(in   ) :: bdat(8)
    integer,                     intent(in   ) :: cdat(8)
    integer,                     intent(in   ) :: nwat
    integer,                     intent(in   ) :: fcst_ntasks
    integer,                     intent(in   ) :: blksz(:)
    character(len=*),            intent(in   ) :: nml_filename
    type(MPI_Comm),              intent(in   ) :: mpicomm
    logical,                     intent(in   ) :: restart
    character(len=:), pointer,   intent(in   ) :: input_nml_file(:)
    character(len=*),            intent(in   ) :: constituent_name(:)
    integer,                     intent(in   ) :: constituent_type(:)
    ! Equivalent Gaussian-grid point count (lonr/latr) and dycore-neutral reference
    ! pressure profile, needed by UGWPv1 under a dycore that has no FV3-style hybrid
    ! sigma-pressure coordinate to derive one from (see
    ! atmos_coupling.F90::ufs_mpas_reference_pressure and atmos_model.F90::atmos_model_init).
    ! Optional: absent for cores/suites that don't need them, in which case
    ! control_initialize's MPAS branch leaves Model%lonr/latr/p_ref at their existing values.
    integer,           optional, intent(in   ) :: gnx
    integer,           optional, intent(in   ) :: gny
    real(kind_phys),   optional, intent(in   ) :: p_ref(:)
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

    ! Locals
    integer :: nb
    integer :: nblks
    integer :: nt
    integer :: nthrds
    integer :: ix


#ifdef _OPENMP
    nthrds = omp_get_max_threads()
#else
    nthrds = 1
#endif

    ! Set control properties (including physics namelist read)
    Model%dycore_active = Model%dycore_mpas
    call Model%init(nml_funit, nml_filename, me, master, 0, levs, real(dt_dyn, kind_phys),   &
         real(dt_phys, kind_phys), 0, bdat, cdat, nwat, constituent_name, constituent_type,  &
         input_nml_file, blksz, restart, mpicomm, fcst_ntasks, nthrds,                       &
         gnx=gnx, gny=gny, p_ref=p_ref)

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
