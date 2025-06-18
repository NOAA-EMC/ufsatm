! ###########################################################################################
!> \file MPAS_init.F90
!>
! ###########################################################################################
module MPAS_init
  use machine,            only : kind_phys
  use ufs_mpas_subdriver, only : MPAS_control_type
  use GFS_typedefs,       only : GFS_control_type, GFS_diag_type, GFS_grid_type, GFS_tbd_type
  use GFS_typedefs,       only : GFS_sfcprop_type, GFS_statein_type, GFS_cldprop_type
  use GFS_typedefs,       only : GFS_radtend_type
  use GFS_typedefs,       only : GFS_coupling_type
  use CCPP_typedefs,      only : GFS_interstitial_type

  implicit none

  public :: MPAS_initialize

contains
  !> #########################################################################################
  !> Procedure to initialize MPAS interface to CCPP Physics.
  !>
  !> #########################################################################################
  subroutine MPAS_initialize (Model, Diag, Grid, Tbd, SfcProp, Statein, CldProp, RadTend,    &
                              Coupling, MPAS, Interstitial)
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
    type(GFS_cldprop_type),      intent(inout) :: Cldprop
    type(GFS_radtend_type),      intent(inout) :: Radtend
    type(GFS_coupling_type),     intent(inout) :: Coupling
    type(MPAS_control_type),     intent(inout) :: MPAS
    type(GFS_interstitial_type), intent(inout) :: Interstitial(:)
    
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

    ! Set control properties (including physics namelist read)
    call Model%init('MPAS', MPAS%nlunit, MPAS%fn_nml, MPAS%me, MPAS%master, MPAS%logunit,    &
         MPAS%levs, real(MPAS%dt_dycore, kind_phys), real(MPAS%dt_phys, kind_phys),          &
         MPAS%iau_offset, MPAS%bdat, MPAS%cdat, MPAS%nwat, MPAS%tracer_names,                &
         MPAS%tracer_types, MPAS%input_nml_file, MPAS%blksz, MPAS%restart, MPAS%mpi_comm,    &
         MPAS%fcst_ntasks, nthrds)

    ! Allocate data containers for physics.
    call Grid%create(Model)
    call Diag%create(Model)
    call Tbd%create(Model)
    call SfcProp%create(Model)
    call Statein%create(Model)
    call Cldprop%create(Model)
    call Radtend%create(Model)
    call Coupling%create(Model)

    ! This logic deals with non-uniform block sizes for CCPP. When non-uniform block sizes
    ! are used, it is required that only the last block has a different (smaller) size than
    ! all other blocks. This is the standard in FV3. If this is the case, set non_uniform_blocks
    ! to .true. and initialize nthreads+1 elements of the interstitial array. The extra element
    ! will be used by the thread that runs over the last, smaller block.
    if (minval(MPAS%blksz)==maxval(MPAS%blksz)) then
       non_uniform_blocks = .false.
    elseif (all(minloc(MPAS%blksz)==(/size(MPAS%blksz)/))) then
       non_uniform_blocks = .true.
    else
       write(0,'(2a)') 'For non-uniform blocksizes, only the last element ', &
                       'in MPAS%blksz can be different from the others'
       stop
    endif

    ! Initialize the Interstitial data type in parallel so that
    ! each thread creates (touches) its Interstitial(nt) first.
    !$OMP parallel do default (shared) &
    !$OMP            schedule (static,1) &
    !$OMP            private  (nt)
    do nt=1,nthrds
       call Interstitial(nt)%create(maxval(MPAS%blksz), Model)
    enddo
    !$OMP end parallel do

    if (non_uniform_blocks) then
       call Interstitial(nthrds+1)%create(MPAS%blksz(nblks), Model)
    end if
    
  end subroutine MPAS_initialize

end module MPAS_init
