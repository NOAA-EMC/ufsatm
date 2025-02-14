! ###########################################################################################
!> \file ufs_mpas_subdriver.F90
! ###########################################################################################
module ufs_mpas_subdriver
  use mpas_kind_types,    only: mpas_real_kind => RKIND
  use mpas_derived_types, only: core_type, domain_type, MPAS_Clock_type
  use mpi_f08,            only: MPI_Comm
  use atm_core_interface

  implicit none

  public :: corelist, domain_ptr
  public :: mpas_init_phase1
  private

  type(core_type),       pointer :: corelist   => null()
  type(domain_type),     pointer :: domain_ptr => null()
  type(MPAS_Clock_type), pointer :: clock      => null()

contains
  ! #########################################################################################
  ! 
  ! #########################################################################################
  subroutine mpas_init_phase1(mpicomm, logUnits)
    use mpas_domain_routines, only : mpas_allocate_domain
    use mpas_framework,       only : mpas_framework_init_phase1
    use atm_core_interface,   only : atm_setup_core, atm_setup_domain

    ! Inputs
    type(MPI_Comm), intent(in) :: mpicomm
    integer, dimension(2), intent(in) :: logUnits
    ! Locals
    integer :: ierr

    allocate(corelist, stat=ierr)
    if( ierr /= 0 ) stop
    nullify(corelist%next)

    allocate(corelist%domainlist, stat=ierr)
    if( ierr /= 0 ) stop
    nullify(corelist%domainlist%next)

    domain_ptr => corelist%domainlist
    domain_ptr%core => corelist

    call mpas_allocate_domain(domain_ptr)
    domain_ptr%domainID = 0

    ! Initialize MPAS infrastructure
    call mpas_framework_init_phase1(domain_ptr%dminfo, external_comm=mpicomm)
    call atm_setup_core(corelist)
    call atm_setup_domain(domain_ptr)
    
  end subroutine mpas_init_phase1
 
end module ufs_mpas_subdriver
