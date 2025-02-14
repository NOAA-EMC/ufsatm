! ###########################################################################################
!> \file MPAS_typedefs.F90
! ########################################################################################### 
module MPAS_typedefs

  use mpi_f08
  use machine, only: kind_phys, kind_dbl_prec, kind_sngl_prec
  implicit none

!> \section arg_table_MPAS_typedefs
!! \htmlinclude MPAS_typedefs.html
!!

  ! MPAS_init_type           !<
  ! MPAS_statein_type        !< prognostic state data in from dycore
  ! MPAS_stateint_type       !< prognostic state or tendencies before call MP.
  ! MPAS_stateout_type       !< prognostic state or tendencies return to dycore

  ! #########################################################################################
  ! MPAS_init_type
  ! #########################################################################################
!! \section arg_table_MPAS_init_type
!! \htmlinclude MPAS_init_type.html
!!
  type MPAS_init_type
  end type MPAS_init_type

  ! #########################################################################################
  ! MPAS_statein_type
  !  Prognostic state variables with layer and level specific data from dycore.
  ! #########################################################################################
!! \section arg_table_MPAS_statein_type
!! \htmlinclude MPAS_statein_type.html
!!
  type MPAS_statein_type
  end type MPAS_statein_type
  
  ! #########################################################################################
  ! MPAS_stateint_type
  !  Prognostic state or tendencies BEFORE calling microphysics.
  ! #########################################################################################
!! \section arg_table_MPAS_stateint_type
!! \htmlinclude MPAS_stateint_type.html
!!
  type MPAS_stateint_type
  end type MPAS_stateint_type
  
  ! #########################################################################################
  ! MPAS_stateout_type
  !  Prognostic state or tendencies after ALL physical parameterizations.
  ! #########################################################################################
!! \section arg_table_MPAS_stateout_type
!! \htmlinclude MPAS_stateout_type.html
!!
  type MPAS_stateout_type
  end type MPAS_stateout_type

  public MPAS_init_type, MPAS_statein_type, MPAS_stateint_type, MPAS_stateout_type

end module MPAS_typedefs
