module MPAS_data

!! \section arg_table_MPAS_data Argument Table
!! \htmlinclude MPAS_data.html
!!
  use MPAS_typedefs, only : MPAS_statein_type
  use MPAS_typedefs, only : MPAS_stateout_type
  implicit none

  private
  
  public MPAS_statein
  public MPAS_stateout

  !------------------------------------------------------!
  ! MPAS data containers.
  !------------------------------------------------------!
  type(MPAS_statein_type),   save, target :: MPAS_statein
  type(MPAS_stateout_type),  save, target :: MPAS_stateout

end module MPAS_data
