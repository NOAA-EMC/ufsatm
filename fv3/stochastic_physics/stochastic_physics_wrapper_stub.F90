module stochastic_physics_wrapper_mod

  use GFS_typedefs,      only: GFS_control_type, GFS_statein_type, GFS_grid_type, GFS_sfcprop_type, GFS_radtend_type, GFS_coupling_type
  use block_control_mod, only: block_control_type

  implicit none

  public stochastic_physics_wrapper
  public stochastic_physics_wrapper_end

contains

  subroutine stochastic_physics_wrapper(GFS_Control, GFS_Statein, GFS_Grid, GFS_Sfcprop, GFS_Radtend, GFS_Coupling, Atm_block, ierr)

    type(GFS_control_type),   intent(inout) :: GFS_Control
    type(GFS_statein_type),   intent(in)    :: GFS_Statein
    type(GFS_grid_type),      intent(in)    :: GFS_Grid
    type(GFS_sfcprop_type),   intent(inout) :: GFS_Sfcprop
    type(GFS_radtend_type),   intent(inout) :: GFS_Radtend
    type(GFS_coupling_type),  intent(inout) :: GFS_Coupling
    type(block_control_type), intent(inout) :: Atm_block
    integer,                  intent(out)   :: ierr

    ierr = 0

  end subroutine stochastic_physics_wrapper

  subroutine stochastic_physics_wrapper_end(GFS_Control)

    type(GFS_control_type), intent(inout) :: GFS_Control

  end subroutine stochastic_physics_wrapper_end

end module stochastic_physics_wrapper_mod
