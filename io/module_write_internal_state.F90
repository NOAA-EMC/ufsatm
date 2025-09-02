!> @file
!> @brief The internal state of the write component.
!> @author Jun Wang @date Feb, 2017

!> @brief The internal state of the write component.
!> @author Jun Wang @date Feb, 2017
module write_internal_state
      use esmf
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type output_grid_info
        integer :: im !< Output grid global I dimension size.
        integer :: jm !< Output grid global J dimension size.
        integer :: lm !< Output grid global L dimension size.
        integer :: i_start !< Output grid lower bound of I dimension on current PE.
        integer :: i_end !< Output grid upper bound of I dimension on current PE.
        integer :: j_start !< Output grid lower bound of J dimension on current PE.
        integer :: j_end !< Output grid upper bound of J dimension on current PE.
        real,dimension(:,:),allocatable  :: lonPtr !< Output grid longitudes.
        real,dimension(:,:),allocatable  :: latPtr !< Output grid latitudes.
        integer,dimension(:),allocatable :: i_start_wrtgrp !< I dimension lower bound of all wrire groups.
        integer,dimension(:),allocatable :: i_end_wrtgrp !< I dimension upper bound of all wrire groups.
        integer,dimension(:),allocatable :: j_start_wrtgrp !< J dimension lower bound of all wrire groups.
        integer,dimension(:),allocatable :: j_end_wrtgrp !< J dimension upper bound of all wrire groups.
        real    :: latse !< Output grid South East corner latitude.
        real    :: latnw !< Output grid North West corner latitude.
        real    :: lonse !< Output grid South East corner longitude.
        real    :: lonnw !< Output grid North West corner longitude.
        real    :: latstart !< Output grid start latitude.
        real    :: latlast !< Output grid last latitude.
        real    :: lonstart !< Output grid start logitude.
        real    :: lonlast !< Output grid last longitude.
      end type output_grid_info

      type wrt_internal_state

!--------------------------------
! pe information and task layout
!--------------------------------
!
      integer :: mype !< MPI rank.
      integer :: petcount !< Number of PEs.
!
!--------------------
!*** grid information
!--------------------
      type(esmf_grid) :: wrtgrid !< ESMF output grid.

      !> Array of output_grid_info for all domains.
      type(output_grid_info) ,dimension(:), allocatable :: out_grid_info 
!
!--------------------------
!*** file bundle for output
!--------------------------
      integer :: FBCount !< Numebr of output ESMF field bundles.
!
!-----------------------------------------------------------------------
!***  THE OUTPUT FILE
!-----------------------------------------------------------------------
!
      integer                 :: num_files !< Number of output files.
!
!-----------------------------------------------------------------------
!***  THE OUTPUT FILE
!-----------------------------------------------------------------------
!
      !> ESMF write field bundles.      
      type(ESMF_FieldBundle),dimension(:),allocatable  :: wrtFB 
!
!-------------------------------------
!***  Times used in history filenames
!-------------------------------------
!
      type(ESMF_Time)         :: io_basetime !< ESMF clock's starting time.
      integer                 :: idate(7) !< Forecast initial time.
      integer                 :: fdate(7) !< Forecast current time.
!
!-----------------------------------------
!***  I/O direction flags (Read or Write)
!-----------------------------------------
!
      logical :: output_history !< True if history output is requested.
!
!-----------------------------------------
!***  POST flags and required variables
!-----------------------------------------
!
      logical                  :: write_dopost !< True if inline post is requested.
      character(80)            :: post_namelist !< File name of the inline post namelist.
!
      real(4)                  :: fhzero !< Hours between clearing of diagnostic buckets.
      integer                  :: ntrac !< Number of tracers.
      integer                  :: ncld !< Number of hydrometeors.
      integer                  :: nsoil !< Number of soil layers.
      integer                  :: imp_physics !< Choice of microphysics scheme.
      integer                  :: landsfcmdl !< Choice of land surface model.
      integer                  :: dtp !< Physics timestep.
      real,dimension(:),allocatable :: ak !< a parameter for sigma pressure level calculations.
      real,dimension(:),allocatable :: bk !< b parameter for sigma pressure level calculations.
!-----------------------------------------------------------------------
!
      end type wrt_internal_state

      !> This state is supported by c pointers but not f90 pointers
      !> therefore we need this wrap.
      type write_wrap
        !> Write grid component internal state.         
        type(wrt_internal_state),pointer :: write_int_state  
      end type write_wrap

!-----------------------------------------------------------
!
      end module write_internal_state
!
!-----------------------------------------------------------
