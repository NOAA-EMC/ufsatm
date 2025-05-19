!> \file fv3atm_oro_io.F90
!! This file defines routines to read orography files for the fv3atm.
module fv3atm_oro_io

  use block_control_mod,  only: block_control_type
  use fms2_io_mod,        only: FmsNetcdfDomainFile_t, &
       register_axis, register_restart_field
  use fv3atm_common_io,   only: get_nx_ny_from_atm
  use GFS_typedefs,       only: GFS_sfcprop_type, GFS_control_type, kind_phys

  implicit none
  private

  public :: Oro_io_data_type, Oro_io_register, Oro_io_copy, Oro_io_final
  public :: Oro_scale_io_data_type, Oro_scale_io_register, Oro_scale_io_copy, Oro_scale_io_final

  !>\defgroup fv3atm_oro_io FV3ATM Orography I/O Module
  !> @{
  !>@ Storage of working arrays for reading orography data.
  type Oro_io_data_type
    character(len=32),    pointer, private, dimension(:)       :: name2 => null()
    real(kind=kind_phys), pointer, private, dimension(:,:,:)   :: var2  => null()
    real(kind=kind_phys), pointer, private, dimension(:,:,:)   :: var3v => null()
    real(kind=kind_phys), pointer, private, dimension(:,:,:)   :: var3s => null()
  contains
    procedure, public :: register => Oro_io_register
    procedure, public :: copy => Oro_io_copy
    final :: Oro_io_final
  end type Oro_io_data_type

  !>@ Storage of working arrays for reading large-scale and small-scale orography data for gravity wave drag schemes.
  type Oro_scale_io_data_type
    character(len=32),    pointer, private, dimension(:)       :: name => null()
    real(kind=kind_phys), pointer, private, dimension(:,:,:)   :: var  => null()
  contains
    procedure, public :: register => Oro_scale_io_register
    procedure, public :: copy => Oro_scale_io_copy
    final :: Oro_scale_io_final
  end type Oro_scale_io_data_type

  !>@ Number of two-dimensional orography fields (excluding large- and small-scale)
  integer, parameter :: nvar_oro_2d  = 19

  !>@ Number of large-scale and small-scale orography fields
  integer, parameter :: nvar_oro_scale = 10

contains

  !>@brief Registers axes and fields for non-quilt restart reading of non-scaled orography variables.
  !> \section oro_io_data_type%register procedure
  !! Calls FMS restart register functions for axes and
  !! variables in the non-scaled orography data. The scaled data is
  !! handled by another function. This includes both 2D and 3D fields.
  subroutine Oro_io_register(oro, Model, Oro_restart, Atm_block)
    implicit none
    class(Oro_io_data_type) :: oro
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Oro_restart
    type(block_control_type), intent(in) :: Atm_block

    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_fr => NULL()
    integer :: nx, ny

    integer :: nvar_vegfr, nvar_soilfr, n, num

    call get_nx_ny_from_atm(Atm_block, nx, ny)

    ! This #define reduces code length by a lot
#define WARN_DISASSOCIATE(name) \
    if(associated(name)) then ; \
      write(0,*) 'Internal error. Called oro%register twice. Will try to keep going anyway.' ; \
      deallocate(name); \
      nullify(name) ; \
    endif

    WARN_DISASSOCIATE(oro%name2)
    WARN_DISASSOCIATE(oro%var2)
    WARN_DISASSOCIATE(oro%var3v)
    WARN_DISASSOCIATE(oro%var3s)
#undef WARN_DISASSOCIATE

    nvar_vegfr  = Model%nvegcat
    nvar_soilfr = Model%nsoilcat

    allocate(oro%name2(nvar_oro_2d))
    allocate(oro%var2(nx,ny,nvar_oro_2d))

    allocate(oro%var3v(nx,ny,nvar_vegfr))
    allocate(oro%var3s(nx,ny,nvar_soilfr))

    oro%var2 = -9999._kind_phys

    num = 1       ; oro%name2(num)  = 'stddev'     ! hprime(ix,1)
    num = num + 1 ; oro%name2(num)  = 'convexity'  ! hprime(ix,2)
    num = num + 1 ; oro%name2(num)  = 'oa1'        ! hprime(ix,3)
    num = num + 1 ; oro%name2(num)  = 'oa2'        ! hprime(ix,4)
    num = num + 1 ; oro%name2(num)  = 'oa3'        ! hprime(ix,5)
    num = num + 1 ; oro%name2(num)  = 'oa4'        ! hprime(ix,6)
    num = num + 1 ; oro%name2(num)  = 'ol1'        ! hprime(ix,7)
    num = num + 1 ; oro%name2(num)  = 'ol2'        ! hprime(ix,8)
    num = num + 1 ; oro%name2(num)  = 'ol3'        ! hprime(ix,9)
    num = num + 1 ; oro%name2(num) = 'ol4'        ! hprime(ix,10)
    num = num + 1 ; oro%name2(num) = 'theta'      ! hprime(ix,11)
    num = num + 1 ; oro%name2(num) = 'gamma'      ! hprime(ix,12)
    num = num + 1 ; oro%name2(num) = 'sigma'      ! hprime(ix,13)
    num = num + 1 ; oro%name2(num) = 'elvmax'     ! hprime(ix,14)
    num = num + 1 ; oro%name2(num) = 'orog_filt'  ! oro
    num = num + 1 ; oro%name2(num) = 'orog_raw'   ! oro_uf
    num = num + 1 ; oro%name2(num) = 'land_frac'  ! land fraction [0:1]
    !--- variables below here are optional
    num = num + 1 ; oro%name2(num) = 'lake_frac'  ! lake fraction [0:1]
    num = num + 1 ; oro%name2(num) = 'lake_depth' ! lake depth(m)

    !--- register axis
    call register_axis( Oro_restart, "lon", 'X' )
    call register_axis( Oro_restart, "lat", 'Y' )
    !--- register the 2D fields
    do n = 1,num
      var2_p => oro%var2(:,:,n)
      if (trim(oro%name2(n)) == 'lake_frac' .or. trim(oro%name2(n)) == 'lake_depth' ) then
        call register_restart_field(Oro_restart, oro%name2(n), var2_p, dimensions=(/'lat','lon'/), is_optional=.true.)
      else
        call register_restart_field(Oro_restart, oro%name2(n), var2_p, dimensions=(/'lat','lon'/))
      endif
    enddo

    !--- register 3D vegetation and soil fractions
    var3_fr => oro%var3v(:,:,:)
    call register_restart_field(Oro_restart, 'vegetation_type_pct', var3_fr, dimensions=(/'num_veg_cat','lat        ','lon        '/) , is_optional=.true.)
    var3_fr => oro%var3s(:,:,:)
    call register_restart_field(Oro_restart, 'soil_type_pct', var3_fr, dimensions=(/'num_soil_cat','lat         ','lon         '/) , is_optional=.true.)

  end subroutine Oro_io_register

  !>@brief Copies orography data from temporary arrays back to Sfcprop grid arrays.
  !> \section oro_io_data_type%copy procedure
  !! After reading the restart, data is on temporary arrays with x-y data storage.
  !! This subroutine copies the x-y fields to Sfcprop's blocked grid storage arrays.
  subroutine Oro_io_copy(oro, Model, Sfcprop, Atm_block)
    implicit none
    class(Oro_io_data_type) :: oro
    type(GFS_control_type),      intent(in) :: Model
    type(GFS_sfcprop_type)                  :: Sfcprop
    type(FmsNetcdfDomainFile_t) :: Oro_restart
    type(block_control_type), intent(in) :: Atm_block

    integer :: i,j,nb,ix,num,im

    !$omp parallel do default(shared) private(i, j, nb, ix, num, im)
    do nb = 1, Atm_block%nblks
      !--- 2D variables
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - Atm_block%isc + 1
        j = Atm_block%index(nb)%jj(ix) - Atm_block%jsc + 1
        im = Model%chunk_begin(nb)+ix-1
        !--- stddev
        !       Sfcprop(nb)%hprim(ix)     = oro%var2(i,j,1)
        !--- hprime(1:14)
        num = 1       ; Sfcprop%hprime(im,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop%hprime(im,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop%hprime(im,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop%hprime(im,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop%hprime(im,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop%hprime(im,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop%hprime(im,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop%hprime(im,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop%hprime(im,num)  = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop%hprime(im,num) = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop%hprime(im,num) = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop%hprime(im,num) = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop%hprime(im,num) = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop%hprime(im,num) = oro%var2(i,j,num)
        !--- oro
        num = num + 1 ; Sfcprop%oro(im)       = oro%var2(i,j,num)
        num = num + 1 ; Sfcprop%oro_uf(im)    = oro%var2(i,j,num)

        Sfcprop%landfrac(im)  = -9999.0
        Sfcprop%lakefrac(im)  = -9999.0

        num = num + 1 ; Sfcprop%landfrac(im)  = oro%var2(i,j,num) !land frac [0:1]
        if (Model%lkm > 0  ) then
          if(oro%var2(i,j,num+1)>Model%lakefrac_threshold .and. &
               oro%var2(i,j,num+2)>Model%lakedepth_threshold) then
            Sfcprop%lakefrac(im)  = oro%var2(i,j,num+1) !lake frac [0:1]
            Sfcprop%lakedepth(im) = oro%var2(i,j,num+2) !lake depth [m]    !YWu
          else
            Sfcprop%lakefrac(im)  = 0
            Sfcprop%lakedepth(im) = -9999
          endif
        else
          Sfcprop%lakefrac(im)  = oro%var2(i,j,num+1) !lake frac [0:1]
          Sfcprop%lakedepth(im) = oro%var2(i,j,num+2) !lake depth [m]    !YWu
        endif
        num = num + 2 ! To account for lakefrac and lakedepth

        Sfcprop%vegtype_frac(im,:)  =  -9999.0
        Sfcprop%soiltype_frac(im,:) =  -9999.0

        Sfcprop%vegtype_frac(im,:)  = oro%var3v(i,j,:) ! vegetation type fractions, [0:1]
        Sfcprop%soiltype_frac(im,:) = oro%var3s(i,j,:) ! soil type fractions, [0:1]

      enddo
    enddo

    !--- deallocate containers and free restart container
    deallocate(oro%name2)
    deallocate(oro%var2)
    deallocate(oro%var3v)
    deallocate(oro%var3s)

    nullify(oro%name2)
    nullify(oro%var2)
    nullify(oro%var3v)
    nullify(oro%var3s)

  end subroutine Oro_io_copy

  !>@brief Destructor for Oro_io_data_type
  subroutine Oro_io_final(oro)
    implicit none
    type(Oro_io_data_type) :: oro

    ! This #define reduces code length by a lot
#define IF_ASSOC_DEALLOC_NULL(var) \
    if(associated(oro%var)) then ; \
      deallocate(oro%var) ; \
      nullify(oro%var) ; \
    endif

    IF_ASSOC_DEALLOC_NULL(name2)
    IF_ASSOC_DEALLOC_NULL(var2)
    IF_ASSOC_DEALLOC_NULL(var3s)
    IF_ASSOC_DEALLOC_NULL(var3v)

#undef IF_ASSOC_DEALLOC_NULL
  end subroutine Oro_io_final

  !>@brief Registers axes and fields for non-quilt restart reading of scaled orography variables.
  !> \section Calls FMS restart register functions for axes and
  !! variables in the large-scale or small-scale orography data. The
  !! scaled data is handled by another function. Each scale needs its
  !! own instance of oro_scale_io_data_type.
  subroutine Oro_scale_io_register(oro_scale, Model, Oro_scale_restart, Atm_block)
    implicit none
    class(Oro_scale_io_data_type) :: oro_scale
    type(GFS_control_type),      intent(in) :: Model
    type(FmsNetcdfDomainFile_t) :: Oro_scale_restart
    type(block_control_type), intent(in) :: Atm_block

    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p  => NULL()
    integer :: num, nx, ny

#define WARN_DISASSOCIATE(name) \
    if(associated(name)) then ; \
      write(0,*) 'Internal error. Called oro_scale%register twice. Will try to keep going anyway.' ; \
      deallocate(name); \
      nullify(name) ; \
    endif

    WARN_DISASSOCIATE(oro_scale%name)
    WARN_DISASSOCIATE(oro_scale%var)
#undef WARN_DISASSOCIATE

    call get_nx_ny_from_atm(Atm_block, nx, ny)

    !--- allocate the various containers needed for orography data
    allocate(oro_scale%name(nvar_oro_scale))
    allocate(oro_scale%var(nx,ny,nvar_oro_scale))

    oro_scale%name(1)  = 'stddev'
    oro_scale%name(2)  = 'convexity'
    oro_scale%name(3)  = 'oa1'
    oro_scale%name(4)  = 'oa2'
    oro_scale%name(5)  = 'oa3'
    oro_scale%name(6)  = 'oa4'
    oro_scale%name(7)  = 'ol1'
    oro_scale%name(8)  = 'ol2'
    oro_scale%name(9)  = 'ol3'
    oro_scale%name(10) = 'ol4'

    call register_axis(Oro_scale_restart, "lon", 'X')
    call register_axis(Oro_scale_restart, "lat", 'Y')

    do num = 1,nvar_oro_scale
      var2_p => oro_scale%var(:,:,num)
      call register_restart_field(Oro_scale_restart, oro_scale%name(num), var2_p, dimensions=(/'lon','lat'/))
    enddo
  end subroutine Oro_scale_io_register

  !>@brief Copies scaled orography data from temporary arrays back to Sfcprop grid arrays.
  !> \section Oro_scale_io_data_type%copy procedure
  !! After reading the restart, data is on temporary arrays with x-y data storage.
  !! This subroutine copies the x-y fields to Sfcprop's blocked grid storage arrays.
  subroutine Oro_scale_io_copy(oro_scale, Model, Sfcprop, Atm_block, first_index)
    implicit none
    class(Oro_scale_io_data_type) :: oro_scale
    type(GFS_control_type),   intent(in) :: Model
    type(GFS_sfcprop_type)               :: Sfcprop
    type(block_control_type), intent(in) :: Atm_block
    integer, intent(in) :: first_index

    integer :: i,j,nb,ix,num,v,im

    !$OMP PARALLEL DO PRIVATE(nb,ix,i,j,v,im)
    do nb = 1, Atm_block%nblks
      !--- 2D variables
      do ix = 1, Atm_block%blksz(nb)
        i = Atm_block%index(nb)%ii(ix) - Atm_block%isc + 1
        j = Atm_block%index(nb)%jj(ix) - Atm_block%jsc + 1
        im = Model%chunk_begin(nb)+ix-1
        do v=1,nvar_oro_scale
          Sfcprop%hprime(im,first_index-1+v)  = oro_scale%var(i,j,v)
        enddo
      enddo
    enddo
  end subroutine Oro_scale_io_copy

  !>@brief Oro_scale_io_data_type destructor
  subroutine Oro_scale_io_final(oro_scale)
    implicit none
    type(Oro_scale_io_data_type) :: oro_scale

#define IF_ASSOC_DEALLOC_NULL(vvarr) \
    if(associated(oro_scale%vvarr)) then ; \
      deallocate(oro_scale%vvarr) ; \
      nullify(oro_scale%vvarr) ; \
    endif

    IF_ASSOC_DEALLOC_NULL(name)
    IF_ASSOC_DEALLOC_NULL(var)

#undef IF_ASSOC_DEALLOC_NULL
  end subroutine Oro_scale_io_final
end module fv3atm_oro_io
!> @}
