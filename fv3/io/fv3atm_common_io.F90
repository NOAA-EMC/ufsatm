!>  \file fv3atm_common_io.F90
!! A set of routines commonly accessed by other io/fv3atm
!! modules. This should not be accessed by other code. Most of the
!! routines in this file copy data between x-y-z arrays and
!! block-decomposed (nb-ix-z) atmosphere arrays.

module fv3atm_common_io
  use GFS_typedefs, only: kind_phys

  implicit none
  private

  public :: copy_from_GFS_Data, copy_to_GFS_Data
  public :: copy_from_GFS_Data_2d_phys2phys, copy_from_GFS_Data_3d_phys2phys, &
       copy_from_GFS_Data_2d_int2phys, copy_from_GFS_Data_3d_int2phys, &
       copy_from_GFS_Data_2d_stack_phys2phys, copy_to_GFS_Data_3d_slice_phys2phys, &
       copy_to_GFS_Data_2d_phys2phys, copy_to_GFS_Data_3d_phys2phys, &
       copy_to_GFS_Data_2d_int2phys, copy_to_GFS_Data_3d_int2phys

  public :: GFS_data_transfer
  public :: GFS_data_transfer_2d_phys2phys, &
       GFS_data_transfer_3d_phys2phys, &
       GFS_data_transfer_2d_int2phys, &
       GFS_data_transfer_3d_int2phys, &
       GFS_data_transfer_3d_slice_phys2phys, &
       GFS_data_transfer_2d_stack_phys2phys

  public :: create_2d_field_and_add_to_bundle
  public :: create_3d_field_and_add_to_bundle
  public :: add_zaxis_to_field

  public :: get_nx_ny_from_atm

#ifdef CCPP_32BIT
  character(len=5), parameter, public :: axis_type = 'float'
#else
  character(len=6), parameter, public :: axis_type = 'double'
#endif

  !>\defgroup fv3atm_common_io FV3ATM Common I/O Utilities Module
  !> @{

  !>@ These subroutines copy data from x-y-z arrays to nb-ix-z grid arrays.
  !! \section copy_from_GFS_Data interface
  !!  There are different combinations of decomposition, copy methods,
  !!  and datatypes. All are combined together into copy_from_GFS_Data
  !!  for convenience
  interface copy_from_GFS_Data
    module procedure copy_from_GFS_Data_2d_phys2phys, &
         copy_from_GFS_Data_3d_phys2phys, &
         copy_from_GFS_Data_2d_int2phys, &
         copy_from_GFS_Data_3d_int2phys, &
         copy_from_GFS_Data_3d_slice_phys2phys, &
         copy_from_GFS_Data_2d_stack_phys2phys
  end interface copy_from_GFS_Data

  !>@ These subroutines copy data from nb-ix-z grid arrays to x-y-z arrays.
  !!  \section copy_to_GFS_Data interface
  !!  There are different combinations of decomposition, copy methods,
  !!  and datatypes. All are combined together into copy_to_GFS_Data
  !!  for convenience
  interface copy_to_GFS_Data
    module procedure copy_to_GFS_Data_2d_phys2phys, &
         copy_to_GFS_Data_3d_phys2phys, &
         copy_to_GFS_Data_2d_int2phys, &
         copy_to_GFS_Data_3d_int2phys, &
         copy_to_GFS_Data_3d_slice_phys2phys, &
         copy_to_GFS_Data_2d_stack_phys2phys
  end interface copy_to_GFS_Data

  !>@brief These subroutines copy data in either direction between nb-ix-z grid arrays and x-y-z arrays.
  !> \section GFS_data_transfer interface functions.
  !! This interface allows a single subroutine to handle both reading
  !! and writing restart files. The direction is controled by the "to"
  !! argument (first argument) which is true when copying from x-y-z
  !! arrays to nb-ix-z arrays.
  !! There are different combinations of decomposition, copy methods,
  !! and datatypes. All are combined together into copy_to_GFS_Data
  !! for convenience
  interface GFS_data_transfer
    module procedure GFS_data_transfer_2d_phys2phys, &
         GFS_data_transfer_3d_phys2phys, &
         GFS_data_transfer_2d_int2phys, &
         GFS_data_transfer_3d_int2phys, &
         GFS_data_transfer_3d_slice_phys2phys, &
         GFS_data_transfer_2d_stack_phys2phys
  end interface GFS_data_transfer

contains

  !>@brief Convenience function to get the x and y dimensions of the grid from Atm_block
  pure subroutine get_nx_ny_from_atm(Atm_block, nx, ny)
    use block_control_mod,  only: block_control_type
    implicit none
    type(block_control_type), intent(in) :: Atm_block
    integer, intent(out), optional :: nx, ny
    integer :: isc, iec, jsc, jec
    if(present(nx)) then
      isc = Atm_block%isc
      iec = Atm_block%iec
      nx  = (iec - isc + 1)
    end if
    if(present(ny)) then
      jsc = Atm_block%jsc
      jec = Atm_block%jec
      ny  = (jec - jsc + 1)
    endif
  end subroutine get_nx_ny_from_atm

  !>@brief copies from the ix-indexed var_block to the 2d x-y real(kind_phys) var2d array
  pure subroutine copy_from_GFS_Data_2d_phys2phys(ii1,jj1,isc,jsc,nt,var2d,var_block)
    implicit none
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc
    integer, intent(inout) :: nt
    real(kind=kind_phys), intent(in) :: var_block(:)
    real(kind=kind_phys), intent(out) :: var2d(:,:,:)
    integer ix

    nt=nt+1
    do ix=1,size(var_block)
      var2d(ii1(ix),jj1(ix),nt) = var_block(ix)
    enddo
  end subroutine copy_from_GFS_Data_2d_phys2phys

  !>@brief copies from the ix-k-indexed var_block to the 3d x-y-z real(kind_phys) var3d array
  pure subroutine copy_from_GFS_Data_3d_phys2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
    implicit none
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc
    integer, intent(inout) :: nt
    real(kind=kind_phys), intent(in) :: var_block(:,:)
    real(kind=kind_phys), intent(out) :: var3d(:,:,:,:)
    integer ix, k

    nt=nt+1
    do k=lbound(var_block,2),ubound(var_block,2)
      do ix=1,size(var_block,1)
        var3d(ii1(ix),jj1(ix),k,nt) = var_block(ix,k)
      enddo
    enddo
  end subroutine copy_from_GFS_Data_3d_phys2phys

  !>@brief copies from the ix-k-indexed var_block to the 3d x-y-z integer var2d array
  pure subroutine copy_from_GFS_Data_2d_int2phys(ii1,jj1,isc,jsc,nt,var2d,var_block)
    implicit none
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc, var_block(:)
    integer, intent(inout) :: nt
    real(kind=kind_phys), intent(out) :: var2d(:,:,:)
    integer ix

    nt=nt+1
    do ix=1,size(var_block)
      var2d(ii1(ix),jj1(ix),nt) = var_block(ix)
    enddo
  end subroutine copy_from_GFS_Data_2d_int2phys

  !>@brief copies a range of levels from the ix-k-indexed var_block to the x-y real(kind_phys) var3d array
  pure subroutine copy_from_GFS_Data_2d_stack_phys2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
    ! For copying phy_f2d and phy_fctd
    implicit none
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc
    integer, intent(inout) :: nt
    real(kind=kind_phys), intent(in) :: var_block(:,:)
    real(kind=kind_phys), intent(out) :: var3d(:,:,:)
    integer ix, k

    nt=nt+1
    do k=lbound(var_block,2),ubound(var_block,2)
      do ix=1,size(var_block,1)
        var3d(ii1(ix),jj1(ix),nt) = var_block(ix,k)
      enddo
    enddo
  end subroutine copy_from_GFS_Data_2d_stack_phys2phys

  !>@brief copies from the ix-k-indexed var_block to the x-y integer var3d array
  pure subroutine copy_from_GFS_Data_3d_int2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
    implicit none
    integer, intent(in) :: ii1(:), jj1(:), var_block(:,:), isc, jsc
    integer, intent(inout) :: nt
    real(kind=kind_phys), intent(out) :: var3d(:,:,:,:)
    integer ix, k

    nt=nt+1
    do k=lbound(var_block,2),ubound(var_block,2)
      do ix=1,size(var_block,1)
        var3d(ii1(ix),jj1(ix),k,nt) = real(var_block(ix,k),kind_phys)
      enddo
    enddo
  end subroutine copy_from_GFS_Data_3d_int2phys

  !>@brief copies a range of levels from from the ix-k-indexed var_block to the x-y-z real(kind_phys) var3d array
  pure subroutine copy_from_GFS_Data_3d_slice_phys2phys(ii1,jj1,isc,jsc,nt,k1,k2,var3d,var_block)
    implicit none
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc, k1, k2
    integer, intent(inout) :: nt
    real(kind=kind_phys), intent(in) :: var_block(:,:)
    real(kind=kind_phys), intent(out) :: var3d(:,:,:,:)
    integer ix, k

    nt=nt+1
    do k=k1,k2
      do ix=1,size(var_block,1)
        var3d(ii1(ix),jj1(ix),k,nt) = var_block(ix,k)
      enddo
    enddo
  end subroutine copy_from_GFS_Data_3d_slice_phys2phys

  !>@brief copies from x-y real(kind_phys) var2d array to the ix-indexed var_block array
  pure subroutine copy_to_GFS_Data_2d_phys2phys(ii1,jj1,isc,jsc,nt,var2d,var_block)
    implicit none
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc
    integer, intent(inout) :: nt
    real(kind=kind_phys), intent(out) :: var_block(:)
    real(kind=kind_phys), intent(in) :: var2d(:,:,:)
    integer ix

    nt=nt+1
    do ix=1,size(var_block)
      var_block(ix) = var2d(ii1(ix),jj1(ix),nt)
    enddo
  end subroutine copy_to_GFS_Data_2d_phys2phys

  !>@brief copies from x-y-z real(kind_phys) var3d array to the ix-k-indexed var_block array
  pure subroutine copy_to_GFS_Data_3d_phys2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
    implicit none
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc
    integer, intent(inout) :: nt
    real(kind=kind_phys), intent(out) :: var_block(:,:)
    real(kind=kind_phys), intent(in) :: var3d(:,:,:,:)
    integer ix, k

    nt=nt+1
    do k=lbound(var_block,2),ubound(var_block,2)
      do ix=1,size(var_block,1)
        var_block(ix,k) = var3d(ii1(ix),jj1(ix),k,nt)
      enddo
    enddo
  end subroutine copy_to_GFS_Data_3d_phys2phys

  !>@brief copies from x-y-z real(kind_phys) var3d array to the ix-k-indexed var_block array
  pure subroutine copy_to_GFS_Data_2d_stack_phys2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
    ! For copying phy_f2d and phy_fctd
    implicit none
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc
    integer, intent(inout) :: nt
    real(kind=kind_phys), intent(out) :: var_block(:,:)
    real(kind=kind_phys), intent(in) :: var3d(:,:,:)
    integer ix, k

    nt=nt+1
    do k=lbound(var_block,2),ubound(var_block,2)
      do ix=1,size(var_block,1)
        var_block(ix,k) = var3d(ii1(ix),jj1(ix),nt)
      enddo
    enddo
  end subroutine copy_to_GFS_Data_2d_stack_phys2phys

  !>@brief copies a range of levels from the x-y-z real(kind_phys) var3d array to the ix-k-indexed var_block array
  pure subroutine copy_to_GFS_Data_3d_slice_phys2phys(ii1,jj1,isc,jsc,nt,k1,k2,var3d,var_block)
    implicit none
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc, k1, k2
    integer, intent(inout) :: nt
    real(kind=kind_phys), intent(out) :: var_block(:,:)
    real(kind=kind_phys), intent(in) :: var3d(:,:,:,:)
    integer ix, k

    nt=nt+1
    do k=k1,k2
      do ix=1,size(var_block,1)
        var_block(ix,k) = var3d(ii1(ix),jj1(ix),k,nt)
      enddo
    enddo
  end subroutine copy_to_GFS_Data_3d_slice_phys2phys

  !>@brief copies from x-y integer var2d array to the ix-indexed var_block array
  pure subroutine copy_to_GFS_Data_2d_int2phys(ii1,jj1,isc,jsc,nt,var2d,var_block)
    implicit none
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc
    integer, intent(inout) :: nt
    integer, intent(out) :: var_block(:)
    real(kind=kind_phys), intent(in) :: var2d(:,:,:)
    integer ix

    nt=nt+1
    do ix=1,size(var_block)
      var_block(ix) = int(var2d(ii1(ix),jj1(ix),nt))
    enddo
  end subroutine copy_to_GFS_Data_2d_int2phys

  !>@brief copies from x-y-z integer var3d array to the ix-k-indexed var_block array
  pure subroutine copy_to_GFS_Data_3d_int2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
    implicit none
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc
    integer, intent(inout) :: nt
    integer, intent(out) :: var_block(:,:)
    real(kind=kind_phys), intent(in) :: var3d(:,:,:,:)
    integer ix

    nt=nt+1
    do ix=1,size(var_block,1)
      var_block(ix,:) = int(var3d(ii1(ix),jj1(ix),:,nt))
    enddo
  end subroutine copy_to_GFS_Data_3d_int2phys

  !>@brief copies between the ix-indexed var_block array and x-y real(kind_phys) var2d array.
  !> \section GFS_data_transfer_2d_phys2phys subroutine from the GFS_data_transfer interface
  !! This is a wrapper around copy_to_GFS_Data and copy_from_GFS_Data routines.
  !! If to=true, then data is copied to var_block (the GFS_Data structures) but if
  !! to=false, it is copied from the var_block arrays. This allows the same subroutine
  !! to both read and write, preventing error-prone code duplication.
  pure subroutine GFS_data_transfer_2d_phys2phys(to,ii1,jj1,isc,jsc,nt,var2d,var_block)
    implicit none
    logical, intent(in) :: to
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc
    integer, intent(inout) :: nt
    real(kind=kind_phys), intent(inout) :: var_block(:)
    real(kind=kind_phys), intent(inout) :: var2d(:,:,:)

    if(to) then
      call copy_to_GFS_Data_2d_phys2phys(ii1,jj1,isc,jsc,nt,var2d,var_block)
    else
      call copy_from_GFS_Data_2d_phys2phys(ii1,jj1,isc,jsc,nt,var2d,var_block)
    end if
  end subroutine GFS_data_transfer_2d_phys2phys

  !>@brief copies between the ix-k-indexed var_block array and x-y-z real(kind_phys) var3d array.
  !> \section GFS_data_transfer_3d_phys2phys subroutine from the GFS_data_transfer interface
  !! This is a wrapper around copy_to_GFS_Data and copy_from_GFS_Data routines.
  !! If to=true, then data is copied to var_block (the GFS_Data structures) but if
  !! to=false, it is copied from the var_block arrays. This allows the same subroutine
  !! to both read and write, preventing error-prone code duplication.
  pure subroutine GFS_data_transfer_3d_phys2phys(to,ii1,jj1,isc,jsc,nt,var3d,var_block)
    implicit none
    logical, intent(in) :: to
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc
    integer, intent(inout) :: nt
    real(kind=kind_phys), intent(inout) :: var_block(:,:)
    real(kind=kind_phys), intent(inout) :: var3d(:,:,:,:)

    if(to) then
      call copy_to_GFS_Data_3d_phys2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
    else
      call copy_from_GFS_Data_3d_phys2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
    endif
  end subroutine GFS_data_transfer_3d_phys2phys

  !>@brief copies a range of levels between the ix-k-indexed var_block array and x-y-z real(kind_phys) var3d array.
  !> \section GFS_data_transfer_3d_slice_phys2phys subroutine from the GFS_data_transfer interface
  !! This is a wrapper around copy_to_GFS_Data and copy_from_GFS_Data routines.
  !! If to=true, then data is copied to var_block (the GFS_Data structures) but if
  !! to=false, it is copied from the var_block arrays. This allows the same subroutine
  !! to both read and write, preventing error-prone code duplication.
  pure subroutine GFS_data_transfer_3d_slice_phys2phys(to,ii1,jj1,isc,jsc,nt,k1,k2,var3d,var_block)
    implicit none
    logical, intent(in) :: to
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc, k1, k2
    integer, intent(inout) :: nt
    real(kind=kind_phys), intent(inout) :: var_block(:,:)
    real(kind=kind_phys), intent(inout) :: var3d(:,:,:,:)

    if(to) then
      call copy_to_GFS_Data_3d_slice_phys2phys(ii1,jj1,isc,jsc,nt,k1,k2,var3d,var_block)
    else
      call copy_from_GFS_Data_3d_slice_phys2phys(ii1,jj1,isc,jsc,nt,k1,k2,var3d,var_block)
    endif
  end subroutine GFS_data_transfer_3d_slice_phys2phys

  !>@brief copies between the ix-indexed var_block array and x-y integer var2d array.
  !> \section GFS_data_transfer_2d_int2phys subroutine from the GFS_data_transfer interface
  !! This is a wrapper around copy_to_GFS_Data and copy_from_GFS_Data routines.
  !! If to=true, then data is copied to var_block (the GFS_Data structures) but if
  !! to=false, it is copied from the var_block arrays. This allows the same subroutine
  !! to both read and write, preventing error-prone code duplication.
  pure subroutine GFS_data_transfer_2d_int2phys(to,ii1,jj1,isc,jsc,nt,var2d,var_block)
    implicit none
    logical, intent(in) :: to
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc
    integer, intent(inout) :: nt
    integer, intent(inout) :: var_block(:)
    real(kind=kind_phys), intent(inout) :: var2d(:,:,:)

    if(to) then
      call copy_to_GFS_Data_2d_int2phys(ii1,jj1,isc,jsc,nt,var2d,var_block)
    else
      call copy_from_GFS_Data_2d_int2phys(ii1,jj1,isc,jsc,nt,var2d,var_block)
    endif
  end subroutine GFS_data_transfer_2d_int2phys

  !>@brief copies between the ix-k-indexed var_block array and x-y-z integer var3d array.
  !> \section GFS_data_transfer_3d_int2phys subroutine from the GFS_data_transfer interface
  !! This is a wrapper around copy_to_GFS_Data and copy_from_GFS_Data routines.
  !! If to=true, then data is copied to var_block (the GFS_Data structures) but if
  !! to=false, it is copied from the var_block arrays. This allows the same subroutine
  !! to both read and write, preventing error-prone code duplication.
  pure subroutine GFS_data_transfer_3d_int2phys(to,ii1,jj1,isc,jsc,nt,var3d,var_block)
    implicit none
    logical, intent(in) :: to
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc
    integer, intent(inout) :: nt
    integer, intent(inout) :: var_block(:,:)
    real(kind=kind_phys), intent(inout) :: var3d(:,:,:,:)

    if(to) then
      call copy_to_GFS_Data_3d_int2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
    else
      call copy_from_GFS_Data_3d_int2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
    endif
  end subroutine GFS_data_transfer_3d_int2phys

  !>@brief copies a range of levels between the ix-k-indexed var_block array and x-y-z real(kind_phys) var3d array.
  !> \section GFS_Data_transfer_2d_stack_phys2phys subroutine from the GFS_data_transfer interface
  !! This is a wrapper around copy_to_GFS_Data and copy_from_GFS_Data routines.
  !! If to=true, then data is copied to var_block (the GFS_Data structures) but if
  !! to=false, it is copied from the var_block arrays. This allows the same subroutine
  !! to both read and write, preventing error-prone code duplication.
  pure subroutine GFS_Data_transfer_2d_stack_phys2phys(to,ii1,jj1,isc,jsc,nt,var3d,var_block)
    ! For copying phy_f2d and phy_fctd
    implicit none
    logical, intent(in) :: to
    integer, intent(in) :: ii1(:), jj1(:), isc, jsc
    integer, intent(inout) :: nt
    real(kind=kind_phys), intent(inout) :: var_block(:,:)
    real(kind=kind_phys), intent(inout) :: var3d(:,:,:)
    integer ix, k

    if(to) then
      call copy_to_GFS_data_2d_stack_phys2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
    else
      call copy_from_GFS_data_2d_stack_phys2phys(ii1,jj1,isc,jsc,nt,var3d,var_block)
    end if
  end subroutine GFS_Data_transfer_2d_stack_phys2phys

  !>@brief adds a 2D restart array to an ESMF bundle for quilting restarts.
  subroutine create_2d_field_and_add_to_bundle(temp_r2d, field_name, outputfile, grid, bundle)

    use esmf

    implicit none

    real(kind_phys), dimension(:,:),   pointer, intent(in)    :: temp_r2d
    character(len=*),                           intent(in)    :: field_name
    character(len=*),                           intent(in)    :: outputfile
    type(ESMF_Grid),                            intent(in)    :: grid
    type(ESMF_FieldBundle),                     intent(inout) :: bundle

    type(ESMF_Field) :: field

    integer :: rc, i

    field = ESMF_FieldCreate(grid, temp_r2d, datacopyflag=ESMF_DATACOPY_REFERENCE, &
         name=trim(field_name), indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, file=__FILE__)) &
         call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", name='output_file', value=trim(outputfile), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_FieldBundleAdd(bundle, (/field/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  end subroutine create_2d_field_and_add_to_bundle

  !>@brief adds a 3D restart array and its vertical axis to an ESMF bundle for quilting restarts.
  subroutine create_3d_field_and_add_to_bundle(temp_r3d, field_name, axis_name, axis_values, outputfile, grid, bundle)

    use esmf

    implicit none

    real(kind_phys), dimension(:,:,:), pointer, intent(in)    :: temp_r3d
    character(len=*),                           intent(in)    :: field_name
    character(len=*),                           intent(in)    :: axis_name
    real(kind_phys), dimension(:),              intent(in)    :: axis_values
    character(len=*),                           intent(in)    :: outputfile
    type(ESMF_Grid),                            intent(in)    :: grid
    type(ESMF_FieldBundle),                     intent(inout) :: bundle

    type(ESMF_Field) :: field

    integer :: rc, i

    field = ESMF_FieldCreate(grid, temp_r3d, datacopyflag=ESMF_DATACOPY_REFERENCE, &
         name=trim(field_name), indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, file=__FILE__)) &
         call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", name='output_file', value=trim(outputfile), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call add_zaxis_to_field(field, axis_name, axis_values)

    call ESMF_FieldBundleAdd(bundle, (/field/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  end subroutine create_3d_field_and_add_to_bundle

  !>@brief adds a vertical axis to an ESMF bundle for quilting restarts.
  subroutine add_zaxis_to_field(field, axis_name, axis_values)

    use esmf

    implicit none

    type(ESMF_Field), intent(inout) :: field
    character(len=*), intent(in)    :: axis_name
    real(kind_phys), dimension(:), intent(in) :: axis_values

    integer :: rc

    call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3", &
         name="ESMF:ungridded_dim_labels", valueList=(/trim(axis_name)/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3-dim", &
         name=trim(axis_name), valueList=axis_values, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeSet(field, convention="NetCDF", purpose="FV3-dim", &
         name=trim(axis_name)//"cartesian_axis", value="Z", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine add_zaxis_to_field

end module fv3atm_common_io
!> @}
