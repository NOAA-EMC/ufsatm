module module_block_harness

  use block_control_mod, only: block_control_type
  use CCPP_data,         only: GFS_control
  use GFS_typedefs,      only: GFS_kind_phys => kind_phys

  implicit none

  private

  type, public :: test_config_type
     integer :: nx
     integer :: ny
     integer :: inpes
     integer :: jnpes
     integer :: blocksize
     character(len=64) :: description
  end type test_config_type

  public :: setup_test_configurations
  public :: initialize_block_control
  public :: cleanup_block_control
  public :: setup_copy2block_state
  public :: packed_index

contains

  subroutine setup_test_configurations(configs)
    type(test_config_type), intent(out) :: configs(:)

    configs(1)%nx = 8
    configs(1)%ny = 8
    configs(1)%inpes = 2
    configs(1)%jnpes = 4
    configs(1)%blocksize = 8
    configs(1)%description = "8x8, inpes=2, jnpes=4, bs=8"

    configs(2)%nx = 8
    configs(2)%ny = 8
    configs(2)%inpes = 2
    configs(2)%jnpes = 2
    configs(2)%blocksize = 16
    configs(2)%description = "8x8, inpes=2, jnpes=2, bs=16"

    configs(3)%nx = 8
    configs(3)%ny = 8
    configs(3)%inpes = 1
    configs(3)%jnpes = 8
    configs(3)%blocksize = 8
    configs(3)%description = "8x8, inpes=1, jnpes=8, bs=8 (linear in Y)"

    configs(4)%nx = 8
    configs(4)%ny = 8
    configs(4)%inpes = 4
    configs(4)%jnpes = 2
    configs(4)%blocksize = 8
    configs(4)%description = "8x8, inpes=4, jnpes=2, bs=8 (uneven both)"

    configs(5)%nx = 16
    configs(5)%ny = 16
    configs(5)%inpes = 1
    configs(5)%jnpes = 4
    configs(5)%blocksize = 16
    configs(5)%description = "16x16, inpes=1, jnpes=4, bs=16 (linear in Y)"

    configs(6)%nx = 16
    configs(6)%ny = 16
    configs(6)%inpes = 3
    configs(6)%jnpes = 3
    configs(6)%blocksize = 16
    configs(6)%description = "16x16, inpes=3, jnpes=3, bs=16 (uneven both)"

    configs(7)%nx = 8
    configs(7)%ny = 8
    configs(7)%inpes = 2
    configs(7)%jnpes = 6
    configs(7)%blocksize = 8
    configs(7)%description = "8x8, inpes=2, jnpes=6, bs=8 (uneven Y division)"

    configs(8)%nx = 16
    configs(8)%ny = 16
    configs(8)%inpes = 2
    configs(8)%jnpes = 8
    configs(8)%blocksize = 16
    configs(8)%description = "16x16, inpes=2, jnpes=8, bs=16 (many blocks)"
  end subroutine setup_test_configurations

  subroutine initialize_block_control(block_ctl, nx, ny, inpes, jnpes, blocksize)
    type(block_control_type), intent(out) :: block_ctl
    integer, intent(in) :: nx, ny, inpes, jnpes, blocksize

    integer :: nblocks, iblock, jblock, block_id
    integer :: i, j, istart, jstart, iend, jend
    integer :: npts, ipt

    nblocks = inpes * jnpes

    allocate(block_ctl%blksz(nblocks))
    allocate(block_ctl%index(nblocks))

    block_ctl%isc = 1
    block_ctl%iec = nx
    block_ctl%jsc = 1
    block_ctl%jec = ny
    block_ctl%nblks = nblocks

    allocate(block_ctl%blkno(block_ctl%isc:block_ctl%iec, block_ctl%jsc:block_ctl%jec))
    allocate(block_ctl%ixp(block_ctl%isc:block_ctl%iec, block_ctl%jsc:block_ctl%jec))
    block_ctl%blkno = 0
    block_ctl%ixp = 0

    block_id = 1
    do jblock = 1, jnpes
       do iblock = 1, inpes
          istart = ((iblock - 1) * nx) / inpes + 1
          iend = (iblock * nx) / inpes
          jstart = ((jblock - 1) * ny) / jnpes + 1
          jend = (jblock * ny) / jnpes

          npts = (iend - istart + 1) * (jend - jstart + 1)
          block_ctl%blksz(block_id) = npts

          allocate(block_ctl%index(block_id)%ii(npts))
          allocate(block_ctl%index(block_id)%jj(npts))

          ipt = 1
          do j = jstart, jend
             do i = istart, iend
                block_ctl%index(block_id)%ii(ipt) = i
                block_ctl%index(block_id)%jj(ipt) = j
                block_ctl%blkno(i, j) = block_id
                block_ctl%ixp(i, j) = ipt
                ipt = ipt + 1
             end do
          end do
          block_id = block_id + 1
       end do
    end do

    print *, "Block Control Initialized (MPI task local domain):"
    print *, "  Local domain: isc=", block_ctl%isc, " iec=", block_ctl%iec, " jsc=", block_ctl%jsc, " jec=", block_ctl%jec
    print *, "  Number of blocks: ", nblocks
    print *, "  Block sizes range from: 1 to ", maxval(block_ctl%blksz)
    if (blocksize <= 0) print *, "  WARNING: non-positive blocksize requested"
  end subroutine initialize_block_control

  subroutine cleanup_block_control(block_ctl)
    type(block_control_type), intent(inout) :: block_ctl
    integer :: i

    if (allocated(block_ctl%blksz)) deallocate(block_ctl%blksz)
    if (allocated(block_ctl%blkno)) deallocate(block_ctl%blkno)
    if (allocated(block_ctl%ixp)) deallocate(block_ctl%ixp)
    if (allocated(block_ctl%index)) then
       do i = 1, size(block_ctl%index)
          if (allocated(block_ctl%index(i)%ii)) deallocate(block_ctl%index(i)%ii)
          if (allocated(block_ctl%index(i)%jj)) deallocate(block_ctl%index(i)%jj)
       end do
       deallocate(block_ctl%index)
    end if
    if (associated(GFS_control%chunk_begin)) deallocate(GFS_control%chunk_begin)
  end subroutine cleanup_block_control

  subroutine setup_copy2block_state(block_ctl)
    type(block_control_type), intent(in) :: block_ctl

    integer :: block_id, offset

    GFS_control%isc = block_ctl%isc
    GFS_control%jsc = block_ctl%jsc
    GFS_control%nx = block_ctl%iec - block_ctl%isc + 1
    GFS_control%ny = block_ctl%jec - block_ctl%jsc + 1
    GFS_control%huge = 9.9692099683868690E30_GFS_kind_phys

    if (associated(GFS_control%chunk_begin)) deallocate(GFS_control%chunk_begin)
    allocate(GFS_control%chunk_begin(block_ctl%nblks))

    offset = 1
    do block_id = 1, block_ctl%nblks
       GFS_control%chunk_begin(block_id) = offset
       offset = offset + block_ctl%blksz(block_id)
    end do
  end subroutine setup_copy2block_state

  integer function packed_index(block_ctl, i, j)
    type(block_control_type), intent(in) :: block_ctl
    integer, intent(in) :: i, j

    packed_index = GFS_control%chunk_begin(block_ctl%blkno(i, j)) + block_ctl%ixp(i, j) - 1
  end function packed_index

end module module_block_harness
