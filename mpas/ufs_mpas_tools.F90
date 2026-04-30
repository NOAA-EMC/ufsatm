!> ###########################################################################################
!> \file ufs_mpas_tools.F90
!>
!>
!> ###########################################################################################
module ufs_mpas_tools
  implicit none

  public
contains
  !> #########################################################################################
  !> Convert one or more values of any intrinsic data types to a character string for pretty
  !> printing.
  !>
  !> If `value` contains more than one element, the elements will be stringified, delimited by
  !> `separator`, then concatenated.
  !> If `value` contains exactly one element, the element will be stringified without using
  !> `separator`.
  !> If `value` contains zero element or is of unsupported data types, an empty character
  !> string is produced.
  !> If `separator` is not supplied, it defaults to ", " (i.e., a comma and a space).
  !> (KCW, 2024-02-04)
  !>
  !> \update: Dustin Swales April 2025 - Modified for use in UWM
  !>
  !> #########################################################################################
  pure function stringify(value, separator)
    use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64

    class(*), intent(in) :: value(:)
    character(*), optional, intent(in) :: separator
    character(:), allocatable :: stringify

    integer, parameter :: sizelimit = 1024

    character(:), allocatable :: buffer, delimiter, format
    character(:), allocatable :: value_c(:)
    integer :: i, n, offset

    if (present(separator)) then
       delimiter = separator
    else
       delimiter = ', '
    end if

    n = min(size(value), sizelimit)

    if (n == 0) then
       stringify = ''

       return
    end if

    select type (value)
    type is (character(*))
       allocate(character(len(value) * n + len(delimiter) * (n - 1)) :: buffer)

       buffer(:) = ''
       offset = 0

       ! Workaround for a bug in GNU Fortran >= 12. This is perhaps the manifestation of GCC Bugzilla Bug 100819.
       ! When a character string array is passed as the actual argument to an unlimited polymorphic dummy argument,
       ! its array index and length parameter are mishandled.
       allocate(character(len(value)) :: value_c(size(value)))

       value_c(:) = value(:)

       do i = 1, n
          if (len(delimiter) > 0 .and. i > 1) then
             buffer(offset + 1:offset + len(delimiter)) = delimiter
             offset = offset + len(delimiter)
          end if

          if (len_trim(adjustl(value_c(i))) > 0) then
             buffer(offset + 1:offset + len_trim(adjustl(value_c(i)))) = trim(adjustl(value_c(i)))
             offset = offset + len_trim(adjustl(value_c(i)))
          end if
       end do

       deallocate(value_c)
    type is (integer(int32))
       allocate(character(11 * n + len(delimiter) * (n - 1)) :: buffer)
       allocate(character(17 + len(delimiter) + floor(log10(real(n))) + 1) :: format)

       write(format, '(a, i0, 3a)') '(ss, ', n, '(i0, :, "', delimiter, '"))'
       write(buffer, format) value
    type is (integer(int64))
       allocate(character(20 * n + len(delimiter) * (n - 1)) :: buffer)
       allocate(character(17 + len(delimiter) + floor(log10(real(n))) + 1) :: format)

       write(format, '(a, i0, 3a)') '(ss, ', n, '(i0, :, "', delimiter, '"))'
       write(buffer, format) value
    type is (logical)
       allocate(character(1 * n + len(delimiter) * (n - 1)) :: buffer)
       allocate(character(13 + len(delimiter) + floor(log10(real(n))) + 1) :: format)

       write(format, '(a, i0, 3a)') '(', n, '(l1, :, "', delimiter, '"))'
       write(buffer, format) value
    type is (real(real32))
       allocate(character(13 * n + len(delimiter) * (n - 1)) :: buffer)

       if (maxval(abs(value)) < 1.0e5_real32) then
          allocate(character(20 + len(delimiter) + floor(log10(real(n))) + 1) :: format)
          write(format, '(a, i0, 3a)') '(ss, ', n, '(f13.6, :, "', delimiter, '"))'
       else
          allocate(character(23 + len(delimiter) + floor(log10(real(n))) + 1) :: format)
          write(format, '(a, i0, 3a)') '(ss, ', n, '(es13.6e2, :, "', delimiter, '"))'
       end if

       write(buffer, format) value
    type is (real(real64))
       allocate(character(13 * n + len(delimiter) * (n - 1)) :: buffer)

       if (maxval(abs(value)) < 1.0e5_real64) then
          allocate(character(20 + len(delimiter) + floor(log10(real(n))) + 1) :: format)
          write(format, '(a, i0, 3a)') '(ss, ', n, '(f13.6, :, "', delimiter, '"))'
       else
          allocate(character(23 + len(delimiter) + floor(log10(real(n))) + 1) :: format)
          write(format, '(a, i0, 3a)') '(ss, ', n, '(es13.6e2, :, "', delimiter, '"))'
       end if

       write(buffer, format) value
    class default
       stringify = ''

       return
    end select

    stringify = trim(buffer)
  end function stringify

  ! ##########################################################################################
  ! Convert <YYYYMMDD> into <YYYY>-<MM>-<DD>
  ! ##########################################################################################
  character(len=10) function date2yyyymmdd (date)
    integer, intent(in) :: date    ! yyyymmdd
    integer             :: year    ! year of yyyy-mm-dd
    integer             :: month   ! month of yyyy-mm-dd
    integer             :: day     ! day of yyyy-mm-dd

    year  = date / 10000
    month = (date - year*10000) / 100
    day   = date - year*10000 - month*100

    write(date2yyyymmdd,80) year, month, day
80  format(i4.4,'-',i2.2,'-',i2.2)

  end function date2yyyymmdd

  ! #########################################################################################
  ! Convert <seconds> into <hours>:<minutes>:<seconds>
  ! #########################################################################################
  character(len=8) function sec2hms (seconds)
    integer, intent(in) :: seconds   ! seconds
    integer             :: hours     ! hours of hh:mm:ss
    integer             :: minutes   ! minutes of hh:mm:ss
    integer             :: secs      ! seconds of hh:mm:ss

    hours   = seconds / 3600
    minutes = (seconds - hours*3600) / 60
    secs    = (seconds - hours*3600 - minutes*60)

    write(sec2hms,80) hours, minutes, secs
80  format(i2.2,':',i2.2,':',i2.2)

  end function sec2hms

  ! #########################################################################################
  ! Convert <integer> into a left justified string.
  ! #########################################################################################
  character(len=10) function int2str(n)
    integer, intent(in) :: n
    write(int2str,'(i0)') n
  end function int2str

  !> #########################################################################################
  !> Convert <logical> as a left justified string.
  !> ######################################################################################### 
  character(len=10) function log2str(n)
    logical, intent(in) :: n
    if (n) then
       write(log2str,'(a4)') 'TRUE'
    else
       write(log2str,'(a4)') 'FALSE'
    endif
  end function log2str

end module ufs_mpas_tools
