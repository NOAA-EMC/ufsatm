module ezxml_mod
  use iso_c_binding
  implicit none
  private

  public :: xml_node
  public :: xml_parse_file
  public :: xml_child
  public :: xml_next
  public :: xml_name
  public :: xml_txt
  public :: xml_attr
  public :: xml_free
  public :: xml_is_valid

  !> Opaque handle to an ezXML node. Copy freely; it just wraps a C pointer.
  type :: xml_node
     type(c_ptr) :: ptr = c_null_ptr
  end type xml_node

  interface

     function c_shim_parse_file(file) bind(C, name="shim_parse_file") result(res)
       import :: c_ptr, c_char
       character(kind=c_char), dimension(*), intent(in) :: file
       type(c_ptr) :: res
     end function c_shim_parse_file

     function c_shim_child(xml, name) bind(C, name="shim_child") result(res)
       import :: c_ptr, c_char
       type(c_ptr), value :: xml
       character(kind=c_char), dimension(*), intent(in) :: name
       type(c_ptr) :: res
     end function c_shim_child

     function c_shim_next(xml) bind(C, name="shim_next") result(res)
       import :: c_ptr
       type(c_ptr), value :: xml
       type(c_ptr) :: res
     end function c_shim_next

     function c_shim_name(xml) bind(C, name="shim_name") result(res)
       import :: c_ptr
       type(c_ptr), value :: xml
       type(c_ptr) :: res
     end function c_shim_name

     function c_shim_txt(xml) bind(C, name="shim_txt") result(res)
       import :: c_ptr
       type(c_ptr), value :: xml
       type(c_ptr) :: res
     end function c_shim_txt

     function c_shim_attr(xml, attrname) bind(C, name="shim_attr") result(res)
       import :: c_ptr, c_char
       type(c_ptr), value :: xml
       character(kind=c_char), dimension(*), intent(in) :: attrname
       type(c_ptr) :: res
     end function c_shim_attr

     subroutine c_shim_free(xml) bind(C, name="shim_free")
       import :: c_ptr
       type(c_ptr), value :: xml
     end subroutine c_shim_free

     function c_strlen(s) bind(C, name="strlen") result(res)
       import :: c_ptr, c_size_t
       type(c_ptr), value :: s
       integer(c_size_t) :: res
     end function c_strlen

  end interface

contains

  !> Fortran string -> null-terminated C string (temporary array).
  function f_to_c_string(f_str) result(c_str)
    character(len=*), intent(in) :: f_str
    character(kind=c_char), dimension(len_trim(f_str)+1) :: c_str
    integer :: i, n
    n = len_trim(f_str)
    do i = 1, n
       c_str(i) = f_str(i:i)
    end do
    c_str(n+1) = c_null_char
  end function f_to_c_string

  !> C string pointer -> allocatable Fortran string. Returns "" for NULL.
  function c_to_f_string(cptr) result(f_str)
    type(c_ptr), intent(in) :: cptr
    character(len=:), allocatable :: f_str
    character(kind=c_char), pointer :: fptr(:)
    integer(c_size_t) :: slen
    integer :: i

    if (.not. c_associated(cptr)) then
       f_str = ""
       return
    end if

    slen = c_strlen(cptr)
    if (slen == 0) then
       f_str = ""
       return
    end if

    call c_f_pointer(cptr, fptr, [int(slen)])
    allocate(character(len=int(slen)) :: f_str)
    do i = 1, int(slen)
       f_str(i:i) = fptr(i)
    end do
  end function c_to_f_string

  !> True if a node handle actually points to something.
  function xml_is_valid(node) result(valid)
    type(xml_node), intent(in) :: node
    logical :: valid
    valid = c_associated(node%ptr)
  end function xml_is_valid

  !> Parse an XML file, returning the root element (invalid handle on failure).
  function xml_parse_file(filename) result(node)
    character(len=*), intent(in) :: filename
    type(xml_node) :: node
    node%ptr = c_shim_parse_file(f_to_c_string(filename))
  end function xml_parse_file

  !> First child element with the given tag name (invalid handle if none).
  function xml_child(node, name) result(child)
    type(xml_node), intent(in) :: node
    character(len=*), intent(in) :: name
    type(xml_node) :: child
    child%ptr = c_shim_child(node%ptr, f_to_c_string(name))
  end function xml_child

  !> Next element with the SAME tag name as `node`.
  !> This is what lets you walk repeated elements:
  !>   book = xml_child(root, "book")
  !>   do while (xml_is_valid(book))
  !>      ...
  !>      book = xml_next(book)
  !>   end do
  function xml_next(node) result(nxt)
    type(xml_node), intent(in) :: node
    type(xml_node) :: nxt
    nxt%ptr = c_shim_next(node%ptr)
  end function xml_next

  !> Tag name of this element.
  function xml_name(node) result(name)
    type(xml_node), intent(in) :: node
    character(len=:), allocatable :: name
    name = c_to_f_string(c_shim_name(node%ptr))
  end function xml_name

  !> Text content directly inside this element.
  function xml_txt(node) result(txt)
    type(xml_node), intent(in) :: node
    character(len=:), allocatable :: txt
    txt = c_to_f_string(c_shim_txt(node%ptr))
  end function xml_txt

  !> Value of a named attribute ("" if absent).
  function xml_attr(node, attrname) result(val)
    type(xml_node), intent(in) :: node
    character(len=*), intent(in) :: attrname
    character(len=:), allocatable :: val
    val = c_to_f_string(c_shim_attr(node%ptr, f_to_c_string(attrname)))
  end function xml_attr

  !> Free an entire tree (call once, on the root you got from xml_parse_file).
  subroutine xml_free(node)
    type(xml_node), intent(inout) :: node
    call c_shim_free(node%ptr)
    node%ptr = c_null_ptr
  end subroutine xml_free

end module ezxml_mod
