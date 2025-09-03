program test_my_new_procedure

    implicit none

    character(20) :: my_input

    my_input = "teststring"

    call my_new_procedure(my_input)

    print*, "Got output: ", trim(my_input)

    if (my_input .ne. "teststring_suffix") then
        stop 1
    endif

end program test_my_new_procedure
