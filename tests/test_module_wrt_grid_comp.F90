program test_module_wrt_grid_comp
    use module_wrt_grid_comp, only: get_outfile, lambert, rtll
    implicit none
  
    call test_get_outfile()
    call test_lambert()
    call test_rtll()

contains
    
    !---------------------------------------------------------------------------
    ! Test get_outfile subroutine
    !---------------------------------------------------------------------------
    subroutine test_get_outfile()
        character(len=128) :: filename(2000,3)
        character(len=128) :: outfile_name(100)
        integer :: noutfile
        character(len=100) :: test_name
        
        ! Test 1: Basic test with unique filenames
        filename = ''
        filename(1,1) = 'file1.nc'
        filename(2,1) = 'file2.nc'
        filename(1,2) = 'file3.nc'
        filename(1,3) = 'file4.nc'
        
        call get_outfile(3, filename, outfile_name, noutfile)
        if ( trim(outfile_name(1)) /= "file1.nc" .or. &
             trim(outfile_name(2)) /= "file2.nc" .or. &
             trim(outfile_name(3)) /= "file3.nc" .or. &
             trim(outfile_name(4)) /= "file4.nc") then
             print *, "Filename trim function not working as intended" 
            stop 1
            
        end if
        
    end subroutine test_get_outfile
    
    !---------------------------------------------------------------------------
    ! Test lambert subroutine
    !---------------------------------------------------------------------------
    subroutine test_lambert()
        real(8) :: stlat1, stlat2, c_lat, c_lon
        real(8) :: glon, glat, x, y
        real(8) :: glon_inv, glat_inv, x_out, y_out
        real(8) :: true_x, true_y
        real(8), parameter :: tol = 1.0e2 ! Difference tolerance
        
        ! Test 1: Forward transformation (glon,glat) -> (x,y)
        ! National mall, Washington DC
        stlat1 = 38.5_8
        stlat2 = 39.5_8
        c_lat = 39.0_8
        c_lon = -77.0_8
        glon = -77.0353_8
        glat = 38.8895_8
        
        true_x = -3055.18
        true_y = -12286.37
        call lambert(stlat1, stlat2, c_lat, c_lon, glon, glat, x, y, 1)

        if ( abs(true_x - x) > tol .or. abs(true_y - y) > tol ) then
          print *, "Forward transformation incorrect"
          stop 2
        end if          
              
        ! Test inverse transformation
        glon_inv = 0.0_8
        glat_inv = 0.0_8
        call lambert(stlat1, stlat2, c_lat, c_lon, glon_inv, glat_inv, x, y, -1)

        if ( (glon - glon_inv) > tol .or. (glat - glat_inv) > tol ) then
          print *, "Results from glon,glat -> lon,lat -> glon,glat not identical"
          stop 3
        end if
        
        ! Test 2: Special case where stlat1 == stlat2
        ! Italy
        stlat1 = 45.5_8
        stlat2 = 45.5_8
        c_lat = 45.5_8
        c_lon = 11.0_8
        glat = 45.4642_8
        glon = 9.1900_8
        
        true_x = -141149.15
        true_y = -2390.66
        
        call lambert(stlat1, stlat2, c_lat, c_lon, glon, glat, x, y, 1)

        if ( (true_x - x) > tol .or. (true_y - y) > tol ) then
           print *, "Issue when stlat1 and stlat2 are equal"
          stop 4
        end if
        
        ! Test 3: Point at projection center
        stlat1 = 38.5_8
        stlat2 = 39.5_8
        c_lat = 39.0_8
        c_lon = -77.0_8
        glon = -77.0_8
        glat = 39.0_8

        call lambert(stlat1, stlat2, c_lat, c_lon, glon, glat, x, y, 1)
      
        if ( abs(x) > 1e-3 .or. abs(y) > 1e-3 ) then
          print *, "Issue when point is at the projection center"
          stop 5
        end if
        
    end subroutine test_lambert
    
    !---------------------------------------------------------------------------
    ! Test rtll subroutine
    !---------------------------------------------------------------------------
    subroutine test_rtll()
        real(8) :: tlmd, tphd, almd, aphd, tlm0d, tph0d
        real(8) :: true_almd, true_aphd
        real(8), parameter :: tol = 1.0e-0
        
        ! Test 1: No rotation
        tlm0d = 0.0_8
        tph0d = 0.0_8
        tlmd = 0.0_8
        tphd = 0.0_8
        true_almd = 0.0_8
        true_aphd = 0.0_8

        call rtll(tlmd, tphd, almd, aphd, tlm0d, tph0d)

        if( abs(true_almd - almd) > tol .or. abs(true_aphd - aphd) > tol) then
           print *, "Coordinates were rotated when no rotation was expected"
          stop 6
        end if
        
        ! Test 2: North pole rotation
        tlm0d = 0.0_8
        tph0d = 90.0_8
        tlmd = 0.0_8
        tphd = 0.0_8
        true_almd = 0.0_8
        true_aphd = 90.0_8
        
        call rtll(tlmd, tphd, almd, aphd, tlm0d, tph0d)

        if( abs(true_almd - almd) > tol .or. abs(true_aphd - aphd) > tol) then
          print *, "Incorrect rotation when 90 deg rotation was expected"
          stop 7
        end if
        
        ! Test 3: Central Europe
        tlm0d = -170.0_8
        tph0d = 40.0_8
        tlmd = 10.0_8
        tphd = 50.0_8
        true_almd = -76.16_8
        true_aphd = 83.57_8
        
        call rtll(tlmd, tphd, almd, aphd, tlm0d, tph0d)

        if( abs(true_almd - almd) > tol .or. abs(true_aphd - aphd) > tol) then
          print *, "Incorrect rotation of coordiantes"
          stop 8
        end if
        
        ! Test 4: Longitude wrapping
        tlm0d = 0.0_8
        tph0d = 45.0_8
        tlmd = 179.0_8
        tphd = 0.0_8
        
        call rtll(tlmd, tphd, almd, aphd, tlm0d, tph0d)

        if( almd > 180 .or. almd < -180) then
          print *, "Incorrect rotation when longitude is at the edge of the domain"
          stop 9
        end if
        
    end subroutine test_rtll
    
end program test_module_wrt_grid_comp
