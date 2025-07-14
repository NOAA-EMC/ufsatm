module module_cdeps_inline
!
!*** this module contains the subroutines for cdeps inline capability 
!
! revision history
!  22 Aug 2023: U. Turuncoglu  Initial development
!
  use ESMF
  use dshr_mod        , only: dshr_pio_init
  use dshr_strdata_mod, only: shr_strdata_type
  use dshr_strdata_mod, only: shr_strdata_init_from_inline
  use dshr_strdata_mod, only: shr_strdata_advance
  use dshr_stream_mod , only: shr_stream_init_from_esmfconfig

  use GFS_typedefs    , only: kp => kind_phys
  use CCPP_data       , only: GFS_control
  use atmos_model_mod , only: setup_inlinedata

  implicit none

  private
  public cdeps_stream_init
  public cdeps_stream_run

  type(ESMF_Grid) :: grid
  type(ESMF_Mesh) :: mesh
  type(ESMF_Field) :: fgrid

  type config
     integer :: year_first
     integer :: year_last
     integer :: year_align
     integer :: offset
     real(kind=8) :: dtlimit
     character(len=ESMF_MAXSTR) :: mesh_filename
     character(len=ESMF_MAXSTR), allocatable :: data_filename(:)
     character(len=ESMF_MAXSTR), allocatable :: fld_list(:)
     character(len=ESMF_MAXSTR), allocatable :: fld_list_model(:)
     character(len=ESMF_MAXSTR) :: mapalgo
     character(len=ESMF_MAXSTR) :: taxmode
     character(len=ESMF_MAXSTR) :: tintalgo
     character(len=ESMF_MAXSTR) :: name
  end type config

  type(config) :: stream ! stream configuration
  type(shr_strdata_type) :: sdat_config
  type(shr_strdata_type), allocatable :: sdat(:) ! input data stream
  real(kind=kp), dimension(:,:), allocatable :: farray !dataptr2d

  integer :: dbug = 0
  integer :: logunit = 6
  real(kind=kp), parameter :: missing_value = 9.99e20_kp
!
  contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

    subroutine cdeps_stream_init(comp, clock, rc)
      ! input/output variables
      type(ESMF_GridComp), intent(in)  :: comp
      type(ESMF_Clock)   , intent(in)  :: clock
      integer            , intent(out) :: rc

      ! local variables
      type(ESMF_ArraySpec) :: arraySpec
      integer :: localPet
      integer :: l, id, nstreams
      integer :: isc, iec, jsc, jec
      character(len=ESMF_MAXSTR) :: streamfilename, stream_name
      character(len=ESMF_MAXSTR), allocatable :: file_list(:), var_list(:,:)

      ! query compontn to retrieve required information
      call ESMF_GridCompGet(comp, grid=grid, localPet=localPet, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      ! create mesh from grid
      mesh = ESMF_MeshCreate(grid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      ! init pio
      call dshr_pio_init(comp, sdat_config, logunit, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      ! read stream configuration file
      streamfilename = 'stream.config' 
      call shr_stream_init_from_esmfconfig(streamfilename, sdat_config%stream, logunit, &
          sdat_config%pio_subsystem, sdat_config%io_type, sdat_config%io_format, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      ! get number of streams
      nstreams = size(sdat_config%stream)

      ! allocate stream data type
      if (.not. allocated(sdat)) allocate(sdat(nstreams))

      ! loop over streams and init
      do id = 1, nstreams
         ! set pio related variables 
         sdat(id)%pio_subsystem => sdat_config%pio_subsystem
         sdat(id)%io_type = sdat_config%io_type
         sdat(id)%io_format = sdat_config%io_format

         ! allocate temporary arrays
         allocate(file_list(sdat_config%stream(id)%nfiles)) 
         allocate(var_list(sdat_config%stream(id)%nvars,2))

         ! fill variables
         do l = 1, sdat_config%stream(id)%nfiles
            file_list(l) = trim(sdat_config%stream(id)%file(l)%name)
         end do
         do l = 1, sdat_config%stream(id)%nvars
            var_list(l,1) = trim(sdat_config%stream(id)%varlist(l)%nameinfile)
            var_list(l,2) = trim(sdat_config%stream(id)%varlist(l)%nameinmodel)
         end do

         ! init stream
         write(stream_name,fmt='(a,i2.2)') 'stream_', id
         call shr_strdata_init_from_inline(sdat(id), &
            my_task=localPet, logunit=logunit, &
            compname = 'cmeps', model_clock=clock, model_mesh=mesh, &
            stream_meshfile=trim(sdat_config%stream(id)%meshfile), &
            stream_filenames=file_list, &
            stream_yearFirst=sdat_config%stream(id)%yearFirst, &
            stream_yearLast=sdat_config%stream(id)%yearLast, &
            stream_yearAlign=sdat_config%stream(id)%yearAlign, &
            stream_fldlistFile=var_list(:,1), &
            stream_fldListModel=var_list(:,2), &
            stream_lev_dimname=trim(sdat_config%stream(id)%lev_dimname), &
            stream_mapalgo=trim(sdat_config%stream(id)%mapalgo), &
            stream_offset=sdat_config%stream(id)%offset, &
            stream_taxmode=trim(sdat_config%stream(id)%taxmode), &
            stream_dtlimit=sdat_config%stream(id)%dtlimit, &
            stream_tintalgo=trim(sdat_config%stream(id)%tInterpAlgo), &
            stream_name=trim(stream_name), &
            stream_src_mask=sdat_config%stream(id)%src_mask_val, &
            stream_dst_mask=sdat_config%stream(id)%dst_mask_val, &
            rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         ! clean memory
         deallocate(file_list)
         deallocate(var_list)
      end do

      ! create temporary field on grid
      isc = GFS_control%isc
      iec = GFS_control%isc+GFS_control%nx-1
      jsc = GFS_control%jsc
      jec = GFS_control%jsc+GFS_control%ny-1
      allocate(farray(isc:iec,jsc:jec))
      fgrid = ESMF_FieldCreate(grid=grid, farray=farray, indexflag=ESMF_INDEX_DELOCAL, name='noname', rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    end subroutine cdeps_stream_init

  !-----------------------------------------------------------------------------

    subroutine cdeps_stream_run(clock, rc)
      ! input/output variables
      type(ESMF_Clock), intent(in)  :: clock
      integer, intent(out) :: rc

      ! local variables
      integer :: item, id, nstreams, nflds
      integer :: curr_ymd, sec
      integer :: year, month, day, hour, minute, second
      character(len=ESMF_MAXSTR) :: filename, istr
      type(ESMF_Time)  :: currTime
      type(ESMF_Field) :: fmesh
      type(ESMF_RouteHandle), save :: rh
      real(kind=kp), dimension(:,:), pointer  :: dataptr2d

      ! query clock
      call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      ! get current time
      call ESMF_TimeGet(currTime, yy=year, mm=month, dd=day, h=hour, m=minute, s=second, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      curr_ymd = abs(year)*10000+month*100+day
      sec = hour*3600+minute*60+second

      ! get number of streams
      nstreams = size(sdat)

      ! loop over streams and get data 
      do id = 1, nstreams
         ! advance cdeps inline
         write(istr,fmt='(a,i2.2)') 'stream_', id
         call shr_strdata_advance(sdat(id), ymd=curr_ymd, tod=sec, logunit=logunit, istr=trim(istr), rc=rc)      
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

         ! get number of fields in FB
         nflds = size(sdat(id)%pstrm(1)%fldlist_model)

         ! loop over fields
         do item = 1, nflds
            ! get field on mesh
            call ESMF_FieldBundleGet(sdat(id)%pstrm(1)%fldbun_model, fieldName=trim(sdat(id)%pstrm(1)%fldlist_model(item)), field=fmesh, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ! create RH and destination field to transfer data from mesh to grid
            if (.not. ESMF_RouteHandleIsCreated(rh, rc=rc)) then
               ! create RH
               call ESMF_FieldRedistStore(fmesh, fgrid, rh, ignoreUnmatchedIndices=.true., rc=rc) 
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            end if

            ! initialize destination field
            call ESMF_FieldFill(fgrid, dataFillScheme="const", const1=missing_value, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ! transfer data from mesh to grid
            call ESMF_FieldRedist(fmesh, fgrid, rh, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ! get field
            call ESMF_FieldGet(fgrid, farrayPtr=dataptr2d, localDE=0, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ! fill internal data structures
            call setup_inlinedata(trim(sdat(id)%pstrm(1)%fldlist_model(item)), dataptr2d, logunit)

            ! diagnostic output
            if (dbug > 0) then
               ! write statistics
               write(logunit,'(A,3g14.7,i8)') '(cdeps inline): '//trim(sdat(id)%pstrm(1)%fldlist_model(item))//' ', &
                  minval(dataptr2d), maxval(dataptr2d), sum(dataptr2d), size(dataptr2d)
            end if

            if (dbug > 5) then
               ! file name
               write(filename, fmt='(a,i4,a1,i2.2,a1,i2.2,a1,i5.5)') trim(sdat(id)%pstrm(1)%fldlist_model(item))//'_', &
                  year, '-', month, '-', day, '-', sec

               ! write field
               if (dbug > 10) then
                  ! write field on mesh to VTK
                  call ESMF_FieldWriteVTK(fmesh, trim(filename), rc=rc)
                  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
               else
                  ! write field on grid to netCDF
                  call ESMF_FieldWrite(fgrid, fileName=trim(filename)//'.nc', rc=rc)
                  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
               end if
            end if
         end do
      end do

    end subroutine cdeps_stream_run

end module module_cdeps_inline
