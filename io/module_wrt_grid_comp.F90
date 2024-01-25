!-----------------------------------------------------------------------
!
    module module_wrt_grid_comp
!
!-----------------------------------------------------------------------
!***  This module includes the functionality of write gridded component.
!-----------------------------------------------------------------------
!***  At initialization step, write grid is defined. The forecast field
!***  bundle is mirrored and output field information inside the field
!***  bundle is used to create ESMF field on the write grid and added in
!***  the output field bundle on write grid component. Also the IO_BaseTime
!***  is set to the initial clock time.
!***  At the run step, output time is set from the write grid comp clock
!***  the ESMF field bundles that contains the data on write grid are
!***  written out through ESMF field bundle write to netcdf files.
!***  The ESMF field bundle write uses parallel write, so if output grid
!***  is cubed sphere grid, the  six tiles file will be written out at
!***  same time.
!-----------------------------------------------------------------------
!***
!***  Revision history
!***
!     Jul 2017:  J. Wang/G. Theurich  - initial code for fv3 write grid component
!     Mar 2018:  S  Moorthi           - changing cfhour to accommodate up to 99999 hours
!     Aug 2019:  J. Wang              - add inline post
!
!---------------------------------------------------------------------------------
!
      use mpi
      use esmf
      use fms

      use write_internal_state
      use module_fv3_io_def,   only : num_pes_fcst,                             &
                                      n_group, num_files,                       &
                                      filename_base, output_grid, output_file,  &
                                      imo,jmo,ichunk2d,jchunk2d,                &
                                      ichunk3d,jchunk3d,kchunk3d,               &
                                      quantize_mode,quantize_nsd,               &
                                      cen_lon, cen_lat,                         &
                                      lon1, lat1, lon2, lat2, dlon, dlat,       &
                                      stdlat1, stdlat2, dx, dy, iau_offset,     &
                                      ideflate, zstandard_level, lflname_fulltime
      use module_write_netcdf, only : write_netcdf
      use module_write_restart_netcdf, only : write_restart_netcdf
      use physcons,            only : pi => con_pi
#ifdef INLINE_POST
      use post_fv3,            only : post_run_fv3
#endif
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------

      private
!
!-----------------------------------------------------------------------
!
!
      integer,save      :: lead_write_task                                !<-- Rank of the first write task in the write group
      integer,save      :: last_write_task                                !<-- Rank of the last write task in the write group
      integer,save      :: ntasks                                         !<-- # of write tasks in the current group
      integer,save      :: itasks, jtasks                                 !<-- # of write tasks in i/j direction in the current group
      integer,save      :: ngrids

      integer,save      :: wrt_mpi_comm                                   !<-- the mpi communicator in the write comp
      integer,save      :: idate(7), start_time(7)
      logical,save      :: write_nsflip
      logical,save      :: change_wrtidate=.false.
      integer,save      :: frestart(999) = -1
      integer,save      :: calendar_type = 3
      logical           :: lprnt
!
!-----------------------------------------------------------------------
!
      type(ESMF_FieldBundle)           :: gridFB
      integer                          :: FBCount
      character(len=esmf_maxstr),allocatable    :: fcstItemNameList(:)
      logical                                   :: top_parent_is_global
!
!-----------------------------------------------------------------------
      REAL(KIND=8)             :: btim,btim0
      REAL(KIND=8),PUBLIC,SAVE :: write_init_tim, write_run_tim
      REAL(KIND=8), parameter  :: radi=180.0d0/pi
!-----------------------------------------------------------------------
!
      public SetServices
!
      interface splat
        module procedure splat4
        module procedure splat8
      end interface splat
!
      type optimizeT
        type(ESMF_State)                   :: state
        type(ESMF_GridComp), allocatable   :: comps(:)
      end type

      contains
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine SetServices(wrt_comp, rc)
        type(ESMF_GridComp)  :: wrt_comp
        integer, intent(out) :: rc

        rc = ESMF_SUCCESS

        call ESMF_GridCompSetEntryPoint(wrt_comp, ESMF_METHOD_INITIALIZE, phase=1, &
                                        userRoutine=wrt_initialize_p1, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridCompSetEntryPoint(wrt_comp, ESMF_METHOD_INITIALIZE, phase=2, &
                                        userRoutine=wrt_initialize_p2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridCompSetEntryPoint(wrt_comp, ESMF_METHOD_INITIALIZE, phase=3, &
                                        userRoutine=wrt_initialize_p3, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridCompSetEntryPoint(wrt_comp, ESMF_METHOD_RUN, &
                                        userRoutine=wrt_run, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridCompSetEntryPoint(wrt_comp, ESMF_METHOD_FINALIZE, &
                                        userRoutine=wrt_finalize, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      end subroutine SetServices
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine wrt_initialize_p1(wrt_comp, imp_state_write, exp_state_write, clock, rc)
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE WRITE GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      type(esmf_GridComp)               :: wrt_comp
      type(ESMF_State)                  :: imp_state_write, exp_state_write
      type(esmf_Clock)                  :: clock
      integer,intent(out)               :: rc
!
!***  LOCAL VARIABLES
!
      TYPE(ESMF_VM)                           :: VM
      type(write_wrap)                        :: WRAP
      type(wrt_internal_state),pointer        :: wrt_int_state

      integer                                 :: tl, i, j, n, k
      integer,dimension(2,6)                  :: decomptile
      integer,dimension(2)                    :: regDecomp !define delayout for the nest grid
      integer                                 :: fieldCount
      integer                                 :: vm_mpi_comm
      character(40)                           :: fieldName
      type(ESMF_Config)                       :: cf, cf_output_grid
      type(ESMF_Info)                         :: info
      type(ESMF_DELayout)                     :: delayout
      type(ESMF_Grid)                         :: fcstGrid
      type(ESMF_Grid), allocatable            :: wrtGrid(:)
      type(ESMF_Grid)                         :: wrtGrid_cubed_sphere
      logical                                 :: create_wrtGrid_cubed_sphere = .true.
      type(ESMF_Grid)                         :: actualWrtGrid
      type(ESMF_Array)                        :: array
      type(ESMF_Field)                        :: field_work, field
      type(ESMF_Decomp_Flag)                  :: decompflagPTile(2,6)

      type(ESMF_StateItem_Flag), allocatable  :: fcstItemTypeList(:)
      type(ESMF_FieldBundle)                  :: fcstFB, fieldbundle, mirrorFB
      type(ESMF_Field),          allocatable  :: fcstField(:)
      type(ESMF_TypeKind_Flag)                :: typekind
      character(len=80),         allocatable  :: fieldnamelist(:)
      integer                                 :: fieldDimCount, gridDimCount, tk, sloc
      integer,                   allocatable  :: petMap(:)
      integer,                   allocatable  :: gridToFieldMap(:)
      integer,                   allocatable  :: ungriddedLBound(:)
      integer,                   allocatable  :: ungriddedUBound(:)
      type(ESMF_StaggerLoc)                   :: staggerloc
      character(len=80)                       :: attName
      character(len=80),         allocatable  :: attNameList(:),attNameList2(:)
      type(ESMF_TypeKind_Flag),  allocatable  :: typekindList(:)
      character(len=80)                       :: valueS
      integer                                 :: valueI4
      real(ESMF_KIND_R4)                      :: valueR4
      real(ESMF_KIND_R8)                      :: valueR8
      logical, allocatable                    :: is_moving(:)
      logical                                 :: isPresent

      integer :: attCount, jidx, idx, noutfile
      character(19)  :: newdate
      character(128) :: FBlist_outfilename(100), outfile_name
      character(128),dimension(:,:), allocatable    :: outfilename
      real(8), dimension(:),         allocatable    :: slat
      real(8), dimension(:),         allocatable    :: lat, lon
      real(ESMF_KIND_R8), dimension(:,:), pointer   :: lonPtr, latPtr
      real(ESMF_KIND_R8)                            :: rot_lon, rot_lat
      real(ESMF_KIND_R8)                            :: geo_lon, geo_lat
      real(ESMF_KIND_R8)                            :: lon1_r8, lat1_r8
      real(ESMF_KIND_R8)                            :: x1, y1, x, y, delat, delon
      type(ESMF_TimeInterval)                       :: IAU_offsetTI

      character(256)                          :: cf_open, cf_close
      character(256)                          :: gridfile
      integer                                 :: num_output_file

      type(ESMF_DistGrid)                     :: acceptorDG, newAcceptorDG
      integer                                 :: grid_id

      logical                    :: history_file_on_native_grid
      character(len=esmf_maxstr) :: output_grid_name
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  initialize the write component timers.
!-----------------------------------------------------------------------
!
      write_init_tim = 0.
      write_run_tim  = 0.
      btim0 = MPI_Wtime()
!
!-----------------------------------------------------------------------
!***  set the write component's internal state.
!-----------------------------------------------------------------------
!
      allocate(wrt_int_state,stat=RC)
      wrap%write_int_state => wrt_int_state
      call ESMF_GridCompSetInternalState(wrt_comp, wrap, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
      call ESMF_VMGetCurrent(vm=VM,rc=RC)
      call ESMF_VMGet(vm=VM, localPet=wrt_int_state%mype,               &
                      petCount=wrt_int_state%petcount,mpiCommunicator=vm_mpi_comm,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call mpi_comm_dup(vm_mpi_comm, wrt_mpi_comm, rc)

      ntasks = wrt_int_state%petcount
      jidx   = wrt_int_state%petcount/6
      lead_write_task = 0
      last_write_task = ntasks -1
      lprnt = lead_write_task == wrt_int_state%mype

      call fms_init(wrt_mpi_comm)

!      print *,'in wrt, lead_write_task=', &
!         lead_write_task,'last_write_task=',last_write_task, &
!         'mype=',wrt_int_state%mype,'jidx=',jidx,' comm=',wrt_mpi_comm
!

!-----------------------------------------------------------------------
!*** get configuration variables
!-----------------------------------------------------------------------
!
      call ESMF_GridCompGet(gridcomp=wrt_comp,config=CF,rc=RC)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out
! variables for post
      call ESMF_ConfigGetAttribute(config=CF,value=wrt_int_state%output_history,default=.true., &
                                   label='output_history:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call ESMF_ConfigGetAttribute(config=CF,value=wrt_int_state%write_dopost,default=.false., &
                                   label='write_dopost:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      call ESMF_ConfigGetAttribute(config=CF,value=write_nsflip,default=.false., &
                                   label='write_nsflip:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      if( wrt_int_state%write_dopost ) then
#ifdef INLINE_POST
        call ESMF_ConfigGetAttribute(config=CF,value=wrt_int_state%post_namelist,default='itag', &
                                     label ='post_namelist:',rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
#else
        rc = ESMF_RC_NOT_IMPL
        print *,'inline post not available on this machine'
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
#endif
      endif

      allocate(output_file(num_files))
      num_output_file = ESMF_ConfigGetLen(config=CF, label ='output_file:',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      if (num_files == num_output_file) then
        call ESMF_ConfigGetAttribute(CF,valueList=output_file,label='output_file:', &
             count=num_files, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        do i = 1, num_files
          if(output_file(i) /= "netcdf" .and. output_file(i) /= "netcdf_parallel") then
            write(0,*)"Only netcdf and netcdf_parallel are allowed for multiple values of output_file"
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
          endif
        enddo
      else if ( num_output_file == 1) then
        call ESMF_ConfigGetAttribute(CF,valuelist=output_file,label='output_file:', count=1, rc=rc)
        output_file(1:num_files) = output_file(1)
      else
        output_file(1:num_files) = 'netcdf'
      endif
      if(lprnt) then
        print *,'num_files=',num_files
        do i=1,num_files
          print *,'num_file=',i,'filename_base= ',trim(filename_base(i)),' output_file= ',trim(output_file(i))
        enddo
      endif

      call ESMF_InfoGetFromHost(imp_state_write, info=info, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_InfoGetAlloc(info, key="is_moving", values=is_moving, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeGet(imp_state_write, convention="NetCDF", purpose="FV3", &
                             name="ngrids", value=ngrids, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeGet(imp_state_write, convention="NetCDF", purpose="FV3", &
                             name="top_parent_is_global", value=top_parent_is_global, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      allocate(wrtGrid(ngrids))

      allocate(output_grid(ngrids))

      allocate(imo(ngrids))
      allocate(jmo(ngrids))

      allocate(cen_lon(ngrids))
      allocate(cen_lat(ngrids))
      allocate(lon1(ngrids))
      allocate(lat1(ngrids))
      allocate(lon2(ngrids))
      allocate(lat2(ngrids))
      allocate(dlon(ngrids))
      allocate(dlat(ngrids))

      allocate(stdlat1(ngrids))
      allocate(stdlat2(ngrids))
      allocate(dx(ngrids))
      allocate(dy(ngrids))

      allocate(ichunk2d(ngrids))
      allocate(jchunk2d(ngrids))
      allocate(ichunk3d(ngrids))
      allocate(jchunk3d(ngrids))
      allocate(kchunk3d(ngrids))
      allocate(ideflate(ngrids))
      allocate(quantize_mode(ngrids))
      allocate(quantize_nsd(ngrids))
      allocate(zstandard_level(ngrids))

      allocate(wrt_int_state%out_grid_info(ngrids))

      do n=1, ngrids

        if (n == 1) then
          ! for top level domain look directly in cf
          cf_output_grid = cf
        else
          ! for nest domains, look under specific section
          write(cf_open,'("<output_grid_",I2.2,">")') n
          write(cf_close,'("</output_grid_",I2.2,">")') n
          cf_output_grid = ESMF_ConfigCreate(cf, openLabel=trim(cf_open), closeLabel=trim(cf_close), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        end if

        call ESMF_ConfigGetAttribute(config=cf_output_grid, value=output_grid(n), label ='output_grid:',rc=rc)
        if (lprnt) then
          print *,'grid_id= ', n, ' output_grid= ', trim(output_grid(n))
        end if

        call ESMF_ConfigGetAttribute(config=CF, value=itasks,default=1,label ='itasks:',rc=rc)
        jtasks = ntasks
        if(itasks > 0 ) jtasks = ntasks/itasks
        if( itasks*jtasks /= ntasks ) then
          itasks = 1
          jtasks = ntasks
        endif

        if (trim(output_grid(n)) == 'gaussian_grid' .or. trim(output_grid(n)) == 'global_latlon') then
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=imo(n), label ='imo:',rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=jmo(n), label ='jmo:',rc=rc)
          if (lprnt) then
            print *,'imo=',imo(n),'jmo=',jmo(n)
          end if
        else if (trim(output_grid(n)) == 'regional_latlon') then
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=lon1(n), label ='lon1:',rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=lat1(n), label ='lat1:',rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=lon2(n), label ='lon2:',rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=lat2(n), label ='lat2:',rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=dlon(n), label ='dlon:',rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=dlat(n), label ='dlat:',rc=rc)
          imo(n) = (lon2(n)-lon1(n))/dlon(n) + 1
          jmo(n) = (lat2(n)-lat1(n))/dlat(n) + 1
          if (lprnt) then
            print *,'lon1=',lon1(n),' lat1=',lat1(n)
            print *,'lon2=',lon2(n),' lat2=',lat2(n)
            print *,'dlon=',dlon(n),' dlat=',dlat(n)
            print *,'imo =',imo(n), ' jmo =',jmo(n)
          end if
        else if (trim(output_grid(n)) == 'rotated_latlon') then
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=cen_lon(n), label ='cen_lon:',rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=cen_lat(n), label ='cen_lat:',rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=lon1(n),    label ='lon1:',   rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=lat1(n),    label ='lat1:',   rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=lon2(n),    label ='lon2:',   rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=lat2(n),    label ='lat2:',   rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=dlon(n),    label ='dlon:',   rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=dlat(n),    label ='dlat:',   rc=rc)
          imo(n) = (lon2(n)-lon1(n))/dlon(n) + 1
          jmo(n) = (lat2(n)-lat1(n))/dlat(n) + 1
          if (lprnt) then
            print *,'cen_lon=',cen_lon(n),' cen_lat=',cen_lat(n)
            print *,'lon1   =',lon1(n),   ' lat1   =',lat1(n)
            print *,'lon2   =',lon2(n),   ' lat2   =',lat2(n)
            print *,'dlon   =',dlon(n),   ' dlat   =',dlat(n)
            print *,'imo    =',imo(n),    ' jmo    =',jmo(n)
          end if
        else if (trim(output_grid(n)) == 'rotated_latlon_moving' .or. &
                 trim(output_grid(n)) == 'regional_latlon_moving') then
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=imo(n),  label ='imo:', rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=jmo(n),  label ='jmo:', rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=dlon(n), label ='dlon:',rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=dlat(n), label ='dlat:',rc=rc)
          if (lprnt) then
            print *,'imo =',imo(n), ' jmo =',jmo(n)
            print *,'dlon=',dlon(n),' dlat=',dlat(n)
          end if
        else if (trim(output_grid(n)) == 'lambert_conformal') then
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=cen_lon(n), label ='cen_lon:',rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=cen_lat(n), label ='cen_lat:',rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=stdlat1(n), label ='stdlat1:',rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=stdlat2(n), label ='stdlat2:',rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=imo(n),     label ='nx:',     rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=jmo(n),     label ='ny:',     rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=lon1(n),    label ='lon1:',   rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=lat1(n),    label ='lat1:',   rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=dx(n),      label ='dx:',     rc=rc)
          call ESMF_ConfigGetAttribute(config=cf_output_grid, value=dy(n),      label ='dy:',     rc=rc)
          if (lprnt) then
            print *,'cen_lon=',cen_lon(n),' cen_lat=',cen_lat(n)
            print *,'stdlat1=',stdlat1(n),' stdlat2=',stdlat2(n)
            print *,'lon1=',lon1(n),' lat1=',lat1(n)
            print *,'nx=',imo(n), ' ny=',jmo(n)
            print *,'dx=',dx(n),' dy=',dy(n)
          endif
        endif ! output_grid

        ! chunksizes for netcdf_parallel
        call ESMF_ConfigGetAttribute(config=CF,value=ichunk2d(n),default=0,label ='ichunk2d:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF,value=jchunk2d(n),default=0,label ='jchunk2d:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF,value=ichunk3d(n),default=0,label ='ichunk3d:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF,value=jchunk3d(n),default=0,label ='jchunk3d:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF,value=kchunk3d(n),default=0,label ='kchunk3d:',rc=rc)

        ! zstandard compression flag
        call ESMF_ConfigGetAttribute(config=CF,value=zstandard_level(n),default=0,label ='zstandard_level:',rc=rc)
        if (zstandard_level(n) < 0) zstandard_level(n)=0

        ! zlib compression flag
        call ESMF_ConfigGetAttribute(config=CF,value=ideflate(n),default=0,label ='ideflate:',rc=rc)
        if (ideflate(n) < 0) ideflate(n)=0

        if (ideflate(n) > 0 .and. zstandard_level(n) > 0) then
           write(0,*)"wrt_initialize_p1: zlib and zstd compression cannot be both enabled at the same time"
           call ESMF_LogWrite("wrt_initialize_p1: zlib and zstd compression cannot be both enabled at the same time",ESMF_LOGMSG_ERROR,rc=RC)
           call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if

        ! quantize_mode and quantize_nsd
        call ESMF_ConfigGetAttribute(config=CF,value=quantize_mode(n),default='quantize_bitgroom',label='quantize_mode:',rc=rc)
        call ESMF_ConfigGetAttribute(config=CF,value=quantize_nsd(n),default=0,label='quantize_nsd:',rc=rc)

        if (.NOT. (trim(quantize_mode(n))=='quantize_bitgroom' &
              .OR. trim(quantize_mode(n))=='quantize_granularbr' &
              .OR. trim(quantize_mode(n))=='quantize_bitround') ) then
           write(0,*)"wrt_initialize_p1: unknown quantize_mode ", trim(quantize_mode(n))
           call ESMF_LogWrite("wrt_initialize_p1: wrt_initialize_p1: unknown quantize_mode "//trim(quantize_mode(n)),ESMF_LOGMSG_ERROR,rc=RC)
           call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if

        if (lprnt) then
            print *,'ideflate=',ideflate(n)
            print *,'quantize_mode=',trim(quantize_mode(n)),' quantize_nsd=',quantize_nsd(n)
            print *,'zstandard_level=',zstandard_level(n)
        end if

        if (cf_output_grid /= cf) then
          ! destroy the temporary config object created for nest domains
          call ESMF_ConfigDestroy(config=cf_output_grid, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        endif

    call ESMF_ConfigGetAttribute(config=CF, value=history_file_on_native_grid, default=.false., &
                                 label='history_file_on_native_grid:', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

#if 1
        if (n == 1 .and. top_parent_is_global .and. history_file_on_native_grid) then
          do tl=1,6
            decomptile(1,tl) = 1
            decomptile(2,tl) = jidx
            decompflagPTile(:,tl) = (/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/)
          enddo
          call ESMF_AttributeGet(imp_state_write, convention="NetCDF", purpose="FV3", &
                                 name="gridfile", value=gridfile, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          wrtGrid_cubed_sphere = ESMF_GridCreateMosaic(filename="INPUT/"//trim(gridfile),                                 &
                                                       regDecompPTile=decomptile,tileFilePath="INPUT/",                   &
                                                       decompflagPTile=decompflagPTile,                                   &
                                                       staggerlocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
                                                       name='wrt_grid', rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          create_wrtGrid_cubed_sphere = .false.
        endif
#endif

        if ( trim(output_grid(n)) == 'cubed_sphere_grid' ) then
          !*** Create cubed sphere grid from file
          if (top_parent_is_global .and. n == 1) then
            do tl=1,6
              decomptile(1,tl) = 1
              decomptile(2,tl) = jidx
              decompflagPTile(:,tl) = (/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/)
            enddo
            call ESMF_AttributeGet(imp_state_write, convention="NetCDF", purpose="FV3", &
                                   name="gridfile", value=gridfile, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            wrtGrid(n) = ESMF_GridCreateMosaic(filename="INPUT/"//trim(gridfile),                              &
                                            regDecompPTile=decomptile,tileFilePath="INPUT/",                   &
                                            decompflagPTile=decompflagPTile,                                   &
                                            staggerlocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
                                            name='wrt_grid', rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          else
            if (top_parent_is_global) then
              write(gridfile,'(A,I2.2,A,I1,A)') 'grid.nest', n, '.tile', n+5, '.nc'
            else
              if (n == 1) then
                gridfile='grid.tile7.halo0.nc'   ! regional top-level parent
              else
                write(gridfile,'(A,I2.2,A,I1,A)') 'grid.nest', n, '.tile', n, '.nc'
              endif
            end if
            regDecomp(1) = 1
            regDecomp(2) = ntasks
            allocate(petMap(ntasks))
            do i=1, ntasks
              petMap(i) = i-1
            enddo
            delayout = ESMF_DELayoutCreate(petMap=petMap, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ! create the nest Grid by reading it from file but use DELayout
            call ESMF_LogWrite("wrtComp: gridfile:"//trim(gridfile),ESMF_LOGMSG_INFO,rc=rc)
            wrtGrid(n) = ESMF_GridCreate(filename="INPUT/"//trim(gridfile),                                    &
                                      fileformat=ESMF_FILEFORMAT_GRIDSPEC, regDecomp=regDecomp,                &
                                      decompflag=(/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/),          &
                                      delayout=delayout, isSphere=.false., indexflag=ESMF_INDEX_DELOCAL,       &
                                      rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            deallocate(petMap)
          endif
        else  ! non 'cubed_sphere_grid'
          if ( trim(output_grid(n)) == 'gaussian_grid') then

            wrtGrid(n) = ESMF_GridCreate1PeriDim(minIndex=(/1,1/),                                        &
                                                 maxIndex=(/imo(n),jmo(n)/), regDecomp=(/itasks,jtasks/), &
                                                 indexflag=ESMF_INDEX_GLOBAL,                             &
                                                 name='wrt_grid',rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_GridAddCoord(wrtGrid(n), staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_GridGetCoord(wrtGrid(n), coordDim=1, farrayPtr=lonPtr, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_GridGetCoord(wrtGrid(n), coordDim=2, farrayPtr=latPtr, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            allocate(slat(jmo(n)), lat(jmo(n)), lon(imo(n)))
            call splat(4, jmo(n), slat)
            if(write_nsflip) then
              do j=1,jmo(n)
                lat(j) = asin(slat(j)) * radi
              enddo
            else
              do j=1,jmo(n)
                lat(jmo(n)-j+1) = asin(slat(j)) * radi
              enddo
            endif
            do j=1,imo(n)
              lon(j) = 360.d0/real(imo(n),8) *real(j-1,8)
            enddo
            do j=lbound(latPtr,2),ubound(latPtr,2)
              do i=lbound(lonPtr,1),ubound(lonPtr,1)
                lonPtr(i,j) = 360.d0/real(imo(n),8) * real(i-1,8)
                latPtr(i,j) = lat(j)
              enddo
            enddo
            lon1(n) = lon(1)
            lon2(n) = lon(imo(n))
            lat1(n) = lat(1)
            lat2(n) = lat(jmo(n))
            dlon(n) =  360.d0/real(imo(n),8)
            dlat(n) =  180.d0/real(jmo(n),8)

            deallocate(slat, lat, lon)

          else if ( trim(output_grid(n)) == 'global_latlon') then
            wrtGrid(n) = ESMF_GridCreate1PeriDim(minIndex=(/1,1/),                                        &
                                                 maxIndex=(/imo(n),jmo(n)/), regDecomp=(/itasks,jtasks/), &
                                                 indexflag=ESMF_INDEX_GLOBAL,                             &
                                                 name='wrt_grid',rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_GridAddCoord(wrtGrid(n), staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_GridGetCoord(wrtGrid(n), coordDim=1, farrayPtr=lonPtr, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_GridGetCoord(wrtGrid(n), coordDim=2, farrayPtr=latPtr, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            allocate(lat(jmo(n)), lon(imo(n)))
            if (mod(jmo(n),2) == 0) then
              ! if jmo even, lats do not include poles and equator
              delat = 180.d0/real(jmo(n),8)
              if(write_nsflip) then
                do j=1,jmo(n)
                  lat(j) = 90.d0 - 0.5*delat - real(j-1,8)*delat
                enddo
              else
                do j=1,jmo(n)
                  lat(j) = -90.d0 + 0.5*delat + real(j-1,8)*delat
                enddo
              endif
            else
              ! if jmo odd, lats include poles and equator
              delat = 180.d0/real(jmo(n)-1,8)
              if(write_nsflip) then
                do j=1,jmo(n)
                  lat(j) = 90.d0 - real(j-1,8)*delat
                enddo
              else
                do j=1,jmo(n)
                  lat(j) = -90.d0 + real(j-1,8)*delat
                enddo
              endif
            endif
            delon = 360.d0/real(imo(n),8)
            do i=1,imo(n)
              lon(i) = real(i-1,8)*delon
            enddo
            do j=lbound(latPtr,2),ubound(latPtr,2)
              do i=lbound(lonPtr,1),ubound(lonPtr,1)
                lonPtr(i,j) = lon(i)
                latPtr(i,j) = lat(j)
              enddo
            enddo
            lon1(n) = lon(1)
            lon2(n) = lon(imo(n))
            lat1(n) = lat(1)
            lat2(n) = lat(jmo(n))
            dlon(n) = delon
            dlat(n) = delat

            deallocate(lat, lon)

          else if ( trim(output_grid(n)) == 'regional_latlon' .or.        &
                    trim(output_grid(n)) == 'regional_latlon_moving' .or. &
                    trim(output_grid(n)) == 'rotated_latlon' .or.         &
                    trim(output_grid(n)) == 'rotated_latlon_moving' .or.  &
                    trim(output_grid(n)) == 'lambert_conformal' ) then

            wrtGrid(n) = ESMF_GridCreateNoPeriDim(minIndex=(/1,1/),                                        &
                                                  maxIndex=(/imo(n),jmo(n)/), regDecomp=(/itasks,jtasks/), &
                                                  indexflag=ESMF_INDEX_GLOBAL,                             &
                                                  name='wrt_grid',rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_GridAddCoord(wrtGrid(n), staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_GridGetCoord(wrtGrid(n), coordDim=1, farrayPtr=lonPtr, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_GridGetCoord(wrtGrid(n), coordDim=2, farrayPtr=latPtr, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            if ( trim(output_grid(n)) == 'regional_latlon' ) then
                do j=lbound(lonPtr,2),ubound(lonPtr,2)
                do i=lbound(lonPtr,1),ubound(lonPtr,1)
                  lonPtr(i,j) = lon1(n) + (lon2(n)-lon1(n))/(imo(n)-1) * (i-1)
                  latPtr(i,j) = lat1(n) + (lat2(n)-lat1(n))/(jmo(n)-1) * (j-1)
                enddo
                enddo
            else if ( trim(output_grid(n)) == 'regional_latlon_moving' ) then
                ! Do not compute lonPtr, latPtr here. Will be done in the run phase
            else if ( trim(output_grid(n)) == 'rotated_latlon' ) then
                do j=lbound(lonPtr,2),ubound(lonPtr,2)
                do i=lbound(lonPtr,1),ubound(lonPtr,1)
                  rot_lon = lon1(n) + (lon2(n)-lon1(n))/(imo(n)-1) * (i-1)
                  rot_lat = lat1(n) + (lat2(n)-lat1(n))/(jmo(n)-1) * (j-1)
                  call rtll(rot_lon, rot_lat, geo_lon, geo_lat, dble(cen_lon(n)), dble(cen_lat(n)))
                  if (geo_lon < 0.0) geo_lon = geo_lon + 360.0
                  lonPtr(i,j) = geo_lon
                  latPtr(i,j) = geo_lat
                enddo
                enddo
                rot_lon = lon1(n)
                rot_lat = lat1(n)
                call rtll(rot_lon, rot_lat, geo_lon, geo_lat, dble(cen_lon(n)), dble(cen_lat(n)))
                if (geo_lon < 0.0) geo_lon = geo_lon + 360.0
                wrt_int_state%out_grid_info(n)%lonstart = geo_lon
                wrt_int_state%out_grid_info(n)%latstart = geo_lat

                rot_lon = lon2(n)
                rot_lat = lat1(n)
                call rtll(rot_lon, rot_lat, geo_lon, geo_lat, dble(cen_lon(n)), dble(cen_lat(n)))
                if (geo_lon < 0.0) geo_lon = geo_lon + 360.0
                wrt_int_state%out_grid_info(n)%lonse = geo_lon
                wrt_int_state%out_grid_info(n)%latse = geo_lat

                rot_lon = lon1(n)
                rot_lat = lat2(n)
                call rtll(rot_lon, rot_lat, geo_lon, geo_lat, dble(cen_lon(n)), dble(cen_lat(n)))
                if (geo_lon < 0.0) geo_lon = geo_lon + 360.0
                wrt_int_state%out_grid_info(n)%lonnw = geo_lon
                wrt_int_state%out_grid_info(n)%latnw = geo_lat

                rot_lon = lon2(n)
                rot_lat = lat2(n)
                call rtll(rot_lon, rot_lat, geo_lon, geo_lat, dble(cen_lon(n)), dble(cen_lat(n)))
                if (geo_lon < 0.0) geo_lon = geo_lon + 360.0
                wrt_int_state%out_grid_info(n)%lonlast = geo_lon
                wrt_int_state%out_grid_info(n)%latlast = geo_lat
            else if ( trim(output_grid(n)) == 'rotated_latlon_moving' ) then
                ! Do not compute lonPtr, latPtr here. Will be done in the run phase
            else if ( trim(output_grid(n)) == 'lambert_conformal' ) then
                lon1_r8 = dble(lon1(n))
                lat1_r8 = dble(lat1(n))
                call lambert(dble(stdlat1(n)),dble(stdlat2(n)),dble(cen_lat(n)),dble(cen_lon(n)), &
                             lon1_r8,lat1_r8,x1,y1, 1)
                do j=lbound(lonPtr,2),ubound(lonPtr,2)
                do i=lbound(lonPtr,1),ubound(lonPtr,1)
                  x = x1 + dx(n) * (i-1)
                  y = y1 + dy(n) * (j-1)
                  call lambert(dble(stdlat1(n)),dble(stdlat2(n)),dble(cen_lat(n)),dble(cen_lon(n)), &
                               geo_lon,geo_lat,x,y,-1)
                  if (geo_lon <0.0) geo_lon = geo_lon + 360.0
                  lonPtr(i,j) = geo_lon
                  latPtr(i,j) = geo_lat
                enddo
                enddo
            endif

          else

            write(0,*)"wrt_initialize_p1: Unknown output_grid ", trim(output_grid(n))
            call ESMF_LogWrite("wrt_initialize_p1: Unknown output_grid "//trim(output_grid(n)),ESMF_LOGMSG_ERROR,rc=RC)
            call ESMF_Finalize(endflag=ESMF_END_ABORT)

          endif

          wrt_int_state%out_grid_info(n)%i_start = lbound(lonPtr,1)
          wrt_int_state%out_grid_info(n)%i_end   = ubound(lonPtr,1)
          wrt_int_state%out_grid_info(n)%j_start = lbound(latPtr,2)
          wrt_int_state%out_grid_info(n)%j_end   = ubound(latPtr,2)

          allocate( wrt_int_state%out_grid_info(n)%i_start_wrtgrp(wrt_int_state%petcount) )
          allocate( wrt_int_state%out_grid_info(n)%i_end_wrtgrp  (wrt_int_state%petcount) )
          allocate( wrt_int_state%out_grid_info(n)%j_start_wrtgrp(wrt_int_state%petcount) )
          allocate( wrt_int_state%out_grid_info(n)%j_end_wrtgrp  (wrt_int_state%petcount) )

          call mpi_allgather(wrt_int_state%out_grid_info(n)%i_start,        1, MPI_INTEGER,    &
                             wrt_int_state%out_grid_info(n)%i_start_wrtgrp, 1, MPI_INTEGER, wrt_mpi_comm, rc)
          call mpi_allgather(wrt_int_state%out_grid_info(n)%i_end,          1, MPI_INTEGER,    &
                             wrt_int_state%out_grid_info(n)%i_end_wrtgrp,   1, MPI_INTEGER, wrt_mpi_comm, rc)
          call mpi_allgather(wrt_int_state%out_grid_info(n)%j_start,        1, MPI_INTEGER,    &
                             wrt_int_state%out_grid_info(n)%j_start_wrtgrp, 1, MPI_INTEGER, wrt_mpi_comm, rc)
          call mpi_allgather(wrt_int_state%out_grid_info(n)%j_end,          1, MPI_INTEGER,    &
                             wrt_int_state%out_grid_info(n)%j_end_wrtgrp,   1, MPI_INTEGER, wrt_mpi_comm, rc)

          allocate( wrt_int_state%out_grid_info(n)%lonPtr(wrt_int_state%out_grid_info(n)%i_start:wrt_int_state%out_grid_info(n)%i_end, &
                                                          wrt_int_state%out_grid_info(n)%j_start:wrt_int_state%out_grid_info(n)%j_end) )
          allocate( wrt_int_state%out_grid_info(n)%latPtr(wrt_int_state%out_grid_info(n)%i_start:wrt_int_state%out_grid_info(n)%i_end, &
                                                          wrt_int_state%out_grid_info(n)%j_start:wrt_int_state%out_grid_info(n)%j_end) )

          if ( trim(output_grid(n)) /= 'regional_latlon_moving' .and. trim(output_grid(n)) /= 'rotated_latlon_moving' ) then
            do j=wrt_int_state%out_grid_info(n)%j_start, wrt_int_state%out_grid_info(n)%j_end
            do i=wrt_int_state%out_grid_info(n)%i_start, wrt_int_state%out_grid_info(n)%i_end
                wrt_int_state%out_grid_info(n)%latPtr(i,j) = latPtr(i,j)
                wrt_int_state%out_grid_info(n)%lonPtr(i,j) = lonPtr(i,j)
            enddo
            enddo
          endif

          wrt_int_state%out_grid_info(n)%im = imo(n)
          wrt_int_state%out_grid_info(n)%jm = jmo(n)

        end if ! non 'cubed_sphere_grid'
      end do !  n = 1, ngrids
!
!-----------------------------------------------------------------------
!***  get write grid component initial time from clock
!-----------------------------------------------------------------------
!
      call ESMF_ClockGet(clock    =CLOCK                                &  !<-- The ESMF Clock
                        ,startTime=wrt_int_state%IO_BASETIME            &  !<-- The Clock's starting time
                        ,rc       =RC)

      call ESMF_TimeGet(time=wrt_int_state%IO_BASETIME,yy=idate(1),mm=idate(2),dd=idate(3), &
                                                        h=idate(4), m=idate(5), s=idate(6),rc=rc)
!     if (lprnt) write(0,*) 'in wrt initial, io_baseline time=',idate,'rc=',rc
      idate(7) = 1
      start_time = idate
      wrt_int_state%idate = idate
      wrt_int_state%fdate = idate
! update IO-BASETIME and idate on write grid comp when IAU is enabled
      if (iau_offset > 0) then
        call ESMF_TimeIntervalSet(IAU_offsetTI, h=iau_offset, rc=rc)
        wrt_int_state%IO_BASETIME = wrt_int_state%IO_BASETIME + IAU_offsetTI
        call ESMF_TimeGet(time=wrt_int_state%IO_BASETIME,yy=idate(1),mm=idate(2),dd=idate(3), &
                                                          h=idate(4), m=idate(5), s=idate(6),rc=rc)
        wrt_int_state%idate = idate
        change_wrtidate = .true.
        if (lprnt) print *,'in wrt initial, with iau, io_baseline time=',idate,'rc=',rc
      endif
!
!--- Look at the incoming FieldBundles in the imp_state_write, and mirror them as 'output_' bundles
!
      call ESMF_StateGet(imp_state_write, itemCount=FBCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      ! if (lprnt) write(0,*)'wrt_initialize_p1: FBCount=',FBCount, ' from imp_state_write'

      allocate(fcstItemNameList(FBCount), fcstItemTypeList(FBCount))
      allocate(outfilename(2000,FBCount))
      outfilename = ''

      call ESMF_StateGet(imp_state_write, itemNameList=fcstItemNameList, &
                         itemTypeList=fcstItemTypeList,                  &
                        !itemorderflag=ESMF_ITEMORDER_ADDORDER,          &
                         rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!loop over all items in the imp_state_write and collect all FieldBundles
      do i=1, FBCount

        if (fcstItemTypeList(i) == ESMF_STATEITEM_FIELDBUNDLE) then

          call ESMF_StateGet(imp_state_write, itemName=fcstItemNameList(i), &
                             fieldbundle=fcstFB, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_AttributeGet(fcstFB, convention="NetCDF", purpose="FV3", &
                                 name="grid_id", value=grid_id, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_AttributeGet(fcstFB, convention="NetCDF", purpose="FV3-nooutput", &
                                 name="frestart", valueList=frestart, isPresent=isPresent, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          if (isPresent) then
            ! if (lprnt) write(0,*)'wrt_initialize_p1: frestart(1:10) = ',frestart(1:10)
            call ESMF_AttributeRemove(fcstFB, convention="NetCDF", purpose="FV3-nooutput", name="frestart", rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          endif


!---  get grid dim count
          ! call ESMF_GridGet(wrtGrid(grid_id), dimCount=gridDimCount, rc=rc)
          ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! create a mirrored 'output_' FieldBundle and add it to importState
          fieldbundle = ESMF_FieldBundleCreate(name="output_"//trim(fcstItemNameList(i)), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_StateAdd(imp_state_write, (/fieldbundle/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! copy the fcstFB Attributes to the 'output_' FieldBundle
          call ESMF_AttributeCopy(fcstFB, fieldbundle, attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! grids in fcstFB for which 'is_moving' is .true. must provide a first level mirror for the Redist() target
          if (is_moving(grid_id)) then

! create a mirrored 'mirror_' FieldBundle and add it to importState
            mirrorFB = ESMF_FieldBundleCreate(name="mirror_"//trim(fcstItemNameList(i)), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_StateAdd(imp_state_write, (/mirrorFB/), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ! copy the fcstFB Attributes to the 'mirror_' FieldBundle
            call ESMF_AttributeCopy(fcstFB, mirrorFB, attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          endif

! deal with all of the Fields inside this fcstFB
          call ESMF_FieldBundleGet(fcstFB, fieldCount=fieldCount, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (fieldCount > 0) then

            call ESMF_FieldBundleGet(fcstFB, grid=fcstGrid, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            allocate(fcstField(fieldCount))
            call ESMF_FieldBundleGet(fcstFB, fieldList=fcstField,     &
                                     itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            if (fcstItemNameList(i)(1:18) == 'cubed_sphere_grid_') then

              if (create_wrtGrid_cubed_sphere) then
                ! create a grid from fcstGrid on forecast grid comp, by rebalancing distgrid to the local PETs
                ! access the acceptor DistGrid
                call ESMF_GridGet(fcstGrid, distgrid=acceptorDG, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                ! rebalance the acceptor DistGrid across the local PETs
                newAcceptorDG = ESMF_DistGridCreate(acceptorDG, balanceflag=.true., rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                wrtGrid_cubed_sphere = ESMF_GridCreate(fcstGrid, newAcceptorDG, copyAttributes=.true., rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

                create_wrtGrid_cubed_sphere = .false.
              end if

              actualWrtGrid = wrtGrid_cubed_sphere
              call ESMF_AttributeSet(fieldbundle, convention="NetCDF", purpose="FV3-nooutput", name="output_grid", value="cubed_sphere_grid", rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            else if (fcstItemNameList(i)(1:8) == 'restart_') then
              ! If this is a 'restart' bundle the actual grid that the output field ('field_work' below) is created on
              ! must be the same grid as forecast grid, not the output grid for this grid_id (wrtGrid(grid_id)).
              ! For 'cubed_sphere_grid' these are the same, but for all other output grids (like Lambert) they are not.

              ! create a grid from fcstGrid on forecast grid comp, by rebalancing distgrid to the local PETs
              ! access the acceptor DistGrid
              call ESMF_GridGet(fcstGrid, distgrid=acceptorDG, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              ! rebalance the acceptor DistGrid across the local PETs
              newAcceptorDG = ESMF_DistGridCreate(acceptorDG, balanceflag=.true., rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              actualWrtGrid = ESMF_GridCreate(fcstGrid, newAcceptorDG, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            else
              actualWrtGrid = wrtGrid(grid_id)
              call ESMF_AttributeSet(fieldbundle, convention="NetCDF", purpose="FV3-nooutput", name="output_grid", value=output_grid(grid_id), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            end if

            do j=1, fieldCount

              call ESMF_FieldGet(fcstField(j), typekind=typekind, dimCount=fieldDimCount, name=fieldName, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              call ESMF_GridGet(actualWrtGrid, dimCount=gridDimCount, rc=rc) ! use actualWrtGrid instead of wrtGrid(grid_id)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              allocate(gridToFieldMap(gridDimCount))
              allocate(ungriddedLBound(fieldDimCount-gridDimCount))
              allocate(ungriddedUBound(fieldDimCount-gridDimCount))

              call ESMF_FieldGet(fcstField(j), gridToFieldMap=gridToFieldMap,                      &
                                 ungriddedLBound=ungriddedLBound, ungriddedUBound=ungriddedUBound, &
                                 staggerloc=staggerloc, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!             if (lprnt) print *,'in wrt,fcstfld,fieldname=',                                         &
!                        trim(fieldname),'fieldDimCount=',fieldDimCount,'gridDimCount=',gridDimCount, &
!                        'gridToFieldMap=',gridToFieldMap,'ungriddedLBound=',ungriddedLBound,         &
!                        'ungriddedUBound=',ungriddedUBound,'rc=',rc

              ! create the output field on output grid
              field_work = ESMF_FieldCreate(actualWrtGrid, typekind, name=fieldName, & ! use actualWrtGrid instead of wrtGrid(grid_id)
                                            staggerloc=staggerloc,             &
                                            gridToFieldMap=gridToFieldMap,     &
                                            ungriddedLBound=ungriddedLBound,   &
                                            ungriddedUBound=ungriddedUBound, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              call ESMF_AttributeCopy(fcstField(j), field_work, attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              ! get output file name
              call ESMF_AttributeGet(fcstField(j), convention="NetCDF", purpose="FV3", &
                                     name="output_file", value=outfile_name, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              call ESMF_LogWrite("bf fcstfield, get output_file "//trim(outfile_name)//" "//trim(fieldName),ESMF_LOGMSG_INFO,rc=RC)
              if (trim(outfile_name) /= '') then
                outfilename(j,i) = trim(outfile_name)
              endif
              call ESMF_LogWrite("af fcstfield, get output_file",ESMF_LOGMSG_INFO,rc=RC)

              ! if (lprnt) print *,' i=',i,' j=',j,' outfilename=',trim(outfilename(j,i))

              ! add the output field to the 'output_' FieldBundle
              call ESMF_FieldBundleAdd(fieldbundle, (/field_work/), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              ! deal with grids for which 'is_moving' is .true.
              if (is_moving(grid_id)) then
                ! create an empty field that will serve as acceptor for GridTransfer of fcstGrid
                field_work = ESMF_FieldEmptyCreate(name=fieldName, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

                ! use attributes to carry information for later FieldEmptyComplete()
                call ESMF_InfoGetFromHost(field_work, info=info, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                tk = typekind ! convert TypeKind_Flag to integer
                call ESMF_InfoSet(info, key="typekind", value=tk, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                sloc = staggerloc ! convert StaggerLoc_Flag to integer
                call ESMF_InfoSet(info, key="staggerloc", value=sloc, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                call ESMF_InfoSet(info, key="gridToFieldMap", values=gridToFieldMap, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                call ESMF_InfoSet(info, key="ungriddedLBound", values=ungriddedLBound, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                call ESMF_InfoSet(info, key="ungriddedUBound", values=ungriddedUBound, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

                ! add to 'mirror_' FieldBundle
                call ESMF_FieldBundleAdd(mirrorFB, (/field_work/), rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              endif

              ! local garbage collection
              deallocate(gridToFieldMap, ungriddedLBound, ungriddedUBound)
            enddo

            call ESMF_AttributeCopy(fcstGrid, actualWrtGrid   , &
                                    attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            deallocate(fcstField)

          endif !if (fieldCount > 0) then

        else  ! anything but a FieldBundle in the state is unexpected here
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
                                msg="Only FieldBundles supported in fcstState.", line=__LINE__, file=__FILE__)
          return
        endif

      enddo !FBCount

      !loop over all items in the imp_state_write and count output FieldBundles
      call get_outfile(FBCount, outfilename, FBlist_outfilename, noutfile)
      wrt_int_state%FBCount = noutfile

      !create output field bundles
      allocate(wrt_int_state%wrtFB(wrt_int_state%FBCount))
      ! if (lprnt) write(0,*)'wrt_initialize_p1: allocated ',wrt_int_state%FBCount, ' wrt_int_state%wrtFB'

      do i=1, wrt_int_state%FBCount

        wrt_int_state%wrtFB(i) = ESMF_FieldBundleCreate(name=trim(FBlist_outfilename(i)), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        ! if (lprnt) write(0,*)'wrt_initialize_p1: created wrtFB ',i, ' with name ', trim(FBlist_outfilename(i))

        ! if (lprnt) write(0,*)'wrt_initialize_p1: loop over ', FBCount, ' forecast bundles'
        do n=1, FBCount

          call ESMF_StateGet(imp_state_write, itemName="output_"//trim(fcstItemNameList(n)), &
                             fieldbundle=fcstFB, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ! if (lprnt) write(0,*)'wrt_initialize_p1: got forecast bundle ', "output_"//trim(fcstItemNameList(n))
          ! if (lprnt) write(0,*)'wrt_initialize_p1: is ', trim(fcstItemNameList(n)), ' == ', trim(FBlist_outfilename(i))

          if (trim_regridmethod_suffix(fcstItemNameList(n)) == trim_regridmethod_suffix(FBlist_outfilename(i))) then

            ! copy the fcstfield bundle Attributes to the output field bundle
            ! if (lprnt) write(0,*)'wrt_initialize_p1: copy atts/fields from ', "output_"//trim(fcstItemNameList(n)), ' to ', trim(FBlist_outfilename(i))
            call ESMF_AttributeCopy(fcstFB,  wrt_int_state%wrtFB(i), &
                                    attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_AttributeGet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                   name="grid_id", value=grid_id, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            ! if (lprnt) write(0,*)'wrt_initialize_p1: got grid_id for wrtFB ', i, ' grid_id =', grid_id, trim(output_grid(grid_id))

            call ESMF_FieldBundleGet(fcstFB, fieldCount=fieldCount, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            allocate(fcstField(fieldCount),fieldnamelist(fieldCount))
            call ESMF_FieldBundleGet(fcstFB, fieldList=fcstField, fieldNameList=fieldnamelist,   &
                                     itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            do j=1, fieldCount

              call ESMF_AttributeGet(fcstField(j),convention="NetCDF", purpose="FV3", &
                                     name='output_file',value=outfile_name, rc=rc)

              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              if( trim(outfile_name) == trim(FBlist_outfilename(i))) then
                call ESMF_FieldBundleAdd(wrt_int_state%wrtFB(i), (/fcstField(j)/), rc=rc)

                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              endif

            enddo ! fieldCount
            deallocate(fcstField, fieldnamelist)

          endif ! index(trim(fcstItemNameList(n)),trim(FBlist_outfilename(i)))

        enddo ! FBCount

        ! add output grid related attributes, only for history files(bundles), skip restart
        if (FBlist_outfilename(i)(1:8) /= 'restart_') then

            call ESMF_AttributeGet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3-nooutput", &
                                   name="output_grid", value=output_grid_name, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_AttributeAdd(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                   attrList=(/"source","grid  "/), rc=rc)
            call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                   name="source", value="FV3GFS", rc=rc)

            if (trim(output_grid_name) == 'cubed_sphere_grid') then

              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="grid", value="cubed_sphere", rc=rc)

            else if (trim(output_grid_name) == 'gaussian_grid') then

              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="grid", value="gaussian", rc=rc)
              call ESMF_AttributeAdd(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     attrList=(/"im","jm"/), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="im", value=imo(grid_id), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="jm", value=jmo(grid_id), rc=rc)

            else if (trim(output_grid_name) == 'regional_latlon'        &
                .or. trim(output_grid_name) == 'regional_latlon_moving' &
                .or. trim(output_grid_name) == 'global_latlon') then

              ! for 'regional_latlon_moving' lon1/2 and lat1/2 will be overwritten in run phase
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="grid", value="latlon", rc=rc)
              call ESMF_AttributeAdd(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     attrList=(/"lon1","lat1","lon2","lat2","dlon","dlat"/), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="dlon", value=dlon(grid_id), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="dlat", value=dlat(grid_id), rc=rc)
              if (trim(output_grid_name) /= 'regional_latlon_moving') then
                call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                       name="lon1", value=lon1(grid_id), rc=rc)
                call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                       name="lat1", value=lat1(grid_id), rc=rc)
                call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                       name="lon2", value=lon2(grid_id), rc=rc)
                call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                       name="lat2", value=lat2(grid_id), rc=rc)
              endif
            else if (trim(output_grid_name) == 'rotated_latlon' &
                .or. trim(output_grid_name) == 'rotated_latlon_moving') then

              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="grid", value="rotated_latlon", rc=rc)
              call ESMF_AttributeAdd(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     attrList=(/"cen_lon",&
                                                "cen_lat",&
                                                "lon1   ",&
                                                "lat1   ",&
                                                "lon2   ",&
                                                "lat2   ",&
                                                "dlon   ",&
                                                "dlat   "/), rc=rc)
              ! for 'rotated_latlon_moving' cen_lon and cen_lat will be overwritten in run phase
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="cen_lon", value=cen_lon(grid_id), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="cen_lat", value=cen_lat(grid_id), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="dlon", value=dlon(grid_id), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="dlat", value=dlat(grid_id), rc=rc)
              if (trim(output_grid_name) /= 'rotated_latlon_moving') then
                call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                       name="lon1", value=lon1(grid_id), rc=rc)
                call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                       name="lat1", value=lat1(grid_id), rc=rc)
                call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                       name="lon2", value=lon2(grid_id), rc=rc)
                call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lat2", value=lat2(grid_id), rc=rc)
              endif
            else if (trim(output_grid_name) == 'lambert_conformal') then

              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="grid", value="lambert_conformal", rc=rc)
              call ESMF_AttributeAdd(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     attrList=(/"cen_lon",&
                                                "cen_lat",&
                                                "stdlat1",&
                                                "stdlat2",&
                                                "nx     ",&
                                                "ny     ",&
                                                "lon1   ",&
                                                "lat1   ",&
                                                "dx     ",&
                                                "dy     "/), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="cen_lon", value=cen_lon(grid_id), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="cen_lat", value=cen_lat(grid_id), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="stdlat1", value=stdlat1(grid_id), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="stdlat2", value=stdlat2(grid_id), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="nx", value=imo(grid_id), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="ny", value=jmo(grid_id), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lat1", value=lat1(grid_id), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="lon1", value=lon1(grid_id), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="dx", value=dx(grid_id), rc=rc)
              call ESMF_AttributeSet(wrt_int_state%wrtFB(i), convention="NetCDF", purpose="FV3", &
                                     name="dy", value=dy(grid_id), rc=rc)

            end if
        end if

      enddo ! end wrt_int_state%FBCount
!
! add time Attribute
! look at the importState attributes and copy those starting with "time"

      call ESMF_AttributeGet(imp_state_write, convention="NetCDF", purpose="FV3", &
                             attnestflag=ESMF_ATTNEST_OFF, count=attCount, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_LogWrite("Write component AttributeGet, attCount ", ESMF_LOGMSG_INFO, rc=rc)

! prepare the lists needed to transfer attributes

      allocate(attNameList(attCount), attNameList2(attCount))
      allocate(typekindList(attCount))

! loop over all the attributes on importState within AttPack
      j = 1
      k = 1
      do i=1, attCount
        call ESMF_AttributeGet(imp_state_write, convention="NetCDF", purpose="FV3",          &
                               attnestflag=ESMF_ATTNEST_OFF, attributeIndex=i, name=attName, &
                               typekind=typekind, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! test for name starting with "time"
        if (index(trim(attName), "time") == 1) then

! add this attribute to the list of transfers
          attNameList(j)  = attName
          typekindList(j) = typekind
          j = j + 1
          if (index(trim(attName), "time:") == 1) then
          ! store names of attributes starting with "time:" for later use
            attNameList2(k) = attName
            k = k+1
          endif
        endif
      enddo

    do n = 1, ngrids
    ! add the transfer attributes from importState to grid
    call ESMF_AttributeAdd(wrtGrid(n), convention="NetCDF", purpose="FV3", &
                           attrList=attNameList(1:j-1), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! add the transfer attributes from importState to special cubed_sphere grid
    if (n == 1 .and. top_parent_is_global .and. history_file_on_native_grid) then
      call ESMF_AttributeAdd(wrtGrid_cubed_sphere, convention="NetCDF", purpose="FV3", &
                             attrList=attNameList(1:j-1), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif

! loop over the added attributes, access the value (only scalar allowed),
! and set them on the grid
    do i=1, j-1
      if (typekindList(i) == ESMF_TYPEKIND_CHARACTER) then
        call ESMF_AttributeGet(imp_state_write,                    &
                               convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueS, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! save calendar_type (as integer) for use in 'coupler.res'
        if (index(trim(attNameList(i)),'time:calendar') > 0) then
          select case( fms_mpp_uppercase(trim(valueS)) )
          case( 'JULIAN' )
              calendar_type = JULIAN
          case( 'GREGORIAN' )
              calendar_type = GREGORIAN
          case( 'NOLEAP' )
              calendar_type = NOLEAP
          case( 'THIRTY_DAY' )
              calendar_type = THIRTY_DAY_MONTHS
          case( 'NO_CALENDAR' )
              calendar_type = NO_CALENDAR
          case default
              call fms_mpp_error ( FATAL, 'fcst_initialize: calendar must be one of '// &
                                      'JULIAN|GREGORIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.' )
          end select
        endif

! update the time:units when idate on write grid component is changed
        if (index(trim(attNameList(i)),'time:units') > 0) then
          if ( change_wrtidate ) then
            idx = index(trim(valueS),' since ')
            if(lprnt) print *,'in write grid comp, time:unit=',trim(valueS)
            write(newdate,'(I4.4,a,I2.2,a,I2.2,a,I2.2,a,I2.2,a,I2.2)') idate(1),'-',   &
              idate(2),'-',idate(3),' ',idate(4),':',idate(5),':',idate(6)
            valueS = valueS(1:idx+6)//newdate
            if(lprnt) print *,'in write grid comp, new time:unit=',trim(valueS)
          endif
        endif
        call ESMF_AttributeSet(wrtGrid(n), convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueS, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if (n == 1 .and. top_parent_is_global .and. history_file_on_native_grid) then
          call ESMF_AttributeSet(wrtGrid_cubed_sphere, convention="NetCDF", purpose="FV3", &
                                 name=trim(attNameList(i)), value=valueS, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        endif

      else if (typekindList(i) == ESMF_TYPEKIND_I4) then
        call ESMF_AttributeGet(imp_state_write,                    &
                               convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueI4, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_AttributeSet(wrtGrid(n), convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueI4, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      else if (typekindList(i) == ESMF_TYPEKIND_R4) then
        call ESMF_AttributeGet(imp_state_write,                    &
                               convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueR4, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_AttributeSet(wrtGrid(n), convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueR4, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      else if (typekindList(i) == ESMF_TYPEKIND_R8) then
        call ESMF_AttributeGet(imp_state_write,                    &
                               convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueR8, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_AttributeSet(wrtGrid(n), convention="NetCDF", purpose="FV3", &
                               name=trim(attNameList(i)), value=valueR8, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
    enddo

! Add special attribute that holds names of "time" related attributes
! for faster access during Run().
    call ESMF_AttributeAdd(wrtGrid(n), convention="NetCDF", purpose="FV3", &
                           attrList=(/"TimeAttributes"/), rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeSet(wrtGrid(n), convention="NetCDF", purpose="FV3", &
                           name="TimeAttributes", valueList=attNameList2(1:k-1), rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return


!
!*** create temporary field bundle for  axes information
! write the Grid coordinate arrays into the output files via temporary FB

      gridFB = ESMF_FieldBundleCreate(rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_GridGetCoord(wrtGrid(n), coordDim=1, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      field = ESMF_FieldCreate(wrtGrid(n), array, name="grid_xt", rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

!add attribute info
! long name
      call ESMF_AttributeAdd(field,convention="NetCDF",purpose="FV3",  &
                             attrList=(/'long_name'/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(field,convention="NetCDF",purpose="FV3",name='long_name', &
                             value="T-cell longitude", rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
! units
      call ESMF_AttributeAdd(field,convention="NetCDF",purpose="FV3",  &
                             attrList=(/'units'/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(field,convention="NetCDF",purpose="FV3",name='units', &
                             value="degrees_E", rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
! cartesian_axis
      call ESMF_AttributeAdd(field,convention="NetCDF",purpose="FV3",  &
                            attrList=(/'cartesian_axis'/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(field,convention="NetCDF",purpose="FV3",name='cartesian_axis', &
                             value="X", rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! add field to bundle
      call ESMF_FieldBundleAdd(gridFB, (/field/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
! get 2nd dimension
      call ESMF_GridGetCoord(wrtGrid(n), coordDim=2, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      field = ESMF_FieldCreate(wrtGrid(n), array, name="grid_yt", rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!add attribute info
! long name
      call ESMF_AttributeAdd(field,convention="NetCDF",purpose="FV3",  &
                             attrList=(/'long_name'/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(field,convention="NetCDF",purpose="FV3",name='long_name', &
                              value="T-cell latitude", rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
! units
      call ESMF_AttributeAdd(field,convention="NetCDF",purpose="FV3",  &
                             attrList=(/'units'/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(field,convention="NetCDF",purpose="FV3",name='units', &
                             value="degrees_N", rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! cartesian_axis
      call ESMF_AttributeAdd(field,convention="NetCDF",purpose="FV3",  &
                             attrList=(/'cartesian_axis'/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_AttributeSet(field,convention="NetCDF",purpose="FV3",name='cartesian_axis', &
                             value="Y", rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_FieldBundleAdd(gridFB, (/field/), rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      end do ! n=1, ngrids

      deallocate(attNameList, attNameList2, typekindList)
!
!     write_init_tim = MPI_Wtime() - btim0
!
!-----------------------------------------------------------------------
!
      end subroutine wrt_initialize_p1
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine wrt_initialize_p2(wrt_comp, imp_state_write, exp_state_write, clock, rc)
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE WRITE GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      type(esmf_GridComp)               :: wrt_comp
      type(ESMF_State)                  :: imp_state_write, exp_state_write
      type(esmf_Clock)                  :: clock
      integer,intent(out)               :: rc
!
!***  LOCAL VARIABLES
      type(ESMF_Info)                         :: info
      logical, allocatable                    :: is_moving(:)
      type(ESMF_StateItem_Flag), allocatable  :: itemTypeList(:)
      character(len=ESMF_MAXSTR),allocatable  :: itemNameList(:)
      integer                                 :: i, j, bundleCount, fieldCount
      type(ESMF_FieldBundle)                  :: mirrorFB
      type(ESMF_Field), allocatable           :: fieldList(:)
      type(ESMF_Grid)                         :: grid
      integer                                 :: sloc
      type(ESMF_StaggerLoc)                   :: staggerloc
      type(ESMF_DistGrid)                     :: acceptorDG, newAcceptorDG
!
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      call ESMF_InfoGetFromHost(imp_state_write, info=info, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_InfoGetAlloc(info, key="is_moving", values=is_moving, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_StateGet(imp_state_write, itemCount=bundleCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      allocate(itemNameList(bundleCount), itemTypeList(bundleCount))

      call ESMF_StateGet(imp_state_write, itemNameList=itemNameList, &
                         itemTypeList=itemTypeList,                  &
                        !itemorderflag=ESMF_ITEMORDER_ADDORDER,      &
                         rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      do i=1, bundleCount

        if (itemTypeList(i) == ESMF_STATEITEM_FIELDBUNDLE) then

          if (index(trim(itemNameList(i)), "mirror_")==1) then
            ! this is a 'mirror_' FieldBundle -> GridTransfer acceptor side
            call ESMF_StateGet(imp_state_write, itemName=trim(itemNameList(i)), fieldbundle=mirrorFB, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            ! access the grid that is passed in from the provider side
            call ESMF_FieldBundleGet(mirrorFB, grid=grid, fieldCount=fieldCount, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            ! access the acceptor DistGrid
            call ESMF_GridGet(grid, distgrid=acceptorDG, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            ! rebalance the acceptor DistGrid across the local PETs
            newAcceptorDG = ESMF_DistGridCreate(acceptorDG, balanceflag=.true., rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            ! create a new Grid on the rebalanced DistGrid
            grid = ESMF_GridCreate(newAcceptorDG, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            ! point all of the acceptor fields to the new acceptor Grid
            allocate(fieldList(fieldCount))
            call ESMF_FieldBundleGet(mirrorFB, fieldList=fieldList, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            do j=1, fieldCount
              ! first access information stored on the field needed for completion
              call ESMF_InfoGetFromHost(fieldList(j), info=info, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              call ESMF_InfoGet(info, key="staggerloc", value=sloc, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              staggerloc = sloc  ! convert integer into StaggerLoc_Flag
              call ESMF_FieldEmptySet(fieldList(j), grid=grid, staggerloc=staggerloc, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            enddo
            ! clean-up
            deallocate(fieldList)
          endif

        else  ! anything but a FieldBundle in the state is unexpected here
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, msg="Only FieldBundles supported in fcstState.", line=__LINE__, file=__FILE__)
          return
        endif

      enddo

!-----------------------------------------------------------------------
!
      end subroutine wrt_initialize_p2
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine wrt_initialize_p3(wrt_comp, imp_state_write, exp_state_write, clock, rc)
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE WRITE GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      type(esmf_GridComp)               :: wrt_comp
      type(ESMF_State)                  :: imp_state_write, exp_state_write
      type(esmf_Clock)                  :: clock
      integer,intent(out)               :: rc
!***  LOCAL VARIABLES
      type(ESMF_Info)                         :: info
      logical, allocatable                    :: is_moving(:)
      type(ESMF_StateItem_Flag), allocatable  :: itemTypeList(:)
      character(len=ESMF_MAXSTR),allocatable  :: itemNameList(:)
      integer                                 :: i, j, bundleCount, fieldCount
      type(ESMF_FieldBundle)                  :: mirrorFB
      type(ESMF_Field), allocatable           :: fieldList(:)
      type(ESMF_TypeKind_Flag)                :: typekind
      integer                                 :: tk
      integer,                   allocatable  :: gridToFieldMap(:)
      integer,                   allocatable  :: ungriddedLBound(:)
      integer,                   allocatable  :: ungriddedUBound(:)

!
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      call ESMF_InfoGetFromHost(imp_state_write, info=info, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_InfoGetAlloc(info, key="is_moving", values=is_moving, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_StateGet(imp_state_write, itemCount=bundleCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      allocate(itemNameList(bundleCount), itemTypeList(bundleCount))

      call ESMF_StateGet(imp_state_write, itemNameList=itemNameList, &
                         itemTypeList=itemTypeList,                  &
                        !itemorderflag=ESMF_ITEMORDER_ADDORDER,      &
                         rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      do i=1, bundleCount

        if (itemTypeList(i) == ESMF_STATEITEM_FIELDBUNDLE) then

          if (index(trim(itemNameList(i)), "mirror_")==1) then
            ! this is a 'mirror_' FieldBundle -> GridTransfer acceptor side
            call ESMF_StateGet(imp_state_write, itemName=trim(itemNameList(i)), fieldbundle=mirrorFB, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            ! finish creating all the mirror Fields
            call ESMF_FieldBundleGet(mirrorFB, fieldCount=fieldCount, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            allocate(fieldList(fieldCount))
            call ESMF_FieldBundleGet(mirrorFB, fieldList=fieldList, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            do j=1, fieldCount
              ! first access information stored on the field needed for completion
              call ESMF_InfoGetFromHost(fieldList(j), info=info, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              call ESMF_InfoGet(info, key="typekind", value=tk, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              typekind = tk  ! convert integer into TypeKind_Flag
              call ESMF_InfoGetAlloc(info, key="gridToFieldMap", values=gridToFieldMap, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              call ESMF_InfoGetAlloc(info, key="ungriddedLBound", values=ungriddedLBound, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              call ESMF_InfoGetAlloc(info, key="ungriddedUBound", values=ungriddedUBound, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              ! now complete the field creation
              call ESMF_FieldEmptyComplete(fieldList(j), typekind=typekind, gridToFieldMap=gridToFieldMap, &
                ungriddedLBound=ungriddedLBound, ungriddedUBound=ungriddedUBound, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            enddo
            ! clean-up
            deallocate(fieldList)
          endif

        else  ! anything but a FieldBundle in the state is unexpected here
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, msg="Only FieldBundles supported in fcstState.", line=__LINE__, file=__FILE__)
          return
        endif

      enddo

!-----------------------------------------------------------------------
!
      end subroutine wrt_initialize_p3
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine wrt_run(wrt_comp, imp_state_write, exp_state_write,clock,rc)
!
!-----------------------------------------------------------------------
!***  the run step for the write gridded component.
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp)            :: wrt_comp
      type(ESMF_State)               :: imp_state_write, exp_state_write
      type(ESMF_Clock)               :: clock
      integer,intent(out)            :: rc
!
!-----------------------------------------------------------------------
!***  local variables
!
      TYPE(ESMF_VM)                         :: VM
      type(ESMF_FieldBundle)                :: file_bundle, mirror_bundle
      type(ESMF_StateItem_Flag)             :: itemType
      type(ESMF_Time)                       :: currtime
      type(ESMF_TimeInterval)               :: io_currtimediff
      type(ESMF_Grid)                       :: fbgrid, wrtGrid
      type(ESMF_State),save                 :: stateGridFB
      type(optimizeT), save                 :: optimize(40)   ! FIXME
      type(ESMF_GridComp), save, allocatable   :: compsGridFB(:)
      type(ESMF_RouteHandle)                :: rh
      type(ESMF_RegridMethod_Flag)          :: regridmethod
      integer                               :: srcTermProcessing
!
      type(write_wrap)                      :: wrap
      type(wrt_internal_state),pointer      :: wrt_int_state
!
      integer                               :: i,j,n,mype,nolog, grid_id, localPet
!
      integer                               :: nf_hours,nf_seconds,nf_minutes
      integer                               :: fcst_seconds
      real(ESMF_KIND_R8)                    :: nfhour
!
      integer                               :: nbdl, cdate(6), ndig, nnnn
      integer                               :: step=1
      integer                               :: out_phase
!
      logical                               :: opened
      logical                               :: lmask_fields
!
      character(esmf_maxstr)                :: filename,compname,wrtFBName,traceString
      character(40)                         :: cfhour, cform
      character(20)                         :: time_iso
      character(15)                         :: time_restart
      character(15)                         :: tile_id
!
      type(ESMF_Grid)                       :: grid
      type(ESMF_Info)                       :: info
      real(ESMF_KIND_R8), allocatable       :: values(:)
      character(160)                        :: msgString
      type(ESMF_Field), allocatable         :: fieldList(:)
      type(ESMF_Array)                      :: coordArray(2)
      type(ESMF_DistGrid)                   :: coordDG
      type(ESMF_DELayout)                   :: coordDL
      integer                               :: fieldCount, deCount, rootPet
      integer                               :: minIndexPTile(2,1), maxIndexPTile(2,1), centerIndex(2)
      integer, allocatable                  :: minIndexPDe(:,:), maxIndexPDe(:,:), petMap(:)
      real(ESMF_KIND_R8), pointer           :: farrayPtr(:,:)
      real(ESMF_KIND_R8)                    :: centerCoord(2)

      integer                               :: ii, jj
      real(ESMF_KIND_R8), pointer           :: lonPtr(:,:), latPtr(:,:)
      real(ESMF_KIND_R8)                    :: rot_lon, rot_lat
      real(ESMF_KIND_R8)                    :: geo_lon, geo_lat
      real(ESMF_KIND_R8), parameter         :: rtod=180.0/pi

      real(kind=8)  :: MPI_Wtime
      real(kind=8)  :: tbeg
      real(kind=8)  :: wbeg,wend

      logical                               :: use_parallel_netcdf
      real, allocatable                     :: output_fh(:)
      logical                               :: is_restart_bundle, restart_written
      integer                               :: tileCount
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      tbeg = MPI_Wtime()
      rc   = esmf_success
!
!-----------------------------------------------------------------------
!***  get the current write grid comp name, and internal state
!
      call ESMF_GridCompGet(wrt_comp, name=compname, localPet=localPet, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
! Provide log message indicating which wrtComp is active
      call ESMF_LogWrite("Write component activated: "//trim(compname), &
                          ESMF_LOGMSG_INFO, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

! access the internal state
      call ESMF_GridCompGetInternalState(wrt_Comp, wrap, rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      wrt_int_state => wrap%write_int_state

      call ESMF_VMGetCurrent(VM,rc=RC)

      mype = wrt_int_state%mype
!    print *,'in wrt run, mype=',mype,'lead_write_task=',lead_write_task

      call ESMF_InfoGetFromHost(imp_state_write, info=info, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      call ESMF_InfoGetAlloc(info, key="output_fh", values=output_fh, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!
!-----------------------------------------------------------------------
!*** get current time and elapsed forecast time

      call ESMF_ClockGet(clock=CLOCK, currTime=currTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_TimeGet(time=currTime,yy=cdate(1),mm=cdate(2),dd=cdate(3), &
                                       h=cdate(4), m=cdate(5), s=cdate(6),rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      wrt_int_state%fdate(7) = 1
      wrt_int_state%fdate(1:6) = cdate(1:6)
      write(time_iso,'(I4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2,"Z")') cdate(1:6)

      io_currtimediff = currtime - wrt_int_state%IO_BASETIME

      call ESMF_TimeIntervalGet(timeinterval=io_currtimediff &
                               ,h_r8=nfhour,h=nf_hours,m=nf_minutes,s=nf_seconds,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (nf_hours < 0) return

      if (lflname_fulltime) then
        ndig = max(log10(nf_hours+0.5)+1., 3.)
        write(cform, '("(I",I1,".",I1,",A1,I2.2,A1,I2.2)")') ndig, ndig
        write(cfhour, cform) nf_hours,'-',nf_minutes,'-',nf_seconds
      else
        ndig = max(log10(nf_hours+0.5)+1., 3.)
        write(cform, '("(I",I1,".",I1,")")') ndig, ndig
        write(cfhour, cform) nf_hours
      endif
!
       if(lprnt) print *,'in wrt run, nfhour=',nfhour,' cfhour=',trim(cfhour)
!
!-----------------------------------------------------------------------
!*** loop on the "output_" FieldBundles, i.e. files that need to write out
!-----------------------------------------------------------------------

      do i=1, FBCount
        call ESMF_StateGet(imp_state_write, itemName="output_"//trim(fcstItemNameList(i)), &
                           fieldbundle=file_bundle, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        ! see whether a "mirror_" FieldBundle exists, i.e. dealing with moving domain that needs updated Regrid() here.
        call ESMF_StateGet(imp_state_write, itemName="mirror_"//trim(fcstItemNameList(i)), &
                           itemType=itemType, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if (itemType == ESMF_STATEITEM_FIELDBUNDLE) then
          ! Regrid() for a moving domain
          call ESMF_LogWrite("Regrid() for moving domain: mirror_"//trim(fcstItemNameList(i)), ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          call ESMF_StateGet(imp_state_write, itemName="mirror_"//trim(fcstItemNameList(i)), &
                             fieldbundle=mirror_bundle, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ! if (fcstItemNameList(i)(1:8) == "restart_" .or. fcstItemNameList(i)(1:18) == 'cubed_sphere_grid_') then
          if (fcstItemNameList(i)(1:8) == "restart_") then
            ! restart output forecast bundles, use Redist instead of Regrid

            call ESMF_FieldBundleRedistStore(mirror_bundle, file_bundle,                &
                                             routehandle=rh, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          else ! not restart bundle

            ! Find the centerCoord of the moving domain
            call ESMF_FieldBundleGet(mirror_bundle, fieldCount=fieldCount, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            allocate(fieldList(fieldCount))
            call ESMF_FieldBundleGet(mirror_bundle, fieldList=fieldList, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            call ESMF_FieldGet(fieldList(1), grid=grid, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            deallocate(fieldList)

            call ESMF_GridGetCoord(grid, coordDim=1, array=coordArray(1), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            call ESMF_GridGetCoord(grid, coordDim=2, array=coordArray(2), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            call ESMF_ArrayGet(coordArray(1), distgrid=coordDG, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            call ESMF_DistGridGet(coordDG, deCount=deCount, minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, &
              delayout=coordDL, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            allocate(petMap(deCount),minIndexPDe(2,deCount), maxIndexPDe(2,deCount))
            call ESMF_DELayoutGet(coordDL, petMap=petMap, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            call ESMF_DistGridGet(coordDG, minIndexPDe=minIndexPDe, maxIndexPDe=maxIndexPDe, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            centerIndex(1) = (maxIndexPTile(1,1)-minIndexPTile(1,1)+1)/2
            centerIndex(2) = (maxIndexPTile(2,1)-minIndexPTile(2,1)+1)/2

  !          write(msgString,*) "Determined centerIndex: ", centerIndex
  !          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_DEBUG, rc=rc)
  !          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            do n=1, deCount
              if (minIndexPDe(1,n)<=centerIndex(1) .and. centerIndex(1)<=maxIndexPDe(1,n) .and. &
                  minIndexPDe(2,n)<=centerIndex(2) .and. centerIndex(2)<=maxIndexPDe(2,n)) then
                ! found the DE that holds the center coordinate
                rootPet = petMap(n)
                if (localPet == rootPet) then
                  ! center DE is on local PET -> fill centerCoord locally
                  call ESMF_ArrayGet(coordArray(1), farrayPtr=farrayPtr, rc=rc)
                  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                  centerCoord(1) = farrayPtr(centerIndex(1)-minIndexPDe(1,n)+1,centerIndex(2)-minIndexPDe(2,n)+1)
                  call ESMF_ArrayGet(coordArray(2), farrayPtr=farrayPtr, rc=rc)
                  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                  centerCoord(2) = farrayPtr(centerIndex(1)-minIndexPDe(1,n)+1,centerIndex(2)-minIndexPDe(2,n)+1)
  !                write(msgString,*) "Found centerCoord: ", centerCoord
  !                call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_DEBUG, rc=rc)
  !                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                endif
                exit
              endif
            enddo

            deallocate(petMap,minIndexPDe,maxIndexPDe)

            call ESMF_VMBroadcast(vm, centerCoord, count=2, rootPet=rootPet, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            write(msgString,*) "All PETs know centerCoord in radians: ", centerCoord
            call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_DEBUG, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ! determine regridmethod
            if (index(fcstItemNameList(i),"_bilinear") >0 )  then
              traceString = "-bilinear"
              regridmethod = ESMF_REGRIDMETHOD_BILINEAR
            else if (index(fcstItemNameList(i),"_patch") >0)  then
              traceString = "-patch"
              regridmethod = ESMF_REGRIDMETHOD_PATCH
            else if (index(fcstItemNameList(i),"_nearest_stod") >0) then
              traceString = "-nearest_stod"
              regridmethod = ESMF_REGRIDMETHOD_NEAREST_STOD
            else if (index(fcstItemNameList(i),"_nearest_dtos") >0) then
              traceString = "-nearest_dtos"
              regridmethod = ESMF_REGRIDMETHOD_NEAREST_DTOS
            else if (index(fcstItemNameList(i),"_conserve") >0) then
              traceString = "-conserve"
              regridmethod = ESMF_REGRIDMETHOD_CONSERVE
            else
              call ESMF_LogSetError(ESMF_RC_ARG_BAD,                          &
                                    msg="Unable to determine regrid method.", &
                                    line=__LINE__, file=__FILE__, rcToReturn=rc)
              return
            endif
            srcTermProcessing = 1 ! have this fixed for bit-for-bit reproducibility
            ! RegridStore()

            ! update output grid coordinates based of fcstgrid center lat/lon
            call ESMF_FieldBundleGet(file_bundle, grid=grid, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            call ESMF_GridGetCoord(grid, coordDim=1, farrayPtr=lonPtr, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            call ESMF_GridGetCoord(grid, coordDim=2, farrayPtr=latPtr, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            call ESMF_AttributeGet(mirror_bundle, convention="NetCDF", purpose="FV3", &
                                   name="grid_id", value=grid_id, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            if (trim(output_grid(grid_id)) == 'regional_latlon_moving' .or. &
                trim(output_grid(grid_id)) == 'rotated_latlon_moving') then
              n = grid_id
              cen_lon(n) = centerCoord(1)*rtod
              cen_lat(n) = centerCoord(2)*rtod
              if (cen_lon(n) > 180.0) cen_lon(n) = cen_lon(n) - 360.0
              cen_lon(n) = NINT(cen_lon(n)*1000.0)/1000.0
              cen_lat(n) = NINT(cen_lat(n)*1000.0)/1000.0
            endif

            if (trim(output_grid(grid_id)) == 'regional_latlon_moving') then
              lon1(n) = cen_lon(n) - 0.5 * (imo(n)-1) * dlon(n)
              lat1(n) = cen_lat(n) - 0.5 * (jmo(n)-1) * dlat(n)
              lon2(n) = cen_lon(n) + 0.5 * (imo(n)-1) * dlon(n)
              lat2(n) = cen_lat(n) + 0.5 * (jmo(n)-1) * dlat(n)
              do jj=lbound(lonPtr,2),ubound(lonPtr,2)
              do ii=lbound(lonPtr,1),ubound(lonPtr,1)
                lonPtr(ii,jj) = lon1(n) + (lon2(n)-lon1(n))/(imo(n)-1) * (ii-1)
                latPtr(ii,jj) = lat1(n) + (lat2(n)-lat1(n))/(jmo(n)-1) * (jj-1)
                wrt_int_state%out_grid_info(n)%latPtr(ii,jj) = latPtr(ii,jj)
                wrt_int_state%out_grid_info(n)%lonPtr(ii,jj) = lonPtr(ii,jj)
              enddo
              enddo
            else if (trim(output_grid(grid_id)) == 'rotated_latlon_moving') then
              lon1(n) = - 0.5 * (imo(n)-1) * dlon(n)
              lat1(n) = - 0.5 * (jmo(n)-1) * dlat(n)
              lon2(n) =   0.5 * (imo(n)-1) * dlon(n)
              lat2(n) =   0.5 * (jmo(n)-1) * dlat(n)
              do jj=lbound(lonPtr,2),ubound(lonPtr,2)
              do ii=lbound(lonPtr,1),ubound(lonPtr,1)
                rot_lon = lon1(n) + (lon2(n)-lon1(n))/(imo(n)-1) * (ii-1)
                rot_lat = lat1(n) + (lat2(n)-lat1(n))/(jmo(n)-1) * (jj-1)
                call rtll(rot_lon, rot_lat, geo_lon, geo_lat, dble(cen_lon(n)), dble(cen_lat(n)))
                if (geo_lon < 0.0) geo_lon = geo_lon + 360.0
                lonPtr(ii,jj) = geo_lon
                latPtr(ii,jj) = geo_lat
                wrt_int_state%out_grid_info(n)%latPtr(ii,jj) = latPtr(ii,jj)
                wrt_int_state%out_grid_info(n)%lonPtr(ii,jj) = lonPtr(ii,jj)
              enddo
              enddo
            endif

            call ESMF_TraceRegionEnter("ESMF_FieldBundleRegridStore()"//trim(traceString), rc=rc)
            call ESMF_FieldBundleRegridStore(mirror_bundle, file_bundle,                &
                                             regridMethod=regridmethod, routehandle=rh, &
                                             unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                                             srcTermProcessing=srcTermProcessing, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            call ESMF_TraceRegionExit("ESMF_FieldBundleRegridStore()"//trim(traceString), rc=rc)

          endif ! fieldbundle restart vs. not restart

          ! Regrid()
          call ESMF_TraceRegionEnter("ESMF_FieldBundleRegrid()"//trim(traceString), rc=rc)
          call ESMF_FieldBundleRegrid(mirror_bundle, file_bundle, &
                                      routehandle=rh, termorderflag=(/ESMF_TERMORDER_SRCSEQ/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          call ESMF_TraceRegionExit("ESMF_FieldBundleRegrid()"//trim(traceString), rc=rc)
          ! RegridRelease()
          call ESMF_FieldBundleRegridRelease(routehandle=rh, noGarbage=.true., rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          ! done
          call ESMF_LogWrite("Done Regrid() for moving domain: mirror_"//trim(fcstItemNameList(i)), ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        endif

        if (fcstItemNameList(i)(1:8) /= "restart_") then
          !recover fields from cartesian vector and sfc pressure
          call recover_fields(file_bundle,rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        end if

      enddo
!
!-----------------------------------------------------------------------
!*** do post
!-----------------------------------------------------------------------
      lmask_fields = .false.
      if( wrt_int_state%write_dopost ) then
#ifdef INLINE_POST
        wbeg = MPI_Wtime()
        do n=1,ngrids

          if (trim(output_grid(n)) /= 'cubed_sphere_grid') then

            if (trim(output_grid(n)) == 'regional_latlon' .or. &
                trim(output_grid(n)) == 'regional_latlon_moving' .or. &
                trim(output_grid(n)) == 'rotated_latlon' .or. &
                trim(output_grid(n)) == 'rotated_latlon_moving' .or. &
                trim(output_grid(n)) == 'lambert_conformal') then

                !mask fields according to sfc pressure, only history bundles
                do nbdl=1, wrt_int_state%FBCount
                  call ESMF_FieldBundleGet(wrt_int_state%wrtFB(nbdl), name=wrtFBName, rc=rc)
                  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) return

                  if (wrtFBName(1:8) == 'restart_') cycle
                  if (wrtFBName(1:18) == 'cubed_sphere_grid_') cycle

                  call mask_fields(wrt_int_state%wrtFB(nbdl),rc)
                  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                enddo
                lmask_fields = .true.
            endif

            call post_run_fv3(wrt_int_state, n, mype, wrt_mpi_comm, lead_write_task, &
                              itasks, jtasks, nf_hours, nf_minutes, nf_seconds)
          else
            rc = ESMF_RC_NOT_IMPL
            print *,'Inline post not available for cubed_sphere_grid'
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return
          endif
        enddo
        wend = MPI_Wtime()
        if (mype == lead_write_task) then
          !** write out inline post log file
          open(newunit=nolog,file='log.atm.inlinepost.f'//trim(cfhour),form='FORMATTED')
          write(nolog,"('completed: fv3atm')")
          write(nolog,"('forecast hour: ',f10.3)") nfhour
          write(nolog,"('valid time: ',6(i4,2x))") wrt_int_state%fdate(1:6)
          close(nolog)
        endif
        if (lprnt) then
          write(*,'(A,F10.5,A,I4.2,A,I2.2)')' actual    inline post time is ',wend-wbeg &
                     ,' at Fcst ',nf_hours,':',nf_minutes
        endif
#else
        rc = ESMF_RC_NOT_IMPL
        print *,'inline post not available on this machine'
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
#endif
      endif
!
!-----------------------------------------------------------------------
! ** now loop through output field bundle
!-----------------------------------------------------------------------

      call ESMF_TimeIntervalGet(timeinterval=io_currtimediff, s=fcst_seconds, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      ! fcst_seconds is number of seconds in io_currtimediff, which is time interval between currenttime and io_basetime.
      ! io_basetime has been adjusted by iau_offset in initialize phase.
      ! Since output_fh and frestart and NOT adjusted by iau_offset, in order to compare
      ! them with fcst_seconds, we must also adjust fcst_seconds by iau_offset
      if (iau_offset > 0) then
        fcst_seconds = fcst_seconds + iau_offset*3600
      endif

      if ( (wrt_int_state%output_history .and. ANY(nint(output_fh(:)*3600.0) == fcst_seconds)) .or. ANY(frestart(:) == fcst_seconds) ) then

        ! if (lprnt) write(0,*)'wrt_run: loop over wrt_int_state%FBCount ',wrt_int_state%FBCount, ' nfhour ',  nfhour, ' cdate ', cdate(1:6)
        two_phase_loop: do out_phase = 1, 2

          restart_written = .false.
          file_loop_all: do nbdl=1, wrt_int_state%FBCount

            call ESMF_FieldBundleGet(wrt_int_state%wrtFB(nbdl), name=wrtFBName, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) return

            is_restart_bundle = .false.
            if (wrtFBName(1:8) == 'restart_') then
              is_restart_bundle = .true.
              if (.not.(ANY(frestart(:) == fcst_seconds))) cycle
            else
              if (.not.(wrt_int_state%output_history .and. ANY(nint(output_fh(:)*3600.0) == fcst_seconds))) cycle
            endif

            if (out_phase == 1 .and. is_restart_bundle) cycle
            if (out_phase == 2 .and. .not.is_restart_bundle) cycle

            ! get grid_id
            call ESMF_AttributeGet(wrt_int_state%wrtFB(nbdl), convention="NetCDF", purpose="FV3", &
                                   name="grid_id", value=grid_id, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ! update lon1/2 and lat1/2 for regional_latlon_moving
            if (trim(output_grid(grid_id)) == 'regional_latlon_moving') then
              call ESMF_AttributeSet(wrt_int_state%wrtFB(nbdl), convention="NetCDF", purpose="FV3", &
                                     name="lon1", value=lon1(grid_id), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              call ESMF_AttributeSet(wrt_int_state%wrtFB(nbdl), convention="NetCDF", purpose="FV3", &
                                     name="lat1", value=lat1(grid_id), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              call ESMF_AttributeSet(wrt_int_state%wrtFB(nbdl), convention="NetCDF", purpose="FV3", &
                                     name="lon2", value=lon2(grid_id), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              call ESMF_AttributeSet(wrt_int_state%wrtFB(nbdl), convention="NetCDF", purpose="FV3", &
                                     name="lat2", value=lat2(grid_id), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            endif

            ! update cen_lon/cen_lat, lon1/2 and lat1/2  for rotated_latlon_moving
            if (trim(output_grid(grid_id)) == 'rotated_latlon_moving') then
              call ESMF_AttributeSet(wrt_int_state%wrtFB(nbdl), convention="NetCDF", purpose="FV3", &
                                     name="cen_lon", value=cen_lon(grid_id), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              call ESMF_AttributeSet(wrt_int_state%wrtFB(nbdl), convention="NetCDF", purpose="FV3", &
                                     name="cen_lat", value=cen_lat(grid_id), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              call ESMF_AttributeSet(wrt_int_state%wrtFB(nbdl), convention="NetCDF", purpose="FV3", &
                                     name="lon1", value=lon1(grid_id), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              call ESMF_AttributeSet(wrt_int_state%wrtFB(nbdl), convention="NetCDF", purpose="FV3", &
                                     name="lat1", value=lat1(grid_id), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              call ESMF_AttributeSet(wrt_int_state%wrtFB(nbdl), convention="NetCDF", purpose="FV3", &
                                     name="lon2", value=lon2(grid_id), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              call ESMF_AttributeSet(wrt_int_state%wrtFB(nbdl), convention="NetCDF", purpose="FV3", &
                                     name="lat2", value=lat2(grid_id), rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            endif

            if(step == 1) then
              file_bundle = wrt_int_state%wrtFB(nbdl)
            endif

            ! FIXME  map nbdl to [1:num_files], only used for output_file
            nnnn = mod(nbdl-1, num_files) + 1

            ! set default chunksizes for netcdf output
            ! (use MPI decomposition size).
            ! if chunksize parameter set to negative value,
            ! netcdf library default is used.
            if (output_file(nnnn)(1:6) == 'netcdf') then
               if (ichunk2d(grid_id) == 0) then
                  if( wrt_int_state%mype == 0 ) &
                    ichunk2d(grid_id) = wrt_int_state%out_grid_info(grid_id)%i_end - wrt_int_state%out_grid_info(grid_id)%i_start + 1
                  call mpi_bcast(ichunk2d(grid_id),1,mpi_integer,0,wrt_mpi_comm,rc)
               endif
               if (jchunk2d(grid_id) == 0) then
                  if( wrt_int_state%mype == 0 ) &
                    jchunk2d(grid_id) = wrt_int_state%out_grid_info(grid_id)%j_end - wrt_int_state%out_grid_info(grid_id)%j_start + 1
                  call mpi_bcast(jchunk2d(grid_id),1,mpi_integer,0,wrt_mpi_comm,rc)
               endif
               if (ichunk3d(grid_id) == 0) then
                  if( wrt_int_state%mype == 0 ) &
                    ichunk3d(grid_id) = wrt_int_state%out_grid_info(grid_id)%i_end - wrt_int_state%out_grid_info(grid_id)%i_start + 1
                  call mpi_bcast(ichunk3d(grid_id),1,mpi_integer,0,wrt_mpi_comm,rc)
               endif
               if (jchunk3d(grid_id) == 0) then
                  if( wrt_int_state%mype == 0 ) &
                    jchunk3d(grid_id) = wrt_int_state%out_grid_info(grid_id)%j_end - wrt_int_state%out_grid_info(grid_id)%j_start + 1
                  call mpi_bcast(jchunk3d(grid_id),1,mpi_integer,0,wrt_mpi_comm,rc)
               endif
               if (kchunk3d(grid_id) == 0 .and. nbdl == 1) then
                  if( wrt_int_state%mype == 0 )  then
                    call ESMF_FieldBundleGet(wrt_int_state%wrtFB(nbdl), grid=wrtGrid)
                    call ESMF_AttributeGet(wrtGrid, convention="NetCDF", purpose="FV3", &
                            attnestflag=ESMF_ATTNEST_OFF, name='pfull', &
                            itemCount=kchunk3d(grid_id), rc=rc)
                  endif
                  call mpi_bcast(kchunk3d(grid_id),1,mpi_integer,0,wrt_mpi_comm,rc)
               endif
               if (lprnt) then
                  print *,'ichunk2d,jchunk2d',ichunk2d(grid_id),jchunk2d(grid_id)
                  print *,'ichunk3d,jchunk3d,kchunk3d',ichunk3d(grid_id),jchunk3d(grid_id),kchunk3d(grid_id)
               endif
            endif

            if (is_restart_bundle) then
              write(time_restart,'(I4,I2.2,I2.2,".",I2.2,I2.2,I2.2)') cdate(1:6)

              ! strip leading 'restart_' from a bundle name and replace it with a directory name 'RESTART/' to create actual file name
              filename = 'RESTART/'//trim(time_restart)//'.'//trim(wrtFBName(9:))//'.nc'

              ! I hate this kind of inconsistencies
              ! If it's a restart bundle and the output grid is not cubed sphere and the output restart file is
              ! from dycore (ie. fv_core, fv_srf_wnd, fv_tracer) append 'tile1' to the end of the file name.
              ! As opposed to physics restart files (phy_data, sfc_data) which do not have 'tile1' appended.
              ! Why can't we have consistent naming?

              if (grid_id > 1) then
                if (top_parent_is_global) then
                  write(tile_id,'(I0)') 6 + grid_id - 1
                else
                  write(tile_id,'(I0)') grid_id
                endif
                filename = 'RESTART/'//trim(time_restart)//'.'//trim(wrtFBName(9:))//'.tile'//trim(tile_id)//'.nc'
              else
                if (.not. top_parent_is_global) then ! non cubed sphere restart bundles
                  if (wrtFBName(9:11) == 'fv_') then ! 'dynamics' restart bundles, append 'tile1'
                    filename = 'RESTART/'//trim(time_restart)//'.'//trim(wrtFBName(9:))//'.tile1'//'.nc'
                  endif
                endif
              endif

            else ! history bundle
              filename = trim(wrtFBName)//'f'//trim(cfhour)//'.nc'
            endif
            if(mype == lead_write_task) print *,'in wrt run,filename= ',nbdl,trim(filename)

  !
  ! set the time Attribute on the grid to carry it into the lower levels
            call ESMF_FieldBundleGet(file_bundle, grid=fbgrid, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_AttributeSet(fbgrid, convention="NetCDF", purpose="FV3", &
                                 name="time", value=nfhour, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            call ESMF_AttributeSet(fbgrid, convention="NetCDF", purpose="FV3", &
                                 name="time_iso", value=trim(time_iso), rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  !*** write out grid bundle:
  ! Provide log message indicating which wrtComp is active
            call ESMF_LogWrite("before Write component before gridFB ", ESMF_LOGMSG_INFO, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            if (trim(output_file(nnnn)) == 'netcdf') then
              use_parallel_netcdf = .false.
            else if (trim(output_file(nnnn)) == 'netcdf_parallel') then
              use_parallel_netcdf = .true.
            else
              call ESMF_LogWrite("wrt_run: Unknown output_file",ESMF_LOGMSG_ERROR,rc=RC)
              call ESMF_Finalize(endflag=ESMF_END_ABORT)
            endif

            wbeg = MPI_Wtime()

            if (is_restart_bundle) then ! restart bundle
              ! restart bundles are always on forecast grid, either cubed sphere or regional/nest

              call ESMF_FieldBundleGet(wrt_int_state%wrtFB(nbdl), grid=grid, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              call ESMF_GridGet(grid, tileCount=tileCount, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              if (tileCount == 6) then ! restart bundle is on cubed sphere
                call ESMFproto_FieldBundleWrite(gridFB, filename=trim(filename),               &
                                                convention="NetCDF", purpose="FV3",            &
                                                status=ESMF_FILESTATUS_REPLACE,                &
                                                state=stateGridFB, comps=compsGridFB,rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

                call ESMFproto_FieldBundleWrite(wrt_int_state%wrtFB(nbdl),                     &
                                                filename=trim(filename), convention="NetCDF",  &
                                                purpose="FV3", status=ESMF_FILESTATUS_OLD,     &
                                                timeslice=step, state=optimize(nbdl)%state,    &
                                                comps=optimize(nbdl)%comps, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              else
                call  write_restart_netcdf(wrt_int_state%wrtFB(nbdl),trim(filename), &
                                    .true., wrt_mpi_comm, wrt_int_state%mype, &
                                    rc)
              endif ! cubed sphere vs. regional/nest write grid

              restart_written = .true.

            else ! history bundle
            if (trim(output_grid(grid_id)) == 'cubed_sphere_grid') then

              if (trim(output_file(nnnn)) == 'netcdf_parallel') then
                call write_netcdf(wrt_int_state%wrtFB(nbdl),trim(filename), &
                                 .true., wrt_mpi_comm,wrt_int_state%mype, &
                                 grid_id,rc)
              else
                call ESMFproto_FieldBundleWrite(gridFB, filename=trim(filename),               &
                                                convention="NetCDF", purpose="FV3",            &
                                                status=ESMF_FILESTATUS_REPLACE,                &
                                                state=stateGridFB, comps=compsGridFB,rc=rc)

                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

                call ESMFproto_FieldBundleWrite(wrt_int_state%wrtFB(nbdl),                     &
                                                filename=trim(filename), convention="NetCDF",  &
                                                purpose="FV3", status=ESMF_FILESTATUS_OLD,     &
                                                timeslice=step, state=optimize(nbdl)%state,    &
                                                comps=optimize(nbdl)%comps, rc=rc)

                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              end if

            else if (trim(output_grid(grid_id)) == 'gaussian_grid' .or. &
                     trim(output_grid(grid_id)) == 'global_latlon') then

              call write_netcdf(wrt_int_state%wrtFB(nbdl),trim(filename), &
                               use_parallel_netcdf, wrt_mpi_comm,wrt_int_state%mype, &
                               grid_id,rc)

            else if (trim(output_grid(grid_id)) == 'regional_latlon' .or.        &
                     trim(output_grid(grid_id)) == 'regional_latlon_moving' .or. &
                     trim(output_grid(grid_id)) == 'rotated_latlon'  .or.        &
                     trim(output_grid(grid_id)) == 'rotated_latlon_moving' .or.  &
                     trim(output_grid(grid_id)) == 'lambert_conformal') then

              !mask fields according to sfc pressure
              if( .not. lmask_fields ) then
                call mask_fields(file_bundle,rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
              endif

              call write_netcdf(wrt_int_state%wrtFB(nbdl),trim(filename), &
                                use_parallel_netcdf, wrt_mpi_comm,wrt_int_state%mype, &
                                grid_id,rc)

            else ! unknown output_grid

              call ESMF_LogWrite("wrt_run: Unknown output_grid",ESMF_LOGMSG_ERROR,rc=RC)
              call ESMF_Finalize(endflag=ESMF_END_ABORT)

            endif
            endif  ! restart or history bundle
            wend = MPI_Wtime()
            if (lprnt) then
              write(*,'(A56,A,F10.5,A,I4.2,A,I2.2,1X,A)') trim(filename),' write time is ',wend-wbeg  &
                     ,' at fcst ',NF_HOURS,':',NF_MINUTES
            endif

          enddo file_loop_all

          if (out_phase == 1 .and. mype == lead_write_task) then
            !** write history log file
            open(newunit=nolog, file='log.atm.f'//trim(cfhour))
            write(nolog,"('completed: fv3atm')")
            write(nolog,"('forecast hour: ',f10.3)") nfhour
            write(nolog,"('valid time: ',6(i4,2x))") wrt_int_state%fdate(1:6)
            close(nolog)
          endif

          if (out_phase == 2 .and. restart_written .and. mype == lead_write_task) then
            !**  write coupler.res log file
            open(newunit=nolog, file='RESTART/'//trim(time_restart)//'.coupler.res')
            write(nolog,"(i6,8x,a)") calendar_type , &
                 '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'
            write(nolog,"(6i6,8x,a)") start_time(1:6), &
                 'Model start time:   year, month, day, hour, minute, second'
            write(nolog,"(6i6,8x,a)") wrt_int_state%fdate(1:6), &
                 'Current model time: year, month, day, hour, minute, second'
            close(nolog)
          endif

        enddo two_phase_loop
      endif ! if ( wrt_int_state%output_history )

      call ESMF_VMBarrier(VM, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      write_run_tim = MPI_Wtime() - tbeg

      IF (lprnt) THEN
        write(*,'(A56,A,F10.5,A,I4.2,A,I2.2,1X,A)')'------- total',' write time is ',write_run_tim &
                 ,' at Fcst ',NF_HOURS,':',NF_MINUTES
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE wrt_run
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      subroutine wrt_finalize(wrt_comp, imp_state_write, exp_state_write, clock, rc)
!
!-----------------------------------------------------------------------
!***  finalize the write gridded component.
!-----------------------------------------------------------------------
!
!***  HISTORY
!       Feb 2017:  J. Wang  - deallocate for fv3
!
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp)            :: wrt_comp
      type(ESMF_State)               :: imp_state_write, exp_state_write
      type(ESMF_Clock)               :: clock
      integer,intent(out)            :: rc
!
!***  local variables
!
      integer                        :: stat
      type(write_wrap)               :: wrap
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  retrieve the write component's esmf internal state(used later for
!***  post finalization)
!-----------------------------------------------------------------------
!
      call ESMF_GridCompGetInternalState(wrt_comp, wrap, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      deallocate(wrap%write_int_state,stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
          msg="Deallocation of internal state memory failed.", &
          line=__LINE__, file=__FILE__)) return

      call fms_end
!
!-----------------------------------------------------------------------
!
    end subroutine wrt_finalize
!
!-----------------------------------------------------------------------
!
   subroutine recover_fields(file_bundle,rc)

     type(ESMF_FieldBundle), intent(in)              :: file_bundle
     integer,                intent(out),   optional :: rc
!
     real, parameter   :: rdgas = 287.04, grav = 9.80
     real, parameter   :: stndrd_atmos_ps = 101325.
     real, parameter   :: stndrd_atmos_lapse = 0.0065

     integer i,j,k,ifld,fieldCount,nstt,nend,fieldDimCount,gridDimCount
     integer istart,iend,jstart,jend,kstart,kend
     logical uPresent, vPresent
     type(ESMF_Grid)  fieldGrid
     type(ESMF_Field)  ufield, vfield
     type(ESMF_TypeKind_Flag) typekind
     character(100) fieldName,uwindname,vwindname
     type(ESMF_Field),   allocatable  :: fcstField(:)
     real(ESMF_KIND_R4), dimension(:,:),     pointer  :: lonr4, latr4
     real(ESMF_KIND_R8), dimension(:,:),     pointer  :: lon, lat
     real(ESMF_KIND_R8), dimension(:,:),     pointer  :: lonloc, latloc
     real(ESMF_KIND_R4), dimension(:,:),     pointer  :: pressfc
     real(ESMF_KIND_R4), dimension(:,:),     pointer  :: uwind2dr4,vwind2dr4
     real(ESMF_KIND_R4), dimension(:,:,:),   pointer  :: uwind3dr4,vwind3dr4
     real(ESMF_KIND_R4), dimension(:,:,:),   pointer  :: cart3dPtr2dr4
     real(ESMF_KIND_R4), dimension(:,:,:,:), pointer  :: cart3dPtr3dr4
     real(ESMF_KIND_R8) :: coslon, sinlon, sinlat

     type(ESMF_Array) :: lon_array, lat_array
!
! get filed count
     call ESMF_FieldBundleGet(file_bundle, fieldCount=fieldCount, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     if (fieldCount == 0) return

     call ESMF_FieldBundleGet(file_bundle, grid=fieldGrid, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     call ESMF_LogWrite("call recover field on wrt comp",ESMF_LOGMSG_INFO,rc=RC)
     call ESMF_GridGet(fieldgrid, dimCount=gridDimCount, rc=rc)

     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     call ESMF_LogWrite("call recover field get coord 1",ESMF_LOGMSG_INFO,rc=RC)

     call ESMF_GridGetCoord(fieldgrid, coordDim=1, array=lon_array, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
     call ESMF_ArrayGet(lon_array, typekind=typekind, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     if (typekind == ESMF_TYPEKIND_R4) then
        call ESMF_GridGetCoord(fieldgrid, coordDim=1, farrayPtr=lonr4, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        allocate(lon(lbound(lonr4,1):ubound(lonr4,1),lbound(lonr4,2):ubound(lonr4,2)))
        lon = lonr4
     else if (typekind == ESMF_TYPEKIND_R8) then
        call ESMF_GridGetCoord(fieldgrid, coordDim=1, farrayPtr=lon, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
     else
        write(0,*)'lon_array unknown typekind'
        rc = 1
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
     endif


     allocate(lonloc(lbound(lon,1):ubound(lon,1),lbound(lon,2):ubound(lon,2)))
     istart = lbound(lon,1)
     iend   = ubound(lon,1)
     jstart = lbound(lon,2)
     jend   = ubound(lon,2)
!$omp parallel do default(none) shared(lon,lonloc,jstart,jend,istart,iend) &
!$omp             private(i,j)
     do j=jstart,jend
      do i=istart,iend
        lonloc(i,j) = lon(i,j) * pi/180.
      enddo
     enddo

     call ESMF_LogWrite("call recover field get coord 2",ESMF_LOGMSG_INFO,rc=RC)

     call ESMF_GridGetCoord(fieldgrid, coordDim=2, array=lat_array, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
     call ESMF_ArrayGet(lat_array, typekind=typekind, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

     if (typekind == ESMF_TYPEKIND_R4) then
        call ESMF_GridGetCoord(fieldgrid, coordDim=2, farrayPtr=latr4, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        allocate(lat(lbound(latr4,1):ubound(latr4,1),lbound(latr4,2):ubound(latr4,2)))
        lat = latr4
     else if (typekind == ESMF_TYPEKIND_R8) then
        call ESMF_GridGetCoord(fieldgrid, coordDim=2, farrayPtr=lat, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
     else
        write(0,*)'lon_array unknown typekind'
        rc = 1
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
     endif

     allocate(latloc(lbound(lat,1):ubound(lat,1),lbound(lat,2):ubound(lat,2)))
     istart = lbound(lat,1)
     iend   = ubound(lat,1)
     jstart = lbound(lat,2)
     jend   = ubound(lat,2)
!$omp parallel do default(none) shared(lat,latloc,jstart,jend,istart,iend) &
!$omp             private(i,j)
     do j=jstart,jend
      do i=istart,iend
        latloc(i,j) = lat(i,j) * pi/180.d0
      enddo
     enddo
!
     allocate(fcstField(fieldCount))
     call ESMF_LogWrite("call recover field get fcstField",ESMF_LOGMSG_INFO,rc=RC)
     call ESMF_FieldBundleGet(file_bundle, fieldList=fcstField, itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
!
     do ifld=1,fieldCount

       call ESMF_LogWrite("call recover field get fieldname, type dimcount",ESMF_LOGMSG_INFO,rc=RC)
       call ESMF_FieldGet(fcstField(ifld),name=fieldName,typekind=typekind,dimCount=fieldDimCount, rc=rc)

! convert back wind
       if(index(trim(fieldName),"vector")>0) then
         nstt = index(trim(fieldName),"wind")+4
         nend = index(trim(fieldName),"vector")-1
         if( nend>nstt ) then
           uwindname = 'u'//fieldName(nstt+1:nend)
           vwindname = 'v'//fieldName(nstt+1:nend)
         else
           uwindname = 'ugrd'
           vwindname = 'vgrd'
         endif
!         print *,'in get 3D vector wind, uwindname=',trim(uwindname),' v=', trim(vwindname),' fieldname=',trim(fieldname)
! get u , v wind
         call ESMF_LogWrite("call recover field get u, v field",ESMF_LOGMSG_INFO,rc=RC)
         call ESMF_FieldBundleGet(file_bundle,trim(uwindname),field=ufield,isPresent=uPresent,rc=rc)
         call ESMF_FieldBundleGet(file_bundle,trim(vwindname),field=vfield,isPresent=vPresent,rc=rc)
         if(.not. uPresent .or. .not.vPresent) then
           rc=990
           print *,' ERROR ,the local wind is not present! rc=', rc
           exit
         endif

! get field data
         if ( typekind == ESMF_TYPEKIND_R4 ) then
           if( fieldDimCount > gridDimCount+1 ) then
             call ESMF_LogWrite("call recover field get 3d card wind farray",ESMF_LOGMSG_INFO,rc=RC)
             call ESMF_FieldGet(fcstField(ifld), localDe=0, farrayPtr=cart3dPtr3dr4, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
             if( ubound(cart3dPtr3dr4,1)-lbound(cart3dPtr3dr4,1)+1/=3) then
               rc=991
               print *,'ERROR, 3D the vector dimension /= 3, rc=',rc
               exit
             endif
             iend   = ubound(cart3dPtr3dr4,2)
             istart = lbound(cart3dPtr3dr4,2)
             iend   = ubound(cart3dPtr3dr4,2)
             jstart = lbound(cart3dPtr3dr4,3)
             jend   = ubound(cart3dPtr3dr4,3)
             kstart = lbound(cart3dPtr3dr4,4)
             kend   = ubound(cart3dPtr3dr4,4)
             call ESMF_FieldGet(ufield, localDe=0, farrayPtr=uwind3dr4,rc=rc)
             call ESMF_FieldGet(vfield, localDe=0, farrayPtr=vwind3dr4,rc=rc)
! update u , v wind
!$omp parallel do default(shared) private(i,j,k,coslon,sinlon,sinlat)
             do k=kstart,kend
!$omp parallel do default(none) shared(uwind3dr4,vwind3dr4,lonloc,latloc,cart3dPtr3dr4,jstart,jend,istart,iend,k) &
!$omp             private(i,j,coslon,sinlon,sinlat)
               do j=jstart, jend
                 do i=istart, iend
                  coslon = cos(lonloc(i,j))
                  sinlon = sin(lonloc(i,j))
                  sinlat = sin(latloc(i,j))
                  uwind3dr4(i,j,k) = cart3dPtr3dr4(1,i,j,k) * coslon           &
                                   + cart3dPtr3dr4(2,i,j,k) * sinlon
                  vwind3dr4(i,j,k) =-cart3dPtr3dr4(1,i,j,k) * sinlat*sinlon    &
                                   + cart3dPtr3dr4(2,i,j,k) * sinlat*coslon    &
                                   + cart3dPtr3dr4(3,i,j,k) * cos(latloc(i,j))
                 enddo
               enddo
             enddo
           else
             call ESMF_FieldGet(fcstField(ifld), localDe=0, farrayPtr=cart3dPtr2dr4, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
             if( ubound(cart3dPtr2dr4,1)-lbound(cart3dPtr2dr4,1)+1 /= 3) then
               rc=991
               write(0,*) 'ERROR, 2D the vector dimension /= 3, rc=',rc
               call ESMF_Finalize(endflag=ESMF_END_ABORT)
             endif
             istart = lbound(cart3dPtr2dr4,2)
             iend   = ubound(cart3dPtr2dr4,2)
             jstart = lbound(cart3dPtr2dr4,3)
             jend   = ubound(cart3dPtr2dr4,3)

             call ESMF_FieldGet(ufield, localDe=0, farrayPtr=uwind2dr4,rc=rc)
             call ESMF_FieldGet(vfield, localDe=0, farrayPtr=vwind2dr4,rc=rc)
              ! update u , v wind
!$omp parallel do default(none) shared(uwind2dr4,vwind2dr4,lonloc,latloc,cart3dPtr2dr4,jstart,jend,istart,iend) &
!$omp             private(i,j,k,coslon,sinlon,sinlat)
             do j=jstart, jend
               do i=istart, iend
                  coslon = cos(lonloc(i,j))
                  sinlon = sin(lonloc(i,j))
                  sinlat = sin(latloc(i,j))
                  uwind2dr4(i,j) = cart3dPtr2dr4(1,i,j) * coslon         &
                                 + cart3dPtr2dr4(2,i,j) * sinlon
                  vwind2dr4(i,j) =-cart3dPtr2dr4(1,i,j) * sinlat*sinlon  &
                                 + cart3dPtr2dr4(2,i,j) * sinlat*coslon  &
                                 + cart3dPtr2dr4(3,i,j) * cos(latloc(i,j))
               enddo
             enddo
           endif
         endif
!     print *,'in 3DCartesian2wind, uwindname=', trim(uwindname),'uPresent =',uPresent, 'vPresent=',vPresent,'fieldDimCount=', &
!       fieldDimCount,'gridDimCount=',gridDimCount

! convert back surface pressure
       else if(index(trim(fieldName),"pressfc")>0) then
         call ESMF_FieldGet(fcstField(ifld),localDe=0, farrayPtr=pressfc, rc=rc)
         istart = lbound(pressfc,1)
         iend   = ubound(pressfc,1)
         jstart = lbound(pressfc,2)
         jend   = ubound(pressfc,2)
!$omp parallel do default(none) shared(pressfc,jstart,jend,istart,iend) private(i,j)
         do j=jstart, jend
           do i=istart, iend
             pressfc(i,j) = pressfc(i,j)**(grav/(rdgas*stndrd_atmos_lapse))*stndrd_atmos_ps
           enddo
         enddo
       endif
     enddo
!
     deallocate(fcstField)
     rc = 0

   end subroutine recover_fields
!
!-----------------------------------------------------------------------
!
   subroutine mask_fields(file_bundle,rc)

     type(ESMF_FieldBundle), intent(in)              :: file_bundle
     integer,                intent(out),   optional :: rc
!
     integer i,j,k,ifld,fieldCount,fieldDimCount,gridDimCount
     integer istart,iend,jstart,jend,kstart,kend
     type(ESMF_Grid)  fieldGrid
     type(ESMF_TypeKind_Flag) typekind
     type(ESMF_TypeKind_Flag) attTypeKind
     character(len=ESMF_MAXSTR) fieldName
     type(ESMF_Field),   allocatable                  :: fcstField(:)
     real(ESMF_KIND_R4), dimension(:,:),     pointer  :: var2dPtr2dr4
     real(ESMF_KIND_R4), dimension(:,:,:),   pointer  :: var3dPtr3dr4
     real(ESMF_KIND_R4), dimension(:,:,:),   pointer  :: vect3dPtr2dr4
     real(ESMF_KIND_R4), dimension(:,:,:,:), pointer  :: vect4dPtr3dr4
     real(ESMF_KIND_R4), dimension(:,:), allocatable  :: maskwrt

     logical :: mvispresent=.false.
     real(ESMF_KIND_R4) :: missing_value_r4=-1.e+10
     real(ESMF_KIND_R8) :: missing_value_r8=9.99e20
     character(len=ESMF_MAXSTR) :: msg

     call ESMF_LogWrite("call mask field on wrt comp",ESMF_LOGMSG_INFO,rc=RC)

! get fieldCount
     call ESMF_FieldBundleGet(file_bundle, fieldCount=fieldCount, &
         grid=fieldGrid, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out
! get gridDimCount
     call ESMF_GridGet(fieldgrid, dimCount=gridDimCount, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return  ! bail out

     allocate(fcstField(fieldCount))
     call ESMF_LogWrite("call mask field get fcstField",ESMF_LOGMSG_INFO,rc=RC)
     call ESMF_FieldBundleGet(file_bundle, fieldList=fcstField, itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)

! generate the maskwrt according to surface pressure
     do ifld=1,fieldCount
       !call ESMF_LogWrite("call mask field get fieldname, type dimcount",ESMF_LOGMSG_INFO,rc=RC)
       call ESMF_FieldGet(fcstField(ifld),name=fieldName,typekind=typekind,dimCount=fieldDimCount, rc=rc)
       !write(msg,*) 'fieldName,typekind,fieldDimCount=',trim(fieldName),typekind,fieldDimCount
       !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
       if (.not. allocated(maskwrt)) then
         if ( typekind == ESMF_TYPEKIND_R4 .and. fieldDimCount == gridDimCount) then
           call ESMF_FieldGet(fcstField(ifld),localDe=0, farrayPtr=var2dPtr2dr4, rc=rc)
           istart = lbound(var2dPtr2dr4,1)
           iend   = ubound(var2dPtr2dr4,1)
           jstart = lbound(var2dPtr2dr4,2)
           jend   = ubound(var2dPtr2dr4,2)
           allocate(maskwrt(istart:iend,jstart:jend))
           maskwrt(istart:iend,jstart:jend)=1.0
         endif
       endif
       if(index(trim(fieldName),"pressfc")>0) then
         call ESMF_FieldGet(fcstField(ifld),localDe=0, farrayPtr=var2dPtr2dr4, rc=rc)
         istart = lbound(var2dPtr2dr4,1)
         iend   = ubound(var2dPtr2dr4,1)
         jstart = lbound(var2dPtr2dr4,2)
         jend   = ubound(var2dPtr2dr4,2)
         if (.not. allocated(maskwrt)) then
           allocate(maskwrt(istart:iend,jstart:jend))
           maskwrt(istart:iend,jstart:jend)=1.0
         endif
!$omp parallel do default(shared) private(i,j)
         do j=jstart, jend
           do i=istart, iend
             if(abs(var2dPtr2dr4(i,j)-0.) < 1.0e-6) maskwrt(i,j)=0.
           enddo
         enddo
         call ESMF_LogWrite("call mask field pressfc found, maskwrt generated",ESMF_LOGMSG_INFO,rc=RC)
         exit
       endif
     enddo

! loop to mask all fields according to maskwrt
     do ifld=1,fieldCount
       !call ESMF_LogWrite("call mask field get fieldname, type dimcount",ESMF_LOGMSG_INFO,rc=RC)
       call ESMF_FieldGet(fcstField(ifld),name=fieldName,typekind=typekind,dimCount=fieldDimCount, rc=rc)
       !write(msg,*) 'fieldName,typekind,fieldDimCount=',trim(fieldName),typekind,fieldDimCount
       !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
       ! For vector fields
       if(index(trim(fieldName),"vector")>0) then
         ! Only work on ESMF_TYPEKIND_R4 fields for now
         if ( typekind == ESMF_TYPEKIND_R4 ) then
           ! 3-d vector fields with 4-d arrays
           if( fieldDimCount > gridDimCount+1 ) then
             !call ESMF_LogWrite("call mask field get vector 3d farray",ESMF_LOGMSG_INFO,rc=RC)
             call ESMF_FieldGet(fcstField(ifld), localDe=0, farrayPtr=vect4dPtr3dr4, rc=rc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                 line=__LINE__, file=__FILE__)) return  ! bail out
             if( ubound(vect4dPtr3dr4,1)-lbound(vect4dPtr3dr4,1)+1/=3 ) then
               rc=991
               write(0,*) 'ERROR, 3D the vector dimension /= 3, rc=',rc
               call ESMF_Finalize(endflag=ESMF_END_ABORT)
             endif
             ! Get the _FillValue from the field attribute if exists
             call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                 name="_FillValue", typekind=attTypeKind, isPresent=mvispresent, rc=rc)
             !write(msg,*) 'fieldName,attTypeKind,isPresent=',trim(fieldName),attTypeKind,mvispresent
             !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
             if ( mvispresent ) then
               if (attTypeKind==ESMF_TYPEKIND_R4) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r4, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r4,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               else if (attTypeKind==ESMF_TYPEKIND_R8) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r8, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r8,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               endif
               istart = lbound(vect4dPtr3dr4,2)
               iend   = ubound(vect4dPtr3dr4,2)
               jstart = lbound(vect4dPtr3dr4,3)
               jend   = ubound(vect4dPtr3dr4,3)
               kstart = lbound(vect4dPtr3dr4,4)
               kend   = ubound(vect4dPtr3dr4,4)
!$omp parallel do default(shared) private(i,j,k)
               do k=kstart,kend
                 do j=jstart, jend
                   do i=istart, iend
                     if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R4) vect4dPtr3dr4(:,i,j,k)=missing_value_r4
                     if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R8) vect4dPtr3dr4(:,i,j,k)=missing_value_r8
                   enddo
                 enddo
               enddo
             endif !mvispresent
           ! 2-d vector fields with 3-d arrays
           else
             call ESMF_FieldGet(fcstField(ifld), localDe=0, farrayPtr=vect3dPtr2dr4, rc=rc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                 line=__LINE__, file=__FILE__)) return  ! bail out
             if( ubound(vect3dPtr2dr4,1)-lbound(vect3dPtr2dr4,1)+1 /= 3 ) then
               rc=991
               write(0,*) 'ERROR, 2D the vector dimension /= 3, rc=',rc
               call ESMF_Finalize(endflag=ESMF_END_ABORT)
             endif
             ! Get the _FillValue from the field attribute if exists
             call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                 name="_FillValue", typekind=attTypeKind, isPresent=mvispresent, rc=rc)
             !write(msg,*) 'fieldName,attTypeKind,isPresent=',trim(fieldName),attTypeKind,mvispresent
             !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
             if ( mvispresent ) then
               if (attTypeKind==ESMF_TYPEKIND_R4) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r4, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r4,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               else if (attTypeKind==ESMF_TYPEKIND_R8) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r8, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r8,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               endif
               istart = lbound(vect3dPtr2dr4,2)
               iend   = ubound(vect3dPtr2dr4,2)
               jstart = lbound(vect3dPtr2dr4,3)
               jend   = ubound(vect3dPtr2dr4,3)
!$omp parallel do default(shared) private(i,j)
               do j=jstart, jend
                 do i=istart, iend
                   if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R4) vect3dPtr2dr4(:,i,j)=missing_value_r4
                   if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R8) vect3dPtr2dr4(:,i,j)=missing_value_r8
                 enddo
               enddo
             endif !mvispresent
           endif
         endif
! For non-vector fields
       else
         ! Only work on ESMF_TYPEKIND_R4 fields for now
         if ( typekind == ESMF_TYPEKIND_R4 ) then
           ! 2-d fields
           if(fieldDimCount == gridDimCount) then
             call ESMF_FieldGet(fcstField(ifld),localDe=0, farrayPtr=var2dPtr2dr4, rc=rc)
             ! Get the _FillValue from the field attribute if exists
             call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                 name="_FillValue", typekind=attTypeKind, isPresent=mvispresent, rc=rc)
             !write(msg,*) 'fieldName,attTypeKind,isPresent=',trim(fieldName),attTypeKind,mvispresent
             !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
             if ( mvispresent ) then
               if (attTypeKind==ESMF_TYPEKIND_R4) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r4, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r4,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               else if (attTypeKind==ESMF_TYPEKIND_R8) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r8, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r8,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               endif
               istart = lbound(var2dPtr2dr4,1)
               iend   = ubound(var2dPtr2dr4,1)
               jstart = lbound(var2dPtr2dr4,2)
               jend   = ubound(var2dPtr2dr4,2)
!$omp parallel do default(shared) private(i,j)
               do j=jstart, jend
                 do i=istart, iend
                   if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R4) var2dPtr2dr4(i,j)=missing_value_r4
                   if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R8) var2dPtr2dr4(i,j)=missing_value_r8
                 enddo
               enddo
             endif !mvispresent
           ! 3-d fields
           else if(fieldDimCount == gridDimCount+1) then
             call ESMF_FieldGet(fcstField(ifld),localDe=0, farrayPtr=var3dPtr3dr4, rc=rc)
             ! Get the _FillValue from the field attribute if exists
             call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                 name="_FillValue", typekind=attTypeKind, isPresent=mvispresent, rc=rc)
             !write(msg,*) 'fieldName,attTypeKind,isPresent=',trim(fieldName),attTypeKind,mvispresent
             !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
             if ( mvispresent ) then
               if (attTypeKind==ESMF_TYPEKIND_R4) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r4, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r4,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               else if (attTypeKind==ESMF_TYPEKIND_R8) then
                 call ESMF_AttributeGet(fcstField(ifld), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=missing_value_r8, isPresent=mvispresent, rc=rc)
                 !write(msg,*) 'fieldName,_FillValue,isPresent=',trim(fieldName),missing_value_r8,mvispresent
                 !call ESMF_LogWrite("call mask field: "//trim(msg),ESMF_LOGMSG_INFO,rc=RC)
               endif
               istart = lbound(var3dPtr3dr4,1)
               iend   = ubound(var3dPtr3dr4,1)
               jstart = lbound(var3dPtr3dr4,2)
               jend   = ubound(var3dPtr3dr4,2)
               kstart = lbound(var3dPtr3dr4,3)
               kend   = ubound(var3dPtr3dr4,3)
!$omp parallel do default(shared) private(i,j,k)
               do k=kstart,kend
                 do j=jstart, jend
                   do i=istart, iend
                     if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R4) var3dPtr3dr4(i,j,k)=missing_value_r4
                     if (maskwrt(i,j)<1.0 .and. attTypeKind==ESMF_TYPEKIND_R8) var3dPtr3dr4(i,j,k)=missing_value_r8
                   enddo
                 enddo
               enddo
             endif !mvispresent
           endif
         endif
       endif
     enddo
!
     deallocate(maskwrt)
     deallocate(fcstField)
     rc = 0

   end subroutine mask_fields
!
!-----------------------------------------------------------------------
!

  subroutine ESMFproto_FieldBundleWrite(fieldbundle, fileName, &
    convention, purpose, status, timeslice, state, comps, rc)
    type(ESMF_FieldBundle),     intent(in)              :: fieldbundle
    character(*),               intent(in)              :: fileName
    character(*),               intent(in),    optional :: convention
    character(*),               intent(in),    optional :: purpose
    type(ESMF_FileStatus_Flag), intent(in),    optional :: status
    integer,                    intent(in),    optional :: timeslice
    type(ESMF_State),           intent(inout), optional :: state
    type(ESMF_GridComp), allocatable, intent(inout), optional :: comps(:)
    integer,                    intent(out),   optional :: rc

    ! Prototype multi-tile implementation for FieldBundleWrite().
    ! Produces as many output files as there are tiles. The naming of the
    ! output files is such that the string in the fileName argument is used
    ! as the basis. If fileName ends with ".nc", then this suffix is replaced
    ! by ".tileN.nc", where "N" is the tile number. If fileName does not
    ! end in ".nc", then ".tileN.nc" will simply be appended.
    !
    ! Restrictions:
    !   - All Fields in the FieldBundle must have the same tileCount

    integer                             :: i, j, ind
    integer                             :: fieldCount, tileCount, itemCount
    type(ESMF_Field), allocatable       :: fieldList(:), tileFieldList(:)
    type(ESMF_Grid)                     :: grid
    type(ESMF_Array)                    :: array
    type(ESMF_DistGrid)                 :: distgrid
    type(ESMF_DELayout)                 :: delayout
    type(ESMF_FieldBundle)              :: wrtTileFB
    type(ESMF_FieldBundle), allocatable :: wrtTileFBList(:)
    character(len=80), allocatable      :: itemNameList(:)
    logical                             :: stateIsEmpty
    character(len=80)                   :: tileFileName
    character(len=80)                   :: statusStr
    integer, allocatable                :: petList(:)
    character(1024)                     :: msgString
    type(ESMF_State), allocatable       :: ioState(:)
    integer                             :: timesliceOpt
    integer                             :: urc

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_VMLogMemInfo("Entering ESMFproto_FieldBundleWrite",rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! query number of fields in fieldbundle
    call ESMF_FieldBundleGet(fieldbundle, fieldCount=fieldCount, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! early successful exit if there are no fields present
    if (fieldCount==0) return
    ! obtain list of fields in the fieldbundle

    allocate(fieldList(fieldCount), tileFieldList(fieldCount))
    call ESMF_FieldBundleGet(fieldbundle, fieldList=fieldList, &
                             itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    ! determine tileCount by looking at first field

    call ESMF_FieldGet(fieldList(1), array=array, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_ArrayGet(array, tileCount=tileCount, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! deal with optional state argument
    stateIsEmpty = .true.
    if (present(state)) then
      if (.not.ESMF_StateIsCreated(state, rc=rc)) then
        state = ESMF_StateCreate(rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
      call ESMF_StateGet(state, itemCount=itemCount, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (itemCount /= 0) then
        stateIsEmpty = .false.
        if (itemCount /= tileCount) then
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
                                msg="Number of items in state must match number of tiles.", &
                                line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
        endif
        allocate(itemNameList(itemCount),wrtTileFBList(itemCount))
        call ESMF_StateGet(state, itemNameList=itemNameList, &
                           itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        do i=1, itemCount
          call ESMF_StateGet(state, itemName=itemNameList(i), &
                             fieldbundle=wrtTileFBList(i), rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        enddo
      endif
    endif

    ! loop over all the tiles and construct a tile specific fieldbundle
    call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() before tileCount-loop",&
      ESMF_LOGMSG_INFO, rc=rc)

    !TODO: remove this once comps is hidden within state
    if (present(comps)) then
      allocate(ioState(tileCount))
      if (.not.allocated(comps)) then
        ! first-write
        allocate(comps(tileCount))
        do i=1, tileCount
          call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() before "// &
                             "ESMFproto_FieldMakeSingleTile() w/ petList", &
                              ESMF_LOGMSG_INFO, rc=rc)
          do j=1, fieldCount
            ! access only tile specific part of field
            call ESMFproto_FieldMakeSingleTile(fieldList(j), tile=i, &
                                               tileField=tileFieldList(j), petList=petList, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          enddo
!          write(msgString, *) petList
!          call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() after "// &
!            "ESMFproto_FieldMakeSingleTile(), petList:"//trim(msgString), &
!            ESMF_LOGMSG_INFO, rc=rc)
          ! create component to handle this tile I/O
          call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() before "// &
                             "tile-component creation", ESMF_LOGMSG_INFO, rc=rc)
          comps(i) = ESMF_GridCompCreate(petList=petList, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          ! convention
          call ESMF_AttributeSet(comps(i), name="convention", &
                                 value=convention, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          ! purpose
          call ESMF_AttributeSet(comps(i), name="purpose", &
                                 value=purpose, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          ! timeslice
          timesliceOpt = -1 ! init
          if (present(timeslice)) timesliceOpt = timeslice
          call ESMF_AttributeSet(comps(i), name="timeslice", &
                                 value=timesliceOpt, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          ! status
          statusStr = "ESMF_FILESTATUS_UNKNOWN" ! default
          if (present(status)) then
            if (status==ESMF_FILESTATUS_UNKNOWN) then
              statusStr="ESMF_FILESTATUS_UNKNOWN" ! default
            else if (status==ESMF_FILESTATUS_NEW) then
              statusStr="ESMF_FILESTATUS_NEW" ! default
            else if (status==ESMF_FILESTATUS_OLD) then
              statusStr="ESMF_FILESTATUS_OLD" ! default
            else if (status==ESMF_FILESTATUS_REPLACE) then
              statusStr="ESMF_FILESTATUS_REPLACE" ! default
            endif
          endif
          call ESMF_AttributeSet(comps(i), name="status", &
                                 value=statusStr, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_GridCompSetServices(comps(i), ioCompSS, userRc=urc, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() after "// &
                             "tile-component creation", ESMF_LOGMSG_INFO, rc=rc)
        enddo
      endif
    endif

    do i=1, tileCount
      if (stateIsEmpty) then
        ! loop over all the fields and add tile specific part to fieldbundle
        call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() before "// &
                           "ESMFproto_FieldMakeSingleTile()", ESMF_LOGMSG_INFO, rc=rc)
        do j=1, fieldCount
          ! access only tile specific part of field
          call ESMFproto_FieldMakeSingleTile(fieldList(j), tile=i, &
                                             tileField=tileFieldList(j), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        enddo
        call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() after "// &
                           "ESMFproto_FieldMakeSingleTile()", ESMF_LOGMSG_INFO, rc=rc)
        ! create tile specific fieldbundle
        wrtTileFB = ESMF_FieldBundleCreate(fieldList=tileFieldList, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        ! ensure global attributes on the fieldbundle are passed on by reference
        call ESMF_AttributeCopy(fieldbundle, wrtTileFB, &
                                attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        ! store this fieldbundle in state if present
        if (present(state)) then
          call ESMF_StateAdd(state, fieldbundleList=(/wrtTileFB/), rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        endif
      else
        ! state brought in existing fieldbundles
        wrtTileFB = wrtTileFBList(i)
      endif
      ! write out the tile specific fieldbundle
      if(tileCount>1) then
        ind=min(index(trim(fileName),".nc",.true.)-1,len_trim(fileName))
        write(tileFileName, "(A,A,I1,A)") fileName(1:ind), ".tile", i, ".nc"
      else
        tileFileName=trim(fileName)
      endif
      if (present(comps)) then
        ioState(i) = ESMF_StateCreate(rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_StateAdd(ioState(i), (/wrtTileFB/), rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_AttributeSet(comps(i), name="tileFileName", value=tileFileName, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_GridCompRun(comps(i), importState=ioState(i), userRc=urc, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_StateDestroy(ioState(i), rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__))  return
      else
        call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() before "//      &
                           "ESMF_FieldBundleWrite(): "//trim(tileFileName), &
                            ESMF_LOGMSG_INFO, rc=rc)
        call ESMF_FieldBundleWrite(fieldbundle=wrtTileFB, fileName=tileFileName,          &
                                   convention=convention, purpose=purpose, status=status, &
                                   timeslice=timeslice, overwrite=.true., rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() after "// &
                           "ESMF_FieldBundleWrite()", ESMF_LOGMSG_INFO, rc=rc)
      endif
      if (.not.present(state)) then
        ! local garbage collection of fields
        do j=1, fieldCount
          call ESMF_FieldGet(tileFieldList(j), array=array, grid=grid, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_FieldDestroy(tileFieldList(j), noGarbage=.true., rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_GridDestroy(grid, noGarbage=.true., rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_ArrayGet(array, distgrid=distgrid, delayout=delayout, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_ArrayDestroy(array, noGarbage=.true., rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_DistGridDestroy(distgrid, noGarbage=.true., rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_DELayoutDestroy(delayout, noGarbage=.true., rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        enddo
        ! destroy tile specific fieldbundle
        call ESMF_FieldBundleDestroy(wrtTileFB, noGarbage=.true., rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
    enddo
    if (present(comps)) then
      deallocate(ioState)
    endif
    call ESMF_LogWrite("In ESMFproto_FieldBundleWrite() after tileCount-loop",&
                        ESMF_LOGMSG_INFO, rc=rc)

    ! deallocate temporary lists
    deallocate(fieldList, tileFieldList)

    call ESMF_VMLogMemInfo("Exiting ESMFproto_FieldBundleWrite",rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine ESMFproto_FieldBundleWrite

  !-----------------------------------------------------------------------------

  subroutine ioCompSS(comp, rc)
    type(ESMF_GridComp)   :: comp
    integer, intent(out)  :: rc

    rc = ESMF_SUCCESS

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
                                    userRoutine=ioCompRun, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine ioCompRun(comp, importState, exportState, clock, rc)
    use netcdf

    type(ESMF_GridComp)   :: comp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    type(ESMF_FieldBundle)           :: wrtTileFB
    character(len=80)                :: tileFileName
    character(len=80)                :: convention
    character(len=80)                :: purpose
    integer                          :: timeslice
    character(len=80)                :: statusStr
    type(ESMF_FileStatus_Flag)       :: status
    character(len=80)                :: itemNameList(1)

    integer                          :: localPet, petCount, i, j, k, ind
    type(ESMF_Grid)                  :: grid
    real(ESMF_KIND_R4), allocatable  :: valueListr4(:)
    real(ESMF_KIND_R8), allocatable  :: valueListr8(:)
    integer                          :: valueCount, fieldCount, udimCount
    character(80),      allocatable  :: udimList(:)
    integer                          :: ncerr, ncid, dimid, varid
    type(ESMF_Field),   allocatable  :: fieldList(:)
    type(ESMF_Field)                 :: field
    logical                          :: isPresent
    real(ESMF_KIND_R8)               :: time
    integer                          :: itemCount, attCount
    character(len=80),  allocatable  :: attNameList(:)
    character(len=80)                :: attName
    type(ESMF_TypeKind_Flag)         :: typekind
    character(len=80)                :: valueS
    integer                          :: valueI4
    real(ESMF_KIND_R4)               :: valueR4
    real(ESMF_KIND_R8)               :: valueR8
    logical                          :: thereAreVerticals
    integer                          :: ch_dimid, timeiso_varid
    character(len=ESMF_MAXSTR)       :: time_iso
    integer                          :: wrt_mpi_comm
    type(ESMF_VM)                    :: vm

    rc = ESMF_SUCCESS

    ! Access the FieldBundle
    call ESMF_StateGet(importState, itemNameList=itemNameList, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_StateGet(importState, itemName=itemNameList(1), &
                       fieldBundle=wrtTileFB, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! Access attributes on the component and use as parameters for Write()
    call ESMF_AttributeGet(comp, name="tileFileName", value=tileFileName, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeGet(comp, name="convention", value=convention, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeGet(comp, name="purpose", value=purpose, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeGet(comp, name="timeslice", value=timeslice, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeGet(comp, name="status", value=statusStr, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (trim(statusStr) == "ESMF_FILESTATUS_UNKNOWN") then
      status = ESMF_FILESTATUS_UNKNOWN
    else if (trim(statusStr) == "ESMF_FILESTATUS_NEW") then
      status = ESMF_FILESTATUS_NEW
    else if (trim(statusStr) == "ESMF_FILESTATUS_OLD") then
      status=ESMF_FILESTATUS_OLD
    else if (trim(statusStr) == "ESMF_FILESTATUS_REPLACE") then
      status = ESMF_FILESTATUS_REPLACE
    endif

    if ( tileFileName(1:7) == 'RESTART' ) then

      ! Write out restart files using write_restart_netcdf, then return from this subroutine
      if (timeslice == 1) then
        call ESMF_GridCompGet(comp, localPet=localPet, petCount=petCount, vm=vm, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        call ESMF_VMGet(vm=vm, mpiCommunicator=wrt_mpi_comm, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        if (petCount > 1) then
          call write_restart_netcdf(wrtTileFB, trim(tileFileName), .true., wrt_mpi_comm, localPet, rc)
        else
          call write_restart_netcdf(wrtTileFB, trim(tileFileName), .false., wrt_mpi_comm, localPet, rc)
        endif

      endif
      return
    endif

    call ESMF_LogWrite("In ioCompRun() before writing to: "// &
                       trim(tileFileName), ESMF_LOGMSG_INFO, rc=rc)

    if (status == ESMF_FILESTATUS_OLD) then
      ! This writes the vertical coordinates and the time dimension into the
      ! file. Doing this before the large data sets are written, assuming that
      ! the first time coming into ioCompRun() with this tileFileName, only
      ! the grid info is written. Second time in, with ESMF_FILESTATUS_OLD,
      ! the large data sets are written. That is when vertical and time info
      ! is also written.
      ! Hoping for better performance because file exists, but is still small.

      call ESMF_GridCompGet(comp, localPet=localPet, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      if (localPet==0) then
        ! do this work only on the root pet
        call ESMF_FieldBundleGet(wrtTileFB, grid=grid, fieldCount=fieldCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        allocate(fieldList(fieldCount))
        call ESMF_FieldBundleGet(wrtTileFB, fieldList=fieldList, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        ! open this tile's NetCDF file
        ncerr = nf90_open(tileFileName, NF90_WRITE, ncid=ncid)

        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        ! loop over all the fields in the bundle and handle their vertical dims

        thereAreVerticals = .false.
        do i=1, fieldCount
          field = fieldList(i)
          call ESMF_AttributeGetAttPack(field, convention="NetCDF", purpose="FV3", &
                                        isPresent=isPresent, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (.not.isPresent) cycle ! field does not have the AttPack
          call ESMF_AttributeGet(field, convention="NetCDF", purpose="FV3", &
                                 name="ESMF:ungridded_dim_labels", isPresent=isPresent, &
                                 itemCount=udimCount, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (udimCount==0 .or. .not.isPresent) cycle ! nothing there to do

          thereAreVerticals = .true.
          allocate(udimList(udimCount))
          call ESMF_AttributeGet(field, convention="NetCDF", purpose="FV3", &
                                 name="ESMF:ungridded_dim_labels", valueList=udimList, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ! loop over all ungridded dimension labels
          do k=1, udimCount
            call write_out_ungridded_dim_atts(dimLabel=trim(udimList(k)), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ! for restart files we store ungridded dimension labels in fields
            call write_out_ungridded_dim_atts_from_field(field, dimLabel=trim(udimList(k)), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          enddo
          deallocate(udimList)
        enddo ! fieldCount
        deallocate(fieldList)
        if (thereAreVerticals) then
          ! see if the vertical_dim_labels attribute exists on the grid, and
          ! if so access it and write out verticals accordingly
          call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                 name="vertical_dim_labels", isPresent=isPresent, &
                                 itemCount=udimCount, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (isPresent .and. (udimCount>0) ) then
            allocate(udimList(udimCount))
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name="vertical_dim_labels", valueList=udimList, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ! loop over all ungridded dimension labels
            do k=1, udimCount
              call write_out_ungridded_dim_atts(dimLabel=trim(udimList(k)), rc=rc)

              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
            enddo
            deallocate(udimList)
          endif
        endif
        ! inquire if NetCDF file already contains the "time" variable
        ncerr = nf90_inq_varid(ncid, "time", varid=varid)
        if (ncerr /= NF90_NOERR) then
          ! the variable does not exist in the NetCDF file yet -> add it
          ! access the "time" attribute on the grid
          call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                 name="time", value=time, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                 name="time_iso", value=time_iso, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          ncerr = nf90_redef(ncid=ncid)
          if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          ncerr = nf90_inq_dimid(ncid, "time", dimid=dimid)
          if (ncerr /= NF90_NOERR) then
            ! "time" dimension does not yet exist, define as unlimited dim
            ncerr = nf90_def_dim(ncid, "time", NF90_UNLIMITED, dimid=dimid)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
          endif
          ncerr = nf90_def_var(ncid, "time", NF90_DOUBLE, &
                               dimids=(/dimid/), varid=varid)
          if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          ncerr = nf90_def_dim(ncid, "nchars", 20, ch_dimid)
          if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
          ncerr = nf90_def_var(ncid, "time_iso", NF90_CHAR, [ch_dimid,dimid], timeiso_varid)
          if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
          ncerr = nf90_put_att(ncid, timeiso_varid, "long_name", "valid time")
          if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
          ncerr = nf90_put_att(ncid, timeiso_varid, "description", "ISO 8601 datetime string")
          if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
          ncerr = nf90_put_att(ncid, timeiso_varid, "_Encoding", "UTF-8")
          if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          ncerr = nf90_enddef(ncid=ncid)
          if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          ncerr = nf90_put_var(ncid, varid, values=time)
          if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          ncerr = nf90_put_var(ncid, timeiso_varid, values=[trim(time_iso)])
          if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          ! loop over all the grid attributes that start with "time:", and
          ! put them on the "time" variable in the NetCDF file

          call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                 name="TimeAttributes", itemCount=itemCount, rc=rc)

          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

          if (itemCount > 0) then
            ncerr = nf90_redef(ncid=ncid)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

            allocate(attNameList(itemCount))
            call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                   name="TimeAttributes", valueList=attNameList, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            do i=1, itemCount
              attName = attNameList(i)
              call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                                     name=trim(attNameList(i)), typekind=typekind, rc=rc)

              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

              if (typekind==ESMF_TYPEKIND_CHARACTER) then
                call ESMF_AttributeGet(grid,                               &
                                       convention="NetCDF", purpose="FV3", &
                                       name=trim(attNameList(i)), value=valueS, rc=rc)

                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

                ncerr = nf90_put_att(ncid, varid, &
                                     trim(attName(6:len(attName))), values=valueS)

                if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

              else if (typekind==ESMF_TYPEKIND_I4) then

                call ESMF_AttributeGet(grid,                               &
                                       convention="NetCDF", purpose="FV3", &
                                       name=trim(attNameList(i)), value=valueI4, rc=rc)

                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                ncerr = nf90_put_att(ncid, varid, &
                                     trim(attName(6:len(attName))), values=valueI4)

                if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

              else if (typekind==ESMF_TYPEKIND_R4) then
                call ESMF_AttributeGet(grid,                               &
                                       convention="NetCDF", purpose="FV3", &
                                       name=trim(attNameList(i)), value=valueR4, rc=rc)

                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

                ncerr = nf90_put_att(ncid, varid, &
                                     trim(attName(6:len(attName))), values=valueR4)

                if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

              else if (typekind==ESMF_TYPEKIND_R8) then
                call ESMF_AttributeGet(grid,                               &
                                       convention="NetCDF", purpose="FV3", &
                                       name=trim(attNameList(i)), value=valueR8, rc=rc)

                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
                ncerr = nf90_put_att(ncid, varid, &
                                     trim(attName(6:len(attName))), values=valueR8)

                if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
              endif
            enddo
            deallocate(attNameList)
            ncerr = nf90_enddef(ncid=ncid)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
          endif
        endif
        ! close the NetCDF file
        ncerr = nf90_close(ncid=ncid)

        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
      call ESMF_LogWrite("In ioCompRun() after "// &
                         "writing vectical and time dimensions.", ESMF_LOGMSG_INFO, rc=rc)
    endif

    !TODO: remove this block once the ESMF_FieldBundleWrite() below allows to
    !TODO: specify the NetCDF 64-bit-offset file format.
    if (status==ESMF_FILESTATUS_REPLACE) then
      ! First time in with this filename (therefore 'replace'), create the
      ! file with 64bit-offset format in order to accommodate larger data
      ! volume.
      call ESMF_GridCompGet(comp, localPet=localPet, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      if (localPet==0) then
        ! only single PET to deal with NetCDF
        ncerr = nf90_create(tileFileName, &
                            cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid)

        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        ncerr = nf90_close(ncid=ncid)

        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
      status = ESMF_FILESTATUS_OLD  ! switch status to 'OLD' to not overwrite
      call ESMF_LogWrite("In ioCompRun() after creating the NetCDF file", &
                         ESMF_LOGMSG_INFO, rc=rc)
    endif

    if (timeslice==-1) then
      call ESMF_FieldBundleWrite(fieldbundle=wrtTileFB, fileName=tileFileName,          &
                                 convention=convention, purpose=purpose, status=status, &
                                 overwrite=.true., rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    else
      call ESMF_FieldBundleWrite(fieldbundle=wrtTileFB, fileName=tileFileName, &
                                 convention=convention, purpose=purpose, status=status, &
                                 timeslice=timeslice, overwrite=.true., rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif
    call ESMF_LogWrite("In ioCompRun() after "// &
                       "ESMF_FieldBundleWrite()", ESMF_LOGMSG_INFO, rc=rc)

  contains

    subroutine write_out_ungridded_dim_atts(dimLabel, rc)
      character(len=*)      :: dimLabel
      integer, intent(out)  :: rc

      logical               :: isPresent

      ! inquire if NetCDF file already contains this ungridded dimension
      ncerr = nf90_inq_varid(ncid, trim(dimLabel), varid=varid)
      if (ncerr == NF90_NOERR) return
      ! the variable does not exist in the NetCDF file yet -> add it
      ! access the undistributed dimension attribute on the grid
      call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                             name=trim(dimLabel), isPresent=isPresent, itemCount=valueCount, typekind=typekind, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (.not.isPresent) return ! nothing there to do

      if( typekind == ESMF_TYPEKIND_R4 ) then
        allocate(valueListr4(valueCount))
        call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                               name=trim(dimLabel), valueList=valueListr4, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      else if ( typekind == ESMF_TYPEKIND_R8) then
        allocate(valueListr8(valueCount))
        call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                              name=trim(dimLabel), valueList=valueListr8, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      else
        write(0,*) 'in write_out_ungridded_dim_atts: ERROR unknown typekind'
      endif
      ! now add it to the NetCDF file
      ncerr = nf90_redef(ncid=ncid)
      if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ncerr = nf90_inq_dimid(ncid, trim(dimLabel), dimid=dimid)
      if (ncerr /= NF90_NOERR) then
        ! dimension does not yet exist, and must be defined
        ncerr = nf90_def_dim(ncid, trim(dimLabel), valueCount, dimid=dimid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
      if( typekind == ESMF_TYPEKIND_R4 ) then
        ncerr = nf90_def_var(ncid, trim(dimLabel), NF90_FLOAT, &
                             dimids=(/dimid/), varid=varid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        ncerr = nf90_enddef(ncid=ncid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        ncerr = nf90_put_var(ncid, varid, values=valueListr4)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        deallocate(valueListr4)
      else if(typekind == ESMF_TYPEKIND_R8) then
        ncerr = nf90_def_var(ncid, trim(dimLabel), NF90_DOUBLE, &
                             dimids=(/dimid/), varid=varid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        ncerr = nf90_enddef(ncid=ncid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        ncerr = nf90_put_var(ncid, varid, values=valueListr8)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        deallocate(valueListr8)
      endif
      ! add attributes to this vertical variable
      call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3", &
                             attnestflag=ESMF_ATTNEST_OFF, count=attCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (attCount>0) then
        ncerr = nf90_redef(ncid=ncid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
      ! loop over all the attributes
      do j=1, attCount
        call ESMF_AttributeGet(grid, convention="NetCDF", purpose="FV3",       &
                               attnestflag=ESMF_ATTNEST_OFF, attributeIndex=j, &
                               name=attName, typekind=typekind, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        ! test for name starting with trim(dimLabel)":"
        if (index(trim(attName), trim(dimLabel)//":") == 1) then
          ind = len(trim(dimLabel)//":")
          ! found a matching attributes
          if (typekind == ESMF_TYPEKIND_CHARACTER) then
            call ESMF_AttributeGet(grid, &
                                   convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=valueS, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ncerr = nf90_put_att(ncid, varid, &
                                 trim(attName(ind+1:len(attName))), values=valueS)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          else if (typekind == ESMF_TYPEKIND_I4) then
            call ESMF_AttributeGet(grid, &
                                   convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=valueI4, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ncerr = nf90_put_att(ncid, varid, &
                                 trim(attName(ind+1:len(attName))), values=valueI4)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          else if (typekind == ESMF_TYPEKIND_R4) then
            call ESMF_AttributeGet(grid, &
                                   convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=valueR4, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ncerr = nf90_put_att(ncid, varid, &
                                 trim(attName(ind+1:len(attName))), values=valueR4)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          else if (typekind == ESMF_TYPEKIND_R8) then
            call ESMF_AttributeGet(grid, &
                                   convention="NetCDF", purpose="FV3", &
                                   name=trim(attName), value=valueR8, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ncerr = nf90_put_att(ncid, varid, &
                                 trim(attName(ind+1:len(attName))), values=valueR8)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
          endif
        endif
      enddo
      if (attCount>0) then
        ncerr = nf90_enddef(ncid=ncid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
    end subroutine write_out_ungridded_dim_atts

    subroutine write_out_ungridded_dim_atts_from_field(field, dimLabel, rc)

      type(ESMF_Field),intent(in) :: field
      character(len=*),intent(in) :: dimLabel
      integer, intent(out)  :: rc

      ! inquire if NetCDF file already contains this ungridded dimension
      ncerr = nf90_inq_varid(ncid, trim(dimLabel), varid=varid)
      if (ncerr == NF90_NOERR) return
      ! the variable does not exist in the NetCDF file yet -> add it
      ! access the undistributed dimension attribute on the grid
      call ESMF_AttributeGet(field, convention="NetCDF", purpose="FV3-dim", &
                             name=trim(dimLabel), itemCount=valueCount, typekind=typekind, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if( typekind == ESMF_TYPEKIND_R4 ) then
        allocate(valueListr4(valueCount))
        call ESMF_AttributeGet(field, convention="NetCDF", purpose="FV3-dim", &
                               name=trim(dimLabel), valueList=valueListr4, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      else if ( typekind == ESMF_TYPEKIND_R8) then
        allocate(valueListr8(valueCount))
        call ESMF_AttributeGet(field, convention="NetCDF", purpose="FV3-dim", &
                              name=trim(dimLabel), valueList=valueListr8, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      else
        write(0,*) 'in write_out_ungridded_dim_atts: ERROR unknown typekind'
      endif
      ! now add it to the NetCDF file
      ncerr = nf90_redef(ncid=ncid)
      if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

      ncerr = nf90_inq_dimid(ncid, trim(dimLabel), dimid=dimid)
      if (ncerr /= NF90_NOERR) then
        ! dimension does not yet exist, and must be defined
        ncerr = nf90_def_dim(ncid, trim(dimLabel), valueCount, dimid=dimid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
      if( typekind == ESMF_TYPEKIND_R4 ) then
        ncerr = nf90_def_var(ncid, trim(dimLabel), NF90_FLOAT, &
                             dimids=(/dimid/), varid=varid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        ncerr = nf90_enddef(ncid=ncid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        ncerr = nf90_put_var(ncid, varid, values=valueListr4)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        deallocate(valueListr4)
      else if(typekind == ESMF_TYPEKIND_R8) then
        ncerr = nf90_def_var(ncid, trim(dimLabel), NF90_DOUBLE, &
                             dimids=(/dimid/), varid=varid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        ncerr = nf90_enddef(ncid=ncid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

        ncerr = nf90_put_var(ncid, varid, values=valueListr8)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        deallocate(valueListr8)
      endif
      ! add attributes to this vertical variable
      call ESMF_AttributeGet(field, convention="NetCDF", purpose="FV3-dim", &
                             attnestflag=ESMF_ATTNEST_OFF, count=attCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (attCount>0) then
        ncerr = nf90_redef(ncid=ncid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
      ! loop over all the attributes
      do j=1, attCount
        call ESMF_AttributeGet(field, convention="NetCDF", purpose="FV3-dim",       &
                               attnestflag=ESMF_ATTNEST_OFF, attributeIndex=j, &
                               name=attName, typekind=typekind, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        ! test for name starting with trim(dimLabel)":"
        if (index(trim(attName), trim(dimLabel)//":") == 1) then
          ind = len(trim(dimLabel)//":")
          ! found a matching attributes
          if (typekind == ESMF_TYPEKIND_CHARACTER) then
            call ESMF_AttributeGet(field, &
                                   convention="NetCDF", purpose="FV3-dim", &
                                   name=trim(attName), value=valueS, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ncerr = nf90_put_att(ncid, varid, &
                                 trim(attName(ind+1:len(attName))), values=valueS)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          else if (typekind == ESMF_TYPEKIND_I4) then
            call ESMF_AttributeGet(field, &
                                   convention="NetCDF", purpose="FV3-dim", &
                                   name=trim(attName), value=valueI4, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ncerr = nf90_put_att(ncid, varid, &
                                 trim(attName(ind+1:len(attName))), values=valueI4)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          else if (typekind == ESMF_TYPEKIND_R4) then
            call ESMF_AttributeGet(field, &
                                   convention="NetCDF", purpose="FV3-dim", &
                                   name=trim(attName), value=valueR4, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ncerr = nf90_put_att(ncid, varid, &
                                 trim(attName(ind+1:len(attName))), values=valueR4)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

          else if (typekind == ESMF_TYPEKIND_R8) then
            call ESMF_AttributeGet(field, &
                                   convention="NetCDF", purpose="FV3-dim", &
                                   name=trim(attName), value=valueR8, rc=rc)

            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

            ncerr = nf90_put_att(ncid, varid, &
                                 trim(attName(ind+1:len(attName))), values=valueR8)

            if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
          endif
        endif
      enddo
      if (attCount>0) then
        ncerr = nf90_enddef(ncid=ncid)
        if (ESMF_LogFoundNetCDFError(ncerr, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      endif
    end subroutine

  end subroutine ioCompRun

  !-----------------------------------------------------------------------------

  subroutine ESMFproto_FieldMakeSingleTile(field, tile, tileField, petList, rc)
    type(ESMF_Field),     intent(in)              :: field
    integer,              intent(in)              :: tile
    type(ESMF_Field),     intent(out)             :: tileField
    integer, allocatable, intent(inout), optional :: petList(:)
    integer,              intent(out),   optional :: rc

    ! Take in a field on a multi-tile grid and return a field that only
    ! references a single tile.

    ! This routine only works with references, no data copies are being
    ! made. The single tile field that is returned points to the original
    ! field allocation.

    ! The original field passed in remains valid.

    type(ESMF_TypeKind_Flag)                :: typekind
    type(ESMF_StaggerLoc)                   :: staggerloc
    type(ESMF_Index_Flag)                   :: indexflag
    type(ESMF_Grid)                         :: grid, tileGrid
    type(ESMF_Array)                        :: array
    type(ESMF_DistGrid)                     :: distgrid
    type(ESMF_DELayout)                     :: delayout
    character(40)                           :: fieldName
    integer                                 :: fieldDimCount, gridDimCount
    integer                                 :: undistDims
    integer                                 :: localDeCount, deCount, tileCount
    integer,                   allocatable  :: gridToFieldMap(:)
    integer,                   allocatable  :: ungriddedLBound(:)
    integer,                   allocatable  :: ungriddedUBound(:)
    integer,                   allocatable  :: localDeToDeMap(:), deToTileMap(:)
    type(ESMF_LocalArray),     allocatable  :: lArrayList(:)
    integer                                 :: i
    integer                                 :: tileDeCount, tileLocalDeCount
    integer,                   allocatable  :: distgridToArrayMap(:)
    integer,                   allocatable  :: undistLBound(:), undistUBound(:)
    integer,                   allocatable  :: petMap(:), tilePetMap(:)
    integer,                   allocatable  :: minIndexPDe(:,:)
    integer,                   allocatable  :: maxIndexPDe(:,:)
    integer,                   allocatable  :: minIndexPTile(:,:)
    integer,                   allocatable  :: maxIndexPTile(:,:)
    integer,                   allocatable  :: deBlockList(:,:,:)

    character(800)                          :: msg

    if (present(rc)) rc = ESMF_SUCCESS

    ! access information from the incoming field
    call ESMF_FieldGet(field, array=array, typekind=typekind, &
                       staggerloc=staggerloc, &
                       dimCount=fieldDimCount, name=fieldName, grid=grid, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! access information from the associated grid
    call ESMF_GridGet(grid, dimCount=gridDimCount, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
#if 0
    write(msg,*) "fieldDimCount=",fieldDimCount,"gridDimCount=",gridDimCount
    call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)
#endif
    ! access list type information from the incoming field
    allocate(gridToFieldMap(gridDimCount))
    undistDims = fieldDimCount-gridDimCount
    if (undistDims < 0) undistDims = 0  ! this supports replicated dimensions
    allocate(ungriddedLBound(undistDims))
    allocate(ungriddedUBound(undistDims))
    call ESMF_FieldGet(field, gridToFieldMap=gridToFieldMap, &
                       ungriddedLBound=ungriddedLBound, ungriddedUBound=ungriddedUBound, &
                       rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! access information from associated array
    call ESMF_ArrayGet(array, distgrid=distgrid, delayout=delayout, &
                       indexflag=indexflag, localDeCount=localDeCount, deCount=deCount, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! access list type information from associated array
    allocate(localDeToDeMap(localDeCount), deToTileMap(deCount))
    allocate(distgridToArrayMap(gridDimCount))
    allocate(undistLBound(undistDims))
    allocate(undistUBound(undistDims))
    call ESMF_ArrayGet(array, tileCount=tileCount,                             &
                       localDeToDeMap=localDeToDeMap, deToTileMap=deToTileMap, &
                       distgridToArrayMap=distgridToArrayMap,                  &
                       undistLBound=undistLBound, undistUBound=undistUBound, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! access list type information from associated distgrid
    allocate(minIndexPDe(gridDimCount,deCount))
    allocate(maxIndexPDe(gridDimCount,deCount))
    allocate(minIndexPTile(gridDimCount,tileCount))
    allocate(maxIndexPTile(gridDimCount,tileCount))
    call ESMF_DistGridGet(distgrid,                                                 &
                          minIndexPDe=minIndexPDe, maxIndexPDe=maxIndexPDe,         &
                          minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, &
                          rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! access list type information from associated delayout
    allocate(petMap(deCount))
    call ESMF_DELayoutGet(delayout, petMap=petMap, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! construct data structures selecting specific tile
    allocate(lArrayList(localDeCount))
    tileLocalDeCount = 0
    do i=1, localDeCount
      if (deToTileMap(localDeToDeMap(i)+1) == tile) then
        ! localDe is on tile
        tileLocalDeCount = tileLocalDeCount + 1
        call ESMF_ArrayGet(array, localDe=i-1, &
                           localarray=lArrayList(tileLocalDeCount), rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      endif
    enddo
    allocate(tilePetMap(deCount))
    allocate(deBlockList(gridDimCount,2,deCount))
    tileDeCount=0
    do i=1, deCount
      if (deToTileMap(i) == tile) then
        ! DE is on tile
        tileDeCount = tileDeCount + 1
        tilePetMap(tileDeCount) = petMap(i)
        deBlockList(:,1,tileDeCount) = minIndexPDe(:,i)
        deBlockList(:,2,tileDeCount) = maxIndexPDe(:,i)
      endif
    enddo
    if (present(petList)) then
      if (.not.allocated(petList)) then
        allocate(petList(tileDeCount))
      else if (size(petList)/=tileDeCount) then
        deallocate(petList)
        allocate(petList(tileDeCount))
      endif
      petList(:) = tilePetMap(1:tileDeCount)
    endif
    ! create DELayout and DistGrid that only contain the single tile
    delayout = ESMF_DELayoutCreate(tilePetMap(1:tileDeCount), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    distgrid = ESMF_DistGridCreate(minIndex=minIndexPTile(:,tile),                    &
                                   maxIndex=maxIndexPTile(:,tile), delayout=delayout, &
                                   deBlockList=deBlockList(:,:,1:tileDeCount), rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! create an Array that only holds tile specific allocations
    if (tileLocalDeCount>0) then
      array = ESMF_ArrayCreate(distgrid, lArrayList(1:tileLocalDeCount),                         &
                               indexflag=indexflag,                                              &
                               distgridToArrayMap=distgridToArrayMap, undistLBound=undistLBound, &
                               undistUBound=undistUBound, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    else
      array = ESMF_ArrayCreate(distgrid, typekind, indexflag=indexflag,                          &
                               distgridToArrayMap=distgridToArrayMap, undistLBound=undistLBound, &
                               undistUBound=undistUBound, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif

    ! create a grid on the new distgrid
    tileGrid = ESMF_GridCreate(distgrid, indexflag=indexflag, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! alias the Attributes on grid level
    call ESMF_AttributeCopy(grid, tileGrid, attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! create the tile specific field from the array
    tileField = ESMF_FieldCreate(tileGrid, array=array, name=fieldName, &
                                 ! staggerloc=staggerloc, &
                                 gridToFieldMap=gridToFieldMap, ungriddedLBound=ungriddedLBound, &
                                 ungriddedUBound=ungriddedUBound, rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! alias the Attributes on field level
    call ESMF_AttributeCopy(field, tileField, attcopy=ESMF_ATTCOPY_REFERENCE, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    ! local garbage collection
    deallocate(localDeToDeMap, deToTileMap)
    deallocate(petMap,         tilePetMap)
    deallocate(minIndexPDe,    maxIndexPDe)
    deallocate(minIndexPTile,  maxIndexPTile)
    deallocate(gridToFieldMap, ungriddedLBound, ungriddedUBound)
    deallocate(deBlockList)

  end subroutine ESMFproto_FieldMakeSingleTile

!
!-----------------------------------------------------------------------
  subroutine splat4(idrt,jmax,aslat)

      implicit none
      integer,intent(in)  :: idrt,jmax
      real(4),intent(out) :: aslat(jmax)
!
      integer,parameter   :: KD=SELECTED_REAL_KIND(15,45)
      real(kind=KD)       :: pk(jmax/2),pkm1(jmax/2),pkm2(jmax/2)
      real(kind=KD)       :: aslatd(jmax/2),sp,spmax,eps=10.d0*epsilon(sp)
      integer,PARAMETER   :: JZ=50
      real(8) bz(jz)
      data bz        / 2.4048255577d0,  5.5200781103d0, &
       8.6537279129d0, 11.7915344391d0, 14.9309177086d0, 18.0710639679d0, &
      21.2116366299d0, 24.3524715308d0, 27.4934791320d0, 30.6346064684d0, &
      33.7758202136d0, 36.9170983537d0, 40.0584257646d0, 43.1997917132d0, &
      46.3411883717d0, 49.4826098974d0, 52.6240518411d0, 55.7655107550d0, &
      58.9069839261d0, 62.0484691902d0, 65.1899648002d0, 68.3314693299d0, &
      71.4729816036d0, 74.6145006437d0, 77.7560256304d0, 80.8975558711d0, &
      84.0390907769d0, 87.1806298436d0, 90.3221726372d0, 93.4637187819d0, &
      96.6052679510d0, 99.7468198587d0, 102.888374254d0, 106.029930916d0, &
      109.171489649d0, 112.313050280d0, 115.454612653d0, 118.596176630d0, &
      121.737742088d0, 124.879308913d0, 128.020877005d0, 131.162446275d0, &
      134.304016638d0, 137.445588020d0, 140.587160352d0, 143.728733573d0, &
      146.870307625d0, 150.011882457d0, 153.153458019d0, 156.295034268d0 /
      real(8)           :: dlt
      integer           :: jhe,jho
!     real(8),parameter :: PI=3.14159265358979d0,C=(1.d0-(2.d0/PI)**2)*0.25d0
      real(8),parameter ::                       C=(1.d0-(2.d0/PI)**2)*0.25d0
      real(8) r
      integer jh,n,j
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  GAUSSIAN LATITUDES
      IF(IDRT.EQ.4) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        R=1.d0/SQRT((JMAX+0.5d0)**2+C)
        DO J=1,MIN(JH,JZ)
          ASLATD(J)=COS(BZ(J)*R)
        ENDDO
        DO J=JZ+1,JH
          ASLATD(J)=COS((BZ(JZ)+(J-JZ)*PI)*R)
        ENDDO
        SPMAX=1.d0
        DO WHILE(SPMAX.GT.EPS)
          SPMAX=0.d0
          DO J=1,JH
            PKM1(J)=1.d0
            PK(J)=ASLATD(J)
          ENDDO
          DO N=2,JMAX
            DO J=1,JH
              PKM2(J)=PKM1(J)
              PKM1(J)=PK(J)
              PK(J)=((2*N-1)*ASLATD(J)*PKM1(J)-(N-1)*PKM2(J))/N
            ENDDO
          ENDDO
          DO J=1,JH
            SP=PK(J)*(1.d0-ASLATD(J)**2)/(JMAX*(PKM1(J)-ASLATD(J)*PK(J)))
            ASLATD(J)=ASLATD(J)-SP
            SPMAX=MAX(SPMAX,ABS(SP))
          ENDDO
        ENDDO
!
        DO J=1,JH
          ASLAT(J)=ASLATD(J)
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  EQUALLY-SPACED LATITUDES INCLUDING POLES
      ELSEIF(IDRT.EQ.0) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        JHO=JHE-1
        DLT=PI/(JMAX-1)
        ASLAT(1)=1.d0
        DO J=2,JH
          ASLAT(J)=COS((J-1)*DLT)
        ENDDO
!
        DO J=1,JH
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  EQUALLY-SPACED LATITUDES EXCLUDING POLES
      ELSEIF(IDRT.EQ.256) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        JHO=JHE
        DLT=PI/JMAX
        ASLAT(1)=1.d0
        DO J=1,JH
          ASLAT(J)=COS((J-0.5)*DLT)
        ENDDO

        DO J=1,JH
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     end subroutine splat4
!----------------------------------------------------------------------
     subroutine splat8(idrt,jmax,aslat)
!$$$
      implicit none
      integer,intent(in)  :: idrt,jmax
      real(8),intent(out) :: aslat(jmax)
!
      integer,parameter   :: KD=SELECTED_REAL_KIND(15,45)
      real(kind=KD)       :: pk(jmax/2),pkm1(jmax/2),pkm2(jmax/2)
      real(kind=KD)       :: aslatd(jmax/2),sp,spmax,eps=10.d0*epsilon(sp)
      integer,parameter   :: jz=50
      real(8) bz(jz)
      data bz        / 2.4048255577d0,  5.5200781103d0, &
       8.6537279129d0, 11.7915344391d0, 14.9309177086d0, 18.0710639679d0, &
      21.2116366299d0, 24.3524715308d0, 27.4934791320d0, 30.6346064684d0, &
      33.7758202136d0, 36.9170983537d0, 40.0584257646d0, 43.1997917132d0, &
      46.3411883717d0, 49.4826098974d0, 52.6240518411d0, 55.7655107550d0, &
      58.9069839261d0, 62.0484691902d0, 65.1899648002d0, 68.3314693299d0, &
      71.4729816036d0, 74.6145006437d0, 77.7560256304d0, 80.8975558711d0, &
      84.0390907769d0, 87.1806298436d0, 90.3221726372d0, 93.4637187819d0, &
      96.6052679510d0, 99.7468198587d0, 102.888374254d0, 106.029930916d0, &
      109.171489649d0, 112.313050280d0, 115.454612653d0, 118.596176630d0, &
      121.737742088d0, 124.879308913d0, 128.020877005d0, 131.162446275d0, &
      134.304016638d0, 137.445588020d0, 140.587160352d0, 143.728733573d0, &
      146.870307625d0, 150.011882457d0, 153.153458019d0, 156.295034268d0 /
      real(8)           :: dlt
      integer(4)        :: jhe,jho
!     real(8),parameter :: PI=3.14159265358979d0,C=(1.d0-(2.d0/PI)**2)*0.25d0
      real(8),parameter ::                       C=(1.d0-(2.d0/PI)**2)*0.25d0
      real(8) r
      integer jh,n,j
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  GAUSSIAN LATITUDES
      IF(IDRT.EQ.4) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        R=1.d0/SQRT((JMAX+0.5d0)**2+C)
        DO J=1,MIN(JH,JZ)
          ASLATD(J)=COS(BZ(J)*R)
        ENDDO
        DO J=JZ+1,JH
          ASLATD(J)=COS((BZ(JZ)+(J-JZ)*PI)*R)
        ENDDO
        SPMAX=1.d0
        DO WHILE(SPMAX.GT.EPS)
          SPMAX=0.d0
          DO J=1,JH
            PKM1(J)=1.d0
            PK(J)=ASLATD(J)
          ENDDO
          DO N=2,JMAX
            DO J=1,JH
              PKM2(J)=PKM1(J)
              PKM1(J)=PK(J)
              PK(J)=((2*N-1)*ASLATD(J)*PKM1(J)-(N-1)*PKM2(J))/N
            ENDDO
          ENDDO
          DO J=1,JH
            SP=PK(J)*(1.d0-ASLATD(J)**2)/(JMAX*(PKM1(J)-ASLATD(J)*PK(J)))
            ASLATD(J)=ASLATD(J)-SP
            SPMAX=MAX(SPMAX,ABS(SP))
          ENDDO
        ENDDO
!
        DO J=1,JH
          ASLAT(J)=ASLATD(J)
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  EQUALLY-SPACED LATITUDES INCLUDING POLES
      ELSEIF(IDRT.EQ.0) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        JHO=JHE-1
        DLT=PI/(JMAX-1)
        ASLAT(1)=1.d0
        DO J=2,JH
          ASLAT(J)=COS((J-1)*DLT)
        ENDDO
        DO J=1,JH
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  EQUALLY-SPACED LATITUDES EXCLUDING POLES
      ELSEIF(IDRT.EQ.256) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        JHO=JHE
        DLT=PI/JMAX
        ASLAT(1)=1.d0
        DO J=1,JH
          ASLAT(J)=COS((J-0.5d0)*DLT)
        ENDDO
!
        DO J=1,JH
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     end subroutine splat8
!
!
   subroutine rtll(tlmd,tphd,almd,aphd,tlm0d,tph0d)
!-------------------------------------------------------------------------------
      real(ESMF_KIND_R8), intent(in) :: tlmd, tphd
      real(ESMF_KIND_R8), intent(out) :: almd, aphd
      real(ESMF_KIND_R8), intent(in) :: tph0d, tlm0d
!-------------------------------------------------------------------------------
!     real(ESMF_KIND_R8), parameter :: pi=3.14159265358979323846
      real(ESMF_KIND_R8), parameter :: dtr=pi/180.0
!
      real(ESMF_KIND_R8) :: tph0, ctph0, stph0, tlm, tph, stph, ctph, ctlm, stlm, aph, cph
      real(ESMF_KIND_R8) :: xx, yy
!-------------------------------------------------------------------------------
!
      tph0=tph0d*dtr
      ctph0=cos(tph0)
      stph0=sin(tph0)
!
      tlm=tlmd*dtr
      tph=tphd*dtr
      stph=sin(tph)
      ctph=cos(tph)
      ctlm=cos(tlm)
      stlm=sin(tlm)
!
      xx=stph0*ctph*ctlm+ctph0*stph
      xx=max(xx,-1.0)
      xx=min(xx, 1.0)
      aph=asin(xx)
      cph=cos(aph)
!
      xx=(ctph0*ctph*ctlm-stph0*stph)/cph
      xx=max(xx,-1.0)
      xx=min(xx, 1.0)
      xx=acos(xx)/dtr
      yy=ctph*stlm/cph
      xx=sign(xx,yy)
      almd=tlm0d+xx

      aphd=aph/dtr
!
      if (almd > 180.0) then
         almd=almd-360.0
      end if
      if (almd < -180.0)  then
         almd=almd+360.0
      end if
!
      return
!
     end subroutine rtll
!
!-----------------------------------------------------------------------
!
     subroutine lambert(stlat1,stlat2,c_lat,c_lon,glon,glat,x,y,inv)

!-------------------------------------------------------------------------------
      real(ESMF_KIND_R8),      intent(in)    :: stlat1,stlat2,c_lat,c_lon
      real(ESMF_KIND_R8),      intent(inout) :: glon, glat
      real(ESMF_KIND_R8),      intent(inout) :: x, y
      integer,                 intent(in)    :: inv
!-------------------------------------------------------------------------------
!     real(ESMF_KIND_R8), parameter :: pi = 3.14159265358979323846
      real(ESMF_KIND_R8), parameter :: dtor = pi/180.0
      real(ESMF_KIND_R8), parameter :: rtod = 180.0/pi
      real(ESMF_KIND_R8), parameter :: a = 6371200.0
!-------------------------------------------------------------------------------
! inv ==  1    (glon,glat) ---> (x,y)
! inv == -1    (x,y) ---> (glon,glat)

      real(ESMF_KIND_R8) :: xp, yp, en, de, rho, rho0, rho2, dlon, theta, dr2
      real(ESMF_KIND_R8) :: h = 1.0

      ! For reference see:
      ! John P. Snyder (1987), Map projections: A working manual (pp. 104-110)
      ! https://doi.org/10.3133/pp1395

      if (stlat1 == stlat2) then
         en = sin(stlat1*dtor)
      else
         en = log(cos(stlat1*dtor)/cos(stlat2*dtor)) / &
              log(tan((45+0.5*stlat2)*dtor)/tan((45+0.5*stlat1)*dtor)) ! (15-3)
      endif
      h = sign(1.0_ESMF_KIND_R8,en)

      de = a*(cos(stlat1*dtor)*tan((45+0.5*stlat1)*dtor)**en)/en       ! (15-2)
      rho0 = de/(tan((45+0.5*c_lat)*dtor)**en)                         ! (15-1a)

      if (inv == 1) then          ! FORWARD TRANSFORMATION
         rho = de/(tan((45+0.5*glat)*dtor)**en)                        ! (15-1)
         dlon = modulo(glon-c_lon+180.0+3600.0,360.0)-180.0
         theta = en*dlon*dtor                                          ! (14-4)
         x = rho*sin(theta)                                            ! (14-1)
         y = rho0-rho*cos(theta)                                       ! (14-2)
      else if (inv == -1) then    ! INVERSE TRANSFORMATION
         xp = h*x;
         yp = h*(rho0-y)
         theta = atan2(xp,yp)                                          ! (14-11)
         glon = c_lon+(theta/en)*rtod                                  ! (14-9)
         glon = modulo(glon+180.0+3600.0,360.0)-180.0
         rho2 = xp*xp+yp*yp                                            ! (14-10)
         if (rho2 == 0.0) then
            glat = h*90.0
         else
            glat = 2.0*atan((de*de/rho2)**(1.0/(2.0*en)))*rtod-90.0    ! (15-5)
         endif
      else
        write (unit=*,fmt=*) " lambert: unknown inv argument"
        return
      end if

      return
     end subroutine lambert
!
!-----------------------------------------------------------------------
!
     subroutine get_outfile(nfl, filename, outfile_name, noutfile)
       integer, intent(in)          :: nfl
       character(*), intent(in)     :: filename(:,:)
       character(*), intent(inout)  :: outfile_name(:)
       integer, intent(inout)       :: noutfile

       integer        :: i,j,n
       logical        :: found
!
       noutfile = 0
       do i=1,nfl

         loopj: do j=1, 2000

           if( trim(filename(j,i)) == '') exit loopj
           if( trim(filename(j,i)) == 'none') cycle

           found = .false.
           loopn: do n=1, noutfile
             if(trim(filename(j,i)) == trim(outfile_name(n))) then
               found = .true.
               exit loopn
             endif
           enddo loopn

           if (.not.found) then
             noutfile = noutfile + 1
             outfile_name(noutfile) = trim(filename(j,i))
!             print *,'in get outfile,noutfile=', noutfile,' outfile=',trim(filename(j,i))
           endif
         enddo loopj
!
       enddo

     end subroutine get_outfile

     pure function trim_regridmethod_suffix(string) result(trimmed_string)
       character(len=*), intent(in) :: string
       character(len=:), allocatable :: trimmed_string

       trimmed_string = trim_suffix(trim(string),  '_bilinear')
       trimmed_string = trim_suffix(trimmed_string,'_patch')
       trimmed_string = trim_suffix(trimmed_string,'_nearest_stod')
       trimmed_string = trim_suffix(trimmed_string,'_nearest_dtos')
       trimmed_string = trim_suffix(trimmed_string,'_conserve')

     end function trim_regridmethod_suffix

     pure function trim_suffix(string, suffix) result(trimmed_string)
       character(len=*), intent(in) :: string, suffix
       character(len=:), allocatable :: trimmed_string
       integer :: suffix_length, string_length

       suffix_length = len(suffix)
       string_length = len(string)

       if (string_length >= suffix_length) then
         if (string(string_length-suffix_length+1:string_length) == suffix) then
           trimmed_string = string(1:string_length-suffix_length)
         else
           trimmed_string = string
         endif
       else
         trimmed_string = string
       endif

     end function trim_suffix
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
    end module  module_wrt_grid_comp
!
!-----------------------------------------------------------------------
