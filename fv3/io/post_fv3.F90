module post_fv3

  use mpi_f08

  use module_fv3_io_def,    only : wrttasks_per_group, fv3atm_output_dir, &
                                   lon1, lat1, lon2, lat2, dlon, dlat,    &
                                   cen_lon, cen_lat, dxin=>dx, dyin=>dy,  &
                                   stdlat1, stdlat2, output_grid
  use write_internal_state, only : wrt_internal_state

  implicit none

  public post_run_fv3

  contains

    subroutine post_run_fv3(wrt_int_state,grid_id,mype,mpicomp,lead_write, &
                            itasks,jtasks,mynfhr,mynfmin,mynfsec)
!
!  revision history:
!     Jul 2019    J. Wang             create interface to run inline post for FV3
!     Sep 2020    J. Dong/J. Wang     create interface to run inline post for FV3-LAM
!     Apr 2021    R. Sun              Added variables for Thomspon MP
!     Apr 2022    W. Meng             1)unify global and regional inline post interfaces
!                                     2)add bug fix for dx/dy computation
!                                     3)add reading pwat from FV3
!                                     4)remove some variable initializations
!                                     5)read max/min 2m T from tmax_max2m/tmin_min2m
!                                       for GFS, and from t02max/min for RRFS
!                                       and  HAFS.
!                                     6)read 3D cloud fraction from cld_amt for GFDL MP,
!                                       and from cldfra for other MPs.
!     Jun 2022    J. Meng             2D decomposition
!     Jul 2022    W. Meng             1)output lat/lon of four corner point for rotated
!                                       lat-lon grid.
!                                     2)read instant model top logwave
!
!-----------------------------------------------------------------------
!*** run post on write grid comp
!-----------------------------------------------------------------------
!
      use ctlblk_mod, only : komax,ifhr,ifmin,modelname,datapd,fld_info, &
                             npset,grib,jsta,  &
                             jend,ista,iend, im, nsoil, filenameflat,numx
      use gridspec_mod, only : maptype, gridtype,latstart,latlast,       &
                               lonstart,lonlast
      use grib2_module, only : gribit2,num_pset,nrecout,first_grbtbl
      use xml_perl_data,only : paramset
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(wrt_internal_state),intent(inout)    :: wrt_int_state
      integer,intent(in)                        :: grid_id
      integer,intent(in)                        :: mype
      type(MPI_Comm),intent(in)                 :: mpicomp
      integer,intent(in)                        :: lead_write
      integer,intent(in)                        :: itasks, jtasks
      integer,intent(in)                        :: mynfhr
      integer,intent(in)                        :: mynfmin
      integer,intent(in)                        :: mynfsec
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      integer              :: n,nwtpg,ierr,i,j,k,its,ite,jts,jte
      integer,allocatable  :: istagrp(:),iendgrp(:),jstagrp(:),jendgrp(:)
      integer,save         :: kpo,kth,kpv
      logical,save         :: first_run=.true.
      logical,save         :: read_postcntrl=.false.
      real(4),dimension(komax),save :: po, th, pv
      character(255)       :: post_fname
      integer,save         :: iostatusD3D=-1
!
!-----------------------------------------------------------------------
!*** set up dimensions
!-----------------------------------------------------------------------
!
      numx = itasks

      call post_getattr_fv3(wrt_int_state, grid_id)

      grib      = "grib2"
      gridtype  = "A"
      nsoil     = wrt_int_state%nsoil
      nwtpg     = wrt_int_state%petcount
      jts       = wrt_int_state%out_grid_info(grid_id)%j_start     !<-- Starting J of this write task's subsection
      jte       = wrt_int_state%out_grid_info(grid_id)%j_end       !<-- Ending J of this write task's subsection
      its       = wrt_int_state%out_grid_info(grid_id)%i_start     !<-- Starting I of this write task's subsection
      ite       = wrt_int_state%out_grid_info(grid_id)%i_end       !<-- Ending I of this write task's subsection

      if(mype==0) print *,'in post_run, numx=',numx,'its=',its,'ite=',ite,'nwtpg=',nwtpg, &
        'jts=',jts,'jte=',jte,'maptype=',maptype,'wrt_int_state%FBCount=',wrt_int_state%FBCount

!
!-----------------------------------------------------------------------
!*** set up fields to run post
!-----------------------------------------------------------------------
!
      if (allocated(jstagrp)) deallocate(jstagrp)
      if (allocated(jendgrp)) deallocate(jendgrp)
      if (allocated(istagrp)) deallocate(istagrp)
      if (allocated(iendgrp)) deallocate(iendgrp)
      allocate(jstagrp(nwtpg),jendgrp(nwtpg))
      allocate(istagrp(nwtpg),iendgrp(nwtpg))
!
      do n=0,nwtpg-1
        jstagrp(n+1) = wrt_int_state%out_grid_info(grid_id)%j_start_wrtgrp(n+1)
        jendgrp(n+1) = wrt_int_state%out_grid_info(grid_id)%j_end_wrtgrp  (n+1)
        istagrp(n+1) = wrt_int_state%out_grid_info(grid_id)%i_start_wrtgrp(n+1)
        iendgrp(n+1) = wrt_int_state%out_grid_info(grid_id)%i_end_wrtgrp  (n+1)
      enddo
      if(mype==0) print *,'in post_run,jstagrp=',jstagrp,'jendgrp=',jendgrp
      if(mype==0) print *,'in post_run,istagrp=',istagrp,'iendgrp=',iendgrp

!-----------------------------------------------------------------------
!*** read namelist for pv,th,po
!-----------------------------------------------------------------------
!
      call read_postnmlt(kpo,kth,kpv,po,th,pv,wrt_int_state%post_namelist)
!
!-----------------------------------------------------------------------
!*** allocate post variables
!-----------------------------------------------------------------------
!
      if(mype==0) print *,'in post_run,be post_alctvars, dim=',wrt_int_state%out_grid_info(grid_id)%im, &
        wrt_int_state%out_grid_info(grid_id)%jm, wrt_int_state%out_grid_info(grid_id)%lm,'mype=',mype,'wrttasks_per_group=', &
        wrttasks_per_group,'lead_write=',lead_write,'jts=',jts,'jte=',jte,   &
        'jstagrp=',jstagrp,'jendgrp=',jendgrp

      call post_alctvars(wrt_int_state%out_grid_info(grid_id)%im, &
                         wrt_int_state%out_grid_info(grid_id)%jm, &
                         wrt_int_state%out_grid_info(grid_id)%lm, &
                         mype,wrttasks_per_group,lead_write, &
                         mpicomp%mpi_val,jts,jte,jstagrp,jendgrp,its,ite,istagrp,iendgrp)
!
!-----------------------------------------------------------------------
!*** read namelist for pv,th,po
!-----------------------------------------------------------------------
!
      first_grbtbl = first_run
      read_postcntrl = .true.
!
!-----------------------------------------------------------------------
!*** fill post variables with values from forecast results
!-----------------------------------------------------------------------
!
      ifhr  = mynfhr
      ifmin = mynfmin
      if (ifhr == 0) ifmin = 0
      if (mype == 0) print *,'bf set_postvars,ifmin=',ifmin,'ifhr=',ifhr

      call set_postvars_fv3(wrt_int_state,grid_id,mype,mpicomp)

      if (read_postcntrl) then
        if (ifhr == 0) then
          filenameflat = 'postxconfig-NT_FH00.txt'
          call read_xml()
        else if(ifhr > 0) then
          filenameflat = 'postxconfig-NT.txt'
          if(associated(paramset)) then
            if(size(paramset)>0) then
              do i=1,size(paramset)
                if (associated(paramset(i)%param)) then
                  if (size(paramset(i)%param)>0) then
                    deallocate(paramset(i)%param)
                    nullify(paramset(i)%param)
                  endif
                endif
              enddo
            endif
            deallocate(paramset)
            nullify(paramset)
          endif
          num_pset = 0
          call read_xml()
          read_postcntrl = .false.
        endif
        if(mype==0) print *,'af read_xml,name=',trim(filenameflat),' ifhr=',ifhr,' num_pset=',num_pset
      endif
!
      do npset = 1, num_pset
        call set_outflds(kth,th,kpv,pv)
        if(allocated(datapd))deallocate(datapd)
        allocate(datapd(ite-its+1,jte-jts+1,nrecout+100))
!$omp parallel do default(none),private(i,j,k),shared(nrecout,jend,jsta,datapd,ista,iend)
        do k=1,nrecout+100
          do j=1,jend+1-jsta
            do i=1,iend+1-ista
              datapd(i,j,k) = 0.
            enddo
          enddo
        enddo
        call get_postfilename(post_fname)
        if (grid_id > 1) then
          write(post_fname, '(A,I2.2)') trim(post_fname)//".nest", grid_id
        endif
        post_fname = trim(fv3atm_output_dir)//trim(post_fname)

        if (mype==0) print *,'post_fname=',trim(post_fname)

        call process(kth,kpv,th(1:kth),pv(1:kpv),iostatusD3D)

        call mpi_barrier(mpicomp,ierr)
        call gribit2(post_fname)
        if(allocated(datapd))deallocate(datapd)
        if(allocated(fld_info))deallocate(fld_info)
      enddo

      if( first_run ) then
         first_run = .false.
      endif
      call post_finalize('grib2')

    end subroutine post_run_fv3
!
!-----------------------------------------------------------------------
!
    subroutine post_getattr_fv3(wrt_int_state,grid_id)
!
      use esmf
      use ctlblk_mod,           only: im, jm, mpi_comm_comp,gdsdegr,spval
      use gridspec_mod,         only: latstart, latlast, lonstart,    &
                                      lonlast, cenlon, cenlat, dxval, &
                                      dyval, truelat2, truelat1,psmapf, &
                                      lonstartv, lonlastv, cenlonv,     &
                                      latstartv, latlastv, cenlatv,     &
                                      latstart_r,latlast_r,lonstart_r,  &
                                      lonlast_r, STANDLON, maptype, gridtype, &
                                      latse,lonse,latnw,lonnw
!
      implicit none
!
      type(wrt_internal_state),intent(inout)    :: wrt_int_state
      integer, intent(in) :: grid_id
!
! local variable
      integer i,j,k,n,kz, attcount, nfb
      integer ni,naryi,nr4,nr8,rc
      integer aklen,varival
      real(4) varr4val
      real(8) varr8val
      character(80) attName, hydrostatics, fldname
      type(ESMF_TypeKind_Flag)           :: typekind
      real(4), dimension(:), allocatable :: ak4,bk4
      real(8), dimension(:), allocatable :: ak8,bk8
      type(ESMF_FieldBundle)             :: fldbundle
      character(128)                     :: wrtFBName
!
      spval = 9.99e20
! field bundle
      do nfb=1, wrt_int_state%FBcount
        fldbundle = wrt_int_state%wrtFB(nfb)

        call ESMF_FieldBundleGet(fldbundle, name=wrtFBName, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) return

        if (wrtFBName(1:8) == 'restart_') cycle
        if (wrtFBName(1:18) == 'cubed_sphere_grid_') cycle

! set grid spec:
!      if(mype==0) print*,'in post_getattr_lam,output_grid=',trim(output_grid(grid_id)),'nfb=',nfb
!      if(mype==0) print*,'in post_getattr_lam, lon1=',lon1,lon2,lat1,lat2,dlon,dlat
      gdsdegr = 1000000.

      if(trim(output_grid(grid_id)) == 'regional_latlon' .or. &
         trim(output_grid(grid_id)) == 'regional_latlon_moving') then
        MAPTYPE=0
        gridtype='A'

        if( lon1(grid_id)<0 ) then
          lonstart = nint((lon1(grid_id)+360.)*gdsdegr)
        else
          lonstart = nint(lon1(grid_id)*gdsdegr)
        endif
        if( lon2(grid_id)<0 ) then
          lonlast = nint((lon2(grid_id)+360.)*gdsdegr)
        else
          lonlast = nint(lon2(grid_id)*gdsdegr)
        endif
        latstart = nint(lat1(grid_id)*gdsdegr)
        latlast  = nint(lat2(grid_id)*gdsdegr)

        dxval = dlon(grid_id)*gdsdegr
        dyval = dlat(grid_id)*gdsdegr

!        if(mype==0) print*,'lonstart,latstart,dyval,dxval', &
!        lonstart,lonlast,latstart,latlast,dyval,dxval

      else if(trim(output_grid(grid_id)) == 'lambert_conformal') then
        MAPTYPE=1
        GRIDTYPE='A'

        if( cen_lon(grid_id)<0 ) then
          cenlon = nint((cen_lon(grid_id)+360.)*gdsdegr)
        else
          cenlon = nint(cen_lon(grid_id)*gdsdegr)
        endif
        cenlat = cen_lat(grid_id)*gdsdegr
        if( lon1(grid_id)<0 ) then
          lonstart = nint((lon1(grid_id)+360.)*gdsdegr)
        else
          lonstart = nint(lon1(grid_id)*gdsdegr)
        endif
        latstart = nint(lat1(grid_id)*gdsdegr)

        truelat1 = nint(stdlat1(grid_id)*gdsdegr)
        truelat2 = nint(stdlat2(grid_id)*gdsdegr)

        if(dxin(grid_id)<spval) then
          dxval = dxin(grid_id)*1.0e3
          dyval = dyin(grid_id)*1.0e3
        else
          dxval = spval
          dyval = spval
        endif

        STANDLON = cenlon
      else if(trim(output_grid(grid_id)) == 'rotated_latlon' .or. &
              trim(output_grid(grid_id)) == 'rotated_latlon_moving') then
        MAPTYPE=207
        GRIDTYPE='A'

        if( cen_lon(grid_id)<0 ) then
          cenlon = nint((cen_lon(grid_id)+360.)*gdsdegr)
        else
          cenlon = nint(cen_lon(grid_id)*gdsdegr)
        endif
        cenlat = cen_lat(grid_id)*gdsdegr
        if( lon1(grid_id)<0 ) then
          lonstart = nint((lon1(grid_id)+360.)*gdsdegr)
        else
          lonstart = nint(lon1(grid_id)*gdsdegr)
        endif
        if( lon2(grid_id)<0 ) then
          lonlast = nint((lon2(grid_id)+360.)*gdsdegr)
        else
          lonlast = nint(lon2(grid_id)*gdsdegr)
        endif
        latstart = nint(lat1(grid_id)*gdsdegr)
        latlast  = nint(lat2(grid_id)*gdsdegr)
        latstart_r = latstart
        lonstart_r = lonstart
        latlast_r = latlast
        lonlast_r = lonlast
        lonstart = nint(wrt_int_state%out_grid_info(grid_id)%lonstart*gdsdegr)
        latstart = nint(wrt_int_state%out_grid_info(grid_id)%latstart*gdsdegr)
        lonse    = nint(wrt_int_state%out_grid_info(grid_id)%lonse*gdsdegr)
        latse    = nint(wrt_int_state%out_grid_info(grid_id)%latse*gdsdegr)
        lonnw    = nint(wrt_int_state%out_grid_info(grid_id)%lonnw*gdsdegr)
        latnw    = nint(wrt_int_state%out_grid_info(grid_id)%latnw*gdsdegr)
        lonlast  = nint(wrt_int_state%out_grid_info(grid_id)%lonlast*gdsdegr)
        latlast  = nint(wrt_int_state%out_grid_info(grid_id)%latlast*gdsdegr)

        if(dlon(grid_id)<spval) then
          dxval = dlon(grid_id)*gdsdegr
          dyval = dlat(grid_id)*gdsdegr
        else
          dxval = spval
          dyval = spval
        endif

!        if(mype==0) print*,'rotated latlon,lonstart,latstart,cenlon,cenlat,dyval,dxval', &
!          lonstart_r,lonlast_r,latstart_r,latlast_r,cenlon,cenlat,dyval,dxval
      else if(trim(output_grid(grid_id)) == 'gaussian_grid') then
        MAPTYPE=4
        gridtype='A'

        if( lon1(grid_id)<0 ) then
          lonstart = nint((lon1(grid_id)+360.)*gdsdegr)
        else
          lonstart = nint(lon1(grid_id)*gdsdegr)
        endif
        if( lon2(grid_id)<0 ) then
          lonlast = nint((lon2(grid_id)+360.)*gdsdegr)
        else
          lonlast = nint(lon2(grid_id)*gdsdegr)
        endif
        latstart = nint(lat1(grid_id)*gdsdegr)
        latlast  = nint(lat2(grid_id)*gdsdegr)

        dxval = dlon(grid_id)*gdsdegr
        dyval = dlat(grid_id)*gdsdegr

      else if(trim(output_grid(grid_id)) == 'global_latlon') then
        MAPTYPE=0
        gridtype='A'

        if( lon1(grid_id)<0 ) then
          lonstart = nint((lon1(grid_id)+360.)*gdsdegr)
        else
          lonstart = nint(lon1(grid_id)*gdsdegr)
        endif
        if( lon2(grid_id)<0 ) then
          lonlast = nint((lon2(grid_id)+360.)*gdsdegr)
        else
          lonlast = nint(lon2(grid_id)*gdsdegr)
        endif
        latstart = nint(lat1(grid_id)*gdsdegr)
        latlast  = nint(lat2(grid_id)*gdsdegr)

        dxval = dlon(grid_id)*gdsdegr
        dyval = dlat(grid_id)*gdsdegr

      endif

! look at the field bundle attributes
      call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
        attnestflag=ESMF_ATTNEST_OFF, Count=attcount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,file=__FILE__))return  ! bail out
!
      aklen=0.
      do i=1, attCount

        call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
          attnestflag=ESMF_ATTNEST_OFF, attributeIndex=i, name=attName, &
          typekind=typekind, itemCount=n,  rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,file=__FILE__))return  ! bail out

        if (typekind==ESMF_TYPEKIND_I4 ) then
          if(n==1) then
            call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
              name=trim(attName), value=varival, rc=rc)
            if (trim(attName) == 'ncnsto') wrt_int_state%ntrac=varival
            if (trim(attName) == 'ncld')   wrt_int_state%ncld=varival
            if (trim(attName) == 'nsoil')  wrt_int_state%nsoil=varival
            if (trim(attName) == 'fhzero')  wrt_int_state%fhzero=varival
            if (trim(attName) == 'imp_physics') wrt_int_state%imp_physics=varival
            if (trim(attName) == 'landsfcmdl') wrt_int_state%landsfcmdl=varival
          endif
        else if (typekind==ESMF_TYPEKIND_R4) then
          if(n==1) then
            call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
              name=trim(attName), value=varr4val, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return  ! bail out
            if (trim(attName) == 'dtp') wrt_int_state%dtp=varr4val
            if (trim(attName) == 'fhzero') wrt_int_state%fhzero=varr4val
!            print *,'in post_fv3, fhzero=',wrt_int_state%fhzero
          else if(n>1) then
            if(trim(attName) =="ak") then
              if(allocated(wrt_int_state%ak)) deallocate(wrt_int_state%ak)
              allocate(wrt_int_state%ak(n))
              call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                name=trim(attName), valueList=wrt_int_state%ak, rc=rc)
              wrt_int_state%out_grid_info(grid_id)%lm = n-1
            else if(trim(attName) =="bk") then
              if(allocated(wrt_int_state%bk)) deallocate(wrt_int_state%bk)
              allocate(wrt_int_state%bk(n))
              call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
                name=trim(attName), valueList=wrt_int_state%bk, rc=rc)
            endif
          endif
        else if (typekind==ESMF_TYPEKIND_R8) then
          if(n==1) then
            call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
              name=trim(attName), value=varr8val, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, file=__FILE__)) return  ! bail out
            if (trim(attName) == 'dtp')   then
               wrt_int_state%dtp=varr8val
            endif
          else if(n>1) then
            if(trim(attName) =="ak") then
              if(allocated(wrt_int_state%ak)) deallocate(wrt_int_state%ak)
              allocate(wrt_int_state%ak(n))
              call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
              name=trim(attName), valueList=wrt_int_state%ak, rc=rc)
              wrt_int_state%out_grid_info(grid_id)%lm = n-1
            else if(trim(attName) =="bk") then
              if(allocated(wrt_int_state%bk)) deallocate(wrt_int_state%bk)
              allocate(wrt_int_state%bk(n))
              call ESMF_AttributeGet(fldbundle, convention="NetCDF", purpose="FV3", &
              name=trim(attName), valueList=wrt_int_state%bk, rc=rc)
            endif
            wrt_int_state%out_grid_info(grid_id)%lm = size(wrt_int_state%ak) - 1
          endif
        endif
!
      enddo
!
      enddo !end nfb
!
    end subroutine post_getattr_fv3
!
!-----------------------------------------------------------------------
!
    subroutine set_postvars_fv3(wrt_int_state,grid_id,mype,mpicomp)
!
!  revision history:
!     Jul 2019    J. Wang      Initial code
!     Apr 2022    W. Meng      Unify set_postvars_gfs and
!                               set_postvars_regional to set_postvars_fv3
!     Apr 2023    W. Meng      Sync RRFS and GFS changes from off-line post
!     Jun 2023    W. Meng      Remove duplicate initialization;
!                              relocate computation of aerosol fields
!
!-----------------------------------------------------------------------
!*** set up post fields from nmint_state
!-----------------------------------------------------------------------
!
      use esmf
      use vrbls4d,     only: dust, smoke, fv3dust, coarsepm, SALT, SUSO, SOOT, &
                             WASO,no3,nh4, PP25, PP10, ebb
      use vrbls3d,     only: t, q, uh, vh, wh, alpint, dpres, zint, zmid, o3,  &
                             qqr, qqs, cwm, qqi, qqw, qqg, qqh, omga, cfr, pmid, &
                             q2, rlwtt, rswtt, tcucn, tcucns, train, el_pbl,   &
                             pint, exch_h, ref_10cm, qqni, qqnr, qqnw, qqnwfa, &
                             qqnifa, effri, effrl, effrs, aextc55, taod5503d,  &
                             duem, dusd, dudp, duwt, dusv, ssem, sssd, ssdp,   &
                             sswt, sssv, bcem, bcsd, bcdp, bcwt, bcsv, ocem,   &
                             ocsd, ocdp, ocwt, ocsv, rhomid
      use vrbls2d,     only: f, pd, sigt4, fis, pblh, ustar, z0, ths, qs, twbs,&
                             qwbs, avgcprate, cprate, avgprec, prec, lspa, sno,&
                             cldefi, th10, q10, tshltr, pshltr, albase,        &
                             avgalbedo, avgtcdc, czen, czmean, mxsnal,landfrac,&
                             radot, cfrach, cfracl, cfracm, avgcfrach, qshltr, &
                             avgcfracl, avgcfracm, cnvcfr, islope, cmc, grnflx,&
                             vegfrc, acfrcv, ncfrcv, acfrst, ncfrst, ssroff,   &
                             bgroff, rlwin,                                    &
                             rlwtoa, cldwork, alwin, alwout, alwtoa, rswin,    &
                             rswinc, rswout, aswin, auvbin, auvbinc, aswout,   &
                             aswtoa, sfcshx, sfclhx, subshx, snopcx, sfcux,    &
                             sfcvx, sfcuvx, gtaux, gtauy, potevp, u10, v10,    &
                             smstav, smstot, ivgtyp, isltyp, sfcevp, sfcexc,   &
                             acsnow, acsnom, sst, thz0, qz0, uz0, vz0, ptop,   &
                             htop, pbot, hbot, ptopl, pbotl, ttopl, ptopm,     &
                             pbotm, ttopm, ptoph, pboth, pblcfr, ttoph, runoff,&
                             tecan, tetran, tedir, twa, sndepac,               &
                             maxtshltr, mintshltr, maxrhshltr, minrhshltr,     &
                             dzice, smcwlt, suntime, fieldcapa, htopd, hbotd,  &
                             htops, hbots, aswintoa, maxqshltr, minqshltr,     &
                             acond, sr, u10h, v10h, avgedir, avgecan,paha,pahi,&
                             avgetrans, avgesnow, avgprec_cont, avgcprate_cont,&
                             avisbeamswin, avisdiffswin, airbeamswin, airdiffswin, &
                             alwoutc, alwtoac, aswoutc, aswtoac, alwinc, aswinc,&
                             avgpotevp, snoavg, ti, si, cuppt, fdnsst,         &
                             w_up_max, w_dn_max, up_heli_max,up_heli_min,      &
                             up_heli_max03,up_heli_min03,rel_vort_max01,       &
                             rel_vort_max, rel_vort_maxhy1, refd_max,          &
                             refdm10c_max, u10max, v10max, wspd10max, sfcuxi,  &
                             sfcvxi, t10m, t10avg, psfcavg, akhsavg, akmsavg,  &
                             albedo, tg, prate_max, pwat, snow_acm, snow_bkt,  &
                             acgraup, graup_bucket, acfrain, frzrn_bucket,     &
                             ltg1_max, ltg2_max, ltg3_max, hwp, albedo,        &
                             aod550,du_aod550,ss_aod550,su_aod550,oc_aod550,   &
                             bc_aod550,maod,                                   &
                             dustpm10, dustcb, bccb, occb, sulfcb, sscb,       &
                             dustallcb, ssallcb, dustpm, sspm, pp25cb, pp10cb, &
                             no3cb, nh4cb, dusmass, ducmass, dusmass25,ducmass25, &
                             snownc, graupelnc, qrmax, hail_maxhailcast,       &
                             smoke_ave,dust_ave,coarsepm_ave,swddif,swddni,    &
                             xlaixy,wspd10umax,wspd10vmax
      use soil,        only: sldpth, sh2o, smc, stc, sllevel
      use masks,       only: lmv, lmh, htm, vtm, gdlat, gdlon, dx, dy, hbm2, sm, sice
      use ctlblk_mod,  only: im, jm, lm, lp1, jsta, jend, jsta_2l, jend_2u, jsta_m,jend_m, &
                             ista, iend, ista_2l, iend_2u, ista_m,iend_m,qmin, &
                             lsm, pt, imp_physics, spval, mpi_comm_comp, gdsdegr,  &
                             tprec, tclod, trdlw, trdsw, tsrfc, tmaxmin, theat, &
                             ardlw, ardsw, asrfc, avrain, avcnvc, iSF_SURFACE_PHYSICS,&
                             td3d, idat, sdat, ifhr, ifmin, dt, nphs, dtq2, pt_tbl, &
                             alsl, spl, ihrst, modelname, nsoil, rdaod, gocart_on,  &
                             gccpp_on, nasa_on, d2d_chem, nbin_ss, nbin_bc, nbin_oc,&
                             nbin_du,nbin_su, nbin_no3, nbin_nh4
      use params_mod,  only: erad, dtr, capa, p1000, small,h1, d608, pi, rd
      use gridspec_mod,only: latstart, latlast, lonstart, lonlast, cenlon, cenlat, &
                             dxval, dyval, truelat2, truelat1, psmapf, cenlat,     &
                             lonstartv, lonlastv, cenlonv, latstartv, latlastv,    &
                             cenlatv,latstart_r,latlast_r,lonstart_r,lonlast_r,    &
                             maptype, gridtype, STANDLON,latse,lonse,latnw,lonnw
      use lookup_mod,  only: thl, plq, ptbl, ttbl, rdq, rdth, rdp, rdthe, pl,   &
                             qs0, sqs, sthe, ttblq, rdpq, rdtheq, stheq, the0q, the0
      use physcons,    only: grav => con_g, fv => con_fvirt, rgas => con_rd,    &
                             eps => con_eps, epsm1 => con_epsm1
      use rqstfld_mod
!
!      use write_internal_state, only: wrt_internal_state
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(wrt_internal_state),intent(in) :: wrt_int_state
      integer,intent(in)                  :: grid_id
      integer,intent(in)                  :: mype
      type(MPI_Comm),intent(in)           :: mpicomp
!
!-----------------------------------------------------------------------
!
      integer i, ip1, j, l, k, n, iret, ibdl, rc, kstart, kend
      integer i1,i2,j1,j2,k1,k2
      integer fieldDimCount,gridDimCount,ncount_field,bundle_grid_id
      integer jdate(8)
      logical foundland, foundice, found, mvispresent
      integer totalLBound3d(3), totalUBound3d(3)
      real(4) rinc(5), fillvalue
      real(8) fillvalue8
      real    tlmh,RADI,TMP,ES,TV,RHOAIR,tem,tstart,dtp
      real, dimension(:),allocatable    :: ak5, bk5
      real(ESMF_KIND_R4),dimension(:,:),pointer    :: arrayr42d
      real(ESMF_KIND_R8),dimension(:,:),pointer    :: arrayr82d
      real(ESMF_KIND_R4),dimension(:,:,:),pointer  :: arrayr43d
      real(ESMF_KIND_R8),dimension(:,:,:),pointer  :: arrayr83d
      real,dimension(:),    allocatable :: slat,qstl
      real,external::FPVSNEW
      real,dimension(:,:),allocatable :: dummy, p2d, t2d, q2d,  qs2d,  &
                             cw2d, cfr2d, snacc_land, snacc_ice,       &
                             acsnom_land, acsnom_ice
      real,dimension(:,:,:),allocatable :: ext550
      character(len=80)              :: fieldname, wrtFBName, flatlon, &
                                        VarName
      type(ESMF_Grid)                :: wrtGrid
      type(ESMF_Field)               :: theField
      type(ESMF_Field), allocatable  :: fcstField(:)
      type(ESMF_TypeKind_Flag)       :: typekind
      type(ESMF_TypeKind_Flag)       :: attTypeKind

!
!-----------------------------------------------------------------------
!***  INTEGER SCALAR/1D HISTORY VARIABLES
!-----------------------------------------------------------------------
!
      imp_physics = wrt_int_state%imp_physics       !set GFS mp physics to 99 for Zhao scheme
      dtp         = wrt_int_state%dtp
      iSF_SURFACE_PHYSICS = wrt_int_state%landsfcmdl
      spval = 9.99e20
!
! nems gfs has zhour defined
      tprec   = wrt_int_state%fhzero
      tclod   = tprec
      trdlw   = tprec
      trdsw   = tprec
      tsrfc   = tprec
      tmaxmin = tprec
      td3d    = tprec
!      if(mype==0)print*,'MP_PHYSICS= ',imp_physics,'tprec=',tprec,'tclod=',tclod, &
!       'dtp=',dtp,'tmaxmin=',tmaxmin,'jsta=',jsta,jend,im,jm

!      write(6,*) 'maptype and gridtype is ', maptype,gridtype
!
!$omp parallel do default(shared),private(i,j)
      do j=jsta,jend
        do  i=ista,iend
          gdlat(i,j) = wrt_int_state%out_grid_info(grid_id)%latPtr(i,j)
          gdlon(i,j) = wrt_int_state%out_grid_info(grid_id)%lonPtr(i,j)
        enddo
      enddo

      call exch(gdlat)
      call exch(gdlon)

!$omp parallel do default(none),private(i,j,ip1), &
!$omp&  shared(jsta,jend_m,im,dx,gdlat,gdlon,dy,ista,iend_m,maptype,dxval,dyval,gdsdegr)
      do j = jsta, jend_m
        do i = ista, iend_m
          ip1 = i + 1
          !if (ip1 > im) ip1 = ip1 - im
          if(maptype==207)then
            dx(i,j)=erad*dxval*dtr/gdsdegr
            dy(i,j)=erad*dyval*dtr/gdsdegr
          else
            dx(i,j) = erad*cos(gdlat(i,j)*dtr)*(gdlon(ip1,j)-gdlon(i,j))*dtr
            dy(i,j) = erad*(gdlat(i,j+1)-gdlat(i,j))*dtr  ! like A*DPH
          endif
        end do
      end do
!
      if(.not. allocated(ak5)) allocate(ak5(lm+1),bk5(lm+1))
      do i=1,lm+1
        ak5(i) = wrt_int_state%ak(i)
        bk5(i) = wrt_int_state%bk(i)
      enddo

!$omp parallel do default(none) private(i,j) shared(jsta,jend,f,gdlat,ista,iend)
      do j=jsta,jend
        do i=ista,iend
          f(I,J) = 1.454441e-4*sin(gdlat(i,j)*dtr)   ! 2*omeg*sin(phi)
        end do
      end do
!
      pt    = ak5(1)

! GFS set up DT to compute accumulated fields, set it to one
      dtq2 = wrt_int_state%dtp
      nphs = 2.
      dt   = dtq2/nphs

      !Allocate for regional models only
      if(modelname=='FV3R') then
        allocate(ext550(ista:iend,jsta:jend,lm))

        do l=1,lm
          do j=jsta,jend
            do i=ista,iend
              ext550(i,j,l)=spval
            end do
          end do
        end do
      endif

        allocate(snacc_ice(ista:iend,jsta:jend))
        allocate(snacc_land(ista:iend,jsta:jend))
        allocate(acsnom_ice(ista:iend,jsta:jend))
        allocate(acsnom_land(ista:iend,jsta:jend))

        do j=jsta,jend
          do i=ista,iend
            snacc_ice(i,j)=spval
            snacc_land(i,j)=spval
            acsnom_ice(i,j)=spval
            acsnom_land(i,j)=spval
          end do
        end do

!
! GFS doesn not yet output soil layer thickness, assign SLDPTH to be the same as nam
      sldpth(1) = 0.10
      sldpth(2) = 0.3
      sldpth(3) = 0.6
      sldpth(4) = 1.0

! set ncfrcv to 1, ncfrst to 1
!$omp parallel do default(none),private(i,j),shared(jsta,jend,spval,ista,iend), &
!$omp& shared(ncfrcv,ncfrst)
      do j=jsta,jend
        do i=ista,iend
          ncfrcv(i,j) = 1.0
          ncfrst(i,j) = 1.0
        enddo
      enddo

! GFS incoming sfc longwave has been averaged over 6 hr bucket, set ARDLW to 1
      ardlw  = 1.0
! GFS incoming sfc longwave has been averaged, set ARDLW to 1
      ardsw = 1.0
! GFS surface flux has been averaged, set  ASRFC to 1
      asrfc = 1.0

! set avrain to 1
      avrain = 1.0
      avcnvc = 1.0
      theat  = 6.0 ! just in case GFS decides to output T tendency

!
! get inital date
      sdat(1)  = wrt_int_state%idate(2)   !month
      sdat(2)  = wrt_int_state%idate(3)   !day
      sdat(3)  = wrt_int_state%idate(1)   !year
      ihrst    = wrt_int_state%idate(4)   !hour

      idat(1)  = wrt_int_state%fdate(2)
      idat(2)  = wrt_int_state%fdate(3)
      idat(3)  = wrt_int_state%fdate(1)
      idat(4)  = wrt_int_state%fdate(4)
      idat(5)  = wrt_int_state%fdate(5)
!
!      if(mype==0) print *,'idat=',idat,'sdat=',sdat,'ihrst=',ihrst
!      CALL W3DIFDAT(JDATE,IDATE,0,RINC)
!
!      if(mype==0)print *,' rinc=',rinc
!      ifhr = nint(rinc(2)+rinc(1)*24.)
!      if(mype==0)print *,' ifhr=',ifhr
!      ifmin = nint(rinc(3))
!      if(ifhr /= nint(fhour))print*,'find wrong Grib file';stop
!      if(mype==0)print*,' in INITPOST ifhr ifmin =',ifhr,ifmin
!
      tstart = 0.
!
!-----------------------------------------------------------------------------
! get post fields
!-----------------------------------------------------------------------------
!
     foundland = .false.
     foundice = .false.

     get_lsmsk: do ibdl=1, wrt_int_state%FBCount

       call ESMF_FieldBundleGet(wrt_int_state%wrtFB(ibdl), name=wrtFBName, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,file=__FILE__)) return

       if (wrtFBName(1:8) == 'restart_') cycle
       if (wrtFBName(1:18) == 'cubed_sphere_grid_') cycle

       call ESMF_AttributeGet(wrt_int_state%wrtFB(ibdl), convention="NetCDF", purpose="FV3", &
                              name="grid_id", value=bundle_grid_id, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

       if (grid_id /= bundle_grid_id) cycle

! find lans sea mask
        found = .false.
        call ESMF_FieldBundleGet(wrt_int_state%wrtFB(ibdl),fieldName='land',isPresent=found, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
!        if(mype==0) print *,'ibdl=',ibdl,'land, found=',found
        if (found) then
          call ESMF_FieldBundleGet(wrt_int_state%wrtFB(ibdl),'land',field=theField, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return  ! bail out
          call ESMF_FieldGet(theField, localDe=0, farrayPtr=arrayr42d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return  ! bail out
          call ESMF_AttributeGet(theField, convention="NetCDF", purpose="FV3", &
                   name='_FillValue', value=fillvalue, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
!           print *,'in post_lam, get land field value,fillvalue=',fillvalue

          !$omp parallel do default(none),private(i,j),shared(jsta,jend,ista,iend,spval,arrayr42d,sm,fillValue)
          do j=jsta, jend
            do i=ista, iend
              if (arrayr42d(i,j) /= spval .and. abs(arrayr42d(i,j)-fillValue)>small ) then
                sm(i,j) = 1.- arrayr42d(i,j)
              else
                sm(i,j) = spval
              endif
            enddo
          enddo
          foundland = .true.
        endif

! find ice fraction
        found = .false.
        call ESMF_FieldBundleGet(wrt_int_state%wrtFB(ibdl),'icec',isPresent=found, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out
!        if(mype==0) print *,'ibdl=',ibdl,'ice, found=',found
        if (found) then
          call ESMF_FieldBundleGet(wrt_int_state%wrtFB(ibdl),'icec',field=theField, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return  ! bail out
          call ESMF_FieldGet(theField, localDe=0, farrayPtr=arrayr42d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return  ! bail out
          call ESMF_AttributeGet(theField, convention="NetCDF", purpose="FV3", &
                   name='_FillValue', value=fillvalue, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
!           if(mype==0) print *,'in post_lam, get icec  field value,fillvalue=',fillvalue

          !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sice,arrayr42d,sm,fillValue)
          do j=jsta, jend
            do i=ista, iend
              sice(i,j) = arrayr42d(i,j)
              if(abs(arrayr42d(i,j)-fillvalue)<small) sice(i,j) = spval
              if (sm(i,j) /= spval .and. sm(i,j) == 0.0) sice(i,j) = 0.0
            enddo
          enddo
          foundice = .true.
        endif

     enddo get_lsmsk
     if (.not.foundland .or. .not.foundice) then
       rc=999
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return  ! bail out
     endif
!     if(mype==0) print *,'after find sm and sice,imp_physics=',imp_physics,'wrt_int_state%FBCount=',wrt_int_state%FBCount
!
     file_loop_all: do ibdl=1, wrt_int_state%FBCount
!
! get grid dimension count
!       if(mype==0) print *,'in setvar, read field, ibdl=',ibdl,'idim=',   &
!         ista,iend,'jdim=',jsta,jend

       call ESMF_FieldBundleGet(wrt_int_state%wrtFB(ibdl), grid=wrtGrid,  &
         fieldCount=ncount_field, name=wrtFBName,rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

       if (wrtFBName(1:8) == 'restart_') cycle
       if (wrtFBName(1:18) == 'cubed_sphere_grid_') cycle

       call ESMF_AttributeGet(wrt_int_state%wrtFB(ibdl), convention="NetCDF", purpose="FV3", &
                              name="grid_id", value=bundle_grid_id, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

       if (grid_id /= bundle_grid_id) cycle

       call ESMF_GridGet(wrtgrid, dimCount=gridDimCount, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

!       if(mype==0) print *,'in setvar, allocate fcstField,ibdl=',ibdl,'count=',ncount_field,'wrtFBname=',trim(wrtFBName)
       allocate(fcstField(ncount_field))
       call ESMF_FieldBundleGet(wrt_int_state%wrtFB(ibdl),           &
         fieldList=fcstField, itemorderflag=ESMF_ITEMORDER_ADDORDER, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return  ! bail out

!       if(mype==0) print *,'in setvar, read field, ibdl=',ibdl, 'nfield=',ncount_field
       do n=1, ncount_field
!
          call ESMF_FieldGet(fcstField(n),typekind=typekind, name=fieldname, &
            dimCount=fieldDimCount,rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return  ! bail out
          if (index(trim(fieldname),"vector") >0) cycle
!
!** for 2D fields
          if (fieldDimCount == 2) then

            if (typekind == ESMF_TYPEKIND_R4) then
              call ESMF_FieldGet(fcstField(n), localDe=0, farrayPtr=arrayr42d, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
              call ESMF_AttributeGet(fcstField(n), convention="NetCDF", purpose="FV3", &
                   name='_FillValue', value=fillValue, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
!              fillValue=9.99E20
!              print *,'in post_lam, get field ',trim(fieldname),' fillValue=',fillValue

            else if (typekind == ESMF_TYPEKIND_R8) then
              call ESMF_FieldGet(fcstField(n), localDe=0, farrayPtr=arrayr82d, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
              call ESMF_AttributeGet(fcstField(n), convention="NetCDF", purpose="FV3", &
                   name='_FillValue', value=fillValue8, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
              fillValue = fillValue8
!              print *,'in post_lam, get field ',trim(fieldname),' fillvalue r82r4=',fillvalue

              allocate( arrayr42d(ista:iend,jsta:jend))
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,arrayr42d,arrayr82d)
              do j=jsta, jend
                do i=ista, iend
                  arrayr42d(i,j) = arrayr82d(i,j)
                enddo
              enddo
            endif

            ! Terrain height (*G later)
            if(trim(fieldname)=='hgtsfc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,fis,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  fis(i,j)=arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) fis(i,j)=spval
                enddo
              enddo
            endif

            ! Surface pressure
!            if(trim(fieldname)=='pressfc') then
!              !$omp parallel do private(i,j)
!              do j=jsta,jend
!                do i=ista, iend
!                  pint(i,j)=arrayr42d(i,j)
!                enddo
!              enddo
!            endif

            ! PBL height using nemsio
            if(trim(fieldname)=='hpbl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,pblh,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  pblh(i,j)=arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillValue) < small) pblh(i,j)=spval
                enddo
              enddo
            endif

            ! Lightning threat index 1
            if(trim(fieldname)=='ltg1_max') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,ltg1_max,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  ltg1_max(i,j)=arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillValue) < small) ltg1_max(i,j)=spval
                enddo
              enddo
            endif

            ! Lightning threat index 2
            if(trim(fieldname)=='ltg2_max') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,ltg2_max,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  ltg2_max(i,j)=arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillValue) < small) ltg2_max(i,j)=spval
                enddo
              enddo
            endif

            ! Lightning threat index 3
            if(trim(fieldname)=='ltg3_max') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,ltg3_max,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  ltg3_max(i,j)=arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillValue) < small) ltg3_max(i,j)=spval
                enddo
              enddo
            endif

            ! Maximum hail diameter (mm) since last output
            if(trim(fieldname)=='hailcast_dhail') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,hail_maxhailcast,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  hail_maxhailcast(i,j)=arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillValue) < small) hail_maxhailcast(i,j)=spval
                enddo
              enddo
            endif

            ! hourly wildfire potential
            if(trim(fieldname)=='hwp_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,hwp,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  hwp(i,j)=arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillValue) < small) hwp(i,j)=spval
                enddo
              enddo
            endif

            !hourly averaged smoke
            if(trim(fieldname)=='smoke_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,smoke_ave,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  smoke_ave(i,j)=arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillValue) < small) smoke_ave(i,j)=spval
                enddo
              enddo
            endif

            !hourly averaged dust
            if(trim(fieldname)=='dust_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,dust_ave,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  dust_ave(i,j)=arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillValue) < small) dust_ave(i,j)=spval
                enddo
              enddo
            endif

            !hourly averaged coarsepm
            if(trim(fieldname)=='coarsepm_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,coarsepm_ave,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  coarsepm_ave(i,j)=arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillValue) < small) coarsepm_ave(i,j)=spval
                enddo
              enddo
            endif

            ! frictional velocity
            if(trim(fieldname)=='fricv') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,ustar,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  ustar(i,j)=arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillValue) < small) ustar(i,j)=spval
                enddo
              enddo
            endif

            ! roughness length
            if(trim(fieldname)=='sfcr') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,z0,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  z0(i,j)=arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillValue) < small) z0(i,j)=spval
                enddo
              enddo
            endif

            ! sfc exchange coeff
            if(trim(fieldname)=='sfexc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,sfcexc,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  sfcexc(i,j)=arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillValue) < small) sfcexc(i,j)=spval
                enddo
              enddo
            endif

            ! aerodynamic conductance
            if(trim(fieldname)=='acond') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,acond,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  acond(i,j)=arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillValue) < small) acond(i,j)=spval
                enddo
              enddo
            endif

            ! surface albedo
            if(trim(fieldname)=='sfalb') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,albedo,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  albedo(i,j)=arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillValue) < small) albedo(i,j)=spval
                enddo
              enddo
            endif

            ! surface potential T
            if(trim(fieldname)=='tmpsfc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,arrayr42d,ths,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval .and. abs(arrayr42d(i,j)-fillValue) > small) then
                    ths(i,j) = arrayr42d(i,j)
                  else
                    ths(i,j) = spval
                  endif
                enddo
              enddo
            endif

            ! foundation temperature
            if(trim(fieldname)=='tref') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,arrayr42d,fdnsst)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval) then
                    fdnsst(i,j) = arrayr42d(i,j)
                  endif
                enddo
              enddo
            endif

            ! convective precip in m per physics time step
            if(trim(fieldname)=='cpratb_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,dtq2,arrayr42d,avgcprate,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval .and. abs(arrayr42d(i,j)-fillValue) > small) then
                    avgcprate(i,j) = arrayr42d(i,j) * (dtq2*0.001)
                  else
                    avgcprate(i,j) = spval
                  endif
                enddo
              enddo
            endif

            ! continuous bucket convective precip in m per physics time step
            if(trim(fieldname)=='cprat_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,dtq2,arrayr42d,avgcprate_cont,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval .and. abs(arrayr42d(i,j)-fillValue) > small) then
                    avgcprate_cont(i,j) = arrayr42d(i,j) * (dtq2*0.001)
                  else
                    avgcprate_cont(i,j) = spval
                  endif
                enddo
              enddo
            endif

            ! time averaged bucketed precip rate
            if(trim(fieldname)=='prateb_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,dtq2,arrayr42d,avgprec,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval .and. abs(arrayr42d(i,j)-fillValue) > small) then
                    avgprec(i,j) = arrayr42d(i,j) * (dtq2*0.001)
                  else
                    avgprec(i,j) = spval
                  endif
                enddo
              enddo
            endif

            ! time averaged continuous precip rate in m per physics time step
            if(trim(fieldname)=='prate_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,dtq2,arrayr42d,avgprec_cont,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval .and.  abs(arrayr42d(i,j)-fillValue) > small) then
                    avgprec_cont(i,j) = arrayr42d(i,j) * (dtq2*0.001)
                  else
                    avgprec_cont(i,j) = spval
                  endif
                enddo
              enddo
            endif

            ! precip rate in m per physics time step
            if(trim(fieldname)=='tprcp') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,dtq2,dtp,arrayr42d,prec,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval .and.  abs(arrayr42d(i,j)-fillValue) > small) then
                    prec(i,j) = arrayr42d(i,j) * (dtq2*0.001) * 1000./dtp
                  else
                    prec(i,j) = spval
                  endif
                enddo
              enddo
            endif

            ! convective precip rate in m per physics time step
            if(trim(fieldname)=='cnvprcp') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,dtq2,dtp,arrayr42d,cprate,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval .and.  abs(arrayr42d(i,j)-fillValue) > small) then
                    cprate(i,j) = max(0.,arrayr42d(i,j)) * (dtq2*0.001) * 1000./dtp
                  else
                    cprate(i,j) = 0.
                  endif
                enddo
              enddo
            endif

            !Accumulated snowfall
            if(trim(fieldname)=='tsnowp') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,snow_acm,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  snow_acm(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) snow_acm(i,j) = spval
                enddo
              enddo
            endif

            !Snowfall bucket
            if(trim(fieldname)=='tsnowpb') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,snow_bkt,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  snow_bkt(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) snow_bkt(i,j) = spval
                enddo
              enddo
            endif

            !Accumulated graupel
            if(trim(fieldname)=='frozr') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,acgraup,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  acgraup(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) acgraup(i,j) = spval
                enddo
              enddo
            endif

            !Graupel bucket
            if(trim(fieldname)=='frozrb') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,graup_bucket,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  graup_bucket(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) graup_bucket(i,j) = spval
                enddo
              enddo
            endif

            !Accumulated freezing rain
            if(trim(fieldname)=='frzr') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,acfrain,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  acfrain(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) acfrain(i,j) = spval
                enddo
              enddo
            endif

            !Freezing rain bucket
            if(trim(fieldname)=='frzrb') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,frzrn_bucket,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  frzrn_bucket(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) frzrn_bucket(i,j) = spval
                enddo
              enddo
            endif

            !time step snow (in m)
            if(trim(fieldname)=='snow') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,snownc,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  snownc(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) snownc(i,j) = spval
                enddo
              enddo
            endif

            !time step graupel (in m)
            if(trim(fieldname)=='graupel') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,graupelnc,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  graupelnc(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) graupelnc(i,j) = spval
                enddo
              enddo
            endif

            ! max hourly surface precipitation rate
            if(trim(fieldname)=='pratemax') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,prate_max,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  prate_max(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) prate_max(i,j) = spval
                enddo
              enddo
            endif

            ! max hourly 1-km agl reflectivity
            if(trim(fieldname)=='refdmax') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,refd_max,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  refd_max(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) refd_max(i,j) = spval
                enddo
              enddo
            endif

            ! max hourly -10C reflectivity
            if(trim(fieldname)=='refdmax263k') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,refdm10c_max,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  refdm10c_max(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) refdm10c_max(i,j) = spval
                enddo
              enddo
            endif

            ! max hourly u comp of 10m agl wind
            if(trim(fieldname)=='u10max') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,u10max,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  u10max(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) u10max(i,j) = spval
                enddo
              enddo
            endif

            ! max hourly v comp of 10m agl wind
            if(trim(fieldname)=='v10max') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,v10max,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  v10max(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) v10max(i,j) = spval
                enddo
              enddo
            endif

            ! max temporal 10m agl wind speed
            if (modelname =='GFS')then
            if(trim(fieldname)=='wind10m_max') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,wspd10max,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  wspd10max(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) wspd10max(i,j) = spval
                enddo
              enddo
            endif
            else
            ! max hourly 10m agl wind speed
            if(trim(fieldname)=='spd10max') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,wspd10max,arrayr42d,sm,fillValue)
              do j=jsta,jend 
                do i=ista, iend
                  wspd10max(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) wspd10max(i,j) = spval
                enddo        
              enddo          
            endif  
            endif !end modelname

            ! u comp of temporal max 10m agl wind speed 
            if(trim(fieldname)=='u10m_max') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,wspd10umax,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  wspd10umax(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) wspd10umax(i,j) = spval
                enddo
              enddo
            endif

            ! v comp of temporal max 10m agl wind speed 
            if(trim(fieldname)=='v10m_max') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,wspd10vmax,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  wspd10vmax(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) wspd10vmax(i,j) = spval
                enddo
              enddo
            endif

            ! inst snow water eqivalent
            if(trim(fieldname)=='weasd') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sno,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  sno(i,j) = arrayr42d(i,j)
                  if (sm(i,j) == 1.0 .and. sice(i,j)==0.)sno(i,j) = spval
                  if (abs(arrayr42d(i,j)-fillValue) < small) sno(i,j) = spval
                enddo
              enddo
            endif

            ! ave snow cover
            if(trim(fieldname)=='snowc_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,snoavg,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  snoavg(i,j) = arrayr42d(i,j)
                  if (sm(i,j)==1.0 .and. sice(i,j)==0.) snoavg(i,j) = spval
                  if (abs(arrayr42d(i,j)-fillValue) < small) snoavg(i,j) = spval
                  if (snoavg(i,j) /= spval) snoavg(i,j) = snoavg(i,j)/100.
                enddo
              enddo
            endif

            ! snow depth in mm
            if(trim(fieldname)=='snod') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,si,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  si(i,j) = arrayr42d(i,j)
                  if (sm(i,j)==1.0 .and. sice(i,j)==0.) si(i,j)=spval
                  if (abs(arrayr42d(i,j)-fillValue) < small) si(i,j)=spval
                  if (si(i,j) /= spval) si(i,j) = si(i,j) * 1000.0
                enddo
              enddo
            endif

            ! 2m potential T (computed later)
            if(trim(fieldname)=='tmp2m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,tshltr,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  tshltr(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) tshltr(i,j) = spval
                enddo
              enddo
            endif

            ! surface potential T
            if(trim(fieldname)=='spfh2m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,qshltr,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  qshltr(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) qshltr(i,j) = spval
                enddo
              enddo
            endif

            ! mid day avg albedo in fraction
            if(trim(fieldname)=='albdo_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgalbedo,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  avgalbedo(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) avgalbedo(i,j) = spval
                  if (avgalbedo(i,j) /= spval) then
                    avgalbedo(i,j) = avgalbedo(i,j) * 0.01
                  endif
                enddo
              enddo
            endif

            ! time averaged column cloud fraction
            if(trim(fieldname)=='tcdc_aveclm') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgtcdc,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  avgtcdc(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) avgtcdc(i,j) = spval
                  if (avgtcdc(i,j) /= spval) then
                    avgtcdc(i,j) = avgtcdc(i,j) * 0.01
                  endif
                enddo
              enddo
            endif

            ! maximum snow albedo in fraction
            if(trim(fieldname)=='snoalb') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,mxsnal,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  mxsnal(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) mxsnal(i,j) = spval
                enddo
              enddo
            endif

            !  land fraction
            if(trim(fieldname)=='lfrac') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,landfrac,arrayr42d,sm)
              do j=jsta,jend
                do i=ista, iend
                  landfrac(i,j) = arrayr42d(i,j)
                  if (sm(i,j) /= 0.0) landfrac(i,j) = spval
                enddo
              enddo
            endif

            ! ave high cloud fraction
            if(trim(fieldname)=='tcdc_avehcl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgcfrach,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  avgcfrach(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) avgcfrach(i,j) = spval
                  if (avgcfrach(i,j) /= spval) then
                    avgcfrach(i,j) = avgcfrach(i,j) * 0.01
                  endif
                enddo
              enddo
            endif

            ! ave low cloud fraction
            if(trim(fieldname)=='tcdc_avelcl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgcfracl,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  avgcfracl(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) avgcfracl(i,j) = spval
                  if (avgcfracl(i,j) /= spval) then
                    avgcfracl(i,j) = avgcfracl(i,j) * 0.01
                  endif
                enddo
              enddo
            endif

            ! ave middle cloud fraction
            if(trim(fieldname)=='tcdc_avemcl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgcfracm,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  avgcfracm(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) avgcfracm(i,j) = spval
                  if (avgcfracm(i,j) /= spval) then
                    avgcfracm(i,j) = avgcfracm(i,j) * 0.01
                  endif
                enddo
              enddo
            endif

            ! inst convective cloud fraction
            if(trim(fieldname)=='tcdccnvcl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,cnvcfr,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  cnvcfr(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) cnvcfr(i,j) = spval
                  if (cnvcfr(i,j) /= spval) then
                    cnvcfr(i,j) = cnvcfr(i,j) * 0.01
                  endif
                enddo
              enddo
            endif

            ! slope type
            if(trim(fieldname)=='sltyp') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,arrayr42d,islope,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) < spval .and. abs(arrayr42d(i,j)-fillValue) > small) then
                    islope(i,j) = nint(arrayr42d(i,j))
                  else
                    islope(i,j) = 0
                  endif
                  if (abs(arrayr42d(i,j)-fillValue) < small) islope(i,j) = 0
                enddo
              enddo
            endif

            ! time averaged column cloud fraction
            if(trim(fieldname)=='cnwat') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,cmc,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  cmc(i,j) = arrayr42d(i,j)
                  if (abs(arrayr42d(i,j)-fillValue) < small) cmc(i,j) = spval
                  if (cmc(i,j) /= spval) cmc(i,j) = cmc(i,j) * 0.001
                  if (sm(i,j) /= 0.0) cmc(i,j) = spval
                enddo
              enddo
            endif

            ! frozen precip fraction
            if(trim(fieldname)=='cpofp') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,arrayr42d,sr,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) /= spval .and.  abs(arrayr42d(i,j)-fillValue) > small) then
                  !set range within (0,1)
                    sr(i,j) = min(1.,max(0.,arrayr42d(i,j)))
                  else
                    sr(i,j) = spval
                  endif
                enddo
              enddo
            endif

            ! sea ice skin temperature
            if(trim(fieldname)=='tisfc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sice,arrayr42d,ti,fillValue)
              do j=jsta,jend
                do i=ista,iend
                  if (arrayr42d(i,j) /= spval .and.  abs(arrayr42d(i,j)-fillValue) > small) then
                    ti(i,j) = arrayr42d(i,j)
                    if (sice(i,j) == spval .or. sice(i,j) == 0.) ti(i,j)=spval
                  else
                    ti(i,j) = spval
                  endif
                enddo
              enddo
            endif

            ! vegetation fraction
            if(trim(fieldname)=='veg') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,vegfrc,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  vegfrc(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) vegfrc(i,j)=spval
                  if (vegfrc(i,j) /= spval) then
                    vegfrc(i,j) = vegfrc(i,j) * 0.01
                  else
                    vegfrc(i,j) = 0.0
                  endif
                  if (sm(i,j) /= 0.0) vegfrc(i,j) = spval
                enddo
              enddo
            endif

            !assign soil depths for RUC LSM, hard wire 9 soil depths here
            !so they aren't missing.
            if (nsoil==9) then
              sllevel(1) = 0.0
              sllevel(2) = 0.01
              sllevel(3) = 0.04
              sllevel(4) = 0.1
              sllevel(5) = 0.3
              sllevel(6) = 0.6
              sllevel(7) = 1.0
              sllevel(8) = 1.6
              sllevel(9) = 3.0
            endif

            ! liquid volumetric soil mpisture in fraction
            if(trim(fieldname)=='soill1') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sh2o,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  sh2o(i,j,1) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sh2o(i,j,1) = spval
                  if (sm(i,j) /= 0.0) sh2o(i,j,1) = spval
                enddo
              enddo
            endif

            ! liquid volumetric soil mpisture in fraction
            if(trim(fieldname)=='soill2') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sh2o,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  sh2o(i,j,2) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sh2o(i,j,2) = spval
                  if (sm(i,j) /= 0.0) sh2o(i,j,2) = spval
                enddo
              enddo
            endif

            ! liquid volumetric soil mpisture in fraction
            if(trim(fieldname)=='soill3') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sh2o,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  sh2o(i,j,3) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sh2o(i,j,3) = spval
                  if (sm(i,j) /= 0.0) sh2o(i,j,3) = spval
                enddo
              enddo
            endif

            ! liquid volumetric soil mpisture in fraction
            if(trim(fieldname)=='soill4') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sh2o,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  sh2o(i,j,4) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sh2o(i,j,4) = spval
                  if (sm(i,j) /= 0.0) sh2o(i,j,4) = spval
                enddo
              enddo
            endif

            if(nsoil==9) then
            ! liquid volumetric soil mpisture in fraction
            if(trim(fieldname)=='soill5') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sh2o,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  sh2o(i,j,5) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sh2o(i,j,5) = spval
                  if (sm(i,j) /= 0.0) sh2o(i,j,5) = spval
                enddo
              enddo
            endif

            ! liquid volumetric soil mpisture in fraction
            if(trim(fieldname)=='soill6') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sh2o,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  sh2o(i,j,6) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sh2o(i,j,6) = spval
                  if (sm(i,j) /= 0.0) sh2o(i,j,6) = spval
                enddo
              enddo
            endif

            ! liquid volumetric soil mpisture in fraction
            if(trim(fieldname)=='soill7') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sh2o,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  sh2o(i,j,7) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sh2o(i,j,7) = spval
                  if (sm(i,j) /= 0.0) sh2o(i,j,7) = spval
                enddo
              enddo
            endif

            ! liquid volumetric soil mpisture in fraction
            if(trim(fieldname)=='soill8') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sh2o,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  sh2o(i,j,8) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sh2o(i,j,8) = spval
                  if (sm(i,j) /= 0.0) sh2o(i,j,8) = spval
                enddo
              enddo
            endif

            ! liquid volumetric soil mpisture in fraction
            if(trim(fieldname)=='soill9') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sh2o,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  sh2o(i,j,9) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sh2o(i,j,9) = spval
                  if (sm(i,j) /= 0.0) sh2o(i,j,9) = spval
                enddo
              enddo
            endif

            endif !nsoil

            ! volumetric soil moisture
            if(trim(fieldname)=='soilw1') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smc,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  smc(i,j,1) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) smc(i,j,1) = spval
                  if (sm(i,j) /= 0.0) smc(i,j,1) = spval
                enddo
              enddo
            endif

            ! volumetric soil moisture
            if(trim(fieldname)=='soilw2') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smc,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  smc(i,j,2) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) smc(i,j,2) = spval
                  if (sm(i,j) /= 0.0) smc(i,j,2) = spval
                enddo
              enddo
            endif

            ! volumetric soil moisture
            if(trim(fieldname)=='soilw3') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smc,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  smc(i,j,3) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) smc(i,j,3) = spval
                  if (sm(i,j) /= 0.0) smc(i,j,3) = spval
                enddo
              enddo
            endif

            ! volumetric soil moisture
            if(trim(fieldname)=='soilw4') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smc,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  smc(i,j,4) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) smc(i,j,4) = spval
                  if (sm(i,j) /= 0.0) smc(i,j,4) = spval
                enddo
              enddo
            endif

            if(nsoil==9) then
            ! volumetric soil moisture
            if(trim(fieldname)=='soilw5') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smc,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  smc(i,j,5) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) smc(i,j,5) = spval
                  if (sm(i,j) /= 0.0) smc(i,j,5) = spval
                enddo
              enddo
            endif

            if(trim(fieldname)=='soilw6') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smc,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  smc(i,j,6) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) smc(i,j,6) = spval
                  if (sm(i,j) /= 0.0) smc(i,j,6) = spval
                enddo
              enddo
            endif

            if(trim(fieldname)=='soilw7') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smc,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  smc(i,j,7) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) smc(i,j,7) = spval
                  if (sm(i,j) /= 0.0) smc(i,j,7) = spval
                enddo
              enddo
            endif

            if(trim(fieldname)=='soilw8') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smc,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  smc(i,j,8) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) smc(i,j,8) = spval
                  if (sm(i,j) /= 0.0) smc(i,j,8) = spval
                enddo
              enddo
            endif

            if(trim(fieldname)=='soilw9') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smc,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  smc(i,j,9) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) smc(i,j,9) = spval
                  if (sm(i,j) /= 0.0) smc(i,j,9) = spval
                enddo
              enddo
            endif

            endif !nsoil

            ! soil temperature
            if(trim(fieldname)=='soilt1') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,stc,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  stc(i,j,1) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) stc(i,j,1) = spval
                  !mask open water areas, combine with sea ice tmp
                  if (sm(i,j) /= 0.0 .and. sice(i,j) ==0.) stc(i,j,1) = spval
                enddo
              enddo
            endif

            ! soil temperature
            if(trim(fieldname)=='soilt2') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,stc,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  stc(i,j,2) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) stc(i,j,2) = spval
                  !mask open water areas, combine with sea ice tmp
                  if (sm(i,j) /= 0.0 .and. sice(i,j) ==0.) stc(i,j,2) = spval
                enddo
              enddo
            endif

            ! soil temperature
            if(trim(fieldname)=='soilt3') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,stc,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  stc(i,j,3) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) stc(i,j,3) = spval
                  !mask open water areas, combine with sea ice tmp
                  if (sm(i,j) /= 0.0 .and. sice(i,j) ==0.) stc(i,j,3) = spval
                enddo
              enddo
            endif

            ! soil temperature
            if(trim(fieldname)=='soilt4') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,stc,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  stc(i,j,4) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) stc(i,j,4) = spval
                  !mask open water areas, combine with sea ice tmp
                  if (sm(i,j) /= 0.0 .and. sice(i,j) ==0.) stc(i,j,4) = spval
                enddo
              enddo
            endif

            if(nsoil==9) then

            ! soil temperature
            if(trim(fieldname)=='soilt5') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,stc,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  stc(i,j,5) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) stc(i,j,5) = spval
                  !mask open water areas, combine with sea ice tmp
                  if (sm(i,j) /= 0.0 .and. sice(i,j) ==0.) stc(i,j,5) = spval
                enddo
              enddo
            endif

            ! soil temperature
            if(trim(fieldname)=='soilt6') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,stc,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  stc(i,j,6) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) stc(i,j,6) = spval
                  !mask open water areas, combine with sea ice tmp
                  if (sm(i,j) /= 0.0 .and. sice(i,j) ==0.) stc(i,j,6) = spval
                enddo
              enddo
            endif

            ! soil temperature
            if(trim(fieldname)=='soilt7') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,stc,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  stc(i,j,7) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) stc(i,j,7) = spval
                  !mask open water areas, combine with sea ice tmp
                  if (sm(i,j) /= 0.0 .and. sice(i,j) ==0.) stc(i,j,7) = spval
                enddo
              enddo
            endif

            ! soil temperature
            if(trim(fieldname)=='soilt8') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,stc,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  stc(i,j,8) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) stc(i,j,8) = spval
                  !mask open water areas, combine with sea ice tmp
                  if (sm(i,j) /= 0.0 .and. sice(i,j) ==0.) stc(i,j,8) = spval
                enddo
              enddo
            endif

            ! soil temperature
            if(trim(fieldname)=='soilt9') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,stc,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  stc(i,j,9) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) stc(i,j,9) = spval
                  !mask open water areas, combine with sea ice tmp
                  if (sm(i,j) /= 0.0 .and. sice(i,j) ==0.) stc(i,j,9) = spval
                enddo
              enddo
            endif

            endif !nsoil

            ! time averaged incoming sfc longwave
            if(trim(fieldname)=='dlwrf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,alwin,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  alwin(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) alwin(i,j) = spval
                enddo
              enddo
            endif

            ! inst incoming sfc longwave
            if(trim(fieldname)=='dlwrf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,rlwin,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  rlwin(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) rlwin(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged outgoing sfc longwave, CLDRAD puts a minus sign
            if(trim(fieldname)=='ulwrf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,alwout,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  alwout(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) alwout(i,j) = spval
                  if (alwout(i,j) /= spval) alwout(i,j) = -alwout(i,j)
                enddo
              enddo
            endif

            ! inst outgoing sfc longwave
            if(trim(fieldname)=='ulwrf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,radot,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  radot(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) radot(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged outgoing model top longwave
            if(trim(fieldname)=='ulwrf_avetoa') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,alwtoa,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  alwtoa(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) alwtoa(i,j) = spval
                enddo
              enddo
            endif

            ! outgoing model top logwave
            if(trim(fieldname)=='ulwrf_toa') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,rlwtoa,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  rlwtoa(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue)<  small) rlwtoa(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged incoming sfc shortwave
            if(trim(fieldname)=='dswrf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,aswin,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  aswin(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) aswin(i,j) = spval
                enddo
              enddo
            endif

            ! inst incoming sfc shortwave
            if(trim(fieldname)=='dswrf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,rswin,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  rswin(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) rswin(i,j) = spval
                enddo
              enddo
            endif

            ! inst incoming clear sky sfc shortwave
            if(trim(fieldname)=='dswrf_clr') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,rswinc,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  rswinc(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) rswinc(i,j) = spval
                enddo
              enddo
            endif

            ! inst incoming direct beam sfc shortwave
            if(trim(fieldname)=='visbmdi') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,swddni,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  swddni(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) swddni(i,j) = spval
                enddo
              enddo
            endif

            ! inst incoming diffuse sfc shortwave
            if(trim(fieldname)=='visdfdi') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,swddif,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  swddif(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) swddif(i,j) = spval
                enddo
              enddo
            endif

            ! leaf area index
            if(trim(fieldname)=='xlaixy') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,xlaixy,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  xlaixy(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) xlaixy(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged incoming sfc uv-b
            if(trim(fieldname)=='duvb_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,auvbin,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  auvbin(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) auvbin(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged incoming sfc clear sky uv-b
            if(trim(fieldname)=='cduvb_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,auvbinc,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  auvbinc(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) auvbinc(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged outgoing sfc shortwave,CLDRAD puts a minus sign
            if(trim(fieldname)=='uswrf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,aswout,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  aswout(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) aswout(i,j) = spval
                  if (aswout(i,j) /= spval) aswout(i,j) = -aswout(i,j)
                enddo
              enddo
            endif

            ! inst outgoing sfc shortwave
            if(trim(fieldname)=='uswrf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,rswout,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  rswout(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) rswout(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged model top incoming shortwave
            if(trim(fieldname)=='dswrf_avetoa') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,aswintoa,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  aswintoa(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) aswintoa(i,j) = spval
                enddo
              enddo
            endif

            ! ime averaged model top outgoing shortwave
            if(trim(fieldname)=='uswrf_avetoa') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,aswtoa,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  aswtoa(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) aswtoa(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged surface sensible heat flux, multiplied by -1 because
            ! wrf model fluxhas reversed sign convention using gfsio
            if(trim(fieldname)=='shtfl_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sfcshx,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  sfcshx(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sfcshx(i,j) = spval
                  if (sfcshx(i,j) /= spval) sfcshx(i,j) = -sfcshx(i,j)
                enddo
              enddo
            endif

            ! inst surface sensible heat flux
            if(trim(fieldname)=='shtfl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,twbs,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  twbs(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) twbs(i,j) = spval
                  if (twbs(i,j) /= spval) twbs(i,j) = -twbs(i,j)
                enddo
              enddo
            endif

            ! time averaged surface latent heat flux, multiplied by -1 because
            ! wrf model flux has reversed sign vonvention using gfsio
            if(trim(fieldname)=='lhtfl_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,sfclhx,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  sfclhx(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sfclhx(i,j) = spval
                  if (sfclhx(i,j) /= spval) sfclhx(i,j) = -sfclhx(i,j)
                enddo
              enddo
            endif

            ! inst surface latent heat flux
            if(trim(fieldname)=='lhtfl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,qwbs,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  qwbs(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) qwbs(i,j) = spval
                  if (qwbs(i,j) /= spval) qwbs(i,j) = -qwbs(i,j)
                enddo
              enddo
            endif

            ! time averaged ground heat flux
            if(trim(fieldname)=='gflux_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,subshx,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  subshx(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) subshx(i,j) = spval
                  if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) subshx(i,j) = spval
                enddo
              enddo
            endif

            ! inst ground heat flux
            if(trim(fieldname)=='gflux') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,grnflx,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  grnflx(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) grnflx(i,j) = spval
                  if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) grnflx(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged zonal momentum flux
            if(trim(fieldname)=='uflx_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,sfcux,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  sfcux(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sfcux(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged meridional momentum flux
            if(trim(fieldname)=='vflx_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,sfcvx,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  sfcvx(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sfcvx(i,j) = spval
                enddo
              enddo
            endif

            ! inst zonal momentum flux
            if(trim(fieldname)=='uflx') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,sfcuxi,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  sfcuxi(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sfcuxi(i,j) = spval
                enddo
              enddo
            endif

            ! inst meridional momentum flux
            if(trim(fieldname)=='vflx') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,sfcvxi,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  sfcvxi(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) sfcvxi(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged zonal gravity wave stress
            if(trim(fieldname)=='u-gwd_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,gtaux,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  gtaux(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) gtaux(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged meridional gravity wave stress
            if(trim(fieldname)=='v-gwd_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,gtauy,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  gtauy(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) gtauy(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged accumulated potential evaporation
            if(trim(fieldname)=='pevpr_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgpotevp,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  avgpotevp(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) avgpotevp(i,j) = spval
                  if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) avgpotevp(i,j) = spval
                enddo
              enddo
            endif

            ! inst potential evaporation
            if(trim(fieldname)=='pevpr') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,potevp,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  potevp(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) potevp(i,j) = spval
                  if (sm(i,j) == 1.0 .and. sice(i,j) ==0.) potevp(i,j) = spval
                enddo
              enddo
            endif

            ! 10 m u
            if(trim(fieldname)=='ugrd10m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,u10,arrayr42d,u10h,spval,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  u10(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) u10(i,j) = spval
                  u10h(i,j) = u10(i,j)
                enddo
              enddo
            endif

            ! 10 m v
            if(trim(fieldname)=='vgrd10m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,v10,arrayr42d,v10h,spval,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  v10(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) v10(i,j) = spval
                  v10h(i,j) = v10(i,j)
                enddo
              enddo
            endif

            ! vegetation type
            if(trim(fieldname)=='vtype') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,arrayr42d,ivgtyp,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) < spval) then
                    ivgtyp(i,j) = nint(arrayr42d(i,j))
                    if( abs(arrayr42d(i,j)-fillValue) < small)  ivgtyp(i,j) = 0
                  else
                    ivgtyp(i,j) = 0
                  endif
                enddo
              enddo
            endif

            ! soil type
            if(trim(fieldname)=='sotyp') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,arrayr42d,isltyp,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  if (arrayr42d(i,j) < spval) then
                    isltyp(i,j) = nint(arrayr42d(i,j))
                    if( abs(arrayr42d(i,j)-fillValue) < small)  isltyp(i,j) = 0
                  else
                    isltyp(i,j) = 0
                  endif
                enddo
              enddo
            endif

            ! wetness
            if(trim(fieldname)=='wetness') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,smstav,arrayr42d,fillvalue,spval)
              do j=jsta,jend
                do i=ista, iend
                  smstav(i,j) = arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillvalue)<small) smstav(i,j) = spval
                enddo
              enddo
            endif

            !sndepac
            if(trim(fieldname)=='snacc_land') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,snacc_land,arrayr42d,fillvalue,spval)
              do j=jsta,jend
                do i=ista, iend
                  snacc_land(i,j) = arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillvalue)<small) snacc_land(i,j) = spval
                enddo
              enddo
            endif
            if(trim(fieldname)=='snacc_ice') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,snacc_ice,arrayr42d,fillvalue,spval)
              do j=jsta,jend
                do i=ista, iend
                  snacc_ice(i,j) = arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillvalue)<small) snacc_ice(i,j) = spval
                enddo
              enddo
            endif

            !snom
            if(trim(fieldname)=='snom_land') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,acsnom_land,arrayr42d,fillvalue,spval)
              do j=jsta,jend
                do i=ista, iend
                  acsnom_land(i,j) = arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillvalue)<small) acsnom_land(i,j) = spval
                enddo
              enddo
            endif
            if(trim(fieldname)=='snom_ice') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,acsnom_ice,arrayr42d,fillvalue,spval)
              do j=jsta,jend
                do i=ista, iend
                  acsnom_ice(i,j) = arrayr42d(i,j)
                  if(abs(arrayr42d(i,j)-fillvalue)<small) acsnom_ice(i,j) = spval
                enddo
              enddo
            endif

            if(rdaod) then
              ! MERRA2 aerosols
              if(trim(fieldname)=='aod550') then
                !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,aod550,arrayr42d,fillValue)
                do j=jsta,jend
                  do i=ista, iend
                    aod550(i,j) = arrayr42d(i,j)
                    if(abs(arrayr42d(i,j)-fillvalue)<small) aod550(i,j) = spval
                  enddo
                enddo
              endif

              if(trim(fieldname)=='du_aod550') then
                !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,du_aod550,arrayr42d,fillValue)
                do j=jsta,jend
                  do i=ista, iend
                    du_aod550(i,j) = arrayr42d(i,j)
                    if(abs(arrayr42d(i,j)-fillvalue)<small) du_aod550(i,j) = spval
                  enddo
                enddo
              endif

              if(trim(fieldname)=='ss_aod550') then
                !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,ss_aod550,arrayr42d,fillValue)
                do j=jsta,jend
                  do i=ista, iend
                    ss_aod550(i,j) = arrayr42d(i,j)
                    if(abs(arrayr42d(i,j)-fillvalue)<small) ss_aod550(i,j) = spval
                  enddo
                enddo
              endif

              if(trim(fieldname)=='su_aod550') then
                !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,su_aod550,arrayr42d,fillValue)
                do j=jsta,jend
                  do i=ista, iend
                    su_aod550(i,j) = arrayr42d(i,j)
                    if(abs(arrayr42d(i,j)-fillvalue)<small) su_aod550(i,j) = spval
                  enddo
                enddo
              endif

              if(trim(fieldname)=='oc_aod550') then
                !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,oc_aod550,arrayr42d,fillValue)
                do j=jsta,jend
                  do i=ista, iend
                    oc_aod550(i,j) = arrayr42d(i,j)
                    if(abs(arrayr42d(i,j)-fillvalue)<small) oc_aod550(i,j) = spval
                  enddo
                enddo
              endif

              if(trim(fieldname)=='bc_aod550') then
                !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,bc_aod550,arrayr42d,fillValue)
                do j=jsta,jend
                  do i=ista, iend
                    bc_aod550(i,j) = arrayr42d(i,j)
                    if(abs(arrayr42d(i,j)-fillvalue)<small) bc_aod550(i,j) = spval
                  enddo
                enddo
              endif

            endif !end rdaod

            if ((gocart_on .or. gccpp_on) .and. d2d_chem) then

              do K = 1, nbin_du
                write(VarName, '(A,I3.3)') 'duem', k
                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,duem,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      duem(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) duem(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo

              do K = 1, nbin_du
                if ( K == 1) VarName='dust1sd'
                if ( K == 2) VarName='dust2sd'
                if ( K == 3) VarName='dust3sd'
                if ( K == 4) VarName='dust4sd'
                if ( K == 5) VarName='dust5sd'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,dusd,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      dusd(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) dusd(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo

              do K = 1, nbin_du
                if ( K == 1) VarName='dust1dp'
                if ( K == 2) VarName='dust2dp'
                if ( K == 3) VarName='dust3dp'
                if ( K == 4) VarName='dust4dp'
                if ( K == 5) VarName='dust5dp'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,dudp,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      dudp(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) dudp(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo

              do K = 1, nbin_du
                if ( K == 1) VarName='dust1wtl'
                if ( K == 2) VarName='dust2wtl'
                if ( K == 3) VarName='dust3wtl'
                if ( K == 4) VarName='dust4wtl'
                if ( K == 5) VarName='dust5wtl'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,duwt,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      duwt(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) duwt(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo

              do K = 1, nbin_du
                if ( K == 1) VarName='dust1wtc'
                if ( K == 2) VarName='dust2wtc'
                if ( K == 3) VarName='dust3wtc'
                if ( K == 4) VarName='dust4wtc'
                if ( K == 5) VarName='dust5wtc'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,dusv,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      dusv(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) dusv(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo

              do K = 1, nbin_ss
                if ( K == 1) VarName='ssem001'
                if ( K == 2) VarName='ssem002'
                if ( K == 3) VarName='ssem003'
                if ( K == 4) VarName='ssem004'
                if ( K == 5) VarName='ssem005'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,ssem,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      ssem(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) ssem(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo

              do K = 1, nbin_ss
                if ( K == 1) VarName='seas1sd'
                if ( K == 2) VarName='seas2sd'
                if ( K == 3) VarName='seas3sd'
                if ( K == 4) VarName='seas4sd'
                if ( K == 5) VarName='seas5sd'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,sssd,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      sssd(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) sssd(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo

              do K = 1, nbin_ss
                if ( K == 1) VarName='seas1dp'
                if ( K == 2) VarName='seas2dp'
                if ( K == 3) VarName='seas3dp'
                if ( K == 4) VarName='seas4dp'
                if ( K == 5) VarName='seas5dp'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,ssdp,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      ssdp(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) ssdp(i,j,K) = spval
                    enddo
                  enddo
                 endif
              enddo

              do K = 1, nbin_ss
                if ( K == 1) VarName='seas1wt'
                if ( K == 2) VarName='seas2wt'
                if ( K == 3) VarName='seas3wt'
                if ( K == 4) VarName='seas4wt'
                if ( K == 5) VarName='seas5wt'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,sswt,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      sswt(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) sswt(i,j,K) = spval
                    enddo
                  enddo
                 endif
              enddo

              do K = 1, nbin_ss
                if ( K == 1) VarName='seas1wtc'
                if ( K == 2) VarName='seas2wtc'
                if ( K == 3) VarName='seas3wtc'
                if ( K == 4) VarName='seas4wtc'
                if ( K == 5) VarName='seas5wtc'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,sssv,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      sssv(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) sssv(i,j,K) = spval
                    enddo
                  enddo
                 endif
              enddo

              do K = 1, nbin_bc
                if ( K == 1) VarName='bceman'
                if ( K == 2) VarName='bcembb'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,bcem,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      bcem(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) bcem(i,j,K) = spval
                    enddo
                  enddo
                 endif
              enddo

              do K = 1, nbin_bc
                if ( K == 1) VarName='bc1sd'
                if ( K == 2) VarName='bc2sd'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,bcsd,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      bcsd(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) bcsd(i,j,K) = spval
                    enddo
                  enddo
                 endif
              enddo

              do K = 1, nbin_bc
                if ( K == 1) VarName='bc1dp'
                if ( K == 2) VarName='bc2dp'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,bcdp,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      bcdp(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) bcdp(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo

              do K = 1, nbin_bc
                if ( K == 1) VarName='bc1wtl'
                if ( K == 2) VarName='bc2wtl'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,bcwt,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      bcwt(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) bcwt(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo

              do K = 1, nbin_bc
                if ( K == 1) VarName='bc1wtc'
                if ( K == 2) VarName='bc2wtc'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,bcsv,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      bcsv(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) bcsv(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo

              do K = 1, nbin_oc
                if ( K == 1) VarName='oceman'
                if ( K == 2) VarName='ocembb'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,ocem,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      ocem(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) ocem(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo

              do K = 1, nbin_oc
                if ( K == 1) VarName='oc1sd'
                if ( K == 2) VarName='oc2sd'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,ocsd,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      ocsd(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) ocsd(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo

              do K = 1, nbin_oc
                if ( K == 1) VarName='oc1dp'
                if ( K == 2) VarName='oc2dp'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,ocdp,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      ocdp(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) ocdp(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo

              do K = 1, nbin_oc
                if ( K == 1) VarName='oc1wtl'
                if ( K == 2) VarName='oc2wtl'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,ocwt,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      ocwt(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) ocwt(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo

              do K = 1, nbin_oc
                if ( K == 1) VarName='oc1wtc'
                if ( K == 2) VarName='oc2wtc'

                if(trim(fieldname)==VarName) then
                  !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,ocsv,arrayr42d,fillvalue)
                  do j=jsta,jend
                    do i=ista, iend
                      ocsv(i,j,K) = arrayr42d(i,j)
                      if( abs(arrayr42d(i,j)-fillValue) < small) ocsv(i,j,K) = spval
                    enddo
                  enddo
                endif
              enddo


              if(trim(fieldname)=='maod') then
                !$omp parallel do default(none) private(i,j,K) shared(jsta,jend,ista,iend,spval,maod,arrayr42d,fillvalue)
                do j=jsta,jend
                  do i=ista, iend
                    maod(i,j) = arrayr42d(i,j)
                    if( abs(arrayr42d(i,j)-fillValue) < small) maod(i,j) = spval
                  enddo
                enddo
              endif

            endif !end gocart_on

            ! inst  cloud top pressure
            if(trim(fieldname)=='prescnvclt') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,ptop,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  ptop(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  ptop(i,j) = spval
                  if(ptop(i,j) <= 0.0) ptop(i,j) = spval
                enddo
              enddo
            endif

            ! inst cloud bottom pressure
            if(trim(fieldname)=='prescnvclb') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,pbot,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  pbot(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) pbot(i,j) = spval
                  if(pbot(i,j) <= 0.0) pbot(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged low cloud top pressure
            if(trim(fieldname)=='pres_avelct') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,ptopl,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  ptopl(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) ptopl(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged low cloud bottom pressure
            if(trim(fieldname)=='pres_avelcb') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,pbotl,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  pbotl(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) pbotl(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged low cloud top temperature
            if(trim(fieldname)=='tmp_avelct') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,ttopl,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  ttopl(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) ttopl(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged middle cloud top pressure
            if(trim(fieldname)=='pres_avemct') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,ptopm,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  ptopm(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) ptopm(i,j) = spval
                enddo
              enddo
            endif
            ! time averaged middle cloud bottom pressure
            if(trim(fieldname)=='pres_avemcb') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,pbotm,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  pbotm(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) pbotm(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged middle cloud top temperature
            if(trim(fieldname)=='tmp_avemct') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,ttopm,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  ttopm(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) ttopm(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged high cloud top pressure
            if(trim(fieldname)=='pres_avehct') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,ptoph,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  ptoph(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) ptoph(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged high cloud bottom pressure
            if(trim(fieldname)=='pres_avehcb') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,pboth,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  pboth(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) pboth(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged high cloud top temperature
            if(trim(fieldname)=='tmp_avehct') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,ttoph,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  ttoph(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) ttoph(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged boundary layer cloud cover
            if(trim(fieldname)=='tcdc_avebndcl') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,pblcfr,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  pblcfr(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  pblcfr(i,j) = spval
                  if (pblcfr(i,j) < spval) pblcfr(i,j) = pblcfr(i,j) * 0.01
                enddo
              enddo
            endif

            ! cloud work function
            if(trim(fieldname)=='cwork_aveclm') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,cldwork,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  cldwork(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  cldwork(i,j) = spval
                enddo
              enddo
            endif

            ! water runoff
            if(trim(fieldname)=='watr_acc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,runoff,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  runoff(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  runoff(i,j) = spval
                  if (sm(i,j) /= 0.0) runoff(i,j) = spval
                enddo
              enddo
            endif

            ! accumulated evaporation of intercepted water
            if(trim(fieldname)=='ecan_acc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,tecan,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  tecan(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) tecan(i,j) = spval
                  if (sm(i,j) /= 0.0) tecan(i,j) = spval
                enddo
              enddo
            endif

            ! accumulated plant transpiration
            if(trim(fieldname)=='etran_acc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,tetran,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  tetran(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) tetran(i,j) = spval
                  if (sm(i,j) /= 0.0) tetran(i,j) = spval
                enddo
              enddo
            endif

            ! accumulated soil surface evaporation
            if(trim(fieldname)=='edir_acc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,tedir,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  tedir(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) tedir(i,j) = spval
                  if (sm(i,j) /= 0.0) tedir(i,j) = spval
                enddo
              enddo
            endif

            ! total water storage in aquifer
            if(trim(fieldname)=='wa_acc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,twa,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  twa(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) twa(i,j) = spval
                  if (sm(i,j) /= 0.0) twa(i,j) = spval
                enddo
              enddo
            endif

            ! shelter max temperature
            if(modelname=='GFS') then
            if(trim(fieldname)=='tmax_max2m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,maxtshltr,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  maxtshltr(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  maxtshltr(i,j) = spval
                enddo
              enddo
            endif
            else
            if(trim(fieldname)=='t02max') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,maxtshltr,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  maxtshltr(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  maxtshltr(i,j) = spval
                enddo
              enddo
            endif
            endif

            ! shelter min temperature
            if(trim(fieldname)=='t02min' .or. trim(fieldname)=='tmin_min2m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,mintshltr,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  mintshltr(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  mintshltr(i,j) = spval
                enddo
              enddo
            endif

            ! shelter max rh
            if(trim(fieldname)=='rh02max') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,maxrhshltr,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  maxrhshltr(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  maxrhshltr(i,j) = spval
                enddo
              enddo
            endif

            ! shelter min rh
            if(trim(fieldname)=='rh02min') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,minrhshltr,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  minrhshltr(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  minrhshltr(i,j) = spval
                enddo
              enddo
            endif

            ! shelter max specific humidity
            if(trim(fieldname)=='spfhmax_max2m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,maxqshltr,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  maxqshltr(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) maxqshltr(i,j) = spval
                enddo
              enddo
            endif

            ! shelter min temperature
            if(trim(fieldname)=='spfhmin_min2m') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,minqshltr,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  minqshltr(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) minqshltr(i,j) = spval
                enddo
              enddo
            endif

            ! ice thickness
            if(trim(fieldname)=='icetk') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,dzice,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  dzice(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  dzice(i,j) = spval
                enddo
              enddo
            endif

            ! wilting point
            if(trim(fieldname)=='wilt') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend, spval,smcwlt,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  smcwlt(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  smcwlt(i,j) = spval
                  if (sm(i,j) /= 0.0) smcwlt(i,j) = spval
                enddo
              enddo
            endif

            ! sunshine duration
            if(trim(fieldname)=='sunsd_acc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,suntime,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  suntime(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  suntime(i,j) = spval
                enddo
              enddo
            endif

            ! field capacity
            if(trim(fieldname)=='fldcp') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,fieldcapa,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  fieldcapa(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  fieldcapa(i,j) = spval
                  if (sm(i,j) /= 0.0) fieldcapa(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged surface visible beam downward solar flux
            if(trim(fieldname)=='vbdsf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,avisbeamswin,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  avisbeamswin(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  avisbeamswin(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged surface visible diffuse downward solar flux
            if(trim(fieldname)=='vddsf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,avisdiffswin,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  avisdiffswin(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  avisdiffswin(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged surface near IR beam downward solar flux
            if(trim(fieldname)=='nbdsf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,airbeamswin,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  airbeamswin(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  airbeamswin(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged surface near IR diffuse downward solar flux
            if(trim(fieldname)=='nddsf_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,airdiffswin,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  airdiffswin(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  airdiffswin(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged surface clear sky outgoing LW
            if(trim(fieldname)=='csulf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,alwoutc,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  alwoutc(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  alwoutc(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged TOA clear sky outgoing LW
            if(trim(fieldname)=='csulftoa') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,alwtoac,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  alwtoac(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  alwtoac(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged surface clear sky outgoing SW
            if(trim(fieldname)=='csusf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,aswoutc,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  aswoutc(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  aswoutc(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged TOA clear sky outgoing SW
            if(trim(fieldname)=='csusftoa') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,aswtoac,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  aswtoac(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small)  aswtoac(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged surface clear sky incoming LW
            if(trim(fieldname)=='csdlf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,alwinc,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  alwinc(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) alwinc(i,j) = spval
                enddo
              enddo
            endif

            ! time averaged surface clear sky incoming SW
            if(trim(fieldname)=='csdsf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,aswinc,arrayr42d,fillValue,spval)
              do j=jsta,jend
                do i=ista, iend
                  aswinc(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) aswinc(i,j) = spval
                enddo
              enddo
            endif

            ! storm runoffs
            if(trim(fieldname)=='ssrun_acc') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,ssroff,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  ssroff(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) ssroff(i,j) = spval
                  if (sm(i,j) /= 0.0) ssroff(i,j) = spval
                enddo
              enddo
            endif

            !  direct soil evaporation
            if(trim(fieldname)=='evbs_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgedir,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  avgedir(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) avgedir(i,j) = spval
                  if (sm(i,j) /= 0.0) avgedir(i,j) = spval
                enddo
              enddo
            endif

            ! canopy water evap
            if(trim(fieldname)=='evcw_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgecan,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  avgecan(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) avgecan(i,j) = spval
                  if (sm(i,j) /= 0.0) avgecan(i,j) = spval
                enddo
              enddo
            endif

            ! AVERAGED PRECIP ADVECTED HEAT FLUX
            if(trim(fieldname)=='pah_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,paha,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  paha(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) paha(i,j) = spval
                  if (sm(i,j) /= 0.0) paha(i,j) = spval
                enddo
              enddo
            endif

            ! instantaneous PRECIP ADVECTED HEAT FLUX
            if(trim(fieldname)=='pahi') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,pahi,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  pahi(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) pahi(i,j) = spval
                  if (sm(i,j) /= 0.0) pahi(i,j) = spval
                enddo
              enddo
            endif

            ! plant transpiration
            if(trim(fieldname)=='trans_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgetrans,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  avgetrans(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) avgetrans(i,j) = spval
                  if (sm(i,j) /= 0.0) avgetrans(i,j) = spval
                enddo
              enddo
            endif

            ! snow sublimation
            if(trim(fieldname)=='sbsno_ave') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,avgesnow,arrayr42d,sm,sice,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  avgesnow(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) avgesnow(i,j) = spval
                  if (sm(i,j)==1.0 .and. sice(i,j)==0.) avgesnow(i,j)=spval
                enddo
              enddo
            endif

            ! total soil moisture
            if(trim(fieldname)=='soilm') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,smstot,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  smstot(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) smstot(i,j) = spval
                  if (sm(i,j) /= 0.0) smstot(i,j) = spval
                enddo
              enddo
            endif

            ! snow phase change heat flux
            if(trim(fieldname)=='snohf') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,snopcx,arrayr42d,sm,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  snopcx(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) snopcx(i,j) = spval
                  if (sm(i,j) /= 0.0) snopcx(i,j) = spval
                enddo
              enddo
            endif


            ! snow phase change heat flux
            if(trim(fieldname)=='pwat') then
              !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,pwat,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  pwat(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) pwat(i,j) = spval
                enddo
              enddo
            endif

            ! model level upvvelmax
            if(trim(fieldname)=='upvvelmax') then
              !$omp parallel do default(none) private(i,j,l) shared(jsta,jend,ista,iend,spval,w_up_max,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  w_up_max(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) w_up_max(i,j) = spval
                enddo
              enddo
            endif

            ! model level dnvvelmax
            if(trim(fieldname)=='dnvvelmax') then
              !$omp parallel do default(none) private(i,j,l) shared(jsta,jend,ista,iend,spval,w_dn_max,arrayr42d,fillValue)
              do j=jsta,jend
                do i=ista, iend
                  w_dn_max(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) w_dn_max(i,j) = spval
                enddo
              enddo
            endif

            ! model level uhmax25
            if(trim(fieldname)=='uhmax25') then
              !$omp parallel do default(none) private(i,j,l) shared(jsta,jend,ista,iend,spval,up_heli_max,arrayr42d,fillvalue)
              do j=jsta,jend
                do i=ista, iend
                  up_heli_max(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) up_heli_max(i,j) = spval
                enddo
              enddo
            endif

            ! model level uhmin25
            if(trim(fieldname)=='uhmin25') then
              !$omp parallel do default(none) private(i,j,l) shared(jsta,jend,ista,iend,spval,up_heli_min,arrayr42d,fillvalue)
              do j=jsta,jend
                do i=ista, iend
                  up_heli_min(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) up_heli_min(i,j) = spval
                enddo
              enddo
            endif

            ! model level uhmax03
            if(trim(fieldname)=='uhmax03') then
              !$omp parallel do default(none) private(i,j,l) shared(jsta,jend,ista,iend,spval,up_heli_max03,arrayr42d,fillvalue)
              do j=jsta,jend
                do i=ista, iend
                  up_heli_max03(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) up_heli_max03(i,j) = spval
                enddo
              enddo
            endif

            ! model level uhmin03
            if(trim(fieldname)=='uhmin03') then
              !$omp parallel do default(none) private(i,j,l) shared(jsta,jend,ista,iend,spval,up_heli_min03,arrayr42d,fillvalue)
              do j=jsta,jend
                do i=ista, iend
                  up_heli_min03(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) up_heli_min03(i,j) = spval
                enddo
              enddo
            endif

            ! model level maxvort01
            if(trim(fieldname)=='maxvort01') then
              !$omp parallel do default(none) private(i,j,l) shared(jsta,jend,ista,iend,spval,rel_vort_max01,arrayr42d,fillvalue)
              do j=jsta,jend
                do i=ista, iend
                  rel_vort_max01(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) rel_vort_max01(i,j) = spval
                enddo
              enddo
            endif

            ! model level maxvort02
            if(trim(fieldname)=='maxvort02') then
              !$omp parallel do default(none) private(i,j,l) shared(jsta,jend,ista,iend,spval,rel_vort_max,arrayr42d,fillvalue)
              do j=jsta,jend
                do i=ista, iend
                  rel_vort_max(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) rel_vort_max(i,j) = spval
                enddo
              enddo
            endif

            ! model level maxvorthy1
            if(trim(fieldname)=='maxvorthy1') then
              !$omp parallel do default(none) private(i,j,l) shared(jsta,jend,ista,iend,spval,rel_vort_maxhy1,arrayr42d,fillvalue)
              do j=jsta,jend
                do i=ista, iend
                  rel_vort_maxhy1(i,j) = arrayr42d(i,j)
                  if( abs(arrayr42d(i,j)-fillValue) < small) rel_vort_maxhy1(i,j) = spval
                enddo
              enddo
            endif


!          else if (fieldDimCount > gridDimCount) then
          else if (fieldDimCount ==3) then
!            print *,'in post_lam, get field value,n=',n,'fieldname=',trim(fieldname)
            if (typekind == ESMF_TYPEKIND_R4) then
              call ESMF_FieldGet(fcstField(n), localDe=0, farrayPtr=arrayr43d, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out

              call ESMF_AttributeGet(fcstField(n), convention="NetCDF", purpose="FV3", &
                 name="_FillValue", typekind=attTypeKind, isPresent=mvispresent, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
              if( mvispresent ) then
                if (attTypeKind==ESMF_TYPEKIND_R4) then
                 call ESMF_AttributeGet(fcstField(n), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=fillvalue, isPresent=mvispresent, rc=rc)
                else
                  call ESMF_AttributeGet(fcstField(n), convention="NetCDF", purpose="FV3", &
                        name="_FillValue", value=fillvalue8, isPresent=mvispresent, rc=rc)
                  fillvalue=fillvalue8
                endif
              endif
!              call ESMF_AttributeGet(fcstField(n), convention="NetCDF", purpose="FV3", &
!                   name='_FillValue', value=fillvalue, rc=rc)
!              print *,'in post_lam, get field value,fillvalue=',fillvalue
            else if (typekind == ESMF_TYPEKIND_R8) then
              call ESMF_FieldGet(fcstField(n), localDe=0, farrayPtr=arrayr83d, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
!              print *,'in post_lam, get field valuer8,n=',n,'fieldname=',trim(fieldname)
              call ESMF_AttributeGet(fcstField(n), convention="NetCDF", purpose="FV3", &
                   name='_FillValue', value=fillvalue8, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, file=__FILE__)) return  ! bail out
!              print *,'in post_lam, get field value,fillvalue8=',fillvalue8

              allocate(arrayr43d(ista:iend,jsta:jend,kstart:kend))
              arrayr43d = 0.
              fillvalue = fillvalue8
              do k=kstart,kend
             !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,k,arrayr43d,arrayr83d)
                do j=jsta,jend
                  do i=ista,iend
                    arrayr43d(i,j,k) = arrayr83d(i,j,k)
                  enddo
                enddo
              enddo
            endif

            ! biomass burning emissions
            if(trim(fieldname)=='ebu_smoke') then
              !$omp parallel do default(none) private(i,j,l) shared(jsta,jend,ista,iend,ebb,arrayr43d,fillValue,spval,lm)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    ebb(i,j,l,1)=arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillValue) < small) ebb(i,j,l,1)=spval
                  enddo
                enddo
              enddo
            endif

            ! model level T
            if(trim(fieldname)=='tmp') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,t,arrayr43d,fillvalue,spval)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    t(i,j,l) = arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue) < small) t(i,j,l) = spval
                  enddo
                enddo
              enddo
!              print *,'in post_lam,tmp 3d=',maxval(t(ista:iend,jsta:jend,1)),minval(t(ista:iend,jsta:jend,1)), &
!                'lm=',maxval(t(ista:iend,jsta:jend,lm)),minval(t(ista:iend,jsta:jend,lm)), &
!                t(ista,jsta,1),arrayr43d(ista,jsta,1),'fillvlaue=',fillvalue

              !! sig4
              !$omp parallel do default(none) private(i,j,tlmh) shared(lm,jsta,jend,ista,iend,t,sigt4,spval)
              do j=jsta,jend
                do i=ista, iend
                  if( t(i,j,lm) /= spval) then
                    tlmh = t(i,j,lm) * t(i,j,lm)
                    sigt4(i,j) = 5.67E-8 * tlmh * tlmh
                  else
                    sigt4(i,j)=spval
                  endif
                enddo
              enddo
            endif

            ! model level spfh
            if(trim(fieldname)=='spfh') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,q,arrayr43d,fillvalue,spval)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    q(i,j,l) = arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) q(i,j,l) = spval
                  enddo
                enddo
              enddo
            endif

            ! model level u wind
            if(trim(fieldname)=='ugrd') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,uh,arrayr43d,fillvalue,spval)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    uh(i,j,l) = arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue) < small) uh(i,j,l) = spval
                  enddo
                enddo
              enddo
!              print *,'in post_lam,uh 3d=',maxval(uh(ista:iend,jsta:jend,1)),minval(uh(ista:iend,jsta:jend,1)), &
!                'lm=',maxval(uh(ista:iend,jsta:jend,lm)),minval(uh(ista:iend,jsta:jend,lm)),'lm=',lm,'dim=',ista,iend,jsta,jend
            endif

            ! model level v wind
            if(trim(fieldname)=='vgrd') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,vh,arrayr43d,fillvalue,spval)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    vh(i,j,l) = arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) vh(i,j,l) = spval
                  enddo
                enddo
              enddo
            endif

            ! model level pressure thinkness
            if(trim(fieldname)=='dpres') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,dpres,arrayr43d,fillvalue,spval)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    dpres(i,j,l) = arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) dpres(i,j,l) = spval
                  enddo
                enddo
              enddo
            endif

            ! model level gh thinkness, model output negative delz
            if(trim(fieldname)=='delz') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,zint,arrayr43d,fillvalue,spval)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    zint(i,j,l) = spval
                    if(abs(arrayr43d(i,j,l)-fillvalue)>small) zint(i,j,l)=-1.*arrayr43d(i,j,l)
                  enddo
                enddo
              enddo
!              print *,'in post_lam,zint 3d=',maxval(zint(ista:iend,jsta:jend,1)),minval(zint(ista:iend,jsta:jend,1)), &
!                'lm=',maxval(zint(ista:iend,jsta:jend,lm)),minval(zint(ista:iend,jsta:jend,lm))
            endif

            ! model level w
            if(trim(fieldname)=='dzdt') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,wh,arrayr43d,fillvalue,spval)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    wh(i,j,l) = arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) wh(i,j,l) = spval
                  enddo
                enddo
              enddo
            endif

            ! model level omga
            if(trim(fieldname)=='omga') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,omga,arrayr43d,fillvalue,spval)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    omga(i,j,l) = arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) omga(i,j,l) = spval
                  enddo
                enddo
              enddo
            endif

            ! soilt
            if(trim(fieldname)=='soilt') then
              !$omp parallel do default(none) private(i,j,l) shared(nsoil,jsta,jend,ista,iend,stc,arrayr43d,sm,sice,fillvalue,spval)
              do l=1,nsoil
                do j=jsta,jend
                  do i=ista, iend
                    stc(i,j,l) = arrayr43d(i,j,l)
                    if( abs(arrayr43d(i,j,l)-fillValue) < small) stc(i,j,l) = spval
                    !mask open water areas, combine with sea ice tmp
                    if (sm(i,j) /= 0.0 .and. sice(i,j) ==0.) stc(i,j,l) = spval
                  enddo
                enddo
              enddo
            endif

            ! soilw
            if(trim(fieldname)=='soilw') then
              !$omp parallel do default(none) private(i,j,l) shared(nsoil,jsta,jend,ista,iend,smc,arrayr43d,sm,fillvalue,spval)
              do l=1,nsoil
                do j=jsta,jend
                  do i=ista, iend
                    smc(i,j,l) = arrayr43d(i,j,l)
                    if( abs(arrayr43d(i,j,l)-fillValue) < small) smc(i,j,l) = spval
                    if (sm(i,j) /= 0.0) smc(i,j,l) = spval
                  enddo
                enddo
              enddo
            endif

            ! soill
            if(trim(fieldname)=='soill') then
              !$omp parallel do default(none) private(i,j,l) shared(nsoil,jsta,jend,ista,iend,sh2o,arrayr43d,sm,fillvalue,spval)
              do l=1,nsoil
                do j=jsta,jend
                  do i=ista, iend
                    sh2o(i,j,l) = arrayr43d(i,j,l)
                    if( abs(arrayr43d(i,j,l)-fillValue) < small) sh2o(i,j,l) = spval
                    if (sm(i,j) /= 0.0) sh2o(i,j,l) = spval
                  enddo
                enddo
              enddo
            endif

            ! model level ozone mixing ratio
#ifdef MULTI_GASES
            if(trim(fieldname)=='spo3') then
#else
            if(trim(fieldname)=='o3mr') then
#endif
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,o3,arrayr43d,fillvalue,spval)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    o3(i,j,l) = arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) o3(i,j,l) = spval
                  enddo
                enddo
              enddo
            endif

! for GFDL MP or Thompson MP
            if (imp_physics == 11 .or. imp_physics == 8) then
              ! model level cloud water mixing ratio
              if(trim(fieldname)=='clwmr') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqw,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqw(i,j,l) = arrayr43d(i,j,l)
                      if(abs(arrayr43d(i,j,l)-fillvalue)<small) qqw(i,j,l)=spval
                    enddo
                  enddo
                enddo
!              print *,'in post_lam,clwmr 3d=',maxval(qqw(ista:iend,jsta:jend,1)),minval(qqw(ista:iend,jsta:jend,1)), &
!                'lm=',maxval(qqw(ista:iend,jsta:jend,lm)),minval(qqw(ista:iend,jsta:jend,lm)), &
!                qqw(ista,jsta,1),arrayr43d(ista,jsta,1),'fillvlaue=',fillvalue
              endif

              ! model level ice mixing ratio
              if(trim(fieldname)=='icmr') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqi,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqi(i,j,l) = arrayr43d(i,j,l)
                      if(abs(arrayr43d(i,j,l)-fillvalue)<small) qqi(i,j,l) = spval
                    enddo
                  enddo
                enddo
              endif

              ! model level rain water mixing ratio
              if(trim(fieldname)=='rwmr') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqr,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqr(i,j,l) = arrayr43d(i,j,l)
                      if(abs(arrayr43d(i,j,l)-fillvalue)<small) qqr(i,j,l) = spval
                    enddo
                  enddo
                enddo
              endif

              ! model level snow mixing ratio
              if(trim(fieldname)=='snmr') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqs,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqs(i,j,l) = arrayr43d(i,j,l)
                      if(abs(arrayr43d(i,j,l)-fillvalue)<small) qqs(i,j,l) = spval
                    enddo
                  enddo
                enddo
!              print *,'in post_lam,snmr 3d=',maxval(qqs(ista:iend,jsta:jend,1)),minval(qqs(ista:iend,jsta:jend,1)), &
!                'lm=',maxval(qqs(ista:iend,jsta:jend,lm)),minval(qqs(ista:iend,jsta:jend,lm)), &
!                qqs(ista,jsta,1),arrayr43d(ista,jsta,1),'fillvlaue=',fillvalue
              endif

              ! model level rain water mixing ratio
              if(trim(fieldname)=='grle') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqg,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqg(i,j,l) = arrayr43d(i,j,l)
                      if(abs(arrayr43d(i,j,l)-fillvalue)<small) qqg(i,j,l) = spval
                    enddo
                  enddo
                enddo
              endif

              ! model level hail mixing ratio
              if(trim(fieldname)=='hail') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqh,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqh(i,j,l) = arrayr43d(i,j,l)
                      if(abs(arrayr43d(i,j,l)-fillvalue)<small) qqh(i,j,l) = spval
                    enddo
                  enddo
                enddo
              endif

              if(imp_physics == 8) then
              ! model level rain water number
              if(trim(fieldname)=='rain_nc') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqnr,arrayr43d,spval,fillvalue)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqnr(i,j,l)=arrayr43d(i,j,l)
                      if(abs(arrayr43d(i,j,l)-fillvalue)<small) qqnr(i,j,l) = spval
                    enddo
                  enddo
                enddo
              endif

              ! model level cloud ice number
              if(trim(fieldname)=='nicp') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqni,arrayr43d,spval,fillvalue)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqni(i,j,l)=arrayr43d(i,j,l)
                      if(abs(arrayr43d(i,j,l)-fillvalue)<small) qqni(i,j,l) = spval
                    enddo
                  enddo
                enddo
              endif

              ! model level cloud water number
              if(trim(fieldname)=='water_nc') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqnw,arrayr43d,spval,fillvalue)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqnw(i,j,l)=arrayr43d(i,j,l)
                      if(abs(arrayr43d(i,j,l)-fillvalue)<small) qqnw(i,j,l) = spval
                    enddo
                  enddo
                enddo
              endif

              ! model level rain number
              if(trim(fieldname)=='nwfa') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqnwfa,arrayr43d,spval,fillvalue)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqnwfa(i,j,l)=arrayr43d(i,j,l)
                      if(abs(arrayr43d(i,j,l)-fillvalue)<small) qqnwfa(i,j,l) = spval
                    enddo
                  enddo
                enddo
              endif

              ! model level rain number
              if(trim(fieldname)=='nifa') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,qqnifa,arrayr43d,spval,fillvalue)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      qqnifa(i,j,l)=arrayr43d(i,j,l)
                      if(abs(arrayr43d(i,j,l)-fillvalue)<small) qqnifa(i,j,l) = spval
                    enddo
                  enddo
                enddo
              endif

              endif !if(imp_physics == 8) then

            endif !if(imp_physics == 11 .or. imp_physics == 8) then

            ! model level ref3d
            if(trim(fieldname)=='ref3D') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,ref_10cm,arrayr43d,fillvalue,spval)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    ref_10cm(i,j,l) = arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) ref_10cm(i,j,l) = spval
                  enddo
                enddo
              enddo
            endif
!              if(mype==0) print *,'in gfs_post, get ref_10cm=',maxval(ref_10cm), minval(ref_10cm)
            if(trim(fieldname)=='refl_10cm') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,ref_10cm,arrayr43d,fillvalue,spval)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    ref_10cm(i,j,l) = arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) ref_10cm(i,j,l) = spval
                  enddo
                enddo
              enddo
!              if(mype==0) print *,'in gfs_post, get ref_10cm=',maxval(ref_10cm), minval(ref_10cm),'ibdl=',ibdl
            endif

            ! model level tke
            if(trim(fieldname)=='qke') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,q2,arrayr43d, fillvalue,spval)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    q2(i,j,l)=arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) then
                      q2(i,j,l) = spval
                    else
                      q2(i,j,l) = q2(i,j,l)/2.0
                    endif
                  enddo
                enddo
!              print *,'in gfs_post, get tke=',maxval(q2(:,:,l)), minval(q2(:,:,l)),'l=',l
              enddo
            endif

            ! model level cloud fraction
            if(imp_physics == 11) then !GFDL MP
              if(trim(fieldname)=='cld_amt') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,cfr,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      cfr(i,j,l) = arrayr43d(i,j,l)
                      if(abs(arrayr43d(i,j,l)-fillvalue)<small) cfr(i,j,l) = spval
                    enddo
                  enddo
                enddo
              endif
            else !Other MP
              ! read cldfra_bl first
              if(trim(fieldname)=='cldfra_bl') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,cfr,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      cfr(i,j,l) = arrayr43d(i,j,l)
                      if(abs(arrayr43d(i,j,l)-fillvalue)<small) cfr(i,j,l) = spval
                    enddo
                  enddo
                enddo
              endif
              if(trim(fieldname)=='cldfra') then
                !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,cfr,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      cfr(i,j,l) = arrayr43d(i,j,l)
                      if(abs(arrayr43d(i,j,l)-fillvalue)<small) cfr(i,j,l) = spval
                    enddo
                  enddo
                enddo
              endif
            endif

            if(modelname=='FV3R') then

            ! model level smoke
            if(trim(fieldname)=='smoke') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,smoke,arrayr43d,spval,fillvalue)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    smoke(i,j,l,1)=arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) smoke(i,j,l,1) = spval
                  enddo
                enddo
              enddo
            endif

            ! model level dust
            if(trim(fieldname)=='dust') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,fv3dust,arrayr43d,spval,fillvalue)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    fv3dust(i,j,l,1)=arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) fv3dust(i,j,l,1) = spval
                  enddo
                enddo
              enddo
            endif

            ! model level ext550 extinction
            if(trim(fieldname)=='ext550') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,ext550,arrayr43d,spval,fillvalue)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    ext550(i,j,l)=arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) ext550(i,j,l) = spval
                  enddo
                enddo
              enddo
            endif

            ! model level coarse dust
            if(trim(fieldname)=='coarsepm') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,coarsepm,arrayr43d,spval,fillvalue)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    coarsepm(i,j,l,1)=arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) coarsepm(i,j,l,1) = spval
                  enddo
                enddo
              enddo
            endif

            endif !end FV3R

            ! Thompson scheme cloud ice effective radius
            if(trim(fieldname)=='cieffr') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,effri,arrayr43d,spval,fillvalue)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    effri(i,j,l)=arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) effri(i,j,l) = spval
                  enddo
                enddo
              enddo
            endif

            ! Thompson scheme cloud water effective radius
            if(trim(fieldname)=='cleffr') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,effrl,arrayr43d,spval,fillvalue)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    effrl(i,j,l)=arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) effrl(i,j,l) = spval
                  enddo
                enddo
              enddo
            endif

            ! Thompson scheme cloud snow effective radius
            if(trim(fieldname)=='cseffr') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,effrs,arrayr43d,spval,fillvalue)
              do l=1,lm
                do j=jsta,jend
                  do i=ista, iend
                    effrs(i,j,l)=arrayr43d(i,j,l)
                    if(abs(arrayr43d(i,j,l)-fillvalue)<small) effrs(i,j,l) = spval
                  enddo
                enddo
              enddo
            endif

            ! read chemical fields
            if(gocart_on .or. gccpp_on .or. nasa_on) then

              if(trim(fieldname)=='dust1') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,dust,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      dust(i,j,l,1) = max(arrayr43d(i,j,l),0.0)
                      if(abs(arrayr43d(i,j,l)-fillvalue) < small) dust(i,j,l,1) = spval
                    enddo
                  enddo
                enddo
              endif

              if(trim(fieldname)=='dust2') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,dust,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      dust(i,j,l,2) = max(arrayr43d(i,j,l),0.0)
                      if(abs(arrayr43d(i,j,l)-fillvalue) < small) dust(i,j,l,2) = spval
                    enddo
                  enddo
                enddo
              endif

              if(trim(fieldname)=='dust3') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,dust,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      dust(i,j,l,3) = max(arrayr43d(i,j,l), 0.0)
                      if(abs(arrayr43d(i,j,l)-fillvalue) < small) dust(i,j,l,3) = spval
                    enddo
                  enddo
                enddo
              endif

              if(trim(fieldname)=='dust4') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,dust,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      dust(i,j,l,4) = max(arrayr43d(i,j,l), 0.0)
                      if(abs(arrayr43d(i,j,l)-fillvalue) < small) dust(i,j,l,4) = spval
                    enddo
                  enddo
                enddo
              endif

              if(trim(fieldname)=='dust5') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,dust,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      dust(i,j,l,5) = max(arrayr43d(i,j,l), 0.0)
                      if(abs(arrayr43d(i,j,l)-fillvalue) < small) dust(i,j,l,5) = spval
                    enddo
                  enddo
                enddo
              endif

              if(trim(fieldname)=='seas1') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,salt,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      salt(i,j,l,1) = max(arrayr43d(i,j,l),0.0)
                      if(abs(arrayr43d(i,j,l)-fillvalue) < small) salt(i,j,l,1) = spval
                    enddo
                  enddo
                enddo
              endif

              if(trim(fieldname)=='seas2') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,salt,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      salt(i,j,l,2) = max(arrayr43d(i,j,l),0.0)
                      if(abs(arrayr43d(i,j,l)-fillvalue) < small) salt(i,j,l,2) = spval
                    enddo
                  enddo
                enddo
              endif

              if(trim(fieldname)=='seas3') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,salt,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      salt(i,j,l,3) = max(arrayr43d(i,j,l),0.0)
                      if(abs(arrayr43d(i,j,l)-fillvalue) < small) salt(i,j,l,3) = spval
                    enddo
                  enddo
                enddo
              endif

              if(trim(fieldname)=='seas4') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,salt,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      salt(i,j,l,4) = max(arrayr43d(i,j,l),0.0)
                      if(abs(arrayr43d(i,j,l)-fillvalue) < small) salt(i,j,l,4) = spval
                    enddo
                  enddo
                enddo
              endif

              if(trim(fieldname)=='seas5') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,salt,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      salt(i,j,l,5) = max(arrayr43d(i,j,l),0.0)
                      if(abs(arrayr43d(i,j,l)-fillvalue) < small) salt(i,j,l,5) = spval
                    enddo
                  enddo
                enddo
              endif

              if(trim(fieldname)=='bc1') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,soot,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      soot(i,j,l,1) = max(arrayr43d(i,j,l),0.0)
                      if(abs(arrayr43d(i,j,l)-fillvalue) < small) soot(i,j,l,1) = spval
                    enddo
                  enddo
                enddo
              endif

              if(trim(fieldname)=='bc2') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,soot,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      soot(i,j,l,2) = max(arrayr43d(i,j,l),0.0)
                      if(abs(arrayr43d(i,j,l)-fillvalue) < small) soot(i,j,l,2) = spval
                    enddo
                  enddo
                enddo
              endif

              if(trim(fieldname)=='oc1') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,waso,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      waso(i,j,l,1) = max(arrayr43d(i,j,l),0.0)
                      if(abs(arrayr43d(i,j,l)-fillvalue) < small) waso(i,j,l,1) = spval
                    enddo
                  enddo
                enddo
              endif

              if(trim(fieldname)=='oc2') then
              !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,waso,arrayr43d,fillvalue,spval)
                do l=1,lm
                  do j=jsta,jend
                    do i=ista, iend
                      waso(i,j,l,2) = max(arrayr43d(i,j,l),0.0)
                      if(abs(arrayr43d(i,j,l)-fillvalue) < small) waso(i,j,l,2) = spval
                    enddo
                  enddo
                enddo
              endif

              if (gocart_on .or. gccpp_on) then
                if(trim(fieldname)=='sulf') then
                 !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,suso,arrayr43d,fillvalue,spval)
                  do l=1,lm
                    do j=jsta,jend
                      do i=ista, iend
                        suso(i,j,l,1) = max(arrayr43d(i,j,l),0.0)
                        if(abs(arrayr43d(i,j,l)-fillvalue) < small) suso(i,j,l,1) = spval
                      enddo
                    enddo
                  enddo
                endif

                if(trim(fieldname)=='pp25') then
                 !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,pp25,arrayr43d,fillvalue,spval)
                  do l=1,lm
                    do j=jsta,jend
                      do i=ista, iend
                        pp25(i,j,l,1) = max(arrayr43d(i,j,l),0.0)
                        if(abs(arrayr43d(i,j,l)-fillvalue) < small) pp25(i,j,l,1) = spval
                      enddo
                    enddo
                  enddo
                endif

                if(trim(fieldname)=='pp10') then
                 !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,pp10,arrayr43d,fillvalue,spval)
                  do l=1,lm
                    do j=jsta,jend
                      do i=ista, iend
                        pp10(i,j,l,1) = max(arrayr43d(i,j,l),0.0)
                        if(abs(arrayr43d(i,j,l)-fillvalue) < small) pp10(i,j,l,1) = spval
                      enddo
                    enddo
                  enddo
                endif

              else if (nasa_on) then
                if(trim(fieldname)=='so4') then
                 !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,suso,arrayr43d,fillvalue,spval)
                  do l=1,lm
                    do j=jsta,jend
                      do i=ista, iend
                        suso(i,j,l,1) = max(arrayr43d(i,j,l),0.0)
                        if(abs(arrayr43d(i,j,l)-fillvalue) < small) suso(i,j,l,1) = spval
                      enddo
                    enddo
                  enddo
                endif

                if(trim(fieldname)=='no3an1') then
                 !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,no3,arrayr43d,fillvalue,spval)
                  do l=1,lm
                    do j=jsta,jend
                      do i=ista, iend
                        no3(i,j,l,1) = max(arrayr43d(i,j,l),0.0)
                        if(abs(arrayr43d(i,j,l)-fillvalue) < small) no3(i,j,l,1) = spval
                      enddo
                    enddo
                  enddo
                endif

                if(trim(fieldname)=='no3an2') then
                 !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,no3,arrayr43d,fillvalue,spval)
                  do l=1,lm
                    do j=jsta,jend
                      do i=ista, iend
                        no3(i,j,l,2) = max(arrayr43d(i,j,l),0.0)
                        if(abs(arrayr43d(i,j,l)-fillvalue) < small) no3(i,j,l,2) = spval
                      enddo
                    enddo
                  enddo
                endif

                if(trim(fieldname)=='no3an3') then
                 !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,no3,arrayr43d,fillvalue,spval)
                  do l=1,lm
                    do j=jsta,jend
                      do i=ista, iend
                        no3(i,j,l,3) = max(arrayr43d(i,j,l),0.0)
                        if(abs(arrayr43d(i,j,l)-fillvalue) < small) no3(i,j,l,3) = spval
                      enddo
                    enddo
                  enddo
                endif

                if(trim(fieldname)=='nh4a') then
                 !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,nh4,arrayr43d,fillvalue,spval)
                  do l=1,lm
                    do j=jsta,jend
                      do i=ista, iend
                        nh4(i,j,l,1) = max(arrayr43d(i,j,l),0.0)
                        if(abs(arrayr43d(i,j,l)-fillvalue) < small) nh4(i,j,l,1) = spval
                      enddo
                    enddo
                  enddo
                endif

                if(trim(fieldname)=='pm25') then
                 !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,pp25,arrayr43d,fillvalue,spval)
                  do l=1,lm
                    do j=jsta,jend
                      do i=ista, iend
                        pp25(i,j,l,1) = max(arrayr43d(i,j,l),0.0)
                        if(abs(arrayr43d(i,j,l)-fillvalue) < small) pp25(i,j,l,1) = spval
                      enddo
                    enddo
                  enddo
                endif

                if(trim(fieldname)=='pm10') then
                 !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,pp10,arrayr43d,fillvalue,spval)
                  do l=1,lm
                    do j=jsta,jend
                      do i=ista, iend
                        pp10(i,j,l,1) = max(arrayr43d(i,j,l),0.0)
                        if(abs(arrayr43d(i,j,l)-fillvalue) < small) pp10(i,j,l,1) = spval
                      enddo
                    enddo
                  enddo
                endif

              endif !nasa_on

            endif !end gocart_on, gccpp_on, nasa_on


!3d fields
          endif

! end loop ncount_field
        enddo

        deallocate(fcstField)

! end file_loop_all
      enddo file_loop_all

! recompute full layer of zint
!$omp parallel do default(none) private(i,j) shared(jsta,jend,lp1,spval,zint,fis,ista,iend)
      do j=jsta,jend
        do i=ista,iend
          if (fis(i,j) /= spval) then
            zint(i,j,lp1) = fis(i,j)
            fis(i,j)      = fis(i,j) * grav
          else
            zint(i,j,lp1) = spval
            fis(i,j)      = spval
          endif
        enddo
      enddo


! compute pint from top down
!$omp parallel do default(none) private(i,j) shared(jsta,jend,ak5,pint,ista,iend)
      do j=jsta,jend
        do i=ista,iend
          pint(i,j,1) = ak5(1)
        end do
      end do

      do l=2,lp1
!$omp parallel do default(none) private(i,j) shared(l,jsta,jend,pint,dpres,spval,ista,iend)
        do j=jsta,jend
          do i=ista,iend
            if(dpres(i,j,l-1) /= spval) then
              pint(i,j,l) = pint(i,j,l-1) + dpres(i,j,l-1)
            else
              pint(i,j,l) = spval
            endif
          enddo
        enddo
      end do

!compute pmid from averaged two layer pint
      do l=lm,1,-1
!$omp parallel do default(none) private(i,j) shared(l,jsta,jend,pmid,pint,spval,ista,iend)
        do j=jsta,jend
          do i=ista,iend
            if(pint(i,j,l+1) /= spval) then
              pmid(i,j,l) = 0.5*(pint(i,j,l)+pint(i,j,l+1))
            else
              pmid(i,j,l) = spval
            endif
          enddo
        enddo
      enddo

!$omp parallel do default(none) private(i,j) shared(jsta,jend,spval,pt,pd,pint,ista,iend)
      do j=jsta,jend
        do i=ista,iend
          pd(i,j)     = spval
          pint(i,j,1) = pt
        end do
      end do
!      print *,'in setvar, pt=',pt,'ak5(lp1)=', ak5(lp1),'ak5(1)=',ak5(1)

! compute alpint
      do l=lp1,1,-1
!$omp parallel do default(none) private(i,j) shared(l,jsta,jend,alpint,pint,spval,ista,iend)
        do j=jsta,jend
          do i=ista,iend
            if(pint(i,j,l) /= spval) then
              alpint(i,j,l) = log(pint(i,j,l))
            else
              alpint(i,j,l) = spval
            endif
          end do
        end do
      end do

      do l=lm,1,-1
!$omp parallel do default(none) private(i,j) shared(l,jsta,jend,omga,wh,dpres,zint,alpint,q,t,spval,ista,iend)
        do j=jsta,jend
          do i=ista,iend
            if(wh(i,j,l) /= spval) then
              if(omga(i,j,l) == spval .and. dpres(i,j,l) /= spval .and. zint(i,j,l) /=spval)  &
                  omga(i,j,l) = (-1.) * wh(i,j,l) * dpres(i,j,l)/zint(i,j,l)
              if(zint(i,j,l+1) /=spval .and. zint(i,j,l) /=spval) &
                  zint(i,j,l) = zint(i,j,l) + zint(i,j,l+1)
            else
              if(zint(i,j,l+1) /=spval .and. t(i,j,l) /= spval .and.  alpint(i,j,l+1) /= spval  &
                             .and. alpint(i,j,l) /=spval .and. q(i,j,l) /= spval) then
                 zint(i,j,l) = zint(i,j,l+1)+(rgas/grav)*t(i,j,l)*(1.+fv*q(i,j,l))*(alpint(i,j,l+1)-alpint(i,j,l))
               else 
                 zint(i,j,l) = spval
              endif
            endif
          enddo
        enddo
      enddo
!      print *,'in post_lam,omga 3d=',maxval(omga(ista:iend,jsta:jend,1)),minval(omga(ista:iend,jsta:jend,1)), &
!           'lm=',maxval(omga(ista:iend,jsta:jend,lm)),minval(omga(ista:iend,jsta:jend,lm))

! compute zmid
      do l=lm,1,-1
!$omp parallel do default(none) private(i,j) shared(l,jsta,jend,zmid,zint,pmid,alpint,spval,ista,iend)
        do j=jsta,jend
          do i=ista,iend
            if( zint(i,j,l+1)/=spval .and. zint(i,j,l)/=spval .and. pmid(i,j,l) /= spval) then
              zmid(i,j,l)=zint(i,j,l+1)+(zint(i,j,l)-zint(i,j,l+1))* &
                    (log(pmid(i,j,l))-alpint(i,j,l+1))/ &
                    (alpint(i,j,l)-alpint(i,j,l+1))
            else
              zmid(i,j,l) = spval
            endif
          end do
        end do
      end do


! surface potential T, and potential T at roughness length
!$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,lp1,sm,ths,sst,thz0,sice,pint)
      do j=jsta,jend
        do i=ista, iend
          !assign sst
          if (sm(i,j) /= 0.0 .and. ths(i,j) /= spval) then
            if (sice(i,j) >= 0.15) then
              sst(i,j) = 271.4
            else
             sst(i,j) = ths(i,j)
            endif
          else
             sst(i,j) = spval
          endif
          if (ths(i,j) /= spval) then
            ths(i,j)  = ths(i,j)* (p1000/pint(i,j,lp1))**capa
            thz0(i,j) = ths(i,j)
          else
            thz0(i,j) = spval
          endif
        enddo
      enddo
!      print *,'in post_gfs,ths=',maxval(ths(1:im,jsta:jend)), &
!          minval(ths(1:im,jsta:jend))

! compute cwm and max qrain in the column to be used later in precip type computation
        do j=jsta,jend
          do i=ista,iend
            qrmax(i,j)=0.
          enddo
        enddo
        do l=1,lm
!$omp parallel do default(none) private(i,j) shared(l,jsta,jend,ista,iend,cwm,qrmax,qqg,qqs,qqr,qqi,qqw,qqh,spval)
          do j=jsta,jend
            do i=ista,iend
              if( qqr(i,j,l) /= spval) then
                cwm(i,j,l) = qqg(i,j,l)+qqs(i,j,l)+qqr(i,j,l)+qqi(i,j,l)+qqw(i,j,l)
                qrmax(i,j)=max(qrmax(i,j),qqr(i,j,l))
                if(qqh(i,j,l) /= spval) then
                  cwm(i,j,l) = cwm(i,j,l)+qqh(i,j,l)
                endif
              else
                cwm(i,j,l) = spval
              endif
            enddo
          enddo
        enddo

! estimate 2m pres and convert t2m to theta
!$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,lm,pshltr,pint,tshltr,spval)
      do j=jsta,jend
        do i=ista, iend
          if( tshltr(i,j) /= spval) then
            pshltr(I,j) = pint(i,j,lm+1)*EXP(-0.068283/tshltr(i,j))
            tshltr(i,j) = tshltr(i,j)*(p1000/pshltr(I,J))**CAPA
          else
            pshltr(I,J) = spval
            tshltr(i,j) = spval
          endif
        enddo
      enddo
!      print *,'in post_gfs,tshltr=',maxval(tshltr(1:im,jsta:jend)), &
!          minval(tshltr(1:im,jsta:jend))

!htop
      do j=jsta,jend
        do i=ista,iend
          htop(i,j) = spval
          if(ptop(i,j) < spval)then
            do l=1,lm
              if(ptop(i,j) <= pmid(i,j,l))then
                htop(i,j)=l
                exit
              end if
            end do
          end if
        end do
      end do

! hbot
      do j=jsta,jend
        do i=ista,iend
          hbot(i,j) = spval
          if(pbot(i,j) < spval)then
            do l=lm,1,-1
              if(pbot(i,j) >= pmid(i,j,l)) then
                hbot(i,j) = l
                exit
              end if
            end do
          end if
        end do
      end do

      if(modelname=='FV3R') then
        ! smoke and dust extinction
        !$omp parallel do default(none) private(i,j,l) shared(lm,jsta,jend,ista,iend,zint,taod5503d,aextc55,ext550,spval)
        do l=1,lm
          do j=jsta,jend
            do i=ista, iend
              if(ext550(i,j,l)<spval) then
                taod5503d(i,j,l)=ext550(i,j,l)
                aextc55(i,j,l)=taod5503d(i,j,l)/(zint(i,j,l)-zint(i,j,l+1))
              else
                taod5503d(i,j,l)=spval
                aextc55(i,j,l)=0.
              endif
            enddo
          enddo
        enddo
        deallocate(ext550)
      endif !end FV3R

        !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,snacc_ice,snacc_land,sndepac)
        do j=jsta,jend
          do i=ista, iend
            if(snacc_land(i,j)<spval) then
              sndepac(i,j) = snacc_land(i,j)
            elseif(snacc_ice(i,j)<spval) then
              sndepac(i,j) = snacc_ice(i,j)
            else
              sndepac(i,j) = spval
            endif
          enddo
        enddo

        !$omp parallel do default(none) private(i,j) shared(jsta,jend,ista,iend,spval,acsnom_ice,acsnom_land,acsnom)
        do j=jsta,jend
          do i=ista, iend
            if(acsnom_land(i,j)<spval) then
              acsnom(i,j) = acsnom_land(i,j)
            elseif(acsnom_ice(i,j)<spval) then
              acsnom(i,j) = acsnom_ice(i,j)
            else
              acsnom(i,j) = spval
            endif
          enddo
        enddo

        deallocate(snacc_ice)
        deallocate(snacc_land)
        deallocate(acsnom_ice)
        deallocate(acsnom_land)


      ! chmical field computation
      if(gocart_on .or. gccpp_on .or. nasa_on) then
        dustcb=0.0
        dustallcb=0.0
        do l=1,lm
          do j=jsta,jend
            do i=ista, iend
              if(dust(i,j,l,1)<spval.and.dust(i,j,l,2)<spval.and. &
                 dust(i,j,l,3)<spval.and.dust(i,j,l,4)<spval.and. &
                 dust(i,j,l,5)<spval) then
                dustcb(i,j)=dustcb(i,j)+ &
                  (dust(i,j,l,1)+0.38*dust(i,j,l,2))* &
                  dpres(i,j,l)/grav
                dustallcb(i,j)=dustallcb(i,j)+ &
                  (dust(i,j,l,1)+dust(i,j,l,2)+ &
                  dust(i,j,l,3)+0.74*dust(i,j,l,4))* &
                  dpres(i,j,l)/grav
              endif
            enddo
          enddo
        enddo

        sscb=0.0
        ssallcb=0.0
        do l=1,lm
          do j=jsta,jend
            do i=ista,iend
              if(salt(i,j,l,1)<spval.and.salt(i,j,l,2)<spval.and. &
                 salt(i,j,l,3)<spval.and.salt(i,j,l,4)<spval.and. &
                 salt(i,j,l,5)<spval) then
            sscb(i,j)=sscb(i,j)+ &
         (salt(i,j,l,1)+salt(i,j,l,2)+0.83*salt(i,j,l,3))*  &
           dpres(i,j,l)/grav


          ssallcb(i,j)=ssallcb(i,j)+ &
         (salt(i,j,l,1)+salt(i,j,l,2)+salt(i,j,l,3)+salt(i,j,l,4))* &
           dpres(i,j,l)/grav
              endif
            enddo
          enddo
        end do

        bccb=0.0
        occb=0.0
        do l=1,lm
          do j=jsta,jend
            do i=ista,iend
              if(soot(i,j,l,1)<spval.and.soot(i,j,l,2)<spval)then
               bccb(i,j)=bccb(i,j)+(soot(i,j,l,1)+soot(i,j,l,2))* &
               dpres(i,j,l)/grav
              endif

              if(waso(i,j,l,1)<spval.and.waso(i,j,l,2)<spval)then
               occb(i,j)=occb(i,j)+ (waso(i,j,l,1)+waso(i,j,l,2))* &
               dpres(i,j,l)/grav
              endif
            enddo
          enddo
        end do

        if(nasa_on) then
          no3cb=0.0
          nh4cb=0.0
          do l=1,lm
            do j=jsta,jend
              do i=ista,iend
                if(no3(i,j,l,1)<spval .and. no3(i,j,l,2)<spval .and. &
                   no3(i,j,l,3)<spval) then
                   no3cb(i,j)=no3cb(i,j)+ (no3(i,j,l,1)+no3(i,j,l,2)+ &
                   no3(i,j,l,3) ) * dpres(i,j,l)/grav
                else
                   no3(i,j,l,1)=0.
                   no3(i,j,l,2)=0.
                   no3(i,j,l,3)=0.
                endif
                if(nh4(i,j,l,1)<spval)then
                   nh4cb(i,j)=nh4cb(i,j)+ nh4(i,j,l,1)* &
                   dpres(i,j,l)/grav
                endif
              enddo
            enddo
          end do
        endif !end nasa_on

        sulfcb=0.0
        pp25cb=0.0
        pp10cb=0.0
        do l=1,lm
          do j=jsta,jend
            do i=ista,iend
              if(suso(i,j,l,1)<spval)then
                sulfcb(i,j)=sulfcb(i,j)+ suso(i,j,l,1)* &
                  dpres(i,j,l)/grav
              endif
              if(pp25(i,j,l,1)<spval)then
                pp25cb(i,j)=pp25cb(i,j)+ pp25(i,j,l,1)* &
                  dpres(i,j,l)/grav
              endif
              if(pp10(i,j,l,1)<spval)then
                pp10cb(i,j)=pp10cb(i,j)+ pp10(i,j,l,1)* &
                  dpres(i,j,l)/grav
              endif
            enddo
          enddo
        enddo ! do loop for l

        do l=1,lm
         do j=jsta,jend
          do i=ista,iend
            tv = t(i,j,l) * (h1+d608*MAX(q(I,J,L),qmin))
            rhomid(i,j,l) = pmid(i,j,l) / (rd*tv)
          enddo
         enddo
        enddo

        l=lm
        do j=jsta,jend
          do i=ista,iend

            dustcb(i,j) = MAX(dustcb(i,j), 0.0)
            dustallcb(i,j) = MAX(dustallcb(i,j), 0.0)
            sscb(i,j) = MAX(sscb(i,j), 0.0)
            ssallcb(i,j) = MAX(ssallcb(i,j), 0.0)
            bccb(i,j) = MAX(bccb(i,j), 0.0)
            occb(i,j) = MAX(occb(i,j), 0.0)
            sulfcb(i,j) = MAX(sulfcb(i,j), 0.0)
            if(nasa_on) then
              no3cb(i,j) = MAX(no3cb(i,j), 0.0)
              nh4cb(i,j) = MAX(nh4cb(i,j), 0.0)
            endif
            pp25cb(i,j) = MAX(pp25cb(i,j), 0.0)
            pp10cb(i,j) = MAX(pp10cb(i,j), 0.0)

           ! Surface PM25 dust and seasalt
           dustpm(i,j)=(dust(i,j,l,1)+0.38*dust(i,j,l,2))*rhomid(i,j,l) !ug/m3
           dustpm10(i,j)=(dust(i,j,l,1)+dust(i,j,l,2)+dust(i,j,l,3)+ &
             0.74*dust(i,j,l,4))*rhomid(i,j,l) !ug/m3
           sspm(i,j)=(salt(i,j,l,1)+salt(i,j,l,2)+ &
             0.83*salt(i,j,l,3))*rhomid(i,j,l)  !ug/m3

           if(gocart_on .or. gccpp_on) then

             !Surface PM10 concentration
             dusmass(i,j)=(dust(i,j,l,1)+dust(i,j,l,2)+dust(i,j,l,3)+ &
               0.74*dust(i,j,l,4)+salt(i,j,l,1)+salt(i,j,l,2)+salt(i,j,l,3)+ &
               salt(i,j,l,4) + soot(i,j,l,1)+soot(i,j,l,2)+waso(i,j,l,1)+ &
               waso(i,j,l,2) +suso(i,j,l,1)+pp25(i,j,l,1)+pp10(i,j,l,1)) &
               *rhomid(i,j,l)  !ug/m3
             !Surface PM25 concentration
             dusmass25(i,j)=(dust(i,j,l,1)+0.38*dust(i,j,l,2)+ &
               salt(i,j,l,1)+salt(i,j,l,2)+0.83*salt(i,j,l,3) + &
               soot(i,j,l,1)+soot(i,j,l,2)+waso(i,j,l,1)+ &
               waso(i,j,l,2) +suso(i,j,l,1)+pp25(i,j,l,1))*rhomid(i,j,l)  !ug/m3

             !PM10 column
             ducmass(i,j)=dustallcb(i,j)+ssallcb(i,j)+bccb(i,j)+ &
               occb(i,j)+sulfcb(i,j)+pp25cb(i,j)+pp10cb(i,j)
             !PM25 column
             ducmass25(i,j)=dustcb(i,j)+sscb(i,j)+bccb(i,j)+occb(i,j) &
               +sulfcb(i,j)+pp25cb(i,j)

           elseif(nasa_on) then
             !Surface PM10 concentration
             dusmass(i,j)=pp10(i,j,l,1)*rhomid(i,j,l)  !ug/m3
             !Surface PM25 concentration
             dusmass25(i,j)=pp25(i,j,l,1)*rhomid(i,j,l)  !ug/m3

             !PM10 column
             ducmass(i,j)=pp10cb(i,j)
             !PM25 column
             ducmass25(i,j)=pp25cb(i,j)

           endif !nasa_on
         end do
       end do

      endif

! generate look up table for lifted parcel calculations
      thl    = 210.
      plq    = 70000.
      pt_tbl = 10000.          ! this is for 100 hPa added by Moorthi

      call table(ptbl,ttbl,pt_tbl,                                     &
                 rdq,rdth,rdp,rdthe,pl,thl,qs0,sqs,sthe,the0)

      call tableq(ttblq,rdpq,rdtheq,plq,thl,stheq,the0q)

      if(mype == 0)then
        write(6,*)'  SPL (POSTED PRESSURE LEVELS) BELOW: '
        write(6,51) (SPL(L),L=1,LSM)
   50   format(14(F4.1,1X))
   51   format(8(F8.1,1X))
      endif
!
!$omp parallel do default(none) private(l) shared(lsm,alsl,spl)
      do l = 1,lsm
         alsl(l) = log(spl(l))
      end do
!
!      print *,'in gfs_post, end ref_10cm=',maxval(ref_10cm), minval(ref_10cm)
!!! above is fv3 change
!
!more fields need to be computed
!

! write lat/lon of the four corner point for rotated lat-lon grid
      if(mype == 0 .and. maptype == 207)then
        write(flatlon,1001)ifhr
        open(112,file=trim(flatlon),form='formatted',status='unknown')
        write(112,1002)latstart/1000,lonstart/1000,&
        latse/1000,lonse/1000,latnw/1000,lonnw/1000, &
        latlast/1000,lonlast/1000
 1001   format('latlons_corners.txt.f',I3.3)
 1002   format(4(I6,I7,X))
        close(112)
      endif
    end subroutine set_postvars_fv3
!
!-----------------------------------------------------------------------
!
end module post_fv3
