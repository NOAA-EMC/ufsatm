!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
    subroutine post_alctvars(imi,jmi,lmi,mype,nwtlpes,lead_write, mpicomp,  &
                             jts,jte,jtsgrp,jtegrp,its,ite,itsgrp,itegrp)
!
!
!   revision history:
!    Jul 2019 Jun Wang: allocate arrays for post processing
!    Feb 2022 J. Meng/B. Cui: create interface to run inline post with post_2d_decomp
!
!-----------------------------------------------------------------------
!*** allocate post variables
!-----------------------------------------------------------------------
!
      use vrbls4d
      use vrbls3d
      use vrbls2d
      use soil
      use masks,      only: lmv, lmh, htm, vtm
      use ctlblk_mod, only: im, jm, lm, im_jm, lp1, grib, gdsdegr, me,      &
                            ioform, jsta, jend, jsta_m, jsta_m2, &
                            jend_m, jend_m2, jvend_2u, jsta_2l, jend_2u, iup, idn, &
                            icnt, idsp, mpi_comm_comp, num_servers,     &
                            modelname, numx, ista, iend, ista_m, ista_m2, &
                            iend_m, iend_m2, ista_2l, iend_2u, &
                            ileft,iright,ileftb,irightb, &
                            icnt2, idsp2,isxa,iexa,jsxa,jexa, &
                            num_procs
!
!-----------------------------------------------------------------------
!
      use mpi
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer,intent(in)            :: imi,jmi,lmi,mype,nwtlpes,mpicomp
      integer,intent(in)            :: lead_write
      integer,intent(in)            :: jts,jte
      integer,intent(in)            :: jtsgrp(nwtlpes),jtegrp(nwtlpes)
      integer,intent(in)            :: its,ite
      integer,intent(in)            :: itsgrp(nwtlpes),itegrp(nwtlpes)
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      integer i,j,l
      integer last_write_task
      integer isumm,isumm2
!
!-----------------------------------------------------------------------
!*** get dims from int_state
!-----------------------------------------------------------------------
!
      im = imi
      jm = jmi
      lm = lmi
      im_jm = im*jm
      lp1 = lm + 1
      grib = 'grib2'
! set ndegr
      gdsdegr = 1000000.
      IOFORM = 'grib'
      me = mype-lead_write
      last_write_task = lead_write+nwtlpes-1
      mpi_comm_comp = mpicomp
      num_procs = nwtlpes
      num_servers = 0
      if(mype==0)print *,'grib=',grib,'ioform=',ioform,'mype=',mype,'me=',me, &
         'lead_write=',lead_write,'last_write_task=',last_write_task, &
         'num_servers=',num_servers,'num_procs=',NUM_PROCS,'gdsdegr=',gdsdegr, &
         'im=',im,jm,lm
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE ARRAYS OF THE POST.
!-----------------------------------------------------------------------
!
      jsta = jts
      jend = jte
      jsta_m  = jsta
      jsta_m2 = jsta
      jend_m  = jend
      jend_m2 = jend
      if ( mype<(lead_write+numx) ) then
         jsta_m  = 2
         jsta_m2 = 3
      end if
      if ( mype>(last_write_task-numx) ) then
         jend_m  = jm - 1
         jend_m2 = jm - 2
      end if
      ista = its
      iend = ite
      ista_m  = ista
      ista_m2 = ista
      iend_m  = iend
      iend_m2 = iend
      if(mod(me,numx)==0)then
        ista_m=2
        ista_m2=3
      end if
      if(mod(me+1,numx)==0)then
        iend_m=im-1
        iend_m2=im-2
      end if
!** neighbors
      iup = mype + 1 - lead_write
      idn = mype - 1 - lead_write
      if ( mype == lead_write ) then
         idn = MPI_PROC_NULL
      end if
      if ( mype == last_write_task ) then
         iup = MPI_PROC_NULL
      end if
!      if(mype==0)print *,'lead_write_task=',lead_write,'last taks=',last_write_task, &
!        'idn=',idn,'iup=',iup,'MPI_PROC_NULL=',MPI_PROC_NULL,'jsta=',jsta,'jend=',jend
!
!     counts, disps for gatherv and scatterv
!
      isumm=0
      isumm2=0
      if(allocated(isxa)) deallocate(isxa)
      if(allocated(jsxa)) deallocate(jsxa)
      if(allocated(iexa)) deallocate(iexa)
      if(allocated(jexa)) deallocate(jexa)
      allocate(isxa(0:num_procs-1) )
      allocate(jsxa(0:num_procs-1) )
      allocate(iexa(0:num_procs-1) )
      allocate(jexa(0:num_procs-1) )

      do i = 1, num_procs
       icnt(i-1) = (jtegrp(i)-jtsgrp(i)+1)*(itegrp(i)-itsgrp(i)+1)
       isxa(i-1) = itsgrp(i)
       iexa(i-1) = itegrp(i)
       jsxa(i-1) = jtsgrp(i)
       jexa(i-1) = jtegrp(i)
       idsp(i-1) = isumm
       isumm=isumm+icnt(i-1)
       if(jtsgrp(i)==1 .or. jtegrp(i)==jm) then
          icnt2(i-1) = (itegrp(i)-itsgrp(i)+1)
       else
          icnt2(i-1) = 0
       endif
       idsp2(i-1)=isumm2
       if(jtsgrp(i)==1 .or. jtegrp(i)==jm) isumm2=isumm2+(itegrp(i)-itsgrp(i)+1)

      enddo
!      write(6,'(a25,i4,16i8)') 'JESSE,me,icnt2,idsp2=',me,icnt2(0:num_procs-1),idsp2(0:num_procs-1)
!
!     extraction limits -- set to two rows
!
      jsta_2l = max(jsta - 2,  1 )
      jend_2u = min(jend + 2, jm )
      if(modelname=='GFS') then
        ista_2l=max(ista-2,0)
        iend_2u=min(iend+2,im+1)
      else
        ista_2l = max(ista - 2,  1 )
        iend_2u = min(iend + 2, im )
      endif
! special for c-grid v
      jvend_2u = min(jend + 2, jm+1 )
!      if(mype==0)print *,'im=',im,'jsta_2l=',jsta_2l,'jend_2u=',jend_2u,'lm=',lm
!      print 901,'GWVX mype/me=',mype,me,'im=',im,'jsta   =',jsta   ,'jend   =',jend   ,'lm=',lm
!      print 901,'GWVX mype/me=',mype,me,'im=',im,'jsta_m =',jsta_m ,'jend_m =',jend_m ,'lm=',lm
!      print 901,'GWVX mype/me=',mype,me,'im=',im,'jsta_2l=',jsta_2l,'jend_2u=',jend_2u,'lm=',lm
!      print 901,'GWVX mype/me=',mype,me,'im=',im,'ista   =',ista   ,'iend   =',iend   ,'lm=',lm
!      print 901,'GWVX mype/me=',mype,me,'im=',im,'ista_m =',ista_m ,'iend_m =',iend_m ,'lm=',lm
!      print 901,'GWVX mype/me=',mype,me,'im=',im,'ista_2l=',ista_2l,'iend_2u=',iend_2u,'lm=',lm
!  901 format(a15,2i4,4(1x,a8,i4))
!       NEW neighbors
      ileft = me - 1
      iright = me + 1
       iup=MPI_PROC_NULL
       idn=MPI_PROC_NULL
!    if(mod(me,numx) .eq. 0) print *,' LEFT POINT',mype,me
!    if(mod(me+1,numx) .eq. 0) print *,' RIGHT  POINT',mype,me
    if(mod(me,numx) .eq. 0) ileft=MPI_PROC_NULL
    if(mod(me,numx) .eq. 0) ileftb=me+numx-1
!    if(mod(me,numx) .eq. 0) print *,' GWVX ILEFTB ',ileftb,mype,me,numx
    if(mod(me+1,numx) .eq. 0 .or. me .eq. num_procs-1)  iright=MPI_PROC_NULL
    if(mod(me+1,numx) .eq. 0 .or. me .eq. num_procs-1)  irightb=me-numx+1
!    if(mod(me+1,numx) .eq. 0 .or. me .eq. num_procs-1)  print *,' GWVX IRIGHTB',irightb,mype,me,numx
    if(me .ge. numx) idn=me-numx
    if(me+1  .le. num_procs-numx) iup=me+numx
!    write(6,'(a12,6i10)') 'GWVX BOUNDS ',me,ileft,iright,iup,idn,num_procs
!
! SETS UP MESSAGE PASSING INFO

      call allocate_all()

!***
! LMH always = LM for sigma-type vert coord
! LMV always = LM for sigma-type vert coord

       do j = jsta_2l, jend_2u
        do i = ista_2l, iend_2u
            lmv ( i, j ) = lm
            lmh ( i, j ) = lm
        end do
       end do
!
! HTM VTM all 1 for sigma-type vert coord

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = ista_2l, iend_2u
            htm ( i, j, l ) = 1.0
            vtm ( i, j, l ) = 1.0
        end do
       end do
      end do
    end subroutine post_alctvars
!
!---------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!---------------------------------------------------------------------
!
  subroutine read_postnmlt(kpo,kth,kpv,po,th,pv,post_namelist)
!
      use ctlblk_mod, only : komax,fileNameD3D,lsm,lsmp1,spl,spldef,  &
                             lsmdef,ALSL,me,d3d_on,gocart_on,hyb_sigp,&
                             pthresh,novegtype,ivegsrc,icu_physics,   &
                             isf_surface_physics,modelname,submodelname,&
                             rdaod,d2d_chem,nasa_on,gccpp_on
      use upp_ifi_mod, only: write_ifi_debug_files
!
!    revision history:
!    Jul 2019 Jun Wang: read post namelist
!
      implicit none
!---
      character (len=*), intent(in) :: post_namelist
      integer,intent(out) :: kpo,kth,kpv
      real(4),dimension(komax),intent(out) :: po,th,pv
      integer :: nlunit
      real :: untcnvt
      logical :: popascal
      integer l,k

      namelist/nampgb/kpo,po,kth,th,kpv,pv,popascal,d3d_on,gocart_on,  &
                      hyb_sigp,write_ifi_debug_files,rdaod,nasa_on,gccpp_on, &
                      d2d_chem
      namelist/model_inputs/modelname,submodelname
!---------------------------------------------------------------------
!
!      print *,'in read_postnmlt'
!
! set default for kpo, kth, th, kpv, pv
      kpo = 0
      po  = 0
      kth = 6
      th  = (/310.,320.,350.,450.,550.,650.,(0.,k=kth+1,komax)/) ! isentropic level to output
      kpv = 8
      pv  = (/0.5,-0.5,1.0,-1.0,1.5,-1.5,2.0,-2.0,(0.,k=kpv+1,komax)/)
      hyb_sigp    = .true.
      d3d_on      = .false.
      gocart_on   = .false.
      popascal    = .false.
      rdaod       = .false.
      nasa_on     = .false.
      gccpp_on    = .false.
      d2d_chem    = .false.
!
      if (me == 0) print *,'post_namelist=',post_namelist
!jw post namelist is using the same file itag as standalone post
      open (newunit=nlunit,file=post_namelist,status='old',action='read')
      read (nlunit,model_inputs)
      read (nlunit,nampgb)
      close (nlunit)
      if (me == 0) then
        print*,'komax,kpo,kth,th,kpv,pv,popascal= ',komax,kpo            &
     &  ,kth,th(1:kth),kpv,pv(1:kpv),popascal,' gocart_on=',gocart_on    &
     &  ,' rdaod=',rdaod
       endif
!
! set up pressure level from POSTGPVARS or DEFAULT
      if(kpo == 0)then
! use default pressure levels
        if (me==0) then
          print*,'using default pressure levels,spldef=',(spldef(l),l=1,lsmdef)
        endif
        lsm = lsmdef
        do l=1,lsm
         spl(l) = spldef(l)
        end do
      else
! use POSTGPVARS
        if (me==0) then
          print*,'using pressure levels from POSTGPVARS'
        endif
        lsm = kpo
        if( .not. popascal ) then
          untcnvt = 100.
        else
          untcnvt = 1.
        endif
        if(po(lsm)<po(1))then ! post logic assumes asscending
          do l=1,lsm
            spl(l) = po(lsm-l+1)*untcnvt
          end do
        else
          do l=1,lsm
            spl(l) = po(l)*untcnvt
          end do
        end if
      end if
      lsmp1 = lsm + 1
      pthresh = 0.000001
      if (me==0) print*,'LSM, SPL = ',lsm,spl(1:lsm),' pthresh=', &
        pthresh
!
! set default novegtype for GFS, need to get this variable from gfs physics
      novegtype = 20
      ivegsrc   = 1
!
! set default CU_PHYSICS for GFS: assigned 4 for SAS
      icu_physics = 4
!
! set default GFS LSM physics to 2 for NOAH
      isf_surface_physics = 2
      if (me==0) print*,'set default value, ivegsrc = ',ivegsrc,'novegtype=',novegtype, &
       'icu_physics=',icu_physics,'isf_surface_physics=',isf_surface_physics

!     COMPUTE DERIVED MAP OUTPUT CONSTANTS.
!$omp parallel do private(l)
      do l = 1,lsm
         alsl(l) = log(spl(l))
      enddo
!
1000  continue

      end subroutine read_postnmlt
!
!---------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!---------------------------------------------------------------------
!
    subroutine post_finalize(post_gribversion)
!
!    revision history:
!    Jul 2019 Jun Wang: finalize post step
!
      use grib2_module, only : grib_info_finalize
!
      character(*),intent(in) :: post_gribversion
!
      IF(trim(post_gribversion)=='grib2') then
         ! call  grib_info_finalize()
      ENDIF
!
      call de_allocate
!
    end subroutine post_finalize

