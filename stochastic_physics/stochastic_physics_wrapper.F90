module stochastic_physics_wrapper_mod

  use machine, only: kind_phys

  implicit none

  ! For stochastic physics pattern generation
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlat
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlon
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: sppt_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: shum_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: skebu_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: skebv_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: sfc_wts
  real(kind=kind_phys), dimension(:,:,:,:), allocatable, save :: spp_wts

  logical, save :: is_initialized = .false.
  integer, save :: lsoil = -999
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: smc
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: stc
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: slc
  !
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: vfrac
  !albedo
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: snoalb
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: alnsf
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: alnwf
  !emissivity
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: semis
  !roughness length for land
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: zorll

  !real(kind=kind_phys), dimension(:,:),   allocatable, save :: stype
  integer, dimension(:,:),   allocatable, save :: stype

  ! For cellular automata
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: sst
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: lmsk
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: lake
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: uwind
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: vwind
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: height
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: dx
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: condition
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca_deep_cpl, ca_turb_cpl, ca_shal_cpl
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca1_cpl, ca2_cpl, ca3_cpl


!----------------
! Public Entities
!----------------
! functions
  public stochastic_physics_wrapper

  contains

!-------------------------------
!  CCPP step
!-------------------------------
  subroutine stochastic_physics_wrapper (GFS_Control, GFS_Statein, GFS_Grid, GFS_Sfcprop, GFS_Radtend, GFS_Coupling, Atm_block, ierr)

#ifdef _OPENMP
    use omp_lib
#endif

    use GFS_typedefs,       only: GFS_control_type, GFS_statein_type, GFS_grid_type, GFS_sfcprop_type, GFS_radtend_type, GFS_coupling_type
    use mpp_mod,            only: FATAL, mpp_error
    use block_control_mod,  only: block_control_type
    use atmosphere_mod,     only: Atm, mygrid

    use stochastic_physics,           only: init_stochastic_physics, run_stochastic_physics
    use cellular_automata_global_mod, only: cellular_automata_global
    use cellular_automata_sgs_mod,    only: cellular_automata_sgs
    use lndp_apply_perts_mod,         only: lndp_apply_perts

    implicit none

    type(GFS_control_type),   intent(inout) :: GFS_Control
    type(GFS_statein_type),   intent(in)    :: GFS_Statein
    type(GFS_grid_type),      intent(in)    :: GFS_Grid
    type(GFS_sfcprop_type),   intent(inout) :: GFS_Sfcprop
    type(GFS_radtend_type),   intent(inout) :: GFS_Radtend
    type(GFS_coupling_type),  intent(inout) :: GFS_Coupling
    type(block_control_type), intent(inout) :: Atm_block
    integer,                  intent(out)   :: ierr

    integer :: nthreads, nb, levs, maxblk, nblks, n, v, ixs, ixe
    logical :: param_update_flag

#ifdef _OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif
    ierr = 0

    levs   = GFS_Control%levs
    maxblk = maxval(GFS_Control%blksz)
    nblks  = Atm_block%nblks

    ! Initialize

    initalize_stochastic_physics: if (.not. is_initialized) then

      if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. (GFS_Control%lndp_type > 0) .OR. GFS_Control%do_spp) then
         allocate(xlat(1:nblks,maxblk))
         allocate(xlon(1:nblks,maxblk))
         call transfer_field_to_stochastics(GFS_Control%blksz, GFS_Grid%xlat, xlat)
         call transfer_field_to_stochastics(GFS_Control%blksz, GFS_Grid%xlon, xlon)
        ! Initialize stochastic physics
        call init_stochastic_physics(levs, GFS_Control%blksz, GFS_Control%dtp, GFS_Control%sppt_amp,                                  &
            GFS_Control%input_nml_file, GFS_Control%fn_nml, GFS_Control%nlunit, xlon, xlat, GFS_Control%do_sppt, GFS_Control%do_shum, &
            GFS_Control%do_skeb, GFS_Control%lndp_type, GFS_Control%n_var_lndp, GFS_Control%use_zmtnblck, GFS_Control%skeb_npass,     &
            GFS_Control%lndp_var_list, GFS_Control%lndp_prt_list,    &
            GFS_Control%n_var_spp, GFS_Control%spp_var_list, GFS_Control%spp_prt_list, GFS_Control%spp_stddev_cutoff, GFS_Control%do_spp,                            &
            GFS_Control%ak, GFS_Control%bk, nthreads, GFS_Control%master, GFS_Control%communicator, ierr)
            if (ierr/=0)  then
                    write(6,*) 'call to init_stochastic_physics failed'
                    return
            endif
      end if
      if (GFS_Control%do_sppt) then
         allocate(sppt_wts(1:nblks,maxblk,1:levs))
      end if
      if (GFS_Control%do_shum) then
         allocate(shum_wts(1:nblks,maxblk,1:levs))
      end if
      if (GFS_Control%do_skeb) then
         allocate(skebu_wts(1:nblks,maxblk,1:levs))
         allocate(skebv_wts(1:nblks,maxblk,1:levs))
      end if
      if ( GFS_Control%do_spp ) then
         allocate(spp_wts(1:nblks,maxblk,1:levs,1:GFS_Control%n_var_spp))
         do n=1,GFS_Control%n_var_spp
           select case (trim(GFS_Control%spp_var_list(n)))
           case('pbl')
             GFS_Control%spp_pbl = 1
           case('sfc')
             GFS_Control%spp_sfc = 1
           case('mp')
             GFS_Control%spp_mp = 7
           case('rad')
             GFS_Control%spp_rad = 1
           case('gwd')
             GFS_Control%spp_gwd = 1
           case('cu_deep')
             GFS_Control%spp_cu_deep = 1
           end select
         end do
      end if
      if ( GFS_Control%lndp_type == 2 ) then
          allocate(sfc_wts(1:nblks,maxblk,1:GFS_Control%n_var_lndp))
          if ( (GFS_Control%lsm == GFS_Control%lsm_noah) .or. (GFS_Control%lsm == GFS_Control%lsm_noahmp)) then
            lsoil = GFS_Control%lsoil
          elseif (GFS_Control%lsm == GFS_Control%lsm_ruc) then
            lsoil = GFS_Control%lsoil_lsm
          endif
          allocate(smc   (1:nblks, maxblk, lsoil))
          do v = 1,GFS_Control%n_var_lndp
            select case (trim(GFS_Control%lndp_var_list(v)))
            case('smc')
              allocate(slc   (1:nblks, maxblk, lsoil))
              allocate(stype (1:nblks, maxblk))
            case('stc')
              allocate(stc   (1:nblks, maxblk, lsoil))
            case('vgf') 
              allocate(vfrac (1:nblks, maxblk))
            case('alb') 
              allocate(alnsf (1:nblks, maxblk))
              allocate(alnwf (1:nblks, maxblk))
            case('sal') 
              allocate(snoalb(1:nblks, maxblk))
            case('emi') 
              allocate(semis (1:nblks, maxblk))
            case('zol') 
              allocate(zorll (1:nblks, maxblk))
            endselect 
          enddo 
      endif


      if ( GFS_Control%lndp_type == 1 ) then ! this scheme sets perts once
         allocate(sfc_wts(1:nblks, maxblk, GFS_Control%n_var_lndp))
         call run_stochastic_physics(levs, GFS_Control%kdt, GFS_Control%fhour, GFS_Control%blksz,       &
                                     sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts,         &
                                     skebv_wts=skebv_wts, sfc_wts=sfc_wts,                              &
                                     spp_wts=spp_wts, nthreads=nthreads)
         ! Copy contiguous data back
         do nb=1,nblks
            GFS_Coupling%sfc_wts(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = sfc_wts(nb,1:GFS_Control%blksz(nb),:)
         end do
         deallocate(sfc_wts)
      end if
      ! Consistency check for cellular automata
      if(GFS_Control%do_ca)then
        if(GFS_Control%ca_sgs)then
           allocate(sst         (1:nblks, maxblk))
           allocate(lmsk        (1:nblks, maxblk))
           allocate(lake        (1:nblks, maxblk))
           allocate(uwind       (1:nblks, maxblk, 1:levs))
           allocate(vwind       (1:nblks, maxblk, 1:levs))
           allocate(height      (1:nblks, maxblk, 1:levs))
           allocate(condition   (1:nblks, maxblk))
           allocate(dx          (1:nblks, maxblk))
           allocate(ca_deep_cpl (1:nblks, maxblk))
           allocate(ca_turb_cpl (1:nblks, maxblk))
           allocate(ca_shal_cpl (1:nblks, maxblk))
        endif
        if(GFS_Control%ca_global)then
          ! Allocate contiguous arrays; no need to copy in (intent out)
          allocate(ca1_cpl (1:nblks, maxblk))
          allocate(ca2_cpl (1:nblks, maxblk))
          allocate(ca3_cpl (1:nblks, maxblk))
        endif
      endif

      is_initialized = .true.

    else initalize_stochastic_physics
      if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. (GFS_Control%lndp_type == 2) .OR. GFS_Control%do_spp) then
         call run_stochastic_physics(levs, GFS_Control%kdt, GFS_Control%fhour, GFS_Control%blksz, &
                                 sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts, skebv_wts=skebv_wts, sfc_wts=sfc_wts, &
                                 spp_wts=spp_wts, nthreads=nthreads)
         ! Copy contiguous data back
         if (GFS_Control%do_sppt) then
            do nb=1,nblks
                GFS_Coupling%sppt_wts(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = sppt_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
         end if
         if (GFS_Control%do_shum) then
            do nb=1,nblks
                GFS_Coupling%shum_wts(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = shum_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
         end if
         if (GFS_Control%do_skeb) then
            do nb=1,nblks
                GFS_Coupling%skebu_wts(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = skebu_wts(nb,1:GFS_Control%blksz(nb),:)
                GFS_Coupling%skebv_wts(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = skebv_wts(nb,1:GFS_Control%blksz(nb),:)
            end do
         end if
         if (GFS_Control%do_spp) then
            do n=1,GFS_Control%n_var_spp
               select case (trim(GFS_Control%spp_var_list(n)))
               case('pbl')
                 do nb=1,Atm_block%nblks
                     GFS_Coupling%spp_wts_pbl(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('sfc')
                 do nb=1,Atm_block%nblks
                     GFS_Coupling%spp_wts_sfc(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('mp')
                 do nb=1,Atm_block%nblks
                     GFS_Coupling%spp_wts_mp(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('gwd')
                 do nb=1,Atm_block%nblks
                     GFS_Coupling%spp_wts_gwd(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('rad')
                 do nb=1,Atm_block%nblks
                     GFS_Coupling%spp_wts_rad(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               case('cu_deep')
                 do nb=1,Atm_block%nblks
                     GFS_Coupling%spp_wts_cu_deep(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = spp_wts(nb,1:GFS_Control%blksz(nb),:,n)
                 end do
               end select
            end do
         end if

         if (GFS_Control%lndp_type == 2) then ! save wts, and apply lndp scheme
             do nb=1,nblks
                GFS_Coupling%sfc_wts(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb),:) = sfc_wts(nb,1:GFS_Control%blksz(nb),:)
             end do
 
             do nb=1,nblks
                ixs = GFS_control%chunk_begin(nb)
                ixe = GFS_control%chunk_end(nb)
                do v = 1,GFS_Control%n_var_lndp
                  ! used to identify locations with land model (=soil) 
                  if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                     smc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Sfcprop%smois(ixs:ixe,1:lsoil)
                  else  ! noah or noah-MP
                     smc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Sfcprop%smc(ixs:ixe,1:lsoil)
                  endif

                  select case (trim(GFS_Control%lndp_var_list(v)))
                  ! DH* is this correct? shouldn't this be slc ?
                  case('smc')
                      ! stype used to fetch soil params
                      stype(nb,1:GFS_Control%blksz(nb))  = GFS_Sfcprop%stype(ixs:ixe)
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                         slc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Sfcprop%sh2o(ixs:ixe,1:lsoil)
                      else  ! noah or noah-MP
                         slc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Sfcprop%slc(ixs:ixe,1:lsoil)
                      endif
                  case('stc')
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                         stc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Sfcprop%tslb(ixs:ixe,1:lsoil)
                      else ! noah or noah-MP 
                         stc(nb,1:GFS_Control%blksz(nb),1:lsoil) = GFS_Sfcprop%stc(ixs:ixe,1:lsoil)
                      endif 
                  case('vgf')
                      if ( (GFS_Control%lsm == GFS_Control%lsm_noahmp) ) then 
                         ! assumes iopt_dveg = 4 (will be checked later)
                         vfrac(nb,1:GFS_Control%blksz(nb))  = GFS_Sfcprop%shdmax(ixs:ixe)
                      else ! ruc or noah-MP
                         vfrac(nb,1:GFS_Control%blksz(nb))  = GFS_Sfcprop%vfrac(ixs:ixe)
                      endif
                  case('alb')
                      alnsf(nb,1:GFS_Control%blksz(nb))  = GFS_Sfcprop%alnsf(ixs:ixe)
                      alnwf(nb,1:GFS_Control%blksz(nb))  = GFS_Sfcprop%alnwf(ixs:ixe)
                  case('sal')
                      snoalb(nb,1:GFS_Control%blksz(nb)) = GFS_Sfcprop%snoalb(ixs:ixe)
                  case('emi')
                      semis(nb,1:GFS_Control%blksz(nb))  = GFS_Radtend%semis(ixs:ixe)
                  case('zol')
                      zorll(nb,1:GFS_Control%blksz(nb))  = GFS_Sfcprop%zorll(ixs:ixe)
                  endselect
              enddo
             enddo


             param_update_flag = .false.
             ! noah and noah-MP treated differently, as global cycle doesn't overwrite shdmax for Noah-MP 
             ! determine whether land paramaters have been over-written to
             ! trigger applying perturbations (logic copied from GFS_driver),
             if ( (GFS_Control%lsm == GFS_Control%lsm_noah) .and. GFS_Control%nscyc >  0) then
                 if (mod(GFS_Control%kdt,GFS_Control%nscyc) == 1 ) then
                   param_update_flag = .true.
                 endif
             endif
             if ( ( GFS_Control%nscyc ==  0 .or. GFS_Control%lsm == GFS_Control%lsm_noahmp) .and. GFS_Control%first_time_step ) then
             ! call once at start of the forecast.
                    param_update_flag = .true.
             endif
              
             call lndp_apply_perts(GFS_Control%blksz, GFS_Control%lsm, GFS_Control%lsm_noah, GFS_Control%lsm_ruc,             &
                               GFS_Control%lsm_noahmp, GFS_Control%iopt_dveg, lsoil, GFS_Control%dtp, GFS_Control%kdt,        &
                               GFS_Control%n_var_lndp, GFS_Control%lndp_var_list, GFS_Control%lndp_prt_list,                  &
                               sfc_wts, xlon, xlat, stype, GFS_Control%pores, GFS_Control%resid,param_update_flag,            &
                               smc, slc, stc, vfrac, alnsf, alnwf, snoalb, semis, zorll, ierr)

             if (ierr/=0)  then
                    write(6,*) 'call to GFS_apply_lndp failed'
                    return
             endif

             do nb=1,nblks
                ixs = GFS_control%chunk_begin(nb)
                ixe = GFS_control%chunk_end(nb)

                do v = 1,GFS_Control%n_var_lndp

                  select case (trim(GFS_Control%lndp_var_list(v)))
                  case('smc')
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                           GFS_Sfcprop%smois(ixs:ixe,1:lsoil) = smc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                           GFS_Sfcprop%sh2o(ixs:ixe,1:lsoil)  = slc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      else  ! noah or noah-MP
                           GFS_Sfcprop%smc(ixs:ixe,1:lsoil) = smc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                           GFS_Sfcprop%slc(ixs:ixe,1:lsoil) = slc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      endif
                  case('stc')
                      if ((GFS_Control%lsm == GFS_Control%lsm_ruc) ) then 
                           GFS_Sfcprop%tslb(ixs:ixe,1:lsoil)  = stc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      else ! noah or noah-MP 
                           GFS_Sfcprop%stc(ixs:ixe,1:lsoil) = stc(nb,1:GFS_Control%blksz(nb),1:lsoil)
                      endif 
                  case('vgf')
                      if ( (GFS_Control%lsm == GFS_Control%lsm_noahmp) ) then 
                        GFS_Sfcprop%shdmax(ixs:ixe)  = vfrac(nb,1:GFS_Control%blksz(nb))
                      else 
                        GFS_Sfcprop%vfrac(ixs:ixe)  = vfrac(nb,1:GFS_Control%blksz(nb))
                      endif
                  case('alb')
                      GFS_Sfcprop%alnsf(ixs:ixe)  = alnsf(nb,1:GFS_Control%blksz(nb))
                      GFS_Sfcprop%alnwf(ixs:ixe)  = alnwf(nb,1:GFS_Control%blksz(nb))
                  case('sal')
                      GFS_Sfcprop%snoalb(ixs:ixe) = snoalb(nb,1:GFS_Control%blksz(nb))
                  case('emi')
                      GFS_Radtend%semis(ixs:ixe)  = semis(nb,1:GFS_Control%blksz(nb))
                  case('zol')
                      GFS_Sfcprop%zorll(ixs:ixe)  = zorll(nb,1:GFS_Control%blksz(nb))
                  end select
                enddo 
            enddo
         endif ! lndp block
      endif ! if do* block

      if (GFS_Control%do_ca) then

       if(GFS_Control%ca_sgs)then
         do nb=1,nblks
             ixs = GFS_control%chunk_begin(nb)
             ixe = GFS_control%chunk_end(nb)
             sst        (nb,1:GFS_Control%blksz(nb)) = GFS_Sfcprop%tsfco(ixs:ixe)
             lmsk       (nb,1:GFS_Control%blksz(nb)) = GFS_Sfcprop%slmsk(ixs:ixe)
             lake       (nb,1:GFS_Control%blksz(nb)) = GFS_Sfcprop%lakefrac(ixs:ixe)
             condition  (nb,1:GFS_Control%blksz(nb)) = GFS_Coupling%condition(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb))
             ca_deep_cpl(nb,1:GFS_Control%blksz(nb)) = GFS_Coupling%ca_deep(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb))
             ca_turb_cpl(nb,1:GFS_Control%blksz(nb)) = GFS_Coupling%ca_turb(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb))
             ca_shal_cpl(nb,1:GFS_Control%blksz(nb)) = GFS_Coupling%ca_shal(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb))
          enddo
          call transfer_field_to_stochastics_3d(GFS_Control%blksz, GFS_Statein%ugrs, uwind)
          call transfer_field_to_stochastics_3d(GFS_Control%blksz, GFS_Statein%vgrs, vwind)
          call transfer_field_to_stochastics_3d(GFS_Control%blksz, GFS_Statein%phil, height)
         call transfer_field_to_stochastics(GFS_Control%blksz, GFS_Grid%dx, dx)
         call cellular_automata_sgs(GFS_Control%kdt,GFS_control%dtp,GFS_control%restart,GFS_Control%first_time_step,              &
            sst,lmsk,lake,uwind,vwind,height,dx,condition,ca_deep_cpl,ca_turb_cpl,ca_shal_cpl, Atm(mygrid)%domain_for_coupler,nblks,      &
            Atm_block%isc,Atm_block%iec,Atm_block%jsc,Atm_block%jec,Atm(mygrid)%npx,Atm(mygrid)%npy, levs,                        &
            GFS_Control%nthresh,GFS_Control%tile_num,GFS_Control%nca,GFS_Control%ncells,GFS_Control%nlives,                       &
            GFS_Control%nfracseed, GFS_Control%nseed,GFS_Control%iseed_ca,GFS_Control%ca_advect,                                  &
            GFS_Control%nspinup,GFS_Control%ca_trigger,Atm_block%blksz(1),GFS_Control%master,GFS_Control%communicator)
         ! Copy contiguous data back as needed
         do nb=1,nblks
             GFS_Coupling%ca_deep(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb)) = ca_deep_cpl (nb,1:GFS_Control%blksz(nb))
             GFS_Coupling%ca_turb(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb)) = ca_turb_cpl (nb,1:GFS_Control%blksz(nb))
             GFS_Coupling%ca_shal(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb)) = ca_shal_cpl (nb,1:GFS_Control%blksz(nb))
         enddo
       endif
       if(GFS_Control%ca_global)then
          call cellular_automata_global(GFS_Control%kdt,GFS_control%restart,GFS_Control%first_time_step,ca1_cpl,ca2_cpl,ca3_cpl,                &
            Atm(mygrid)%domain_for_coupler, nblks,Atm_block%isc,Atm_block%iec,Atm_block%jsc,Atm_block%jec,Atm(mygrid)%npx,Atm(mygrid)%npy,levs, &
            GFS_Control%nca_g,GFS_Control%ncells_g,GFS_Control%nlives_g,GFS_Control%nfracseed,GFS_Control%nseed_g,                              &
            GFS_Control%iseed_ca,GFS_control%tile_num,GFS_Control%ca_smooth,GFS_Control%nspinup,Atm_block%blksz(1),                             &
            GFS_Control%nsmooth,GFS_Control%ca_amplitude,GFS_Control%master,GFS_Control%communicator)
          ! Copy contiguous data back
          do nb=1,nblks
             GFS_Coupling%ca1(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb)) = ca1_cpl(nb,1:GFS_Control%blksz(nb))
             GFS_Coupling%ca2(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb)) = ca2_cpl(nb,1:GFS_Control%blksz(nb))
             GFS_Coupling%ca3(GFS_Control%chunk_begin(nb):GFS_Control%chunk_end(nb)) = ca3_cpl(nb,1:GFS_Control%blksz(nb))
          enddo
       endif

      endif !do_ca

    endif initalize_stochastic_physics

  contains

    subroutine transfer_field_to_stochastics(blksz, data_in, data_out)

      integer, dimension(:), intent(in) :: blksz
      real(kind=kind_phys), dimension(:), intent(in) :: data_in
      real(kind=kind_phys), dimension(:,:), intent(out) :: data_out
      integer :: i, nb, ni

      nb = 1
      ni = 1
      do i=1,size(data_in)
        if (ni>blksz(nb)) then
          nb = nb+1
          ni = 1
        end if
        data_out(nb,ni) = data_in(i)
        ni =  ni+1
      end do

    end subroutine transfer_field_to_stochastics

    subroutine transfer_field_to_stochastics_3d(blksz, data_in, data_out)

      integer, dimension(:), intent(in) :: blksz
      real(kind=kind_phys), dimension(:,:), intent(in) :: data_in
      real(kind=kind_phys), dimension(:,:,:), intent(out) :: data_out
      integer :: j

      do j=1,size(data_in, dim=2)
         call transfer_field_to_stochastics(blksz, data_in(:,j), data_out(:,:,j))
      end do

    end subroutine transfer_field_to_stochastics_3d

    subroutine transfer_field_from_stochastics(blksz, data_in, data_out)

      integer, dimension(:), intent(in) :: blksz
      real(kind=kind_phys), dimension(:,:), intent(in) :: data_in
      real(kind=kind_phys), dimension(:), intent(out) :: data_out
      integer :: i, nb, ni

      nb = 1
      ni = 1
      do i=1,size(data_out)
        if (ni>blksz(nb)) then
          nb = nb+1
          ni = 1
        end if
        data_out(i)= data_in(nb,ni)
        ni =  ni+1
      end do

    end subroutine transfer_field_from_stochastics

  end subroutine stochastic_physics_wrapper


  subroutine stochastic_physics_wrapper_end (GFS_Control)

  use GFS_typedefs,       only: GFS_control_type
  use stochastic_physics, only: finalize_stochastic_physics

  implicit none

  type(GFS_control_type),   intent(inout) :: GFS_Control

  if (GFS_Control%do_sppt .OR. GFS_Control%do_shum .OR. GFS_Control%do_skeb .OR. (GFS_Control%lndp_type > 0) .OR. GFS_Control%do_spp) then
      if (allocated(xlat)) deallocate(xlat)
      if (allocated(xlon)) deallocate(xlon)
      if (GFS_Control%do_sppt) then
         if (allocated(sppt_wts)) deallocate(sppt_wts)
      end if
      if (GFS_Control%do_shum) then
         if (allocated(shum_wts)) deallocate(shum_wts)
      end if
      if (GFS_Control%do_skeb) then
         if (allocated(skebu_wts)) deallocate(skebu_wts)
         if (allocated(skebv_wts)) deallocate(skebv_wts)
      end if
      if (GFS_Control%do_spp) then
         if (allocated(spp_wts)) deallocate(spp_wts)
      end if
      if ( GFS_Control%lndp_type == 2 ) then
         lsoil = -999
         if (allocated(sfc_wts)) deallocate(sfc_wts)
      end if
      if (GFS_Control%lndp_type == 2) then
          if (allocated(smc))    deallocate(smc)
          if (allocated(slc))    deallocate(slc)
          if (allocated(stc))    deallocate(stc)
          if (allocated(stype))  deallocate(stype)
          if (allocated(vfrac))  deallocate(vfrac)
          if (allocated(snoalb)) deallocate(snoalb)
          if (allocated(alnsf))  deallocate(alnsf)
          if (allocated(alnwf))  deallocate(alnwf)
          if (allocated(semis))  deallocate(semis)
          if (allocated(zorll))  deallocate(zorll)
      endif
      call finalize_stochastic_physics()
   endif
   if(GFS_Control%do_ca)then
        if(GFS_Control%ca_sgs)then
           deallocate(sst         )
           deallocate(lmsk        )
           deallocate(lake        )
           deallocate(uwind       )
           deallocate(vwind       )
           deallocate(height      )
           deallocate(dx          )
           deallocate(condition   )
           deallocate(ca_deep_cpl )
           deallocate(ca_turb_cpl )
           deallocate(ca_shal_cpl )
        endif
        if(GFS_Control%ca_global)then
            deallocate(ca1_cpl )
            deallocate(ca2_cpl )
            deallocate(ca3_cpl )
        endif
   endif
  end subroutine stochastic_physics_wrapper_end

end module stochastic_physics_wrapper_mod
