! This program provides unit testing for the subroutines in io/post_nems_routines.F90
!
! Alex Richert, 11 Jan 2024
program test_post_nems_routines

  use ctlblk_mod, only : komax,hyb_sigp,d3d_on,gocart_on, &
   rdaod,nasa_on,gccpp_on,d2d_chem,modelname,submodelname, lsm

  implicit none

  character (len=*), parameter :: post_namelist_empty="data/post_namelist_empty.nml"
  character (len=*), parameter :: post_namelist="data/post_namelist.nml"
  integer :: kpo,kth,kpv
  real(4),dimension(komax) :: po,th,pv
  logical :: popascal
  real, parameter :: tini=tiny(1.0)

  ! Verify default settings by using empty nml file
  call read_postnmlt(kpo,kth,kpv,po,th,pv,trim(post_namelist_empty))
  if (kpo.ne.0) stop 1
  if (kth.ne.6) stop 2
  if (kpv.ne.8) stop 3
  if (any(po.ne.0.0)) stop 4
  if (any(abs(th(1:6)-(/310.,320.,350.,450.,550.,650./)).gt.tini)) stop 5
  if (any(abs(pv(1:8)-(/0.5,-0.5,1.0,-1.0,1.5,-1.5,2.0,-2.0/)).gt.tini)) stop 6
  if (.not.hyb_sigp) stop 7
  if (d3d_on) stop 8
  if (gocart_on) stop 9
  if (lsm.ne.46) stop 10 ! 'lsm' is determined by 'popascal'
  if (rdaod) stop 11
  if (nasa_on) stop 12
  if (gccpp_on) stop 13
  if (d2d_chem) stop 14

  ! Now use fully populated nml file
  call read_postnmlt(kpo,kth,kpv,po,th,pv,trim(post_namelist))
  if (kpo.ne.5) stop 101
  if (kth.ne.7) stop 102
  if (kpv.ne.9) stop 103
  if (po(1).ne.0.5) stop 104
  if (any(po(2:komax).ne.1.0)) stop 104
  if (any(abs(th(1:7)-(/1.,2.,3.,4.,5.,6.,7./)).gt.tini)) stop 105
  if (any(abs(pv(1:9)-(/11.,12.,13.,14.,15.,16.,17.,18.,19./)).gt.tini)) stop 106
  if (hyb_sigp) stop 107
  if (.not.d3d_on) stop 108
  if (.not.gocart_on) stop 109
  if (lsm.ne.5) stop 110 ! 'lsm' is determined by 'popascal'
  if (.not.rdaod) stop 111
  if (.not.nasa_on) stop 112
  if (.not.gccpp_on) stop 113
  if (.not.d2d_chem) stop 114
  if (trim(modelname).ne."DMMY") stop 115
  if (trim(submodelname).ne."SUBM") stop 116

end program test_post_nems_routines
