cmake_minimum_required(VERSION 3.19)

project(CCPP-UFS
        LANGUAGES C CXX Fortran)
set(PROJECT "CCPP-UFS")

#------------------------------------------------------------------------------
#
# Create list of SCHEME_FILES, HOST_FILES, and SUITE_FILES.
# Paths should be relative to CMAKE_SOURCE_DIR (this file's directory).
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# SCHEME_FILES
#------------------------------------------------------------------------------
set(SCHEME_FILES_DIR "${CMAKE_SOURCE_DIR}/${PARENT_DIR}/ccpp/physics/physics")
set(SCHEME_METADATA_FILES
  # time-vary
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmg_setup.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_time_vary_pre.fv3.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_phys_time_vary.fv3.meta"
  # radiation
  "${SCHEME_FILES_DIR}/Radiation/RRTMG/radsw_param.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmg_pre.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_radiation_surface.meta"
  "${SCHEME_FILES_DIR}/Radiation/RRTMG/rad_sw_pre.meta"
  "${SCHEME_FILES_DIR}/Radiation/RRTMG/radsw_main.meta"
  "${SCHEME_FILES_DIR}/Radiation/RRTMG/rrtmg_sw_post.meta"
  "${SCHEME_FILES_DIR}/Radiation/RRTMG/radlw_main.meta"
  "${SCHEME_FILES_DIR}/Radiation/RRTMG/rrtmg_lw_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_radiation_post.meta"
  # physics
  "${SCHEME_FILES_DIR}/tools/get_prs_fv3.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_1.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_generic_pre.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_composites_pre.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/dcyc2t3.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_composites_inter.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_2.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_loop_control_part1.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_diff.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/scm_sfc_flux_spec.meta"
  "${SCHEME_FILES_DIR}/SFC_Models/Land/Noahmp/lnd_iau_mod.meta"
  "${SCHEME_FILES_DIR}/SFC_Models/Land/sfc_land.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_nst_pre.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_nst.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_nst_post.meta"
  "${SCHEME_FILES_DIR}/SFC_Models/Land/Noahmp/noahmpdrv.meta"
  "${SCHEME_FILES_DIR}/SFC_Models/SeaIce/CICE/sfc_sice.meta"
  "${SCHEME_FILES_DIR}/SFC_Models/SeaIce/CICE/sfc_cice.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_loop_control_part2.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_composites_post.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_diag.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_diag_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_generic_post.meta"
  "${SCHEME_FILES_DIR}/PBL/SATMEDMF/satmedmfvdifq.meta"
  "${SCHEME_FILES_DIR}/PBL/SATMEDMF/satmedmfvdif.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_PBL_generic_pre.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_PBL_generic_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_GWD_generic_pre.meta"
  "${SCHEME_FILES_DIR}/GWD/unified_ugwp.meta"
  "${SCHEME_FILES_DIR}/GWD/unified_ugwp_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_GWD_generic_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_photochemistry.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_photochemistry_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_stochastics.meta"
  "${SCHEME_FILES_DIR}/tools/get_phi_fv3.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_3.meta"
  "${SCHEME_FILES_DIR}/CONV/SAMF/samfdeepcnv.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_DCNV_generic_post.meta"
  "${SCHEME_FILES_DIR}/CONV/SAMF/samfshalcnv.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_SCNV_generic_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_4.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/cnvc90.meta"
  "${SCHEME_FILES_DIR}/MP/Thompson/mp_thompson.meta"
  "${SCHEME_FILES_DIR}/MP/Thompson/mp_thompson_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_MP_generic_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/maximum_hourly_diagnostics.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_physics_post.meta"
  "${SCHEME_FILES_DIR}/CONV/C3/CCPP/cu_c3_driver_ccpp.meta"
  "${SCHEME_FILES_DIR}/CONV/C3/CCPP/cu_c3_driver_post.meta"
  "${SCHEME_FILES_DIR}/CONV/C3/CCPP/cu_c3_driver_pre.meta"
  "${SCHEME_FILES_DIR}/GWD/ugwpv1_gsldrag.meta"
  "${SCHEME_FILES_DIR}/GWD/ugwpv1_gsldrag_post.meta"
  "${SCHEME_FILES_DIR}/SFC_Models/Ocean/UFS/sfc_ocean.meta"
  "${SCHEME_FILES_DIR}/MP/TEMPO/mp_tempo_post.meta"
  "${SCHEME_FILES_DIR}/MP/TEMPO/mp_tempo.meta"
  "${SCHEME_FILES_DIR}/MP/TEMPO/mp_tempo_condensation.meta"
  "${SCHEME_FILES_DIR}/MP/TEMPO/tempo_v3/src/module_mp_tempo_cfgs.meta"
  "${SCHEME_FILES_DIR}/SFC_Models/Land/Noah/lsm_noah.meta"
  "${SCHEME_FILES_DIR}/GWD/cires_ugwp.meta"
  "${SCHEME_FILES_DIR}/GWD/cires_ugwp_post.meta"
  "${SCHEME_FILES_DIR}/MP/GFDL/v1_2019/gfdl_cloud_microphys.meta"
  "${SCHEME_FILES_DIR}/MP/GFDL/fv_sat_adj.meta"
  "${SCHEME_FILES_DIR}/MP/GFDL/v3_2022/gfdl_cloud_microphys_v3.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_cloud_overlap.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_setup.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_pre.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_cloud_mp.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_cloud_diagnostics.meta"
  "${SCHEME_FILES_DIR}/Radiation/RRTMGP/rrtmgp_aerosol_optics.meta"
  "${SCHEME_FILES_DIR}/Radiation/RRTMGP/rrtmgp_sw_main.meta"
  "${SCHEME_FILES_DIR}/Radiation/RRTMGP/rrtmgp_lw_main.meta"
  "${SCHEME_FILES_DIR}/PBL/HEDMF/hedmf.meta"
  "${SCHEME_FILES_DIR}/MP/Ferrier_Aligo/mp_fer_hires.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/GFDL/gfdl_sfc_layer.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/sgscloud_radpre.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/sgscloud_radpost.meta"
  "${SCHEME_FILES_DIR}/PBL/MYNN_EDMF/mynnedmf_wrapper.meta"
  "${SCHEME_FILES_DIR}/PBL/MYNN_EDMF/mynnedmf_wrapper_post.meta"
  "${SCHEME_FILES_DIR}/CONV/Grell_Freitas/cu_gf_driver_pre.meta"
  "${SCHEME_FILES_DIR}/CONV/Grell_Freitas/cu_gf_driver.meta"
  "${SCHEME_FILES_DIR}/CONV/Grell_Freitas/cu_gf_driver_post.meta"
  "${SCHEME_FILES_DIR}/MP/NSSL/mp_nssl.meta"
  "${SCHEME_FILES_DIR}/SFC_Models/Land/RUC/lsm_ruc.meta"
  "${SCHEME_FILES_DIR}/smoke_dust/rrfs_smoke_wrapper.meta"
  "${SCHEME_FILES_DIR}/smoke_dust/rrfs_smoke_postpbl.meta"
  "${SCHEME_FILES_DIR}/GWD/drag_suite.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/MYNN/MYNN/CCPP/module_sf_mynnsfc_driver.meta"
  "${SCHEME_FILES_DIR}/SFC_Models/Lake/CLM/clm_lake.meta"
  "${SCHEME_FILES_DIR}/SFC_Models/Lake/Flake/flake_driver.meta"
  "${SCHEME_FILES_DIR}/CONV/Chikira_Sugiyama/cs_conv_pre.meta"
  "${SCHEME_FILES_DIR}/CONV/Chikira_Sugiyama/cs_conv.meta"
  "${SCHEME_FILES_DIR}/CONV/Chikira_Sugiyama/cs_conv_post.meta"
  "${SCHEME_FILES_DIR}/MP/Morrison_Gettelman/m_micro_pre.meta"
  "${SCHEME_FILES_DIR}/MP/Morrison_Gettelman/m_micro.meta"
  "${SCHEME_FILES_DIR}/MP/Morrison_Gettelman/m_micro_post.meta"
  "${SCHEME_FILES_DIR}/CONV/Chikira_Sugiyama/cs_conv_aw_adj.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/MYJ/myjsfc_wrapper.meta"
  "${SCHEME_FILES_DIR}/PBL/MYJ/myjpbl_wrapper.meta"
  "${SCHEME_FILES_DIR}/CONV/nTiedtke/cu_ntiedtke.meta"
  "${SCHEME_FILES_DIR}/CONV/nTiedtke/cu_ntiedtke_pre.meta"
  "${SCHEME_FILES_DIR}/CONV/nTiedtke/cu_ntiedtke_post.meta"
  "${SCHEME_FILES_DIR}/PBL/saYSU/shinhongvdif.meta"
  "${SCHEME_FILES_DIR}/PBL/YSU/ysuvdif.meta"
  "${SCHEME_FILES_DIR}/PBL/SATMEDMF/canopy_driver.meta"
)

# Dycore (FV3) specific interstitial files.
if (FV3)
  list(APPEND SCHEME_METADATA_FILES
    "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_rad_time_vary.fv3.meta"
  )
endif()

# Dycore (MPAS) specific interstitial files.
if (MPAS)
  list(APPEND SCHEME_METADATA_FILES
    "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_rad_time_vary.mpas.meta"
  )
endif()

# Create list of scheme source files from metadata files
ccpp_source_files(SCHEME_FORTRAN_FILES ${SCHEME_METADATA_FILES})

#------------------------------------------------------------------------------
# HOST_FILES
#------------------------------------------------------------------------------
set(HOST_FILES_DIR "${CMAKE_SOURCE_DIR}/${PARENT_DIR}/ccpp")
set(HOST_METADATA_FILES
   "${HOST_FILES_DIR}/data/CCPP_typedefs.meta"
   "${HOST_FILES_DIR}/data/GFS_typedefs.meta"
   "${HOST_FILES_DIR}/CCPP_driver.meta"
   "${HOST_FILES_DIR}/data/CCPP_data.meta"
   "physics/physics/hooks/physcons.meta"
)

# Create list of host source files from metadata files
ccpp_source_files(HOST_FORTRAN_FILES ${HOST_METADATA_FILES})
message(INFO "${HOST_METADATA_FILES} --> ${HOST_FORTRAN_FILES}")

set(HOST "ufs")

if (FV3)
  message(STATUS "Build CCPP interface to FV3 dynamical core")
  set(CCPP_TARGET fv3ccpp)
  add_definitions(-DFV3)
endif()

if (MPAS)
  message(STATUS "Build CCPP interface to MPAS dynamical core")
  set(CCPP_TARGET mpasccpp)
  add_definitions(-DMPAS)
endif()

#------------------------------------------------------------------------------
# SUITE_FILES
#------------------------------------------------------------------------------
set(SUITE_FILES_DIR "${CMAKE_SOURCE_DIR}/${PARENT_DIR}/ccpp/suites")
if(DEFINED CCPP_SUITES)
  # Take comma seperated string containing SDF names <CCPP_SUITES>,
  # append <SUITE_FILES_DIR> as file prefix, append ".xml" file suffix,
  # and replace "," with ";"
  set(SUITE_FILES ${SUITE_FILES_DIR}/suite_${CCPP_SUITES})
  string(REPLACE "," ".xml;${SUITE_FILES_DIR}/suite_" SUITE_FILES "${SUITE_FILES}")
  set(SUITE_FILES "${SUITE_FILES}.xml")
else()
  set(SUITE_FILES
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v15_thompson_mynn_lam3km.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v15p2.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v16.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v16_csawmg.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v16_flake.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v16_fv3wam.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v16_gfdlmpv3.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v17_coupled_p8.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v17_coupled_p8_c3.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v17_coupled_p8_sfcocn.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v17_coupled_p8_ugwpv1.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v17_p8.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v17_p8_c3.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v17_p8_mynn.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v17_p8_rrtmgp.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v17_p8_ugwpv1.xml"
    "${SUITE_FILES_DIR}/suite_FV3_GFS_v17_p8_ugwpv1_tempo.xml"
    "${SUITE_FILES_DIR}/suite_FV3_HAFS_v1_gfdlmp_tedmf.xml"
    "${SUITE_FILES_DIR}/suite_FV3_HAFS_v1_gfdlmp_tedmf_nonsst.xml"
    "${SUITE_FILES_DIR}/suite_FV3_HAFS_v1_gfdlmp_tedmf_ntiedtke.xml"
    "${SUITE_FILES_DIR}/suite_FV3_HAFS_v1_gfdlmpv3_tedmf.xml"
    "${SUITE_FILES_DIR}/suite_FV3_HAFS_v1_thompson.xml"
    "${SUITE_FILES_DIR}/suite_FV3_HAFS_v1_thompson_nonsst.xml"
    "${SUITE_FILES_DIR}/suite_FV3_HAFS_v1_thompson_tedmf_gfdlsf.xml"
    "${SUITE_FILES_DIR}/suite_FV3_HAFS_v2.xml"
    "${SUITE_FILES_DIR}/suite_FV3_HAFS_v2_coupled.xml"
    "${SUITE_FILES_DIR}/suite_FV3_HRRR.xml"
    "${SUITE_FILES_DIR}/suite_FV3_HRRR_c3.xml"
    "${SUITE_FILES_DIR}/suite_FV3_HRRR_gf.xml"
    "${SUITE_FILES_DIR}/suite_FV3_HRRR_gf_nogwd.xml"
    "${SUITE_FILES_DIR}/suite_FV3_RAP.xml"
    "${SUITE_FILES_DIR}/suite_FV3_RAP_cires_ugwp.xml"
    "${SUITE_FILES_DIR}/suite_FV3_RAP_clm_lake.xml"
    "${SUITE_FILES_DIR}/suite_FV3_RAP_flake.xml"
    "${SUITE_FILES_DIR}/suite_FV3_RAP_noah.xml"
    "${SUITE_FILES_DIR}/suite_FV3_RAP_noah_sfcdiff_cires_ugwp.xml"
    "${SUITE_FILES_DIR}/suite_FV3_RAP_sfcdiff.xml"
    "${SUITE_FILES_DIR}/suite_FV3_RAP_unified_ugwp.xml"
    "${SUITE_FILES_DIR}/suite_FV3_RRFS_v1beta.xml"
    "${SUITE_FILES_DIR}/suite_FV3_RRFS_v1nssl.xml"
    "${SUITE_FILES_DIR}/suite_FV3_WoFS_v0.xml"
    "${SUITE_FILES_DIR}/suite_FV3_coupled_lowres.xml"
    "${SUITE_FILES_DIR}/suite_FV3_global_nest_v1.xml"
    "${SUITE_FILES_DIR}/suite_FV3_ideal_mp_nssl.xml"
    "${SUITE_FILES_DIR}/suite_FV3_ideal_pbl_mp_nssl.xml"
    "${SUITE_FILES_DIR}/suite_FV3_lowres.xml"
    "${SUITE_FILES_DIR}/suite_MPAS_RRFS.xml"
    "${SUITE_FILES_DIR}/suite_RRFS_sas.xml"
    "${SUITE_FILES_DIR}/suite_RRFS_sas_nogwd.xml"
    "${SUITE_FILES_DIR}/suite_RRFSens_phy1.xml"
    "${SUITE_FILES_DIR}/suite_RRFSens_phy2.xml"
    "${SUITE_FILES_DIR}/suite_RRFSens_phy3.xml"
    "${SUITE_FILES_DIR}/suite_RRFSens_phy4.xml"
    "${SUITE_FILES_DIR}/suite_RRFSens_phy5.xml"
  )
endif()
message(STATUS "Compiling the following SDFs: ${SUITE_FILES}")

set(CCPP_VERBOSITY "2" CACHE STRING "Verbosity level of output (default: 0)")

# By default, generated caps go in ccpp subdir
set(OUTPUT_ROOT "${CMAKE_CURRENT_BINARY_DIR}/ccpp")

#------------------------------------------------------------------------------
# Run ccpp_validator
#------------------------------------------------------------------------------
ccpp_validator(VERBOSITY ${CCPP_VERBOSITY}
               SOURCE_FILES ${SCHEME_FORTRAN_FILES}
               METADATA_FILES ${SCHEME_METADATA_FILES}
               TYPE "SCHEME")
ccpp_validator(VERBOSITY ${CCPP_VERBOSITY}
               SOURCE_FILES ${HOST_FORTRAN_FILES}
               METADATA_FILES ${HOST_METADATA_FILES}
               TYPE "HOST")

if(CCPP_32BIT)
  message(STATUS "Compile CCPP slow physics with 32-bit precision")
  set(CCPP_KINDS "kind_phys=REAL32")
endif()

#------------------------------------------------------------------------------
# Run ccpp_capgen_ng
#------------------------------------------------------------------------------
ccpp_capgen(VERBOSITY ${CCPP_VERBOSITY}
            HOSTFILES ${HOST_METADATA_FILES}
            SCHEMEFILES ${SCHEME_METADATA_FILES}
            SUITES ${SUITE_FILES}
            HOST_NAME ${HOST}
            OUTPUT_ROOT "${OUTPUT_ROOT}"
            KIND_SPECS "${CCPP_KINDS}")

# Retrieve the list of Fortran files required for test host from datatable.xml;
# this includes capgen-generated files and dependencies inferred from metadata
ccpp_datafile(DATATABLE "${OUTPUT_ROOT}/datatable.xml"
              REPORT_NAME "--dependencies")
set(CAPGEN_DEPENDENCIES ${CCPP_FILES})
ccpp_datafile(DATATABLE "${OUTPUT_ROOT}/datatable.xml"
              REPORT_NAME "--scheme-files")
set(SCHEME_FORTRAN_FILES ${CCPP_FILES})
ccpp_datafile(DATATABLE "${OUTPUT_ROOT}/datatable.xml"
              REPORT_NAME "--capgen-files")
set(CAPGEN_FILES ${CCPP_FILES})

message(STATUS "List of capgen dependencies: ${CAPGEN_DEPENDENCIES}")
message(STATUS "List of filtered scheme files: ${SCHEME_FORTRAN_FILES}")
message(STATUS "List of capgen-generated files: ${CAPGEN_FILES}")

#------------------------------------------------------------------------------
# Set compile options for either 32/64bit physics build
#------------------------------------------------------------------------------
if(CCPP_32BIT)
  message(STATUS "Compile CCPP slow physics with 32-bit precision")
  add_definitions(-DSINGLE_PREC)
  add_definitions(-DRTE_USE_SP)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set_source_files_properties(
      ${CAPGEN_DEPENDENCIES} ${SCHEME_FORTRAN_FILES} ${HOST_FORTRAN_FILES} ${CAPGEN_FILES}
      PROPERTIES
      COMPILE_OPTIONS "-real-size;32"
    )
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set_source_files_properties(
      ${CAPGEN_DEPENDENCIES} ${SCHEME_FORTRAN_FILES} ${HOST_FORTRAN_FILES} ${CAPGEN_FILES}
      PROPERTIES
      COMPILE_OPTIONS "-fno-default-real-8;-fdefault-double-8"
    )
  endif()
else(CCPP_32BIT)
  message(STATUS "Compile CCPP slow physics with 64-bit precision")
  remove_definitions(-DSINGLE_PREC)
  remove_definitions(-DRTE_USE_SP)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set_source_files_properties(
      ${CAPGEN_DEPENDENCIES} ${SCHEME_FORTRAN_FILES} ${HOST_FORTRAN_FILES} ${CAPGEN_FILES}
      PROPERTIES
      COMPILE_OPTIONS "-real-size;64"
    )
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set_source_files_properties(
      ${CAPGEN_DEPENDENCIES} ${SCHEME_FORTRAN_FILES} ${HOST_FORTRAN_FILES} ${CAPGEN_FILES}
      PROPERTIES
      COMPILE_OPTIONS "-fdefault-real-8;-fdefault-double-8"
    )
  endif()
endif(CCPP_32BIT)

#------------------------------------------------------------------------------
# Set compiler options for 32bit dynamics build	(only fv_sat_adj.F90)
#------------------------------------------------------------------------------
set(SCHEMES_DYNAMICS ${SCHEME_FILES_DIR}/MP/GFDL/fv_sat_adj.F90)

if(32BIT)
  message(STATUS "Compile CCPP fast physics with 32-bit precision")
  add_definitions(-DOVERLOAD_R4)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set_property(SOURCE ${SCHEMES_DYNAMICS} APPEND PROPERTY COMPILE_OPTIONS -real-size 32)
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set_property(SOURCE ${SCHEMES_DYNAMICS} APPEND PROPERTY COMPILE_OPTIONS -fno-default-real-8 -fdefault-double-8)
  endif()
else()
  message(STATUS "Compile CCPP fast physics with 64-bit precision")
  remove_definitions(-DOVERLOAD_R8)
  remove_definitions(-DOVERLOAD_R4)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set_property(SOURCE ${SCHEMES_DYNAMICS} APPEND PROPERTY COMPILE_OPTIONS -real-size 64)
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set_property(SOURCE ${SCHEMES_DYNAMICS} APPEND PROPERTY COMPILE_OPTIONS -fdefault-real-8 -fdefault-double-8)
  endif()
endif()

#------------------------------------------------------------------------------
# Disable openmp for rte-rrtmgp source files
#------------------------------------------------------------------------------
message(STATUS "Disabling OpenMP for rte-rrtmgp radiation")
set(RRTMGP_SOURCES
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/gas-optics/mo_gas_optics.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/gas-optics/mo_gas_optics_constants.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/gas-optics/mo_gas_concentrations.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/gas-optics/mo_gas_optics_util_string.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rrtmgp-frontend/mo_gas_optics_rrtmgp.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rrtmgp-kernels/mo_gas_optics_rrtmgp_kernels.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rrtmgp-frontend/mo_cloud_optics_rrtmgp.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/extensions/mo_rrtmgp_clr_all_sky.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/extensions/mo_fluxes_byband.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/extensions/solar_variability/mo_solar_variability.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/extensions/mo_heating_rates.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/extensions/mo_fluxes_bygpoint.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/extensions/mo_compute_bc.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/extensions/mo_cloud_sampling.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rte-frontend/mo_rte_util_array_validation.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rte-frontend/mo_rte_config.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rte-frontend/mo_source_functions.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rte-frontend/mo_rte_sw.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rte-frontend/mo_fluxes.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rte-frontend/mo_rte_lw.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rte-kernels/mo_rte_util_array.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rte-kernels/mo_rte_solver_kernels.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rte-kernels/mo_optical_props_kernels.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rte-kernels/mo_fluxes_broadband_kernels.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rte-frontend/mo_rte_kind.F90
    ${SCHEME_FILES_DIR}/Radiation/RRTMGP/rte-rrtmgp/rte-frontend/mo_optical_props.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set_property(SOURCE ${RRTMGP_SOURCES} APPEND PROPERTY COMPILE_OPTIONS -fno-openmp)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set_property(SOURCE ${RRTMGP_SOURCES} APPEND PROPERTY COMPILE_OPTIONS -qno-openmp)
else()
  message(FATAL_ERROR "Flags to disable OpenMP not defined for Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")
endif()

#------------------------------------------------------------------------------
# Set flags for 32-bit rte-rrtmgp radiation
#------------------------------------------------------------------------------
if(RRTMGP_32BIT)
  message(STATUS "Compile CCPP RTE-RRTMGP with 32-bit precision")
  add_definitions(-DRTE_USE_SP)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set_property(SOURCE ${RRTMGP_SOURCES} APPEND PROPERTY COMPILE_OPTIONS -fno-default-real-8 -fdefault-double-8)
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set_property(SOURCE ${RRTMGP_SOURCES} APPEND PROPERTY COMPILE_OPTIONS -real-size 32)
  else()
    message(FATAL_ERROR "rte-rrtmgp 32-bit not supported for Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")
  endif()
endif()

#------------------------------------------------------------------------------
# Set flags for SAMF Convection (32-bit option)
#------------------------------------------------------------------------------
option(SAMFCNV_32BIT "Compile SAMF Convection with 32-bit precision" OFF)

if(SAMFCNV_32BIT)
  message(STATUS "Compile CCPP SAMFCNV with 32-bit precision")

  set(SAMF_NAMES
    mo_conv_kind.F90
    samfdeepcnv.f
    samfshalcnv.f
    progomega_calc.f90
    progsigma_calc.f90
  )
  set(SAMF_FORTRAN_FILES)

  foreach(NAME IN LISTS SAMF_NAMES)
    set(TMP_LIST ${SCHEME_FORTRAN_FILES})
    list(FILTER TMP_LIST INCLUDE REGEX ".*${NAME}$")
    list(APPEND SAMF_FORTRAN_FILES ${TMP_LIST})
  endforeach()

  if(SAMF_FORTRAN_FILES)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
      set_property(SOURCE ${SAMF_FORTRAN_FILES} APPEND PROPERTY COMPILE_OPTIONS "-real-size;32;-DSAMFCNV_USE_SP")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      set_property(SOURCE ${SAMF_FORTRAN_FILES} APPEND PROPERTY COMPILE_OPTIONS "-fno-default-real-8;-fdefault-double-8;-DSAMFCNV_USE_SP")
    endif()
  endif()
endif()

#------------------------------------------------------------------------------
# For CAPS, remove bounds checks with Intel Classic (ifort), Intel LLVM (ifx),
# and GNU (gfortran) to avoid the compilation being killed or taking forever
#------------------------------------------------------------------------------
message(STATUS "Removing bounds checks for autogenerated capgen files")
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM")
  set_property(SOURCE ${CAPGEN_FILES} APPEND PROPERTY COMPILE_OPTIONS -check nobounds)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set_property(SOURCE ${CAPGEN_FILES} APPEND PROPERTY COMPILE_OPTIONS -fcheck=no-bounds)
endif()

#------------------------------------------------------------------------------
# Cap the optimization level of the (large) autogenerated files  at -O2: at -O3
# they take a very long time to compile for negligible runtime benefit.
# Appended after the inherited Release flags, so the trailing -O2 wins.
#------------------------------------------------------------------------------
if (CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "Reducing optimization to -O2 for autogenerated capgen files")
  set_property(SOURCE ${CAPGEN_FILES} APPEND PROPERTY COMPILE_OPTIONS -O0)
endif()

#------------------------------------------------------------------------------
# GFS_typedefs.F90 is a very large host file that cannot be optimized with
# Intel in Release mode without enormous compile times; build it at -O0.
#------------------------------------------------------------------------------
if (CMAKE_BUILD_TYPE STREQUAL "Release" AND (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM"))
  message(STATUS "Disabling optimization (-O0) for GFS_typedefs.F90 with Intel")
  set_property(SOURCE ${HOST_FILES_DIR}/data/GFS_typedefs.F90 APPEND PROPERTY COMPILE_OPTIONS -O0)
endif()

#------------------------------------------------------------------------------
# Schemes that Intel mis-optimizes at -O2/-O3 in Release mode (runtime SIGSEGV);
# build them at -O1. GNU is unaffected, so this is Intel/IntelLLVM only.
#------------------------------------------------------------------------------
if (CMAKE_BUILD_TYPE STREQUAL "Release" AND (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM"))
  message(STATUS "Lowering optimization (-O1) for select schemes with Intel")
  set_property(SOURCE
    ${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/gcycle.F90
    ${SCHEME_FILES_DIR}/MP/NSSL/module_mp_nssl_2mom.F90
    ${SCHEME_FILES_DIR}/PBL/MYNN_EDMF/mynnedmf_wrapper.F90
    APPEND PROPERTY COMPILE_OPTIONS -O1)
endif()

#------------------------------------------------------------------------------
# Build CCPP_TARGET
#------------------------------------------------------------------------------
set(BUILD_SHARED_LIBS OFF)
add_library(${CCPP_TARGET}
  ${CAPGEN_DEPENDENCIES}
  ${SCHEME_FORTRAN_FILES}
  ${HOST_FORTRAN_FILES}
  ${CAPGEN_FILES}
)

# UFSFV3 requires FMS
if(FV3)
  target_link_libraries(${CCPP_TARGET} fms)
endif()

target_link_libraries(${CCPP_TARGET} MPI::MPI_Fortran)
if(OPENMP)
  target_link_libraries(${CCPP_TARGET} OpenMP::OpenMP_Fortran)
endif()
target_link_libraries(${CCPP_TARGET} NetCDF::NetCDF_Fortran)
if(TARGET bacio::bacio)
  set(BACIO_TARGET bacio::bacio)
elseif(TARGET bacio::bacio_4)
  set(BACIO_TARGET bacio::bacio_4)
else()
  message(FATAL_ERROR
    "bacio found but expected target not exported (bacio::bacio or bacio::bacio_4)")
endif()
target_link_libraries(${CCPP_TARGET} ${BACIO_TARGET})
if(ip_FOUND)
  target_link_libraries(${CCPP_TARGET} ip::ip_d)
else()
  target_link_libraries(${CCPP_TARGET} sp::sp_d)
endif()
target_link_libraries(${CCPP_TARGET} w3emc::w3emc_d)

set_target_properties(${CCPP_TARGET} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
target_include_directories(${CCPP_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>)

#------------------------------------------------------------------------------
# Install
#------------------------------------------------------------------------------
install(
  TARGETS ${CCPP_TARGET}
  EXPORT ${CCPP_TARGET}-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX})

install(EXPORT ${CCPP_TARGET}-config
  DESTINATION lib/cmake)
