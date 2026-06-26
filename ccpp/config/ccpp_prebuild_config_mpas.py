#!/usr/bin/env python

# CCPP prebuild config for MPAS - Model for Prediction Across Scales


###############################################################################
# Definitions                                                                 #
###############################################################################

HOST_MODEL_IDENTIFIER = "MPAS"

# Add all files with metadata tables on the host model side and in CCPP,
# relative to basedir = top-level directory of host model. This includes
# kind and type definitions used in CCPP physics. Also add any internal
# dependencies of these files to the list.
VARIABLE_DEFINITION_FILES = [
    # actual variable definition files
    'framework/src/ccpp_types.F90',
    'physics/physics/hooks/machine.F',
    'physics/physics/Radiation/RRTMG/radsw_param.f',
    'physics/physics/Radiation/RRTMG/radlw_param.f',
    'physics/physics/photochem/module_ozphys.F90',
    'physics/physics/MP/TEMPO/tempo_v3/src/module_mp_tempo_cfgs.F90',    
    'physics/physics/photochem/module_h2ophys.F90',
    'physics/physics/SFC_Models/Land/Noahmp/lnd_iau_mod.F90',
    '../ccpp/data/CCPP_typedefs.F90',
    '../ccpp/data/GFS_typedefs.F90',
    '../ccpp/data/MPAS_typedefs.F90',
    '../ccpp/data/CCPP_data.F90'
    ]

TYPEDEFS_NEW_METADATA = {
    'ccpp_types' : {
        'ccpp_t' : 'cdata',
        'MPI_Comm' : '',
        'ccpp_types' : '',
        },
    'machine' : {
        'machine' : '',
        },
    'module_radlw_parameters' : {
        'module_radsw_parameters' : '',
        },
    'module_radlw_parameters' : {
        'module_radlw_parameters' : '',
        },
    'module_ozphys' : {
        'module_ozphys' : '',
        'ty_ozphys'     : '',
        },
    'module_mp_tempo_cfgs' : {
        'module_mp_tempo_cfgs' : '',
        'ty_tempo_cfgs'     : '',
        },
    'module_h2ophys' : {
        'module_h2ophys' : '',
        'ty_h2ophys'     : '',
        },
    'land_iau_mod' : {
        'land_iau_mod' : '',
        'land_iau_external_data_type' : '',
        'land_iau_state_type' : '',
        'land_iau_control_type' : '',
        },
    'CCPP_typedefs' : {
        'GFS_interstitial_type' : 'GFS_Interstitial(cdata%thrd_no)',
        'GFDL_interstitial_type' : 'GFDL_interstitial',
        'CCPP_typedefs' : '',
        },
    'CCPP_data' : {
        'CCPP_data' : '',
        },
    'MPAS_typedefs' : {
        'MPAS_typedefs'         : '',
        },
    'GFS_typedefs' : {
        'GFS_control_type'      : 'GFS_Control',
        'GFS_statein_type'      : 'GFS_Statein',
        'GFS_stateout_type'     : 'GFS_Stateout',
        'GFS_grid_type'         : 'GFS_Grid',
        'GFS_tbd_type'          : 'GFS_Tbd',
        'GFS_cldprop_type'      : 'GFS_Cldprop',
        'GFS_sfcprop_type'      : 'GFS_Sfcprop',
        'GFS_radtend_type'      : 'GFS_Radtend',
        'GFS_coupling_type'     : 'GFS_Coupling',
        'GFS_diag_type'         : 'GFS_Intdiag',
        'GFS_typedefs' : '',
        },
    }

# Add all physics scheme files relative to basedir
SCHEME_FILES = [
    # Relative path to source (from where ccpp_prebuild.py is called) : [ list of physics sets in which scheme may be called ];
    # current restrictions are that each scheme can only belong to one physics set, and all schemes within one group in the
    # suite definition file have to belong to the same physics set
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_phys_time_vary.fv3.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rad_time_vary.mpas.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_time_vary_pre.fv3.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_radiation_surface.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_radiation_post.F90',
    'physics/physics/MP/TEMPO/mp_tempo.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmg_pre.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmg_setup.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_setup.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_pre.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_cloud_diagnostics.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_cloud_mp.F90',
    'physics/physics/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmgp_cloud_overlap.F90',
    'physics/physics/Radiation/RRTMG/radsw_main.F90',
    'physics/physics/Radiation/RRTMG/radlw_main.F90',
    'physics/physics/Radiation/RRTMG/rrtmg_lw_post.F90',
    'physics/physics/Radiation/RRTMG/rrtmg_sw_post.F90',
    'physics/physics/Radiation/RRTMG/rad_sw_pre.F90',
    'physics/physics/Radiation/RRTMGP/rrtmgp_aerosol_optics.F90',
    'physics/physics/Radiation/RRTMGP/rrtmgp_lw_main.F90',
    'physics/physics/Radiation/RRTMGP/rrtmgp_sw_main.F90',
    'physics/physics/photochem/module_h2ophys.F90',
    'physics/physics/photochem/module_ozphys.F90'
]

# Default build dir, relative to current working directory,
# if not specified as command-line argument
DEFAULT_BUILD_DIR = 'build'

# Auto-generated makefile/cmakefile snippets that contain all type definitions
TYPEDEFS_MAKEFILE   = '{build_dir}/physics/CCPP_TYPEDEFS.mk'
TYPEDEFS_CMAKEFILE  = '{build_dir}/physics/CCPP_TYPEDEFS.cmake'
TYPEDEFS_SOURCEFILE = '{build_dir}/physics/CCPP_TYPEDEFS.sh'

# Auto-generated makefile/cmakefile snippets that contain all schemes
SCHEMES_MAKEFILE   = '{build_dir}/physics/CCPP_SCHEMES.mk'
SCHEMES_CMAKEFILE  = '{build_dir}/physics/CCPP_SCHEMES.cmake'
SCHEMES_SOURCEFILE = '{build_dir}/physics/CCPP_SCHEMES.sh'

# Auto-generated makefile/cmakefile snippets that contain all caps
CAPS_MAKEFILE   = '{build_dir}/physics/CCPP_CAPS.mk'
CAPS_CMAKEFILE  = '{build_dir}/physics/CCPP_CAPS.cmake'
CAPS_SOURCEFILE = '{build_dir}/physics/CCPP_CAPS.sh'

# Directory where to put all auto-generated physics caps
CAPS_DIR = '{build_dir}/physics'

# Directory where the suite definition files are stored
SUITES_DIR = '../ccpp/suites'

# Directory where to write static API to
STATIC_API_DIR = '{build_dir}/physics'
STATIC_API_CMAKEFILE = '{build_dir}/physics/CCPP_STATIC_API.cmake'
STATIC_API_SOURCEFILE = '{build_dir}/physics/CCPP_STATIC_API.sh'

# Directory for writing HTML pages generated from metadata files
# used by metadata2html.py for generating scientific documentation
METADATA_HTML_OUTPUT_DIR = '{build_dir}/physics/physics/docs'

# HTML document containing the model-defined CCPP variables
HTML_VARTABLE_FILE = '{build_dir}/physics/CCPP_VARIABLES_MPAS.html'

# LaTeX document containing the provided vs requested CCPP variables
LATEX_VARTABLE_FILE = '{build_dir}/framework/doc/DevelopersGuide/CCPP_VARIABLES_MPAS.tex'
