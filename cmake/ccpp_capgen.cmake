# CMake function to collect source files referenced by metadata files
#
# Usage:
#   ccpp_source_files(OUT_VAR ${META_FILES})
#
# Example:
#   ccpp_source_files(CCPP_SRCS ${META_FILES})
#
# Result:
#   OUT_VAR contains all matching source files
#
function(ccpp_source_files OUT_VAR)
  set(result)
  foreach(meta_file IN LISTS ARGN)
    # Convert meta_file to an absolute path first
    get_filename_component(
        meta_file
        "${meta_file}"
        ABSOLUTE
    )
    # Then extract the directory
    get_filename_component(
        meta_dir
        "${meta_file}"
        DIRECTORY
    )
    
    # Read metadata file
    file(READ "${meta_file}" meta_contents)
    
    # Extract: source_path = something
    string(REGEX MATCH
        "source_path[ \t]*=[ \t]*([^\r\n]+)"
        _match
        "${meta_contents}"
    )
    
    # If no source_path exists, use directory of metadata file
    if(NOT _match)
        set(source_dir ${meta_dir})
    else()
      # Extract just the path
      string(REGEX REPLACE
          "source_path[ \t]*=[ \t]*"
          ""
          source_path
          "${_match}"
      )
      string(STRIP "${source_path}" source_path)
      # Absolute directory containing source
      get_filename_component(
        source_dir
        "${meta_dir}/${source_path}"
        ABSOLUTE
      )
    endif()
    
    # Try supported Fortran suffixes
    foreach(ext
      f
      F
      f90
      F90
      F95
      f95
    )
      string(REGEX REPLACE
          "\\.meta$"
          ".${ext}"
          candidate
          "${meta_file}"
      )
      if(EXISTS "${candidate}")
          list(APPEND result "${candidate}")
          break()
      endif()
    endforeach()
  endforeach()
  
  # Remove duplicates and return
  list(REMOVE_DUPLICATES result)
  set(${OUT_VAR}
      ${result}
      PARENT_SCOPE)
endfunction()


# CMake wrapper for ccpp_validator.py
#
# SOURCE_FILES   - CMake list of Fortran source files
# METADATA_FILES - CMake list of corresponding metadata files
# TYPE           - "HOST" or "SCHEME" (case-insensitive)
function(ccpp_validator)
  set(optionalArgs)
  set(oneValueArgs VERBOSITY TYPE)
  set(multi_value_keywords SOURCE_FILES METADATA_FILES)
  cmake_parse_arguments(arg "${optionalArgs}" "${oneValueArgs}" "${multi_value_keywords}" ${ARGN})

  # Error if script file not found.
  set(CCPP_VALIDATOR_CMD_LIST "${CMAKE_SOURCE_DIR}/${PARENT_DIR}/ccpp/framework/capgen/ccpp_validator.py")
  if(NOT EXISTS ${CCPP_VALIDATOR_CMD_LIST})
    message(FATAL_ERROR "function(ccpp_validator): Could not find ccpp_validator.py.  Looked for ${CCPP_VALIDATOR_CMD_LIST}.")
  endif()

  # Interpret parsed arguments
  if(NOT DEFINED arg_SOURCE_FILES)
    message(FATAL_ERROR "function(ccpp_validator): SOURCE_FILES not set.")
  endif()
  list(JOIN arg_SOURCE_FILES "," SOURCE_FILES_SEPARATED)
  list(APPEND CCPP_VALIDATOR_CMD_LIST "--source-files" "${SOURCE_FILES_SEPARATED}")

  if(NOT DEFINED arg_METADATA_FILES)
    message(FATAL_ERROR "function(ccpp_validator): METADATA_FILES not set.")
  endif()
  list(JOIN arg_METADATA_FILES "," METADATA_FILES_SEPARATED)

  if(NOT DEFINED arg_TYPE)
    message(FATAL_ERROR "function(ccpp_validator): TYPE must be HOST or SCHEME")
  endif()
  string(TOUPPER "${arg_TYPE}" _type)
  if(NOT (_type MATCHES "^(HOST|SCHEME)$"))
    message(FATAL_ERROR "function(ccpp_validator): TYPE must be HOST or SCHEME")
  endif()

  if(_type MATCHES "^HOST$")
    list(APPEND CCPP_VALIDATOR_CMD_LIST "--host-files" "${METADATA_FILES_SEPARATED}")
  endif()
  if(_type MATCHES "^SCHEME$")
    list(APPEND CCPP_VALIDATOR_CMD_LIST "--scheme-files" "${METADATA_FILES_SEPARATED}")
  endif()

  if(DEFINED arg_VERBOSITY)
    string(REPEAT "--verbose " ${arg_VERBOSITY} VERBOSE_PARAMS_SEPARATED)
    separate_arguments(VERBOSE_PARAMS UNIX_COMMAND "${VERBOSE_PARAMS_SEPARATED}")
    list(APPEND CCPP_VALIDATOR_CMD_LIST ${VERBOSE_PARAMS})
  endif()

  # DH* 20260513 TEMPORARY: add --legacy-mode to allow parsing
  # of horizontal_loop_extent for scheme run phase metadata
  list(APPEND CCPP_VALIDATOR_CMD_LIST "--legacy-mode")
  # *DH

  message(STATUS "Running ccpp_validator.py from ${CMAKE_CURRENT_SOURCE_DIR}")

  unset(VALIDATOR_OUT)
  execute_process(COMMAND ${CCPP_VALIDATOR_CMD_LIST}
                  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                  OUTPUT_VARIABLE VALIDATOR_OUT
                  ERROR_VARIABLE VALIDATOR_OUT
                  RESULT_VARIABLE RES
                  COMMAND_ECHO STDOUT)
  if(RES EQUAL 0)
    message(STATUS "ccpp-validator stdout: ${VALIDATOR_OUT}")
  else()
    message(STATUS "${CCPP_VALIDATOR_CMD_LIST} FAILED: result = ${RES}")
    message(STATUS "ccpp-validator stdout: ${VALIDATOR_OUT}")
    message(FATAL_ERROR "Validation of source files failed, abort.")
  endif()

endfunction()


# CMake wrapper for ccpp_capgen.py
#
# TRACE          - ON/OFF (Default: OFF) - Add --trace flag to capgen call
# HOST_NAME      - String name of host (drives <host>_ccpp_cap.F90 filename
#                  and module name; required)
# OUTPUT_ROOT    - String path to put generated caps
# VERBOSITY      - Number of --verbose flags to pass to capgen
# HOSTFILES      - CMake list of host metadata filenames
# SCHEMEFILES    - CMake list of scheme metadata files
# SUITES         - CMake list of suite xml files
# KIND_SPECS     - Comma-separated kind mappings, e.g. "kind_phys=REAL32" or
#                  "kind_phys=my_mod:kind_r4,kind_dyn=REAL64".  Each pair is
#                  forwarded as `--kind-type <pair>` to capgen (see the
#                  capgen docstring for the `<name>=[<module>:]<spec>`
#                  grammar; bare ISO specs default to iso_fortran_env).
function(ccpp_capgen)
  set(optionalArgs TRACE)
  set(oneValueArgs HOST_NAME OUTPUT_ROOT VERBOSITY KIND_SPECS)
  set(multi_value_keywords HOSTFILES SCHEMEFILES SUITES)

  cmake_parse_arguments(arg "${optionalArgs}" "${oneValueArgs}" "${multi_value_keywords}" ${ARGN})

  # Error if script file not found.
  set(CCPP_CAPGEN_CMD_LIST "${CMAKE_SOURCE_DIR}/${PARENT_DIR}/ccpp/framework/capgen/ccpp_capgen.py")
  if(NOT EXISTS ${CCPP_CAPGEN_CMD_LIST})
    message(FATAL_ERROR "function(ccpp_capgen): Could not find ccpp_capgen.py.  Looked for ${CCPP_CAPGEN_CMD_LIST}.")
  endif()

  # Interpret parsed arguments
  if(NOT DEFINED arg_HOSTFILES)
    message(FATAL_ERROR "function(ccpp_capgen): HOSTFILES not set.")
  endif()
  list(JOIN arg_HOSTFILES "," HOSTFILES_SEPARATED)
  list(APPEND CCPP_CAPGEN_CMD_LIST "--host-files" "${HOSTFILES_SEPARATED}")

  if(NOT DEFINED arg_SCHEMEFILES)
    message(FATAL_ERROR "function(ccpp_capgen): SCHEMEFILES not set.")
  endif()
  list(JOIN arg_SCHEMEFILES "," SCHEMEFILES_SEPARATED)
  list(APPEND CCPP_CAPGEN_CMD_LIST "--scheme-files" "${SCHEMEFILES_SEPARATED}")

  if(NOT DEFINED arg_SUITES)
    message(FATAL_ERROR "function(ccpp_capgen): SUITES not set.")
  endif()
  list(JOIN arg_SUITES "," SUITES_SEPARATED)
  list(APPEND CCPP_CAPGEN_CMD_LIST "--suites" "${SUITES_SEPARATED}")

  if(NOT DEFINED arg_HOST_NAME)
    message(FATAL_ERROR "function(ccpp_capgen): HOSTNAME not set.")
  endif()
  list(APPEND CCPP_CAPGEN_CMD_LIST "--host-name" "${arg_HOST_NAME}")

  if(NOT DEFINED arg_OUTPUT_ROOT)
    message(FATAL_ERROR "function(ccpp_capgen): OUTPUT_ROOT not set.")
  endif()
  #file(MAKE_DIRECTORY "${arg_OUTPUT_ROOT}")
  list(APPEND CCPP_CAPGEN_CMD_LIST "--output-root" "${arg_OUTPUT_ROOT}")

  if(DEFINED arg_VERBOSITY)
    string(REPEAT "--verbose " ${arg_VERBOSITY} VERBOSE_PARAMS_SEPARATED)
    separate_arguments(VERBOSE_PARAMS UNIX_COMMAND "${VERBOSE_PARAMS_SEPARATED}")
    list(APPEND CCPP_CAPGEN_CMD_LIST ${VERBOSE_PARAMS})
  endif()

  if(DEFINED arg_KIND_SPECS)
    # Accept either a comma-separated string ("kind_phys=REAL64,kind_dyn=REAL32")
    # or a CMake list of pairs.  Each pair becomes a separate
    # `--kind-type <pair>` argv pair so capgen's argparse sees one
    # `--kind-type` per pair (the flag is `action='append'`).
    string(REPLACE "," ";" KIND_SPEC_LIST "${arg_KIND_SPECS}")
    foreach(pair IN LISTS KIND_SPEC_LIST)
      list(APPEND CCPP_CAPGEN_CMD_LIST "--kind-type" "${pair}")
    endforeach()
  endif()

  if(arg_TRACE)
    list(APPEND CCPP_CAPGEN_CMD_LIST "--trace")
  endif()

  # DH* 20260513 TEMPORARY: add --legacy-mode to allow parsing
  # of horizontal_loop_extent for scheme run phase metadata
  list(APPEND CCPP_CAPGEN_CMD_LIST "--legacy-mode")
  # *DH

  # DH* 20260521 TEMPORARY: add --gfs-dim-aliases so that
  # capgen considers specific dimension names as equal
  list(APPEND CCPP_CAPGEN_CMD_LIST "--gfs-dim-aliases")
  # *DH

  # DH* 20260514 TEMPORARY: add --no-host-introspection to
  # suppress writing lists of variables etc to ccpp_static_api
  list(APPEND CCPP_CAPGEN_CMD_LIST "--no-host-introspection")
  # *DH

  message(STATUS "Running ccpp_capgen.py from ${CMAKE_CURRENT_SOURCE_DIR}")

  # Unset CAPGEN_OUT to prevent incorrect output on subsequent ccpp_capgen(...) calls
  unset(CAPGEN_OUT)
  execute_process(COMMAND ${CCPP_CAPGEN_CMD_LIST}
                  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                  OUTPUT_VARIABLE CAPGEN_OUT
                  ERROR_VARIABLE CAPGEN_OUT
                  RESULT_VARIABLE RES
                  COMMAND_ECHO STDOUT)

  message(STATUS "ccpp-capgen stdout: ${CAPGEN_OUT}")

endfunction()


# CMake wrapper for ccpp_datafile.py
#
# DATATABLE   - Path to generated datatable.xml file
# REPORT_NAME - String report name to get list of generated files form capgen (typically --ccpp-files)
function(ccpp_datafile)
  set(mandatoryArgs DATATABLE REPORT_NAME)
  cmake_parse_arguments(arg "" "${mandatoryArgs}" "" ${ARGN})

  set(CCPP_DATAFILE_CMD "${CMAKE_SOURCE_DIR}/${PARENT_DIR}/ccpp/framework/capgen/ccpp_datafile.py")

  if(NOT EXISTS ${CCPP_DATAFILE_CMD})
    message(FATAL_ERROR "function(ccpp_datafile): Could not find ccpp_datafile.py.  Looked for ${CCPP_DATAFILE_CMD}.")
  endif()

  if(NOT DEFINED arg_REPORT_NAME)
    message(FATAL_ERROR "function(ccpp_datafile): REPORT_NAME not set.")
  endif()
  list(APPEND CCPP_DATAFILE_CMD "${arg_REPORT_NAME}")

  if(NOT DEFINED arg_DATATABLE)
    message(FATAL_ERROR "function(ccpp_datafile): DATATABLE not set.")
  endif()
  list(APPEND CCPP_DATAFILE_CMD "${arg_DATATABLE}")

  message(STATUS "Running ccpp_datafile from ${CMAKE_CURRENT_SOURCE_DIR}")

  # Unset CCPP_FILES to prevent incorrect output on subsequent ccpp_datafile(...) calls
  unset(CCPP_FILES)
  execute_process(COMMAND ${CCPP_DATAFILE_CMD}
                  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                  OUTPUT_VARIABLE CCPP_FILES
                  RESULT_VARIABLE RES
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE
                  COMMAND_ECHO STDOUT)
  #message(STATUS "CCPP_FILES = ${CCPP_FILES}")
  if(RES EQUAL 0)
    message(STATUS "CCPP files retrieved")
  else()
    message(FATAL_ERROR "CCPP file retrieval FAILED: result = ${RES}")
  endif()
  if(CCPP_FILES)
    # Convert "," separated list from python back to ";" separated list for CMake
    string(REPLACE "," ";" CCPP_FILES ${CCPP_FILES})
  endif()
  set(CCPP_FILES "${CCPP_FILES}" PARENT_SCOPE)
endfunction()
