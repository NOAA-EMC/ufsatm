
# ufsatm

This repository contains a driver and key subcomponents of the
atmospheric component of the NOAA's [Unified Forecast System
(UFS)](https://ufscommunity.org/) weather model.

The subcomponents include:

- The Finite-Volume Cubed-Sphere (FV3) dynamical core, originally
   from the [Geophysical Fluid Dynamics Laboratory](https://www.gfdl.noaa.gov/).

- The Model for Prediction Across Scales - Atmosphere (MPAS-A)	dynamical
   core developed by [NSF-NCAR Mesoscale and Micrometeorology Laboratory (MMM)](https://www.mmm.ucar.edu).
   - [MPAS Model](https://github.com/ufs-community/MPAS-Model)
   - [MPAS documentation](https://www.mmm.ucar.edu/models/mpas)

- The Common Community Physics Package (CCPP) supported by the
   [Developmental Testbed Center
   (DTC)](https://dtcenter.org/community-code/common-community-physics-package-ccpp),
   including:
   - [CCPP Framework](https://github.com/NCAR/ccpp-framework).
   - [CCPP Physics](https://github.com/NCAR/ccpp-physics)

- wrapper code to call [UFS stochastic
   physics](https://stochastic-physics.readthedocs.io/en/latest/)

- The io code handles netCDF I/O.

- The cpl coupler code connects the different components and allows
   them to communicate.

## Prerequisites

This package requires the following
[NCEPLIBS](https://github.com/NOAA-EMC/NCEPLIBS) packages:
 - [NCEPLIBS-w3emc](https://github.com/NOAA-EMC/NCEPLIBS-w3emc)
 - [NCEPLIBS-bacio](https://github.com/NOAA-EMC/NCEPLIBS-bacio)
 - [NCEPLIBS-nemsio](https://github.com/NOAA-EMC/NCEPLIBS-nemsio)
 - [NCEPLIBS-sp](https://github.com/NOAA-EMC/NCEPLIBS-sp)

If the INLINE_POST cmake variable is set, the upp library will be
needed:
 - [Unified Post Processing Library](https://github.com/NOAA-EMC/EMC_post)

This package also requires the following external packages:
 - [netcdf-c Library](https://github.com/Unidata/netcdf-c)
 - [netcdf-fortran Library](https://github.com/Unidata/netcdf-fortran)
 - [ESMF](https://github.com/esmf-org/esmf)
 - [GFDL's Flexible Modeling System](https://github.com/NOAA-GFDL/FMS)

## Obtaining ufsatm

To obtain ufsatm, clone the git repository, and update the submodules:

```
git clone https://github.com/NOAA-EMC/ufsatm.git
cd ufsatm
git submodule update --init --recursive
```

## Disclaimer

The United States Department of Commerce (DOC) GitHub project code is
provided on an "as is" basis and the user assumes responsibility for
its use. DOC has relinquished control of the information and no longer
has responsibility to protect the integrity, confidentiality, or
availability of the information. Any claims against the Department of
Commerce stemming from the use of its GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial
products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of
Commerce. The Department of Commerce seal and logo, or the seal and
logo of a DOC bureau, shall not be used in any manner to imply
endorsement of any commercial product or activity by DOC or the United
States Government.

