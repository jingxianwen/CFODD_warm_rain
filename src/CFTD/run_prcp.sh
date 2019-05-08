#!/bin/bash
### remove exacutable file
if [ -f "prcp.out" ]
then
        rm ./prcp.out
fi

### compile and run
gfortran-7 -o prcp.out -fconvert=big-endian prog_miroc_prcp_pdf.f90
#gfortran-7 -o a.out -fconvert=big-endian prog_miroc_cfodd_subcol_adiabatic_aerosol_categ.f90
#gfortran-7 -o a.out -fconvert=big-endian prog_miroc_cfodd_subcol_adiabatic.f90
#gfortran-7 -o a.out -fconvert=big-endian prog_miroc_cfodd_subcol.f90
if [ -f "prcp.out" ]
then
        ./prcp.out
fi

### remove exacutable file
if [ -f "prcp.out" ]
then
        rm ./prcp.out
fi
