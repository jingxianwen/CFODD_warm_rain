#!/bin/bash
if [ -f "b.out" ]
then
        rm ./b.out
fi
gfortran-7 -o b.out -fconvert=big-endian prog_miroc_cftd_subcol.f90
#gfortran-7 -o a.out -fconvert=big-endian prog_miroc_cfodd_subcol_adiabatic_aerosol_categ.f90
#gfortran-7 -o a.out -fconvert=big-endian prog_miroc_cfodd_subcol_adiabatic.f90
#gfortran-7 -o a.out -fconvert=big-endian prog_miroc_cfodd_subcol.f90
if [ -f "b.out" ]
then
        ./b.out
fi