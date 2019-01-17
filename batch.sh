#! /bin/bash

sed -e '646,647d;649,650d;653,654d;656,657d;666,890d;897,905d;958,1041d;12d'  mpm/microphysics.f90 > mpm/microphysics.f90.new

mv mpm/microphysics.f90 mpm/microphysics.f90.bak
mv mpm/microphysics.f90.new mpm/microphysics.f90

make
#make NETCDFLIB=-L /lib64/ NETCDFMOD=/usr/lib64/gfortran/modules/ FFLAGSOMP='-O3 -fbounds-check -g -o '

mv mpm/microphysics.f90.bak mpm/microphysics.f90


