#!/usr/bin/env bash

FC=gfortran-8
FFLAGS="-O0 -g -fimplicit-none  -Wall  -Wextra  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace"
FFLAGS+="  -Wno-unused-function  -Wno-unused-dummy-argument"

rm -f *.mod *.o a.out
${FC} ${FFLAGS} main.f90

./a.out

