#!/bin/sh -l
gfortran eno.f90 -o eno
./eno
gnuplot "plotOutput.p"
