#!/bin/sh -l
gfortran pol.f90 -o pol
./pol
gnuplot "plotOutput.p"
