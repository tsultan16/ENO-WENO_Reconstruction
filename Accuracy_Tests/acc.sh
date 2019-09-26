#!/bin/sh -l
gfortran accuracy.f90 -o accuracy
./accuracy
gnuplot "plotErr.p"
