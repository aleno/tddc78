#!/bin/bash
#SBATCH -N 2
#SBATCH -t 00:10:00
mpprun ./blurmain 50 im1.ppm om1.ppm

