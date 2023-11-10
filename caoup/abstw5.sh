#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -M mmi728227@gmail.com
#$ -m ea
#$ -V
#                      
#$ -q all.q@mikan
./ana5.out
