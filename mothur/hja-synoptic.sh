#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=249gb,walltime=240:00:00
#PBS -M wisnoski@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/hja-synoptic
module load gcc/9.1.0
/N/u/wisnoski/Carbonate/mothur/mothur /gpfs/home/w/i/wisnoski/Carbonate/GitHub/hja-synoptic/mothur/hja-synoptic.batch
