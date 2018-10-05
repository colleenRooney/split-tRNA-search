#! /bin/bash
cd /home/rm/vfp
g++ -o vfp main.cpp fasta.cpp matrix.cpp footprint.cpp
echo "compilation finished... press <Enter>"
read a

