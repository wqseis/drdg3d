#!/bin/bash

set -e
set -x

bin_dir='../../bin'

${bin_dir}/exe_get_neigh mesh.nc
mkdir -p data
${bin_dir}/exe_part_mesh mesh.nc
mpirun --oversubscribe -np 4 ${bin_dir}/exe_solver 
