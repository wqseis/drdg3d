

m_scripts = '../../matlab';

% set the directory of metis (version 5)
metis_dir = '/usr/local/bin';
%metis_dir = '/home/wzhang/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-9.4.0/metis-5.1.0-g3hkcjmmfvljvil6u6okitnlxozpw3mn/bin';

addpath(genpath(m_scripts));

myconstants;

MAX_NUM_RECV_FAULT = 100;
MAX_NUM_RECV_FREE = 100;
