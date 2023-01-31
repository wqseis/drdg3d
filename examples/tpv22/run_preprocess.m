
% read mesh 
% and build mesh.nc
conf_mesh;

% modify media (vp,vs,rho) in mesh.nc
conf_media;

% modify initial stress and 
% frictions (Tx0,Ty0,Tz0,mu_s,mu_d,Dc,C0, etc.) in mesh.nc
conf_stress;

% mesh partition using METIS
conf_part;

% build mesh connectivity
%system('../../bin/exe_get_neigh mesh.nc');

%system('mkdir -p data');

% write mesh database in parallel
% results will be written in
% data/meshVar000000.nc
% data/meshVar000001.nc
% ...
% data/meshVar00000N.nc (N=nproc)
% system('../../bin/exe_part_mesh mesh.nc');
