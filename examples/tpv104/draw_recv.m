clc
clear
%close all
addmypath

coord = [0 0 -7.5];
varnm = 'ratem';
bc = BC_FAULT;


% bc = BC_FREE;
% coord = [6 12 0];
% varnm = 'Vy';

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

[ t, v, coord1 ] = extract_seismo( data_dir, nproc, coord, bc, varnm );

figure
plot(t,v)
title(varnm)
xlabel('Time (sec)')