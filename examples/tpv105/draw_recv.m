clc
clear
%close all
addmypath

bc = BC_FAULT;
coord = [0,-9,-7.5];
varnm = 'state';

% bc = BC_FREE;
% coord = [-6 -12 0];
% varnm = 'Vz';

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

[ t, v, coord1 ] = extract_seismo( data_dir, nproc, coord, bc, varnm );
%v = -v;
u = cumtrapz(t,v);

figure
plot(t,v);
%xlim([0 15])
title(varnm)
xlabel('Time (sec)')