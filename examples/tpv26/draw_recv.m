clc
clear
%close all
addmypath

id = 3;
varnm = 'ratem';
%varnm = 'all';

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

%[ t, v, coord1 ] = extract_seismo( data_dir, nproc, coord, bc, varnm );
[ t, v, bc, coord, nor ] = extract_seismo_from_id( data_dir, nproc, id, varnm );

figure
plot(t,v)
title(varnm)
xlabel('Time (sec)')