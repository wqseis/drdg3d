clc
clear
%close all
addmypath

varnm = 'rate';
id = 1;

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

[ t, v, bc, coord, nor ] = extract_seismo_from_id( data_dir, nproc, id, varnm);

figure
plot(t,v,'x-')
title([varnm,': ',num2str(coord')])
xlabel('Time (sec)')
axis tight
