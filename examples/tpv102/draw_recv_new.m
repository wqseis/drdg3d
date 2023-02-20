clc
clear
%close all
addmypath

varnm = 'Vy';
id = 21;
varnm = 'state';
id = 2;

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

[ t, v, bc, coord, nor ] = extract_seismo_from_id( data_dir, nproc, id, varnm, [-1 0 0]);
bc

if strcmp(varnm,'state')
    f0 = 0.6;
    b = 0.012;
    Dc = 0.02;
    V0 = 1e-6;
    %v = Dc/V0*exp((v-f0)/b);
    %v = log10(v);
    v = log10(Dc/V0)+log10(exp(1))*(v-f0)/b;
end

figure
plot(t,v,'-')
title([varnm,': ',num2str(coord')])
xlabel('Time (sec)')
axis tight
