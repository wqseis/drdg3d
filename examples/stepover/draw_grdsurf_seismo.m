clc
clear
close all

addmypath;

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

varnm = 'Vy';
 

x0 = 15;
y0 = 5;
z0 = 0; % on the ground surface

[v,t] = extract_seismo_from_snap_grdsurf(data_dir,nproc,[x0,y0,z0],varnm);

nt = length(v);
 
h=figure;
%set(h,'Visible','off');
plot(t,v,'linewidth',1)
xlabel('T (sec)')
%saveas(h,'seismo.png')
