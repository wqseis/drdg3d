clc
clear
close all

addmypath;

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

data_dir = 'data';

varnm = 'rate';
%varnm = 'slip';
%varnm = 'stress';
%varnm = 'sigma';

x0 = 0; % on the fault surface
y0 = 5;
z0 = -7.5;
%z0 = 0;

%y0 = 1.9;
%z0 = 1.9;
% y0 = 5.6+0.0;
% z0 = 2.3-0.0;

[v,t] = extract_seismo_from_snap(data_dir,nproc,[x0,y0,z0],varnm);

nt = length(v);
 
h=figure;
%set(h,'Visible','off');
plot(t,v,'linewidth',1)
xlabel('T (sec)')
%saveas(h,'seismo.png')
