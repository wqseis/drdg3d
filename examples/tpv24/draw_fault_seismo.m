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
%varnm = 'tau';
varnm = 'sigma';

x0 = 0; % on the fault surface
y0 = 0;
z0 = -7.5;
%z0 = 0;
ang = 0;
strike = -8; dip = 10;
x0 = strike*sind(ang);
y0 = strike*cosd(ang);
z0 = -dip;
%coord = [0 -2 -10];

%y0 = 1.9;
%z0 = 1.9;
% y0 = 5.6+0.0;
% z0 = 2.3-0.0;

[v,t] = extract_seismo_from_snap(data_dir,nproc,[x0,y0,z0],varnm);

nt = length(v);

%b22 =  0.926793;
%b33 =  1.073206;
%b23 = -0.169029;
%
%Pf = 9.8*dip;
%Sn = -2.67*9.8*dip;
%Sn0=b33*(Sn+Pf)
%
%if strcmp(varnm,'sigma')
%v = v-v(1)+Sn0;
%end
 
h=figure;
%set(h,'Visible','off');
plot(t,v,'linewidth',1)
% hold on
% plot(t,v1,'linewidth',1)
xlabel('T (sec)')
%saveas(h,'seismo.png')
