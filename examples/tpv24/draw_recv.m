clc
clear
%close all
addmypath

bc = BC_FAULT;
varnm = 'ratem';

strike = 2; dip =   10  ;
ang = 30;
coord = [strike*sind(ang) strike*cosd(ang) -dip];
%coord = [0 -2 -10];

% bc = BC_FREE;
% 
% coord = [0.6 2 0];
% coord = [2.3 8 0];
% coord = [3 -2 0];
% coord = [4.2  2 0];
% coord = [7.6 8 0];
% coord = [-3 -2 0];
% coord = [-3 2 0];
% coord = [-3 8 0];


par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

[ t, v, coord1 ] = extract_seismo_all( data_dir, nproc, coord, bc);

figure
if bc == BC_FAULT
subplot(322)
plot(t,v(1,:));title('Vm')
subplot(324)
plot(t,v(2,:));title('Vl')

subplot(321)
plot(t,v(3,:));title('Tm')
subplot(323)
plot(t,v(4,:));title('Tl')
subplot(325)
plot(t,v(5,:));title('Tn')
xlabel('Time (sec)')
end

if bc == BC_FREE
subplot(312)
plot(t,-v(1,:));title('Vx')
subplot(311)
plot(t,v(2,:));title('Vy')
subplot(313)
plot(t,-v(3,:));title('Vz')
xlabel('Time (sec)')
end