% clc
% clear
% close all

addmypath;

par = ReadYaml('parameters.yaml');

nproc = par.nproc;
data_dir = par.data_dir;

x = gather_fault_var( data_dir, nproc, 'x' );
y = gather_fault_var( data_dir, nproc, 'y' );
z = gather_fault_var( data_dir, nproc, 'z' );
nx = gather_fault_var( data_dir, nproc, 'nx' );
t = gather_fault_var( data_dir, nproc, 'ruptime' );
fault_id = gather_fault_var1d( data_dir, nproc, 'fault_id' );

idx = find(mean(nx(:,:))>0 & fault_id == 1);
x = x(:,idx);
y = y(:,idx);
z = z(:,idx);
t = t(:,idx);

tri = get_face_connect(x);
figure
tricontour(tri,y(:),z(:),t(:),0:0.5:100,'r')
axis equal

%%
fid = fopen('cplot_tpv102.txt','wt');
fprintf(fid,'# problem = TPV102\n');
fprintf(fid,'# author = Wenqiang Zhang\n');
fprintf(fid,'# date = 2023/02/13\n');
fprintf(fid,'# code = DRDG3D\n');
fprintf(fid,'# code_version = 0.0\n');
fprintf(fid,'# element_size = 200 m on fault, Order 4\n');
fprintf(fid,'# Column #1 = horizontal coordinate, distance along strike (m)\n');
fprintf(fid,'# Column #2 = vertical coordinate, distance down-dip (m)\n');
fprintf(fid,'# Column #3 = rupture time (s)\n');
fprintf(fid,'# The line below lists the names of the data fields:\n');
fprintf(fid,'j        k       t\n');

t(t<0)=1e9;
dat = [y(:)*1e3,z(:)*-1e3,t(:)];
fprintf(fid,'%14.6e %14.6e %14.6e\n',dat');
fclose(fid);
