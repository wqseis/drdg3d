% clc
% clear
% close all

addmypath;

par = ReadYaml('parameters.yaml');

nproc = par.nproc;
data_dir = par.data_dir;

%data_dir = '/home/wqzhang/scratch/data_tpv24_200m_O4_err'
%data_dir = 'data'

for branch = 0:1

x = gather_fault_var( data_dir, nproc, 'x' );
y = gather_fault_var( data_dir, nproc, 'y' );
z = gather_fault_var( data_dir, nproc, 'z' );
nx = gather_fault_var( data_dir, nproc, 'nx' );
t = gather_fault_var( data_dir, nproc, 'ruptime' );
fault_id = gather_fault_var1d( data_dir, nproc, 'fault_id' );

branch

if branch
idx = find(mean(nx(:,:))>0 & (fault_id ==1) );
else
idx = find(mean(nx(:,:))>0 & (fault_id ==3 | fault_id == 2) );
end 
x = x(:,idx);
y = y(:,idx);
z = z(:,idx);
t = t(:,idx);


if 0
tri = get_face_connect(x);
figure
tricontour(tri,y(:),z(:),t(:),0:0.5:100,'r')
end

%%
if branch
fid = fopen('cplot_branch.txt','wt');
y=y/cosd(30);
else
fid = fopen('cplot_main.txt','wt');
end
fprintf(fid,'# problem = TPV24\n');
fprintf(fid,'# author = Wenqiang Zhang\n');
fprintf(fid,'# date = 2023/02/18\n');
fprintf(fid,'# code = DRDG3D\n');
fprintf(fid,'# code_version = 0.0\n');
fprintf(fid,'# element_size = 200 m on fault, Order 6\n');
fprintf(fid,'# Column #1 = horizontal coordinate, distance along strike (m)\n');
fprintf(fid,'# Column #2 = vertical coordinate, distance down-dip (m)\n');
fprintf(fid,'# Column #3 = rupture time (s)\n');
fprintf(fid,'# The line below lists the names of the data fields:\n');
fprintf(fid,'j        k       t\n');

t(t<0)=1e9;
dat = [y(:)*1e3,z(:)*-1e3,t(:)];
fprintf(fid,'%14.6e %14.6e %14.6e\n',dat');
fclose(fid);

end
