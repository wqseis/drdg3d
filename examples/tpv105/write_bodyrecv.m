clc
clear
%close all
addmypath

body =   9;
strike =  0;
dip = 0;


bc = BC_FREE;
coord = [-body strike dip];

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

[ t, v, coord1 ] = extract_seismo_all( data_dir, nproc, coord, bc);

v(7:end,:) = [];
% note that our coordinate convention differs from the benchmark website
% normal direction
v(1,:) = -v(1,:);
v(4,:) = -v(4,:);
% vertical direction
v(3,:) = -v(3,:);
v(6,:) = -v(6,:);


figure
for i = 1:3
subplot(3,1,i)
plot(t,v(i+0,:))
end

fnm = sprintf('body%03dst%03ddp%03d.txt',body*10,strike*10,dip*10);
fid = fopen(fnm,'wt');
fprintf(fid,'# problem = TPV105-3D\n');
fprintf(fid,'# author = Wenqiang Zhang\n');
fprintf(fid,'# date = 2023/02/01\n');
fprintf(fid,'# code = DRDG3D\n');
fprintf(fid,'# code_version = 0.0\n');
fprintf(fid,'# element_size = 200 m on fault, O4\n');
fprintf(fid,'# time_step = %g\n',t(2)-t(1));
fprintf(fid,'# num_time_steps = %d\n',length(t));
fprintf(fid,'# location= 6 km off fault, -12 km along strike, 0 km down-normal\n');
fprintf(fid,'# Column #2 = horizontal displacement (m)\n');
fprintf(fid,'# Column #3 = horizontal velocity (m/s)\n');
fprintf(fid,'# Column #4 = vertical displacement (m)\n');
fprintf(fid,'# Column #5 = vertical velocity (m/s)\n');
fprintf(fid,'# Column #6 = normal displacement (m)\n');
fprintf(fid,'# Column #7 = normal velocity (m/s)\n');
fprintf(fid,'#\n');
fprintf(fid,'# The line below lists the names of the data fields:\n');
fprintf(fid,'t h-disp h-vel v-disp v-vel n-disp n-vel\n');
fprintf(fid,'#\n');

v = v([5,2,6,3,4,1],:);
dat = [t,v'];
fprintf(fid,'%20.12e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n',dat');
fclose(fid);
