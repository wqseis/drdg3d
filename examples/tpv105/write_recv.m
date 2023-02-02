clc
clear
%close all
addmypath

strike = 18;
dip = 7.5;

bc = BC_FAULT;
coord = [0, strike,-dip];

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

[ t, v, coord1 ] = extract_seismo_all( data_dir, nproc, coord, bc);
v = v([6,1,3,7,2,4,5,8,9,10],:);
v(7,:) = -v(7,:);
v(9,:) = v(9,:) + 483.15;
v(10,:) = v(10,:) + 80;

fnm = sprintf('faultst%03ddp%03d.txt',strike*10,dip*10);
fid = fopen(fnm,'wt');
fprintf(fid,'# problem = TPV105-3D\n');
fprintf(fid,'# author = Wenqiang Zhang\n');
fprintf(fid,'# date = 2023/02/01\n');
fprintf(fid,'# code = DRDG3D\n');
fprintf(fid,'# code_version = 0.0\n');
fprintf(fid,'# element_size = 200 m on fault, O4\n');
fprintf(fid,'# time_step = %g\n',t(2)-t(1));
fprintf(fid,'# num_time_steps = %d\n',length(t));
fprintf(fid,'# location = on fault, %g km along strike, %g km down-dip\n',strike,dip);
fprintf(fid,'# Column #1 = time (s)\n');
fprintf(fid,'# Column #2 = horizontal slip (m)\n');
fprintf(fid,'# Column #3 = horizontal slip rate (m/s)\n');
fprintf(fid,'# Column #4 = horizontal shear stress (MPa)\n');
fprintf(fid,'# Column #5 = vertical slip (m)\n');
fprintf(fid,'# Column #6 = vertical slip rate (m/s)\n');
fprintf(fid,'# Column #7 = vertical shear stress (MPa)\n');
fprintf(fid,'# Column #8 = effective normal stress (MPa)\n');
fprintf(fid,'# Column #9 = state variable psi (dimensionless)\n');
fprintf(fid,'# Column #10 = temperature (K)\n');
fprintf(fid,'# Column #11 = pore pressure (MPa)\n');
fprintf(fid,'# The line below lists the names of the data fields:\n');
fprintf(fid,'t h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress psi temperature pressure\n');
fprintf(fid,'# Here is the time-series data.\n');
dat = [t,v'];
fprintf(fid,'%20.12e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n',dat');
fclose(fid);
