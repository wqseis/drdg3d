clc
clear
%close all
addmypath


par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

for id = 1:9

[ t, v, bc, coord, nor ] = extract_seismo_from_id( data_dir, nproc, id, 'all');
v = v([6,1,3,7,2,4,5,8,9,10],:);
%v(7,:) = -v(7,:);
v = v(1:8,:);
f0 = 0.6;
b = 0.012;
Dc = 0.02;
V0 = 1e-6;
%v = Dc/V0*exp((v-f0)/b);
%v = log10(v);
v1 = v(8,:);
v1 = log10(Dc/V0)+log10(exp(1))*(v1-f0)/b;
v(8,:) = v1;

strike = coord(2);
dip = -coord(3);

if (strike>-1e-6)
    fmt1='%03d';
else
    fmt1='%04d';
end
if (dip>-1e-6)
    fmt2='%03d';
else
    fmt2='%04d';
end
fmt = ['faultst',fmt1,'dp',fmt2];
fnm = sprintf(fmt,strike*10,dip*10)
fid = fopen(fnm,'wt');
fprintf(fid,'# problem = TPV102\n');
fprintf(fid,'# author = Wenqiang Zhang\n');
fprintf(fid,'# date = 2023/02/13\n');
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
fprintf(fid,'# Column #9 = log10(theta)\n');
fprintf(fid,'# The line below lists the names of the data fields:\n');
fprintf(fid,'t h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress log-theta\n');
fprintf(fid,'# Here is the time-series data.\n');
dat = [t,v'];
fprintf(fid,'%20.12e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n',dat');
fclose(fid);
end
