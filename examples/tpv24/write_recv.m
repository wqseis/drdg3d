clc
clear
%close all

addmypath

for id = 1:14

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

[ t, v, bc, coord, nor ] = extract_seismo_from_id( data_dir, nproc, id, 'all');
v = v([6,1,3,7,2,4,5,8,9,10],:);
%v(7,:) = -v(7,:);
v = v(1:7,:);
bc
coord
nor
if(bc==100)
  strike = round(coord(2)/cosd(30));
  dip = -(coord(3));
  prefix = 'branch';
else
  strike = coord(2);
  dip = -coord(3);
  prefix = 'fault';
end

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

fmt = [prefix,'st',fmt1,'dp',fmt2];
fnm = sprintf(fmt,strike*10,dip*10);
fid = fopen(fnm,'wt');
fprintf(fid,'# problem = TPV24\n');
fprintf(fid,'# author = Wenqiang Zhang\n');
fprintf(fid,'# date = 2023/02/07\n');
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
fprintf(fid,'# The line below lists the names of the data fields:\n');
fprintf(fid,'t h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress\n');
fprintf(fid,'# Here is the time-series data.\n');
dat = [t,v'];
fprintf(fid,'%21.13e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n',dat');
fclose(fid);

end
