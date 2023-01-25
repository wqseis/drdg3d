%clc
%clear
%close all

addmypath;

par = ReadYaml('parameters.yaml');
nproc = par.nproc;

fnm_out = 'mesh.nc';

elem = ncread(fnm_out,'elem');
Nelem = size(elem,2);

if nproc > 1
    part = metis_part(elem,nproc,metis_dir) + 1;
else
    part = zeros(Nelem,1)+1;
end

disp('modified part ...')

ncid = netcdf.open(fnm_out,'WRITE');
var1 = netcdf.inqVarID(ncid,'part');
netcdf.putVar(ncid,var1,part);
netcdf.close(ncid);
