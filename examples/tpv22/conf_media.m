%clc
%clear
%close all

fnm_out = 'mesh.nc';

node = ncread(fnm_out,'node');
elem = ncread(fnm_out,'elem');

x_tet = reshape(node(1,elem),size(elem));
xc = mean(x_tet); % center of tetrahedron

ncid = netcdf.open(fnm_out,'WRITE');

dimid_elem = netcdf.inqDimID(ncid,'Nelem');
[~,num_elem] = netcdf.inqDim(ncid,dimid_elem);

var1 = netcdf.inqVarID(ncid,'rho');
var2 = netcdf.inqVarID(ncid,'vp');
var3 = netcdf.inqVarID(ncid,'vs');

vp = zeros(num_elem,1);
vs = zeros(num_elem,1);
rho = zeros(num_elem,1);
vp(:) = 6;
vs(:) = 3.464;
rho(:) = 2.67;

if 0
% for tpv6, bimaterial fault
idx = find(xc<0);
vp(idx) = 6/1.6;
vs(idx) = 3.464/1.6;
rho(idx) = 2.67/1.2;
end


netcdf.putVar(ncid,var1,rho);
netcdf.putVar(ncid,var2,vp);
netcdf.putVar(ncid,var3,vs);

netcdf.close(ncid);

