function [ v ] = gather_grdsurf_var( data_dir, nproc, varnm )

v = [];
for iproc = 0:nproc-1
    fnm = sprintf('%s/grdsurf_mpi%06d.nc',data_dir,iproc);
    if (exist(fnm,'file'))
        v1 = ncread(fnm,varnm);
        v = [v,v1];   
    end
end

end
