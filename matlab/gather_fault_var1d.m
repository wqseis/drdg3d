function [ v ] = gather_fault_var1d( data_dir, nproc, varnm )

v = [];
for iproc = 0:nproc-1
    fnm = sprintf('%s/fault_mpi%06d.nc',data_dir,iproc);
    if (exist(fnm,'file'))
        v1 = ncread(fnm,varnm);
        v = [v,reshape(v1,1,[])];   
    end
end

end
