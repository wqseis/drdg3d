function [ v ] = gather_fault_var( data_dir, nproc, varnm )

v = [];
for iproc = 0:nproc-1
    fnm = sprintf('%s/fault_mpi%06d.nc',data_dir,iproc);
    if (exist(fnm,'file'))
        v1 = ncread(fnm,varnm);
        %v1 = squeeze(v1);
        %size(v1)
        %disp(fnm)
        v = [v,v1];   
    end
end

end