function [ t ] = gather_fault_time( data_dir, nproc )

v = [];
for iproc = 0:nproc-1
    fnm = sprintf('%s/fault_mpi%06d.nc',data_dir,iproc);
    if (exist(fnm,'file'))
        t = ncread(fnm, 'time');
        return
    end
end

end
