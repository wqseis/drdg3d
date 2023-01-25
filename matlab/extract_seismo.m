function [ t, v, coord1 ] = extract_seismo( data_dir, nproc, coord, bc, varnm )

myconstants;

rmin = 1e30;
ridx = 0;
mpi_id = 0;
for i = 0:nproc-1
    fnm = [data_dir, '/recv_mpi',num2str(i,'%06d'),'.nc'];
    if(exist(fnm,'file'))
        c = ncread(fnm,'coord');
        %[~,nr] = size(c);
        
        r = sqrt((c(1,:)-coord(1)).^2+...
            (c(2,:)-coord(2)).^2+...
            (c(3,:)-coord(3)).^2);
        [val,idx] = min(r(:));
    end
    if(val<rmin)
        rmin = val;
        ridx = idx;
        mpi_id = i;
        
    end
end

fnm = ['data/recv_mpi',num2str(mpi_id,'%06d'),'.nc'];

if bc == BC_FAULT
    switch varnm
        case 'rate'
            ivar = 1;
        case 'stress';
            ivar = 2;
        case 'sigma';
            ivar = 3;
        case 'slip';
            ivar = 4;
        otherwise
            ivar = 1;
    end
end
if bc == BC_FREE
    switch varnm
        case 'Vx'
            ivar = 1;
        case 'Vy'
            ivar = 2;
        case 'Vz'
            ivar = 3;
        otherwise
            ivar = 1;
    end
end

t = ncread(fnm,'time');
v = ncread(fnm,'var',[ridx,ivar,1],[1,1,Inf]);
%bc = ncread(fnm,'bctype');
coord1 = ncread(fnm,'coord',[1,ridx],[3,1]);

%v1 = v(1,1,:);
v = squeeze(v);
end
