function [ t, v, coord1 ] = extract_seismo_all( data_dir, nproc, coord, bc)

myconstants;

rmin = 1e30;
ridx = 0;
mpi_id = 0;
for i = 0:nproc-1
    fnm = [data_dir, '/recv_mpi',num2str(i,'%06d'),'.nc'];
    %disp(fnm)
    if(exist(fnm,'file'))
        %disp(fnm)
        c = ncread(fnm,'coord');
        %[~,nr] = size(c);
        
        r = sqrt(...
            (c(1,:)-coord(1)).^2+...
            (c(2,:)-coord(2)).^2+...
            (c(3,:)-coord(3)).^2 );
        [val,idx] = min(r(:));

        if(val<rmin)
            rmin = val;
            ridx = idx;
            mpi_id = i;
        end
    end
end

fnm = [data_dir, '/recv_mpi',num2str(mpi_id,'%06d'),'.nc'];

t = ncread(fnm,'time');
v = ncread(fnm,'var',[ridx,1,1],[1,Inf,Inf]);
%bc = ncread(fnm,'bctype');
coord1 = ncread(fnm,'coord',[1,ridx],[3,1]);

%v1 = v(1,1,:);
v = squeeze(v);

if (bc >= BC_FAULT)
%v = v([6,1,3,7,2,4,5,8,9,10],:);
end
if (bc == BC_FREE)
    for i = 1:3
        v(i+3,:) = cumtrapz(t,v(i,:));
    end
end
end
