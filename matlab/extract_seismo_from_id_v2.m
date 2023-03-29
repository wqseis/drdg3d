function [ t, v, bc, coord, normal ] = extract_seismo_from_id( data_dir, nproc, id1, varnm, normal1 )

myconstants;

if (nargin < 5)
  normal1 = [-1 0 0]; % positive side of the fault surface
end

ridx = 0;
mpi_id = 0;
for i = 0:nproc-1
    fnm = [data_dir, '/recv_mpi',num2str(i,'%06d'),'.nc'];
    %disp(fnm)
    if(exist(fnm,'file'))
        disp(fnm)
        id = ncread(fnm,'id');
        nor = ncread(fnm,'normal');
        ndotn = zeros(size(id));

        for j = 1:length(id)
            ndotn(j) = dot(nor(:,j),normal1);
        end

        %idx = find(id==id1 & ndotn > 0.707);
        idx = find(id==id1 & ndotn > 0.01);

        if (~isempty(idx))
            ridx = idx(1);
            mpi_id = i;
            break;
        end
    end
end

%idx
%ridx
fnm = [data_dir, '/recv_mpi',num2str(mpi_id,'%06d'),'.nc'];
t = ncread(fnm,'time');
coord = ncread(fnm,'coord',[1,ridx],[3,1]);
bc = ncread(fnm,'bctype',[ridx],[1]);
normal = ncread(fnm,'normal',[1,ridx],[3,1]);

nvar = 1;

if bc >= BC_FAULT
    switch varnm
        case 'ratem'
            ivar = 1;
        case 'ratel'
            ivar = 2;
        case 'taum'
            ivar = 3;
        case 'taul'
            ivar = 4;
        case 'sigma'
            ivar = 5;
        case 'slipm'
            ivar = 6;
        case 'slipl'
            ivar = 7;
        case 'state'
            ivar = 8;
        case 'TP_T'
            ivar = 9;
        case 'TP_P'
            ivar = 10;
        case 'Vx'
            ivar = 11;
        case 'Vy'
            ivar = 12;
        case 'Vz'
            ivar = 13;
        case 'all'
            ivar = 1; nvar = 13;
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
        case 'all'
            ivar = 1; nvar = 6;
        otherwise
            ivar = 1;
    end
end

v = ncread(fnm,'var',[ridx,ivar,1],[1,nvar,Inf]);
v = squeeze(v);
if (bc == BC_FREE && strcmp(varnm,'all'))
    for i = 1:3
        v(i+3,:) = cumtrapz(t,v(i,:));
    end
end

if (bc >= BC_FAULT && strcmp(varnm,'all'))
    v1 = zeros(size(v,1)+3,size(v,2));
    v1(1:13,:) = v;
    for i = 1:3
        v1(i+13,:) = cumtrapz(t,v(i+10,:));
    end
    v = v1;
end

end