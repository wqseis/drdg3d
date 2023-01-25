clc
clear

addmypath;

fnm = 'mesh.nc';

node = ncread(fnm,'node');
elem = ncread(fnm,'elem');
fault2wave = ncread(fnm,'fault2wave');
neighbor = ncread(fnm,'neighbor');
bctype = ncread(fnm,'bctype');


nelem = size(elem,2);
nfault_elem = length(fault2wave);

r = zeros(nfault_elem,1);


for ief = 1:nfault_elem

    ie = fault2wave(ief);
    A = node(:,elem(1,ie));
    B = node(:,elem(2,ie));
    C = node(:,elem(3,ie));
    D = node(:,elem(4,ie));

    vol1 = tet_vol(A,B,C,D);

    % find neighor
    for is = 1:4
        if( bctype(is,ie)>=BC_FAULT )
            ien = neighbor(is,ie);
        end
    end

    A = node(:,elem(1,ien));
    B = node(:,elem(2,ien));
    C = node(:,elem(3,ien));
    D = node(:,elem(4,ien));

    vol2 = tet_vol(A,B,C,D);

    v1=max(vol1,vol2);
    v2=min(vol1,vol2);

    r(ief)=v1/v2;

end

hist(r)
xlabel('volume ratio')
ylabel('number of elements')
