function [bctype] = set_bctype_from_nodes(elem,fnodes,id)

if (nargin <3 )
    id = -3;
end

Nelem = size(elem,2); % elem(4,nelem);
Nnode = max(elem(:));

bctype = zeros(4,Nelem);

if(isempty(bctype))
    return
end

FToV = [ ...
    1,2,3;
    1,2,4;
    2,3,4;
    1,3,4];

node_flag = zeros(Nnode,1);
%for i = 1:length(fnodes)
%    j = fnodes(i);
%    node_flag(j) = -1;
%end
node_flag(fnodes)=-1;

fn=node_flag(elem(FToV',:));
fn=reshape(fn,[3,4,Nelem]);
bctype(:,:) = sum(fn,1);

bctype1 = zeros(size(bctype));
bctype1(bctype==-3) = id;

bctype = bctype1;

%bctype = bctype / (-3);
%bctype = bctype * id;
%for ie = 1:Nelem
%    for is = 1:4
%        fnode = elem(FToV(is,1:3),ie);
%        bctype(is,ie) = sum(node_flag(fnode));
%    end
%end
 
end
