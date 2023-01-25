function [ftype] = set_fluxtype_from_faultnodes(elem,fnodes)

fnodes = reshape(fnodes,[],1);
fnodes = unique(sort(fnodes));

Nelem = size(elem,2); % elem(4,nelem);
Nnode = max(elem(:));

ftype = zeros(4,Nelem);

FToV = [ ...
    1,2,3;
    1,2,4;
    2,3,4;
    1,3,4];

node_flag = zeros(Nnode,1);
node_flag(fnodes)=-1;
% for i = 1:length(fnodes)
%    j = fnodes(i);
%    node_flag(j) = -1;
% end


fn=node_flag(elem(FToV',:));
fn=reshape(fn,[3,4,Nelem]);
ftype(:,:) = sum(fn,1);

ftype1=zeros(size(ftype(:)));
ftype1(ftype(:)==-1)=1;
ftype1(ftype(:)==-2)=1;
ftype=reshape(ftype1,size(ftype));


% bcytype(:,:)=0;
% for ie = 1:Nelem
%    for is = 1:4
%        fnode = elem(FToV(is,1:3),ie);
%        ftype(is,ie) = sum(node_flag(fnode));
%    end
% end
 
end
