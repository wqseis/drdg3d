function check_edge_length(fnm)

if(nargin < 1)
    fnm = 'mesh.nc';
end

elem = ncread(fnm,'elem');
node = ncread(fnm,'node');

idx = [1,2;1,3;1,4;2,3;2,4;3,4];

velem = reshape(node(:,elem),[3,size(elem)]);
v1 = reshape(velem(:,idx(:,1),:),3,[]);
v2 = reshape(velem(:,idx(:,2),:),3,[]);
vv = v1-v2;
len = sqrt(vv(1,:).^2+vv(2,:).^2+vv(3,:).^2);
fprintf('edge len = %g ~ %g\n',min(len),max(len));
figure
hist(len)
title('Edge length (km)')
end

