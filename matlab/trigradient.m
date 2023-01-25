function [zx,zy] = trigradient(x,y,z,t,opt)
%   [ZX,ZY] = TRIGRADIENT(X,Y,Z,T) returns an approximation to the gradient
%   of Z defined on the irregular mesh with delaunay triangulation T. 
%   ZX corresponds to the partial derivative dZ/dX, and ZY corresponds to 
%   the partial derivative  dZ/dY. 
%
%   If Z is a matrix of size MxN then ZX and ZY will be computed for each
%   of the N columns of Z.
%
%   [ZX,ZY] = TRIGRADIENT(X,Y,Z) returns an approximation to the gradient
%   of Z for the scattered data (X,Y). 
%
%   [ZX,ZY] = TRIGRADIENT(X,Y,Z,'face') returns the value of the of the
%   gradient for each triangular face assuming z is a piece-wise planar
%   function. ZX(1) is an approximation to the gradient in the x direction
%   of the function Z for the triangle defined by the nodes x(T(1,:)),
%   y(T(1,:)).
%

%   Trigradient uses a first order approximation. Z is assumed to be a 
%   piece-wise planar function over the triangulated region. The partial
%   derivative of each plane (ax + by + c = z) is a in the x direction and 
%   b in the y direction. A given node in the triangulation is associated 
%   with k triangles. The partial derivative at each node is the weighted 
%   sum of each of the partial derivatives of the triangles associated with 
%   that node. The weightings are the areas of each triangle divided by the
%   total area of all triangles associated with that node.
%   
%
%   Examples:
%       [x y] = meshgrid(-2:.1:2,-2:.1:2);
%       x(2:2:end,:)=x(2:2:end,:)+.1/2;     
%       z = x.*exp(-x.^2 - y.^2);         
%       [zx,zy] = trigradient(x(:),y(:),z(:));       
%       dzdx = (1./x-2*x).*z;
%       dzdy = -2*y.*z;       
%       figure(1); hold on
%       quiver(x(:),y(:),dzdx(:),dzdy(:),'r')
%       quiver(x(:),y(:),zx,zy,'k')      
%       contour(x,y,z,10)
%       legend('Exact Gradient','Trigradient')
%
%   See also DELAUNAY, DELAUNAYTRI, TRIMESH

%   Copyright 2013 
%   Mick Warehime
%   $Version 2.1$ $Date: 2013/04/05$

% find the delaunay triangulation of the scaterred data if none provided
if nargin<4; t = delaunay(x,y); opt ='normal';end
p = [x,y];

% catch for no option given.
if nargin<5; opt='normal'; end

% number of triangles and points
nt = size(t,1); np = size(p,1);

% determine the coefficients for the plane above each triangle
% [x1 y1 1; x2 y2 1; x3 y3 1][a; b; c] = [z1; z2; z3] 
C = sparse(repmat(1:3*nt,3,1)',kron(reshape(1:3*nt,3,nt)',[1; 1; 1]),[p(t',:),ones(3*nt,1)])\z(t',:);

% return the derivatives of the triangles if the 'face' opt is selected
if strcmp(opt,'face')
    zx = C(1:3:end,:);
    zy = C(2:3:end,:);
    return
end

% calculate triangle areas using the cross product
areas = repmat(.5*abs((p(t(:,2),1)-p(t(:,1),1)).*(p(t(:,3),2)-p(t(:,1),2))-(p(t(:,3),1)-p(t(:,1),1)).*(p(t(:,2),2)-p(t(:,1),2))),1,size(z,2));

% map matrix from triangles to nodes
M = sparse(repmat((1:nt)',3,1),t,1,nt,np);

% weight the partial derivative (a or b coefficient) by the area of each
% triangle and sum for a given node. divide by the total area 
zx = (M'*(C(1:3:end,:).*areas))./(M'*areas);
zy = (M'*(C(2:3:end,:).*areas))./(M'*areas);

