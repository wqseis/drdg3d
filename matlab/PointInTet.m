function flag = PointInTet(v1,v2,v3,v4,p)

% flag = PointInTet(v1,v2,v3,v4,p)
% https://people.sc.fsu.edu/~jburkardt/presentations/cg_lab_barycentric_tetrahedrons.pdf

% flag = ...
%     SameSide(v1,v2,v3,v4,p) && ...
%     SameSide(v2,v3,v4,v1,p) && ...
%     SameSide(v3,v4,v1,v2,p) && ...
%     SameSide(v4,v1,v2,v3,p);
[c1,c2,c3,c4] = barycentricPointInTet(v1,v2,v3,v4,p);
flag = (c1>=0 && c2>=0 && c3>=0 && c4>=0);

% flag2 = ...
%     DistPointToTri(v1,v2,v3,p) < 1e-3 || ...
%     PointInTri(v2,v3,v4,p) || ...
%     PointInTri(v3,v4,v1,p) || ...
%     PointInTri(v4,v1,v2,p);
% 
% flag = flag || flag2;

end

function flag = SameSide(v1,v2,v3,v4,p)

n = cross(v2-v1, v3-v1);
dotV4 = dot(n, v4-v1);
dotP = dot(n, p-v1);

flag = (mysign(dotV4)==mysign(dotP));

end

function d = SignedDistance(v1,v2,v3,p)

n = cross(v2-v1, v3-v1);
d = dot(n, p-v1);

end

function dist = DistPointToTri(v1,v2,v3,p)
n = cross(v2-v1, v3-v1);
dist = dot(n, p-v1);
end

function y = mysign(x)
y = sign(x);
if(y==0)
    y=1;
end
end
