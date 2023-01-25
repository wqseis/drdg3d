function d = SignedDistance(v1,v2,v3,p)

n = cross(v2-v1, v3-v1);
d = dot(n, p-v1);

end
