function [c1,c2,c3,c4] = barycentricPointInTet(v1,v2,v3,v4,p)
c1 = SignedDistance(v2,v3,v4,p) / ...
    SignedDistance(v2,v3,v4,v1);
c2 = SignedDistance(v1,v3,v4,p) / ...
    SignedDistance(v1,v3,v4,v2);
c3 = SignedDistance(v1,v2,v4,p) / ...
    SignedDistance(v1,v2,v4,v3);
c4 = SignedDistance(v1,v2,v3,p) / ...
    SignedDistance(v1,v2,v3,v4);
end