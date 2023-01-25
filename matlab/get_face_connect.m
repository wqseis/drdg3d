function [tri1] = get_face_connect(x)

[Nfp,Nelem]=size(x);
Order = fix(sqrt(2*Nfp-1/4)-1/2);
tris = tris_per_face(Order+1);
ntris = size(tris,1);
tri1 = zeros(ntris*Nelem,3);
for i = 1:Nelem
    tri = tris + (i-1)*Nfp;
    tri1(((i-1)*ntris+1):i*ntris,:) = tri;
end
end

