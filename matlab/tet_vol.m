function vol = tet_vol(A,B,C,D)

vol = dot(B-A,cross(C-A,D-A))/6;

end