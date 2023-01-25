function [Tn] = rotate_tensor(xo,yo,zo,xn,yn,zn,To)
% To: old tensor; Tn: new tensor;
% e.g., unit vector xo = [0,0,1];
% reference: https://slideplayer.com/slide/14088614/

R = [dot(xo,xn),dot(xo,yn),dot(xo,zn);
    dot(yo,xn),dot(yo,yn),dot(yo,zn);
    dot(zo,xn),dot(zo,yn),dot(zo,zn)];

Tn = R'*To*R;
