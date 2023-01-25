function [nx,ny,nz] = calcu_normal_vectors(elem,node)
%[nx,ny,nz] = calcu_normal_vectors(elem,node)

% FtoV = [...
%     1,2,3;
%     1,2,4;
%     2,3,4;
%     1,3,4];
% point outward
FtoV = [...
    1,3,2;
    1,2,4;
    2,3,4;
    1,4,3;
    ];

nelem = size(elem,2);

a =   (node(2,elem(FtoV(:,2),:))-node(2,elem(FtoV(:,1),:))) ...
    .*(node(3,elem(FtoV(:,3),:))-node(3,elem(FtoV(:,1),:))) ...
    - (node(2,elem(FtoV(:,3),:))-node(2,elem(FtoV(:,1),:))) ...
    .*(node(3,elem(FtoV(:,2),:))-node(3,elem(FtoV(:,1),:)));
b =   (node(1,elem(FtoV(:,3),:))-node(1,elem(FtoV(:,1),:))) ...
    .*(node(3,elem(FtoV(:,2),:))-node(3,elem(FtoV(:,1),:))) ...
    - (node(1,elem(FtoV(:,2),:))-node(1,elem(FtoV(:,1),:))) ...
    .*(node(3,elem(FtoV(:,3),:))-node(3,elem(FtoV(:,1),:)));
c =   (node(1,elem(FtoV(:,2),:))-node(1,elem(FtoV(:,1),:))) ...
    .*(node(2,elem(FtoV(:,3),:))-node(2,elem(FtoV(:,1),:))) ...
    - (node(2,elem(FtoV(:,2),:))-node(2,elem(FtoV(:,1),:))) ...
    .*(node(1,elem(FtoV(:,3),:))-node(1,elem(FtoV(:,1),:)));

a = reshape(a,4,nelem);
b = reshape(b,4,nelem);
c = reshape(c,4,nelem);
Fscale = sqrt(a.*a+b.*b+c.*c);
nx = a./Fscale;
ny = b./Fscale;
nz = c./Fscale;

end 
