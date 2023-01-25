function fault_tris = get_fault_connectivity(elem,fault2wave,bctype)

BC_FAULT = 100;

FToV = [ ...
    1,2,3;
    1,2,4;
    2,3,4;
    1,3,4];
 
nfault_elem = length(fault2wave);

% counting fault faces ...')

num_fault_face = 0;

for ief = 1:nfault_elem
    ie = fault2wave(ief);
    for is = 1:4
        if(bctype(is,ie)>=BC_FAULT)
            num_fault_face = num_fault_face + 1;
        end
    end
end

nfault = num_fault_face;
tri_dense = zeros(nfault,3);
 
k = 1;
% for ie = 1:nelem
%     ief = wave2fault(ie);
for ief = 1:nfault_elem
    ie = fault2wave(ief);
    for is = 1:4
        if (bctype(is,ie)>=BC_FAULT)
            %tri_dense(k,:) = [1,2,3] + 3*(k-1); 
            tri_dense(k,:) = elem(ie,FToV(is,:));
            k = k + 1;
        end
    end
end


fault_tris = tri_dense;

end

