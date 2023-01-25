function [wave2fault,fault2wave] = set_fault_tags(bctype)
BC_FAULT = 100;
[~,Nelem] = size(bctype);
%k = 0;
%kf = 0;
%ke = 0;
nfault_elem = 0;
nfault_face = 0;
wave2fault = zeros(Nelem,1);
fault2wave = zeros(Nelem,1);
for ie = 1:Nelem
   isfault = false;
   for is = 1:4
       if (bctype(is,ie)>=BC_FAULT)
           % = k + 1;
           %kf = kf + 1;
           isfault = true;

           nfault_face = nfault_face + 1;
           %fault_elem(k) = ie;
           %fault_face(k) = is;
           %fault_face_elem(1,k) = is;
           %fault_face_elem(2,k) = ie;
           %wave2fault(ie) = k;
           %fault2wave(k) = ie;
       end
   end
   if (isfault)
       %ke = ke + 1;
       nfault_elem = nfault_elem + 1;
       wave2fault(ie) = nfault_elem;
       fault2wave(nfault_elem) = ie;       
   end
end 
fault2wave = fault2wave(1:nfault_elem);
end

