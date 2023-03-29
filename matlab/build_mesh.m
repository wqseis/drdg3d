function msh = build_mesh(fnm_in,tet_id,faultsurf_id,freesurf_id)

myconstants;
%num_fault_segments = 3;
num_fault_segments = length(faultsurf_id);
num_free_segments = length(freesurf_id);

%fnm_in = 'stepover3.exo';
%fnm_out = 'mesh.nc';

coord = ncread(fnm_in, 'coord');
%coord = coord * 1e-3; % Add by xxt@20230224, for the present mesh is in m, but not km.

num_tets = length(tet_id);

connect = [];
for i = 1:num_tets
    connect1 = ncread(fnm_in,['connect',num2str(tet_id(i))]);
    connect = [connect,connect1];
end
    

%connect2 = ncread(fnm,'connect51');
%connect3 = ncread(fnm,'connect52');
%connect = [connect1,connect2,connect3];
%connect = connect1;

% if there are > 2 connects
% connect = [connect1,connect2,...,connectn];


% find fault and free surface in *.exo file:
% multiple fault surfaces:
cf = cell(num_fault_segments);
for i = 1:num_fault_segments
    cf{i} = ncread(fnm_in,['connect',num2str(faultsurf_id(i))]); 
%cf{2} = ncread(fnm_in,'connect28');
%cf{3} = ncread(fnm_in,'connect35');
end

if num_free_segments > 0
% free surface

fr = [];
for i = 1:num_free_segments
    fr1 = ncread(fnm_in,['connect',num2str(freesurf_id(i))]);
    fr = [fr,fr1];
end

%fr = ncread(fnm_in,['connect',num2str(freesurf_id)]);


freenodes = unique(sort(fr(:)));

else
    freenodes = [];

end

%fault = unique(sort(ft(:)));
fault = [];
for i = 1:num_fault_segments
    f1 = unique(sort(cf{i}(:)));
    fault = [fault;f1];
end
fault = unique(sort(fault(:)));
 
elem = double(connect');
node = coord;

%x = node(:,1);
%y = node(:,2);
%z = node(:,3);

Nelem = size(elem,1);
Nfault = size(fault,1);
Nnode = size(node,1);

fnodes = fault(:);
fnodes = unique(sort(fnodes));
Nfnodes = length(fnodes);

msh.Nelem = Nelem;
msh.Nnode = Nnode;
msh.Nnode_fault = length(fnodes);
msh.elem = elem';
msh.node = node';
msh.node_fault = fnodes;

msh.neighbor  = zeros(4,Nelem);
msh.face      = zeros(4,Nelem);
msh.direction = zeros(4,Nelem);
msh.bctype    = zeros(4,Nelem);
msh.elemtype  = zeros(1,Nelem);
msh.vp        = zeros(1,Nelem);
msh.vs        = zeros(1,Nelem);
msh.rho       = zeros(1,Nelem);
msh.part      = zeros(1,Nelem);

disp('set bctype ...')
%tic

if num_fault_segments >= 1
% label fault surface as: 100,101,102,...
fnode1 = unique(sort(cf{1}(:)));
bc1 = set_bctype_from_nodes(msh.elem,fnode1,BC_FAULT);
for i = 2:num_fault_segments
    fnode1 = unique(sort(cf{i}(:)));
    bc2 = set_bctype_from_nodes(msh.elem,fnode1,BC_FAULT+(i-1));
    bc1 = bc1+bc2;
end
else
    bc1 = [];
end

%if num_free_segments > 0
% label free surface
bc2 = set_bctype_from_nodes(msh.elem,freenodes,BC_FREE);

%end

if ~isempty(bc1)
    msh.bctype = bc1+bc2;
else
    msh.bctype = bc2;
end

fnodes = fault(:);
fnodes = unique(sort(fnodes));
%size(fnodes)
msh.fluxtype = set_fluxtype_from_faultnodes(msh.elem,fnodes);

%toc

% set fault to wave indexing
[wave2fault,fault2wave] = set_fault_tags(msh.bctype);
msh.nfault_elem = length(fault2wave);
msh.wave2fault = wave2fault;
msh.fault2wave = fault2wave;
 
end
