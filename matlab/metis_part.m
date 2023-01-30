function part = metis_part(elem,np,metis_dir)

if nargin < 3
    metis_dir = '/usr/local/bin';
    metis_dir = '/home/wzhang/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-9.4.0/metis-5.1.0-g3hkcjmmfvljvil6u6okitnlxozpw3mn/bin';
end

metis_version = 5; % choose metis version 4 or 5

% https://people.sc.fsu.edu/~jburkardt/data/metis_mesh/metis_mesh.html
% ASCII
% a mesh of N elements is stored in a file of N+1 lines;
% the first line lists the number of elements, and a code that indicates the type of element;
% 1: 2D triangles;
% 2: 3D tetrahedrons;
% 3: 3D hexahedrons ("bricks" with 6 sides and 8 vertices);
% 4: 2D quadrilaterals;
% each subsequent line lists the vertices of one element;
% for triangular and tetrahedral elements, the vertices may be listed in any order; quadrilateral and hexahedral elements require that the vertices be listed in a particular order;
% comment lines begin with a "%" sign;

nelem = size(elem,2);
fid = fopen('tmp.met','w');
if (metis_version == 4)
    fprintf(fid,'%d %d\n',nelem,2);
elseif (metis_version == 5)
    % in metis-5, the second number in the first line is weight
    % (can be omitted)
    fprintf(fid,'%d\n',nelem);
end

%fprintf(fid,'%d %d %d %d\n',elem');
fprintf(fid,'%d %d %d %d\n',elem);
fclose(fid);

if (metis_version == 4)
    %metis_dir = '../metis-4.0.3';
    cmd = sprintf('%s/mesh2dual tmp.met',metis_dir);
    system(cmd);
    cmd = sprintf('%s/pmetis tmp.met.dgraph %d',metis_dir,np);
    system(cmd);
    part = load(['tmp.met.dgraph.part.',num2str(np)]);
elseif (metis_version == 5)
    %metis_dir = '/usr/local/bin';
    %metis_dir = '/home/wzhang/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-9.4.0/metis-5.1.0-g3hkcjmmfvljvil6u6okitnlxozpw3mn/bin';

    exe = strtrim(fileread('../../path.m2gmetis'));
    %cmd = sprintf('%s/m2gmetis tmp.met tmp.met.dgraph',metis_dir);
    cmd = sprintf('%s tmp.met tmp.met.dgraph',exe);
    disp(cmd);
    system(cmd);

    exe = strtrim(fileread('../../path.gpmetis'));
    %cmd = sprintf('%s/gpmetis tmp.met.dgraph %d',mesh_dir,np);
    cmd = sprintf('%s tmp.met.dgraph %d',exe,np);
    disp(cmd);
    system(cmd);

    part = load(['tmp.met.dgraph.part.',num2str(np)]);
end

% part number may start from 0 or 1 depend on the version, so here I force
% it start from 0
part = part - min(part); % start from 0

% clean up
system(['rm -rf tmp.met tmp.met.dgraph tmp.met.dgraph.part.',num2str(np)]);
end
