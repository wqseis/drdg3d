% clc
clear
% close all

addmypath;

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

varnm = 'PGVh';

x = gather_grdsurf_var( data_dir, nproc, 'x' );
y = gather_grdsurf_var( data_dir, nproc, 'y' );
z = gather_grdsurf_var( data_dir, nproc, 'z' );
v = gather_grdsurf_var( data_dir, nproc, varnm );

tri = get_face_connect(x);

h=figure;
trisurf(tri,x(:),y(:),z(:),log10(v(:)+1e-3))
view(2)
shading interp
axis image
colormap jet
xlabel('X (km)')
ylabel('Y (km)')
colorbar
 
title([varnm])

if 0
    print('-dpng','-r300',...
        ['grd_',varnm,'.png'])
end

