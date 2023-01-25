% clc
clear
% close all

addmypath;

par = ReadYaml('parameters.yaml');
nproc = par.nproc;
data_dir = par.data_dir;

varnm = 'Uy';

x = gather_grdsurf_var( data_dir, nproc, 'x' );
y = gather_grdsurf_var( data_dir, nproc, 'y' );
z = gather_grdsurf_var( data_dir, nproc, 'z' );

tri1 = get_face_connect(x);

h=figure;
%set(h,'Visible','off')
%for it = 10:10:100
for it =  1 : 1 : 100000
  [v,t] = gather_grdsurf_snap(data_dir,nproc,varnm,it);
  trisurf(tri1,x(:),y(:),z(:),v(:))
  view(2)
  shading interp
  axis image
  xlabel('X (km)')
  ylabel('Y (km)')
  colorbar
  vmax=max(abs(v(:)));
  vmax=max(vmax,1e-16);
  %vmax=.1;
  caxis([-1 1]*vmax/2)
  
  colormap  coolwarm_ext
  %caxis([-1 1]*1)
  title([varnm, ' T(',num2str(it),') = ',num2str(t),' s'])
  pause(0.002)

  if 0
     print('-dpng','-r300',...
            ['grd_',varnm,'_it_',num2str(it,'%06d'),'.png'])
  end

  %saveas(h,'grdsurf_snap.png')
end
