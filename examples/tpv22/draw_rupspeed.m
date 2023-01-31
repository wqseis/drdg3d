clc
clear
close all

addmypath;

par = ReadYaml('parameters.yaml');

nproc = par.nproc;
data_dir = par.data_dir;

x = gather_fault_var( data_dir, nproc, 'x' );
y = gather_fault_var( data_dir, nproc, 'y' );
z = gather_fault_var( data_dir, nproc, 'z' );
nx = gather_fault_var( data_dir, nproc, 'nx' );
t = gather_fault_var( data_dir, nproc, 'ruptime' );
fault_id = gather_fault_var1d( data_dir, nproc, 'fault_id' );

idx = find(mean(nx(:,:))>0 & fault_id == 3);
x = x(:,idx);
y = y(:,idx);
z = z(:,idx);
t = t(:,idx);

% r = sqrt((y-6).^2+(z+8.5).^2);
% r = r/3.464;
% %t = r + 0.01*randn(size(t));

x = mean(x);
y = mean(y);
z = mean(z);
t = mean(t);

%tri = get_face_connect(x);

[zx,zy] = trigradient(y(:),z(:),t(:) );
v=sqrt(zx.^2+zy.^2);
%v = abs(zx);
v = 1./(v+1e-16);

v = v/3.464;

yy=linspace(min(y(:)),max(y(:)),200);
zz=linspace(min(z(:)),max(z(:)),100);
[yy,zz]=meshgrid(yy,zz);

F = scatteredInterpolant(y(:),z(:),v(:));
vv = F(yy,zz);
F = scatteredInterpolant(y(:),z(:),t(:));
tt = F(yy,zz);
figure
%trisurf(tri,x(:),y(:),z(:),v(:) )
%scatter3(x,y,z,60,v,'filled')
pcolor(yy,zz,vv)
hold on
contour(yy,zz,tt,0:0.5:100,'color','k')
hold off
colorbar
axis image
shading interp
colormap jet
caxis([0 2])