function [p,t] = check_surface_id(fnm,surface_id,ifplot)
% check_surface_id(filename,surface_id,ifplot)
% check_surface_id('filename.exo',1)
% check_surface_id('filename.exo',1,0)
% [p,t] = check_surface_id('filename.exo',1)

if(nargin<3)
    ifplot = 1;
end
% if (nargin<2)
%     surface_id = 1;
% end
% if (nargin<1)
%     fnm = 'stepover3.exo';
% end

t = ncread(fnm,['connect',num2str(surface_id)])';

if(size(t,2)==3)
    fprintf('there are %d triangles\n',size(t,1));
elseif (size(t,2)==4)
    fprintf('there are %d tetrahedrons\n',size(t,1));
else
    error(['connect',num2str(surface_id),' is not a surface! check: ncdisp(''',fnm,''')'])
end

p = ncread(fnm,'coord');

x=p(t,1);
y=p(t,2);
z=p(t,3);
x1=min(x(:));
y1=min(y(:));
z1=min(z(:));
x2=max(x(:));
y2=max(y(:));
z2=max(z(:));

fprintf('x = %g ~ %g\n',x1,x2);
fprintf('y = %g ~ %g\n',y1,y2);
fprintf('z = %g ~ %g\n',z1,z2);

if ifplot
    figure

    if (size(t,2)==4)
        tet = double(t);
        X = double(p);
        trep = triangulation(tet, X);
        [tri,xf] = freeBoundary(trep);
        trisurf(tri, xf(:,1),xf(:,2),xf(:,3), ...
            'FaceColor', 'cyan', 'FaceAlpha', .9);
        title(['connnect ',num2str(surface_id)])
    end
    if (size(t,2)==3)
        trisurf(t,p(:,1),p(:,2),p(:,3))
        title(['Surface ',num2str(surface_id)])
    end
    xlabel('X (km)')
    ylabel('Y (km)')
    zlabel('Z (km)')
end
end
