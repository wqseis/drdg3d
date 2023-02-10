%clc
clear
close all

addmypath;

check_stress = 1;

fnm_out = 'mesh.nc';

node = ncread(fnm_out,'node');
elem = ncread(fnm_out,'elem');
nelem = size(elem,2);
nnode = size(node,2);

fault2wave = ncread(fnm_out,'fault2wave');
wave2fault = ncread(fnm_out,'wave2fault');

nfault_elem = length(fault2wave);

[nx,ny,nz] = calcu_normal_vectors(elem,node);

Tx = zeros(3,4,nfault_elem);
Ty = zeros(3,4,nfault_elem);
Tz = zeros(3,4,nfault_elem);
dTx = zeros(3,4,nfault_elem);
dTy = zeros(3,4,nfault_elem);
dTz = zeros(3,4,nfault_elem);
mu_s = zeros(3,4,nfault_elem);
mu_d = zeros(3,4,nfault_elem);
Dc = zeros(3,4,nfault_elem);
C0 = zeros(3,4,nfault_elem);
temp = zeros(3,4,nfault_elem);

bctype = ncread(fnm_out,'bctype');

disp('counting fault faces ...')
tic
num_fault_face = 0;
for ie = 1:nelem
    for is = 1:4
        if(bctype(is,ie)>=BC_FAULT)
            num_fault_face = num_fault_face + 1;
        end
    end
end
toc
%
% disp('assigning fault triangulates ...')
% tic
% tri_fault = zeros(num_fault_face,3);
% k = 1;
% for ie = 1:nelem
%     for is = 1:4
%         if (bctype(is,ie)==BC_FAULT)
%             tri =  (elem(FtoV(is,:),ie));
%             tri_fault(k,:) = tri;
%             k = k + 1;
%             %fprintf('%g %g %g\n',nx(is,ie),ny(is,ie),nz(is,ie));
%         end
%     end
% end
% toc

mu_s1 = 0.7;
mu_d1 = 0.5;
%C01 = 0.2;
C01 = 0;
Dc1 = 0.5;

dip_angle = 15;
vec_x0 = [1,0,0];
vec_y0 = [0,1,0];
vec_z0 = [0,0,1];
vec_n0 = [sind(dip_angle),0,cosd(dip_angle)];
vec_m0 = [0,1,0];
vec_l0 = cross(vec_n0,vec_m0);
stress0_bak0 = [-120,0,70;
    0,0,0;
    70 0 0];
stress0_bak = rotate_tensor(vec_n0,vec_m0,vec_l0,vec_x0,vec_y0,vec_z0,stress0_bak0);
stress0_asp0 = [-120,0,85;
    0 0 0;
    85 0 0];
stress0_asp = rotate_tensor(vec_n0,vec_m0,vec_l0,vec_x0,vec_y0,vec_z0,stress0_asp0); 


for ie = 1:nelem
    ief = wave2fault(ie);
    for is = 1:4
        if (bctype(is,ie)>=BC_FAULT)
            xc = mean(node(1,elem(FtoV(is,:),ie)));
            yc = mean(node(2,elem(FtoV(is,:),ie)));
            zc = mean(node(3,elem(FtoV(is,:),ie)));

            for i = 1:3
                j = FtoV(is,i);
                x = node(1,elem(j,ie));
                y = node(2,elem(j,ie));
                z = node(3,elem(j,ie));

                T1 = stress0_bak;
                if (abs(yc-0)<=1.5 && ...
                        abs(zc+7.5*sind(dip_angle))<=1.5*sind(dip_angle))
                    T1 = stress0_asp;
                end

                % Coordinate transforms
                % rotations of the coordinate system 
%                 % while the object is held constant. 
%                  angle = -15;
%                 sina=sind(angle);cosa=cosd(angle);
%                 Rot= [ cosa sina 0;
%                       -sina cosa 0;
%                           0    0 1];
% 
%                  Rot= [ cosa   0   sina;
%                         0      0    0;
%                         -sina  1    cosa];
%
%                 T2=Rot*T1*transpose(Rot);

                vec_n = [nx(is,ie);ny(is,ie);nz(is,ie)];

                % linear 0 to 70 MPa
                T2 = T1 * (x/(15*cosd(dip_angle)));
                %T2 = T1;

                traction0 = T2 * vec_n;

                tn=dot(traction0,vec_n);
                ts_vec=traction0-tn*vec_n;
                ts=norm(ts_vec);
                tp = mu_s1*(-tn); % peak stress
                td = mu_d1*(-tn);

                Tx(i,is,ief) = traction0(1);
                Ty(i,is,ief) = traction0(2);
                Tz(i,is,ief) = traction0(3);

                temp(i,is,ief) = ts;%(tp-ts)./(ts-td);% (ts-td)/(tp-td);
                mu_s(i,is,ief) = mu_s1;
                mu_d(i,is,ief) = mu_d1;
                Dc(i,is,ief) = Dc1;
                C0(i,is,ief) = C01;
            end
        end
    end
end


nfault = num_fault_face;
tri_dense = zeros(nfault,3);
v_dense = zeros(nfault,3);
x_dense = zeros(nfault,3);
y_dense = zeros(nfault,3);
z_dense = zeros(nfault,3);
nx_dense = zeros(nfault,1);
k = 1;
for ie = 1:nelem
    ief = wave2fault(ie);
    for is = 1:4
        if (bctype(is,ie)>=BC_FAULT)
            tri_dense(k,:) = [1,2,3] + 3*(k-1);
            x_dense(k,:) = node(1,elem(FtoV(is,1:3),ie));
            y_dense(k,:) = node(2,elem(FtoV(is,1:3),ie));
            z_dense(k,:) = node(3,elem(FtoV(is,1:3),ie));
            nx_dense(k,:) = nx(is,ie);
             if (nx(is,ie)>0)
                 v_dense(k,:) = temp(1:3,is,ief);
             else
                 v_dense(k,:) = nan;
             end
            k = k + 1;
        end
    end
end

ncid = netcdf.open(fnm_out,'WRITE');

var1 = netcdf.inqVarID(ncid,'Tx0');
var2 = netcdf.inqVarID(ncid,'Ty0');
var3 = netcdf.inqVarID(ncid,'Tz0');
var4 = netcdf.inqVarID(ncid,'mu_s');
var5 = netcdf.inqVarID(ncid,'mu_d');
var6 = netcdf.inqVarID(ncid,'Dc');
var7 = netcdf.inqVarID(ncid,'C0');
var8 = netcdf.inqVarID(ncid,'dTx0');
var9 = netcdf.inqVarID(ncid,'dTy0');
var10 = netcdf.inqVarID(ncid,'dTz0');

netcdf.putVar(ncid,var1,Tx);
netcdf.putVar(ncid,var2,Ty);
netcdf.putVar(ncid,var3,Tz);
netcdf.putVar(ncid,var4,mu_s);
netcdf.putVar(ncid,var5,mu_d);
netcdf.putVar(ncid,var6,Dc);
netcdf.putVar(ncid,var7,C0);
netcdf.putVar(ncid,var8,dTx);
netcdf.putVar(ncid,var9,dTy);
netcdf.putVar(ncid,var10,dTz);

netcdf.close(ncid);

idx = find(nx_dense>1e-6);
if check_stress
    trisurf(tri_dense(idx,:),x_dense',y_dense',z_dense',v_dense',...
        'FaceColor','interp')
    colorbar
    axis equal
    view(30,30)
    axis image
    shading interp
end
