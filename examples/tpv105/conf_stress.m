%clc
clear
close all

%addpath('../matlab');
addmypath;

check_stress = 1;

%BC_FAULT = 100;

%FtoV = [1,3,2;1,2,4;2,3,4;1,4,3];% point outward


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
a = zeros(3,4,nfault_elem);
b = zeros(3,4,nfault_elem);
state = zeros(3,4,nfault_elem);
Vw = zeros(3,4,nfault_elem);
TP_hy = zeros(3,4,nfault_elem);

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

disp('assigning fault triangulates ...')
tic
tri_fault = zeros(num_fault_face,3);
k = 1;
for ie = 1:nelem
    for is = 1:4
        if (bctype(is,ie)>=BC_FAULT)
            tri =  (elem(FtoV(is,:),ie));
            tri_fault(k,:) = tri;
            k = k + 1;
            %fprintf('%g %g %g\n',nx(is,ie),ny(is,ie),nz(is,ie));
        end
    end
end
toc

W = 15e0; w = 3e0;

W = 4; w = 3;
W = 15; w= 3;

V0 = 1e-6;
Vini = 1e-16;
% By = Bfunc(y, W, w);
% Bz = Bfunc(z+7.5e3, W/2, w);
% B = (1-By.*Bz);
% a = 0.01+0.01*B;
 
mu_s1 = 0.677;
mu_d1 = 0.525;
C01 = 0;
Dc1 = 0.4;


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
                
%                 depth = -z;
%                 if(depth<0.1)
%                     depth=0.1;
%                 end
%                 
%                 
%                 if(depth<20)
%                     sigma_v = depth/20*100;
%                 else
%                     sigma_v = 100;
%                 end
%                 
%                 sigma_v = 20;
%           
%                 
%                 sigma_v =  -sigma_v;
%                 sigma_H = 1*sigma_v;
%                 sigma_h = 0.5*sigma_v;
%                 
%                 T1 = diag([sigma_H,sigma_h,sigma_v]);
                
                sigma_ini = max(1.67*9.8*z,-45);
                % small sigma is prone to instablity
                sigma_ini = min(sigma_ini,-1e-3);

                %sigma_ini = -45;
                tau_ini = 0.41 * sigma_ini;
                
                T1 = [-120 -45 0;
                    -45 0 0;
                    0 0 0];
                T1 = [sigma_ini tau_ini,0;
                    tau_ini,0,0;
                    0 0 0];
%                 if (abs(yc-0)<=1.5+1e-6 & abs(zc-0)<=1.5+1e-6)
%                     T1 = [-120 -100 0;
%                     -100 0 0;
%                     0 0 0];
%                     
%                 end

             
                % Coordinate transforms
                % rotations of the coordinate system 
%                 % while the object is held constant. 
%                 angle = 0;
%                 sina=sind(angle);cosa=cosd(angle);
%                 Rot= [ cosa sina 0;
%                       -sina cosa 0;
%                           0    0 1];
%                 
%                 T2=Rot*T1*transpose(Rot);
%                 
                vec_n = [nx(is,ie);ny(is,ie);nz(is,ie)];
                
                
                
                traction0 = T1 * vec_n;
                
                tn=dot(traction0,vec_n);
                ts_vec=traction0-tn*vec_n;
                ts=norm(ts_vec);
                
                % add stress asperity to nucleate
                y0 = -4;
                z0 = -7.5;
                %z0 = -1.5;
                %z0 = -5;
                %z0=0;
                dist = sqrt((y-y0)^2+ ((z-z0 )/sind(90))^2);
                R1=2;
                R2=3;

                R = 1.5e0;
                %R = 3;
                r = dist;
                if(r<R)
                    Fr = exp(r^2/(r^2-R^2));
                else
                    Fr = 0;
                end
%{
                if(dist<R1)
                    asp_ratio=1;
                elseif(dist<R2)
                    asp_ratio=exp(-2*(dist-R1)^2/(R2-R1)^2);
                else
                    asp_ratio=0;
                end
                
                ts_vec_one = ts_vec/max(ts,1e-9);
                
                mu_s1 = 0.677;
                mu_d1 = 0.525;
              %  C01 = 0.2;
                C01 = 0;
                Dc1 = 0.4;
                
                ts_vec_amp = asp_ratio*((...
                    1.001*(-tn)*mu_s1+C01)-ts)+ts;
%}
                ts0 = ts; % save background intial stress
                ts_vec_amp = ts + Fr * 0;
                ts_vec_one = ts_vec/max(ts,1e-9);
                ts_vec = ts_vec_amp * ts_vec_one;
                ts=norm(ts_vec);
                traction0 = ts_vec+tn*vec_n;

                Tx(i,is,ief) = traction0(1);
                Ty(i,is,ief) = traction0(2);
                Tz(i,is,ief) = traction0(3);

                ts_vec_amp = Fr * 50;
                ts_vec_one = ts_vec/max(ts,1e-9);
                ts_vec = ts_vec_amp * ts_vec_one;
                ts=norm(ts_vec);
                traction0 = ts_vec+0*vec_n;

                dTx(i,is,ief) = traction0(1);
                dTy(i,is,ief) = traction0(2);
                dTz(i,is,ief) = traction0(3);
                
                mu_s(i,is,ief) = mu_s1;
                mu_d(i,is,ief) = mu_d1;
                Dc(i,is,ief) = Dc1;
                C0(i,is,ief) = C01;

                By = Bfunc(y-0, W, w);
                %Bz = Bfunc(z+7.5, W/2, w);
                Bz = Bfunc_free( 0-z,W,w);
                Bz3 = Bfunc(0-z, W, w);
                B = (1-By.*Bz);
                BB = (1-By.*Bz3);
                %a1 = 0.01+0.01*B;
                a1 = 0.01*(1+B);
                b1 = 0.014;
                Vw1 = 0.1+0.9*B;
                hy1 = 4e-4+1*BB;

                state1 = a1.*log(2.0*V0/Vini*sinh(ts0./abs(tn)./a1));
                
                a(i,is,ief) = a1;
                b(i,is,ief) = b1;
                Dc(i,is,ief) = 0.4;
                state(i,is,ief) = state1;
                Vw(i,is,ief) = Vw1;
                TP_hy(i,is,ief) = hy1;

                temp(i,is,ief) = Tx(i,is,ief)+dTx(i,is,ief);
                %temp(i,is,ief) = (tp-ts)./(ts-td);
                %temp(i,is,ief) = (ts-td)/(tp-td);
                %temp(i,is,ief) = hy1;
                %temp(i,is,ief) = ts0+ts;
                %temp(i,is,ief) = Vw1;
                %temp(i,is,ief) = a1;
                %temp(i,is,ief) = state1;
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
var13 = netcdf.inqVarID(ncid,'dTx0');
var14 = netcdf.inqVarID(ncid,'dTy0');
var15 = netcdf.inqVarID(ncid,'dTz0');
var4 = netcdf.inqVarID(ncid,'mu_s');
var5 = netcdf.inqVarID(ncid,'mu_d');
var6 = netcdf.inqVarID(ncid,'Dc');
var7 = netcdf.inqVarID(ncid,'C0');
var8 = netcdf.inqVarID(ncid,'a');
var9 = netcdf.inqVarID(ncid,'b');
var10 = netcdf.inqVarID(ncid,'Vw');
var11 = netcdf.inqVarID(ncid,'state');
var12 = netcdf.inqVarID(ncid,'TP_hy');

netcdf.putVar(ncid,var1,Tx);
netcdf.putVar(ncid,var2,Ty);
netcdf.putVar(ncid,var3,Tz);
netcdf.putVar(ncid,var13,dTx);
netcdf.putVar(ncid,var14,dTy);
netcdf.putVar(ncid,var15,dTz);
netcdf.putVar(ncid,var4,mu_s);
netcdf.putVar(ncid,var5,mu_d);
netcdf.putVar(ncid,var6,Dc);
netcdf.putVar(ncid,var7,C0);
netcdf.putVar(ncid,var8,a);
netcdf.putVar(ncid,var9,b);
netcdf.putVar(ncid,var10,Vw);
netcdf.putVar(ncid,var11,state);
netcdf.putVar(ncid,var12,TP_hy);

netcdf.close(ncid);
%%
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
