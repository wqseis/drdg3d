%clc
clear
%close all

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
a = zeros(3,4,nfault_elem);
b = zeros(3,4,nfault_elem);
state = zeros(3,4,nfault_elem);
Vw = zeros(3,4,nfault_elem);
TP_hy = zeros(3,4,nfault_elem);

temp = zeros(3,4,nfault_elem);

bctype = ncread(fnm_out,'bctype');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu_s1 = 0.18; % not used
mu_d1 = 0.12;
C01 = 0.0;
Dc1 = 0.4;
W = 15; w = 3;
V0 = 1e-6;
Vini = 1e-16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for ief = 1:nfault_elem
    ie = fault2wave(ief);
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

                T1 = [-120 -40 0;
                    -40 0 0;
                    0 0 0];

                vec_n = [nx(is,ie);ny(is,ie);nz(is,ie)];

                traction0 = T1 * vec_n;

                tn = dot(traction0,vec_n);
                ts_vec = traction0-vec_n*tn;
                ts = norm(ts_vec);
                tp = mu_s1*(-tn);
                td = mu_d1*(-tn);

                % add stress asperity to nucleate
                y0=0;
                z0=-7.5;
                %z0=0;
                dist = sqrt((y-y0)^2+ ((z-z0 )/sind(90))^2);
                R1=2;
                R2=3;

                R = 3e0;
                r = dist;
                if(r<R)
                    Fr = exp(r^2/(r^2-R^2));
                else
                    Fr = 0;
                end

                if(dist<R1)
                    asp_ratio=1;
                elseif(dist<R2)
                    asp_ratio=exp(-2*(dist-R1)^2/(R2-R1)^2);
                else
                    asp_ratio=0;
                end

                ts_vec_one = ts_vec/max(ts,1e-9);

                ts_vec_amp = asp_ratio*((...
                    1.001*(-tn)*mu_s1+C01)-ts)+ts;

                ts0 = ts;
                ts_vec_amp = ts + Fr * 45;

                ts_vec = ts_vec_amp * ts_vec_one;
                ts=norm(ts_vec);

                tp = mu_s1*(-tn); % peak stress
                td = mu_d1*(-tn);

                traction0 = ts_vec+tn*vec_n;

                Tx(i,is,ief) = traction0(1);
                Ty(i,is,ief) = traction0(2);
                Tz(i,is,ief) = traction0(3);

                mu_s(i,is,ief) = mu_s1;
                mu_d(i,is,ief) = mu_d1;
                Dc(i,is,ief) = Dc1;
                C0(i,is,ief) = C01;

                By = Bfunc(y-0, W, w);
                Bz = Bfunc(z+7.5, W/2, w);
                B = (1-By.*Bz);
                %a1 = 0.01+0.01*B;
                a1 = 0.01*(1+B);
                b1 = 0.014;
                Vw1 = 0.1+0.9*B;

                state1 = a1.*log(2.0*V0/Vini*sinh(ts0./abs(tn)./a1));

                a(i,is,ief) = a1;
                b(i,is,ief) = b1;
                Dc(i,is,ief) = 0.4;
                state(i,is,ief) = state1;
                Vw(i,is,ief) = Vw1;


                % for checking
                %temp(i,is,ief) = (tp-ts)./(ts-td);
                temp(i,is,ief) = (ts-td)./((tp-td)+1e-6);
                temp(i,is,ief) = -tn;
            end

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


%% checking ...
% The following code is unnecessary, but useful for checking
% initial values on the fault surface

% counting fault faces ...
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
v_dense = zeros(nfault,3);
x_dense = zeros(nfault,3);
y_dense = zeros(nfault,3);
z_dense = zeros(nfault,3);
nx_dense = zeros(nfault,1);
k = 1;

for ief = 1:nfault_elem
    ie = fault2wave(ief);
    for is = 1:4
        if (bctype(is,ie)>=BC_FAULT )
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

idx = find(nx_dense>1e-6); % extract the minus side of the fault surface
%idx = 1:size(nx_dense,1);
if check_stress
    figure
    trisurf(tri_dense(idx,:),x_dense',y_dense',z_dense',v_dense',...
        'FaceColor','interp')
    colorbar
    %view(30,30)
    axis image
    shading interp
end
