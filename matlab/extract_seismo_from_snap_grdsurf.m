function [v,t,x,y,z,idx2,idx3] = extract_seismo_from_snap_grdsurf(data_dir,nproc,loc,varnm)
x0=loc(1);
y0=loc(2);
z0=loc(3);
rmin = 1e30;
for iproc = 0:nproc-1
    fnm = sprintf('%s/grdsurf_mpi%06d.nc',data_dir,iproc);
    if (exist(fnm,'file'))
        x1 = ncread(fnm,'x');
        y1 = ncread(fnm,'y');
        z1 = ncread(fnm,'z');

        [Nfp,Nelem]=size(x1);
        
        r = sqrt( (x1(:)-x0).^2 + ( y1(:)-y0 ).^2+ ( z1(:) - z0 ) .^2 );
        [r1,idx1]=min(r(:));

        if (r1<rmin)
            rmin = r1;
            mpi_id = iproc;
            nflt = Nelem;
            
            idx = idx1;
            
            idx2 = mod(idx-1,Nfp)+1;
            idx3 = (idx-idx2)/Nfp+1;
            
            %fprintf("%d,y = %g, %g\n",iproc,y1(idx),y1(idx2,idx3))
            %fprintf("%d,z = %g, %g\n",iproc,z1(idx),z1(idx2,idx3))
        end
    end
end
clear x1 y1 z1 nx1
 
%idx2 = mod(idx-1,Nfp)+1;
%idx1 = (idx-idx2)/Nfp+1;
 
for iproc = mpi_id
    fnm = sprintf('%s/grdsurf_mpi%06d.nc',data_dir,iproc);
    
    x1 = ncread(fnm,'x');
    y1 = ncread(fnm,'y');
    z1 = ncread(fnm,'z');
    
    v = ncread(fnm, varnm ,[idx2,idx3,1],[1 1 Inf],[1 1 1]);
    t = ncread(fnm, 'time');
    
    fprintf('%d,x = %g, %g\n',iproc,x1(idx),x1(idx2,idx3))
    fprintf('%d,y = %g, %g\n',iproc,y1(idx),y1(idx2,idx3))
    fprintf('%d,z = %g, %g\n',iproc,z1(idx),z1(idx2,idx3))
end

v = squeeze(v);
t = squeeze(t);
 
end
