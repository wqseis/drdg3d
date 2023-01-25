function write_mesh_nc(fnm,mesh)

ncid = netcdf.create(fnm,'CLOBBER');
dimid_node = netcdf.defDim(ncid,'Nnode',mesh.Nnode);
dimid_elem = netcdf.defDim(ncid,'Nelem',mesh.Nelem);
%dimid_face_fault = netcdf.defDim(ncid,'num_face_fault',mesh.num_face_fault);
dimid_fault_elem = netcdf.defDim(ncid,'nfault_elem',mesh.nfault_elem);
%dimid_nfault_face = netcdf.defDim(ncid,'nfault_face',mesh.nfault_face);
dimid_node_fault = netcdf.defDim(ncid,'Nnode_fault',mesh.Nnode_fault);

dimid2 = netcdf.defDim(ncid,'two',2);
dimid3 = netcdf.defDim(ncid,'three',3);
dimid4 = netcdf.defDim(ncid,'four',4);

varid1 = netcdf.defVar(ncid,'node','NC_DOUBLE',[dimid3,dimid_node]);
varid2 = netcdf.defVar(ncid,'elem','NC_INT',[dimid4,dimid_elem]);
varid21 = netcdf.defVar(ncid,'node_fault','NC_INT',[dimid_node_fault]);
varid3 = netcdf.defVar(ncid,'neighbor','NC_INT',[dimid4,dimid_elem]);
varid4 = netcdf.defVar(ncid,'face','NC_INT',[dimid4,dimid_elem]);
varid5 = netcdf.defVar(ncid,'direction','NC_INT',[dimid4,dimid_elem]);
varid6 = netcdf.defVar(ncid,'bctype','NC_INT',[dimid4,dimid_elem]);
varid7 = netcdf.defVar(ncid,'elemtype','NC_INT',[dimid_elem]);
varid8 = netcdf.defVar(ncid,'rho','NC_DOUBLE',[dimid_elem]);
varid9 = netcdf.defVar(ncid,'vp','NC_DOUBLE',[dimid_elem]);
varid10 = netcdf.defVar(ncid,'vs','NC_DOUBLE',[dimid_elem]);
varid11 = netcdf.defVar(ncid,'part','NC_INT',[dimid_elem]);
if 1
% fault congfiguration
varid31 = netcdf.defVar(ncid,'wave2fault','NC_INT',[dimid_elem]);
varid32 = netcdf.defVar(ncid,'fault2wave','NC_INT',[dimid_fault_elem]);
varid12 = netcdf.defVar(ncid,'Tx0','NC_DOUBLE' ,[dimid3,dimid4,dimid_fault_elem]);
varid13 = netcdf.defVar(ncid,'Ty0','NC_DOUBLE' ,[dimid3,dimid4,dimid_fault_elem]);
varid14 = netcdf.defVar(ncid,'Tz0','NC_DOUBLE' ,[dimid3,dimid4,dimid_fault_elem]);
varid15 = netcdf.defVar(ncid,'mu_s','NC_DOUBLE',[dimid3,dimid4,dimid_fault_elem]);
varid16 = netcdf.defVar(ncid,'mu_d','NC_DOUBLE',[dimid3,dimid4,dimid_fault_elem]);
varid17 = netcdf.defVar(ncid,'Dc','NC_DOUBLE'  ,[dimid3,dimid4,dimid_fault_elem]);
varid18 = netcdf.defVar(ncid,'C0','NC_DOUBLE'  ,[dimid3,dimid4,dimid_fault_elem]);
end
netcdf.endDef(ncid);

netcdf.putVar(ncid,varid1,mesh.node);
netcdf.putVar(ncid,varid2,mesh.elem);
netcdf.putVar(ncid,varid21,mesh.node_fault);
netcdf.putVar(ncid,varid3,mesh.neighbor);
netcdf.putVar(ncid,varid4,mesh.face);
netcdf.putVar(ncid,varid5,mesh.direction);
netcdf.putVar(ncid,varid6,mesh.bctype);
netcdf.putVar(ncid,varid7,mesh.elemtype);
netcdf.putVar(ncid,varid8,mesh.rho);
netcdf.putVar(ncid,varid9,mesh.vp);
netcdf.putVar(ncid,varid10,mesh.vs);
netcdf.putVar(ncid,varid11,mesh.part);
netcdf.putVar(ncid,varid31,mesh.wave2fault);
netcdf.putVar(ncid,varid32,mesh.fault2wave);
netcdf.close(ncid);
